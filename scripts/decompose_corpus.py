"""Decompose a corpus of real PROTACs and build the component library.

This is the validation step. The slot machine is only interesting if the
decomposition it rests on is right, and the only way to know that is to run it
over every published degrader we can get hold of and count.

Input is a SMILES file (``<smiles>\t<id>`` per line); ``scripts/prepare_protac8k.py``
writes one from the PROTAC-8K structure set. Output is
``data/components.json``: every warhead, linker and E3 ligand observed, with
the number of degraders each appears in.

Usage::

    python scripts/decompose_corpus.py data/protac8k/protac.smi
"""

from __future__ import annotations

import json
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rdkit import Chem, RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

from lasprotacs.assemble import assemble, same_molecule  # noqa: E402
from lasprotacs.decompose import decompose  # noqa: E402
from lasprotacs.e3_library import identify_e3  # noqa: E402


def load(path: Path) -> dict[str, str]:
    """Canonical SMILES -> first id seen, de-duplicated."""
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t") if "\t" in line else line.split()
        smi = parts[0]
        name = parts[1] if len(parts) > 1 else f"mol{len(out)}"
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        out.setdefault(Chem.MolToSmiles(mol), name)
    return out


def main() -> None:
    src = Path(sys.argv[1] if len(sys.argv) > 1 else "data/protac8k/protac.smi")
    if not src.exists():
        sys.exit(f"no corpus at {src} — run scripts/prepare_protac8k.py first")

    t0 = time.time()
    corpus = load(src)
    print(f"corpus              : {len(corpus)} unique structures from {src}")

    rows = []
    no_e3 = cut_failed = 0
    strict = cip_only = 0
    for smi, name in corpus.items():
        try:
            d = decompose(smi)
        except Exception:
            d = None
        if d is None:
            mol = Chem.MolFromSmiles(smi)
            if mol is None or identify_e3(mol) is None:
                no_e3 += 1
            else:
                cut_failed += 1
            continue
        back = assemble(d.warhead, d.linker, d.e3)
        if back == smi:
            strict += 1
            rt = "exact"
        elif same_molecule(back, smi):
            cip_only += 1
            rt = "stereo-relabelled"
        else:
            rt = "FAILED"
        rows.append({"id": name, "smiles": smi, "roundtrip": rt, **d.as_dict()})

    n = len(rows)
    broken = [r for r in rows if r["roundtrip"] == "FAILED"]
    fam = Counter(r["e3_family"] for r in rows)

    print(f"  no E3 ligand found: {no_e3}")
    print(f"  E3 found, no clean cut: {cut_failed}")
    print(f"decomposed          : {n}  ({100 * n / len(corpus):.1f}%)")
    print(f"  by ligase         : {dict(fam)}")
    print("round trip (decompose -> reassemble):")
    print(f"  identical SMILES  : {strict}  ({100 * strict / n:.1f}%)")
    print(f"  same molecule, ring-stereo relabelled: {cip_only}")
    print(f"  genuinely wrong   : {len(broken)}")

    lib = {}
    for field, label in (("warhead", "warheads"), ("linker", "linkers"), ("e3", "e3s")):
        counts = Counter(r[field] for r in rows)
        items = []
        for smi, k in counts.most_common():
            mol = Chem.MolFromSmiles(smi)
            entry = {"smiles": smi, "n": k, "heavy": mol.GetNumHeavyAtoms()}
            if field == "e3":
                fams = Counter(
                    r["e3_family"] for r in rows if r[field] == smi
                ).most_common(1)
                entry["family"] = fams[0][0] if fams else "?"
            items.append(entry)
        lib[label] = items
        cover10 = sum(i["n"] for i in items[:10]) / n
        print(f"{label:9s}: {len(items):5d} unique   top-10 cover {100 * cover10:.0f}%")

    space = len(lib["warheads"]) * len(lib["linkers"]) * len(lib["e3s"])
    print(f"\nformal combination space: {space:,}  ({space:.2e})")

    lib["source"] = str(src)
    lib["degraders"] = n
    lib["space"] = space
    out = ROOT / "data" / "components.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(lib, separators=(",", ":")))
    (ROOT / "data" / "decomposed.json").write_text(json.dumps(rows, separators=(",", ":")))
    print(f"wrote {out}  ({out.stat().st_size / 1024:.0f} KB)   [{time.time() - t0:.0f}s]")


if __name__ == "__main__":
    main()
