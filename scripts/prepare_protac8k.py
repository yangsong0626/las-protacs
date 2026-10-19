"""Turn the PROTAC-8K structure set into a SMILES corpus.

PROTAC-8K ships 3D structures as mol2, one file per degrader. Two things are
worth knowing before you point RDKit at them:

* **RDKit's mol2 reader rejects most of them.** The files carry Tripos ``ar``
  bond types, and on fused heteroaromatics RDKit's kekulisation fails with
  ``KekulizeException``. Lenient sanitisation flags do not rescue it. On a
  sample of three files, one parsed. Open Babel reads all 6,113.
* **So Open Babel does the reading and RDKit does the chemistry.** Each
  structure goes mol2 -> SMILES (Open Babel) -> canonical SMILES (RDKit).
  Round-tripped through canonical SMILES, 6,113 of 6,113 survive.

Usage::

    python scripts/prepare_protac8k.py ~/Downloads/PROTAC-8K/protac \\
        -o data/protac8k/protac.smi

Accepts a directory of ``*.mol2``, or a single concatenated multi-record mol2
(faster, and what you want if the files live on another machine).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    from openbabel import pybel
except ImportError:  # pragma: no cover
    sys.exit("needs Open Babel: pip install openbabel-wheel")

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
pybel.ob.obErrorLog.SetOutputLevel(0)


def records(src: Path):
    """Yield (title, obabel molecule) from a file or a directory of mol2."""
    if src.is_dir():
        files = sorted(
            src.glob("*.mol2"),
            key=lambda p: (len(p.stem), p.stem),  # 1, 2, ... 10 rather than 1, 10
        )
        for f in files:
            try:
                yield f.stem, next(pybel.readfile("mol2", str(f)))
            except Exception:
                yield f.stem, None
    else:
        for i, mol in enumerate(pybel.readfile("mol2", str(src))):
            yield (mol.title.strip() or f"mol{i}"), mol


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("src", type=Path, help="directory of *.mol2, or one multi-record mol2")
    ap.add_argument("-o", "--out", type=Path, default=Path("data/protac8k/protac.smi"))
    args = ap.parse_args()
    if not args.src.exists():
        sys.exit(f"not found: {args.src}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    n = ok = ob_fail = rd_fail = 0
    seen: set[str] = set()
    with args.out.open("w") as fh:
        for title, mol in records(args.src):
            n += 1
            if mol is None:
                ob_fail += 1
                continue
            try:
                parts = mol.write("smi").split()
            except Exception:
                parts = []
            if not parts:
                ob_fail += 1
                continue
            rd = Chem.MolFromSmiles(parts[0])
            if rd is None:
                rd_fail += 1
                continue
            smi = Chem.MolToSmiles(rd)
            ok += 1
            if smi in seen:
                continue
            seen.add(smi)
            fh.write(f"{smi}\t{title}\n")

    print(f"read            : {n}")
    print(f"  Open Babel failed : {ob_fail}")
    print(f"  RDKit rejected    : {rd_fail}")
    print(f"converted       : {ok} ({100 * ok / max(n, 1):.1f}%)")
    print(f"unique written  : {len(seen)} -> {args.out}")


if __name__ == "__main__":
    main()
