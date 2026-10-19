"""Build the reel the web demo spins.

Components are the ones decomposed out of PROTAC-8K by
``scripts/decompose_corpus.py``, taken in order of how often they actually
appear in reported degraders. The demo shows the head of that distribution —
enough to spin, few enough to precompute every product — while the funnel it
reports is measured over the whole library by ``scripts/enumerate_space.py``.

Names are assigned by hand from structure, for the warheads and E3 ligands
whose chemotype is unambiguous. They are *not* part of the dataset, which
ships structures only. Anything not confidently identifiable keeps a neutral
label rather than a guess; a wrong drug name is worse than no name.
"""

from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rdkit import Chem, RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

from lasprotacs.assemble import assemble  # noqa: E402
from lasprotacs.rules import attachment_ok, attachment_reason, new_alerts  # noqa: E402
from lasprotacs.score import developability, flags, properties  # noqa: E402

N_WARHEADS, N_LINKERS, N_E3 = 14, 22, 8

# Hand assignments, keyed by canonical SMILES of the decomposed fragment.
# Only chemotypes that identify a compound or a well-defined inhibitor class
# on sight.
WARHEAD_NAMES = {
    "Nc1ncnc2c1c(-c1ccc(Oc3ccccc3)cc1)nn2[*:1]": ("Ibrutinib", "BTK"),
    "O=c1[nH]nc(Cc2ccc(F)c([*:1])c2)c2ccccc12": ("Olaparib (phthalazinone)", "PARP1"),
    "CC(=O)c1c(C)c2cnc(Nc3ccc([*:1])cn3)nc2n(C2CCCC2)c1=O": ("Palbociclib", "CDK4/6"),
    "Cc1sc2c(c1C)C(c1ccc(Cl)cc1)=N[C@@H]([*:1])c1nnc(C)n1-2": ("JQ1", "BRD4"),
    "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1[*:1]": ("4-anilinoquinazoline", "EGFR"),
    "Cc1cc(Nc2ncc(Cl)c(Nc3ccccc3S(=O)(=O)C(C)C)n2)c(OC(C)C)cc1[*:1]": (
        "Ceritinib", "ALK"),
    "COc1cc([*:1])ccc1Nc1ncc(Cl)c(Nc2ccccc2P(C)(C)=O)n1": ("Brigatinib", "ALK"),
    "Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc([*:1])n1": ("Dasatinib", "ABL/SRC"),
    "Fc1cc(I)ccc1Nc1c([*:1])ccc(F)c1F": ("Diarylamine MEK inhibitor", "MEK1/2"),
    "CCN(c1cc(-c2ccc([*:1])cc2)cc(C(=O)NCc2c(C)cc(C)[nH]c2=O)c1C)C1CCOCC1": (
        "Tazemetostat-like", "EZH2"),
    "O=C(NO)c1cccc(C(=O)NCCc2ccc([*:1])cc2)c1": ("Hydroxamic acid", "HDAC"),
}

E3_NAMES = {
    "O=C1CC[C@@H](N2C(=O)c3cccc([*:2])c3C2=O)C(=O)N1": "Pomalidomide 4-substituted",
    "O=C1CC[C@H](N2C(=O)c3cccc([*:2])c3C2=O)C(=O)N1": "Pomalidomide 4-subst. (R)",
    "O=C1CC[C@@H](N2C(=O)c3ccc([*:2])cc3C2=O)C(=O)N1": "Thalidomide 5-substituted",
    "O=C1CC[C@H](N2C(=O)c3ccc([*:2])cc3C2=O)C(=O)N1": "Thalidomide 5-subst. (R)",
    "O=C1CC[C@@H](N2Cc3c(cccc3[*:2])C2=O)C(=O)N1": "Lenalidomide 4-substituted",
    "O=C1CC[C@H](N2Cc3c(cccc3[*:2])C2=O)C(=O)N1": "Lenalidomide 4-subst. (R)",
    "O=C1CC[C@@H](N2Cc3cc([*:2])ccc3C2=O)C(=O)N1": "Lenalidomide 5-substituted",
    "Cc1ncsc1-c1ccc(CNC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](N[*:2])C(C)(C)C)cc1": (
        "VH032"),
    "Cc1ncsc1-c1ccc([C@H](C)NC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](N[*:2])C(C)(C)C)cc1": (
        "VH101 (α-methyl)"),
    "Cc1ncsc1-c1ccc(CNC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@H](N[*:2])C(C)(C)C)cc1": (
        "VH032 epimer"),
}


def linker_name(smiles: str) -> tuple[str, str]:
    """Describe a linker from its structure: chemists read these at a glance."""
    mol = Chem.MolFromSmiles(smiles)
    heavy = mol.GetNumHeavyAtoms()
    n_o = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "O" and not any(
        b.GetBondType() == Chem.BondType.DOUBLE for b in a.GetBonds()))
    amide = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[NX3]"))
    ring = mol.GetRingInfo().NumRings()
    triazole = mol.HasSubstructMatch(Chem.MolFromSmarts("c1cn[nH0]n1"))

    if triazole:
        return f"triazole, {heavy} atoms", "click"
    if ring:
        return f"ring spacer, {heavy} atoms", "rigid"
    if n_o >= 2:
        tag = "PEG-amide" if amide else "PEG"
        return f"PEG{n_o}{', amide' if amide else ''} ({heavy} atoms)", tag
    if amide:
        return f"alkyl amide, {heavy} atoms", "alkyl-amide"
    return f"alkyl, {heavy} atoms", "alkyl"


def take(items, n, namer) -> list[dict]:
    out = []
    for c in items[:n]:
        name, tag = namer(c)
        out.append({
            "name": name, "tag": tag, "smiles": c["smiles"],
            "n": c["n"], "heavy": c["heavy"],
        })
    return out


def main() -> None:
    lib_path = ROOT / "data" / "components.json"
    if not lib_path.exists():
        sys.exit("no data/components.json — run scripts/decompose_corpus.py first")
    lib = json.loads(lib_path.read_text())

    def wh_namer(c):
        name, target = WARHEAD_NAMES.get(c["smiles"], (None, None))
        if name:
            return f"{name}", target
        return f"warhead ({c['heavy']} atoms)", "unassigned"

    def e3_namer(c):
        return E3_NAMES.get(c["smiles"], f"{c['family']} ligand"), c["family"]

    warheads = take(lib["warheads"], N_WARHEADS, wh_namer)
    linkers = take(lib["linkers"], N_LINKERS, lambda c: linker_name(c["smiles"]))
    e3s = take(lib["e3s"], N_E3, e3_namer)

    named = sum(1 for w in warheads if w["tag"] != "unassigned")
    print(f"reel: {len(warheads)} warheads ({named} named), "
          f"{len(linkers)} linkers, {len(e3s)} E3 ligands")

    raw = len(warheads) * len(linkers) * len(e3s)
    combos: dict[str, dict] = {}
    rejects: dict[str, str] = {}
    n_attach = n_built = n_clean = 0
    for (wi, w), (li, l), (ei, e) in product(
        enumerate(warheads), enumerate(linkers), enumerate(e3s)
    ):
        if not attachment_ok(w["smiles"], l["smiles"], e["smiles"]):
            rejects[f"{wi}-{li}-{ei}"] = attachment_reason(
                w["smiles"], l["smiles"], e["smiles"]
            )
            continue
        n_attach += 1
        smi = assemble(w["smiles"], l["smiles"], e["smiles"])
        if smi is None:
            continue
        n_built += 1
        props = properties(smi)
        if props is None:
            continue
        alerts = new_alerts(smi, [w["smiles"], l["smiles"], e["smiles"]])
        if not alerts:
            n_clean += 1
        combos[f"{wi}-{li}-{ei}"] = {
            "smiles": smi, "mw": props.mw, "clogp": props.clogp,
            "tpsa": props.tpsa, "hbd": props.hbd, "rotb": props.rotb,
            "dev": developability(props), "alerts": alerts, "flags": flags(props),
        }

    print(f"  raw {raw} -> attachment {n_attach} -> assembled {n_built} "
          f"-> no new alerts {n_clean}")

    funnel_path = ROOT / "data" / "funnel.json"
    full = json.loads(funnel_path.read_text()) if funnel_path.exists() else {}

    payload = {
        "warheads": warheads, "linkers": linkers, "e3s": e3s, "combos": combos,
        "rejects": rejects,
        "funnel": {"raw": raw, "attachment": n_attach,
                   "assembled": n_built, "clean": n_clean},
        "library": {
            "warheads": len(lib["warheads"]), "linkers": len(lib["linkers"]),
            "e3s": len(lib["e3s"]), "degraders": lib["degraders"],
            "space": lib["space"],
        },
        "full_funnel": full,
    }
    out = ROOT / "docs" / "reel.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, separators=(",", ":")))
    print(f"wrote {out}  ({out.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
