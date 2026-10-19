"""Check the hand-assigned names against the compounds they claim.

The dataset ships structures, not names. The demo labels its reels, so those
labels are a claim this repo is making on its own and they should be checked
rather than trusted.

Two kinds of claim, checked two different ways:

* **A compound name** ("Ibrutinib", "Palbociclib") claims the fragment is the
  recognisable core of that drug. Cap the attachment point with hydrogen and
  the result must be a substructure of the published parent. It will usually be
  a proper substructure: a degrader replaces part of the drug with the linker,
  which is the whole point, so an exact formula match is the wrong test.
* **A chemotype name** ("4-anilinoquinazoline", "Hydroxamic acid") claims only
  that a diagnostic motif is present, so it is checked with a SMARTS.

Run it after changing any name in ``build_demo_data.py``::

    python scripts/validate_components.py
"""

from __future__ import annotations

import sys
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_demo_data import E3_NAMES, WARHEAD_NAMES  # noqa: E402

# Published parents, drawn independently of the fragments they are checked
# against.
PARENTS = {
    "Ibrutinib":
        "C=CC(=O)N1CCC[C@@H](C1)n1nc(-c2ccc(Oc3ccccc3)cc2)c2c(N)ncnc21",
    # 2-{[5-(piperazin-1-yl)pyridin-2-yl]amino}: the piperazine sits on the
    # 5-position, para-ish to the ring nitrogen, not meta.  Drawing it wrong
    # here is what this script is for.
    "Palbociclib":
        "CC(=O)c1c(C)c2cnc(Nc3ccc(N4CCNCC4)cn3)nc2n(C2CCCC2)c1=O",
    "JQ1":
        "Cc1sc2c(c1C)C(c1ccc(Cl)cc1)=N[C@@H](CC(=O)O)c1nnc(C)n1-2",
    "Ceritinib":
        "Cc1cc(Nc2ncc(Cl)c(Nc3ccccc3S(=O)(=O)C(C)C)n2)c(OC(C)C)cc1C1CCNCC1",
    "Brigatinib":
        "COc1cc(N2CCC(N3CCN(C)CC3)CC2)ccc1Nc1ncc(Cl)c(Nc2ccccc2P(C)(C)=O)n1",
    "Dasatinib":
        "Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1",
    "Pomalidomide": "Nc1cccc2c1C(=O)N(C1CCC(=O)NC1=O)C2=O",
    "Thalidomide": "O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1",
    "Lenalidomide": "Nc1cccc2c1CN(C1CCC(=O)NC1=O)C2=O",
    "VH032":
        "Cc1ncsc1-c1ccc(CNC(=O)[C@@H]2C[C@@H](O)CN2C(=O)"
        "[C@@H](NC(C)=O)C(C)(C)C)cc1",
    "VH101":
        "Cc1ncsc1-c1ccc([C@H](C)NC(=O)[C@@H]2C[C@@H](O)CN2C(=O)"
        "[C@@H](NC(C)=O)C(C)(C)C)cc1",
}

# Chemotype claims: the motif that has to be there for the label to be fair.
CHEMOTYPES = {
    "4-anilinoquinazoline": "[#7]c1ncnc2ccccc12",
    "Olaparib (phthalazinone)": "O=c1[nH]nc(C)c2ccccc12",
    "Hydroxamic acid": "[OX2H1][NX3H1][CX3]=O",
    "Diarylamine MEK inhibitor": "c1ccccc1[NX3H1]c1ccccc1",
    "Tazemetostat-like": "O=c1[nH]c(C)cc(C)c1C[NX3H1]C=O",
}


def capped(fragment: str) -> Chem.Mol | None:
    """Replace the attachment dummy with hydrogen to recover a real molecule."""
    mol = Chem.MolFromSmiles(fragment)
    if mol is None:
        return None
    rw = Chem.RWMol(mol)
    for atom in [a for a in rw.GetAtoms() if a.GetAtomicNum() == 0]:
        rw.ReplaceAtom(atom.GetIdx(), Chem.Atom(1))
    try:
        out = Chem.RemoveHs(rw.GetMol())
        Chem.SanitizeMol(out)
    except Exception:
        return None
    return out


def flat(mol: Chem.Mol) -> Chem.Mol:
    """Stereochemistry is drawn inconsistently across sources; ignore it."""
    m = Chem.Mol(mol)
    Chem.RemoveStereochemistry(m)
    return m


def check(label: str, fragment: str) -> tuple[bool, str]:
    frag = capped(fragment)
    if frag is None:
        return False, "fragment will not parse"

    base = label.split(" (")[0].replace("-like", "").strip()
    base = base.replace(" 4-substituted", "").replace(" 5-substituted", "")
    base = base.replace(" 4-subst. (R)", "").replace(" 5-subst. (R)", "")
    base = base.replace(" epimer", "").replace(" (α-methyl)", "").strip()

    if label in CHEMOTYPES:
        patt = Chem.MolFromSmarts(CHEMOTYPES[label])
        if patt is None:
            return False, "bad chemotype SMARTS"
        ok = frag.HasSubstructMatch(patt)
        return ok, "chemotype present" if ok else "chemotype ABSENT"

    if base in PARENTS:
        parent = Chem.MolFromSmiles(PARENTS[base])
        if parent is None:
            return False, "parent will not parse"
        ok = flat(parent).HasSubstructMatch(flat(frag), useChirality=False)
        n_f, n_p = frag.GetNumHeavyAtoms(), parent.GetNumHeavyAtoms()
        if ok:
            return True, f"core of {base} ({n_f}/{n_p} heavy atoms)"
        return False, f"NOT a substructure of {base}"

    return True, "no claim to check (neutral label)"


def main() -> None:
    bad = 0
    # WARHEAD_NAMES maps smiles -> (name, target); E3_NAMES maps smiles -> name.
    for fragment, value in list(WARHEAD_NAMES.items()) + list(E3_NAMES.items()):
        label = value[0] if isinstance(value, tuple) else value
        ok, note = check(label, fragment)
        print(f"  {'ok ' if ok else 'FAIL'}  {label:34s} {note}")
        if not ok:
            bad += 1
    total = len(WARHEAD_NAMES) + len(E3_NAMES)
    print(f"\n{total - bad}/{total} name claims check out")
    if bad:
        sys.exit(1)


if __name__ == "__main__":
    main()
