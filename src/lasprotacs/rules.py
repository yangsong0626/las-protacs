"""Chemistry rules for which pieces may be joined, and what must be rejected.

A naive cross-product of warheads x linkers x E3 ligands produces bonds no
chemist would draw.  Joining a VHL ligand (exit vector: an amine nitrogen) to
a piperazine linker (exit vector: a ring nitrogen) makes a hydrazine.  Every
such product is a wasted slot in the output, and a reviewer notices
immediately.

So the reel is filtered twice:

1. ``attachment_ok`` — before assembly, on the pair of atoms about to be
   bonded.  Cheap, and removes most of the nonsense.
2. ``structure_alerts`` — after assembly, on the whole molecule, for motifs
   that only appear once the pieces meet.
"""

from __future__ import annotations

from rdkit import Chem

HETERO = {7, 8, 16}  # N, O, S

# Motifs that are unstable, synthetically awkward, or usually unwanted.
#
# These are counted on the product *and* on the pieces it came from, and only
# the difference is reported — a warhead that already contains a sulfonamide
# or an acetal is the chemist's business, not the recombiner's.  Oxidised
# nitrogen and sulfur (nitro, N-oxide, sulfonamide, sulfinamide) are excluded
# explicitly so they are never mistaken for hydrazines.
_NOT_OXIDISED = "!$([#7,#16]~[OX1]);!$([#7,#16][#16X4]);!$([#7,#16][#16X3]);!$([#7+]);!$([#8-])"
_ALERTS: list[tuple[str, str]] = [
    ("hetero-hetero single bond", f"[NX3,OX2,SX2;{_NOT_OXIDISED}]-[NX3,OX2,SX2;{_NOT_OXIDISED}]"),
    ("anhydride", "[CX3](=O)[OX2][CX3](=O)"),
    ("N,N-acetal", "[NX3][CX4H1,CX4H2][NX3]"),
    ("N,O-acetal", "[NX3][CX4H1,CX4H2][OX2]"),
    ("O,O-acetal (acyclic)", "[OX2;R0][CX4H1,CX4H2;R0][OX2;R0]"),
    ("imine (acyclic)", "[CX3;R0]=[NX2;R0][#6]"),
    ("gem-diol", "[OX2H1][CX4][OX2H1]"),
]
ALERTS = [(name, Chem.MolFromSmarts(sma)) for name, sma in _ALERTS]


def attachment_atoms(smiles: str) -> dict[int, int]:
    """Map attachment label -> atomic number of the atom carrying the dummy."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    out: dict[int, int] = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        nbrs = atom.GetNeighbors()
        if not nbrs:
            continue
        out[atom.GetAtomMapNum()] = nbrs[0].GetAtomicNum()
    return out


def bond_ok(a: int, b: int) -> bool:
    """Would a bond between elements ``a`` and ``b`` be a sensible one to make?

    At least one partner must be carbon.  Heteroatom-heteroatom single bonds
    formed by recombination are hydrazines, hydroxylamines, peroxides and
    sulfenamides; none belong in a degrader design unless deliberately placed.
    """
    if a == 0 or b == 0:
        return False
    return a == 6 or b == 6


def attachment_ok(warhead: str, linker: str, e3: str) -> bool:
    wh = attachment_atoms(warhead)
    lk = attachment_atoms(linker)
    e3a = attachment_atoms(e3)
    if 1 not in wh or 2 not in e3a or 1 not in lk or 2 not in lk:
        return False
    return bond_ok(wh[1], lk[1]) and bond_ok(e3a[2], lk[2])


SYMBOL = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 15: "P", 16: "S"}

# What you actually get when two heteroatoms are bonded together.
_PRODUCT = {
    frozenset("NN"): "a hydrazine",
    frozenset("NO"): "a hydroxylamine",
    frozenset("OO"): "a peroxide",
    frozenset("NS"): "a sulfenamide",
    frozenset("OS"): "a sulfenate ester",
    frozenset("SS"): "a disulfide",
}


def attachment_reason(warhead: str, linker: str, e3: str) -> str | None:
    """Why this combination was refused, in words, or None if it was not.

    The UI shows this verbatim, so it names the bond and what it would make.
    """
    wh = attachment_atoms(warhead)
    lk = attachment_atoms(linker)
    e3a = attachment_atoms(e3)
    if 1 not in wh or 2 not in e3a or 1 not in lk or 2 not in lk:
        return "a piece is missing an attachment point"
    for a, b, side in ((wh[1], lk[1], "warhead-linker"),
                       (e3a[2], lk[2], "linker-E3")):
        if bond_ok(a, b):
            continue
        sa, sb = SYMBOL.get(a, "?"), SYMBOL.get(b, "?")
        made = _PRODUCT.get(frozenset(sa + sb), "a heteroatom-heteroatom bond")
        return f"{side} join is {sa}–{sb}: that makes {made}"
    return None


def _alert_counts(mol_or_smiles) -> dict[str, int]:
    mol = (
        Chem.MolFromSmiles(mol_or_smiles)
        if isinstance(mol_or_smiles, str)
        else mol_or_smiles
    )
    if mol is None:
        return {"unparseable": 1}
    counts: dict[str, int] = {}
    for name, patt in ALERTS:
        if patt is None:
            continue
        n = len(mol.GetSubstructMatches(patt))
        if n:
            counts[name] = n
    return counts


def structure_alerts(mol_or_smiles) -> list[str]:
    """Names of unwanted motifs present anywhere in the molecule."""
    return sorted(_alert_counts(mol_or_smiles))


def new_alerts(product, pieces) -> list[str]:
    """Motifs the product has that its pieces did not — i.e. created by the join.

    A sulfonamide or acetal that arrived inside a published warhead is not a
    defect of the recombination, so only the increase is reported.
    """
    prod = _alert_counts(product)
    if "unparseable" in prod:
        return ["unparseable"]
    before: dict[str, int] = {}
    for p in pieces:
        for k, v in _alert_counts(p).items():
            before[k] = before.get(k, 0) + v
    return sorted(k for k, v in prod.items() if v > before.get(k, 0))
