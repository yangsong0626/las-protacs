"""Split a PROTAC into warhead / linker / E3-ligand with labelled exit vectors.

The cut is *anchored*, not searched.  A free search over bond pairs will
happily shave the tert-leucine off a VHL ligand or hand the whole alkyl chain
to the E3 arm, because a short linker always looks more "linker-like" than a
long one.  So instead:

1. Match the whole E3 ligand (``e3_library.ligand_smarts``) and cut the one
   bond that leaves it.  The ligand boundary is chemistry, not a free
   parameter.
2. Walk outwards from that exit vector through chain atoms and saturated
   rings.  Stop at the first aromatic system — that is where the warhead's
   scaffold begins.  Everything walked through is the linker.

Attachment points come back as atom-mapped dummies so the pieces reassemble
with ``Chem.molzip``:

    warhead  [*:1]-...        linker  [*:1]-...-[*:2]        e3  [*:2]-...
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass

from rdkit import Chem

from .e3_library import identify_e3

MIN_WARHEAD_HEAVY = 6


@dataclass
class Decomposition:
    warhead: str
    linker: str
    e3: str
    e3_family: str
    e3_anchor: str
    linker_heavy: int

    def as_dict(self) -> dict:
        return {
            "warhead": self.warhead,
            "linker": self.linker,
            "e3": self.e3,
            "e3_family": self.e3_family,
            "e3_anchor": self.e3_anchor,
            "linker_heavy": self.linker_heavy,
        }


def _heavy(mol) -> int:
    return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)


def _exit_bonds(mol, atoms: set[int]) -> list[tuple[int, int, int]]:
    """Bonds leaving ``atoms``, as (bond_idx, inside_atom, outside_atom)."""
    out = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in atoms and j not in atoms:
            out.append((b.GetIdx(), i, j))
        elif j in atoms and i not in atoms:
            out.append((b.GetIdx(), j, i))
    return out


MAX_SPACER_SUBSTITUENTS = 2


def _scaffold_atoms(mol) -> set[int]:
    """Atoms belonging to a ring system that reads as a ligand scaffold.

    A scaffold is a ring system that is aromatic, *fused* (two rings sharing a
    bond), or *decorated*.  Everything else is a spacer ring — piperazine,
    azetidine, a cyclohexyl — which linkers are made of and the walk passes
    straight through.

    Three distinctions earn their keep, and each was put here by a molecule
    that broke without it:

    * **Aromatic or fused.** JQ1 attaches through a saturated diazepine fused
      to a triazole and a thiophene, so an aromatic-only test walks into the
      middle of the warhead.
    * **Fused means sharing a bond, not an atom.** Spiro linkers share exactly
      one atom and are still spacers.
    * **Decorated means more than two exocyclic bonds.** A piperazine carrying
      a chain at each nitrogen is a spacer: you go in one side and out the
      other.  The same ring carrying four substituents is a core.  Nutlin-class
      MDM2 warheads hang on a non-aromatic imidazoline with four exocyclic
      bonds; without this test the walk crosses it and emerges in three
      directions at once.
    """
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    systems: list[set[int]] = []
    for r in rings:
        merged = set(r)
        # Fused = sharing a bond, i.e. at least two atoms.  One shared atom is
        # a spiro junction and keeps the rings separate.
        for s in [s for s in systems if len(s & r) >= 2]:
            merged |= s
            systems.remove(s)
        systems.append(merged)

    out: set[int] = set()
    for s in systems:
        if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in s):
            out |= s
            continue
        if sum(1 for r in rings if r <= s) >= 2:  # fused
            out |= s
            continue
        exocyclic = sum(
            1
            for i in s
            for nb in mol.GetAtomWithIdx(i).GetNeighbors()
            if nb.GetIdx() not in s and nb.GetAtomicNum() > 1
        )
        if exocyclic > MAX_SPACER_SUBSTITUENTS:
            out |= s
    return out


def _walk_linker(mol, start: int, blocked: set[int], scaffold: set[int]) -> set[int]:
    """Atoms reachable from ``start`` without entering a ligand scaffold."""
    seen: set[int] = set()
    q = deque([start])
    while q:
        idx = q.popleft()
        if idx in seen or idx in blocked or idx in scaffold:
            continue
        seen.add(idx)
        for nb in mol.GetAtomWithIdx(idx).GetNeighbors():
            n = nb.GetIdx()
            if n not in seen and n not in blocked and n not in scaffold:
                q.append(n)
    return seen


def _frag_with_dummies(mol, bond_idxs: list[int], labels) -> list[Chem.Mol]:
    frag = Chem.FragmentOnBonds(mol, bond_idxs, addDummies=True, dummyLabels=labels)
    return list(Chem.GetMolFrags(frag, asMols=True, sanitizeFrags=True))


def _relabel(frag) -> str:
    rw = Chem.RWMol(frag)
    for a in rw.GetAtoms():
        if a.GetAtomicNum() == 0:
            a.SetAtomMapNum(a.GetIsotope())
            a.SetIsotope(0)
    return Chem.MolToSmiles(rw.GetMol())


def decompose(smiles: str) -> Decomposition | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fam = identify_e3(mol)
    if fam is None:
        return None

    matches = fam.ligand_atoms(mol)
    if not matches:
        return None
    scaffold = _scaffold_atoms(mol)

    best: tuple[int, tuple] | None = None
    for e3_atoms in matches:
        exits = _exit_bonds(mol, e3_atoms)
        # A ligand with no exit vector is the free ligand, not a PROTAC arm.
        for e3_bond, _inside, outside in exits:
            linker_atoms = _walk_linker(mol, outside, blocked=e3_atoms, scaffold=scaffold)
            if not linker_atoms:
                continue
            # Where the linker meets the warhead scaffold.  Only bonds leaving
            # the *linker* count.  Substituents hanging off the recognised E3
            # ligand — the acyl cap on a VHL tert-leucine, the second amide of
            # a hydroxyproline — leave the ligand, not the linker, and counting
            # them here made every capped VHL degrader look branched.
            wh_bonds = [
                b for b, _i, o in _exit_bonds(mol, linker_atoms)
                if o not in e3_atoms
            ]
            if len(wh_bonds) != 1:
                continue  # branched or no clean warhead boundary
            wh_bond = wh_bonds[0]
            if wh_bond == e3_bond:
                continue

            try:
                pieces = _frag_with_dummies(
                    mol, [wh_bond, e3_bond], [(1, 1), (2, 2)]
                )
            except Exception:
                continue
            if len(pieces) != 3:
                continue

            linkers = [
                p for p in pieces
                if sum(a.GetAtomicNum() == 0 for a in p.GetAtoms()) == 2
            ]
            if len(linkers) != 1:
                continue
            linker = linkers[0]
            ends = [p for p in pieces if p is not linker]
            carries = [e for e in ends if e.HasSubstructMatch(fam.ligand())]
            if len(carries) != 1:
                continue
            e3_frag = carries[0]
            warhead = ends[0] if ends[1] is e3_frag else ends[1]
            if _heavy(warhead) < MIN_WARHEAD_HEAVY:
                continue

            # Prefer the split that keeps the most atoms in the warhead:
            # a bigger recognised ligand means a tighter, better-anchored cut.
            key = _heavy(warhead)
            if best is None or key > best[0]:
                best = (key, (warhead, linker, e3_frag))

    if best is None:
        return None
    _, (warhead, linker, e3_frag) = best

    # Convention: warhead carries :1, E3 carries :2.
    wh_label = next(
        a.GetIsotope() for a in warhead.GetAtoms() if a.GetAtomicNum() == 0
    )
    if wh_label != 1:
        for frag in (warhead, linker, e3_frag):
            for a in frag.GetAtoms():
                if a.GetAtomicNum() == 0:
                    a.SetIsotope(1 if a.GetIsotope() == 2 else 2)

    return Decomposition(
        warhead=_relabel(warhead),
        linker=_relabel(linker),
        e3=_relabel(e3_frag),
        e3_family=fam.name,
        e3_anchor=fam.anchor,
        linker_heavy=_heavy(linker),
    )
