"""Recognition patterns for the E3 ligase ligands used in PROTACs.

Two patterns per family:

``core_smarts``    small, highly specific — answers "which ligase is this?"
``ligand_smarts`` the whole ligand as a medicinal chemist would draw it —
                  this is what anchors the cut, so the linker starts where the
                  recognised ligand stops rather than wherever a search happens
                  to land.

Aromatic/aliphatic atoms in a SMARTS match regardless of what else hangs off
them, so one pattern covers a whole substituted series.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem


@dataclass
class E3Family:
    name: str  # ligase recruited
    anchor: str  # chemotype label
    core_smarts: str
    ligand_smarts: str
    note: str = ""
    _core: object = field(default=None, repr=False)
    _ligand: object = field(default=None, repr=False)

    def core(self):
        if self._core is None:
            self._core = Chem.MolFromSmarts(self.core_smarts)
        return self._core

    def ligand(self):
        if self._ligand is None:
            self._ligand = Chem.MolFromSmarts(self.ligand_smarts)
        return self._ligand

    def matches(self, mol) -> bool:
        c = self.core()
        return c is not None and mol.HasSubstructMatch(c)

    def ligand_atoms(self, mol) -> list[set[int]]:
        """Every distinct whole-ligand match, largest first."""
        lig = self.ligand()
        if lig is None:
            return []
        return [set(m) for m in mol.GetSubstructMatches(lig)]


E3_FAMILIES: list[E3Family] = [
    E3Family(
        name="CRBN",
        anchor="phthalimide-glutarimide",
        # Glutarimide = piperidine-2,6-dione; the ring NH donates to CRBN.
        core_smarts="O=C1CCCC(=O)[NX3H1]1",
        # Thalidomide / pomalidomide class: glutarimide + phthalimide.
        ligand_smarts="O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)[NX3H1]1",
        note="IMiD class, phthalimide variant",
    ),
    E3Family(
        name="CRBN",
        anchor="isoindolinone-glutarimide",
        core_smarts="O=C1CCCC(=O)[NX3H1]1",
        # Lenalidomide class: glutarimide + isoindolinone (one carbonyl).
        ligand_smarts="O=C1CCC(N2Cc3ccccc3C2=O)C(=O)[NX3H1]1",
        note="IMiD class, isoindolinone variant",
    ),
    E3Family(
        name="VHL",
        anchor="VH032",
        core_smarts="[OX2H1][CH1]1C[CH1](C(=O)[NX3])N([CX3]=O)C1",
        # Full VH032: methylthiazole-benzyl amide + hydroxyproline +
        # tert-leucine.  The tert-leucine nitrogen is the usual exit vector.
        ligand_smarts="Cc1ncsc1-c1ccc(CNC(=O)C2CC(O)CN2C(=O)C([NX3])C(C)(C)C)cc1",
        note="VH032 / VHL-ligand class",
    ),
    E3Family(
        name="VHL",
        anchor="hydroxyproline-generic",
        core_smarts="[OX2H1][CH1]1C[CH1](C(=O)[NX3])N([CX3]=O)C1",
        # Fallback for VHL ligands with a non-standard cap on the benzyl side.
        ligand_smarts="[OX2H1]C1CC(C(=O)[NX3])N(C(=O)C([NX3])C(C)(C)C)C1",
        note="VHL ligand without the methylthiazole recognised",
    ),
    E3Family(
        name="IAP",
        anchor="AVPI-mimetic",
        core_smarts="[CH3][NX3H1][CH1]([CH3])C(=O)[NX3H1][CH1]([CH1]([CH3])[CH3])C(=O)",
        ligand_smarts="[CH3][NX3H1][CH1]([CH3])C(=O)N1[CH1](CCC1)C(=O)[NX3]",
        note="SMAC-mimetic class",
    ),
    E3Family(
        name="MDM2",
        anchor="nutlin",
        core_smarts="Clc1ccc(cc1)C1=N[CH1]([CH1](N1)c1ccc(Cl)cc1)",
        ligand_smarts="Clc1ccc(cc1)C1=N[CH1]([CH1](N1)c1ccc(Cl)cc1)",
        note="nutlin / idasanutlin class",
    ),
]


def identify_e3(mol) -> E3Family | None:
    """First family whose core is present AND whose whole ligand matches."""
    fallback = None
    for fam in E3_FAMILIES:
        if not fam.matches(mol):
            continue
        if fam.ligand_atoms(mol):
            return fam
        if fallback is None:
            fallback = fam
    return fallback
