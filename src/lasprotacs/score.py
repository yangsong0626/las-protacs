"""Property calculation and triage for assembled PROTACs.

PROTACs live in "beyond rule of 5" space, so Lipinski thresholds are the wrong
yardstick — applying them would reject every clinical degrader, including the
ones that reached patients.  The bands here come from the published property
distributions of degraders that showed cellular activity, and they are used to
*rank*, not to hard-filter.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

# Soft bands for degraders in cells.  Outside these a molecule is not wrong,
# just further from where most active degraders sit.
BANDS = {
    "mw": (700.0, 1100.0),
    "clogp": (2.0, 6.5),
    "tpsa": (90.0, 180.0),
    "hbd": (0, 4),
    "rotb": (6, 20),
}


@dataclass
class Properties:
    mw: float
    clogp: float
    tpsa: float
    hbd: int
    hba: int
    rotb: int
    rings: int
    heavy: int

    def as_dict(self) -> dict:
        return asdict(self)


def properties(mol_or_smiles) -> Properties | None:
    mol = (
        Chem.MolFromSmiles(mol_or_smiles)
        if isinstance(mol_or_smiles, str)
        else mol_or_smiles
    )
    if mol is None:
        return None
    return Properties(
        mw=round(Descriptors.MolWt(mol), 1),
        clogp=round(Crippen.MolLogP(mol), 2),
        tpsa=round(rdMolDescriptors.CalcTPSA(mol), 1),
        hbd=rdMolDescriptors.CalcNumHBD(mol),
        hba=rdMolDescriptors.CalcNumHBA(mol),
        rotb=rdMolDescriptors.CalcNumRotatableBonds(mol),
        rings=rdMolDescriptors.CalcNumRings(mol),
        heavy=mol.GetNumHeavyAtoms(),
    )


def _band_score(value: float, lo: float, hi: float, slack: float) -> float:
    """1.0 inside the band, tapering to 0 one ``slack`` outside it."""
    if lo <= value <= hi:
        return 1.0
    d = lo - value if value < lo else value - hi
    return max(0.0, 1.0 - d / slack)


def developability(props: Properties) -> float:
    """0-1 summary of how close a molecule sits to active-degrader space."""
    parts = [
        _band_score(props.mw, *BANDS["mw"], slack=250.0),
        _band_score(props.clogp, *BANDS["clogp"], slack=2.5),
        _band_score(props.tpsa, *BANDS["tpsa"], slack=60.0),
        _band_score(props.hbd, *BANDS["hbd"], slack=3.0),
        _band_score(props.rotb, *BANDS["rotb"], slack=8.0),
    ]
    return round(sum(parts) / len(parts), 3)


def flags(props: Properties) -> list[str]:
    """Human-readable reasons a molecule sits outside the usual bands."""
    out = []
    if props.mw > BANDS["mw"][1]:
        out.append(f"MW {props.mw:.0f} above the usual degrader range")
    if props.mw < BANDS["mw"][0]:
        out.append(f"MW {props.mw:.0f} small for a degrader")
    if props.clogp > BANDS["clogp"][1]:
        out.append(f"cLogP {props.clogp:.1f} high — solubility risk")
    if props.clogp < BANDS["clogp"][0]:
        out.append(f"cLogP {props.clogp:.1f} low — permeability risk")
    if props.tpsa > BANDS["tpsa"][1]:
        out.append(f"TPSA {props.tpsa:.0f} high — permeability risk")
    if props.rotb > BANDS["rotb"][1]:
        out.append(f"{props.rotb} rotatable bonds — entropic cost on ternary formation")
    if props.hbd > BANDS["hbd"][1]:
        out.append(f"{props.hbd} H-bond donors — permeability risk")
    return out
