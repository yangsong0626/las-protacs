"""Decompose -> reassemble must reproduce the original molecule exactly."""

import sys
from pathlib import Path

import pytest
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lasprotacs.assemble import assemble, same_molecule
from lasprotacs.decompose import _scaffold_atoms, decompose
from lasprotacs.rules import (
    attachment_ok,
    attachment_reason,
    bond_ok,
    new_alerts,
    repair_attachment,
)

JQ1 = "CC1=C(C)C2=C(S1)N3C(C)=NN=C3[C@@H](CC(=O)N{R})N=C2C4=CC=C(Cl)C=C4"

PROTACS = {
    "vhl_hdac": (
        "Cc1ncsc1-c1ccc(CNC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](NC(=O)COCCCCCCC"
        "CCOCC(=O)Nc2ccc(C(=O)Nc3ccc(F)cc3N)cc2)C(C)(C)C)cc1"
    ),
    "dbet_like": JQ1.format(R="CCCCCC(=O)Nc5cccc6c5C(=O)N(C7CCC(=O)NC7=O)C6=O"),
    "mz1_like": JQ1.format(
        R="CCOCCOCCC(=O)N[C@H](C(=O)N5C[C@H](O)C[C@H]5C(=O)NCc6ccc(-c7scnc7C)cc6)C(C)(C)C"
    ),
    "piperazine_crbn": (
        "O=C1CCC(N2C(=O)c3cccc(N4CCN(CCC(=O)Nc5ccc(-c6ccncc6)cc5)CC4)c3C2=O)C(=O)N1"
    ),
}

EXPECTED_FAMILY = {
    "vhl_hdac": "VHL",
    "dbet_like": "CRBN",
    "mz1_like": "VHL",
    "piperazine_crbn": "CRBN",
}


@pytest.mark.parametrize("name", sorted(PROTACS))
def test_roundtrip(name):
    smiles = PROTACS[name]
    d = decompose(smiles)
    assert d is not None, f"{name} did not decompose"
    assert d.e3_family == EXPECTED_FAMILY[name]
    rebuilt = assemble(d.warhead, d.linker, d.e3)
    assert rebuilt == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def test_linker_keeps_tert_leucine_with_the_vhl_ligand():
    """The classic failure mode: a free search shaves tBu off the VHL ligand."""
    d = decompose(PROTACS["mz1_like"])
    assert "C(C)(C)C" in d.e3, "tert-leucine must stay on the E3 ligand"
    assert "C(C)(C)C" not in d.linker


def test_piperazine_stays_in_the_linker():
    d = decompose(PROTACS["piperazine_crbn"])
    assert "N1CCN" in d.linker.replace("[*:1]", "").replace("[*:2]", "")


def test_heteroatom_bonds_are_refused():
    assert bond_ok(6, 7) and bond_ok(6, 6)      # C-N, C-C fine
    assert not bond_ok(7, 7)                     # hydrazine
    assert not bond_ok(8, 8)                     # peroxide


def test_attachment_rule_blocks_amine_to_amine():
    vhl = "Cc1ncsc1-c1ccc(CNC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](N[*:2])C(C)(C)C)cc1"
    piperazine = "O=C(CCN1CCN([*:2])CC1)N[*:1]"
    warhead = "Nc1cc(F)ccc1NC(=O)c1ccc([*:1])cc1"
    assert not attachment_ok(warhead, piperazine, vhl)


def test_alerts_ignore_motifs_inherited_from_pieces():
    warhead = "CS(=O)(=O)Nc1ccccc1[*:1]"       # sulfonamide already present
    linker = "O=C(C[*:1])NCC[*:2]"
    e3 = "O=C1CCC(N2C(=O)c3cccc([*:2])c3C2=O)C(=O)N1"
    product = assemble(warhead, linker, e3)
    assert product is not None
    assert new_alerts(product, [warhead, linker, e3]) == []


def test_rejection_names_the_bond_and_what_it_makes():
    vhl = "Cc1ncsc1-c1ccc(CNC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](N[*:2])C(C)(C)C)cc1"
    piperazine = "O=C(CCN1CCN([*:2])CC1)N[*:1]"
    warhead = "Nc1cc(F)ccc1NC(=O)c1ccc([*:1])cc1"
    why = attachment_reason(warhead, piperazine, vhl)
    assert "N–N" in why and "hydrazine" in why
    assert attachment_reason(warhead, "O=C(C[*:1])NCC[*:2]", vhl) is None


# ---------------------------------------------------------------------------
# The three bugs below were invisible on a hand-picked set of four degraders
# and only showed up running the full PROTAC-8K corpus. Each test pins one.
# ---------------------------------------------------------------------------

CAPPED_VHL = (
    "COc1cc(-c2cn(C)c(=O)c3cnccc23)cc(OC)c1CN1CCN(CCOCCOCCOc2cc(-c3scnc3C)ccc2CNC"
    "(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](NC(=O)C2(F)CC2)C(C)(C)C)CC1"
)

NUTLIN_CRBN = (
    "CCOc1cc(C(C)(C)C)ccc1C1=N[C@@H](c2ccc(Cl)cc2)[C@H](c2ccc(Cl)cc2)N1C(=O)N1CCN"
    "(C(=O)CCC#Cc2cccc3c2CN([C@@H]2CCC(=O)NC2=O)C3=O)CC1"
)

RING_STEREO = (
    "N#CC1(CNc2cccc(-c3cc(N[C@H]4CC[C@H](NCCCCCCNc5cccc6c5C(=O)N([C@H]5CCC(=O)NC5"
    "=O)C6=O)CC4)ncc3Cl)n2)CCOCC1"
)


def test_a_cap_on_the_e3_ligand_is_not_a_warhead_boundary():
    """VHL ligands are often acylated on the tert-leucine nitrogen.

    That cap leaves the *ligand*, not the linker. Counting it as a second
    warhead boundary made every capped VHL degrader look branched and dropped
    758 structures on the floor.
    """
    d = decompose(CAPPED_VHL)
    assert d is not None, "capped VHL degrader must still decompose"
    assert d.e3_family == "VHL"
    assert "C2(F)CC2" in d.e3 or "C(F)" in d.e3, "the acyl cap belongs to the ligand"
    assert same_molecule(assemble(d.warhead, d.linker, d.e3), CAPPED_VHL)


def test_a_decorated_non_aromatic_ring_is_a_scaffold_not_a_spacer():
    """Nutlin-class warheads hang on a non-aromatic imidazoline.

    It carries four exocyclic bonds, so it is a core, not something the walk
    should pass through on its way out of the linker.
    """
    d = decompose(NUTLIN_CRBN)
    assert d is not None
    assert d.e3_family == "CRBN"
    # Both chlorophenyls belong to the warhead, not the linker.
    assert d.warhead.count("Cl") == 2
    assert "Cl" not in d.linker
    assert same_molecule(assemble(d.warhead, d.linker, d.e3), NUTLIN_CRBN)


def test_spiro_rings_are_spacers_but_fused_rings_are_not():
    spiro = Chem.MolFromSmiles("C1CCC2(CC1)CCNCC2")       # shares one atom
    fused = Chem.MolFromSmiles("C1CCC2CCCCC2C1")          # shares a bond
    assert _scaffold_atoms(spiro) == set()
    assert _scaffold_atoms(fused) != set()


def test_ring_stereo_relabelling_is_not_a_round_trip_failure():
    """1,4-disubstituted cyclohexanes rewrite both @ tags when reassembled.

    The compound is unchanged — every CIP label matches — so comparing
    canonical SMILES as strings reports a failure that is not one.
    """
    d = decompose(RING_STEREO)
    assert d is not None
    rebuilt = assemble(d.warhead, d.linker, d.e3)
    assert rebuilt != Chem.MolToSmiles(Chem.MolFromSmiles(RING_STEREO))
    assert same_molecule(rebuilt, RING_STEREO)


def test_same_molecule_still_separates_genuine_stereoisomers():
    assert not same_molecule("C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O")
    assert same_molecule("C[C@H](N)C(=O)O", "C[C@H](N)C(=O)O")
    assert not same_molecule("CCO", "CCN")


# ---------------------------------------------------------------------------
# Repair: a refused join is usually one atom from a legal one.
# ---------------------------------------------------------------------------

VHL_LIGAND = (
    "Cc1ncsc1-c1ccc(CNC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](N[*:2])C(C)(C)C)cc1"
)
PIPERAZINE = "O=C(CCN1CCN([*:2])CC1)N[*:1]"
ANILIDE = "Nc1cc(F)ccc1NC(=O)c1ccc([*:1])cc1"


def test_repair_turns_a_hydrazine_join_into_an_amine():
    assert not attachment_ok(ANILIDE, PIPERAZINE, VHL_LIGAND)
    fix = repair_attachment(ANILIDE, PIPERAZINE, VHL_LIGAND)
    assert fix is not None
    # Piperazine becomes a piperidine: the ring nitrogen at the exit vector is
    # now carbon, so the bond formed is C-N rather than N-N.
    assert attachment_ok(ANILIDE, fix["linker"], VHL_LIGAND)
    assert attachment_reason(ANILIDE, fix["linker"], VHL_LIGAND) is None
    product = assemble(ANILIDE, fix["linker"], VHL_LIGAND)
    assert product is not None and Chem.MolFromSmiles(product) is not None


def test_repair_edits_the_linker_and_never_the_ligands():
    """The ligands' exit vectors are pharmacophore; only the linker is fair game."""
    fix = repair_attachment(ANILIDE, PIPERAZINE, VHL_LIGAND)
    product = assemble(ANILIDE, fix["linker"], VHL_LIGAND)
    mol = Chem.MolFromSmiles(product)
    for ligand in (ANILIDE, VHL_LIGAND):
        q = Chem.MolFromSmiles(ligand.replace("[*:1]", "[H]").replace("[*:2]", "[H]"))
        q = Chem.RemoveHs(q)
        assert mol.HasSubstructMatch(q), "a published ligand was altered by the repair"


def test_repair_declines_when_nothing_is_wrong():
    good = "O=C(C[*:1])NCC[*:2]"
    assert attachment_ok(ANILIDE, good, VHL_LIGAND)
    assert repair_attachment(ANILIDE, good, VHL_LIGAND) is None


def test_repair_is_refused_when_it_would_break_the_ring():
    """Mutating the nitrogen out of a pyrrole leaves a ring that cannot exist.

    The saturated control matters: it shows the None above is the guard firing,
    not the input failing to parse.
    """
    from lasprotacs.rules import _mutate_attachment_to_carbon as mutate

    assert Chem.MolFromSmiles("c1ccn([*:2])c1") is not None
    assert mutate("c1ccn([*:2])c1", 2) is None
    assert mutate("C1CCN([*:2])CC1", 2) == ("C1CCC([*:2])CC1", "N→C")
