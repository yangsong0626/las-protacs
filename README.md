# Las Protacs

**A slot machine for PROTACs.** Pull three reels — a warhead, a linker, an E3
ligase ligand — and watch a degrader assemble. Then watch most of what you
assemble get thrown away.

The house edge is chemistry: of the 165 million molecules these components
formally spell out, a third are bonds no chemist would make.

**[Try it →](https://yangsong0626.github.io/las-protacs/)**

![Las Protacs: spinning a warhead, linker and E3 ligand into a degrader, and a
pairing refused before assembly](docs/demo.gif)

*Three spins. The middle one is refused before anything is built: joining that
linker to that E3 ligand would make a hydroxylamine, so no molecule is drawn
and nothing is scored.*

A PROTAC is three pieces bolted together, which makes it look like a
combinatorics problem. Decompose every degrader in a public structure set and
the pieces recombine into **1.65 × 10⁸** formal molecules. That number is a
trap: a third of them are not molecules a chemist would draw, and the
arithmetic cannot tell you which third.

This repo is the smallest honest version of that pipeline — decompose real
degraders into their parts, recombine the parts, and apply the filters that
decide what is worth a slot in a synthesis queue. Everything below is measured
on a public dataset, not asserted.

## The library is decomposed, not invented

| | |
| --- | --- |
| degraders read (PROTAC-8K) | 6,107 unique structures |
| decomposed into warhead / linker / E3 | **5,181 (84.8%)** |
| reassembled to the original molecule | **5,181 — zero failures** |
| warheads / linkers / E3 ligands recovered | 990 / 2,829 / 59 |
| formal combinations they span | 165,241,890 |

The recovered components rank the way a medicinal chemist would expect:
pomalidomide and VH032 are the two most-used E3 ligands (2,059 and 878
degraders), and the most-used warheads are ibrutinib, an olaparib-type
phthalazinone, palbociclib and JQ1. Nothing in the pipeline knows those names —
they fall out of the decomposition.

## The funnel

Enumeration is the easy part. The interesting part is what each filter removes.

| stage | what it does | survives |
| --- | --- | --- |
| raw cross-product | every warhead × linker × E3 | 165,241,890 |
| attachment chemistry | the bond about to be formed must be one a chemist would make | **105,765,932** (64.0%) |
| structure alerts | motifs *created by the join*, not inherited from the pieces | ~95,500,000 (90.3% ± 0.4) |
| property bands | ranks what is left against cell-active degrader space | ranked |

**The attachment count is exact, and nothing was enumerated to get it.**
Whether a bond may be formed depends only on the two elements being joined, so
the rule factorises: group the components by the element at their attachment
point and the survivor count is a sum over element pairs. 1.65 × 10⁸
combinations, filtered in milliseconds, no molecule built.

The alert and property stages depend on the assembled molecule, so those are
estimated from a uniform random sample of 20,000 with a Wilson interval, and
reported as an interval rather than a point.

The first filter is the one people skip. Joining a VHL ligand (exit vector: an
amine nitrogen) to a piperazine linker (exit vector: a ring nitrogen) gives a
hydrazine. **Thirty-six percent of this space is nonsense before a single
property is computed** — 59 million molecules that a naive enumeration would
have shipped.

## Method

### Decomposition is anchored, not searched

Given a PROTAC, find the warhead, the linker and the E3 ligand.

The obvious approach — try every pair of cut bonds and score the result — fails
in a specific and instructive way. A short linker always looks more
"linker-like" than a long one, so the optimiser shaves the *tert*-leucine off a
VHL ligand and hands the alkyl chain to the E3 arm. Both pieces look plausible
in isolation. Both are wrong.

So the cut is anchored to chemistry instead:

1. Match the **whole E3 ligand** (not just its recognition core) and cut the
   single bond that leaves it. The ligand boundary is a known quantity, not a
   free parameter.
2. Walk outwards from that exit vector through chain atoms and *spacer* rings.
   Stop at the first **scaffold** ring system — that is where the warhead
   begins.

A ring is a scaffold if it is aromatic, fused, or decorated with more than two
exocyclic bonds; otherwise it is a spacer the walk passes through. All three
tests were put there by a molecule that broke without them, and the last two
only surfaced at corpus scale.

### What the corpus found that four hand-picked molecules did not

The decomposition passed a four-degrader round-trip test for a while. Running
it over 6,107 real structures moved it from 71% to 84.8%, by exposing three
bugs that a curated test set is structurally incapable of showing:

- **A cap on the E3 ligand is not a warhead boundary.** VHL ligands are
  routinely acylated on the *tert*-leucine nitrogen. That cap leaves the
  *ligand*; the original code counted bonds leaving the ligand *and* the
  linker, so every capped VHL degrader looked branched and was dropped. **758
  structures.**
- **A decorated non-aromatic ring is a core, not a spacer.** Nutlin-class MDM2
  warheads hang on an imidazoline carrying four exocyclic bonds. Treating any
  non-aromatic ring as a spacer walked straight through it and emerged in three
  directions at once. Hence the exocyclic-bond test — and the matching
  refinement that *fused* means sharing a bond, so spiro linkers stay spacers.
- **219 "round-trip failures" were not failures.** An `@`/`@@` tag is written
  relative to neighbour ordering, so fragmenting and rezipping a
  1,4-disubstituted cyclohexane — a common rigidifying motif — flips both tags
  together and leaves the compound untouched. Comparing CIP labels instead of
  SMILES strings turns 219 spurious mismatches into zero real ones. The
  reported 100% round trip depends on getting this comparison right, so
  `same_molecule()` is tested against genuine stereoisomers too.

Each is pinned by a test named after the molecule that exposed it.

### Assembly

Pieces carry atom-mapped dummies (`[*:1]` toward the warhead, `[*:2]` toward
the E3 ligand) and are fused with `Chem.molzip`. Every bond formed this way is
a bond that appeared in a real degrader at least once.

### Rules

`attachment_ok` runs **before** assembly, on the pair of atoms about to be
bonded: at least one must be carbon. That single rule removes hydrazines,
hydroxylamines, peroxides and sulfenamides in one pass, cheaply.
`attachment_reason` says which bond was refused and what it would have made —
the demo prints it verbatim rather than saying "invalid".

`new_alerts` runs **after** assembly and reports only motifs the product has
that its pieces did not. This distinction matters: a warhead that already
contains a sulfonamide or a cyclic acetal is the medicinal chemist's business.
Flagging it would be the tool second-guessing a published molecule.

### Scoring

PROTACs live beyond the rule of 5, so Lipinski thresholds are the wrong
yardstick — applying them rejects every clinical degrader. The bands here
(MW 700–1100, cLogP 2–6.5, TPSA 90–180, 6–20 rotatable bonds) come from the
property distributions of degraders with reported cellular activity, and they
**rank** rather than gate.

## What this does not do

The limitations are the point of the exercise, so they get their own section.

- **No ternary complex geometry.** This is the big one. Degradation depends on
  whether the warhead and E3 ligand can hold two proteins in a productive
  orientation long enough to transfer ubiquitin. Nothing here models that.
  Two molecules with identical 2D properties can differ completely in whether
  they degrade anything. A 2D filter tells you what is worth *modelling*, not
  what will work.
- **No trained degradation model.** The property bands are literature
  heuristics. Training on measured DC50 is the obvious next step and the data
  exists ([PROTAC-DB 3.0](https://doi.org/10.1093/nar/gkae768) reports DC50 for
  ~1,300 compounds), but cross-paper DC50 is not a single comparable quantity:
  it varies with cell line, timepoint, detection method and normalisation, and
  the resulting inter-laboratory variance puts a ceiling on how well any model
  can generalise. A model trained on the pooled table without accounting for
  that would report an accuracy it does not have. PROTAC-8K, used here, carries
  structures and labels but no SMILES-level activity table.
- **15.2% of the corpus does not decompose**, and the residue is not random:
  - ~390 structures use an E3 ligand the library does not recognise. Coverage
    is CRBN, VHL, IAP and MDM2; the long tail (DCAF15, DCAF16, RNF114, KEAP1)
    is not.
  - ~100 have **no linker at all** — the E3 ligand is bonded straight onto the
    warhead scaffold, biaryl or alkyne. A three-piece model has nowhere to put
    them.
  - ~190 have **branched or peptidic linkers** that reach the warhead at more
    than one point. Linear topology only; macrocyclic and bivalent degraders
    are out too.
  - ~12 have a warhead with no ring system at all (free hydroxamic acids), so
    there is no scaffold for the walk to stop at.
- **Names are not data.** The dataset ships structures only. The warhead and E3
  names in the demo were assigned by hand from structure, and only where the
  chemotype is unambiguous; the rest carry a neutral label, because a wrong
  drug name is worse than no name. Because they are a claim this repo makes on
  its own, `scripts/validate_components.py` checks every one of them against an
  independently drawn parent — a compound name must be a substructure of that
  drug, a chemotype name must match a diagnostic SMARTS. It earned its keep
  immediately by catching a misdrawn piperazine in the palbociclib reference.
- **No synthetic accessibility.** A combination can pass every filter here and
  still take a chemist six weeks.

## Running it

```bash
pip install -e .

# 1. mol2 -> SMILES.  RDKit's mol2 reader fails on most of these files
#    (Tripos `ar` bonds, KekulizeException); Open Babel reads all 6,113.
python scripts/prepare_protac8k.py ~/Downloads/PROTAC-8K/protac

# 2. decompose the corpus and build the component library
python scripts/decompose_corpus.py data/protac8k/protac.smi

# 3. measure the funnel over the whole 1.65e8 space
python scripts/enumerate_space.py --sample 20000

# 4. build the reel the demo spins, and serve it
python scripts/build_demo_data.py
python -m http.server --directory docs    # open localhost:8000
```

```python
from lasprotacs.decompose import decompose
from lasprotacs.assemble import assemble, same_molecule
from lasprotacs.rules import attachment_reason
from lasprotacs.score import properties, developability

d = decompose("Cc1sc2c(c1C)C(c1ccc(Cl)cc1)=N[C@H](CC(=O)NCCCCCC(=O)"
              "Nc1cccc3c1C(=O)N(C1CCC(=O)NC1=O)C3=O)c1nnc(C)n1-2")
print(d.e3_family, d.linker)        # CRBN  O=C(C[*:1])NCCCCCC(=O)N[*:2]

smi = assemble(d.warhead, d.linker, d.e3)
print(developability(properties(smi)))

vhl = "Cc1ncsc1-c1ccc(CNC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](N[*:2])C(C)(C)C)cc1"
print(attachment_reason("Nc1ccc([*:1])cc1", "O=C(CCN1CCN([*:2])CC1)N[*:1]", vhl))
# linker-E3 join is N–N: that makes a hydrazine
```

## Data and provenance

Structures come from the **PROTAC-8K** set
([Zenodo](https://doi.org/10.5281/zenodo.7278277)), read with
[Open Babel](https://openbabel.org/) and canonicalised with
[RDKit](https://www.rdkit.org/). The pipeline also reads
[PROTAC-DB 3.0](http://cadd.zju.edu.cn/protacdb/) exports, which ship their own
warhead / linker / E3 decomposition alongside activity data — useful as an
independent check on the decomposition, and as the route to a trained model.

RDKit.js is vendored under `docs/vendor/` rather than loaded from a CDN, so the
demo still draws molecules on a network that blocks one.

This is an independent implementation built on public data and published
methods. See [`docs/references.md`](docs/references.md).

## License

MIT.
