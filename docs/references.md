# References

## Data

- **PROTAC-8K** — 3D structures of PROTAC ternary-complex components (PROTAC,
  target ligand, E3 ligase ligand and both pockets), released with the
  DeepPROTACs degradation-prediction work.
  <https://doi.org/10.5281/zenodo.7278277>. This is the corpus decomposed here:
  6,113 mol2 records, 6,107 unique structures. It ships structures and
  degradation labels, but no SMILES-level activity table.
- Li S. et al. **DeepPROTACs is a deep learning-based targeted degradation
  predictor for PROTACs.** *Nature Communications* 13, 7133 (2022).
  <https://doi.org/10.1038/s41467-022-34807-3> — the model PROTAC-8K was built
  for, and the source of its train/test construction.
- Weng G. et al. **PROTAC-DB 3.0: an updated database of PROTACs with extended
  pharmacokinetic parameters.** *Nucleic Acids Research* 53, D1510 (2025).
  <https://doi.org/10.1093/nar/gkae768> — 6,111 PROTACs, 569 warheads,
  2,753 linkers, 107 E3 ligands, DC50 for ~1,300 compounds. Ships its own
  warhead / linker / E3 split, which makes it an independent check on the
  decomposition here as well as the route to a trained model.

## Method background

- Landrum G. **RDKit: open-source cheminformatics.** <https://www.rdkit.org/>
  — decomposition, `molzip` assembly, descriptors, CIP labelling.
- O'Boyle N.M. et al. **Open Babel: An open chemical toolbox.** *Journal of
  Cheminformatics* 3, 33 (2011). <https://doi.org/10.1186/1758-2946-3-33> —
  used to read the Tripos mol2 files, which RDKit's reader rejects on
  kekulisation of fused heteroaromatics.
- Ribes S., Nittinger E., Tyrchan C., Mercado R. **Modeling PROTAC degradation
  activity with machine learning.** *Artificial Intelligence in the Life
  Sciences* 6, 100104 (2024). <https://doi.org/10.1016/j.ailsci.2024.100104>
  — curation choices for pooled PROTAC activity data.
- **A Comprehensive Review of Emerging Approaches in Machine Learning for
  De Novo PROTAC Design.** arXiv:2406.16681 — survey of generative approaches;
  most published work designs the linker with the two ligands fixed, rather
  than varying all three positions.

## On the limits of pooled degradation data

Reported DC50 and Dmax are not a single comparable quantity across
publications: they depend on cell line, timepoint, detection method and
normalisation. Any model trained on a pooled table inherits that variance,
which bounds how well it can generalise to a new target or ligase. This is why
this repo ranks on properties and chemistry rules and does not yet claim a
degradation prediction.
