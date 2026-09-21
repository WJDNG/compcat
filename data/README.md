# Data

## raw/

The archived output of the screening run, unmodified apart from gzip
compression. `SHA256SUMS.txt` covers every file.

| File | Rows | What it is |
|---|---|---|
| `all_energy_results.csv.gz` | 646,000 | the screen: one E_arch per configuration, 17 adsorbates, 12 systems |
| `df_ads.csv.gz` | 38,000 | the CH4-resolved subset, eight surface positions for CoNiPdRhRu |
| `df_ads_alignn_results.csv.gz` | 24,000 | post-processed CH4 table carrying the second composition record |
| `clean_slab_alignn.csv.gz` | 7,000 | clean-slab formation-energy screen, 1,000 decorations per system |
| `cacecolani_negative_proxy.csv.gz` | 62 | configurations the archived screen listed as having a negative proxy |
| `alignn_test_predictions.csv.gz` | 5,572 | held-out predictions of the formation-energy surrogate |
| `alignn_config.json` | — | identifies the surrogate checkpoint |
| `adsorbate_summary.csv` | 17 | distribution statistics per adsorbate |
| `wasserstein_matrix.csv`, `jensen_shannon_matrix.csv` | 17 × 17 | distributional distances |

### Columns of the energy tables

`slab_index`, `element_combo`, `adsorbate`, `position`, `energy`, `traj_file`,
`input_file`, `output_file`.

Two fields need care. `energy` is the archived quantity written `E_arch` in the
paper; it is not an adsorption energy and carries no recorded reference state.
`element_combo` is stored character-separated, so `Al22Co21Cu19Fe23Ni23` appears
as `A-l-2-2-C-o-2-1-...`; a parser reading element symbols from it recovers
single letters — C from Co or Cu, F from Fe, H from Hf — which is the origin of
the spurious species in the archived composition and entropy columns. Use the
directory name inside `traj_file` instead, and see `scripts/08` for the test of
what that name is worth.

## derived/

Everything the scripts compute. Safe to delete; rerunning the scripts rebuilds
it. The three files from the reference calculations
(`gasphase_energies.json`, `ni111_reference.json`, `bcc_facet_control.json`) are
kept in the repository because recomputing them takes hours and needs the
optional dependencies.

## Licence

CC BY 4.0. Attribute the accompanying paper.
