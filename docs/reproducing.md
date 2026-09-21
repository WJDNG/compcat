# Where every reported number comes from

Each row names a result in the paper and the script that produces it. The public
package intentionally omits the internal `figures/` and `manuscript/`
directories; rows that name those paths are provenance references to the
non-public submission bundle and cannot be run from this repository alone.

## Main text

| Result | Script | Output |
|---|---|---|
| Table 1 — gas-phase bond-dissociation energies, MAE 0.25 eV | `01_gasphase_references.py` | `data/derived/gasphase_energies.json` |
| Table 2 — CaCeCoLaNi proxy, regenerated entropy | `06_configurational_entropy.py` | `data/derived/entropy_regenerated.csv` |
| Table 3 — proxy across seven clean-slab systems | archived directly | `data/raw/clean_slab_alignn.csv.gz` |
| Table 4 — archived CH4 statistics by system | `05_elemental_correlations.py` | console summary |
| Table 5 — Ni(111) reference, r = 0.95, MAE 0.77 eV on CHx | `02_ni111_benchmark.py` | `data/derived/ni111_reference.json` |
| Tables 6, 7 — site geometries and layer composition | archived directly | `data/raw/all_energy_results.csv.gz` |
| Table 8 — leading CoNiPdRhRu compositions | `09_candidate_descriptors.py` | `data/derived/descriptors_conipdrhru.csv` |
| Figure 3 / Section 2.6 — R² = 0.50–0.65, upper 95% bound 0.674, site spread 0.30 eV | `04_site_representativeness.py` | `data/derived/site_statistics.csv` |
| Figure 4 / Section 2.7 — neighbour census, Ni and Co enrichment 0.98 | `07_neighbour_normalisation.py` | `data/derived/neighbour_normalisation.csv` |
| Figure 5 — descriptor maps | `figures/fig05_desc.py` | `figures/fig05_desc.png` |
| Figure 6 / Section 3.5 — pooled max |r| = 0.73, within-system median 0.10 | `05_elemental_correlations.py` | `data/derived/elemental_correlations.csv` |
| Figure 7 / Section 2.9 — reference benchmark, affine relation | `02_ni111_benchmark.py`, `figures/fig_ref.py` | `figures/fig_reference.png` |
| Figure 10 — DMR energetics | `figures/fig_energy.py` | `figures/fig11_energy.png` |
| Figure 11 / Section 3.5 — candidates, Pd +3.6 pp, Co −1.9 pp, split half r = 0.79 | `09_candidate_descriptors.py`, `figures/fig_cand.py` | `figures/fig_cand.png` |
| Section 3.5 — composition-record test | `08_composition_record_test.py` | `data/derived/composition_record_test.csv` |
| Section 3.7 — bcc facet control, 0.83 Å vs 2.03 Å, CN 7.0 vs 10.0 | `03_bcc_facet_control.py` | `data/derived/bcc_facet_control.json` |
| Section 3.7 — retention 93.0 / 92.9 / 61.0 / 29.5 / 9.1 % | `05_elemental_correlations.py` | console summary |
| Section 2.3 — ALIGNN test MAE 0.033, RMSE 0.070, n = 5,572 | archived directly | `data/raw/alignn_test_predictions.csv.gz` |

## Supporting Information

| Result | Script |
|---|---|
| Table S1, Section S9 — gas-phase energies and coordinates | `01_gasphase_references.py` |
| Table S2 — archived-value distribution statistics | `data/raw/adsorbate_summary.csv` |
| Tables S3, S4 — pooled and within-system correlations | `05_elemental_correlations.py` |
| Table S5 — site statistics with bootstrap intervals | `04_site_representativeness.py` |
| Table S6 — slab geometry and retention, with the (110) control | `03_bcc_facet_control.py` |
| Tables S7, S8 — within-system spread and regenerated entropy | `06_configurational_entropy.py` |
| Figure S3 — distributional distances | `figures/fig12_dist.py` |

## Numbers that cannot be reproduced from this archive

Stated here so that nobody looks for them.

- **The reference state of the archived values.** Never recorded by the
  pipeline. Script 02 determines it externally; that is a reconstruction, not a
  recovery.
- **Which composition record the generator used.** Script 08 identifies one on
  the evidence; the archive does not state it.
- **The Graphormer pre-relaxation stage.** No checkpoint, configuration, code or
  log survives.
- **The cause of individual non-physical outliers.** No convergence flags, no
  trajectories.
- **The (7,6) grid position.** Never computed; it is absent from the trajectory
  filenames themselves, so the omission occurred at job submission.
- **The per-system CoNiPdRhRu file.** Lost after the combined table was
  assembled; the 272,000 rows here were recovered from the combined table, which
  is why their trajectory paths carry no system prefix.
