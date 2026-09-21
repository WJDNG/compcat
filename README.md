# compcat

Reproduction code and data for **"Reference Diagnostics and Sampling Design for
Machine-Learning Screening of High-Entropy Alloy Catalysts: 646,000
Configurations for Dry Methane Reforming."**

Most numerical results reported in the paper are produced by a script in
`scripts/` from data in `data/`. Each script names, in its docstring, the table
or numbered claim it reproduces. `docs/reproducing.md` gives the full mapping
and records which figure-generation sources belong to the non-public submission
bundle.

## What is here

The study screens twelve five-element high-entropy alloy systems for dry methane
reforming with a machine-learning surrogate, and sets out the diagnostics that
establish what such a screen can support. The repository holds the archived
screening output, the independent reference calculations performed against it,
and the analysis that connects the two.

## Layout

```
data/raw/        the archived screening output, gzipped, with checksums
data/derived/    everything the scripts compute from it
src/compcat/     a small package: archive I/O, the window, descriptors, statistics
scripts/         numbered reproduction scripts, run in any order
docs/            the mapping from every reported number to the script that makes it
```

The public package intentionally omits the internal `figures/` and `manuscript/`
directories. This keeps the repository limited to the shared analysis code,
archived data and reproduction notes while preserving the provenance map in
`docs/reproducing.md`.

## Quick start

```bash
git clone <repository-url> && cd compcat
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

python scripts/07_neighbour_normalisation.py     # seconds, no heavy dependencies
python scripts/06_configurational_entropy.py     # seconds
python scripts/04_site_representativeness.py     # ~1 min
python scripts/05_elemental_correlations.py      # ~2 min
python scripts/08_composition_record_test.py     # ~2 min
python scripts/09_candidate_descriptors.py       # ~3 min
```

Scripts 01–03 recompute the reference calculations rather than the analysis and
need `pyscf` (01) or `ase` and `mace-torch` (02, 03). Their outputs are included
in `data/derived/`, so the analysis scripts run without them.

## Data provenance

`data/raw/` is the archived output of the original screening run, unmodified
apart from compression. `SHA256SUMS.txt` covers every file. The ALIGNN
checkpoint used for the formation-energy surrogate is the upstream JARVIS-DFT
model and is not redistributed here; `alignn_config.json` identifies it
(dataset `dft_3d`, target `formation_energy_peratom`, 300 epochs, split
80/10/10, seed 123, code version `9835fe0d4b313e2522034ff39f0ebdbfecde99a2`),
and `alignn_test_predictions.csv.gz` carries its held-out predictions
(MAE 0.033 eV/atom, RMSE 0.070 eV/atom, n = 5,572).

The Graphormer pre-relaxation stage described in the workflow is not represented
in the archive by any checkpoint, configuration or code, and cannot be
reproduced from what survives. It supplies starting geometries only; every
analysed number is read from the energy tables produced after the subsequent
BFGS optimization.

## Licence

Code under MIT (`LICENSE`). Data under CC BY 4.0 (`data/LICENSE`). If you use
either, please cite the paper (`CITATION.cff`).
