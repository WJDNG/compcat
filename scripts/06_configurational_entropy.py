"""Configurational entropy regenerated over the metal sublattice.

The archived entropies were evaluated over every species the parser returned,
including adsorbate atoms and species that cannot be present. The archived mean
exceeds the ideal five-component limit and is therefore impossible.

Produces: data/derived/entropy_regenerated.csv
Reproduces: Table 2, Table S8 and the 62 -> 11 result of Section 3.2.
"""
import csv
import gzip
import math
import re
import statistics as st

import _paths
from compcat.io import ARCHIVE, label_composition

KB = 8.617333262e-5      # eV / (atom K)
T = 1273.0
IDEAL = math.log(5)      # five-component ideal limit, in units of kB

rows = []
with gzip.open(ARCHIVE / "cacecolani_negative_proxy.csv.gz", "rt", newline="") as fh:
    for row in csv.DictReader(fh):
        label = re.search(r"outputs/([A-Za-z0-9]+)-FCC/", row["traj_file"]).group(1)
        x = list(label_composition(label).values())
        s = -KB * sum(xi * math.log(xi) for xi in x if xi > 0)
        ef = float(row["formation_energy_eV_per_atom"])
        sa = float(row["S_config_eV_per_atom_K"])
        rows.append({"composition_label": label,
                     "E_form": round(ef, 4),
                     "S_archived_kB": round(sa / KB, 3),
                     "S_metal_sublattice_kB": round(s / KB, 3),
                     "S_metal_sublattice_1e4": round(s * 1e4, 4),
                     "TS_1273K": round(T * s, 4),
                     "proxy_archived": round(ef - T * sa, 4),
                     "proxy_regenerated": round(ef - T * s, 4)})

rows.sort(key=lambda r: r["proxy_regenerated"])
with open(_paths.DERIVED / "entropy_regenerated.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

neg_arch = sum(r["proxy_archived"] < 0 for r in rows)
neg_new = sum(r["proxy_regenerated"] < 0 for r in rows)
shift = st.mean(r["proxy_regenerated"] - r["proxy_archived"] for r in rows)
print(f"configurations listed as having a negative proxy: {neg_arch}")
print(f"after regeneration over the metal sublattice:     {neg_new}")
print(f"every proxy rises by {shift:+.3f} eV/atom")
print(f"mean archived entropy    {st.mean(r['S_archived_kB'] for r in rows):.3f} kB"
      f"  (ideal five-component limit {IDEAL:.3f} kB - the archived value exceeds it)")
print(f"mean regenerated entropy {st.mean(r['S_metal_sublattice_kB'] for r in rows):.3f} kB")
