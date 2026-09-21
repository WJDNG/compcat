"""Elemental correlations, pooled and within alloy family.

Produces: data/derived/elemental_correlations.csv
Reproduces: Figure 6, Tables S3 and S4, and the Simpson's paradox result of
Section 3.5 (largest pooled |r| = 0.73; within-system median 0.10, 86% below 0.2).
"""
import csv
import statistics as st
from collections import defaultdict


def _mean(seq):
    """statistics.mean uses exact arithmetic and is far too slow here."""
    seq = list(seq)
    return sum(seq) / len(seq)

import _paths
from compcat import load_configurations, label_composition, pearson, SYSTEMS
from compcat.window import WINDOW

ADSORBATES = ["CH4", "CH3", "CH2", "CH", "C", "H", "O", "OH", "CO", "CO2"]
MIN_N = 150

# One configuration is one observation. Averaging over positions first would
# discard the site-to-site variation that Section 2.6 shows to be real, so the
# correlations are evaluated on the configurations themselves.
configs = load_configurations(lattice="fcc", window=WINDOW)
per_system = defaultdict(lambda: defaultdict(list))
for system, cfg in configs.items():
    for (label, _pos), ads in cfg.items():
        for a, e in ads.items():
            per_system[system][a].append((label, e))

out = []
# --- pooled ---
elements = sorted({e for s in configs for e in SYSTEMS[s][0]})
for a in ADSORBATES:
    xs, ys = defaultdict(list), []
    for system in configs:
        for lab, e in per_system[system].get(a, []):
            comp = label_composition(lab)
            for el in elements:
                xs[el].append(comp.get(el, 0.0))
            ys.append(e)
    for el in elements:
        if len(ys) >= MIN_N and len(set(xs[el])) > 1:
            out.append({"scope": "pooled", "system": "", "element": el,
                        "adsorbate": a, "n": len(ys),
                        "r": round(pearson(xs[el], ys), 4)})
# --- within system ---
for system, tables in per_system.items():
    for a in ADSORBATES:
        obs = tables.get(a, [])
        if len(obs) < MIN_N:
            continue
        ys = [e for _l, e in obs]
        for el in SYSTEMS[system][0]:
            xs = [label_composition(l).get(el, 0.0) for l, _e in obs]
            out.append({"scope": "within", "system": system, "element": el,
                        "adsorbate": a, "n": len(obs),
                        "r": round(pearson(xs, ys), 4)})

with open(_paths.DERIVED / "elemental_correlations.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["scope", "system", "element", "adsorbate", "n", "r"])
    w.writeheader(); w.writerows(out)

pooled = [abs(r["r"]) for r in out if r["scope"] == "pooled"]
within = [abs(r["r"]) for r in out if r["scope"] == "within"]
print(f"pooled: n = {len(pooled)}, largest |r| = {max(pooled):.3f}")
print(f"within: n = {len(within)}, median |r| = {st.median(within):.3f}, "
      f"{100 * sum(v < 0.2 for v in within) / len(within):.0f}% below 0.2")
for system in sorted({r['system'] for r in out if r['scope'] == 'within'}):
    v = [abs(r["r"]) for r in out if r["system"] == system]
    print(f"  {system:12s} n {len(v):3d}  max {max(v):.3f}  "
          f"above 0.2: {100 * sum(x > 0.2 for x in v) / len(v):.0f}%")
