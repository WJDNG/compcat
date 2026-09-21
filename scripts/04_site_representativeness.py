"""Single-site versus multi-site sampling on CoNiPdRhRu(111).

Produces: data/derived/site_statistics.csv
Reproduces: Table S5, Figure 3, and the R2 = 0.50-0.65 result of Section 2.6,
with bootstrap intervals (largest upper bound 0.674).
"""
import csv
import statistics as st
from collections import defaultdict

import _paths
from compcat import load_energy_table, pearson

try:                      # the bootstrap is the only step that wants numpy
    import numpy as np
except ImportError:
    np = None

def _bootstrap_r2(xs, ys, n_resamples=400, seed=0):
    """Percentile interval on R2. Uses numpy when available; the pure-Python
    fallback is correct but slow, so it runs fewer resamples."""
    if np is None:
        import random
        rng = random.Random(seed)
        n = len(xs)
        draws = []
        for _ in range(min(n_resamples, 100)):
            idx = [rng.randrange(n) for _ in range(n)]
            draws.append(pearson([xs[i] for i in idx], [ys[i] for i in idx]) ** 2)
        draws.sort()
        k = int(0.025 * len(draws))
        return draws[k], draws[-(k + 1)]
    x = np.asarray(xs); y = np.asarray(ys); n = x.size
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_resamples, n))
    bx, by = x[idx], y[idx]
    bx = bx - bx.mean(axis=1, keepdims=True)
    by = by - by.mean(axis=1, keepdims=True)
    r = (bx * by).sum(axis=1) / np.sqrt((bx ** 2).sum(axis=1) * (by ** 2).sum(axis=1))
    draws = np.sort(r ** 2)
    k = int(0.025 * n_resamples)
    return float(draws[k]), float(draws[-(k + 1)])


by_slab = defaultdict(dict)
for row in load_energy_table("df_ads.csv.gz"):
    if row["system"] != "CoNiPdRhRu" or row["adsorbate"] != "CH4":
        continue
    by_slab[row["slab_index"]][row["position"]] = row["energy"]

full = {k: v for k, v in by_slab.items() if len(v) == 8}
positions = sorted({p for v in full.values() for p in v})
print(f"slabs carrying all eight positions: {len(full)}")

# statistics.mean is exact-arithmetic and far too slow for this many calls
def _mean(seq):
    seq = list(seq)
    return sum(seq) / len(seq)


slabs = list(full.values())
totals = [sum(v.values()) for v in slabs]

rows = []
for p in positions:
    xs = [(t - v[p]) / 7.0 for t, v in zip(totals, slabs)]   # mean of the other seven
    ys = [v[p] for v in slabs]
    r = pearson(xs, ys)
    mx, my = _mean(xs), _mean(ys)
    slope = (sum((a - mx) * (b - my) for a, b in zip(xs, ys))
             / sum((a - mx) ** 2 for a in xs))
    lo, hi = _bootstrap_r2(xs, ys, n_resamples=400, seed=0)
    rows.append({"position": p, "n": len(ys), "mean_Earch": round(my, 4),
                 "r": round(r, 4), "R2": round(r * r, 4),
                 "R2_ci_low": round(lo, 4), "R2_ci_high": round(hi, 4),
                 "slope": round(slope, 3)})

rows.sort(key=lambda r: -r["R2"])
with open(_paths.DERIVED / "site_statistics.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

means = [r["mean_Earch"] for r in rows]
print(f"R2 range {min(r['R2'] for r in rows):.3f}-{max(r['R2'] for r in rows):.3f}; "
      f"largest upper 95% bound {max(r['R2_ci_high'] for r in rows):.3f}")
print(f"slopes {min(r['slope'] for r in rows):.2f}-{max(r['slope'] for r in rows):.2f}; "
      f"spread of site means {max(means) - min(means):.3f} eV")
