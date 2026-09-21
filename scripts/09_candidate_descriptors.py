"""Composition-resolved candidates from offset-invariant descriptors.

Produces: data/derived/descriptors_<system>.csv
Reproduces: Table 8, Figure 11 and Section 3.5, including the split-half
reproducibility test and the bootstrap null band on the elemental enrichment.

The compositional result is conditional on the record identified by
08_composition_record_test.py; both records are evaluated here.
"""
import csv
import statistics as st

import _paths
from compcat import (load_configurations, configuration_descriptors,
                     composite_score, label_composition, pearson,
                     bootstrap_null, split_half, DESCRIPTORS)
from compcat.io import load_parsed_compositions
from compcat.window import WINDOW

SYSTEMS = ["CoNiPdRhRu", "AgAuCuPdPt"]
configs = load_configurations(systems=set(SYSTEMS), window=WINDOW)

for system in SYSTEMS:
    cfg = configs.get(system, {})
    rows = composite_score(configuration_descriptors(cfg))
    rows.sort(key=lambda r: -r["z"])
    print(f"\n=== {system}: {len(rows)} compositions with the complete descriptor set ===")
    for name in DESCRIPTORS:
        v = [r[name] for r in rows]
        print(f"  {name:7s} mean {st.mean(v):+6.2f}  sd {st.pstdev(v):5.2f}")
    print("  top five")
    for r in rows[:5]:
        print("   ", r["composition_label"],
              " ".join(f"{n} {r[n]:+.2f}" for n in DESCRIPTORS), f"z {r['z']:+.2f}")

    with open(_paths.DERIVED / f"descriptors_{system.lower()}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["composition_label", "n_positions",
                                           *DESCRIPTORS, "z"])
        w.writeheader(); w.writerows(rows)

    # --- elemental enrichment of the leading group, under both records ---
    elements = sorted(label_composition(rows[0]["composition_label"]))
    records = [("directory_label", label_composition)]
    if system == "CoNiPdRhRu":
        parsed = load_parsed_compositions(system)
        rows_p = [r for r in rows if r["composition_label"] in parsed]
        records.append(("parsed_field", lambda lab: parsed[lab]))
    for record, getter in records:
        pool = [r for r in rows if record == "directory_label"
                or r["composition_label"] in parsed]
        top = pool[:30]
        print(f"  enrichment of the top 30 - {record}")
        for e in elements:
            base = st.mean([getter(r["composition_label"]).get(e, 0.0) for r in pool]) * 100
            obs = st.mean([getter(r["composition_label"]).get(e, 0.0) for r in top]) * 100 - base
            lo, hi = bootstrap_null(
                [getter(r["composition_label"]).get(e, 0.0) for r in pool],
                lambda s: st.mean(s) * 100 - base, 30, n_resamples=2000, seed=1)
            r_el = pearson([getter(r["composition_label"]).get(e, 0.0) for r in pool],
                           [r["z"] for r in pool])
            flag = "significant" if (obs < lo or obs > hi) else "-"
            print(f"    {e:3s} {obs:+5.2f} pp   null [{lo:+.2f}, {hi:+.2f}]   "
                  f"r {r_el:+.3f}   {flag}")

# --- split-half reproducibility of the CoNiPdRhRu ranking ---
cfg = configs["CoNiPdRhRu"]
positions = sorted({p for _lab, p in cfg})
half_a = {k: v for k, v in cfg.items() if k[1] in positions[:4]}
half_b = {k: v for k, v in cfg.items() if k[1] in positions[4:]}
sa = {r["composition_label"]: r["z"] for r in composite_score(configuration_descriptors(half_a))}
sb = {r["composition_label"]: r["z"] for r in composite_score(configuration_descriptors(half_b))}
r, shared, expected = split_half(sorted(set(sa) & set(sb)), sa, sb, top=100)
print(f"\nsplit half: r = {r:.3f}; {shared} of the top 100 shared "
      f"against {expected:.0f} expected by chance")
