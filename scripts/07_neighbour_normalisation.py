"""Normalising the neighbour census against the sampling design.

The raw census returns Ni and Co in first place. Both occur in the same eight of
the twelve systems, so their expected counts are identical under any weighting
of the systems whatsoever, and the exact tie is the null expectation rather than
an anomaly. Putting a number on the enrichment requires assuming the systems are
equally represented in the census, which the archive does not record; the
agreement of eight independent elements with that null is the evidence for it.

Produces: data/derived/neighbour_normalisation.csv
Reproduces: Section 2.7 and the caption of Figure 4.
"""
import csv
from collections import Counter

import _paths
from compcat import SYSTEMS

# counts read off the archived neighbour census (2,742 neighbours in total)
OBSERVED = {"Ni": 360, "Co": 360, "Zn": 242, "Cu": 194, "Zr": 181,
            "Pd": 115, "Nb": 112, "Al": 110}
TOTAL = 2742

occurrences = Counter()
for elements, _lattice in SYSTEMS.values():
    occurrences.update(elements)
slots = len(SYSTEMS) * 5

rows = []
for e, n in sorted(OBSERVED.items(), key=lambda kv: -kv[1]):
    expected = 100 * occurrences[e] / slots
    observed = 100 * n / TOTAL
    rows.append({"element": e, "n_systems": occurrences[e],
                 "expected_pct": round(expected, 2),
                 "observed_pct": round(observed, 2),
                 "enrichment": round(observed / expected, 3)})

with open(_paths.DERIVED / "neighbour_normalisation.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

print(f"{'element':8s} {'systems':>8s} {'expected':>9s} {'observed':>9s} {'enrichment':>11s}")
for r in rows:
    print(f"{r['element']:8s} {r['n_systems']:8d} {r['expected_pct']:9.1f} "
          f"{r['observed_pct']:9.1f} {r['enrichment']:11.2f}")

co = {k for k, (els, _) in SYSTEMS.items() if "Co" in els}
ni = {k for k, (els, _) in SYSTEMS.items() if "Ni" in els}
print(f"\nCo and Ni occur in the same systems: {co == ni} ({len(co)} of {len(SYSTEMS)})")
print("their expected counts are therefore equal under any weighting of the systems")
