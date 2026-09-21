"""Which of the two composition records the archive supports.

The archive holds two records of the slab composition: the trajectory directory
label, and a composition field in the post-processed tables. This script
measures how far apart they are and correlates each against E_arch. One record
produces a chemically ordered pattern, the other produces nothing; the paper
takes the first as the operative composition and says so explicitly.

Produces: data/derived/composition_record_test.csv
Reproduces: the provenance paragraph of Section 3.5.
"""
import csv
import statistics as st
from collections import defaultdict


def _mean(seq):
    """statistics.mean uses exact arithmetic and is far too slow here."""
    seq = list(seq)
    return sum(seq) / len(seq)

import _paths
from compcat import (load_configurations, label_composition, pearson)
from compcat.io import load_parsed_compositions
from compcat.window import WINDOW

SYSTEM = "CoNiPdRhRu"
ELEMENTS = ["Co", "Ni", "Pd", "Rh", "Ru"]
ADSORBATES = ["CH3", "OH", "CO", "O", "C", "H"]

parsed = load_parsed_compositions(SYSTEM)
configs = load_configurations(systems={SYSTEM}, window=WINDOW)[SYSTEM]

# --- how far apart are the two records ---
dev = defaultdict(list)
for label in parsed:
    a, b = parsed[label], label_composition(label)
    for e in b:
        dev[e].append(100 * (a.get(e, 0.0) - b[e]))
print("deviation of the parsed field from the directory label, percentage points")
for e in ELEMENTS:
    print(f"  {e:3s} mean {st.mean(dev[e]):+5.2f}   mean |dev| "
          f"{st.mean(map(abs, dev[e])):5.2f}   max {max(map(abs, dev[e])):5.2f}")

# --- correlate each record against E_arch ---
values = defaultdict(lambda: defaultdict(list))
for (label, _pos), ads in configs.items():
    for a, e in ads.items():
        values[a][label].append(e)

out = []
for record, getter in (("directory_label", label_composition),
                       ("parsed_field", lambda lab: parsed[lab])):
    print(f"\nPearson r of element fraction with E_arch - {record}")
    print("      " + " ".join(f"{a:>7s}" for a in ADSORBATES))
    for el in ELEMENTS:
        line = []
        for a in ADSORBATES:
            labels = [l for l in values[a] if l in parsed]
            r = pearson([getter(l).get(el, 0.0) for l in labels],
                        [_mean(values[a][l]) for l in labels])
            line.append(r)
            out.append({"record": record, "element": el, "adsorbate": a, "r": round(r, 4)})
        print(f"{el:5s} " + " ".join(f"{v:+7.3f}" for v in line))

with open(_paths.DERIVED / "composition_record_test.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["record", "element", "adsorbate", "r"])
    w.writeheader(); w.writerows(out)
