"""Offset-invariant descriptors.

The archived quantity carries a reference state that is undetermined to within a
per-system constant. That constant cancels exactly in any difference taken
between two adsorbates on the same configuration, so descriptors built as such
differences are quantitative within a system even though the values themselves
are not. A residual scale factor survives the cancellation and is common to all
configurations of one system, which is why these descriptors order compositions
within a system but do not compare systems.
"""
import statistics as st
from collections import defaultdict

# name -> (adsorbate A, adsorbate B); the descriptor is E_arch(A) - E_arch(B)
DESCRIPTORS = {
    "D_act":  ("CH3", "CH4"),   # first C-H activation step
    "D_coke": ("C",   "CH3"),   # full dehydrogenation to atomic carbon
    "D_rem":  ("CO",  "C"),     # removal of that carbon as CO
    "D_CO2":  ("CO2", "CO"),    # CO2 activation step
}

# sign that makes a larger contribution more favourable for dry methane reforming
SIGNS = {"D_act": -1, "D_coke": +1, "D_rem": -1, "D_CO2": +1}


def configuration_descriptors(configurations, names=None):
    """Average each descriptor over the sampled positions of every composition.

    `configurations` is the mapping returned by io.load_configurations for one
    system. Averaging over positions applies the multi-site rule: a single
    surface site accounts for no more than about two thirds of the variance of
    its own local ensemble.
    """
    names = names or list(DESCRIPTORS)
    per = defaultdict(lambda: defaultdict(list))
    for (label, _position), adsorbates in configurations.items():
        for name in names:
            a, b = DESCRIPTORS[name]
            if a in adsorbates and b in adsorbates:
                per[label][name].append(adsorbates[a] - adsorbates[b])
    rows = []
    for label, d in per.items():
        if all(name in d and d[name] for name in names):
            rec = {"composition_label": label,
                   "n_positions": max(len(v) for v in d.values())}
            rec.update({name: st.mean(d[name]) for name in names})
            rows.append(rec)
    return rows


def composite_score(rows, names=None, signs=None):
    """Equally weighted z score over the descriptors, signed so that larger is
    more favourable. Added to each row in place as `z`, and returned."""
    names = names or list(DESCRIPTORS)
    signs = signs or SIGNS
    stats = {n: (st.mean([r[n] for r in rows]), st.pstdev([r[n] for r in rows]))
             for n in names}
    for r in rows:
        r["z"] = sum(signs[n] * (r[n] - stats[n][0]) / stats[n][1]
                     for n in names) / len(names)
    return rows
