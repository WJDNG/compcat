"""The physical window applied to the archived values.

The archived tables contain entries lying orders of magnitude outside any
physically meaningful range. The archive records no convergence flag and retains
no trajectories, so the cause of an individual outlier cannot be recovered; the
window is applied as a filter on non-physical values and nothing is inferred
about why a given entry falls outside it.
"""

WINDOW = (-10.0, 5.0)   # eV, on the archived quantity E_arch


def in_window(value, window=WINDOW):
    lo, hi = window
    return lo <= value <= hi


def retention(values, window=WINDOW):
    """Fraction of `values` inside the window."""
    values = list(values)
    if not values:
        return 0.0
    return sum(in_window(v, window) for v in values) / len(values)
