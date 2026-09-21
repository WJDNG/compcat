"""The statistics used in the paper: Pearson r, bootstrap intervals, split-half.

Written with the standard library only so that the reproduction scripts run
without a scientific stack if necessary.
"""
import math
import random
import statistics as st


def pearson(x, y):
    """Pearson correlation. Returns nan when either series is constant."""
    x, y = list(x), list(y)
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("series must be the same length and hold at least two points")
    mx, my = st.mean(x), st.mean(y)
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx == 0 or syy == 0:
        return float("nan")
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(sxx * syy)


def bootstrap_ci(sample, statistic, n_resamples=2000, level=0.95, seed=0):
    """Percentile bootstrap interval for `statistic` computed on `sample`."""
    rng = random.Random(seed)
    n = len(sample)
    draws = []
    for _ in range(n_resamples):
        draws.append(statistic([sample[rng.randrange(n)] for _ in range(n)]))
    draws.sort()
    k = int((1 - level) / 2 * n_resamples)
    return draws[k], draws[-(k + 1)]


def bootstrap_null(population, statistic, size, n_resamples=2000, level=0.95, seed=0):
    """Null band for `statistic` evaluated on `size` members drawn at random
    from `population`. Used to test whether the top-scoring group differs from a
    random group of the same size."""
    rng = random.Random(seed)
    n = len(population)
    draws = []
    for _ in range(n_resamples):
        draws.append(statistic([population[rng.randrange(n)] for _ in range(size)]))
    draws.sort()
    k = int((1 - level) / 2 * n_resamples)
    return draws[k], draws[-(k + 1)]


def split_half(keys, score_a, score_b, top=100):
    """Reproducibility of a ranking across two independent halves of the data.

    Returns (r, shared, expected_by_chance) where `shared` counts the members
    common to the two top-`top` lists.
    """
    common = [k for k in keys if k in score_a and k in score_b]
    r = pearson([score_a[k] for k in common], [score_b[k] for k in common])
    ta = set(sorted(common, key=lambda k: -score_a[k])[:top])
    tb = set(sorted(common, key=lambda k: -score_b[k])[:top])
    return r, len(ta & tb), top * top / len(common)
