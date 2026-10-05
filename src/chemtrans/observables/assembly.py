"""Assembly-complexity observables for polymer strings.

Molecular assembly index (string version): the minimum number of join
steps in a concat-DAG producing the string, where each distinct
previously-built block may be reused at no additional cost.

    a(s) = min over binary parse DAGs of |internal nodes|

Formally computed as a shortest path over "needed block" sets:

    f(N) = 0                                          if N empty
         = 1 + min_{w in N, splits u+v=w} f(N \\ {w} u {u,v} \\ leaves)

where N is the set of distinct blocks still to be built (leaves are the
monomers, free). Shared substrings are counted once — this is the reuse
credit that distinguishes assembly index from length-1. Exact for short
strings via memoized search; for len > `exact_limit` an explicit upper
bound (len-1) is returned with `estimator="bound"` so downstream code
can never silently treat a bound as a measurement.

Signature criterion (per plan): a species only counts toward the
assembly-complexity signature if it is BOTH high-assembly (a > null
quantile for its length) AND high-abundance (mass fraction above the
calibrated abundance floor). Either condition alone is diagnostic only.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

_ALPHABET_SIZE = 2  # default; overridden per network where available


def _leaf_strings(alphabet_size: int) -> set[str]:
    return {format(i, f"0{1}b")[-1] for i in range(alphabet_size)}


def assembly_index(s: str, exact_limit: int = 12,
                   node_budget: int = 200_000) -> tuple[int, str]:
    """Return (index, estimator) where estimator is 'exact' or 'bound'."""
    if len(s) <= 1:
        return 0, "exact"
    if len(s) > exact_limit:
        return len(s) - 1, "bound"

    chars = set(s)
    budget = [node_budget]

    @lru_cache(maxsize=None)
    def f(needed: frozenset) -> int:
        if budget[0] <= 0:
            raise _BudgetExceeded
        budget[0] -= 1
        if not needed:
            return 0
        best = None
        for w in needed:
            for k in range(1, len(w)):
                u, v = w[:k], w[k:]
                nxt = (needed - {w})
                if len(u) > 1:
                    nxt = nxt | {u}
                if len(v) > 1:
                    nxt = nxt | {v}
                cost = 1 + f(frozenset(nxt))
                if best is None or cost < best:
                    best = cost
                    if best == 1:      # can't beat 1 + 0
                        break
        return best

    try:
        return f(frozenset({s})), "exact"
    except _BudgetExceeded:
        return len(s) - 1, "bound"
    finally:
        f.cache_clear()


class _BudgetExceeded(Exception):
    pass


def species_assembly(net, species_index: int,
                     **kw) -> tuple[int, str]:
    """Assembly index of species i's string label."""
    return assembly_index(net.labels[species_index], **kw)


def null_ma_distribution(rng: np.random.Generator, length: int,
                         alphabet_size: int = 2,
                         n_samples: int = 500,
                         **kw) -> np.ndarray:
    """Null MA distribution for random strings of the given length.

    The null is "random string of same length over same alphabet" — the
    threshold is therefore length-conditional, not a fixed cutoff.
    """
    out = np.empty(n_samples)
    for i in range(n_samples):
        s = "".join(rng.choice(list("01"[:alphabet_size]),
                               size=length))
        out[i], _ = assembly_index(s, **kw)
    return out


def ma_threshold(null_samples: np.ndarray, q: float = 0.95) -> float:
    """Null-calibrated assembly threshold: q-quantile of the null."""
    return float(np.quantile(null_samples, q))


def assembly_signature(
    net,
    x: np.ndarray,
    ma_thresholds: dict[int, float],
    abundance_floor: float,
    exact_limit: int = 12,
) -> dict:
    """Joint assembly+abundance signature over non-food species.

    ma_thresholds: {length: null-calibrated threshold}
    abundance_floor: minimum monomer-mass fraction of total non-food
        mass for "high abundance" (calibrated against shuffled/
        unstructured nulls upstream).

    Returns per-species records plus the signature summary:
      n_signature, mass_signature (fraction of non-food mass that is
      jointly high-MA and high-abundance).
    """
    comp_mass = net.comp.sum(axis=0)
    nonfood = ~net.is_food
    total = float(comp_mass[nonfood] @ x[nonfood])
    records = []
    n_sig = 0
    mass_sig = 0.0
    for i in range(net.n_species):
        if not nonfood[i]:
            continue
        lab = net.labels[i]
        a, est = assembly_index(lab, exact_limit=exact_limit)
        frac = float(comp_mass[i] * x[i] / total) if total > 0 else 0.0
        thr = ma_thresholds.get(len(lab))
        high_ma = est == "exact" and thr is not None and a > thr
        high_ab = frac >= abundance_floor
        sig = bool(high_ma and high_ab)
        n_sig += sig
        mass_sig += frac * sig
        records.append({
            "species": i, "label": lab, "assembly": a, "estimator": est,
            "mass_fraction": frac, "high_assembly": high_ma,
            "high_abundance": high_ab, "signature": sig,
        })
    return {"species": records, "n_signature": n_sig,
            "mass_signature": mass_sig}
