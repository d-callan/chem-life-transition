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

Signature criterion: a species counts toward the assembly-complexity
signature iff its EXACT assembly index is >= `a_cut` AND its absolute
monomer-mass concentration is >= `mass_floor` — see
`assembly_signature_absolute`. `a_cut` is either explicit, derived from
the pool's maximal exact class, or derived from a measured abiotic null
(`derive_a_cut`): under a sequence-blind null each length-L string is
equiprobable, so the expected copy number of any specific string is
V * M_L^null * |Sigma|^{-L} with M_L^null ~ c*r^L; the cut is set where
that expectation falls below one copy.
"""

from __future__ import annotations

from functools import cache

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

    budget = [node_budget]

    @cache
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
    return assembly_index(net.species[species_index], **kw)


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


def pool_ma_distribution(net, exact_limit: int = 12) -> dict:
    """Exact assembly-index distribution over the non-food species pool.

    Returns {a_max_exact, counts: {index: n_species}}. Counts only exact
    estimates; species returning "bound" are tallied separately.
    """
    counts: dict[int, int] = {}
    n_bound = 0
    amax = 0
    for i in range(net.n_species):
        if net.is_food[i]:
            continue
        a, est = assembly_index(net.species[i], exact_limit=exact_limit)
        if est == "exact":
            counts[a] = counts.get(a, 0) + 1
            amax = max(amax, a)
        else:
            n_bound += 1
    return {"a_max_exact": amax, "counts": {str(k): v
            for k, v in sorted(counts.items())},
            "n_bound": n_bound}


def assembly_signature_absolute(
    net,
    x: np.ndarray,
    a_cut: int | None = None,
    mass_floor: float = 0.0,
    exact_limit: int = 12,
) -> dict:
    """Absolute-threshold assembly signature over non-food species.

    A species is a signature hit iff its EXACT assembly index is
    >= `a_cut` AND its absolute monomer-mass concentration is
    >= `mass_floor`. `a_cut=None` derives the pool's maximal exact
    index (top complexity class). `mass_floor` is an absolute
    concentration, NOT a fraction — relative fractions spike when total
    polymer mass collapses.

    Returns the full (a, mass) spectrum plus:
      a_cut, mass_floor, pool_ma_max, n_signature,
      mass_signature_abs (absolute monomer-mass concentration of
      signature species), mass_signature_frac (diagnostic only).
    """
    comp_mass = net.comp.sum(axis=0)
    nonfood = ~net.is_food
    total = float(comp_mass[nonfood] @ x[nonfood])
    dist = pool_ma_distribution(net, exact_limit=exact_limit)
    cut = dist["a_max_exact"] if a_cut is None else int(a_cut)
    records, n_sig, mass_sig, mass_sig_frac = [], 0, 0.0, 0.0
    for i in range(net.n_species):
        if not nonfood[i]:
            continue
        lab = net.species[i]
        a, est = assembly_index(lab, exact_limit=exact_limit)
        mass_abs = float(comp_mass[i] * x[i])
        frac = mass_abs / total if total > 0 else 0.0
        high_ma = est == "exact" and a >= cut
        high_ab = mass_abs >= mass_floor
        sig = bool(high_ma and high_ab)
        n_sig += sig
        mass_sig += mass_abs * sig
        mass_sig_frac += frac * sig
        records.append({
            "species": i, "label": lab, "assembly": a, "estimator": est,
            "mass_abs": mass_abs, "mass_fraction": frac,
            "high_assembly": high_ma, "high_abundance": high_ab,
            "signature": sig,
        })
    return {"species": records, "a_cut": cut, "mass_floor": mass_floor,
            "pool_ma_max": dist["a_max_exact"],
            "pool_ma_dist": dist["counts"], "n_pool_bound": dist["n_bound"],
            "n_signature": n_sig, "mass_signature_abs": mass_sig,
            "mass_signature_frac": mass_sig_frac}


def length_mass_spectrum(net, x: np.ndarray) -> dict[int, float]:
    """{length: absolute nonfood monomer-mass concentration} — the
    length-resolved spectrum needed for the null fit. Requires
    net.structure (species labels)."""
    comp_mass = net.comp.sum(axis=0)
    out: dict[int, float] = {}
    for i in range(net.n_species):
        if net.is_food[i]:
            continue
        lab = net.structure.get(i) if net.structure is not None else None
        if not isinstance(lab, str):
            continue
        L = len(lab)
        out[L] = out.get(L, 0.0) + float(comp_mass[i] * x[i])
    return dict(sorted(out.items()))


def derive_a_cut(length_masses: dict[int, float], volume: float,
                 alphabet_size: int = 2, copies_floor: float = 1.0,
                 ma_to_length: int = 2, max_L: int = 64) -> dict:
    """Null-derived assembly cut (the toy Walker analogue).

    Under a sequence-blind abiotic null, each length-L string is
    equiprobable (P = alphabet_size^{-L}), so the expected copy number of
    any specific string is

        E[n_s] = volume * M_L^null * alphabet_size^{-L} / L

    (M_L is a monomer-MASS concentration; one length-L molecule carries
    L monomers, so per-sequence molecule concentration is M_L / (L *
    |Sigma|^L)).

    with M_L^null ~ c * r^L fit by log-linear regression on the measured
    length spectrum.  L* is the first length with E < copies_floor, and

        a_cut = L* - ma_to_length

    (binary strings empirically reach max MA ~ L-2; the offset converts
    the length threshold into an index threshold).

    Returns a derivation record suitable for a campaign manifest:
    c, r, L_star, a_cut, per-length expected copies, plus `extrapolated`
    (L* beyond the measured range) and `reason` when the fit is
    impossible. The cut is only as good as the null: if the null becomes
    sequence-dependent, alphabet_size^{-L} must be replaced by the
    measured P_null(s) and this derivation redone.
    """
    Ls = np.array(sorted(length_masses), dtype=float)
    Ms = np.array([length_masses[int(L)] for L in Ls], dtype=float)
    pos = Ms > 0
    rec = {"volume": float(volume), "alphabet_size": int(alphabet_size),
           "copies_floor": float(copies_floor),
           "ma_to_length": int(ma_to_length),
           "measured_spectrum": {str(int(L)): float(m)
                                 for L, m in zip(Ls, Ms)},
           "n_lengths_fit": int(pos.sum())}
    if pos.sum() < 2:
        rec.update({"a_cut": None, "L_star": None,
                    "reason": "fewer than two nonzero length classes in "
                              "the null spectrum"})
        return rec
    slope, intercept = np.polyfit(Ls[pos], np.log(Ms[pos]), 1)
    r, c = float(np.exp(slope)), float(np.exp(intercept))
    rec["c"], rec["r"] = c, r

    # expected specific-copy number vs length; first L below the floor
    expected = {}
    L_star = None
    for L in range(1, max_L + 1):
        e = volume * c * r ** L * alphabet_size ** (-L) / float(L)
        expected[str(L)] = float(e)
        if e < copies_floor:
            L_star = L
            break
    if L_star is None:
        rec.update({"a_cut": None, "L_star": None,
                    "expected_copies": expected,
                    "reason": "expected copies never fall below the "
                              "floor within max_L"})
        return rec
    a_cut = int(L_star - ma_to_length)
    rec.update({
        "L_star": int(L_star), "a_cut": a_cut,
        "expected_copies": expected,
        "extrapolated": bool(L_star > Ls.max()),
        # a_cut <= 0: the null is so dilute that no specific string is
        # resolvable at this volume — the signature is vacuous, not
        # permissive. Increase V or the abiotic rate before interpreting.
        "vacuous": bool(a_cut <= 0),
        "derivation": ("L* = first L with V*c*r^L*|Sigma|^-L/L < "
                       "copies_floor; a_cut = L* - ma_to_length"),
    })
    return rec
