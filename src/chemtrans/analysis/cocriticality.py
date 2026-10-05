"""Co-criticality: do independently-defined life-observables share λc?

Pipeline:
  1. extract per-network (lam, O_k) series from life_vector components
  2. per-observable critical-point estimation:
       continuous -> fit_all_models model selection, λc from the
                     selected model (none if null/ramp wins)
       discrete   -> changepoint: earliest λ where P(O > 0) >= 0.5
                     across networks, with binomial uncertainty
  3. paired network-joint bootstrap: resample WHOLE networks, refit
     every observable on the same resample, take differences
     D_jk = λc_j - λc_k. Pairwise verdicts:

       "indistinguishable" : CI95(D_jk) contains 0 — failure to
                             reject equality, NOT proof (report width)
       "separated"         : CI95(D_jk) excludes 0
       "incomparable"      : one or both observables lack a defined
                             λc in >= `min_defined_frac` of resamples
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .critical_points import fit_all_models

CONTINUOUS = ("O_raf", "O_maint", "O_persist", "O_select", "O_assembly")
DISCRETE = ("O_hered",)

# map selected model -> (params key holding λc, transition class)
_LC_MAP = {"hinge": ("lc", "onset"),
           "itrans": ("lc", "onset_imperfect"),
           "fold": ("lf", "saddle_node"),
           "sigmoid": ("lc", "crossover")}


@dataclass
class LcEstimate:
    lc: float | None
    transition_class: str      # onset/onset_imperfect/saddle_node/
                               # crossover/none/discrete_changepoint
    detail: dict


def group_by_network(records: list[dict]) -> dict[str, list[dict]]:
    """records: run dicts carrying a 'life_vector' or raw blocks."""
    by_net: dict[str, list[dict]] = {}
    for r in records:
        by_net.setdefault(r["net_hash"], []).append(r)
    return by_net


def extract_series(by_net: dict[str, list[dict]],
                   obs: str) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """{net_hash: (lam, y)} keeping only records where O_k is defined.

    Requires each record's life_vector (call observables.life_vector
    first and attach as rec["life_vector"]) and a scalar control value
    at rec["control"]["lam"].
    """
    out = {}
    for h, recs in by_net.items():
        lam, y = [], []
        for r in recs:
            lv = r.get("life_vector", {})
            comp = lv.get(obs, {})
            if comp.get("defined") and comp["value"] is not None:
                lam.append(float(r["control"]["lam"]))
                y.append(float(comp["value"]))
        if len(lam) >= 3:
            order = np.argsort(lam)
            out[h] = (np.asarray(lam)[order], np.asarray(y)[order])
    return out


def estimate_lc_continuous(lam: np.ndarray, y: np.ndarray,
                           seed: int = 0) -> LcEstimate:
    fit = fit_all_models(lam, y, seed=seed)
    sel, params = fit["selected"], fit["params"]
    if sel in _LC_MAP:
        key, cls = _LC_MAP[sel]
        return LcEstimate(float(params[key]), cls,
                          {"selected": sel, "aicc": fit["aicc"],
                           "delta_aicc": fit["delta_aicc"]})
    return LcEstimate(None, "none" if sel == "null" else "crossover",
                      {"selected": sel, "aicc": fit["aicc"]})


def estimate_lc_discrete(by_net_series: dict, active_pred,
                         threshold: float = 0.5) -> LcEstimate:
    """Changepoint for a discrete observable (e.g. O_hered > 0).

    λc = smallest λ where the fraction of networks active >=
    `threshold`. Per-point Clopper-Pearson-ish uncertainty reported as
    the λ-range over which the active fraction straddles `threshold`.
    """
    grid = sorted({v for lam, _ in by_net_series.values() for v in lam})
    fracs = []
    for lv in grid:
        act = sum(1 for lam, y in by_net_series.values()
                  for l0, y0 in zip(lam, y) if l0 == lv and active_pred(y0))
        tot = sum(1 for lam, _ in by_net_series.values()
                  for l0 in lam if l0 == lv)
        fracs.append(act / tot if tot else np.nan)
    fracs = np.asarray(fracs)
    hits = [grid[i] for i, f in enumerate(fracs)
            if np.isfinite(f) and f >= threshold]
    if not hits:
        return LcEstimate(None, "none",
                          {"active_fraction": dict(zip(grid, fracs.tolist()))})
    lc = hits[0]
    # uncertainty: neighboring grid values where frac < threshold
    below = [g for g, f in zip(grid, fracs)
             if g < lc and np.isfinite(f)]
    above = [g for g, f in zip(grid, fracs)
             if g > lc and np.isfinite(f) and f < threshold]
    return LcEstimate(lc, "discrete_changepoint",
                      {"lc_lo": below[-1] if below else None,
                       "lc_hi": above[0] if above else None,
                       "active_fraction": dict(zip(grid,
                                                   fracs.tolist()))})


def _lc_for_obs(by_net: dict[str, list[dict]], obs: str,
                seed: int) -> LcEstimate:
    series = extract_series(by_net, obs)
    if not series:
        return LcEstimate(None, "undefined", {})
    if obs in DISCRETE:
        return estimate_lc_discrete(series, lambda v: v > 0)
    lam = np.concatenate([s[0] for s in series.values()])
    y = np.concatenate([s[1] for s in series.values()])
    return estimate_lc_continuous(lam, y, seed=seed)


def coincidence_report(by_net: dict[str, list[dict]],
                       observables: tuple[str, ...] = CONTINUOUS + DISCRETE,
                       n_boot: int = 400, seed: int = 0,
                       min_defined_frac: float = 0.5) -> dict:
    """Full pairwise coincidence matrix over observables.

    Bootstrap unit = whole network (all its λ-points and all
    observables resampled together). λc is refit per resample.
    """
    rng = np.random.default_rng(seed)
    hashes = sorted(by_net)
    n = len(hashes)

    point_ests = {o: _lc_for_obs(by_net, o, seed) for o in observables}

    boot_lc: dict[str, list] = {o: [] for o in observables}
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        resampled = {f"{hashes[i]}#{b}": by_net[hashes[i]]
                     for b, i in enumerate(idx)}
        s = int(rng.integers(2**31))
        for o in observables:
            boot_lc[o].append(_lc_for_obs(resampled, o, s).lc)

    matrix = {}
    for j, oj in enumerate(observables):
        for ok in observables[j + 1:]:
            a = np.asarray([v if v is not None else np.nan
                            for v in boot_lc[oj]])
            b = np.asarray([v if v is not None else np.nan
                            for v in boot_lc[ok]])
            both = np.isfinite(a) & np.isfinite(b)
            if both.mean() < min_defined_frac:
                verdict, d = "incomparable", np.array([np.nan])
            else:
                d = a[both] - b[both]
                lo, hi = np.quantile(d, [0.025, 0.975])
                verdict = ("separated" if lo > 0 or hi < 0
                           else "indistinguishable")
            matrix[f"{oj}|{ok}"] = {
                "verdict": verdict,
                "frac_both_defined": float(both.mean()),
                "ci95_width": (float(np.quantile(d, 0.975)
                                     - np.quantile(d, 0.025))
                               if np.isfinite(d).any() else None),
                "d_ci95": ([float(np.quantile(d, 0.025)),
                            float(np.quantile(d, 0.975))]
                           if np.isfinite(d).any() else None),
            }
    return {
        "point_estimates": {
            o: {"lc": e.lc, "class": e.transition_class,
                "detail": e.detail}
            for o, e in point_ests.items()},
        "pairwise": matrix,
        "n_boot": n_boot, "n_networks": n,
        "note": ("'indistinguishable' = failure to reject equality; "
                 "compare ci95_width to the sweep resolution"),
    }
