"""Critical-point / transition-shape inference.

Fits candidate response curves to order-parameter-vs-control data and
selects by corrected AIC. Continuous fits only — no grid-snapped
critical points, no "first index where O jumps" estimators.

Models (functions of control lam):

  null:    O = c0
  ramp:    O = c0 + c1*lam   (k=2; an explicit lam0 is non-identifiable)
  sigmoid: O = c0 + A / (1 + exp(-(lam - lc)/w))
  hinge:   O = c0 + A * max(0, lam - lc)                (continuous onset)
  fold:    O = c0           (lam < lf)
           O = c1 + A*sqrt(lam - lf)   (lam >= lf)      (saddle-node jump)
  itrans:  O = c0 + (A/2)[(lam-lc) + sqrt((lam-lc)^2 + 4*c_field)]
           imperfect/transcritical bifurcation with rounding field
           c_field ~ eps; mean-field exponents beta=1, delta=2.

Selection reports the full AICc table — a flat/null win IS a valid
outcome (honest null). OSCILLATORY is a named outcome of dynamics,
not a fit model; trajectories labeled oscillatory must be excluded
from these fits by the caller and reported separately.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class FitResult:
    model: str
    params: dict
    rss: float
    aicc: float
    n: int
    predict: callable


def _model_preds():
    return {
        "null": (lambda p, L: np.full_like(L, p["c0"]),
                 {"c0": None},
                 [{"c0": m} for m in np.linspace(0, 1, 3)]),
        "ramp": (lambda p, L: p["c0"] + p["c1"] * L,
                 {"c0": None, "c1": None},
                 [{"c0": 0.1, "c1": s} for s in (0.5, 2.0)]),
        "sigmoid": (lambda p, L: p["c0"] + p["A"] /
                    (1.0 + np.exp(np.clip(-(L - p["lc"]) / p["w"],
                                          -60.0, 60.0))),
                    {"c0": None, "A": None, "lc": None, "w": None},
                    [{"c0": 0.0, "A": 0.8, "lc": None, "w": w}
                     for w in (0.02, 0.1, 0.3)]),
        "hinge": (lambda p, L: p["c0"] + p["A"] * np.maximum(
                      0.0, L - p["lc"]),
                  {"c0": None, "A": None, "lc": None},
                  [{"c0": 0.02, "A": a, "lc": None}
                   for a in (0.5, 2.0)]),
        "fold": (lambda p, L: np.where(
                     L < p["lf"], p["c0"],
                     p["c1"] + p["A"] * np.sqrt(
                         np.maximum(0.0, L - p["lf"]))),
                 {"c0": None, "c1": None, "A": None, "lf": None},
                 [{"c0": 0.05, "c1": 0.6, "A": 0.8, "lf": None}]),
        "itrans": (lambda p, L: p["c0"] + 0.5 * p["A"] * (
                       (L - p["lc"]) + np.sqrt(np.maximum(
                           (L - p["lc"]) ** 2 + 4.0 * p["c_field"],
                           0.0))),
                   {"c0": None, "A": None, "lc": None,
                    "c_field": None},
                   [{"c0": 0.02, "A": 1.0, "lc": None,
                     "c_field": cf} for cf in (1e-4, 1e-2, 0.1)]),
        # mirrored (offset) variants: the same transition shapes
        # opening toward decreasing L. Without these, a monotone
        # decreasing curve (order dying as lam grows) can only be fit
        # by the direction-agnostic sigmoid/ramp — a structural AICc
        # bias, not evidence for crossover physics. Same dof per pair,
        # so the penalty is symmetric.
        "hinge_rev": (lambda p, L: p["c0"] + p["A"] * np.maximum(
                          0.0, p["lc"] - L),
                      {"c0": None, "A": None, "lc": None},
                      [{"c0": 0.02, "A": a, "lc": None}
                       for a in (0.5, 2.0)]),
        "fold_rev": (lambda p, L: np.where(
                         L > p["lf"], p["c0"],
                         p["c1"] + p["A"] * np.sqrt(
                             np.maximum(0.0, p["lf"] - L))),
                     {"c0": None, "c1": None, "A": None, "lf": None},
                     [{"c0": 0.05, "c1": 0.6, "A": 0.8, "lf": None}]),
        "itrans_rev": (lambda p, L: p["c0"] + 0.5 * p["A"] * (
                           (p["lc"] - L) + np.sqrt(np.maximum(
                               (p["lc"] - L) ** 2 + 4.0 * p["c_field"],
                               0.0))),
                       {"c0": None, "A": None, "lc": None,
                        "c_field": None},
                       [{"c0": 0.02, "A": 1.0, "lc": None,
                         "c_field": cf} for cf in (1e-4, 1e-2, 0.1)]),
    }


def _free_keys(model: str) -> list[str]:
    return list(_model_preds()[model][1].keys())


def _bounds(model: str, lam: np.ndarray):
    """Physical parameter bounds for trust-region fitting, in
    _free_keys(model) order. Critical points may lie up to 0.5 outside
    the sampled window (onset between grid points) but not arbitrarily
    far; rounding fields and widths are nonnegative."""
    lo = float(lam.min()) - 0.5
    hi = float(lam.max()) + 0.5
    return {
        "ramp": ([-np.inf, -np.inf], [np.inf, np.inf]),
        "sigmoid": ([-np.inf, -np.inf, lo, 1e-4],
                    [np.inf, np.inf, hi, 10.0]),
        "hinge": ([-np.inf, 0.0, lo], [np.inf, np.inf, hi]),
        "fold": ([-np.inf, -np.inf, 0.0, lo],
                 [np.inf, np.inf, np.inf, hi]),
        "itrans": ([-np.inf, 0.0, lo, 0.0],
                   [np.inf, np.inf, hi, 100.0]),
        "hinge_rev": ([-np.inf, 0.0, lo], [np.inf, np.inf, hi]),
        "fold_rev": ([-np.inf, -np.inf, 0.0, lo],
                     [np.inf, np.inf, np.inf, hi]),
        "itrans_rev": ([-np.inf, 0.0, lo, 0.0],
                       [np.inf, np.inf, hi, 100.0]),
    }[model]


def _fit_one(model: str, lam: np.ndarray, y: np.ndarray,
             weights: np.ndarray | None = None,
             seed: int = 0) -> FitResult:
    from scipy.optimize import least_squares
    pred_fn, _, starts = _model_preds()[model]
    keys = _free_keys(model)
    w = np.ones_like(y) if weights is None else np.asarray(weights, float)

    def unpack(v):
        return dict(zip(keys, v))

    best = None
    # expand starts: an unspecified breakpoint (lc/lf) is tried at
    # several quantiles of lam — the non-smooth landscape around the
    # breakpoint makes bounded TRF far more start-sensitive than LM was
    import itertools
    start_grid: list[np.ndarray] = []
    for start in starts:
        per_key = []
        for k in keys:
            if k in ("lc", "lf"):
                per_key.append(
                    [start[k]] if start.get(k) is not None
                    else list(np.quantile(lam, (0.25, 0.5, 0.75))))
            elif k == "c0" and start.get(k) is None:
                per_key.append([float(np.min(y))])
            else:
                per_key.append([start.get(k, 0.1)])
        for combo in itertools.product(*per_key):
            start_grid.append(np.array(combo, float))
    for p0 in start_grid:
        if model == "null":
            theta = np.array([np.average(y, weights=w)])
            resid = (y - theta[0]) * np.sqrt(w)
            rss = float(resid @ resid)
            cand = (theta, rss)
        else:
            def r(v):
                return (pred_fn(unpack(v), lam) - y) * np.sqrt(w)
            try:
                lb, ub = _bounds(model, lam)
                sol = least_squares(
                    r, np.clip(p0, lb, ub), bounds=(lb, ub),
                    method="trf", max_nfev=4000)
            except Exception:  # noqa: BLE001,S112 - a failed start skips to the next grid point
                continue
            rss = float(sol.fun @ sol.fun)
            cand = (sol.x, rss)
        if best is None or cand[1] < best[1]:
            best = cand
    if best is None:
        return FitResult(model, {}, np.inf, np.inf, len(y),
                         lambda L: np.full_like(np.asarray(L, float),
                                                np.nan))
    theta, rss = best
    params = unpack(theta)
    if model == "null":
        params = {"c0": float(theta[0])}
    k = len(keys) if model != "null" else 1
    n = len(y)
    aicc = (n * np.log(max(rss / n, 1e-300))
            + 2 * k + 2 * k * (k + 1) / max(n - k - 1, 1))
    return FitResult(model, params, rss, float(aicc), n,
                     lambda L, fn=pred_fn, pp=params:
                     fn(pp, np.asarray(L, float)))


def fit_all_models(lam: np.ndarray, y: np.ndarray,
                   weights: np.ndarray | None = None,
                   models: tuple[str, ...] | None = None,
                   seed: int = 0) -> dict:
    """Fit every candidate; return AICc table + selected model."""
    lam = np.asarray(lam, float)
    y = np.asarray(y, float)
    ok = np.isfinite(lam) & np.isfinite(y)
    lam, y = lam[ok], y[ok]
    if weights is not None:
        weights = np.asarray(weights, float)[ok]
    names = list(models) if models else list(_model_preds().keys())
    fits = {m: _fit_one(m, lam, y, weights, seed) for m in names}
    order = sorted(fits, key=lambda m: fits[m].aicc)
    best = fits[order[0]]
    return {
        "selected": order[0],
        "params": best.params,
        "aicc": {m: fits[m].aicc for m in order},
        "rss": {m: fits[m].rss for m in order},
        "n": len(y),
        "delta_aicc": {m: fits[m].aicc - best.aicc for m in order},
        "_fits": fits,
    }


# models that carry a critical-point parameter, and its name
LC_PARAM = {"hinge": "lc", "itrans": "lc", "fold": "lf",
            "sigmoid": "lc",
            "hinge_rev": "lc", "itrans_rev": "lc", "fold_rev": "lf"}


def lc_from_fit(fit: dict) -> float | None:
    """Critical-point estimate from a fit_all_models result.

    Returns the selected model's breakpoint parameter, or None when
    the winner carries no changepoint (null/ramp)."""
    key = LC_PARAM.get(fit["selected"])
    if key is None:
        return None
    v = fit["params"].get(key)
    return float(v) if v is not None else None


def network_joint_bootstrap(
    records: list[dict],
    lam_key: str, y_key: str,
    n_boot: int = 400, seed: int = 0,
) -> dict:
    """Bootstrap over NETWORK records, preserving paired observables.

    A record is one network's (possibly multi-lambda) measurement set.
    Resampling whole records preserves the nested trajectory structure
    and the pairing between lambda and O — unlike resampling points.

    ANTI-POOLING: records are never concatenated into one scatter —
    each resampled record is fit separately, so the output is the
    dispersion of *per-family* critical points plus the across-family
    spread that a pooled fit would hide. ``selection_freq`` counts
    per-record model selections (denominator n_boot * n records drawn,
    i.e. n_boot * n).
    """
    rng = np.random.default_rng(seed)
    n = len(records)
    sel_counts: dict[str, int] = {}
    fam_lc: dict[int, list[float]] = {}
    spread_samples: list[float] = []
    n_fits = 0
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        lcs = []
        for i in idx:
            lam_i = np.atleast_1d(records[i][lam_key])
            y_i = np.atleast_1d(records[i][y_key])
            res = fit_all_models(lam_i, y_i,
                                 seed=int(rng.integers(2**31)))
            n_fits += 1
            sel_counts[res["selected"]] = \
                sel_counts.get(res["selected"], 0) + 1
            lc_i = lc_from_fit(res)
            if lc_i is not None:
                lcs.append(lc_i)
                fam_lc.setdefault(int(i), []).append(lc_i)
        uniq = sorted(set(lcs))
        if len(uniq) >= 2:
            spread_samples.append(uniq[-1] - uniq[0])
    fam_summ = {}
    for i, lcs in fam_lc.items():
        a = np.asarray(lcs)
        fam_summ[str(i)] = {
            "median": float(np.median(a)),
            "ci95": [float(q) for q in np.quantile(a, [0.025, 0.975])],
            "n": len(lcs)}
    return {
        "selection_freq": {m: c / max(n_fits, 1)
                           for m, c in sel_counts.items()},
        "selection_counts": sel_counts,
        "n_boot": n_boot,
        "n_records": n,
        "n_fits": n_fits,
        "lc_per_family": fam_summ,
        "lc_cross_family_spread": (
            {"median": float(np.median(spread_samples)),
             "ci95": [float(q) for q in
                      np.quantile(spread_samples, [0.025, 0.975])],
             "n": len(spread_samples)}
            if spread_samples else None),
        "note": ("per-record fits only; no pooled curve. "
                 "lc_cross_family_spread measures the heterogeneity a "
                 "pooled fit would hide."),
    }
