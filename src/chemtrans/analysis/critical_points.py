"""Critical-point / transition-shape inference.

Fits candidate response curves to order-parameter-vs-control data and
selects by corrected AIC. Continuous fits only — no grid-snapped
critical points, no "first index where O jumps" estimators.

Models (functions of control lam):

  null:    O = c0
  ramp:    O = c0 + c1*(lam - lam0)
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
        "ramp": (lambda p, L: p["c0"] + p["c1"] * (L - p["lam0"]),
                 {"c0": None, "c1": None, "lam0": None},
                 [{"c0": 0.1, "c1": s, "lam0": 0.0}
                  for s in (0.5, 2.0)]),
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
    }


def _free_keys(model: str) -> list[str]:
    return list(_model_preds()[model][1].keys())


def _fit_one(model: str, lam: np.ndarray, y: np.ndarray,
             weights: np.ndarray | None = None,
             seed: int = 0) -> FitResult:
    from scipy.optimize import least_squares
    pred_fn, proto, starts = _model_preds()[model]
    keys = _free_keys(model)
    lam_mid = float(np.median(lam))
    rng = np.random.default_rng(seed)
    w = np.ones_like(y) if weights is None else np.asarray(weights, float)

    def unpack(v):
        return dict(zip(keys, v))

    best = None
    for start in starts:
        p0 = []
        for k in keys:
            if k in ("lc", "lf", "lam0"):
                p0.append(lam_mid if start.get(k) is None
                          else start[k])
            elif k == "c0" and start.get(k) is None:
                p0.append(float(np.min(y)))
            else:
                p0.append(start.get(k, 0.1))
        p0 = np.array(p0, float)
        if model == "null":
            theta = np.array([np.average(y, weights=w)])
            resid = (y - theta[0]) * np.sqrt(w)
            rss = float(resid @ resid)
            cand = (theta, rss)
        else:
            def r(v):
                return (pred_fn(unpack(v), lam) - y) * np.sqrt(w)
            try:
                sol = least_squares(r, p0, method="lm", max_nfev=4000)
            except Exception:
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
        "n": int(len(y)),
        "delta_aicc": {m: fits[m].aicc - best.aicc for m in order},
        "_fits": fits,
    }


def network_joint_bootstrap(
    records: list[dict],
    lam_key: str, y_key: str,
    n_boot: int = 400, seed: int = 0,
) -> dict:
    """Bootstrap over NETWORK records, preserving paired observables.

    A record is one network's (possibly multi-lambda) measurement set.
    Resampling whole records preserves the nested trajectory structure
    and the pairing between lambda and O — unlike resampling points.
    Returns bootstrap distribution of selected-model labels and of each
    model's parameters; CIs are only as meaningful as n_records.
    """
    rng = np.random.default_rng(seed)
    n = len(records)
    sel = []
    par_samples: dict[str, list[dict]] = {}
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        lam = np.concatenate([np.atleast_1d(records[i][lam_key])
                              for i in idx])
        y = np.concatenate([np.atleast_1d(records[i][y_key])
                            for i in idx])
        res = fit_all_models(lam, y, seed=int(rng.integers(2**31)))
        sel.append(res["selected"])
        par_samples.setdefault(res["selected"], []).append(res["params"])
    labels, counts = np.unique(sel, return_counts=True)
    return {
        "selection_freq": {l: int(c) for l, c in zip(labels, counts)},
        "n_boot": n_boot,
        "n_records": n,
        "param_samples": par_samples,
    }
