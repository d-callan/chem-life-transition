"""Jacobian spectral diagnostics at certified fixed points.

For CSTR, the linearization is  J = S dv/dx - D*I  (net.jacobian).
The leading nonzero eigenvalue's real part measures local relaxation
rate; critical slowing down appears as Re(lambda_max) -> 0 near a
bifurcation. Serial-dilution stability is instead Floquet (multipliers
of the cycle map) — see ode.sd_fixed_point_multipliers.
"""

from __future__ import annotations

import numpy as np

from ..networks.core import ReactionNetwork


def jacobian_eigs(net: ReactionNetwork, x: np.ndarray, dilution: float,
                  eps: float | None = None) -> np.ndarray:
    """Eigenvalues of the CSTR Jacobian at x, sorted by real part desc."""
    J = net.jacobian(np.maximum(np.asarray(x, float), 0.0), dilution, eps)
    w = np.linalg.eigvals(J)
    return w[np.argsort(-w.real)]


def leading_rate(net: ReactionNetwork, x: np.ndarray, dilution: float,
                 eps: float | None = None) -> float:
    """Re(lambda) of the leading eigenvalue (<0 for stable, ->0 critical)."""
    return float(jacobian_eigs(net, x, dilution, eps)[0].real)


def relaxation_rate_series(records: list[dict], lam_key: str = "lambda",
                           x_key: str = "x_eq") -> np.ndarray:
    """Leading relaxation rate vs control parameter, from run records.

    Each record needs net (or net handle resolved by caller), lambda value,
    and a certified x_eq. Returns array of (lam, re_lam1, n_neg_eigs).
    """
    out = []
    for rec in records:
        net = rec["net"]
        x = np.asarray(rec[x_key], float)
        w = jacobian_eigs(net, x, rec["dilution"], rec.get("eps"))
        out.append((float(rec[lam_key]), float(w[0].real),
                    int(np.sum(w.real > -1e-9))))
    return np.array(out)


def fit_slowing_down(lam: np.ndarray, re1: np.ndarray,
                     lam_c: float) -> dict:
    """Fit |Re lam1| ~ A * |lam - lam_c|^nu on the stable branch.

    Mean-field expectation nu = 1. Returns params + fit diagnostics;
    caller decides whether the window is justified.
    """
    lam = np.asarray(lam, float)
    re1 = np.asarray(re1, float)
    mask = (lam > lam_c) & (re1 < 0)
    if mask.sum() < 3:
        return {"status": "insufficient_stable_points", "n": int(mask.sum())}
    x = np.log(np.abs(lam[mask] - lam_c))
    y = np.log(np.abs(re1[mask]))
    A = np.polyfit(x, y, 1)
    pred = np.polyval(A, x)
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return {
        "status": "ok",
        "nu": float(A[0]),
        "log_A": float(A[1]),
        "r2": 1.0 - ss_res / max(ss_tot, 1e-300),
        "n": int(mask.sum()),
        "theory_nu_meanfield": 1.0,
    }
