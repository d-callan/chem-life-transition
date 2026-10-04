"""Deterministic mass-action dynamics under a chemostat.

Model: dx/dt = S @ a(x) + feed - dilution * x

- feed: constant influx of food species (energy/matter flux knob)
- dilution: uniform outflow; the "mortality" pressure that organizational
  closure must overcome to persist

This is the workhorse for lambda sweeps: cheap enough to run thousands of
networks per sweep point.
"""

import numpy as np
from scipy.integrate import solve_ivp

from ..networks.core import ReactionNetwork


def simulate(
    net: ReactionNetwork,
    t_span: tuple[float, float],
    feed_rate: float = 1.0,
    dilution_rate: float = 0.1,
    x0: np.ndarray | None = None,
    n_points: int = 2000,
    steady_check: bool = True,
    x_cap: float = 1e5,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate the chemostat ODE.

    Returns (t, X) with X shape (n_points, n_species). Rows after a
    blowup event are NaN — check `steady` via the warning or treat
    np.isnan(X[-1]).any() as a runaway (unbounded autocatalytic growth;
    a genuine dynamical regime, not an integration failure).

    feed_rate: influx per food species (total feed = feed_rate * n_food).
    x0 defaults to zeros + food at feed_rate/dilution_rate (its chemostat
    fixed point), so food species start near their equilibrium.
    x_cap: terminal blowup detector — mass-action autocatalysis has no
    resource saturation, so supercritical networks can diverge; dilution
    cannot catch quadratic growth. We stop rather than grind the solver.
    """
    S = net.stoichiometry()
    n = net.n_species

    feed = np.zeros(n)
    feed[list(net.food)] = feed_rate

    if x0 is None:
        x0 = feed / max(dilution_rate, 1e-12)
    x0 = np.asarray(x0, dtype=float)

    def rhs(t, x):
        x = np.maximum(x, 0.0)
        return S @ net.propensity(x) + feed - dilution_rate * x

    def blowup(t, x):
        return x_cap - x.max()
    blowup.terminal = True
    blowup.direction = -1

    t_eval = np.linspace(t_span[0], t_span[1], n_points)
    sol = solve_ivp(
        rhs, t_span, x0, method="LSODA", t_eval=t_eval,
        rtol=1e-7, atol=1e-10, events=blowup,
    )
    if not sol.success:
        raise RuntimeError(f"ODE integration failed: {sol.message}")

    t_out, X = sol.t, sol.y.T
    if sol.t_events[0].size:  # blew up: pad remainder with NaN
        t_end = sol.t_events[0][-1]
        keep = t_eval <= t_end
        t_out, X = t_eval.copy(), np.full((n_points, n), np.nan)
        m = min(int(keep.sum()), sol.y.shape[1])
        X[:m] = sol.y.T[:m]
        import warnings
        warnings.warn(
            f"runaway: concentration exceeded x_cap={x_cap:g} at t={t_end:.2f}"
        )
    elif steady_check:
        _warn_if_not_steady(net, X[-1], feed, dilution_rate)

    return t_out, X


def _warn_if_not_steady(
    net: ReactionNetwork, x_end: np.ndarray, feed: np.ndarray, dilution: float
) -> None:
    import warnings

    drift = np.linalg.norm(net.stoichiometry() @ net.propensity(x_end)
                           + feed - dilution * x_end)
    scale = max(np.linalg.norm(feed), 1e-12)
    if drift / scale > 1e-4:
        warnings.warn(
            f"final state may not be steady (|dx/dt|/|feed| = {drift/scale:.2e}); "
            "consider a longer t_span"
        )
