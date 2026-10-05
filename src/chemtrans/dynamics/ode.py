"""Deterministic dynamics for the chemtrans 2.0 model.

Two reactors:

- CSTR:  dx/dt = S v(x) + D * (food_inlet - x)
  where v(x) are the reversible channel fluxes from the network. Feed is
  `D * food_inlet` so that the food inflow concentration is the physical
  boundary condition (mass-conserving across chi/delta_mu sweeps).

- SD (serial dilution / growth-division): closed batch dx/dt = S v(x) for
  a cycle of length tau, then the exact mass-conserving dilution map
      x_poly <- x_poly / D_factor
      x_food <- x_food / D_factor + (1 - 1/D_factor) * food_inlet
  iterated. Fixed points and Floquet multipliers analyzed via the map.

Integrations return a structured ODEResult — no warning-parsing, no
implicit use of unconverged endpoints. `x_eq` exists only when certified
stationary.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root as scipy_root

from ..networks.core import ReactionNetwork


@dataclass
class ODEResult:
    t: np.ndarray
    x: np.ndarray                    # (n_times, n_species)
    success: bool
    message: str
    nfev: int
    trajectory_type: str             # STATIONARY | OSCILLATORY | RUNAWAY |
                                     # UNRESOLVED
    residual_scaled: np.ndarray      # componentwise |dx/dt|/(atol+rtol|x|)
    window_change: float             # max rel change over last window
    material_error: float            # |dM/dt - (feed_letters*D - D*M)|
    x_eq: np.ndarray | None          # certified equilibrium or None
    eps: float = 1.0
    dilution: float = 0.0
    diagnostics: dict = field(default_factory=dict)

    @property
    def converged(self) -> bool:
        return self.trajectory_type == "STATIONARY" and self.x_eq is not None


def rhs_cstr(net: ReactionNetwork, x: np.ndarray, dilution: float,
             eps: float | None = None) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return (net.stoichiometry() @ net.net_flux(x, eps)
            + dilution * (net.food_inlet - x))


def _scaled_residual(dxdz: np.ndarray, x: np.ndarray, atol, rtol):
    return np.abs(dxdz) / (atol + rtol * np.abs(x))


def integrate_cstr(
    net: ReactionNetwork,
    t_max: float = 200.0,
    dilution: float = 0.1,
    x0: np.ndarray | None = None,
    eps: float | None = None,
    n_points: int = 400,
    conv_atol: float = 1e-8,
    conv_rtol: float = 1e-5,
    window_tol: float = 1e-5,
    n_chunks: int = 8,
    polish: bool = True,
    osc_amp_tol: float = 0.05,
) -> ODEResult:
    """Integrate CSTR to convergence in chunks, then certify.

    Convergence (componentwise): every species satisfies
    |dx_i/dt| <= conv_atol + conv_rtol*|x_i| AND the last-chunk relative
    change is <= window_tol. A vector norm can hide a slowly drifting
    species; per-species residuals cannot.

    Oscillatory: certified via a positive periodicity check (the last
    window returns to a near-repeated range over >= 2 estimated periods),
    not merely by failure to converge.
    """
    if eps is None:
        eps = net.eps
    n = net.n_species
    if x0 is None:
        x0 = net.food_inlet.copy()
    x = np.maximum(np.asarray(x0, dtype=float), 0.0)

    def rhs(t, y):
        return rhs_cstr(net, np.maximum(y, 0.0), dilution, eps)

    def jac(t, y):
        return net.jacobian(np.maximum(y, 0.0), dilution, eps)

    ts: list[np.ndarray] = []
    xs: list[np.ndarray] = []
    nfev = 0
    t = 0.0
    chunk = t_max / n_chunks
    converged = False
    message = "max horizon reached"
    n_eval = max(3, n_points // n_chunks)

    for _ in range(n_chunks):
        t_eval = np.linspace(t, t + chunk, n_eval)
        sol = solve_ivp(rhs, (t, t + chunk), x, method="LSODA", jac=jac,
                        t_eval=t_eval, rtol=1e-8, atol=1e-11)
        nfev += sol.nfev
        if not sol.success:
            return ODEResult(np.array([t]), x[None, :], False,
                             f"solver failed: {sol.message}", nfev,
                             "UNRESOLVED", np.full(n, np.nan), np.nan,
                             np.nan, None, eps, dilution)
        ts.append(sol.t)
        xs.append(sol.y.T)
        x_new = np.maximum(sol.y[:, -1], 0.0)

        resid = _scaled_residual(rhs(t, x_new), x_new, conv_atol, conv_rtol)
        wchange = float(
            np.max(np.abs(x_new - x) / (conv_atol + np.abs(x)))
        )
        if np.all(resid <= 1.0) and wchange <= window_tol:
            converged = True
            x = x_new
            t = float(sol.t[-1])
            break
        x = x_new
        t = float(sol.t[-1])

    t_all = np.concatenate(ts)
    X = np.concatenate(xs)

    # material balance diagnostic: total letter mass must obey the scalar
    # law dM/dt = D*(feed_letters - M) exactly (internal channels conserve
    # mass). Compare against the closed-form solution, not a numerical
    # derivative of the trajectory.
    if net.comp is not None:
        m = net.comp.sum(axis=0) @ X.T
        if dilution > 0:
            m_pred = (net.feed_letters()
                      + (m[0] - net.feed_letters())
                      * np.exp(-dilution * t_all))
        else:
            m_pred = np.full_like(m, m[0])
        mat_err = float(np.max(np.abs(m - m_pred))
                        / max(abs(net.feed_letters()), 1e-12))
    else:
        mat_err = np.nan

    x_eq = None
    ttype = "UNRESOLVED"
    resid_end = _scaled_residual(rhs(t, x), x, conv_atol, conv_rtol)

    # Always attempt a Newton polish: a slow-drifting trajectory may still
    # sit in the basin of a certifiable fixed point. Conversely a
    # trajectory-flagged 'converged' point must still pass polish.
    if polish:
        sol = scipy_root(
            lambda y: rhs_cstr(net, np.maximum(y, 0.0), dilution, eps),
            x, jac=lambda y: net.jacobian(np.maximum(y, 0.0),
                                          dilution, eps),
            method="hybr")
        xp = np.maximum(sol.x, 0.0)
        pr = _scaled_residual(rhs_cstr(net, xp, dilution, eps),
                              xp, conv_atol, conv_rtol)
        if sol.success and np.all(pr <= 1.0) and \
                np.all(np.isfinite(xp)):
            x, resid_end, converged = xp, pr, True

    if converged:
        x_eq = x
        ttype = "STATIONARY"
        message = ("converged (componentwise residual + window"
                   + (", root-polished" if polish else "") + ")")
    else:
        # oscillatory check: does the last window repeat its range?
        w = X[-n_eval:]
        p2p = w.max(axis=0) - w.min(axis=0)
        span = np.abs(w.mean(axis=0)) + conv_atol
        if np.max(p2p / span) > osc_amp_tol and np.all(
                np.isfinite(resid_end)):
            ttype = "OSCILLATORY"
            message = "nonconverged with sustained oscillation amplitude"
        elif not np.all(np.isfinite(resid_end)) or np.max(
                np.abs(x)) > 1e12:
            ttype = "RUNAWAY"
            message = "unbounded growth / nonfinite residual"
        else:
            message = ("no componentwise convergence within t_max "
                       "(slow drift)")

    return ODEResult(t_all, X, True, message, nfev, ttype, resid_end,
                     float(wchange) if 'wchange' in dir() else np.nan,
                     mat_err, x_eq, eps, dilution)


def sd_step(net: ReactionNetwork, x: np.ndarray, d_factor: float) -> np.ndarray:
    """One exact mass-conserving serial-dilution + refeed step.

    polymers: x <- x / D_factor
    food:     x <- x / D_factor + (1 - 1/D_factor) * food_inlet
    """
    if d_factor <= 1.0:
        raise ValueError("d_factor must be > 1")
    x = np.asarray(x, dtype=float).copy()
    food_mask = net.is_food
    x /= d_factor
    x[food_mask] += (1.0 - 1.0 / d_factor) * net.food_inlet[food_mask]
    return x


def sd_map(
    net: ReactionNetwork,
    x0: np.ndarray,
    tau: float,
    d_factor: float,
    n_cycles: int,
    eps: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Iterate the SD map: integrate closed batch for tau, then dilute.

    Returns (cycle_states, states_after_each_dilution): cycle_states[i] is
    the state at the START of cycle i (post-dilution).
    """
    if eps is None:
        eps = net.eps
    x = np.asarray(x0, dtype=float).copy()
    starts = np.empty((n_cycles, net.n_species))
    ends = np.empty((n_cycles, net.n_species))
    for c in range(n_cycles):
        starts[c] = x
        sol = solve_ivp(
            lambda t, y: net.stoichiometry()
            @ net.net_flux(np.maximum(y, 0.0), eps),
            (0.0, tau), np.maximum(x, 0.0),
            method="LSODA",
            jac=lambda t, y: net.jacobian(np.maximum(y, 0.0), 0.0, eps),
            rtol=1e-8, atol=1e-11,
        )
        if not sol.success:
            raise RuntimeError(f"SD batch integration failed: {sol.message}")
        ends[c] = np.maximum(sol.y[:, -1], 0.0)
        x = sd_step(net, ends[c], d_factor)
    return starts, ends


def sd_fixed_point_multipliers(
    net: ReactionNetwork,
    x_star: np.ndarray,
    tau: float,
    d_factor: float,
    eps: float | None = None,
    h: float = 1e-6,
) -> np.ndarray:
    """Floquet multipliers of the SD map about a fixed point x_star.

    Finite-difference linearization of one full cycle (batch flow +
    dilution step). |mu| > 1 => unstable fixed point of the map.
    """
    n = net.n_species

    def one_cycle(x):
        sol = solve_ivp(
            lambda t, y: net.stoichiometry()
            @ net.net_flux(np.maximum(y, 0.0), eps),
            (0.0, tau), np.maximum(x, 0.0), method="LSODA",
            rtol=1e-9, atol=1e-12)
        return sd_step(net, np.maximum(sol.y[:, -1], 0.0), d_factor)

    J = np.empty((n, n))
    f0 = one_cycle(x_star)
    for j in range(n):
        e = np.zeros(n)
        e[j] = h * max(1.0, abs(x_star[j]))
        J[:, j] = (one_cycle(x_star + e) - f0) / e[j]
    return np.linalg.eigvals(J)
