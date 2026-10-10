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
import scipy.sparse as sp_sparse
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


def _krylov_equilibrium(net: ReactionNetwork, x0: np.ndarray,
                        dilution: float, eps: float,
                        conv_atol: float, conv_rtol: float,
                        max_restarts: int = 2,
                        ) -> tuple[np.ndarray, np.ndarray, int, bool, str]:
    """Equilibrium seek for sparse networks: trust-region Newton with
    Krylov inner solves (least_squares trf/lsmr) — needs only sparse
    Jacobian products, never factors the Jacobian.

    Sparse LU is unusable at L>=12: catalytic hubs make LU fill-in
    ~150x (measured: 44 s/factorization, 45M nnz at L=12 vs 0.4 s at
    L=10). lsmr solves the trust-region subproblem by matrix-vector
    products (~ms each).

    Returns (x, scaled_residual, n_rhs_evals, converged, message).
    No trajectory is produced — oscillation/runaway classification is
    unavailable on this path; status is STATIONARY or UNRESOLVED.
    """
    from scipy.optimize import least_squares

    nfev = 0

    def rhs(y):
        return rhs_cstr(net, y, dilution, eps)

    def jac(y):
        return net.jacobian(y, dilution, eps)

    def resid_ok(x):
        # a root in the negative orthant is a root of the polynomial
        # dynamics but not a physical equilibrium — reject it
        r = _scaled_residual(rhs(x), x, conv_atol, conv_rtol)
        return r, bool(np.all(r <= 1.0) and np.all(np.isfinite(x))
                       and x.min() >= -conv_atol)

    x = np.maximum(np.asarray(x0, float), 0.0)
    best_x, best_r = x, resid_ok(x)[0]
    rng = np.random.default_rng(0)

    # Phase 1: pseudo-transient continuation. Follow the physical
    # transient with implicit Euler steps (x_{k+1} = x_k + dt*f) solved
    # by damped Newton with an lsmr inner solve on (I/dt - J) — one
    # Jacobian build per Newton iteration, matvecs only after that.
    # A step is accepted when the *step equation* is solved loosely
    # (||F(y)|| <= 0.2*||f(x)||); errors are damped by continuation
    # itself. dt grows geometrically; converges only to states
    # reachable from x0 — the same basin a trajectory would find.
    # Runs BEFORE the direct solve: Newton polish is for roots, not
    # for finding basins.
    from scipy.sparse.linalg import lsmr

    def ptc(x_start, dt0=0.1, steps=200, budget=300):
        n = len(x_start)
        I = sp_sparse.eye(n, format='csc')
        x = np.maximum(x_start, 0.0)
        dt, spent = dt0, 0
        for _ in range(steps):
            _, ok = resid_ok(x)
            if ok:
                return x, spent
            f0 = np.linalg.norm(rhs(x), np.inf)
            step_ok = False
            for _retry in range(8):
                y = x.copy()
                J = None
                for it in range(10):
                    F = (y - x) / dt - rhs(y)
                    if np.linalg.norm(F, np.inf) <= 0.2 * f0 + 1e-10:
                        break
                    if J is None or it % 3 == 0:
                        spent += 1
                        J = jac(y)
                    d = lsmr(I / dt - J, -F, atol=1e-6, btol=1e-6,
                             maxiter=200)[0]
                    y = np.maximum(y + d, 0.0)
                    if not np.all(np.isfinite(y)):
                        break
                if np.all(np.isfinite(y)) and \
                        np.linalg.norm((y - x) / dt - rhs(y), np.inf) \
                        <= 0.25 * f0 + 1e-10:
                    x = y
                    step_ok = True
                    dt = min(dt * 3.0, 1e12)
                    break
                dt /= 5.0
                if dt < 1e-10 or spent >= budget:
                    break
            if not step_ok or spent >= budget:
                break
        return x, spent

    x_ptc, spent_ptc = ptc(x)
    r, ok = resid_ok(x_ptc)
    nfev += spent_ptc
    if np.nanmax(r) < np.nanmax(best_r):
        best_x, best_r = x_ptc, r
    if ok:
        x = x_ptc
    else:
        x = np.maximum(best_x, 0.0)

    # Phase 2: Newton polish (trf/lsmr) with jittered restarts.
    for attempt in range(1 + max_restarts):
        sol = least_squares(rhs, x, jac=jac, method="trf",
                            tr_solver="lsmr",
                            xtol=1e-11, ftol=1e-11, gtol=1e-11,
                            max_nfev=300)
        nfev += int(sol.nfev)
        # no clip on the polish result: resid_ok rejects negative-
        # orthant roots rather than projecting them into existence
        xp = sol.x
        r, ok = resid_ok(xp)
        if np.nanmax(r) < np.nanmax(best_r):
            best_x, best_r = xp, r
        if ok:
            return xp, r, nfev, True, \
                f"PTC + trf/lsmr polish (attempt {attempt + 1})"
        # perturb restart: sparse multiplicative jitter escapes bad
        # basins without destroying the current scale
        x = np.maximum(sol.x * (1 + 0.1 * rng.standard_normal(len(x))),
                       0.0) if np.all(np.isfinite(sol.x)) \
            else np.maximum(x0, 0.0)
    return best_x, best_r, nfev, False, \
        (f"no equilibrium found (PTC + direct; "
         f"resid_max={np.nanmax(best_r):.2e})")


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
    neg_tol: float = 1e-6,
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

    if net.use_sparse and n > 3000:
        # Sparse-LU implicit steppers are infeasible at this scale
        # (catalytic hubs cause ~150x LU fill-in; measured 44 s/LU at
        # L=12). Trust-region Krylov equilibrium seek instead; no
        # trajectory => no oscillation classification.
        xk, resid, nfev_k, ok, msg = _krylov_equilibrium(
            net, x, dilution, eps, conv_atol, conv_rtol)
        ttype = "STATIONARY" if ok else "UNRESOLVED"
        # steady-state material balance: rhs=0 implies
        # D(feed - M) + 0 = 0 is embedded in the residual; report the
        # letter-mass residual separately as the balance error
        if net.comp is not None:
            dm = net.comp.sum(axis=0) @ rhs_cstr(
                net, xk, dilution, eps)
            feed = max(abs(net.feed_letters()), 1e-12)
            mat_err = float(np.max(np.abs(dm)) / feed)
        else:
            mat_err = np.nan
        return ODEResult(np.array([0.0]), xk[None, :], ok,
                         f"krylov equilibrium seek: {msg}",
                         nfev_k, ttype, resid,
                         np.nan, mat_err,
                         xk if ok else None, eps, dilution)

    def rhs(t, y):
        # vector field evaluated at the TRUE state — never clipped
        # (Invariant 4). The rate law is polynomial in x, so it is
        # well-defined for negative y; sustained negativity is a real
        # solver excursion and is flagged below, not hidden.
        return rhs_cstr(net, y, dilution, eps)

    def jac(t, y):
        return net.jacobian(y, dilution, eps)

    ts: list[np.ndarray] = []
    xs: list[np.ndarray] = []
    nfev = 0
    t = 0.0
    chunk = t_max / n_chunks
    converged = False
    message = "max horizon reached"
    n_eval = max(3, n_points // n_chunks)
    # LSODA wraps ODEPACK (dense Jacobians only). In sparse mode use
    # BDF, which factors the sparse Jacobian by sparse LU — feasible
    # at n <= ~3000 (L=10: ~0.4 s/LU); larger nets take the
    # equilibrium-seek path above.
    method = "BDF" if net.use_sparse else "LSODA"

    for _ in range(n_chunks):
        t_eval = np.linspace(t, t + chunk, n_eval)
        sol = solve_ivp(rhs, (t, t + chunk), x, method=method, jac=jac,
                        t_eval=t_eval, rtol=1e-8, atol=1e-11)
        nfev += sol.nfev
        if not sol.success:
            # strict-JSON-safe: zeros/None, never NaN — the caller must
            # be able to serialize this record as a failed run
            return ODEResult(np.array([t]), x[None, :], False,
                             f"solver failed: {sol.message}", nfev,
                             "UNRESOLVED", np.zeros(n), None,
                             None, None, eps, dilution)
        ts.append(sol.t)
        xs.append(sol.y.T)
        # raw endpoint — no positivity clipping (Invariant 4): a
        # significantly negative state is a real unphysical excursion
        # and must surface, not be silently repaired to zero
        x_new = sol.y[:, -1]
        if x_new.min() < -neg_tol:
            x = x_new
            t = float(sol.t[-1])
            break

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
    # trajectory-flagged 'converged' point must still pass polish — if
    # polish was requested and failed, the window-flagged endpoint is NOT
    # certified (revoke convergence rather than emit an unpolished root).
    polished = False
    polish_note = ""
    if polish:
        if net.use_sparse:
            # trust-region Newton with Krylov inner solves — uses the
            # sparse Jacobian by products, never factors it
            from scipy.optimize import least_squares
            sol_ls = least_squares(
                lambda y: rhs_cstr(net, y, dilution, eps),
                x,
                jac=lambda y: net.jacobian(y, dilution, eps),
                method="trf", tr_solver="lsmr",
                xtol=1e-11, ftol=1e-11, gtol=1e-11, max_nfev=300)
            xp = sol_ls.x
            polish_ok = bool(sol_ls.success)   # lsmr can exit unconverged
        else:
            sol = scipy_root(
                lambda y: rhs_cstr(net, y, dilution, eps),
                x, jac=lambda y: net.jacobian(y, dilution, eps),
                method="hybr")
            xp, polish_ok = sol.x, bool(sol.success)
        pr = _scaled_residual(rhs_cstr(net, xp, dilution, eps),
                              xp, conv_atol, conv_rtol)
        # a polished root with significantly negative components is a
        # root of the polynomial dynamics, not a physical equilibrium —
        # reject rather than clip it into existence
        polished = bool(polish_ok and np.all(pr <= 1.0)
                        and np.all(np.isfinite(xp))
                        and xp.min() >= -neg_tol)
        if polished:
            x, resid_end, converged = xp, pr, True
        else:
            polish_note = ("trajectory passed window test but root "
                           "polish failed" if converged else
                           "root polish failed; trajectory endpoint "
                           "not certified")
            converged = False

    x_min = float(X.min()) if X.size else float(x.min())
    if converged and x.min() >= -neg_tol:
        x_eq = x
        ttype = "STATIONARY"
        message = ("converged (componentwise residual + window"
                   + (", root-polished" if polished else "") + ")")
        if x_min < -neg_tol:
            # certified root, but the trajectory crossed unphysical
            # (negative) territory to reach it — surface the excursion
            # in message + diagnostics rather than hiding it
            message += (f" [transient negative excursion "
                        f"x_min={x_min:.2e}]")
    elif x.min() < -neg_tol:
        ttype = "UNRESOLVED"
        message = (f"negative concentration excursion "
                   f"(x_min={x_min:.2e} <= -neg_tol={neg_tol:g}): "
                   f"trajectory or endpoint left the physical "
                   f"orthant; state not certified")
    elif (not np.all(np.isfinite(resid_end))
          or not np.all(np.isfinite(X))
          or np.max(np.abs(x)) > 1e10):
        # runaway is checked FIRST: a diverging trajectory has large
        # last-window amplitude by construction, so ordering it after
        # the oscillation check mislabels divergence as oscillation
        ttype = "RUNAWAY"
        message = "unbounded growth / nonfinite residual"
    else:
        # positive periodicity check: sustained amplitude AND a
        # near-repeated range between window halves AND low correlation
        # with time (a monotone slow drift fails both halves)
        w = X[-n_eval:]
        p2p = w.max(axis=0) - w.min(axis=0)
        span = np.abs(w.mean(axis=0)) + conv_atol
        amp_rel = float(np.max(p2p / span))
        half = len(w) // 2
        w1, w2 = w[:half], w[half:2 * half]
        mean_diff = float(np.max(
            np.abs(w1.mean(axis=0) - w2.mean(axis=0)) / span))
        t_w = t_all[-n_eval:]
        t_c = t_w - t_w.mean()
        drift_corr = 0.0
        for j in range(net.n_species):
            if p2p[j] / span[j] <= osc_amp_tol:
                continue
            sj = np.std(w[:, j])
            if sj <= 0:
                continue
            drift_corr = max(drift_corr, abs(float(
                np.corrcoef(t_c, w[:, j])[0, 1])))
        if amp_rel > osc_amp_tol and mean_diff < 0.05 \
                and drift_corr < 0.8:
            ttype = "OSCILLATORY"
            message = "nonconverged with periodic oscillation"
        else:
            ttype = "UNRESOLVED"
            if polish_note:
                message = polish_note
            else:
                message = ("no componentwise convergence within t_max "
                           "(slow drift)")

    return ODEResult(t_all, X, True, message, nfev, ttype, resid_end,
                     float(wchange) if 'wchange' in dir() else np.nan,
                     mat_err, x_eq, eps, dilution,
                     diagnostics={
                         "x_min": x_min,
                         "neg_tol": neg_tol,
                         "neg_excursion": bool(x_min < -neg_tol)})


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
    neg_tol: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray, bool, str]:
    """Iterate the SD map: integrate closed batch for tau, then dilute.

    Returns (starts, ends, ok, message): starts[i]/ends[i] are the
    post-dilution / end-of-batch states of cycle i. On solver failure,
    returns the cycles completed so far with ok=False — callers must
    emit UNRESOLVED records rather than crash (mirrors the CSTR
    solver-failure contract). Batch endpoints significantly negative
    (x < -neg_tol) are real unphysical excursions — reported as
    failure, never clipped to zero (Invariant 4)."""
    if eps is None:
        eps = net.eps
    x = np.asarray(x0, dtype=float).copy()
    starts = np.empty((n_cycles, net.n_species))
    ends = np.empty((n_cycles, net.n_species))
    method = "Radau" if net.use_sparse else "LSODA"
    for c in range(n_cycles):
        starts[c] = x
        sol = solve_ivp(
            lambda t, y: net.stoichiometry() @ net.net_flux(y, eps),
            (0.0, tau), np.maximum(x, 0.0),
            method=method,
            jac=lambda t, y: net.jacobian(y, 0.0, eps),
            rtol=1e-8, atol=1e-11,
        )
        if not sol.success:
            return (starts[:c], ends[:c], False,
                    (f"SD batch integration failed at cycle {c}: "
                     f"{sol.message}"))
        x_end = sol.y[:, -1]
        if x_end.min() < -neg_tol:
            return (starts[:c], ends[:c], False,
                    (f"SD batch reached unphysical state at cycle {c} "
                     f"(x_min={float(x_end.min()):.2e})"))
        ends[c] = x_end
        x = sd_step(net, ends[c], d_factor)
    return starts, ends, True, "ok"


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
    failed = {"flag": False}

    def one_cycle(x):
        sol = solve_ivp(
            lambda t, y: net.stoichiometry() @ net.net_flux(y, eps),
            (0.0, tau), np.maximum(x, 0.0), method="LSODA",
            rtol=1e-9, atol=1e-12)
        if not sol.success or sol.y[:, -1].min() < -1e-6:
            failed["flag"] = True
            return sd_step(net, np.maximum(x, 0.0), d_factor)
        return sd_step(net, sol.y[:, -1], d_factor)

    J = np.empty((n, n))
    f0 = one_cycle(x_star)
    for j in range(n):
        if failed["flag"]:
            return None
        e = np.zeros(n)
        e[j] = h * max(1.0, abs(x_star[j]))
        J[:, j] = (one_cycle(x_star + e) - f0) / e[j]
    if failed["flag"]:
        return None
    return np.linalg.eigvals(J)
