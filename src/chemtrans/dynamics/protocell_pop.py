"""Multilevel protocell population dynamics.

A population of N semi-permeable vesicles shares an external nutrient
bath.  Each vesicle j carries concentration vector x_j (same species
ordering as the ReactionNetwork) and volume v_j.  Three physical
couplings connect internal chemistry to compartment fate:

1. Selective permeation.  Species of length 1 (monomers, fuel, waste)
   exchange with the bath at rate P_fast; dimers at P_dimer; polymers
   of length >= 3 cannot cross the membrane (P = 0).  Flux per unit
   internal volume is  J_i = P_i * (x_bath_i - x_i) * A_j / v_j  with
   A_j ~ v_j^{2/3} (spherical vesicle), i.e. J_i ~ v_j^{-1/3}.

2. Osmotic growth / lipid stealing.  Internal osmotic pressure is the
   retained-polymer concentration  Pi_j = sum_{L>=3} x_i.  Vesicles
   above the population mean steal lipids from those below:

       dv_j/dt = mu * v_j * (Pi_j - Pi_bar) - d_leak * v_j

   Growth dilutes internal concentrations: dx_i/dt gets a
   -x_i * (dv_j/dt)/v_j term (product rule on m_i = v_j x_i).

3. Discrete lifecycle events checked each chunk:
   - lysis:     v_j <= v_lysis   -> vesicle bursts, contents lost.
   - division:  v_j >= v_div     -> two daughters at v0 each; species
     molecule counts n_i = round(v_j x_i omega) are partitioned by
     binomial(p=0.5) sampling.  Carrying capacity n_max is enforced by
     flushing the lowest-pressure vesicles.

Deterministic chunks are integrated with adaptive LSODA on the stacked
state [x_flat, v]; events are applied between chunks.  All bookkeeping
is honest: lysis/division are the only places species counts change
non-continuously, and both are recorded in the event log.

Retained polymers leak slowly (p_polymer): real fatty-acid membranes
are imperfect, and without any polymer outflow a vesicle's internal
concentration has no bounded steady state — synthesis would exceed
hydrolysis forever, and even frequent division only splits molecule
counts without lowering concentration.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from ..networks.core import ReactionNetwork


@dataclass
class PopResult:
    """Final population state plus diagnostics."""
    x_cells: np.ndarray          # (n_alive, n_species)
    v_cells: np.ndarray          # (n_alive,)
    n_divisions: int
    n_lysis: int
    n_births: int                # = n_divisions (each adds one cell)
    n_initial: int
    t: np.ndarray                # chunk times recorded
    n_alive_hist: list           # population size over time
    pi_hist: list                # mean osmotic pressure over time
    status: str                  # "ok" | "EXTINCT" | "UNRESOLVED"
    message: str = ""


class ProtocellPopulation:
    """Vectorized population of competing vesicles on one network."""

    def __init__(self, net: ReactionNetwork, n_cells: int = 40,
                 v0: float = 1.0, v_lysis: float | None = None,
                 v_div: float | None = None, mu_growth: float = 1.0,
                 d_leak: float = 0.0, p_fast: float = 1.0,
                 p_dimer: float = 0.2, p_polymer: float = 0.01,
                 omega: float = 1000.0,
                 n_max: int | None = None,
                 chunk: float = 2.0, seed: int = 0,
                 m_init: int = 40):
        self.net = net
        self.v0 = float(v0)
        self.v_lysis = (0.4 * v0) if v_lysis is None else float(v_lysis)
        self.v_div = (2.0 * v0) if v_div is None else float(v_div)
        self.mu = float(mu_growth)
        self.d_leak = float(d_leak)
        self.omega = float(omega)
        self.n_max = int(n_max) if n_max else 4 * int(n_cells)
        self.chunk = float(chunk)
        self.rng = np.random.default_rng(seed)
        if self.chunk <= 0:
            raise ValueError("chunk must be > 0")
        if self.v_lysis >= self.v0 or self.v_div <= self.v0:
            raise ValueError("require v_lysis < v0 < v_div")
        if self.mu < 0 or self.d_leak < 0:
            raise ValueError("mu_growth and d_leak must be >= 0")

        S = net.stoichiometry()
        self._S = np.asarray(S.todense()
                             if hasattr(S, "todense") else S)
        L = np.asarray(net.polymer_lengths, dtype=float)
        # permeability: small molecules fast, dimers slow, polymers
        # strongly retained but leaky (imperfect membrane bounds the
        # internal steady state)
        self._P = np.where(L <= 1, p_fast,
                           np.where(L == 2, p_dimer, p_polymer))
        self._retained = L >= 3
        # external bath concentrations: the network's feed boundary
        self._x_bath = np.maximum(np.asarray(net.food_inlet, float), 0.0)

        # --- encapsulate: each vesicle traps m_init molecules drawn
        # uniformly over the non-food polymer pool (a random thin slice
        # of an abiotic soup), then gets fresh bath concentrations.
        n = net.n_species
        pool = np.where(~net.is_food)[0]
        self.x = np.zeros((n_cells, n))
        self.x[:, :] = self._x_bath[None, :]
        self.v = np.full(n_cells, self.v0)
        p = np.full(len(pool), 1.0 / len(pool))
        for j in range(n_cells):
            counts = self.rng.multinomial(m_init, p)
            self.x[j, pool] += counts / (self.v0 * self.omega)
        self.n_divisions = 0
        self.n_lysis = 0
        self.n_initial = n_cells

    # ------------------------------------------------ dynamics ----
    def _rhs(self, x: np.ndarray, v: np.ndarray):
        """Return (dx/dt, dv/dt) for all live cells."""
        n = x.shape[0]
        dx = np.empty_like(x)
        for j in range(n):
            dx[j] = self._S @ self.net.net_flux(x[j])
        # permeation ~ A/v ~ v^{-1/3}; normalized to rate P at v=v0
        shape = (self.v0 / np.maximum(v, 1e-12)) ** (1.0 / 3.0)
        dx += (self._P[None, :] * (self._x_bath[None, :] - x)
               * shape[:, None])
        # osmotic pressure + competitive growth
        pi = x[:, self._retained].sum(axis=1)
        pi_bar = pi.mean() if n else 0.0
        dv = self.mu * v * (pi - pi_bar) - self.d_leak * v
        # growth dilutes / shrinkage concentrates internal species
        dx -= (dv / np.maximum(v, 1e-12))[:, None] * x
        return dx, dv

    def _rhs_cell(self, j: int, pi_bar: float):
        """Curried (t, y) -> dy/dt for cell j with the population-mean
        pressure frozen (Strang-style operator splitting: the
        mean-field coupling in dv is held constant over a chunk, which
        is second-order-consistent when the chunk is short relative to
        the growth timescale)."""
        P, xb, S, net = self._P, self._x_bath, self._S, self.net
        retained, v0, mu, d_leak = (self._retained, self.v0,
                                    self.mu, self.d_leak)

        def f(t, y):
            x, v = y[:-1], max(y[-1], 1e-12)
            dx = S @ net.net_flux(x)
            shape = (v0 / v) ** (1.0 / 3.0)
            dx = dx + P * (xb - x) * shape
            pi = x[retained].sum()
            dv = mu * v * (pi - pi_bar) - d_leak * v
            dx = dx - (dv / v) * x
            out = np.empty_like(y)
            out[:-1], out[-1] = dx, dv
            return out
        return f

    def _integrate_chunk(self, t_span: tuple[float, float]) -> bool:
        """Advance each cell independently over t_span with LSODA and
        the mean-field coupling frozen at the chunk-start value.
        Returns False when a cell's chunk cannot be certified."""
        pi = self.x[:, self._retained].sum(axis=1)
        pi_bar = float(pi.mean()) if pi.size else 0.0
        xs, vs = [], []
        for j in range(len(self.v)):
            y0 = np.concatenate([self.x[j], [self.v[j]]])
            sol = solve_ivp(self._rhs_cell(j, pi_bar), t_span, y0,
                            method="LSODA", rtol=1e-7, atol=1e-10)
            if not sol.success:
                return False
            y = sol.y[:, -1]
            xs.append(y[:-1])
            vs.append(y[-1])
        self.x = np.asarray(xs)
        self.v = np.asarray(vs)
        # clip only sub-roundoff solver excursions
        self.x[(self.x < 0) & (self.x > -1e-9)] = 0.0
        self.v = np.maximum(self.v, 0.0)
        return True

    # -------------------------------------------------- events ----
    def _events(self):
        # lysis first: below-threshold cells burst
        dead = self.v <= self.v_lysis
        if np.any(dead):
            self.n_lysis += int(dead.sum())
            keep = ~dead
            self.x, self.v = self.x[keep], self.v[keep]
        if len(self.v) == 0:
            return
        # division: above-threshold cells split binomially
        split = self.v >= self.v_div
        if np.any(split):
            idx = np.where(split)[0]
            new_x, new_v = [], []
            for j in idx:
                n_i = np.floor(self.x[j] * self.v[j] * self.omega
                               + 0.5).astype(int)
                n_a = self.rng.binomial(np.maximum(n_i, 0), 0.5)
                xa = n_a / (self.v0 * self.omega)
                xb = (n_i - n_a) / (self.v0 * self.omega)
                # residual (sub-molecule) mass stays with daughter b
                resid = np.maximum(
                    self.x[j] - n_i / (self.v[j] * self.omega), 0.0)
                xb += resid
                self.x[j] = xa
                self.v[j] = self.v0
                new_x.append(xb)
                new_v.append(self.v0)
                self.n_divisions += 1
            self.x = np.vstack([self.x] + [np.asarray(new_x)])
            self.v = np.concatenate([self.v, np.asarray(new_v)])
        # carrying capacity: flush lowest-pressure cells
        if len(self.v) > self.n_max:
            pi = self.x[:, self._retained].sum(axis=1)
            keep = np.argsort(pi)[::-1][: self.n_max]
            n_drop = len(self.v) - self.n_max
            self.n_lysis += n_drop
            self.x, self.v = self.x[keep], self.v[keep]

    # ---------------------------------------------------- driver --
    def run(self, t_max: float = 500.0) -> PopResult:
        t_hist, n_hist, pi_hist = [0.0], [len(self.v)], [
            float(self.x[:, self._retained].sum(axis=1).mean())
            if len(self.v) else 0.0]
        t, status, msg = 0.0, "ok", ""
        while t < t_max - 1e-12:
            dt_chunk = min(self.chunk, t_max - t)
            if not self._integrate_chunk((t, t + dt_chunk)):
                status, msg = "UNRESOLVED", "solver failed on chunk"
                return self._result(t, t_hist, n_hist, pi_hist,
                                    status, msg)
            t += dt_chunk
            if not np.all(np.isfinite(self.x)):
                status, msg = "UNRESOLVED", "non-finite state"
                return self._result(t, t_hist, n_hist, pi_hist,
                                    status, msg)
            self._events()
            if len(self.v) == 0:
                status, msg = "EXTINCT", "population lysed out"
                return self._result(t, t_hist, n_hist, pi_hist,
                                    status, msg)
            t_hist.append(t)
            n_hist.append(len(self.v))
            pi_hist.append(float(
                self.x[:, self._retained].sum(axis=1).mean()))
        return self._result(t, t_hist, n_hist, pi_hist, status, msg)

    def _result(self, t, t_hist, n_hist, pi_hist, status, msg):
        return PopResult(
            x_cells=self.x.copy(), v_cells=self.v.copy(),
            n_divisions=self.n_divisions, n_lysis=self.n_lysis,
            n_births=self.n_divisions, n_initial=self.n_initial,
            t=np.asarray(t_hist), n_alive_hist=n_hist,
            pi_hist=pi_hist, status=status, message=msg)


# ------------------------------------------------- observables ----
def speciation_index(x_cells: np.ndarray, retained_mask: np.ndarray,
                     floor: float = 1e-12) -> float | None:
    """Mean pairwise cosine distance between cells' retained-polymer
    sequence profiles.  0 = identical commons; ->1 = specialized
    non-overlapping lineages.  None when <2 live cells or the state
    is non-finite (UNRESOLVED run)."""
    n = x_cells.shape[0]
    if n < 2 or not np.all(np.isfinite(x_cells)):
        return None
    p = x_cells[:, retained_mask]
    norms = np.linalg.norm(p, axis=1)
    norms = np.maximum(norms, floor)
    sim = (p @ p.T) / np.outer(norms, norms)
    iu = np.triu_indices(n, k=1)
    return float((1.0 - sim[iu]).mean())


def per_cell_sparsity(x_cells: np.ndarray,
                      length_masks: dict[int, np.ndarray],
                      floor: float = 1e-12) -> dict[str, float | None]:
    """Mean effective diversity and Gini-style sparsity per length pool
    across live cells: sparsity = 1 - N_eff/n_pool where
    N_eff = exp(H), H = -sum p ln p."""
    out: dict[str, float | None] = {}
    if x_cells.shape[0] == 0 or not np.all(np.isfinite(x_cells)):
        return out
    for L, mask in length_masks.items():
        sub = x_cells[:, mask]
        mass = sub.sum(axis=1)
        ok = mass > floor
        if not np.any(ok):
            out[f"L{L}_n_eff"] = None
            out[f"L{L}_sparsity"] = None
            continue
        p = np.clip(sub[ok] / mass[ok, None], 0.0, None)
        with np.errstate(divide="ignore", invalid="ignore"):
            H = -np.sum(np.where(p > 0, p * np.log(p), 0.0), axis=1)
        n_eff = np.exp(H)
        n_pool = mask.sum()
        out[f"L{L}_n_eff"] = float(n_eff.mean())
        out[f"L{L}_sparsity"] = float(
            (1.0 - n_eff / n_pool).mean())
    return out


def signature_cells(x_cells: np.ndarray, ma: np.ndarray,
                    a_cut: int, mass_floor: float) -> dict:
    """High-assembly abundance signature across the population.

    A cell 'has signature' when at least one species with molecular
    assembly index >= a_cut reaches concentration >= mass_floor."""
    hi = ma >= a_cut
    if (not np.any(hi) or x_cells.shape[0] == 0
            or not np.all(np.isfinite(x_cells))):
        return {"n_signature_cells": 0, "max_signature_mass": 0.0}
    sig_mass = x_cells[:, hi].max(axis=1)
    return {
        "n_signature_cells": int((sig_mass >= mass_floor).sum()),
        "max_signature_mass": float(sig_mass.max()),
    }
