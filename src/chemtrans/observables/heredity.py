"""Compositional heredity: multistability and serial transfer.

Two distinct phenomena, explicitly not conflated:

- MULTISTABILITY (Gate B prerequisite): the deterministic system has >1
  certified attractor under identical parameters, distinguished by
  composition (mass-fraction vector over non-food species). Detected by
  integrating from varied initial conditions and clustering certified
  x_eq values.

- COMPOSITIONAL TRANSMISSION (E5): after serial-dilution shocks the
  post-dilution composition returns to the SAME attractor's signature —
  measured as similarity between the pre-shock attractor signature and
  the composition at the end of the next cycle. This is compositional
  heredity, not sequence heredity (the model has no copying mechanism).

"Survival of original molecules" is never a heredity measure: counts
decay by dilution; what must persist is the relative composition.
"""

from __future__ import annotations

import numpy as np

from ..dynamics import ode as ode_mod
from ..networks.core import ReactionNetwork


def attractor_signature(net: ReactionNetwork, x: np.ndarray) -> np.ndarray:
    """Normalized non-food monomer-mass-fraction vector.

    Composition is measured in MONOMER MASS units, not molecule counts,
    consistent with the mass_raf order parameter.
    """
    if net.comp is None:
        w = np.where(net.is_food, 0.0, 1.0)
    else:
        w = np.where(net.is_food, 0.0, net.comp.sum(axis=0))
    v = w * np.maximum(np.asarray(x, float), 0.0)
    s = v.sum()
    return v / s if s > 0 else v


def compositional_distance(sig_a: np.ndarray, sig_b: np.ndarray) -> float:
    """L1 distance between signatures (bounded [0,2] -> normalized [0,1])."""
    return float(0.5 * np.abs(np.asarray(sig_a) - np.asarray(sig_b)).sum())


def cluster_attractors(
    net: ReactionNetwork,
    equilibria: list[np.ndarray],
    tol: float = 0.05,
    min_mass: float = 0.0,
) -> dict:
    """Greedy clustering of certified equilibria by signature distance.

    `min_mass` filters out equilibria whose total non-food mass is below
    the threshold — without it, near-dead endpoints (~1e-35 residue)
    produce normalized "ghost" signatures that mimic organized
    attractors. Returns {n_attractors, labels, centers, members,
    kept_idx}; labels are indexed into the *kept* equilibria.
    """
    if min_mass > 0:
        w = np.where(net.is_food, 0.0, 1.0)
        if net.comp is not None:
            w = w * net.comp.sum(axis=0)
        masses = [float(w @ np.maximum(x, 0.0)) for x in equilibria]
        kept = [x for x, m in zip(equilibria, masses) if m >= min_mass]
    else:
        kept = list(equilibria)
    sigs = [attractor_signature(net, x) for x in kept]
    labels = np.full(len(sigs), -1, dtype=int)
    centers: list[np.ndarray] = []
    for i, s in enumerate(sigs):
        for c, ctr in enumerate(centers):
            if compositional_distance(s, ctr) <= tol:
                labels[i] = c
                break
        else:
            labels[i] = len(centers)
            centers.append(s.copy())
    return {
        "n_attractors": len(centers),
        "labels": labels,
        "centers": centers,
        "members": {c: int((labels == c).sum()) for c in range(len(centers))},
        "n_input": len(equilibria),
        "n_kept": len(kept),
    }


def multistability_scan(
    net: ReactionNetwork,
    dilution: float,
    ic_scale: float = 1e-3,
    n_ic: int = 12,
    seed: int = 0,
    tol: float = 0.05,
    min_mass: float = 0.0,
    **integ_kw,
) -> dict:
    """Integrate from varied initial conditions; cluster attractors.

    ICs: food-only, food+small-random-polymer, and perturbed food-only.
    Returns clustering plus per-IC status so inconclusive runs are
    visible rather than silently dropped.
    """
    rng = np.random.default_rng(seed)
    nonfood = np.where(~net.is_food)[0]
    ics = [net.food_inlet.copy()]
    for _ in range(n_ic - 1):
        x = net.food_inlet.copy()
        pick = rng.choice(nonfood, size=min(8, len(nonfood)),
                          replace=False)
        x[pick] = ic_scale * rng.random(len(pick))
        ics.append(x)

    eqs, labels_ok, status = [], [], []
    for k, x0 in enumerate(ics):
        r = ode_mod.integrate_cstr(net, x0=x0, dilution=dilution,
                                   **integ_kw)
        status.append(r.trajectory_type)
        if r.converged:
            eqs.append(r.x_eq)
            labels_ok.append(k)
    cl = cluster_attractors(net, eqs, tol, min_mass) if eqs else {
        "n_attractors": 0, "labels": np.array([]), "centers": [],
        "members": {}, "n_input": 0, "n_kept": 0}
    return {
        "n_attractors": cl["n_attractors"],
        "n_certified": len(eqs),
        "n_kept": cl["n_kept"],
        "n_ics": len(ics),
        "ic_status": status,
        "converged_ic": labels_ok,
        "cluster_labels": cl["labels"].tolist(),
        "center_distances": _pairwise(cl["centers"]),
        "status": ("ok" if len(eqs) >= 2 else
                   "insufficient_certified_equilibria"),
    }


def _pairwise(centers: list[np.ndarray]) -> list[list[float]]:
    return [[compositional_distance(a, b) for b in centers]
            for a in centers]


def serial_transfer_test(
    net: ReactionNetwork,
    x_organized: np.ndarray,
    tau: float,
    d_factor: float,
    n_cycles: int = 8,
    eps: float | None = None,
) -> dict:
    """Compositional transmission through serial dilution.

    Start at an organized state; iterate the SD map. Transmission is
    assessed by the signature distance between the pre-shock attractor
    and the state at the END of each subsequent batch (post-recovery),
    NOT by raw molecule counts. Returns per-cycle recovery distances and
    a pass/fail against a caller-set threshold upstream.
    """
    if eps is None:
        eps = net.eps
    sig_ref = attractor_signature(net, x_organized)
    x = np.asarray(x_organized, float).copy()

    w = np.where(net.is_food, 0.0, 1.0)
    if net.comp is not None:
        w = w * net.comp.sum(axis=0)

    dist_post_dilution = []
    dist_post_recovery = []
    mass_post_recovery = []
    for _ in range(n_cycles):
        sol = ode_mod.sd_map(net, x, tau, d_factor, 1, eps)
        batch_end = sol[1][0]                    # end of batch, pre-dilution
        x_next = ode_mod.sd_step(net, batch_end, d_factor)
        dist_post_dilution.append(
            compositional_distance(sig_ref,
                                   attractor_signature(net, x_next)))
        dist_post_recovery.append(
            compositional_distance(sig_ref,
                                   attractor_signature(net, batch_end)))
        mass_post_recovery.append(float(w @ np.maximum(batch_end, 0.0)))
        x = x_next
    return {
        "dist_post_dilution": dist_post_dilution,
        "dist_post_recovery": dist_post_recovery,
        "mass_post_recovery": mass_post_recovery,
        "n_cycles": n_cycles,
        "d_factor": d_factor,
        "tau": tau,
    }
