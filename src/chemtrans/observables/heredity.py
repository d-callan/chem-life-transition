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


def compositional_distance(sig_a: np.ndarray, sig_b: np.ndarray
                           ) -> float | None:
    """L1 distance between signatures (bounded [0,2] -> normalized [0,1]).

    Extinction is categorical, not compositional: a zero signature
    (no nonfood mass) has no composition to compare. Two zero
    signatures are the same categorical state (distance 0.0); a zero
    vs nonzero signature returns None — an extinct endpoint is not
    "half as far" from a lineage as an orthogonal rival, it is a
    qualitatively different outcome. Callers must handle None.
    """
    a = np.asarray(sig_a, float)
    b = np.asarray(sig_b, float)
    za, zb = a.sum() <= 0, b.sum() <= 0
    if za and zb:
        return 0.0
    if za or zb:
        return None
    return float(0.5 * np.abs(a - b).sum())


def cluster_attractors(
    net: ReactionNetwork,
    equilibria: list[np.ndarray],
    tol: float = 0.05,
    min_mass: float = 0.0,
    mass_ratio_tol: float | None = None,
) -> dict:
    """Greedy clustering of certified equilibria by signature distance.

    `min_mass` filters out equilibria whose total non-food mass is below
    the threshold — without it, near-dead endpoints (~1e-35 residue)
    produce normalized "ghost" signatures that mimic organized
    attractors. Returns {n_attractors, labels, centers, members,
    kept_idx}; labels are indexed into the *kept* equilibria.

    `mass_ratio_tol`: when set, two equilibria join a cluster only if
    BOTH compositional distance <= tol AND their nonfood masses are
    within a `mass_ratio_tol` factor of the cluster center's mass.
    Without it, composition-normalized signatures make amplitude
    bistability (high/low branches of the same composition) invisible.
    """
    w = np.where(net.is_food, 0.0, 1.0)
    if net.comp is not None:
        w = w * net.comp.sum(axis=0)
    masses = [float(w @ np.maximum(x, 0.0)) for x in equilibria]
    if min_mass > 0:
        kept_idx = [i for i, m in enumerate(masses) if m >= min_mass]
        kept = [equilibria[i] for i in kept_idx]
        masses = [masses[i] for i in kept_idx]
    else:
        kept = list(equilibria)
    sigs = [attractor_signature(net, x) for x in kept]
    n = len(sigs)
    # Order-invariant clustering: connected components of the
    # within-tolerance graph (single linkage). The previous greedy
    # first-seen-center scheme depended on IC ordering — permuting the
    # input could merge or split clusters, flipping O_hered on
    # identical data. An edge requires composition within tol AND
    # (when mass_ratio_tol is set) mass ratio within tolerance, so
    # amplitude bistability still separates. Single linkage can merge
    # through chains — that direction is conservative for
    # multistability claims (undercounts rather than fabricates).
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for i in range(n):
        for j in range(i + 1, n):
            d_ij = compositional_distance(sigs[i], sigs[j])
            ok = d_ij is not None and d_ij <= tol
            if ok and mass_ratio_tol is not None:
                ok = (abs(np.log(max(masses[i], 1e-300)
                                 / max(masses[j], 1e-300)))
                      <= np.log(mass_ratio_tol))
            if ok:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[rj] = ri
    comps: dict[int, list[int]] = {}
    for i in range(n):
        comps.setdefault(find(i), []).append(i)
    ordered = sorted(comps.values(), key=lambda mem: mem[0])
    labels = np.full(n, -1, dtype=int)
    centers: list[np.ndarray] = []
    center_mass: list[float] = []
    for c, members in enumerate(ordered):
        for i in members:
            labels[i] = c
        # medoid signature: the member minimizing its worst distance
        # to the rest of the component — representative, not first-seen
        def _worst(i, members=members):
            return max(
                (d if (d := compositional_distance(sigs[i], sigs[j]))
                 is not None else np.inf)
                for j in members)
        med = min(members, key=_worst)
        centers.append(sigs[med].copy())
        center_mass.append(float(np.mean([masses[i] for i in members])))
    return {
        "n_attractors": len(centers),
        "labels": labels,
        "centers": centers,
        "center_masses": center_mass,
        "members": {c: int((labels == c).sum()) for c in range(len(centers))},
        "n_input": len(equilibria),
        "n_kept": len(kept),
        "linkage": "single",
    }


def multistability_scan(
    net: ReactionNetwork,
    dilution: float,
    ic_scale: float = 1e-3,
    n_ic: int = 12,
    seed: int = 0,
    tol: float = 0.05,
    min_mass: float = 0.0,
    n_interior: int = 4,
    **integ_kw,
) -> dict:
    """Integrate from varied initial conditions; cluster attractors.

    ICs: food-only, food+small-random-polymer boundary draws, plus
    `n_interior` interior Dirichlet draws carrying organized-scale
    nonfood mass — boundary-only banks structurally cannot reach
    interior attractors, so n_attractors is reported as a lower bound.
    Returns clustering plus per-IC status so inconclusive runs are
    visible rather than silently dropped.
    """
    rng = np.random.default_rng(seed)
    nonfood = np.where(~net.is_food)[0]
    ics = [net.food_inlet.copy()]
    kinds = ["food"]
    for _ in range(n_ic - 1):
        x = net.food_inlet.copy()
        pick = rng.choice(nonfood, size=min(8, len(nonfood)),
                          replace=False)
        x[pick] = ic_scale * rng.random(len(pick))
        ics.append(x)
        kinds.append("boundary")
    cm_all = np.maximum(net.comp.sum(axis=0), 1e-12)
    feed = max(float(net.feed_letters()), 0.1)
    for _ in range(n_interior):
        x = net.food_inlet.copy()
        k = int(rng.integers(1, min(16, len(nonfood)) + 1))
        pick = rng.choice(nonfood, size=k, replace=False)
        w = rng.dirichlet(np.ones(k))
        x[pick] = w * feed * 10 ** rng.uniform(-1.0, 0.3) / cm_all[pick]
        ics.append(x)
        kinds.append("interior")

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
        "n_attractors_semantics": "lower_bound_on_ic_bank",
        "n_certified": len(eqs),
        "n_kept": cl["n_kept"],
        "n_ics": len(ics),
        "ic_kinds": kinds,
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
    failure = None
    for _ in range(n_cycles):
        _, ends, ok, msg = ode_mod.sd_map(net, x, tau, d_factor, 1, eps)
        if not ok:
            failure = msg
            break
        batch_end = ends[0]                      # end of batch, pre-dilution
        x_next = ode_mod.sd_step(net, batch_end, d_factor)
        dist_post_dilution.append(
            compositional_distance(sig_ref,
                                   attractor_signature(net, x_next)))
        dist_post_recovery.append(
            compositional_distance(sig_ref,
                                   attractor_signature(net, batch_end)))
        mass_post_recovery.append(float(w @ np.maximum(batch_end, 0.0)))
        x = x_next
    out = {
        "dist_post_dilution": dist_post_dilution,
        "dist_post_recovery": dist_post_recovery,
        "mass_post_recovery": mass_post_recovery,
        "n_cycles": n_cycles,
        "n_cycles_completed": len(dist_post_dilution),
        "d_factor": d_factor,
        "tau": tau,
    }
    if failure is not None:
        out["status"] = "UNRESOLVED"
        out["message"] = failure
    return out
