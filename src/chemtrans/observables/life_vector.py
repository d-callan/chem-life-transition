"""The life-vector: six independently-defined observables, one record.

Each component is a separate physical question about the same run.
Components return {"value": v, "defined": bool, "kind": str} — a
component is `defined=False` (with a `reason`) when its measurement
context is absent or failed; it is NEVER silently set to zero.

  O_raf      structural maxRAF coverage     kind="structural"
             value = raf_fraction_channels
  O_maint    certified stationary nonfood mass  kind="continuous"
             value = M_poly at x_eq (NOT mass_raf, which saturates)
  O_persist  finite-volume organized dwell      kind="continuous"
             value = mean dwell_organized over SSA reps (V-tagged)
  O_hered    attractor diversity                kind="discrete"
             value = n_organized_attractors - 1 (min_mass-floored)
  O_select   mixture-competition fixation bias  kind="probability"
             value = |P(A wins) - 0.5| per mixture protocol
             CONDITIONAL: only meaningful inside the multistable
             regime (O_hered > 0); never fit as a global curve
  O_assembly joint high-MA + high-abundance     kind="continuous"
             value = mass_signature_abs (absolute monomer-mass
             concentration of a_cut + mass_floor signature species)

Units are deliberately NOT normalized to [0,1] where that would hide
the discrete-vs-continuous distinction (O_hered stays an integer).
"""

from __future__ import annotations

import numpy as np

from ..networks.core import ReactionNetwork


def _nonfood_mass(net: ReactionNetwork, x: np.ndarray) -> float:
    w = np.where(net.is_food, 0.0, 1.0)
    if net.comp is not None:
        w = w * net.comp.sum(axis=0)
    return float(w @ np.maximum(np.asarray(x, float), 0.0))


def _undefined(reason: str) -> dict:
    return {"value": None, "defined": False, "reason": reason}


def _basin_block(ms: dict | None, key: str) -> dict | None:
    """Basin-conditioned observable: per-attractor values + declared
    scalar reductions.

    O_maint/O_assembly's `value` reads whichever basin the reactor's
    policy seed landed in; `multistability.per_attractor_obs` (written
    by measure_run) evaluates the same observable on EVERY certified
    attractor center. Reductions: `max` (best-organized basin) and
    `basin_weighted_mean` (weighted by IC count per cluster — a
    basin-volume proxy from the multi-IC scan).
    """
    pao = (ms or {}).get("per_attractor_obs")
    if not pao:
        return None
    members = ms.get("members") or {}
    w = np.array([float(members.get(str(o["attractor"]),
                                    members.get(o["attractor"], 1)))
                  for o in pao])
    w = w / w.sum() if w.sum() > 0 else np.full(len(pao), 1 / len(pao))
    vals = np.array([o[key] for o in pao], dtype=float)
    return {"per_attractor": {str(o["attractor"]): float(v)
                              for o, v in zip(pao, vals)},
            "reductions": {"max": float(vals.max()),
                           "basin_weighted_mean": float(w @ vals)},
            "weight_rule": "ic_count_per_cluster"}


def life_vector(rec: dict, net: ReactionNetwork | None = None) -> dict:
    """Extract the six-observable vector from a measure_run record.

    `net` is only required for O_maint/O_assembly fallbacks that need
    the composition matrix; records already carry precomputed values.
    """
    out: dict = {}
    ms = rec.get("multistability")

    # ---- O_raf: structural RAF coverage -------------------------
    st = rec.get("structural")
    if st is not None:
        out["O_raf"] = {"value": st["raf_fraction_channels"],
                        "defined": True, "kind": "structural",
                        "n_raf_species": st["n_raf_species"]}
    else:
        out["O_raf"] = _undefined("no structural block")

    # ---- O_maint: stationary nonfood mass -----------------------
    eq = rec.get("equilibrium")
    det = rec.get("det", {})
    if eq is not None and det.get("status") == "STATIONARY":
        if "m_poly" in eq:
            mpoly = eq["m_poly"]
        elif net is not None and "x_eq" in eq:
            mpoly = _nonfood_mass(net, eq["x_eq"])
        else:
            mpoly = None
        out["O_maint"] = ({"value": mpoly, "defined": mpoly is not None,
                           "kind": "continuous",
                           **({"basin": b}
                              if (b := _basin_block(ms, "m_poly"))
                              else {})}
                          if mpoly is not None else
                          _undefined("equilibrium present but no mass"))
    elif det.get("status") in ("OSCILLATORY", "UNRESOLVED", "RUNAWAY"):
        out["O_maint"] = _undefined(f"nonstationary:{det['status']}")
    else:
        out["O_maint"] = _undefined("no certified equilibrium")

    # ---- O_persist: SSA organized dwell -------------------------
    sto = rec.get("stochastic")
    if sto and sto.get("replicates"):
        dwells = [r["dwell_organized"] for r in sto["replicates"]]
        n_surv = sum(1 for r in sto["replicates"]
                     if r.get("t_extinction") is None)
        out["O_persist"] = {
            "value": float(np.mean(dwells)), "defined": True,
            "kind": "continuous", "volume": sto["volume"],
            "n_survive": n_surv, "n_reps": len(dwells)}
    else:
        out["O_persist"] = _undefined("no stochastic block")

    # ---- O_hered: attractor diversity ---------------------------
    if ms is not None:
        out["O_hered"] = {
            "value": max(0, ms["n_attractors"] - 1), "defined": True,
            "kind": "discrete", "n_attractors": ms["n_attractors"],
            "n_attractors_total": ms.get("n_attractors_total"),
            "n_kept": ms.get("n_kept")}
    else:
        out["O_hered"] = _undefined("no multistability block")

    # ---- O_select: competition fixation bias --------------------
    comp = rec.get("competition")
    if comp is not None and comp.get("completed"):
        pa = comp.get("p_a_wins")
        if pa is not None:
            out["O_select"] = {
                "value": abs(pa - 0.5), "defined": True,
                "kind": "probability", "scope": "conditional",
                "condition": "O_hered > 0",
                "n_reps": comp["n_reps"],
                "init_protocol": comp.get("init_protocol")}
        else:
            out["O_select"] = _undefined("competition had no winner "
                                         "labels")
    elif comp is not None:
        out["O_select"] = _undefined(
            f"competition incomplete: {comp.get('reason','?')}")
    else:
        out["O_select"] = _undefined("no competition block")

    # ---- O_assembly: absolute-MA + absolute-mass signature --------
    if eq is not None and "assembly_signature" in eq:
        sig = eq["assembly_signature"]
        if "mass_signature_abs" in sig:
            # new absolute-cut record: value is absolute monomer-mass
            # concentration in signature species
            out["O_assembly"] = {
                "value": sig["mass_signature_abs"],
                "defined": True, "kind": "continuous",
                "a_cut": sig["a_cut"], "mass_floor": sig["mass_floor"],
                "n_signature": sig["n_signature"],
                "pool_ma_max": sig["pool_ma_max"],
                **({"basin": b}
                   if (b := _basin_block(ms, "mass_signature_abs"))
                   else {})}
        else:
            out["O_assembly"] = {
                "value": sig.get("mass_signature"),
                "defined": True, "kind": "continuous",
                "legacy_fraction": True,
                "n_signature": sig.get("n_signature")}
    else:
        out["O_assembly"] = _undefined("no assembly block")

    return out
