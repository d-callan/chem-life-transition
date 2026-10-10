"""measure_run: one normalized measurement task -> one strict-JSON record.

A "task spec" combines a network generator spec with reactor and
stochastic settings. The record is self-describing: net_hash, all
generator/reactor parameters, status (stationary/oscillatory/runaway/
unresolved/truncated/complete), observables, and diagnostics.

Conventions:
- Order parameter: mass_raf (monomer mass in maxRAF species / total
  non-food monomer mass). phi_* are diagnostics.
- Only STATIONARY results carry equilibrium observables under
  `equilibrium`; nonconverged runs report endpoint/trajectory observables
  under `final_state` with status flags — never conflated.
- Seed policy: integer seeds are derived deterministically from
  seed_key strings so task reordering preserves results.
"""

from __future__ import annotations

import hashlib
import time

import numpy as np

from .dynamics import gillespie, ode
from .networks import benchmarks as bench
from .networks.core import ReactionNetwork
from .networks.generators import polymer_network
from .observables import assembly, closure, heredity, spectral
from .provenance import run_id, sha256_json, strict_jsonable


def derive_seed(seed_key: str) -> int:
    return int.from_bytes(
        hashlib.sha256(seed_key.encode()).digest()[:8], "little")


def build_network(gen_spec: dict) -> tuple[ReactionNetwork, dict]:
    """Dispatch on gen_spec['kind']: 'polymer' | 'transcritical' |
    'schlogl'. Returns (net, analytic-or-empty dict). Unknown keys are
    filtered out by signature so injected keys like 'seed' are safe."""
    import inspect
    kind = gen_spec.get("kind", "polymer")
    fns = {"polymer": polymer_network,
           "transcritical": bench.transcritical_motif,
           "schlogl": bench.schlogl}
    if kind not in fns:
        raise ValueError(f"unknown generator kind {kind!r}")
    fn = fns[kind]
    params = {k: v for k, v in gen_spec.items()
              if k in inspect.signature(fn).parameters}
    if kind == "transcritical":
        return fn(**params)               # returns (net, analytic)
    return fn(**params), {}               # polymer, schlogl


def mass_raf_series(net: ReactionNetwork, n: np.ndarray,
                    volume: float) -> np.ndarray:
    """mass_raf along a recorded SSA count trajectory (rows = times)."""
    core = np.zeros(net.n_species, bool)
    for i in closure.raf_species(net):
        core[i] = True
    comp_mass = net.comp.sum(axis=0)
    nonfood = ~net.is_food
    conc = n.astype(float) / volume
    core_mass = conc[:, core] @ comp_mass[core]
    tot_mass = conc[:, nonfood] @ comp_mass[nonfood]
    return np.where(tot_mass > 0, core_mass / np.maximum(tot_mass, 1e-300),
                    0.0)


def measure_run(task: dict) -> dict:
    """Execute one measurement task; return strict-JSON-able record."""
    t0 = time.time()
    gen_spec = task["generator"]
    seed_key = str(task.get("seed_key", "default"))
    gen_spec = dict(gen_spec)
    gen_spec.setdefault("seed", derive_seed(seed_key + "|gen"))

    net, analytic = build_network(gen_spec)
    reactor = task.get("reactor", {})
    rtype = reactor.get("type", "cstr")
    dilution = float(reactor.get("dilution", 0.1))
    eps = task.get("eps", net.eps)

    rec = {
        "generator": gen_spec,
        "reactor": {**reactor, "type": rtype},
        "eps": eps,
        "seed_key": seed_key,
        "net_hash": net.hash,
        # family = the latent topology: the resolved generator seed plus
        # pool-defining params. Sweeping thresholds/energies (p_cat,
        # delta_mu, kappa, ...) masks the same latent draws and keeps
        # the family; changing the pool (alphabet_size, max_len,
        # ligation) does not.
        "family_id": seed_key + "|" + sha256_json({
            "kind": gen_spec.get("kind", "polymer"),
            "seed": gen_spec["seed"],
            "alphabet_size": gen_spec.get("alphabet_size"),
            "max_len": gen_spec.get("max_len"),
            "ligation": gen_spec.get("ligation", False),
        })[:12],
        "n_species": net.n_species,
        "n_channels": net.n_channels,
        "n_catalyzed_channels": int(sum(bool(c.catalysts)
                                        for c in net.channels)),
    }
    # run identity covers every analysis-affecting task key (assembly,
    # competition, continuation, ...) — not just gen/reactor/eps/
    # stochastic — so a task that differs only in an analysis block can
    # never be silently skipped as "already done" on resume
    task_identity = {k: v for k, v in task.items()
                     if k not in ("seed_key", "name")}
    task_identity["generator"] = gen_spec   # resolved seed included
    rec["run_id"] = run_id(net.hash, task_identity, seed_key)

    # ---- structural block: always emitted (cheap; O_raf context) --
    core = closure.raf_reactions(net)
    core_sp = closure.raf_species(net)
    n_elong = sum(c.mechanism == "elongation" for c in net.channels)
    rec["structural"] = {
        "n_raf_channels": len(core),
        "n_raf_species": len(core_sp),
        "n_catalyzed": int(sum(bool(c.catalysts)
                               for c in net.channels)),
        "n_elong_channels": n_elong,
        "realized_p_cat": (sum(bool(c.catalysts) for c in net.channels
                               if c.mechanism == "elongation")
                           / max(n_elong, 1)),
        "n_kinetic_accessible": len(closure.kinetic_accessible(net)),
        "raf_fraction_channels": len(core) / max(net.n_channels, 1),
        "raf_labels": sorted(net.species[i] for i in core_sp),
    }

    # ---- deterministic ----
    if "x0" in reactor:
        x0_det = np.asarray(reactor["x0"], float)
    elif reactor.get("x0_policy") == "dense":
        x0_det = net.food_inlet.copy()
        x0_det[~net.is_food] = float(reactor.get("x0_level", 0.15))
    elif reactor.get("x0_policy") == "food_only":
        x0_det = net.food_inlet.copy()
    elif str(reactor.get("x0_policy", "")).startswith("lineage"):
        # lineage:<letter> or lineage_a/lineage_b -> seed the species
        # group dominated by that letter (same rule as _ic_bank).
        pol = str(reactor["x0_policy"])
        letter = (pol.split(":", 1)[1] if ":" in pol
                  else pol.rsplit("_", 1)[-1]).upper()
        x0_det = _lineage_seed(net, seed_key, letter,
                               level=float(reactor.get("x0_level", 0.08)))
    else:
        x0_det = None
    if rtype == "structural":
        # E1: structural landscape only — RAF topology without dynamics.
        rec["det"] = {"type": "structural", "status": "NA"}
    elif rtype == "cstr":
        det = ode.integrate_cstr(
            net, t_max=float(reactor.get("t_max", 200.0)),
            dilution=dilution, eps=eps, x0=x0_det,
            n_points=int(reactor.get("n_points", 400)))
        rec["det"] = {
            "type": rtype,
            "status": det.trajectory_type,
            "success": det.success,
            "message": det.message,
            "nfev": det.nfev,
            "residual_scaled_max": (
                float(np.nanmax(det.residual_scaled))
                if np.any(np.isfinite(det.residual_scaled)) else None),
            "window_change": (float(det.window_change)
                              if det.window_change is not None
                              and np.isfinite(det.window_change)
                              else None),
            "material_error": (float(det.material_error)
                               if det.material_error is not None
                               and np.isfinite(det.material_error)
                               else None),
            "x_min": (float(det.diagnostics["x_min"])
                      if det.diagnostics
                      and np.isfinite(det.diagnostics.get(
                          "x_min", np.nan))
                      else None),
            "neg_excursion": (bool(det.diagnostics["neg_excursion"])
                              if det.diagnostics else None),
        }
        if det.converged:
            x = det.x_eq
            obs = closure.phi_decomposition(net, x, eps)
            eigs, spec_conv = spectral.jacobian_eigs(net, x, dilution, eps)
            rec["equilibrium"] = {
                "mass_raf": obs["mass_raf"],
                "phi_gross": obs["phi_gross"],
                "phi_cat": obs["phi_cat"],
                "phi_net": obs["phi_net"],
                "denominators": obs["denominators"],
                "n_raf_channels": obs["n_raf_channels"],
                "leading_eig_re": (float(eigs[0].real)
                                 if spec_conv and len(eigs) else None),
                "spectral_converged": bool(spec_conv),
                "scaled_resid_max": float(np.max(det.residual_scaled)),
                "nonzero_mass_species": int(np.sum(x > 1e-10)),
                "m_poly": _nonfood_mass(net, x),
                "assembly_signature": _assembly_sig(net, x, task),
                "x_eq": x,
            }
            a_spec = (task.get("assembly") or {})
            if a_spec.get("derive_cut"):
                vol = float(a_spec.get(
                    "volume",
                    (task.get("stochastic") or {}).get("volume", 1.0)))
                rec["assembly_calibration"] = assembly.derive_a_cut(
                    assembly.length_mass_spectrum(net, x), vol,
                    alphabet_size=int(net.params.get("alphabet_size", 2)),
                    copies_floor=float(a_spec.get("copies_floor", 1.0)))
        else:
            x = det.x[-1]
            obs = closure.phi_decomposition(net, x, eps)
            rec["final_state"] = {
                "mass_raf": obs["mass_raf"],
                "osc_amplitude": _osc_amplitude(det),
            }
    elif rtype == "sd":
        tau = float(reactor["tau"])
        d_factor = float(reactor["d_factor"])
        n_cycles = int(reactor.get("n_cycles", 40))
        x0 = (x0_det if x0_det is not None
              else net.food_inlet.copy())
        starts, ends, sd_ok, sd_msg = ode.sd_map(
            net, x0, tau, d_factor, n_cycles, eps)
        if not sd_ok:
            # graceful solver-failure record — same contract as the CSTR
            # path (strict-JSON-safe status, no crash, partial cycles kept)
            rec["det"] = {
                "type": "sd", "tau": tau, "d_factor": d_factor,
                "n_cycles": n_cycles,
                "n_cycles_completed": len(starts),
                "status": "UNRESOLVED",
                "message": sd_msg,
            }
        else:
            # fixed-point check: CONSECUTIVE start drift over the last
            # cycles. Comparing starts[-1] to starts[-5] aliases
            # period-2/4 orbits onto fixed points (same phase evaluates
            # drift = 0).
            k_tail = min(4, n_cycles - 1)
            drift = float(max(
                np.max(np.abs(starts[c] - starts[c - 1])
                       / (np.abs(starts[c - 1]) + 1e-12))
                for c in range(n_cycles - 1, n_cycles - 1 - k_tail, -1)))
            stationary = drift < 1e-4
            # propagation diagnostics: is the nonfood composition a
            # persistent lineage signature across dilution cycles, or
            # does it wander? comp_dist_cycle[i] = d_L1(sig_i, sig_{i-1})
            # on end-of-cycle states; a stable lineage -> ~0
            # asymptotically.
            end_sigs = [heredity.attractor_signature(net, e)
                        for e in ends]
            comp_dist = [None] + [
                heredity.compositional_distance(end_sigs[i],
                                                end_sigs[i - 1])
                for i in range(1, n_cycles)]
            seed_sig = heredity.attractor_signature(net, x0)
            dist_seed = heredity.compositional_distance(end_sigs[-1],
                                                        seed_sig)
            # compact lineage record: top-8 species of the final signature
            fin = end_sigs[-1]
            top = np.argsort(fin)[::-1][:8]
            top = top[fin[top] > 0]
            lin = [{"species": str(net.structure.get(i, i)),
                    "frac": float(fin[i])} for i in top]
            rec["det"] = {"type": "sd", "tau": tau,
                          "d_factor": d_factor,
                          "n_cycles": n_cycles,
                          "cycle_start_drift": drift,
                          "m_poly_cycles": [
                              float(_nonfood_mass(net, e)) for e in ends],
                          "comp_dist_cycle": comp_dist,
                          "comp_dist_to_seed": dist_seed,
                          "final_signature_top": lin,
                          "status": "STATIONARY" if stationary else
                                    "UNRESOLVED"}
            if stationary:
                if net.n_species <= 600:
                    mults = ode.sd_fixed_point_multipliers(
                        net, starts[-1], tau, d_factor, eps)
                    if mults is not None:
                        floquet = {
                            "floquet_max_abs": float(
                                np.max(np.abs(mults))),
                            "stable": bool(np.all(np.abs(mults) < 1.0))}
                    else:
                        floquet = {"floquet_max_abs": None,
                                   "stable": None,
                                   "note": "floquet solve failed"}
                else:
                    # FD map Jacobian costs n integrations — infeasible
                    # at sparse scale; report as unverified rather than
                    # silently skipping
                    floquet = {"floquet_max_abs": None,
                               "stable": None,
                               "note": "floquet not computed at n>600"}
                obs = closure.phi_decomposition(net, ends[-1], eps)
                rec["equilibrium"] = {
                    "mass_raf": obs["mass_raf"],
                    "phi_gross": obs["phi_gross"],
                    "phi_cat": obs["phi_cat"],
                    "phi_net": obs["phi_net"],
                    **floquet,
                    "m_poly": _nonfood_mass(net, ends[-1]),
                    "x_eq": ends[-1],
                }
            else:
                obs = closure.phi_decomposition(net, ends[-1], eps)
                rec["final_state"] = {"mass_raf": obs["mass_raf"]}
    else:
        raise ValueError(f"unknown reactor type {rtype!r}")

    # ---- multistability scan (O_hered context) -------------------
    if reactor.get("multistability"):
        mspec = reactor["multistability"]
        min_mass = float(mspec.get("min_mass", 0.01))
        n_ic = int(mspec.get("n_ic", 16))
        ms = _multistability_block(net, dilution, eps, seed_key,
                                   n_ic=n_ic, min_mass=min_mass,
                                   t_max=float(mspec.get(
                                       "t_max", reactor.get(
                                           "t_max", 200.0))))
        _add_per_attractor_obs(net, ms, task)
        rec["multistability"] = ms
    else:
        rec["multistability"] = None

    # ---- continuation / hysteresis sweep -------------------------
    if reactor.get("continuation"):
        rec["continuation"] = _continuation_block(
            net, reactor, dilution, eps, x0_det)
    else:
        rec["continuation"] = None

    # ---- competition SSA (O_select context) ----------------------
    if reactor.get("competition"):
        cspec = reactor["competition"]
        ms_block = rec["multistability"]
        if ms_block is None:
            ms_block = _multistability_block(
                net, dilution, eps, seed_key,
                n_ic=int(cspec.get("n_ic", 16)),
                min_mass=float(cspec.get("min_mass", 0.01)),
                t_max=float(cspec.get("t_max", 200.0)))
            _add_per_attractor_obs(net, ms_block, task)
            rec["multistability"] = ms_block
        rec["competition"] = _competition_block(
            net, dilution, eps, seed_key, ms_block, cspec)
    else:
        rec["competition"] = None

    # ---- stochastic ----
    stoch = task.get("stochastic", {})
    if stoch.get("enabled"):
        volume = float(stoch.get("volume", 100.0))
        replicates = int(stoch.get("replicates", 3))
        t_max = float(stoch.get("t_max", reactor.get("t_max", 200.0)))
        theta = float(stoch.get("theta_organized", 0.1))
        n0 = None
        init_policy = "food"
        if stoch.get("x0_from") == "reactor" and x0_det is not None:
            n0 = gillespie.counts_from_concentration(
                np.asarray(x0_det, dtype=float), volume)
            init_policy = "reactor_x0_round"
        reps = []
        nonfood = ~net.is_food
        comp_mass = (net.comp.sum(axis=0) if net.comp is not None
                     else np.ones(net.n_species))
        mass_floor = float(stoch.get("mass_floor", 0.01))
        hazard_mode = str(stoch.get("hazard_mode", "participants"))
        n_seed_dropped = None
        if n0 is not None:
            n_seed_dropped = int(np.sum(
                (np.asarray(x0_det) > 0) & (n0 == 0)))
        for r_i in range(replicates):
            seed = derive_seed(f"{seed_key}|ssa|{r_i}")
            res = gillespie.simulate(
                net, t_max=t_max, dilution=dilution, volume=volume,
                eps=eps, seed=seed, n0=n0,
                max_events=int(stoch.get("max_events", 2_000_000)),
                record_dt=float(stoch.get("record_dt", t_max / 1000)),
                hazard_mode=hazard_mode)
            o = mass_raf_series(net, res.n, volume)
            nf = res.n[:, nonfood].sum(axis=1)
            # absolute nonfood monomer-mass concentration
            nf_mass = res.n[:, nonfood] @ comp_mass[nonfood] / volume
            # organized = high RAF fraction AND above absolute mass floor
            organized = (o > theta) & (nf_mass >= mass_floor)
            # censor at censor_time: padded records must not count
            mask = res.t <= res.censor_time + 1e-12
            idx_org = np.nonzero(organized & mask)[0]
            # quasi-extinction: first dip below the absolute mass floor
            # after first organization (meaningful at eps>=0; the strict
            # absorbing hit nf==0 is only possible at eps=0). Both are
            # masked to post-organization — an absorbing hit at t=0 on a
            # food-only start is vacuous, not an extinction.
            if len(idx_org):
                post = mask & (np.arange(len(nf_mass)) >= idx_org[0])
                idx_abs = np.nonzero(post & (nf == 0))[0]
                idx_ext = np.nonzero(post & (nf_mass < mass_floor))[0]
            else:
                idx_abs = idx_ext = np.array([], dtype=int)
            reps.append({
                "seed": seed,
                "completed": res.completed,
                "truncated": res.truncated,
                "censor_time": res.censor_time,
                "events": res.events,
                "O_final": float(o[-1]),
                "O_max": float(o.max()),
                "t_first_organized": (float(res.t[idx_org[0]])
                                      if len(idx_org) else None),
                "t_extinction": (float(res.t[idx_ext[0]])
                                 if len(idx_ext) else None),
                "t_absorbing": (float(res.t[idx_abs[0]])
                                if len(idx_abs) else None),
                # dwell as a fraction of the intended horizon t_max:
                # a run truncated while organized at t < t_max is NOT
                # 100% dwell — the unobserved remainder is not credited
                "dwell_organized": min(1.0, float(
                    organized[mask].sum()
                    * res.diagnostics["record_dt"]
                    / max(t_max, 1e-12))),
                "dwell_untruncated": bool(not res.truncated),
            })
        rec["stochastic"] = {
            "volume": volume, "replicates": reps,
            "n_replicates": len(reps),
            "theta_organized": theta,
            "mass_floor": mass_floor,
            "hazard_mode": hazard_mode,
            "n_seed_dropped": n_seed_dropped,
            "init_policy": init_policy,
        }
    else:
        rec["stochastic"] = None

    if analytic:
        out = {}
        for k, v in analytic.items():
            if callable(v):
                try:
                    out[k] = float(v(dilution))
                except TypeError:
                    out[k] = "callable"
            else:
                out[k] = v
        rec["analytic"] = out
    rec["wall_s"] = time.time() - t0
    return strict_jsonable(rec)


def _osc_amplitude(det) -> float:
    w = det.x[-max(3, len(det.x) // 8):]
    span = np.abs(w.mean(axis=0)) + 1e-12
    return float(np.max((w.max(axis=0) - w.min(axis=0)) / span))


def _nonfood_mass(net: ReactionNetwork, x: np.ndarray) -> float:
    w = np.where(net.is_food, 0.0, 1.0)
    if net.comp is not None:
        w = w * net.comp.sum(axis=0)
    return float(w @ np.maximum(np.asarray(x, float), 0.0))


def _assembly_sig(net: ReactionNetwork, x: np.ndarray,
                  task: dict) -> dict:
    """Absolute-threshold assembly signature (derived cut, absolute mass
    floor). task['assembly'] may set a_cut / mass_floor overrides."""
    a = (task or {}).get("assembly", {})
    return assembly.assembly_signature_absolute(
        net, np.asarray(x, float),
        a_cut=a.get("a_cut"),
        mass_floor=float(a.get("mass_floor", 0.01)))


def _continuation_block(net: ReactionNetwork, reactor: dict,
                        dilution: float, eps: float,
                        x0_det: np.ndarray | None) -> dict:
    """Forward/reverse warm-start sweep for hysteresis measurement.

    reactor["continuation"] = {"param": "dilution"|"eps",
                               "values": [...], "mass_floor": float}
    Each point integrates from the previous point's endpoint; the rev
    branch starts from the fwd branch's final state (true continuation).
    Generator-param axes are not supported here — the network is fixed
    inside a task; use a campaign sweep for those.

    A "transition" is the first λ where organized status (m_poly >
    mass_floor) flips relative to the branch's starting state. Folds
    located this way are SPINODALS — the coexistence point needs the
    stochastic/quasi-potential analysis; report width, not a single
    'transition point'.
    """
    cont = reactor["continuation"]
    param = cont["param"]
    if param not in ("dilution", "eps"):
        raise ValueError(
            f"continuation param {param!r} unsupported "
            "(only 'dilution' and 'eps' are integrator-level axes)")
    values = [float(v) for v in cont["values"]]
    floor = float(cont.get("mass_floor", 0.01))
    t_max = float(cont.get("t_max", reactor.get("t_max", 200.0)))

    branches = {}
    rev_reseeded = False
    x = (np.asarray(x0_det, float).copy() if x0_det is not None
         else net.food_inlet.copy())
    for direction, seq in (("fwd", values), ("rev", values[::-1])):
        # At eps == 0 a dead fwd endpoint is absorbing: uncatalyzed
        # synthesis is zero, so warm-starting the rev branch from it
        # can never re-ignite and hysteresis_width is always None.
        # Reseed the rev branch once from the organized reactor seed
        # (or a dense nonfood seed if none) to probe the lower
        # spinodal. The eps actually used at the first rev point is
        # seq[0] when the swept parameter IS eps.
        eps_start = seq[0] if param == "eps" else eps
        if (direction == "rev" and eps_start == 0.0
                and _nonfood_mass(net, x) < floor):
            x = (np.asarray(x0_det, float).copy()
                 if x0_det is not None else net.food_inlet.copy())
            if _nonfood_mass(net, x) < floor:
                x = net.food_inlet.copy()
                x[~net.is_food] = 0.15
            rev_reseeded = True
        rows = []
        for lam in seq:
            kw = {"dilution": dilution, "eps": eps}
            kw[param] = lam
            r = ode.integrate_cstr(net, t_max=t_max, x0=x, **kw)
            x = r.x_eq if r.converged else r.x[-1]
            rows.append({"lam": lam,
                         "status": r.trajectory_type,
                         "m_poly": _nonfood_mass(net, x)})
        branches[direction] = rows

    out = {"param": param, "values": values, "mass_floor": floor,
           "branches": branches, "rev_reseeded": rev_reseeded}
    transitions = {}
    for direction in ("fwd", "rev"):
        rows = branches[direction]
        org = [row["m_poly"] > floor for row in rows]
        base = org[0]
        flips = [i for i, o in enumerate(org) if o != base]
        if flips:
            i = flips[0]
            transitions[direction] = {
                "lam": rows[i]["lam"],
                "lc_lo": rows[i - 1]["lam"] if i else None,
                "branch_started_organized": base}
        else:
            transitions[direction] = None
    out["transitions"] = transitions
    tf, tr = transitions["fwd"], transitions["rev"]

    # Genuine hysteresis requires the two branches to occupy DIFFERENT
    # states at the same control value. A single-valued response curve
    # that crosses mass_floor always produces one fwd flip and one rev
    # flip separated by exactly one grid step — that is a threshold
    # bracket, not a fold pair. Verify branch divergence at common lams.
    fwd_by_lam = {row["lam"]: row["m_poly"] for row in branches["fwd"]}
    rev_by_lam = {row["lam"]: row["m_poly"] for row in branches["rev"]}
    common_lams = sorted(set(fwd_by_lam) & set(rev_by_lam))
    state_diffs = {lam: abs(fwd_by_lam[lam] - rev_by_lam[lam])
                   for lam in common_lams}
    max_branch_diff = max(state_diffs.values()) if state_diffs else 0.0
    bistable_lams = [lam for lam, d in state_diffs.items()
                     if d > 2.0 * floor]
    out["max_branch_diff"] = float(max_branch_diff)
    out["n_bistable_points"] = len(bistable_lams)

    if tf and tr and bistable_lams:
        # fold pair located AND branch states actually diverge ->
        # hysteresis window between the spinodals
        out["hysteresis_width"] = abs(tf["lam"] - tr["lam"])
        out["spinodal_fwd"], out["spinodal_rev"] = tf["lam"], tr["lam"]
        out["bistable_window"] = [float(min(bistable_lams)),
                                  float(max(bistable_lams))]
    elif tf and tr:
        # both branches crossed the floor, but states coincide at every
        # shared lam -> single-valued response, NOT hysteresis
        out["hysteresis_width"] = None
        out["spinodal_fwd"] = out["spinodal_rev"] = None
        out["note"] = ("no genuine branch bifurcation detected "
                       "(single-valued response)")
    else:
        out["hysteresis_width"] = None
        out["spinodal_fwd"] = out["spinodal_rev"] = None
        out["note"] = ("no fold pair located within the sweep; widen "
                       "values or check mass_floor")
    return out


def _lineage_seed(net: ReactionNetwork, seed_key: str, letter: str,
                  level: float = 0.08) -> np.ndarray:
    """Species group dominated by `letter` (>=2/3 composition) at
    `level` — the same lineage-biased rule as _ic_bank, usable as a
    standalone reactor x0_policy ("lineage:a"). Falls back to a dense
    seed if the letter group is empty."""
    rng = np.random.default_rng(
        derive_seed(seed_key + f"|lineage:{letter}"))
    nonfood = np.where(~net.is_food)[0]
    x = net.food_inlet.copy()
    grp = np.array(
        [i for i in nonfood
         if isinstance(net.structure.get(i), str)
         and len(net.structure[i]) > 1
         and net.structure[i].count(letter)
         / len(net.structure[i]) >= 2 / 3],
        dtype=int)
    if len(grp):
        x[grp] = level * (0.5 + 0.5 * rng.random(len(grp)))
    else:
        x[nonfood] = level
    return x


def _ic_bank(net: ReactionNetwork, seed_key: str, n_ic: int,
             ic_scales=(1e-3, 5e-2, 0.5)) -> list[np.ndarray]:
    """Deterministic IC bank: food-only + dense + lineage-biased +
    sparse random draws.

    Sparse seeds can sit below the p>=2 separatrix; the bank mixes
    scales so both basin structures are probed. Lineage-biased seeds
    (one per letter, species where that letter is >= 2/3 of the
    composition) probe lineage-basin boundaries that random sparse
    seeds can miss — they are how p_info's heredity-protecting claim
    gets tested rather than assumed.
    """
    rng = np.random.default_rng(derive_seed(seed_key + "|icbank"))
    nonfood = np.where(~net.is_food)[0]
    ics = [net.food_inlet.copy()]
    dense = net.food_inlet.copy()
    dense[nonfood] = 0.15
    ics.append(dense)
    if net.structure is not None:
        letters = sorted({s for s in net.structure.values()
                          if isinstance(s, str) and len(s) == 1})
        for letter in letters:
            grp = np.array(
                [i for i in nonfood
                 if isinstance(net.structure.get(i), str)
                 and len(net.structure[i]) > 1
                 and net.structure[i].count(letter)
                 / len(net.structure[i]) >= 2 / 3],
                dtype=int)
            if len(grp):
                x = net.food_inlet.copy()
                x[grp] = 0.08 * (0.5 + 0.5 * rng.random(len(grp)))
                ics.append(x)
    for k in range(max(0, n_ic - len(ics))):
        x = net.food_inlet.copy()
        scale = ic_scales[k % len(ic_scales)]
        size = max(2, int(rng.integers(2, min(12, len(nonfood) + 1))))
        pick = rng.choice(nonfood, size=size, replace=False)
        x[pick] = scale * (0.1 + 0.9 * rng.random(size))
        ics.append(x)
    return ics


def _multistability_block(net: ReactionNetwork, dilution: float,
                          eps: float, seed_key: str, *,
                          n_ic: int, min_mass: float,
                          t_max: float) -> dict:
    """Certified multi-IC scan -> attractor clustering + member x_eqs.

    Returns the clustering plus a representative certified equilibrium
    vector per organized attractor (used to seed competition SSA).
    """
    ics = _ic_bank(net, seed_key, n_ic)
    eqs, status, organized_ic = [], [], []
    for k, x0 in enumerate(ics):
        r = ode.integrate_cstr(net, x0=x0, dilution=dilution, eps=eps,
                               t_max=t_max)
        status.append(r.trajectory_type)
        if r.converged:
            eqs.append(r.x_eq)
            organized_ic.append(k)
    # stability certification: a "STATIONARY" endpoint with a positive
    # leading Jacobian eigenvalue is an unstable fixed point or slow
    # transient, not an attractor — exclude it (ghost class #2).
    stable_eqs, stable_ic = [], []
    for x, k in zip(eqs, organized_ic):
        if _nonfood_mass(net, x) < min_mass:
            continue
        eigs, spec_conv = spectral.jacobian_eigs(net, x, dilution, eps)
        # a partially-converged spectrum may be missing the leading
        # unstable eigenvalue — do not certify stability on it
        if spec_conv and eigs[0].real < 0:
            stable_eqs.append(x)
            stable_ic.append(k)
    cl = heredity.cluster_attractors(net, stable_eqs, tol=0.05,
                                     min_mass=min_mass)
    # joint composition+mass clustering: amplitude bistability (high/low
    # branches with the same composition) is invisible to normalized
    # signatures — report the joint count separately
    cl_total = heredity.cluster_attractors(net, stable_eqs, tol=0.05,
                                           min_mass=min_mass,
                                           mass_ratio_tol=1.5)
    # representative member x_eq per organized cluster
    centers_x = {}
    for c in range(cl["n_attractors"]):
        idx = [i for i, l in enumerate(cl["labels"]) if l == c]
        if idx:
            centers_x[c] = stable_eqs[idx[0]]
    return {
        "n_attractors": cl["n_attractors"],
        "n_attractors_total": cl_total["n_attractors"],
        "n_certified": len(eqs),
        "n_stable": len(stable_eqs),
        "n_kept": cl["n_kept"],
        "n_ics": len(ics),
        "ic_status": status,
        "organized_ic": stable_ic,
        "cluster_labels": cl["labels"].tolist(),
        "members": cl["members"],
        "center_distances": heredity._pairwise(cl["centers"]),
        "min_mass": min_mass,
        "attractor_x0": {str(c): x.tolist() for c, x in centers_x.items()},
        "attractor_masses": {str(c): float(_nonfood_mass(net, x))
                             for c, x in centers_x.items()},
        "attractor_top_species": {
            str(c): net.species[int(np.argmax(
                heredity.attractor_signature(net, x)))]
            for c, x in centers_x.items()},
    }


def _add_per_attractor_obs(net: ReactionNetwork, ms: dict,
                           task: dict) -> None:
    """Basin-conditioned observables: evaluate m_poly and the assembly
    signature on EACH certified attractor center, not just whichever
    basin the reactor's policy seed landed in. Closes the
    basin-conditioning gap where O_maint/O_assembly read only the
    primary integration endpoint while O_hered uses the full IC scan.
    """
    obs = []
    for cid, x0 in sorted(ms.get("attractor_x0", {}).items(),
                          key=lambda kv: int(kv[0])):
        x = np.asarray(x0, float)
        sig = _assembly_sig(net, x, task)
        obs.append({"attractor": int(cid),
                    "m_poly": float(_nonfood_mass(net, x)),
                    "mass_signature_abs": sig["mass_signature_abs"],
                    "mass_signature_frac": sig["mass_signature_frac"],
                    "n_signature": sig["n_signature"]})
    ms["per_attractor_obs"] = obs


def _competition_block(net: ReactionNetwork, dilution: float,
                       eps: float, seed_key: str,
                       ms_block: dict | None, cspec: dict) -> dict:
    """Mixture-competition SSA between the first two organized
    attractors. Two initializations per plan:

      raw:     n0 ~ V * (x_A + x_B)/2              (natural basins)
      eq_mass: n0 ~ V * M_t * (s_A + s_B)/2        (equal monomer mass)

    Winner per rep = nearest attractor signature at the recorded final
    state, labeled DEAD if final nonfood mass < min_mass floor.
    """
    if ms_block is None or ms_block["n_attractors"] < 2:
        return {"completed": False,
                "reason": "fewer_than_two_organized_attractors"}
    min_mass = float(cspec.get("min_mass", ms_block["min_mass"]))
    volume = float(cspec.get("volume", 16.0))
    n_reps = int(cspec.get("replicates", 30))
    t_max = float(cspec.get("t_max", 400.0))
    tol = float(cspec.get("tol", 0.05))
    hazard_mode = str(cspec.get("hazard_mode", "participants"))
    tie_margin = float(cspec.get("tie_margin", tol))

    # caller may nominate the competing attractor pair (cluster ids as
    # strings); default is the first two certified attractors
    c_pair = [str(c) for c in cspec.get("attractor_pair", ("0", "1"))]
    if any(c not in ms_block["attractor_x0"] for c in c_pair):
        return {"completed": False,
                "reason": "requested_attractor_pair_not_certified"}
    xA = np.asarray(ms_block["attractor_x0"][c_pair[0]], float)
    xB = np.asarray(ms_block["attractor_x0"][c_pair[1]], float)
    sigs = [heredity.attractor_signature(net, xA),
            heredity.attractor_signature(net, xB)]
    m_target = min(_nonfood_mass(net, xA), _nonfood_mass(net, xB))
    # signatures are monomer-mass fractions; dividing by monomer mass
    # per molecule converts the equal-mass target into molecular
    # concentrations — counts_from_concentration takes molecule counts
    comp_mass = (net.comp.sum(axis=0) if net.comp is not None
                 else np.ones(net.n_species))

    out = {"completed": True, "volume": volume, "n_reps": n_reps,
           "attractor_pair": c_pair,
           "attractor_labels": ms_block["attractor_top_species"],
           "hazard_mode": hazard_mode,
           "protocols": {}}
    for pname, n0f in (
            ("raw", lambda V: gillespie.counts_from_concentration(
                0.5 * (xA + xB), V)),
            ("eq_mass", lambda V: gillespie.counts_from_concentration(
                0.5 * m_target * (sigs[0] + sigs[1])
                / np.maximum(comp_mass, 1.0)
                + net.food_inlet, V))):
        n0 = n0f(volume)
        counts = {"A": 0, "B": 0, "dead": 0, "other": 0,
                  "undecided": 0}
        rep_rows = []
        for rep in range(n_reps):
            seed = derive_seed(f"{seed_key}|comp|{pname}|{rep}")
            res = gillespie.simulate(
                net, t_max=t_max, dilution=dilution, volume=volume,
                n0=n0, eps=eps, seed=seed,
                max_events=int(cspec.get("max_events", 4_000_000)),
                record_dt=t_max / 500,
                hazard_mode=hazard_mode)
            xf = res.final_counts / volume
            if _nonfood_mass(net, xf) < min_mass:
                label = "dead"
            else:
                sf = heredity.attractor_signature(net, xf)
                d = [heredity.compositional_distance(sf, s)
                     for s in sigs]
                if min(d) > tol * 4:
                    label = "other"
                elif abs(d[0] - d[1]) <= tie_margin:
                    # terminal composition sits (near-)equidistant
                    # between the attractors — an undecided replicate,
                    # not a coin-flip win for whichever lineage was
                    # indexed first (was: d[0] <= d[1] -> "A")
                    label = "undecided"
                else:
                    label = "A" if d[0] < d[1] else "B"
            counts[label] += 1
            rep_rows.append({"seed": seed, "winner": label,
                             "truncated": res.truncated,
                             "censor_time": res.censor_time})
        n_ab = counts["A"] + counts["B"]
        n_total = len(rep_rows)
        out["protocols"][pname] = {
            "counts": counts, "n_decided": n_ab,
            "p_a_wins": (counts["A"] / n_ab if n_ab else None),
            "p_extinct": counts["dead"] / n_total,
            "p_other": counts["other"] / n_total,
            "p_undecided": counts["undecided"] / n_total,
            "tie_margin": tie_margin,
            "reps": rep_rows}
    # headline: raw protocol asymmetry (None if undecided)
    raw = out["protocols"]["raw"]
    out["p_a_wins"] = raw["p_a_wins"]
    out["init_protocol"] = "raw"
    return out
