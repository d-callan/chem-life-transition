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

from .networks.core import ReactionNetwork
from .networks.generators import polymer_network
from .networks import benchmarks as bench
from .dynamics import gillespie, ode
from .observables import assembly, closure, heredity, spectral
from .provenance import run_id, strict_jsonable


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
        "n_species": net.n_species,
        "n_channels": net.n_channels,
        "n_catalyzed_channels": int(sum(bool(c.catalysts)
                                        for c in net.channels)),
    }
    rec["run_id"] = run_id(net.hash, {"gen": gen_spec, "reactor": reactor,
                                      "eps": eps,
                                      "stochastic": task.get("stochastic")},
                           seed_key)

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
    if rtype == "structural":
        # E1: structural landscape only — RAF topology without dynamics.
        rec["det"] = {"type": "structural", "status": "NA"}
    elif rtype == "cstr":
        det = ode.integrate_cstr(
            net, t_max=float(reactor.get("t_max", 200.0)),
            dilution=dilution, eps=eps,
            x0=np.asarray(reactor["x0"]) if "x0" in reactor else None,
            n_points=int(reactor.get("n_points", 400)))
        rec["det"] = {
            "type": rtype,
            "status": det.trajectory_type,
            "success": det.success,
            "message": det.message,
            "nfev": det.nfev,
            "residual_scaled_max": float(np.nanmax(det.residual_scaled)),
            "window_change": det.window_change,
            "material_error": det.material_error,
        }
        if det.converged:
            x = det.x_eq
            obs = closure.phi_decomposition(net, x, eps)
            eigs = spectral.jacobian_eigs(net, x, dilution, eps)
            rec["equilibrium"] = {
                "mass_raf": obs["mass_raf"],
                "phi_gross": obs["phi_gross"],
                "phi_cat": obs["phi_cat"],
                "phi_net": obs["phi_net"],
                "denominators": obs["denominators"],
                "n_raf_channels": obs["n_raf_channels"],
                "leading_eig_re": float(eigs[0].real),
                "nonzero_mass_species": int(np.sum(x > 1e-10)),
                "m_poly": _nonfood_mass(net, x),
                "assembly_signature": _assembly_sig(net, x, seed_key),
                "x_eq": x,
            }
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
        x0 = (np.asarray(reactor["x0"]) if "x0" in reactor
              else net.food_inlet.copy())
        starts, ends = ode.sd_map(net, x0, tau, d_factor, n_cycles, eps)
        # fixed-point check: start drift over last 5 cycles
        drift = float(np.max(np.abs(starts[-1] - starts[-5])
                             / (np.abs(starts[-5]) + 1e-12)))
        stationary = drift < 1e-4
        rec["det"] = {"type": "sd", "tau": tau, "d_factor": d_factor,
                      "n_cycles": n_cycles, "cycle_start_drift": drift,
                      "status": "STATIONARY" if stationary else
                                "UNRESOLVED"}
        if stationary:
            mults = ode.sd_fixed_point_multipliers(net, starts[-1], tau,
                                                   d_factor, eps)
            obs = closure.phi_decomposition(net, ends[-1], eps)
            rec["equilibrium"] = {
                "mass_raf": obs["mass_raf"],
                "phi_gross": obs["phi_gross"],
                "phi_cat": obs["phi_cat"],
                "phi_net": obs["phi_net"],
                "floquet_max_abs": float(np.max(np.abs(mults))),
                "stable": bool(np.all(np.abs(mults) < 1.0)),
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
        rec["multistability"] = ms
    else:
        rec["multistability"] = None

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
        if stoch.get("x0_from") == "reactor" and "x0" in reactor:
            n0 = gillespie.counts_from_concentration(
                np.asarray(reactor["x0"], dtype=float), volume)
            init_policy = "reactor_x0_round"
        reps = []
        nonfood = ~net.is_food
        for r_i in range(replicates):
            seed = derive_seed(f"{seed_key}|ssa|{r_i}")
            res = gillespie.simulate(
                net, t_max=t_max, dilution=dilution, volume=volume,
                eps=eps, seed=seed, n0=n0,
                max_events=int(stoch.get("max_events", 2_000_000)),
                record_dt=float(stoch.get("record_dt", t_max / 1000)))
            o = mass_raf_series(net, res.n, volume)
            nf = res.n[:, nonfood].sum(axis=1)
            above = o > theta
            idx_org = np.nonzero(above)[0]
            if nf[0] > 0:
                # seeded start: first time all nonfood counts hit zero
                idx_ext = np.nonzero(nf == 0)[0]
            elif len(idx_org):
                # unseeded: only post-organization extinction is meaningful
                idx_ext = np.nonzero(nf[idx_org[0]:] == 0)[0] + idx_org[0]
            else:
                idx_ext = np.array([], dtype=int)
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
                "dwell_organized": float(above.mean()),
            })
        rec["stochastic"] = {
            "volume": volume, "replicates": reps,
            "n_replicates": len(reps),
            "theta_organized": theta,
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
                  seed_key: str) -> dict:
    """Null-calibrated joint assembly+abundance signature at x."""
    rng = np.random.default_rng(derive_seed(seed_key + "|manull"))
    lengths = sorted({len(s) for s in net.structure.values()})
    thr = {L: assembly.ma_threshold(
        assembly.null_ma_distribution(
            rng, L, alphabet_size=int(net.params.get(
                "alphabet_size", 2))))
        for L in lengths}
    return assembly.assembly_signature(net, np.asarray(x, float),
                                       thr, abundance_floor=0.05)


def _ic_bank(net: ReactionNetwork, seed_key: str, n_ic: int,
             ic_scales=(1e-3, 5e-2, 0.5)) -> list[np.ndarray]:
    """Deterministic IC bank: food-only + dense + sparse random draws.

    Sparse seeds can sit below the p>=2 separatrix; the bank mixes
    scales so both basin structures are probed.
    """
    rng = np.random.default_rng(derive_seed(seed_key + "|icbank"))
    nonfood = np.where(~net.is_food)[0]
    ics = [net.food_inlet.copy()]
    dense = net.food_inlet.copy()
    dense[nonfood] = 0.15
    ics.append(dense)
    for k in range(n_ic - 2):
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
    cl = heredity.cluster_attractors(net, eqs, tol=0.05,
                                     min_mass=min_mass)
    # representative member x_eq per organized cluster
    centers_x = {}
    kept = [e for e in eqs if _nonfood_mass(net, e) >= min_mass]
    for c in range(cl["n_attractors"]):
        idx = [i for i, l in enumerate(cl["labels"]) if l == c]
        if idx:
            centers_x[c] = kept[idx[0]]
    return {
        "n_attractors": cl["n_attractors"],
        "n_certified": len(eqs),
        "n_kept": cl["n_kept"],
        "n_ics": len(ics),
        "ic_status": status,
        "organized_ic": organized_ic,
        "cluster_labels": cl["labels"].tolist(),
        "members": cl["members"],
        "center_distances": heredity._pairwise(cl["centers"]),
        "min_mass": min_mass,
        "attractor_x0": {str(c): x.tolist() for c, x in centers_x.items()},
        "attractor_top_species": {
            str(c): net.species[int(np.argmax(
                heredity.attractor_signature(net, x)))]
            for c, x in centers_x.items()},
    }


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

    xA = np.asarray(ms_block["attractor_x0"]["0"], float)
    xB = np.asarray(ms_block["attractor_x0"]["1"], float)
    sigs = [heredity.attractor_signature(net, xA),
            heredity.attractor_signature(net, xB)]
    m_target = min(_nonfood_mass(net, xA), _nonfood_mass(net, xB))

    out = {"completed": True, "volume": volume, "n_reps": n_reps,
           "attractor_labels": ms_block["attractor_top_species"],
           "protocols": {}}
    for pname, n0f in (
            ("raw", lambda V: gillespie.counts_from_concentration(
                0.5 * (xA + xB), V)),
            ("eq_mass", lambda V: gillespie.counts_from_concentration(
                0.5 * m_target * (sigs[0] + sigs[1])
                + net.food_inlet, V))):
        n0 = n0f(volume)
        counts = {"A": 0, "B": 0, "dead": 0, "other": 0}
        rep_rows = []
        for rep in range(n_reps):
            seed = derive_seed(f"{seed_key}|comp|{pname}|{rep}")
            res = gillespie.simulate(
                net, t_max=t_max, dilution=dilution, volume=volume,
                n0=n0, eps=eps, seed=seed,
                max_events=int(cspec.get("max_events", 4_000_000)),
                record_dt=t_max / 500)
            xf = res.final_counts / volume
            if _nonfood_mass(net, xf) < min_mass:
                label = "dead"
            else:
                sf = heredity.attractor_signature(net, xf)
                d = [heredity.compositional_distance(sf, s)
                     for s in sigs]
                label = ("A" if d[0] <= d[1] else "B") \
                    if min(d) <= tol * 4 else "other"
            counts[label] += 1
            rep_rows.append({"seed": seed, "winner": label,
                             "truncated": res.truncated,
                             "censor_time": res.censor_time})
        n_ab = counts["A"] + counts["B"]
        out["protocols"][pname] = {
            "counts": counts, "n_decided": n_ab,
            "p_a_wins": (counts["A"] / n_ab if n_ab else None),
            "reps": rep_rows}
    # headline: raw protocol asymmetry (None if undecided)
    raw = out["protocols"]["raw"]
    out["p_a_wins"] = raw["p_a_wins"]
    out["init_protocol"] = "raw"
    return out
