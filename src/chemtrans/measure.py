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
from .observables import closure, spectral
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
                                      "eps": eps}, seed_key)

    # ---- deterministic ----
    if rtype == "cstr":
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
                "x_eq": ends[-1],
            }
        else:
            obs = closure.phi_decomposition(net, ends[-1], eps)
            rec["final_state"] = {"mass_raf": obs["mass_raf"]}
    else:
        raise ValueError(f"unknown reactor type {rtype!r}")

    # ---- stochastic ----
    stoch = task.get("stochastic", {})
    if stoch.get("enabled"):
        volume = float(stoch.get("volume", 100.0))
        replicates = int(stoch.get("replicates", 3))
        t_max = float(stoch.get("t_max", reactor.get("t_max", 200.0)))
        reps = []
        for r_i in range(replicates):
            seed = derive_seed(f"{seed_key}|ssa|{r_i}")
            res = gillespie.simulate(
                net, t_max=t_max, dilution=dilution, volume=volume,
                eps=eps, seed=seed,
                max_events=int(stoch.get("max_events", 2_000_000)),
                record_dt=float(stoch.get("record_dt", t_max / 1000)))
            o = mass_raf_series(net, res.n, volume)
            reps.append({
                "seed": seed,
                "completed": res.completed,
                "truncated": res.truncated,
                "censor_time": res.censor_time,
                "events": res.events,
                "O_final": float(o[-1]),
                "O_max": float(o.max()),
            })
        rec["stochastic"] = {
            "volume": volume, "replicates": reps,
            "n_replicates": len(reps),
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
