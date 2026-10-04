"""Shared run-and-measure helper for experiments.

One call = generate nothing, simulate one network under given chemostat
conditions, compute all milestone-1 order parameters. Every experiment
script should use this so the observable definitions and the output schema
stay consistent across sweeps, hysteresis ramps, and finite-size runs.
"""

import warnings

import numpy as np

from .dynamics.ode import simulate
from .observables.assembly import assembly_distribution, path_assembly_index
from .observables.closure import (
    autocatalytic_activity,
    organizational_closure,
    raf_core,
)
from .observables.spectral import spectral_radius


def measure_run(
    net,
    feed_rate: float,
    dilution_rate: float,
    t_max: float = 150.0,
    x0: np.ndarray | None = None,
    with_persistence: bool = False,
    persistence_seed: int | None = None,
) -> dict:
    """Simulate net to steady state and return the order-parameter record.

    The returned dict always contains rho_K plus run status flags; on a
    successful integration it also contains phi, A, n_core_species, the
    assembly summary, and x_ss (final state, for continuation).
    """
    row = {"rho_K": spectral_radius(net), "runaway": False}
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _, X = simulate(
                net, (0.0, t_max),
                feed_rate=feed_rate, dilution_rate=dilution_rate, x0=x0,
            )
            msgs = [str(x.message) for x in w]
            row["runaway"] = any("runaway" in m for m in msgs)
            row["steady_state_reached"] = not any(
                "may not be steady" in m or "runaway" in m for m in msgs
            )
    except RuntimeError as e:
        row["steady_state_reached"] = False
        row["error"] = str(e)
        return row

    x_ss = X[~np.isnan(X).any(axis=1)][-1]
    flux = net.propensity(x_ss)
    A_idx = path_assembly_index(net)
    edges, hist = assembly_distribution(net, x_ss)
    finite = A_idx[np.isfinite(A_idx)]

    row.update({
        "phi": organizational_closure(net, flux),
        "A": autocatalytic_activity(net, x_ss),
        "n_core_species": len(raf_core(net)),
        "assembly": {
            "mean": float(np.average(
                finite, weights=x_ss[np.isfinite(A_idx)] + 1e-30
            )) if len(finite) else None,
            "max": float(finite.max()) if len(finite) else None,
            "tail_mass": float(
                hist[len(hist) // 2:].sum() / max(hist.sum(), 1e-30)
            ) if len(hist) > 1 else 0.0,
            "hist_edges": edges.tolist(),
            "hist": hist.tolist(),
        },
        "x_ss": x_ss.tolist(),
    })

    if with_persistence and not row["runaway"]:
        from .observables.persistence import (
            evolvability_proxy,
            reproductive_number,
            self_maintenance,
        )
        row["S"] = self_maintenance(net, x_ss, feed_rate, dilution_rate)
        row["R0"] = reproductive_number(net, x_ss, dilution_rate)
        row["E"] = evolvability_proxy(net, seed=persistence_seed)
    return row
