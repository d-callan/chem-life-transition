"""Stochastic simulation (Gillespie SSA) for small networks.

Deterministic ODEs can't capture extinction: an autocatalytic core at
concentration 1e-9 is "alive" in the ODE but dead in reality. SSA is needed
to measure genuine persistence vs. numerical persistence, and to see
demographic noise near the transition.

Intended for SMALL networks (hundreds of reactions): each event rescores
all propensities, so cost is O(n_reactions) per firing — fine for
random_crn toy nets, impractical for polymer_ligation (~6000 rxns and
thousands of events/second). For big nets use tau-leaping (not yet
implemented).
"""

import numpy as np

from ..networks.core import ReactionNetwork


def simulate(
    net: ReactionNetwork,
    t_max: float,
    feed_rate: float = 1.0,
    dilution_rate: float = 0.1,
    volume: float = 100.0,
    x0: np.ndarray | None = None,
    seed: int | None = None,
    max_events: int = 5_000_000,
    record_dt: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Run SSA with food influx and dilution as Poisson/death processes.

    Concentrations are counts / volume; reaction propensity a_r evaluated
    by net.propensity at x = n/V then multiplied by V (event/time).
    Food influx: feed_rate*V molecules/time per food species.
    Dilution: each molecule leaves at dilution_rate (death process).

    Returns (t, X): t are record times, X molecule counts per species.
    If the trajectory hits max_events the tail is truncated — treat
    long organized trajectories as 'running until cap'.
    """
    rng = np.random.default_rng(seed)
    n_species = net.n_species
    food_idx = np.array(sorted(net.food), dtype=int)
    influx_rate = feed_rate * volume

    if x0 is None:
        n = np.zeros(n_species)
        n[food_idx] = round(feed_rate * volume / max(dilution_rate, 1e-12))
    else:
        n = np.floor(np.asarray(x0, dtype=float)).copy()
    n = np.maximum(n, 0.0)

    if record_dt is None:
        record_dt = t_max / 2000.0

    ts, xs = [0.0], [n.copy()]
    t = 0.0
    t_next_record = record_dt
    events = 0
    while t < t_max and events < max_events:
        x = n / volume
        a = net.propensity(x) * volume
        # process list: [reactions | food influx | dilution per species]
        death = dilution_rate * n
        a_in = influx_rate * np.ones(len(food_idx))
        total = a.sum() + a_in.sum() + death.sum()
        if total <= 0.0:
            break  # nothing can happen (empty reactor, no feed)
        t += float(rng.exponential(1.0 / total))
        if t > t_max:
            break
        u = float(rng.random()) * total
        a_tot = a.sum()
        if u < a_tot:
            r = int(np.searchsorted(np.cumsum(a), u, side="right"))
            rxn = net.reactions[r]
            for i, nu in rxn.reactants.items():
                n[i] -= nu
            for i, nu in rxn.products.items():
                n[i] += nu
        elif u < a_tot + a_in.sum():
            j = int((u - a_tot) // max(influx_rate, 1e-12))
            n[food_idx[min(j, len(food_idx) - 1)]] += 1.0
        else:
            i = int(np.searchsorted(np.cumsum(death),
                                    u - a_tot - a_in.sum(), side="right"))
            n[min(i, n_species - 1)] -= 1.0
        np.maximum(n, 0.0, out=n)
        events += 1
        while t >= t_next_record and t_next_record <= t_max:
            ts.append(t_next_record)
            xs.append(n.copy())
            t_next_record += record_dt
    return np.asarray(ts), np.asarray(xs)
