"""Persistence and self-maintenance S(lambda), heredity/evolvability H/E.

self_maintenance
    Deepest dilution shock the organized state survives: re-simulate from
    x_ss * (1 - s) for increasing s and find where Phi fails to recover.
    The dynamical definition of self-maintenance — does the RAF core
    re-establish after losing most of its material?
reproductive_number
    R0 of the organized state via a next-generation matrix:
    G[c, j] = (flux of j catalyzed by c) / (dilution * x_c) — expected
    "offspring" species-production events per unit lifetime of c.
    R0 = rho(G): >1 means the catalytic population grows.
evolvability_proxy
    Fraction of single-edge catalyst mutations that change the RAF core —
    does heritable phenotypic variation exist at all?
"""

import numpy as np

from ..dynamics.ode import simulate
from ..networks.core import Reaction, ReactionNetwork
from .closure import organizational_closure, raf_core, raf_reactions


def self_maintenance(
    net: ReactionNetwork,
    x_ss: np.ndarray,
    feed_rate: float = 1.0,
    dilution_rate: float = 0.1,
    t_recover: float = 100.0,
    shocks: tuple = (0.5, 0.7, 0.9, 0.99, 0.999),
) -> float:
    """S(lambda): deepest fractional dilution shock survived, in [0,1].

    0.0 = dies to any tested shock; ~1.0 = recovers from near-total loss.
    """
    flux0 = net.propensity(x_ss)
    phi0 = organizational_closure(net, flux0)
    if phi0 <= 0.0:
        return 0.0
    survived = 0.0
    for s in shocks:
        try:
            _, X = simulate(net, (0.0, t_recover), feed_rate=feed_rate,
                            dilution_rate=dilution_rate,
                            x0=x_ss * (1.0 - s), steady_check=False)
        except RuntimeError:
            break
        last = X[~np.isnan(X).any(axis=1)]
        if len(last) == 0:
            break
        phi_post = organizational_closure(net, net.propensity(last[-1]))
        if phi_post >= 0.5 * phi0:
            survived = s
        else:
            break
    return float(survived)


def reproductive_number(
    net: ReactionNetwork, x: np.ndarray, dilution_rate: float = 0.1
) -> float:
    """R0 = spectral radius of the next-generation matrix at state x."""
    x = np.maximum(np.asarray(x, dtype=float), 1e-12)
    G = np.zeros((net.n_species, net.n_species))
    for rxn in net.reactions:
        if rxn.k_sat is None:
            mass = 1.0
            for i, nu in rxn.reactants.items():
                mass *= x[i] ** nu
        else:
            mass = 1.0
            for i, nu in rxn.reactants.items():
                mass *= (x[i] / (rxn.k_sat + x[i])) ** nu
        for c in rxn.catalysts:
            for j in rxn.products:
                # per-capita production rate of j catalyzed by c, times the
                # dilution lifetime 1/dilution_rate
                if rxn.k_sat is None:
                    G[c, j] += rxn.cat_rate * mass / dilution_rate
                else:
                    G[c, j] += (rxn.cat_rate * mass
                                / (rxn.k_sat + x[c]) / dilution_rate)
    if not np.any(G):
        return 0.0
    return float(max(abs(np.linalg.eigvals(G))))


def evolvability_proxy(
    net: ReactionNetwork, n_mut: int = 20, seed: int | None = None
) -> float:
    """Fraction of single-catalyst-edge mutations that alter the RAF core.

    Heredity requires that organization admits variants; if every mutation
    leaves the core identical (or destroys it entirely is counted as a
    change), there is heritable structure to select on.
    """
    rng = np.random.default_rng(seed)
    core0 = raf_reactions(net)
    cat_rxns = [r for r, rxn in enumerate(net.reactions) if rxn.catalysts]
    if not core0 or not cat_rxns:
        return 0.0
    changed = 0
    valid = 0
    for _ in range(n_mut):
        reactions = [
            Reaction(dict(r.reactants), dict(r.products), list(r.catalysts),
                     r.rate, r.cat_rate)
            for r in net.reactions
        ]
        r = int(rng.choice(cat_rxns))
        reactions[r].catalysts = [int(rng.integers(net.n_species))]
        mutant = ReactionNetwork(net.species, reactions, set(net.food),
                                 dict(net.structure))
        if raf_reactions(mutant):
            valid += 1
            if raf_reactions(mutant) != core0:
                changed += 1
    return float(changed / valid) if valid else 0.0
