"""Organizational closure Phi — the primary candidate order parameter.

Phi = (chemical flux inside self-maintaining organization) / (total flux)

Two pieces:
1. Structural: find the self-sustaining core of the network — the RAF set.
   We use the strict definition: a set of catalyzed reactions where every
   reaction's reactants are food-or-produced-within-the-set AND every
   reaction's catalyst is produced within the set (not externally
   supplied — an imported catalyst isn't self-maintaining).
   Computed by iterative pruning to a fixed point.
2. Dynamical: at steady state, what fraction of total reaction flux runs
   through that core?

Also the autocatalytic-activity observable A(lambda) for the
transition-discrimination experiment.
"""

import numpy as np

from ..networks.core import ReactionNetwork


def producible_closure(net: ReactionNetwork) -> tuple[set[int], set[int]]:
    """Species producible from food, and the reactions that run.

    A reaction is usable when all reactants are available AND (it is
    spontaneous (rate>0) or at least one catalyst is available).
    Grow pass to fixed point. Includes non-catalytic chemistry — this is
    "what can run", not "what is self-maintaining".
    """
    available = set(net.food)
    used: set[int] = set()
    changed = True
    while changed:
        changed = False
        for r, rxn in enumerate(net.reactions):
            if r in used or not set(rxn.reactants) <= available:
                continue
            if rxn.rate == 0.0 and not set(rxn.catalysts) & available:
                continue
            used.add(r)
            available.update(rxn.products)
            changed = True
    return available, used


def raf_reactions(net: ReactionNetwork) -> set[int]:
    """Indices of reactions in the maximal RAF (self-maintaining) core.

    Candidate set = catalyzed reactions only. Iterate: drop any reaction
    whose reactants are not food-union-core-products, or whose catalysts
    are not produced by the core itself. Fixed point is the RAF core.
    """
    core = {r for r, rxn in enumerate(net.reactions) if rxn.catalysts}
    while True:
        produced = set()
        for r in core:
            produced.update(net.reactions[r].products)
        avail = set(net.food) | produced
        pruned = {
            r
            for r in core
            if set(net.reactions[r].reactants) <= avail
            and set(net.reactions[r].catalysts) <= produced
        }
        if pruned == core:
            return core
        core = pruned


def raf_core(net: ReactionNetwork) -> set[int]:
    """Species indices belonging to the maximal reflexively autocatalytic set."""
    species: set[int] = set()
    for r in raf_reactions(net):
        species.update(net.reactions[r].products)
        species.update(net.reactions[r].reactants)
        species.update(net.reactions[r].catalysts)
    return species - net.food


def organizational_closure(net: ReactionNetwork, flux: np.ndarray) -> float:
    """Phi: fraction of reaction flux inside the RAF core.

    flux: per-reaction rates at the state of interest, i.e.
    net.propensity(x_ss).
    """
    total = float(np.sum(flux))
    if total <= 0.0:
        return 0.0
    core = raf_reactions(net)
    return float(np.sum(flux[list(core)])) / total if core else 0.0


def autocatalytic_activity(net: ReactionNetwork, x: np.ndarray) -> float:
    """A(lambda): fraction of production flux attributable to catalysis."""
    cat_flux = 0.0
    total_flux = 0.0
    x = np.maximum(np.asarray(x, dtype=float), 0.0)
    for rxn in net.reactions:
        mass = 1.0
        for i, nu in rxn.reactants.items():
            mass *= x[i] ** nu
        cat = rxn.cat_rate * mass * sum(x[c] for c in rxn.catalysts)
        cat_flux += cat
        total_flux += rxn.rate * mass + cat
    return float(cat_flux / total_flux) if total_flux > 0.0 else 0.0
