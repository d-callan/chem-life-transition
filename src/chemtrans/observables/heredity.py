"""Heritable-information observable H(lambda).

Only meaningful for generators with structured species (polymer_ligation):
a catalyst that SHARES sequence with its product — product ⊆ catalyst
(template stores the motif), catalyst ⊆ product (product copies the
motif), or equal (self-replication) — is informational. That is the
first place where "information" is causally relevant to chemistry rather
than just being complicated molecules.

Two measures:

template_fidelity
    Fraction of catalysis edges where product ⊆ catalyst sequence —
    the topological existence of template relations.
heredity_index (H)
    Flux-weighted version at a given state: fraction of catalyzed
    production flux that copies catalyst sequence information. This is
    the observable intended to decouple from RAF core size — a big
    structureless RAF gives high Phi but low H.
"""

import numpy as np

from ..networks.core import ReactionNetwork


def _is_template_match(net: ReactionNetwork, cat: int, prod: int) -> bool:
    """Catalyst and product share sequence: product ⊆ catalyst (product is
    a piece of the template) or catalyst ⊆ product (the product copies a
    catalyst motif). Equality is the self-replicating template."""
    c, p = net.structure.get(cat), net.structure.get(prod)
    return bool(c and p and (p in c or c in p))


def template_fidelity(net: ReactionNetwork) -> float:
    """Fraction of catalyst edges that are template matches (topological)."""
    hits = tot = 0
    for rxn in net.reactions:
        for c in rxn.catalysts:
            tot += 1
            if any(_is_template_match(net, c, j) for j in rxn.products):
                hits += 1
    return hits / tot if tot else 0.0


def heredity_index(net: ReactionNetwork, x: np.ndarray) -> float:
    """H(lambda): catalyzed flux through template-matching reactions /
    total catalyzed flux, at state x."""
    x = np.maximum(np.asarray(x, dtype=float), 0.0)
    info = tot = 0.0
    for rxn in net.reactions:
        if rxn.k_sat is None:
            mass = 1.0
            for i, nu in rxn.reactants.items():
                mass *= x[i] ** nu
            cat = rxn.cat_rate * mass * sum(x[c] for c in rxn.catalysts)
        else:
            mass = 1.0
            for i, nu in rxn.reactants.items():
                mass *= (x[i] / (rxn.k_sat + x[i])) ** nu
            cat = (rxn.cat_rate * mass
                   * sum(x[c] / (rxn.k_sat + x[c]) for c in rxn.catalysts))
        if cat <= 0:
            continue
        tot += cat
        if any(_is_template_match(net, c, j)
               for c in rxn.catalysts for j in rxn.products):
            info += cat
    return info / tot if tot > 0 else 0.0
