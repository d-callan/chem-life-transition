"""Tunable chemical reaction network generators.

These are the "fiddle knobs". Every generator must return a ReactionNetwork
AND record its full parameter set in the run metadata — results are
sensitive to generator assumptions, so provenance is first-class.

Planned generators
------------------
random_crn
    Erdos-Renyi-flavored: n_species, n_reactions sampled with random
    reactant/product sets; each reaction gets a random catalyst with prob
    p_cat. Simplest baseline.
raf_biased
    Generates networks with a tunable fraction of reactions arranged into
    reflexively autocatalytic sets (Kauffman-style), to probe whether
    seeded closure changes the transition structure.
polymer_ligation
    Species are strings over a small alphabet; reactions are
    ligation (A + B -> AB) and cleavage (AB -> A + B), with optional
    template catalysis. Gives a REAL assembly index for free and is the
    most physically motivated model for the information-chemistry link.
"""

import numpy as np

from .core import Reaction, ReactionNetwork


def random_crn(
    n_species: int,
    n_reactions: int,
    p_cat: float,
    n_food: int = 5,
    max_reactants: int = 2,
    max_products: int = 2,
    p_cat_required: float = 0.5,
    cat_boost: float = 20.0,
    seed: int | None = None,
) -> ReactionNetwork:
    """Random CRN with tunable catalytic connectivity.

    p_cat is the primary component of the control parameter lambda:
    raising it raises the spectral radius of the catalysis matrix K.

    p_cat_required: of the reactions that get a catalyst, this fraction is
    catalyst-REQUIRED (rate=0, cannot fire spontaneously). The rest are
    catalyst-accelerated. This knob interpolates between "catalysis is a
    speedup" (0) and "catalysis is a precondition for new chemistry" (1) —
    the latter is the physically interesting regime for RAF formation.

    Food species are indices 0..n_food-1 and are supplied by the chemostat.
    Rates are lognormal(sigma=0.5) around 1.0 so the network is
    heterogeneous rather than a mean-field caricature.
    """
    rng = np.random.default_rng(seed)
    species = [f"m{i}" for i in range(n_species)]
    food = set(range(n_food))

    reactions: list[Reaction] = []
    seen: set[tuple] = set()
    attempts = 0
    while len(reactions) < n_reactions and attempts < n_reactions * 100:
        attempts += 1
        n_r = int(rng.integers(1, max_reactants + 1))
        n_p = int(rng.integers(1, max_products + 1))
        reactants = rng.choice(n_species, size=n_r, replace=False)
        products = rng.choice(n_species, size=n_p, replace=False)
        if set(products.tolist()) <= set(reactants.tolist()):
            continue
        key = (tuple(sorted(reactants.tolist())), tuple(sorted(products.tolist())))
        if key in seen:
            continue
        seen.add(key)

        base_rate = float(np.exp(rng.normal(0.0, 0.5)))
        rxn = Reaction(
            reactants={int(i): 1 for i in reactants},
            products={int(i): 1 for i in products},
            rate=base_rate,
            cat_rate=base_rate * cat_boost,
        )
        if rng.random() < p_cat:
            cat = int(rng.integers(n_species))
            if cat not in rxn.reactants:
                rxn.catalysts = [cat]
                if rng.random() < p_cat_required:
                    rxn.rate = 0.0
        reactions.append(rxn)

    return ReactionNetwork(species=species, reactions=reactions, food=food)


def raf_biased(
    n_species: int,
    n_reactions: int,
    closure_fraction: float,
    seed: int | None = None,
) -> ReactionNetwork:
    """CRN biased toward containing reflexively autocatalytic (RAF) sets."""
    raise NotImplementedError


def polymer_ligation(
    alphabet_size: int,
    max_length: int,
    p_template: float,
    n_food_monomers: int = 2,
    seed: int | None = None,
) -> ReactionNetwork:
    """Ligation/cleavage polymer chemistry with template catalysis.

    Species get structure[species] = string, enabling true assembly-index
    computation in observables.assembly.
    """
    raise NotImplementedError
