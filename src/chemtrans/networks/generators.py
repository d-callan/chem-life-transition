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

from .core import ReactionNetwork


def random_crn(
    n_species: int,
    n_reactions: int,
    p_cat: float,
    n_food: int = 5,
    max_reactants: int = 2,
    max_products: int = 2,
    seed: int | None = None,
) -> ReactionNetwork:
    """Random CRN with tunable catalytic connectivity.

    p_cat is the primary component of the control parameter lambda:
    raising it raises the spectral radius of the catalysis matrix K.
    """
    raise NotImplementedError


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
