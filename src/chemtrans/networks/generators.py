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
    exclusive_catalysis: bool = False,
    k_sat: float | None = None,
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

    exclusive_catalysis: if True, species produced by catalyst-required
    reactions are given NO spontaneous production route — spontaneous and
    catalyst-accelerated reactions are resampled so their products avoid
    those species entirely. This creates genuine nucleation barriers (the
    RAF core cannot bootstrap from food alone) and makes dynamical
    bistability possible; see runs/NOTEBOOK.md.

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
        # disallow reactant->product overlap: A->A+B is a spontaneous
        # replicator creating matter ex nihilo and diverging under mass
        # action. Amplification belongs exclusively in the catalysis
        # channel (that's the thing being studied).
        if set(products.tolist()) & set(reactants.tolist()):
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
            k_sat=k_sat,
        )
        if rng.random() < p_cat:
            cat = int(rng.integers(n_species))
            if cat not in rxn.reactants:
                rxn.catalysts = [cat]
                if rng.random() < p_cat_required:
                    rxn.rate = 0.0
        reactions.append(rxn)

    if exclusive_catalysis:
        excl: set[int] = set()
        for rxn in reactions:
            if rxn.rate == 0.0 and rxn.catalysts:
                excl.update(rxn.products)
        excl -= food
        allowed = [i for i in range(n_species) if i not in excl]
        excl_list = sorted(excl)
        kept = []
        for rxn in reactions:
            if rxn.rate > 0.0 and set(rxn.products) & excl:
                if not allowed:
                    continue  # no spontaneous-eligible product exists; drop
                n_p = len(rxn.products)
                rxn.products = {
                    int(i): 1
                    for i in rng.choice(
                        allowed, size=n_p, replace=len(allowed) < n_p
                    )
                }
            kept.append(rxn)
        reactions = kept
        # Nucleation barrier also requires the CATALYSTS of required
        # reactions to be exclusive species — otherwise a spontaneous
        # route to a catalyst ignites the whole core anyway
        # (see runs/NOTEBOOK.md, M2/M3 findings).
        if excl_list:
            for rxn in reactions:
                if rxn.rate == 0.0 and rxn.catalysts:
                    rxn.catalysts = [
                        int(rng.choice(excl_list)) for _ in rxn.catalysts
                    ]

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
    alphabet_size: int = 3,
    max_length: int = 9,
    n_ligations: int = 150,
    p_cat: float = 0.3,
    p_info: float = 0.7,
    p_selfrep: float = 0.3,
    p_cat_required: float = 0.5,
    cat_boost: float = 20.0,
    cleavage_rate: float = 0.02,
    ligation_rate: float = 0.05,
    k_sat: float | None = None,
    food_len: int = 1,
    seed: int | None = None,
) -> ReactionNetwork:
    """Ligation/cleavage polymer chemistry with template catalysis.

    Species are strings over `alphabet_size` letters; structure[idx] holds
    the string, enabling true assembly-index computation.

    Reactions:
      - ligation  a + b -> ab   (rate ligation_rate, possibly template-
        catalyzed; rate=0 for the p_cat_required fraction)
      - cleavage  ab -> a + b   (spontaneous hydrolysis, cleavage_rate)

    Catalysis: each ligation gets a catalyst with prob p_cat. Of those,
    fraction p_info are informational — the product itself with prob
    p_selfrep (minimal self-template replicator; this is the only motif
    that closes small RAF cores) or a species sharing the product's
    sequence (super-/substring). The rest are random species. p_cat is
    thus the connectivity knob, p_info the information-channel strength,
    p_selfrep the replicator density — they decouple, which is what lets
    H vary independently of Phi in the discrimination test.

    Strings of length <= food_len are food (supplied by the chemostat);
    food_len=1 gives monomers only, food_len=2 a "rich broth" of short
    oligomers that makes small RAFs much more accessible.
    """
    rng = np.random.default_rng(seed)
    letters = [chr(65 + i) for i in range(alphabet_size)]
    species_set = list(letters)
    index = {s: i for i, s in enumerate(species_set)}
    structure: dict[int, str] = {i: s for s, i in index.items()}

    def add_species(s: str) -> int:
        if s not in index:
            index[s] = len(species_set)
            species_set.append(s)
            structure[index[s]] = s
        return index[s]

    # pass 1: grow the species pool via random ligations. Pair sampling
    # is biased toward shorter strings (~1/len): short oligomers are more
    # abundant/mobile, and this is what makes monomer-level self-
    # replicators (a+b -> ab cat ab with a,b in food) common enough to
    # form RAF cores.
    lig_pairs: list[tuple[str, str, str]] = []
    seen: set[tuple[str, str]] = set()
    attempts = 0
    while len(lig_pairs) < n_ligations and attempts < n_ligations * 200:
        attempts += 1
        w = 1.0 / np.array([len(s) for s in species_set])
        a, b = rng.choice(species_set, size=2, p=w / w.sum())
        prod = a + b
        if len(prod) > max_length or (a, b) in seen:
            continue
        seen.add((a, b))
        lig_pairs.append((a, b, prod))
        add_species(prod)

    # pass 2: close the pool under cleavage — every substring of a
    # sampled product is itself a species (hydrolysis produces it).
    for s in list(species_set):
        for k in range(1, len(s)):
            for st in range(len(s) - k + 1):
                add_species(s[st:st + k])

    # pass 3: for every product string, EVERY split a+b is a ligation
    # route (and its reverse a cleavage). Catalysis is per product:
    # template matching is a property of what is being assembled, so all
    # routes to a self-templated product share the self-rep catalyst.
    reactions: list[Reaction] = []
    for prod in species_set:
        if len(prod) < 2:
            continue
        routes = [(prod[:k], prod[k:]) for k in range(1, len(prod))]
        # catalyst identity for this product (shared by all its routes)
        catalyst: int | None = None
        if rng.random() < p_cat:
            if rng.random() < p_info:
                # informational: self-template (cat == prod, the minimal
                # replicator RAF) with prob p_selfrep, else a species
                # sharing the product's sequence (super-/substring)
                if rng.random() < p_selfrep:
                    catalyst = index[prod]
                else:
                    kin = [s for s in species_set
                           if s != prod and (prod in s or s in prod)]
                    catalyst = index[kin[int(rng.integers(len(kin)))]] \
                        if kin else index[prod]
            else:
                catalyst = index[
                    species_set[int(rng.integers(len(species_set)))]]
        required = catalyst is not None and rng.random() < p_cat_required
        for a, b in routes:
            base = ligation_rate * float(np.exp(rng.normal(0.0, 0.5)))
            lig = Reaction(
                reactants={index[a]: 1, index[b]: 1},
                products={index[prod]: 1},
                rate=0.0 if required else base,
                cat_rate=base * cat_boost,
                k_sat=k_sat,
            )
            if catalyst is not None and catalyst not in lig.reactants:
                lig.catalysts = [catalyst]
            reactions.append(lig)
            reactions.append(Reaction(
                reactants={index[prod]: 1},
                products={index[a]: 1, index[b]: 1},
                rate=cleavage_rate * float(np.exp(rng.normal(0.0, 0.5))),
                k_sat=k_sat,
            ))

    return ReactionNetwork(
        species=species_set,
        reactions=reactions,
        food={i for i, s in structure.items() if len(s) <= food_len},
        structure=structure,
    )
