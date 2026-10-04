"""Core data structure for chemical reaction networks.

A ReactionNetwork is the directed hypergraph G = (M, R, C) from the
conversation: molecules M, reactions R, catalytic relationships C.

Conventions
-----------
- Each reaction maps multisets of reactants -> products (dict species -> stoich).
- Catalysis is stored separately from stoichiometry: a catalyst accelerates a
  reaction without being consumed. This keeps the catalysis matrix K clean,
  which is what the spectral-radius threshold rho(K) is computed from.
- "Food" species are externally supplied (chemostat feed) and are the main
  lever for the energy-flux part of the control parameter lambda.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Reaction:
    """One reaction in the network.

    reactants / products: {species_index: stoichiometric coefficient}
    catalysts: species indices that accelerate this reaction (not consumed)
    rate: mass-action rate constant (uncatalyzed)
    cat_rate: catalytic rate constant (per unit catalyst, mass-action)
    k_sat: None = pure mass-action catalysis (unbounded — a closed
        autocatalytic core diverges under a chemostat). Finite value =
        Michaelis-Menten saturation: the catalyzed rate becomes
        cat_rate * prod_i[x_i/(k_sat+x_i)] * sum_c[x_c/(k_sat+x_c)],
        bounded by ~cat_rate, giving the organized phase a real steady
        state.

    Convention: rate == 0 with a nonempty catalysts list means the reaction
    REQUIRES a catalyst — it cannot fire spontaneously. The RAF/closure
    analysis uses this to decide which reactions are conditionally enabled.
    """

    reactants: dict[int, int]
    products: dict[int, int]
    catalysts: list[int] = field(default_factory=list)
    rate: float = 1.0
    cat_rate: float = 10.0
    k_sat: float | None = None


@dataclass
class ReactionNetwork:
    """A chemical reaction network with explicit catalysis."""

    species: list[str]
    reactions: list[Reaction]
    food: set[int] = field(default_factory=set)
    # optional structural info per species (e.g. polymer string) used by
    # observables.assembly to compute a real assembly index
    structure: dict[int, str] = field(default_factory=dict)
    # cached vectorized propensity arrays (False = use slow loop)
    _prop_cache: object = field(
        default=None, init=False, repr=False, compare=False)

    @property
    def n_species(self) -> int:
        return len(self.species)

    @property
    def n_reactions(self) -> int:
        return len(self.reactions)

    def stoichiometry(self) -> np.ndarray:
        """S[i, r] = net production of species i by reaction r."""
        S = np.zeros((self.n_species, self.n_reactions))
        for r, rxn in enumerate(self.reactions):
            for i, nu in rxn.products.items():
                S[i, r] += nu
            for i, nu in rxn.reactants.items():
                S[i, r] -= nu
        return S

    def catalysis_matrix(self) -> np.ndarray:
        """K[i, j] = catalytic strength of species i on the production of j.

        Built by attributing each catalyzed reaction's output to its
        catalysts. The spectral radius rho(K) is the candidate threshold
        variable: rho(K) < 1 catalytic activity dies out, rho(K) > 1 permits
        self-sustaining amplification.
        """
        K = np.zeros((self.n_species, self.n_species))
        for rxn in self.reactions:
            for cat in rxn.catalysts:
                for j, nu in rxn.products.items():
                    K[cat, j] += rxn.cat_rate * nu
        return K

    def catalysis_graph(self):
        """Directed graph of catalytic edges (for RAF/closure analysis)."""
        import networkx as nx

        g = nx.DiGraph()
        g.add_nodes_from(range(self.n_species))
        for rxn in self.reactions:
            for cat in rxn.catalysts:
                for j in rxn.products:
                    g.add_edge(cat, j, reaction=self.reactions.index(rxn))
        return g

    def propensity(self, x: np.ndarray) -> np.ndarray:
        """Propensity of each reaction at concentration vector x.

        Uncatalyzed term: rate * prod(x_i^nu_i).
        Catalyzed term: mass-action (k_sat None) or Michaelis-Menten
        saturating in both substrates and catalyst (k_sat finite).
        """
        x = np.maximum(x, 0.0)
        fa = self._fast_arrays()
        if fa is not None:
            R, rate, cat_rate, C, k_sat = fa
            if k_sat is None:
                xp = np.append(x, 1.0)
                mass = xp[R[:, 0]] * xp[R[:, 1]]
                cat = C @ x
            else:
                xs = x / (k_sat + x)
                xsp = np.append(xs, 1.0)
                mass = xsp[R[:, 0]] * xsp[R[:, 1]]
                cat = C @ xs
            return np.maximum(rate * mass + cat_rate * mass * cat, 0.0)
        a = np.zeros(self.n_reactions)
        for r, rxn in enumerate(self.reactions):
            if rxn.k_sat is None:
                mass = 1.0
                for i, nu in rxn.reactants.items():
                    mass *= x[i] ** nu
                cat = sum(x[c] for c in rxn.catalysts)
            else:
                # Michaelis-Menten: substrate AND catalyst saturate —
                # every rate is bounded, so steady states always exist
                mass = 1.0
                for i, nu in rxn.reactants.items():
                    mass *= (x[i] / (rxn.k_sat + x[i])) ** nu
                cat = sum(x[c] / (rxn.k_sat + x[c]) for c in rxn.catalysts)
            a[r] = rxn.rate * mass + rxn.cat_rate * mass * cat
        return np.maximum(a, 0.0)

    def _fast_arrays(self):
        """Vectorized propensity arrays, or None if any reaction needs
        the general path (>2 reactants, non-unit stoich, mixed k_sat).
        Cached on the instance; networks are treated as immutable."""
        if self._prop_cache is not None:
            return self._prop_cache or None
        k_sats = {r.k_sat for r in self.reactions}
        ok = len(k_sats) == 1 and all(
            len(r.reactants) <= 2
            and all(nu == 1 for nu in r.reactants.values())
            for r in self.reactions
        )
        if not ok:
            self._prop_cache = False
            return None
        from scipy import sparse

        n = self.n_species
        R = np.full((self.n_reactions, 2), n)  # sentinel -> appended 1.0
        C = sparse.dok_matrix((self.n_reactions, n))
        rate = np.empty(self.n_reactions)
        cat_rate = np.empty(self.n_reactions)
        for r, rxn in enumerate(self.reactions):
            for k, s in enumerate(rxn.reactants):
                R[r, k] = s
            for c in rxn.catalysts:
                C[r, c] += 1.0
            rate[r] = rxn.rate
            cat_rate[r] = rxn.cat_rate
        self._prop_cache = (R, rate, cat_rate, C.tocsr(), k_sats.pop())
        return self._prop_cache
