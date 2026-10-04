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
    cat_rate: additional rate constant per unit catalyst concentration

    Convention: rate == 0 with a nonempty catalysts list means the reaction
    REQUIRES a catalyst — it cannot fire spontaneously. The RAF/closure
    analysis uses this to decide which reactions are conditionally enabled.
    """

    reactants: dict[int, int]
    products: dict[int, int]
    catalysts: list[int] = field(default_factory=list)
    rate: float = 1.0
    cat_rate: float = 10.0


@dataclass
class ReactionNetwork:
    """A chemical reaction network with explicit catalysis."""

    species: list[str]
    reactions: list[Reaction]
    food: set[int] = field(default_factory=set)
    # optional structural info per species (e.g. polymer string) used by
    # observables.assembly to compute a real assembly index
    structure: dict[int, str] = field(default_factory=dict)

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
        """Mass-action propensity of each reaction at concentration vector x.

        a_r = rate * prod(x_i^nu_i) + cat_rate * prod(x_i^nu_i) * sum(x_cat)
        """
        a = np.zeros(self.n_reactions)
        for r, rxn in enumerate(self.reactions):
            mass = 1.0
            for i, nu in rxn.reactants.items():
                mass *= x[i] ** nu
            cat_conc = sum(x[c] for c in rxn.catalysts)
            a[r] = mass * (rxn.rate + rxn.cat_rate * cat_conc)
        return np.maximum(a, 0.0)
