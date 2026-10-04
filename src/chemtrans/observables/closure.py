"""Organizational closure Phi — the primary candidate order parameter.

Phi = (chemical flux inside self-maintaining organization) / (total flux)

Two pieces:
1. Structural: find the self-sustaining core of the network — reactions
   whose reactants and catalysts are all produced within the set
   (RAF-set detection; networkx over the catalysis + stoichiometry graph).
2. Dynamical: at steady state, what fraction of total reaction flux runs
   through that core?

Also the autocatalytic-activity observable A(lambda) for the
transition-discrimination experiment.
"""

import numpy as np

from ..networks.core import ReactionNetwork


def raf_core(net: ReactionNetwork) -> set[int]:
    """Species indices belonging to the maximal reflexively autocatalytic set."""
    raise NotImplementedError


def organizational_closure(net: ReactionNetwork, flux: np.ndarray) -> float:
    """Phi: fraction of steady-state reaction flux inside the RAF core.

    flux: per-reaction steady-state rates from dynamics.simulate.
    """
    raise NotImplementedError


def autocatalytic_activity(net: ReactionNetwork, x: np.ndarray) -> float:
    """A(lambda): fraction of production flux attributable to catalysis."""
    raise NotImplementedError
