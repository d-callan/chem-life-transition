"""Assembly-index observables C(lambda).

Two variants:

- path_assembly: cheap proxy — shortest reaction-path length from food
  species to each molecule in THIS network. Well-defined for any
  ReactionNetwork.
- true assembly index (assembly theory): minimum joining steps over
  molecular substructure. Only defined for generators that attach real
  structure (e.g. polymer_ligation strings).

The key observable is the DISTRIBUTION P(A), not the mean: the
conversation predicts the high-assembly tail gets populated and becomes
persistent only above the transition.
"""

import numpy as np

from ..networks.core import ReactionNetwork


def path_assembly_index(net: ReactionNetwork) -> np.ndarray:
    """Shortest-path assembly index for every species (inf if unreachable)."""
    raise NotImplementedError


def molecular_assembly_index(net: ReactionNetwork) -> np.ndarray:
    """True assembly index from species structure (requires structure field)."""
    raise NotImplementedError


def assembly_distribution(
    net: ReactionNetwork, abundance: np.ndarray, true_index: bool = False
) -> tuple[np.ndarray, np.ndarray]:
    """Abundance-weighted P(A). Returns (bin_edges, probabilities)."""
    raise NotImplementedError
