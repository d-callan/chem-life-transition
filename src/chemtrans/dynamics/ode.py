"""Deterministic mass-action dynamics under a chemostat.

Model: dx/dt = S @ a(x) + feed - dilution * x

- feed: constant influx of food species (energy/matter flux knob)
- dilution: uniform outflow; the "mortality" pressure that organizational
  closure must overcome to persist

This is the workhorse for lambda sweeps: cheap enough to run thousands of
networks per sweep point.
"""

import numpy as np

from ..networks.core import ReactionNetwork


def simulate(
    net: ReactionNetwork,
    t_span: tuple[float, float],
    feed_rate: float = 1.0,
    dilution_rate: float = 0.1,
    x0: np.ndarray | None = None,
    n_points: int = 2000,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate the chemostat ODE.

    Returns (t, X) with X shape (n_points, n_species).
    """
    raise NotImplementedError
