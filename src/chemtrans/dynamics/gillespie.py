"""Stochastic simulation (Gillespie SSA) for small networks.

Deterministic ODEs can't capture extinction: an autocatalytic core at
concentration 1e-9 is "alive" in the ODE but dead in reality. SSA is needed
to measure genuine persistence vs. numerical persistence, and to see
demographic noise near the transition.
"""

import numpy as np

from ..networks.core import ReactionNetwork


def simulate(
    net: ReactionNetwork,
    t_max: float,
    feed_rate: float = 1.0,
    dilution_rate: float = 0.1,
    volume: float = 1.0,
    x0: np.ndarray | None = None,
    max_events: int = 10_000_000,
) -> tuple[np.ndarray, np.ndarray]:
    """Run SSA with food influx and dilution as Poisson/death processes.

    Returns (t, X) event history; X in molecule counts.
    """
    raise NotImplementedError
