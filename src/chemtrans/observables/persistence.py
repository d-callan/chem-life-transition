"""Persistence and self-maintenance S(lambda), heredity/evolvability H/E.

persistence_time
    How long the organized state survives dilution shocks / noise. The
    dynamical definition of self-maintenance: does the RAF core reform
    after perturbation?
reproductive_number
    R0 for the organized state: expected copies of the autocatalytic core
    produced per decay — the Model-A T4 / evolvability threshold probe.
evolvability_proxy
    Variance of offsprings' core composition under mutation of network
    edges (does heritable variation exist at all?).
"""

import numpy as np

from ..networks.core import ReactionNetwork


def persistence_time(
    net: ReactionNetwork, shock_strength: float = 0.9, n_trials: int = 20
) -> float:
    """Median time for the RAF core to reform (or die) after a dilution shock."""
    raise NotImplementedError


def reproductive_number(net: ReactionNetwork, x: np.ndarray) -> float:
    """R0 of the organized state at steady state."""
    raise NotImplementedError
