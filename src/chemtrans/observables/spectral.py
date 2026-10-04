"""Spectral properties of the catalysis matrix K.

rho(K) (spectral radius) is the static threshold candidate:
  rho(K) < 1 -> catalytic perturbations die out (subcritical)
  rho(K) > 1 -> self-sustaining amplification possible (supercritical)

Comparing rho(K)=1 against the *dynamically observed* critical point
lambda_c is itself a useful test of how well linear theory predicts the
nonlinear transition.
"""

import numpy as np

from ..networks.core import ReactionNetwork


def spectral_radius(net: ReactionNetwork) -> float:
    """Largest-magnitude eigenvalue of the catalysis matrix K."""
    raise NotImplementedError


def spectral_gap(net: ReactionNetwork) -> float:
    """Gap between leading and next eigenvalue — sharpens near transitions."""
    raise NotImplementedError
