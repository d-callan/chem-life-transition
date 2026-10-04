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


def _eigvals_sorted(net: ReactionNetwork) -> np.ndarray:
    K = net.catalysis_matrix()
    if K.size == 0:
        return np.array([0.0])
    ev = np.abs(np.linalg.eigvals(K))
    return np.sort(ev)[::-1]


def spectral_radius(net: ReactionNetwork) -> float:
    """Largest-magnitude eigenvalue of the catalysis matrix K."""
    return float(_eigvals_sorted(net)[0])


def spectral_gap(net: ReactionNetwork) -> float:
    """Gap between leading and next eigenvalue — sharpens near transitions."""
    ev = _eigvals_sorted(net)
    return float(ev[0] - ev[1]) if len(ev) > 1 else 0.0
