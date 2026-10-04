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


def jacobian_spectral_abscissa(
    net: ReactionNetwork,
    x: np.ndarray,
    feed_rate: float,
    dilution_rate: float,
    h: float = 1e-5,
    alive_threshold: float = 1e-6,
) -> float:
    """Dominant eigenvalue alpha of the RESIDENT chemostat Jacobian at x.

    The Jacobian is restricted to species with x > alive_threshold.
    This matters: the positivity clamp (x >= 0) makes naive finite
    differences at extinct species produce spurious unstable directions
    — a state can be a stable fixed point on its support while the full
    Jacobian reports alpha > 0. alpha < 0: locally stable attractor;
    alpha -> 0 approaching a bifurcation (critical slowing down).
    """
    from scipy import sparse

    S = sparse.csr_matrix(net.stoichiometry())
    feed = np.zeros(net.n_species)
    feed[list(net.food)] = feed_rate
    x = np.maximum(np.asarray(x, dtype=float), 0.0)
    alive = np.where(x > alive_threshold)[0]
    if len(alive) == 0:
        return 0.0

    def f(z):
        z = np.maximum(z, 0.0)
        return S @ net.propensity(z) + feed - dilution_rate * z

    J = np.empty((len(alive), len(alive)))
    for k, j in enumerate(alive):
        e = np.zeros(net.n_species)
        e[j] = h * max(1.0, x[j])
        J[:, k] = (f(x + e) - f(x - e))[alive] / (2.0 * e[j])
    return float(np.linalg.eigvals(J).real.max())


def invasion_pressure(
    net: ReactionNetwork,
    x: np.ndarray,
    feed_rate: float,
    dilution_rate: float,
    alive_threshold: float = 1e-6,
) -> float:
    """Max net growth rate among extinct species at state x.

    > 0 means the current community could be invaded by an absent species
    IF it could nucleate; combined with a catalyst-required barrier this
    is the signature of bistability (organized attractor coexists with a
    dead attractor). <= 0: the community is invasion-proof.
    """
    from scipy import sparse

    S = sparse.csr_matrix(net.stoichiometry())
    feed = np.zeros(net.n_species)
    feed[list(net.food)] = feed_rate
    x = np.maximum(np.asarray(x, dtype=float), 0.0)
    dead = np.where(x <= alive_threshold)[0]
    if len(dead) == 0:
        return 0.0
    drift = S @ net.propensity(x) + feed - dilution_rate * x
    return float(drift[dead].max())
