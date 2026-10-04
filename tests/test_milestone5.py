"""M5: Gillespie SSA correctness + 2-D manifold machinery."""

import numpy as np

from chemtrans.dynamics.gillespie import simulate as ssa
from chemtrans.networks.core import Reaction, ReactionNetwork


def test_ssa_food_only_poisson():
    """Food species with no reactions: N ~ Poisson(feed*V/dilution)."""
    net = ReactionNetwork(species=["F", "X"], reactions=[], food={0})
    _, X = ssa(net, t_max=200, feed_rate=1.0, dilution_rate=0.1,
               volume=100, seed=0)
    mean = X[10:, 0].mean()
    assert 700 < mean < 1300      # expected 1000, loose band
    assert X[-1, 1] == 0          # non-food species never appears


def test_ssa_catalyst_required_never_nucleates():
    """rate=0 reaction without catalyst present cannot fire in SSA —
    the chicken-and-egg barrier is exact at finite copy number."""
    net = ReactionNetwork(
        species=["F", "A"],
        reactions=[Reaction(reactants={0: 1}, products={1: 1},
                            rate=0.0, cat_rate=5.0, catalysts=[1])],
        food={0})
    _, X = ssa(net, t_max=50, feed_rate=1.0, dilution_rate=0.1,
               volume=50, seed=1)
    assert X[:, 1].max() == 0


def test_ssa_nucleation_and_growth():
    """A self-replicator channel with small spontaneous rate nucleates
    and reaches nonzero copy number (extinction-vs-ODE regime)."""
    net = ReactionNetwork(
        species=["F", "A"],
        reactions=[Reaction(reactants={0: 1}, products={1: 1},
                            rate=0.05, cat_rate=3.0, catalysts=[1],
                            k_sat=1.0)],
        food={0})
    _, X = ssa(net, t_max=100, feed_rate=1.0, dilution_rate=0.1,
               volume=50, seed=2)
    assert X[:, 1].max() > 50


def test_nu2_reactant_fast_path_matches_slow():
    """A+A->AA (nu=2) must consume two copies; vectorized propensity
    must equal the slow loop on mixed-stoichiometry networks."""
    net = ReactionNetwork(
        species=["A", "AA"],
        reactions=[Reaction(reactants={0: 2}, products={1: 1},
                            rate=1.0, cat_rate=0.0)],
        food={0})
    x = np.array([3.0, 1.0])
    fast = net.propensity(x)
    net._prop_cache = False
    slow = net.propensity(x)
    assert np.allclose(fast, slow)
    assert fast[0] == 9.0        # 1.0 * 3^2
