import numpy as np

from chemtrans.networks.core import Reaction, ReactionNetwork


def tiny_net() -> ReactionNetwork:
    # A + B -> C, catalyzed by C itself (minimal autocatalytic loop)
    return ReactionNetwork(
        species=["A", "B", "C"],
        reactions=[Reaction(reactants={0: 1, 1: 1}, products={2: 1}, catalysts=[2])],
        food={0, 1},
    )


def test_stoichiometry():
    S = tiny_net().stoichiometry()
    assert S.shape == (3, 1)
    assert list(S[:, 0]) == [-1, -1, 1]


def test_catalysis_matrix_is_autocatalytic():
    K = tiny_net().catalysis_matrix()
    assert K.shape == (3, 3)
    # C catalyzes its own production -> nonzero diagonal element
    assert K[2, 2] > 0
    # spectral radius > 0 for an autocatalytic network
    assert max(abs(np.linalg.eigvals(K))) > 0


def test_propensity_positive_and_catalyst_boost():
    net = tiny_net()
    x = np.array([1.0, 1.0, 0.0])
    a_no_cat = net.propensity(x)
    x = np.array([1.0, 1.0, 2.0])
    a_cat = net.propensity(x)
    assert a_cat[0] > a_no_cat[0] > 0
