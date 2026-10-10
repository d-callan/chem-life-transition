"""Channel/network invariants: stoichiometry, thermodynamic consistency,
conservation, hashing."""

import numpy as np
import pytest

from chemtrans.networks.core import Channel, ReactionNetwork


def mini_net():
    """A + B <=> C catalyzed by C, food={A,B}."""
    ch = Channel(reactants=((0, 1), (1, 1)), products=((2, 1),),
                 catalysts=(2,), k_f=2.0, k_r=0.5, synth=True,
                 mechanism="elongation")
    return ReactionNetwork(
        ("A", "B", "C"), [ch], (0, 1),
        comp=np.array([[1.0, 1.0, 2.0]]),
        g=np.zeros(3), food_inlet=np.array([1.0, 1.0, 0.0]),
        kappa=10.0, eps=1e-3)


def test_stoichiometry_signs():
    net = mini_net()
    S = net.stoichiometry()
    assert S.shape == (3, 1)
    assert list(S[:, 0]) == [-1.0, -1.0, 1.0]


def test_reverse_identity_and_conservation():
    net = mini_net()
    S = net.stoichiometry()
    # forward and reverse of the same channel are the same column
    assert np.allclose(net.comp @ S, 0.0)  # letter mass conserved


def test_local_detailed_balance():
    """k_f/k_r must equal exp(-dg) for every channel."""
    g = np.array([0.4, -0.1])
    k_eq = np.exp(-(g[1] - g[0]))          # nu = (-1, +1)
    ch = Channel(reactants=((0, 1),), products=((1, 1),),
                 catalysts=(), k_f=3.0, k_r=3.0 / k_eq, synth=True)
    net = ReactionNetwork(
        ("A", "B"), [ch], (0,), comp=np.ones((1, 2)),
        g=g, food_inlet=np.array([1.0, 0.0]))
    assert net.check_detailed_balance()["ok"]


def test_inconsistent_g_detected():
    """A channel whose k_f/k_r contradicts exp(-dg) must be caught."""
    ch = Channel(reactants=((0, 1),), products=((1, 1),),
                 catalysts=(), k_f=5.0, k_r=1.0, synth=True)
    net = ReactionNetwork(
        ("A", "B"), [ch], (0,), comp=np.ones((1, 2)),
        g=np.array([0.0, 0.0]), food_inlet=np.array([1.0, 0.0]))
    assert not net.check_detailed_balance()["ok"]


def test_undeclared_irreversible_fails_detailed_balance():
    """k_r=0 with finite g cannot satisfy detailed balance — it must
    fail loudly, not be silently skipped."""
    ch = Channel(reactants=((0, 1),), products=((1, 1),),
                 catalysts=(), k_f=1.0, k_r=0.0, synth=True)
    net = ReactionNetwork(
        ("A", "B"), [ch], (0,), comp=np.ones((1, 2)),
        g=np.array([0.0, 0.0]), food_inlet=np.array([1.0, 0.0]))
    db = net.check_detailed_balance()
    assert not db["ok"]
    row = db["channels"][0]
    assert row["error"] == "irreversible_channel"
    assert row["declared_driven"] is False
    assert "rel_err" not in row


def test_declared_driven_irreversible_reported_not_failed():
    """A deliberately one-way driven limit (declared driven=True)
    passes validation but is still reported as irreversible — never
    silently treated as satisfying detailed balance."""
    ch = Channel(reactants=((0, 1),), products=((1, 1),),
                 catalysts=(), k_f=1.0, k_r=0.0, synth=True, driven=True)
    net = ReactionNetwork(
        ("A", "B"), [ch], (0,), comp=np.ones((1, 2)),
        g=np.array([0.0, 0.0]), food_inlet=np.array([1.0, 0.0]))
    db = net.check_detailed_balance()
    assert db["ok"]
    row = db["channels"][0]
    assert row["error"] == "irreversible_channel"
    assert row["declared_driven"] is True


def test_synth_dir_validated():
    ch = Channel(reactants=((0, 1),), products=((1, 1),),
                 synth_dir="sideways")
    with pytest.raises(ValueError):
        ReactionNetwork(("A", "B"), [ch], (0,),
                        comp=np.ones((1, 2)),
                        g=np.zeros(2),
                        food_inlet=np.array([1.0, 0.0]))


def test_channel_new_fields_serialize():
    """synth_dir and driven must round-trip through to_dict/from_dict."""
    ch = Channel(reactants=((0, 1),), products=((1, 1),), catalysts=(),
                 k_f=1.0, k_r=0.5, synth=True, synth_dir="rev",
                 driven=True)
    net = ReactionNetwork(
        ("A", "B"), [ch], (0,), comp=np.ones((1, 2)),
        g=np.array([0.0, -0.3]), food_inlet=np.array([1.0, 0.0]))
    net2 = ReactionNetwork.from_dict(net.to_dict())
    assert net2.channels[0].synth_dir == "rev"
    assert net2.channels[0].driven is True
    assert net2.hash == net.hash
    # without_catalysis must preserve them too
    net3 = net.without_catalysis()
    assert net3.channels[0].synth_dir == "rev"
    assert net3.channels[0].driven is True


def test_catalysis_matrix():
    net = mini_net()
    K = net.catalysis_matrix()
    assert K.shape == (1, 3)
    assert K[0, 2] == 1.0


def test_net_flux_factorizes():
    net = mini_net()
    x = np.array([0.4, 0.4, 0.3])
    eps = 0.01
    f, b = net.rates(x, eps)
    # catalyzed boost applied to both directions
    boost = eps + net.kappa * x[2]
    assert np.isclose(f[0], 2.0 * 0.4 * 0.4 * boost)
    assert np.isclose(b[0], 0.5 * 0.3 * boost)


def test_hash_deterministic_and_sensitive():
    a = mini_net()
    b = mini_net()
    assert a.hash == b.hash
    ch = Channel(reactants=((0, 1), (1, 1)), products=((2, 1),),
                 catalysts=(2,), k_f=2.5, k_r=0.5, synth=True)
    c = ReactionNetwork(("A", "B", "C"), [ch], (0, 1),
                        comp=np.array([[1., 1., 2.]]), g=np.zeros(3),
                        food_inlet=np.array([1., 1., 0.]))
    assert c.hash != a.hash
