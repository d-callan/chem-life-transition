"""Generator invariants: mass conservation, thermodynamic consistency,
non-food catalysts, absorbing limit, deterministic hashing."""

import numpy as np

from chemtrans.networks.generators import polymer_network


def test_letter_conservation_all_channels():
    net = polymer_network(alphabet_size=2, max_len=3, seed=0)
    rep = net.check_letter_conservation()
    assert rep["ok"], rep["violations"]


def test_detailed_balance_by_construction():
    for dm in (0.0, 2.0, 5.0):
        net = polymer_network(alphabet_size=2, max_len=3, seed=1,
                              delta_mu=dm)
        assert net.check_detailed_balance()["ok"], dm


def test_catalysts_are_nonfood_polymers():
    net = polymer_network(alphabet_size=2, max_len=3, seed=2, p_cat=1.0,
                          min_cat_len=2)
    for ch in net.channels:
        for c in ch.catalysts:
            assert c not in net.food
            assert "*" not in net.species[c]
            assert len(net.species[c]) >= 2


def test_p_cat_zero_means_no_catalysis():
    net = polymer_network(alphabet_size=2, max_len=3, seed=3, p_cat=0.0)
    assert all(not ch.catalysts for ch in net.channels)


def test_p_selfrep_assigns_own_product():
    net = polymer_network(alphabet_size=2, max_len=3, seed=4, p_cat=1.0,
                          p_selfrep=1.0, min_cat_len=2)
    selfrep = [ch for ch in net.channels
               if ch.mechanism == "elongation" and ch.catalysts]
    assert selfrep
    for ch in selfrep:
        prod = {i for i, _ in ch.products}
        assert set(ch.catalysts) <= prod


def test_feed_mass_conserving_across_chi():
    a = polymer_network(alphabet_size=2, max_len=2, chi=0.0, seed=5)
    b = polymer_network(alphabet_size=2, max_len=2, chi=1.0, seed=5)
    assert np.isclose(a.feed_letters(), b.feed_letters())


def test_eps_zero_absorbing():
    """eps=0 + non-food catalysts => the zero-polymer state is absorbing
    for polymer synthesis (decay channels may still run)."""
    net = polymer_network(alphabet_size=2, max_len=3, seed=6, eps=0.0,
                          p_cat=0.5)
    x = net.food_inlet.copy()
    f, b = net.rates(x, eps=0.0)
    for r, ch in enumerate(net.channels):
        if ch.synth:
            assert f[r] == 0.0 and b[r] == 0.0


def test_deterministic_hash():
    a = polymer_network(alphabet_size=2, max_len=3, seed=7)
    b = polymer_network(alphabet_size=2, max_len=3, seed=7)
    c = polymer_network(alphabet_size=2, max_len=3, seed=8)
    assert a.hash == b.hash
    assert a.hash != c.hash


def test_roundtrip_serialization():
    from chemtrans.networks.core import ReactionNetwork
    net = polymer_network(alphabet_size=2, max_len=3, seed=9)
    back = ReactionNetwork.from_dict(net.to_dict())
    assert back.hash == net.hash
