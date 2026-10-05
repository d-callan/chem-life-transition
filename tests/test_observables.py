"""Observables: RAF fixed point, mass_raf, assembly index, hysteretic
region state, survival estimators, compositional signatures."""

import numpy as np

from chemtrans.networks.core import Channel, ReactionNetwork
from chemtrans.networks.generators import polymer_network
from chemtrans.observables import (
    assembly, closure, heredity, persistence)


def raf_net():
    """F -> P1 catalyzed by P2;  F + P1 -> P2 catalyzed by P1 (a 2-channel
    food-generated RAF), plus a food-disconnected catalyzed channel that
    must NOT be in the RAF."""
    # 0 F (food), 1 P1, 2 P2, 3 X, 4 X2, 5 X3 — food-disconnected branch
    chans = [
        Channel(reactants=((0, 1),), products=((1, 1),),
                catalysts=(2,), k_f=1.0, k_r=1.0, synth=True),
        Channel(reactants=((0, 1), (1, 1)), products=((2, 1),),
                catalysts=(1,), k_f=1.0, k_r=1.0, synth=True),
        Channel(reactants=((3, 1),), products=((4, 1),),
                catalysts=(5,), k_f=1.0, k_r=1.0, synth=True),
    ]
    species = ["F", "P1", "P2", "X", "X2", "X3"]
    net = ReactionNetwork(species, chans, (0,),
                          food_inlet=np.array([1, 0, 0, 0, 0, 0]),
                          comp=np.ones((1, 6)), g=np.zeros(6),
                          kappa=5.0, eps=0.0)
    return net


def test_maxraf_fixed_point():
    net = raf_net()
    core = closure.raf_reactions(net)
    assert core == {0, 1}
    sp = closure.raf_species(net)
    assert sp == {1, 2}


def test_mass_raf_counts_only_core():
    net = raf_net()
    x = np.array([0.5, 0.2, 0.3, 0.9, 0.0, 0.0])
    # RAF species P1, P2 mass = 0.5 ; non-food total = 0.5+0.9 = 1.4
    m = closure._mass_raf(net, x)
    assert np.isclose(m, 0.5 / 1.4)


def test_assembly_index_values():
    assert assembly.assembly_index("A") == (0, "exact")
    a_ab, est = assembly.assembly_index("AB")
    assert (a_ab, est) == (1, "exact")
    a_abab, _ = assembly.assembly_index("ABAB")
    assert a_abab == 2          # AB built once, reused
    a_rand, _ = assembly.assembly_index("ABBA")
    assert a_rand == 3          # no reuse possible


def test_region_state_hysteresis():
    t = np.linspace(0, 10, 11)
    o = np.array([0.0, 0.5, 0.7, 0.5, 0.45, 0.3, 0.2, 0.65, 0.5, 0.55,
                  0.2])
    st = persistence.region_state(t, o, theta_hi=0.6, theta_lo=0.4)
    # on at idx2 (0.7), stays on through band (0.45), off at 0.3 (idx5)
    assert st[2] and st[4] and not st[5]
    assert not st[0] and not st[-1]


def test_first_passage_dwell():
    t = np.arange(20.0)
    o = np.zeros(20)
    o[5:8] = 0.9            # too short for dwell=5
    o[10:16] = 0.9          # qualifies
    out = persistence.first_passage_times([(t, o)], theta_hi=0.6,
                                          dwell=4.0)
    assert out[0][1] and out[0][0] == 14.0


def test_kaplan_meier_known():
    # 4 subjects: events at 2,4; censor at 3
    times = np.array([2.0, 3.0, 4.0, 6.0])
    events = np.array([True, False, True, False])
    tk, S, se = persistence.kaplan_meier(times, events)
    assert np.isclose(S[1], 1 - 1 / 4)
    # at t=4 only 2 are still at risk (the t=3 censor left the risk set)
    assert np.isclose(S[-1], 0.75 * (1 - 1 / 2))


def test_signature_is_compositional():
    net = polymer_network(alphabet_size=2, max_len=3, seed=0)
    x = np.zeros(net.n_species)
    x[0] = 0.5
    for i in range(net.n_species):
        if net.species[i] == "AAB":
            x[i] = 0.2
        if net.species[i] == "BBB":
            x[i] = 0.3
    sig = heredity.attractor_signature(net, x)
    assert np.isclose(sig.sum(), 1.0)
    same = heredity.compositional_distance(sig, sig)
    assert same == 0.0
