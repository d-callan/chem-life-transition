"""Observables: RAF fixed point, mass_raf, assembly index, hysteretic
region state, survival estimators, compositional signatures."""

import numpy as np
import pytest

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
    # entry is credited at the first crossing, not end of dwell window
    assert out[0][1] and out[0][0] == 10.0


def test_first_passage_honors_censor_time():
    # grid extends to t_max but the run was truncated at t=5 — censor
    # time must be 5, not the padded grid endpoint 19
    t = np.arange(20.0)
    o = np.zeros(20)
    out = persistence.first_passage_times([(t, o, 5.0)],
                                          theta_hi=0.6, dwell=4.0)
    assert out[0] == (5.0, False)
    # state above threshold only after the censor time -> still censored
    o[10:] = 0.9
    out = persistence.first_passage_times([(t, o, 5.0)],
                                          theta_hi=0.6, dwell=4.0)
    assert out[0] == (5.0, False)


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


def _max_ma_label(net):
    """Label of a highest-exact-MA species in the pool."""
    best, ba = None, -1
    for i in range(net.n_species):
        if net.is_food[i]:
            continue
        a, est = assembly.assembly_index(net.species[i])
        if est == "exact" and a > ba:
            best, ba = net.species[i], a
    return best, ba


def test_assembly_signature_absolute_positive_control():
    """All nonfood mass in the maximal-MA species -> signature fires."""
    net = polymer_network(alphabet_size=2, max_len=4, seed=0)
    lab, amax = _max_ma_label(net)
    x = np.zeros(net.n_species)
    x[0] = 0.5                              # food
    i = net.species.index(lab)
    x[i] = 0.2
    sig = assembly.assembly_signature_absolute(
        net, x, a_cut=None, mass_floor=0.01)
    assert sig["a_cut"] == amax
    assert sig["n_signature"] == 1
    assert sig["mass_signature_abs"] == 0.2 * len(lab)
    # full spectrum emitted
    nonfood = sum(1 for i in range(net.n_species)
                  if not net.is_food[i])
    assert len(sig["species"]) == nonfood


def test_assembly_signature_absolute_mass_floor():
    """Trace mass below the absolute floor must not sign even when it
    is a large FRACTION of a collapsing pool."""
    net = polymer_network(alphabet_size=2, max_len=4, seed=0)
    lab, amax = _max_ma_label(net)
    x = np.zeros(net.n_species)
    x[0] = 0.5
    i = net.species.index(lab)
    x[i] = 1e-6                             # 100% of nonfood mass, tiny abs
    sig = assembly.assembly_signature_absolute(
        net, x, a_cut=amax, mass_floor=0.01)
    assert sig["n_signature"] == 0
    assert sig["mass_signature_abs"] == 0.0
    # but the fraction diagnostic still sees it as "high" — proving the
    # floor is what vetoes, not the fraction
    rec = [r for r in sig["species"] if r["label"] == lab][0]
    assert rec["high_assembly"] and not rec["high_abundance"]
    assert rec["mass_fraction"] == 1.0


def test_assembly_signature_explicit_cut():
    net = polymer_network(alphabet_size=2, max_len=4, seed=0)
    x = np.zeros(net.n_species)
    x[0] = 0.5
    for i in range(net.n_species):
        if not net.is_food[i]:
            x[i] = 0.1                      # all species abundant
    sig = assembly.assembly_signature_absolute(
        net, x, a_cut=1, mass_floor=0.01)
    assert sig["a_cut"] == 1
    assert sig["n_signature"] == sig["n_signature"] > 0
    # bound estimators never qualify
    for r in sig["species"]:
        if r["estimator"] == "bound":
            assert not r["high_assembly"]


def test_derive_a_cut_synthetic():
    """Known null: M_L = 10*0.9^L, V=100, binary ->
    E[L]=1000*0.45^L/L (per-molecule mass correction) crosses 1 at
    L*=7, a_cut=5, extrapolated."""
    lm = {L: 10 * 0.9 ** L for L in range(1, 7)}
    rec = assembly.derive_a_cut(lm, volume=100.0, alphabet_size=2)
    assert rec["c"] == pytest.approx(10.0, rel=1e-6)
    assert rec["r"] == pytest.approx(0.9, rel=1e-6)
    assert rec["L_star"] == 7
    assert rec["a_cut"] == 5
    assert rec["extrapolated"] and not rec["vacuous"]


def test_derive_a_cut_vacuous():
    """A null too dilute to resolve any specific string must flag
    vacuous rather than return a permissive cut."""
    lm = {2: 1e-4, 3: 1e-7, 4: 1e-10}
    rec = assembly.derive_a_cut(lm, volume=10.0, alphabet_size=2)
    assert rec["a_cut"] <= 0 and rec["vacuous"]


def test_derive_a_cut_insufficient_data():
    rec = assembly.derive_a_cut({2: 0.0, 3: 1e-9}, volume=10.0)
    assert rec["a_cut"] is None and "reason" in rec
