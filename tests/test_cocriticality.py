"""Co-criticality apparatus: life_vector schema, multistability and
competition task blocks, ghost-signature regression, and the paired
coincidence analysis on synthetic data."""

import numpy as np
import pytest

from chemtrans.analysis import cocriticality as cc
from chemtrans.measure import measure_run
from chemtrans.networks.generators import polymer_network
from chemtrans.observables import heredity
from chemtrans.observables.life_vector import life_vector

P2_GEN = {"kind": "polymer", "max_len": 4, "alphabet_size": 2,
          "p_cat": 0.5, "p_info": 0.0, "p_selfrep": 0.3,
          "delta_mu": 4.0, "chi": 0.7, "eps": 0.001,
          "kappa": 500.0, "cat_power": 2}


def test_ghost_signature_regression():
    """Dead endpoints (~1e-35 mass) must not normalize into
    attractors when min_mass is enforced."""
    net = polymer_network(max_len=4, alphabet_size=2, p_cat=0.5,
                          seed=7, cat_power=2, kappa=500.0)
    dead = [net.food_inlet + 1e-35 * np.random.default_rng(i).random(
        net.n_species) for i in range(4)]
    cl = heredity.cluster_attractors(net, dead, tol=0.05,
                                     min_mass=0.01)
    assert cl["n_attractors"] == 0 and cl["n_kept"] == 0


def test_life_vector_undefined_semantics():
    # structural-only record: dynamical observables undefined, not 0
    net = polymer_network(max_len=4, alphabet_size=2, p_cat=0.3,
                          seed=3)
    rec = {"structural": {"raf_fraction_channels": 0.2,
                          "n_raf_species": 5}}
    lv = life_vector(rec)
    assert lv["O_raf"]["defined"] and lv["O_raf"]["value"] == 0.2
    for k in ("O_maint", "O_persist", "O_hered", "O_select",
              "O_assembly"):
        assert not lv[k]["defined"] and lv[k]["value"] is None


def test_life_vector_full_record():
    task = {"generator": dict(P2_GEN), "eps": 0.0,
            "seed_key": "e1land|net3",
            "reactor": {"type": "cstr", "dilution": 0.5,
                        "t_max": 400.0,
                        "multistability": {"n_ic": 16,
                                           "min_mass": 0.01,
                                           "t_max": 600.0}}}
    rec = measure_run(task)
    lv = life_vector(rec)
    assert lv["O_raf"]["defined"]
    assert lv["O_maint"]["defined"]          # dead eq is still defined
    assert lv["O_hered"]["kind"] == "discrete"
    assert isinstance(lv["O_hered"]["value"], int)
    assert not lv["O_persist"]["defined"]    # no stochastic block
    assert not lv["O_select"]["defined"]     # competition not requested


def test_multistability_deterministic():
    """Same seed_key -> identical attractor count and centers."""
    task = {"generator": dict(P2_GEN), "eps": 0.0,
            "seed_key": "e1land|net3",
            "reactor": {"type": "cstr", "dilution": 0.5,
                        "multistability": {"n_ic": 16,
                                           "min_mass": 0.01,
                                           "t_max": 600.0}}}
    a = measure_run(task)["multistability"]
    b = measure_run(task)["multistability"]
    assert a["n_attractors"] == b["n_attractors"] >= 2
    assert a["attractor_top_species"] == b["attractor_top_species"]
    assert a["min_mass"] == 0.01


def test_competition_undefined_when_monostable():
    """p=1 net has a single attractor -> competition must report
    incomplete with an explicit reason, never a silent zero."""
    gen = dict(P2_GEN, cat_power=1, kappa=20.0)
    task = {"generator": gen, "eps": 0.0, "seed_key": "e1land|net3",
            "reactor": {"type": "cstr", "dilution": 0.5,
                        "multistability": {"n_ic": 8,
                                           "min_mass": 0.01,
                                           "t_max": 400.0},
                        "competition": {"replicates": 2, "volume": 4,
                                        "t_max": 50.0}}}
    comp = measure_run(task)["competition"]
    assert comp["completed"] is False
    assert comp["reason"] == "fewer_than_two_organized_attractors"


def test_competition_mixture_init_equal_mass():
    """eq_mass protocol seeds equal monomer mass per lineage."""
    from chemtrans.measure import (_multistability_block, _nonfood_mass,
                                   derive_seed)
    from chemtrans.dynamics import gillespie
    net = polymer_network(**{k: v for k, v in P2_GEN.items()
                             if k != "kind"},
                          seed=derive_seed("e1land|net3|gen"))
    ms = _multistability_block(net, 0.5, 0.0, "e1land|net3",
                               n_ic=16, min_mass=0.01, t_max=600.0)
    assert ms["n_attractors"] >= 2
    xA = np.array(ms["attractor_x0"]["0"])
    xB = np.array(ms["attractor_x0"]["1"])
    sigA = heredity.attractor_signature(net, xA)
    sigB = heredity.attractor_signature(net, xB)
    m_t = min(_nonfood_mass(net, xA), _nonfood_mass(net, xB))
    n0 = gillespie.counts_from_concentration(
        0.5 * m_t * (sigA + sigB) + net.food_inlet, 64.0)
    xf = n0 / 64.0
    # each lineage contributes exactly m_t/2 of nonfood mass
    sig = heredity.attractor_signature(net, xf)
    dA = heredity.compositional_distance(sig, sigA)
    dB = heredity.compositional_distance(sig, sigB)
    assert dA == pytest.approx(dB, abs=0.05)


def _synthetic_net_series(lc, n_nets=4, n_pts=11, span=(0.5, 1.6)):
    """Fake records: hinge data with critical point lc."""
    recs = []
    rng = np.random.default_rng(0)
    for j in range(n_nets):
        for lam in np.linspace(*span, n_pts):
            recs.append({
                "net_hash": f"net{j}",
                "control": {"lam": float(lam)},
                "life_vector": {
                    "O_maint": {"defined": True,
                                "value": float(
                                    max(0, lam - lc)
                                    + 0.01 * rng.normal())},
                    "O_hered": {"defined": True,
                                "value": int(lam >= lc)},
                }})
    return recs


def test_coincidence_matrix_synthetic():
    """Coincident observables -> coincident_within_tol; separated ones
    -> separated. Verdicts come from per-family paired diffs."""
    recs = _synthetic_net_series(lc=1.0)
    for r in recs:  # add a second observable shifted by +0.4
        lam = r["control"]["lam"]
        r["life_vector"]["O_persist"] = {
            "defined": True,
            "value": float(max(0, lam - 1.4))}
    by_net = cc.group_by_network(recs)
    rep = cc.coincidence_report(
        by_net, observables=("O_maint", "O_persist", "O_hered"),
        n_boot=60, seed=1)
    assert rep["pairwise"]["O_maint|O_persist"]["verdict"] == \
        "separated"
    assert rep["pairwise"]["O_maint|O_hered"]["verdict"] == \
        "coincident_within_tol"
    # per-family paired diffs are reported exactly
    pd = rep["pairwise"]["O_maint|O_hered"]["paired_diffs"]
    assert pd and all(abs(d) <= rep["tol_pair"] for d in pd.values())
    pe = rep["point_estimates"]
    assert pe["O_maint"]["lc"] == pytest.approx(1.0, abs=0.15)
    assert pe["O_hered"]["lc"] == pytest.approx(1.0, abs=0.15)
    # per-family estimates exist for every family
    assert set(rep["family_estimates"]["O_maint"]) == set(by_net)


def test_conditional_observable_never_fit():
    """O_select is conditional on multistability: even when defined it
    is never fit as a global lambda_c curve."""
    recs = _synthetic_net_series(lc=1.0)
    for r in recs:
        lam = r["control"]["lam"]
        r["life_vector"]["O_select"] = {
            "defined": True, "value": float(abs(lam - 1.0) > 0.2)}
    rep = cc.coincidence_report(
        cc.group_by_network(recs),
        observables=("O_maint", "O_select"), n_boot=20, seed=0)
    assert rep["point_estimates"]["O_select"]["lc"] is None
    assert rep["family_estimates"]["O_select"] == {}
    assert rep["pairwise"]["O_maint|O_select"]["verdict"] == \
        "incomparable"
