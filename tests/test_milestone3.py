import numpy as np
import pandas as pd

from chemtrans.analysis.critical_points import (
    align_critical_points,
    estimate_critical_point,
    susceptibility,
)
from chemtrans.networks.core import Reaction, ReactionNetwork
from chemtrans.networks.generators import random_crn
from chemtrans.observables.persistence import (
    evolvability_proxy,
    reproductive_number,
    self_maintenance,
)


def test_exclusive_catalysis_blocks_spontaneous_routes():
    net = random_crn(40, 150, p_cat=0.8, p_cat_required=1.0,
                     exclusive_catalysis=True, seed=3)
    excl = set()
    for rxn in net.reactions:
        if rxn.rate == 0.0 and rxn.catalysts:
            excl.update(rxn.products)
    excl -= net.food
    for rxn in net.reactions:
        if rxn.rate > 0.0:
            assert not (set(rxn.products) & excl), \
                "spontaneous reaction producing exclusive species"
        if rxn.rate == 0.0 and rxn.catalysts:
            assert set(rxn.catalysts) <= excl, \
                "required catalyst outside exclusive set leaks a nucleation route"


def _catalyzed_net() -> ReactionNetwork:
    # F -> A cat A (self-sustaining once seeded), plus spontaneous F -> G
    return ReactionNetwork(
        species=["F", "A", "G"],
        reactions=[
            Reaction({0: 1}, {1: 1}, catalysts=[1], rate=0.0, cat_rate=5.0),
            Reaction({0: 1}, {2: 1}, rate=1.0),
        ],
        food={0},
    )


def test_reproductive_number():
    net = _catalyzed_net()
    x = np.array([5.0, 2.0, 1.0])
    r0 = reproductive_number(net, x, dilution_rate=0.1)
    assert r0 > 0
    # no catalysis at all -> R0 = 0
    net0 = ReactionNetwork(species=["F", "A"],
                           reactions=[Reaction({0: 1}, {1: 1})], food={0})
    assert reproductive_number(net0, np.array([5.0, 1.0])) == 0.0


def test_self_maintenance_bounds():
    net = _catalyzed_net()
    # organized steady state: A present at moderate concentration
    x_ss = np.array([5.0, 3.0, 1.0])
    s = self_maintenance(net, x_ss, feed_rate=1.0, dilution_rate=0.1,
                         t_recover=60.0, shocks=(0.5, 0.9))
    assert 0.0 <= s <= 1.0
    # no RAF at all -> S = 0
    net0 = ReactionNetwork(species=["F", "A"],
                           reactions=[Reaction({0: 1}, {1: 1})], food={0})
    assert self_maintenance(net0, np.array([5.0, 1.0]),
                            shocks=(0.5,)) == 0.0


def test_evolvability_proxy_bounds():
    net = random_crn(40, 200, p_cat=0.4, max_reactants=1,
                     n_food=8, seed=11)
    e = evolvability_proxy(net, n_mut=10, seed=0)
    assert 0.0 <= e <= 1.0


def _synthetic_df(stagger: bool) -> pd.DataFrame:
    """Ensemble with two observables showing a percolation-like transition:
    each rep independently 'flips on' with a logistic probability rising
    through lambda_c — so ensemble variance peaks AT lambda_c, like our
    real RAF data. staggered=True puts the two observables at different
    critical lambdas."""
    rng = np.random.default_rng(0)
    lambdas = np.linspace(0, 1, 9)
    lcs = (0.3, 0.7) if stagger else (0.5, 0.5)
    rows = []
    for lam in lambdas:
        for rep in range(12):
            rows.append({
                "p_cat": lam,
                "o1": 0.9 if rng.random() < 1/(1+np.exp(-(lam-lcs[0])/0.03))
                      else 0.05,
                "o2": 0.9 if rng.random() < 1/(1+np.exp(-(lam-lcs[1])/0.03))
                      else 0.05,
            })
    return pd.DataFrame(rows)


def test_critical_point_estimation():
    df = _synthetic_df(stagger=True)
    lams = np.sort(df.p_cat.unique())
    mat = np.array([df.loc[df.p_cat == l, "o1"].to_numpy() for l in lams])
    lc = estimate_critical_point(lams, mat)
    assert abs(lc - 0.3) < 0.15
    chi = susceptibility(mat)
    assert chi.shape == (9,)
    assert np.nanargmax(chi) == np.argmin(np.abs(lams - lc))


def test_alignment_detects_stagger_vs_colocation():
    df_stag = _synthetic_df(stagger=True)
    _, pair_s = align_critical_points(df_stag, ["o1", "o2"],
                                      n_bootstrap=50, seed=1)
    # lambda_c's separated by 0.4 >> grid step 0.125 -> not aligned
    assert pair_s.iloc[0]["aligned"] < 0.3

    df_same = _synthetic_df(stagger=False)
    _, pair_c = align_critical_points(df_same, ["o1", "o2"],
                                      n_bootstrap=50, seed=1)
    assert pair_c.iloc[0]["aligned"] > 0.7
