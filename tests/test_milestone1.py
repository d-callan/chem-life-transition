import numpy as np

from chemtrans.dynamics.ode import simulate
from chemtrans.networks.core import Reaction, ReactionNetwork
from chemtrans.networks.generators import random_crn
from chemtrans.observables.assembly import (
    assembly_distribution,
    path_assembly_index,
)
from chemtrans.observables.closure import (
    autocatalytic_activity,
    organizational_closure,
    raf_core,
    raf_reactions,
)
from chemtrans.observables.spectral import spectral_gap, spectral_radius


def test_random_crn_reproducible_and_connectivity():
    a = random_crn(40, 80, p_cat=0.3, seed=7)
    b = random_crn(40, 80, p_cat=0.3, seed=7)
    assert a.stoichiometry().tolist() == b.stoichiometry().tolist()

    lo = random_crn(60, 120, p_cat=0.0, seed=1)
    hi = random_crn(60, 120, p_cat=0.9, seed=1)
    n_cat = lambda net: sum(bool(r.catalysts) for r in net.reactions)
    assert n_cat(lo) == 0
    assert n_cat(hi) > n_cat(lo)
    assert spectral_radius(hi) > spectral_radius(lo)


def test_ode_food_only_chemostat():
    # no reactions: food sits at feed/dilution, everything else zero
    net = ReactionNetwork(species=["F", "X"], reactions=[], food={0})
    _, X = simulate(net, (0, 100), feed_rate=1.0, dilution_rate=0.1,
                    x0=np.zeros(2))
    assert np.isclose(X[-1, 0], 10.0, rtol=1e-3)
    assert X[-1, 1] == 0.0


def test_ode_nonnegative_and_runs():
    net = random_crn(30, 50, p_cat=0.4, seed=3)
    _, X = simulate(net, (0, 100), feed_rate=1.0, dilution_rate=0.2)
    assert X.shape[1] == 30
    assert X.min() > -1e-6


def _raf_net() -> ReactionNetwork:
    # food F. C+D -> E (required catalyst E' = species 4, produced by core)
    # Build a 3-reaction autocatalytic triangle:
    #   F -> A  cat A (required)      A makes A... simple RAF:
    #   F -> A  cat B                 B required-catalyst made inside
    #   F -> B  cat A
    # Both reactions need the other species => neither is self-bootstrapping
    # unless seed exists. Simpler RAF: F -> A cat A (self-loop required cat).
    return ReactionNetwork(
        species=["F", "A", "G"],
        reactions=[
            # self-catalyzed, catalyst-required: classic minimal RAF
            Reaction({0: 1}, {1: 1}, catalysts=[1], rate=0.0, cat_rate=5.0),
            # spontaneous side reaction, not part of RAF (no catalyst)
            Reaction({0: 1}, {2: 1}, rate=1.0),
        ],
        food={0},
    )


def test_raf_core_finds_self_maintaining_set():
    net = _raf_net()
    core_rxns = raf_reactions(net)
    assert core_rxns == {0}          # only the catalyzed self-loop
    assert raf_core(net) == {1}      # species A, not food, not G


def test_raf_empty_when_catalyst_not_internal():
    # catalyst is a food species => externally supplied, not self-maintaining
    net = ReactionNetwork(
        species=["F", "A"],
        reactions=[Reaction({0: 1}, {1: 1}, catalysts=[0], rate=0.0)],
        food={0},
    )
    assert raf_reactions(net) == set()


def test_organizational_closure_between_0_and_1():
    net = _raf_net()
    x = np.array([5.0, 2.0, 1.0])
    phi = organizational_closure(net, net.propensity(x))
    assert 0.0 < phi <= 1.0


def test_autocatalytic_activity():
    net = _raf_net()
    assert autocatalytic_activity(net, np.array([5.0, 0.0, 0.0])) == 0.0
    assert autocatalytic_activity(net, np.array([5.0, 2.0, 0.0])) > 0.0


def test_path_assembly_index():
    # F -> A -> B chain: A has index 1, B index 2, unreachable inf
    net = ReactionNetwork(
        species=["F", "A", "B", "Z"],
        reactions=[
            Reaction({0: 1}, {1: 1}),
            Reaction({1: 1}, {2: 1}),
        ],
        food={0},
    )
    A = path_assembly_index(net)
    assert A[0] == 0 and A[1] == 1 and A[2] == 2 and not np.isfinite(A[3])


def test_assembly_distribution_normalized():
    net = random_crn(30, 60, p_cat=0.3, seed=5)
    edges, hist = assembly_distribution(net, np.ones(30))
    assert np.isclose(np.sum(hist * np.diff(edges)), 1.0, atol=0.2)


def test_spectral_gap_nonneg():
    net = random_crn(20, 40, p_cat=0.5, seed=2)
    assert spectral_gap(net) >= 0.0
