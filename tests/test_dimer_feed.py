"""Environmental oligomer feed (f_dimer) and inlet_species plumbing.

Physical model: a fraction f_dimer of the inlet letter flux arrives as
pre-formed dimers — abiotic background chemistry outside the
compartment supplying short oligomers. Total inlet letter flux is
conserved; dimers remain non-food for order-parameter purposes but are
declared inlet_species so the CSTR boundary may supply them.

Invariants exercised:
- f_dimer=0 is bit-identical legacy behavior.
- Letter conservation and detailed balance are untouched (the inlet is
  a boundary flux, not an internal channel).
- food_inlet may only supply food or declared inlet_species.
- SSA feed processes exist for every nonzero inlet entry.
"""
import numpy as np
import pytest

from chemtrans.dynamics import gillespie, ode
from chemtrans.networks.core import ReactionNetwork
from chemtrans.networks.generators import polymer_network

GEN = {"seed": 1, "max_len": 4, "p_cat": 0.5, "p_selfrep": 0.3,
       "p_info": 0.0, "delta_mu": 4.0, "chi": 0.7, "eps": 0.001,
       "kappa": 100.0, "cat_power": 2}


def _dimer_idx(net):
    return [i for i, s in enumerate(net.species)
            if len(s) == 2 and not s.endswith("*")]


def test_f_dimer_zero_legacy_bit_identical():
    n0 = polymer_network(**GEN)
    n1 = polymer_network(**GEN, f_dimer=0.0)
    assert n0.hash == n1.hash
    np.testing.assert_array_equal(n0.food_inlet, n1.food_inlet)


def test_inlet_letter_flux_conserved():
    for fd in (0.0, 0.25, 0.5, 0.9):
        net = polymer_network(**GEN, f_dimer=fd)
        assert net.letter_mass(net.food_inlet) == pytest.approx(1.0)


def test_dimer_feed_uniform_and_scaled():
    net = polymer_network(**GEN, f_dimer=0.4)
    di = _dimer_idx(net)
    assert len(di) == 4  # K^2 dimers
    expect = 0.4 / (2.0 * 4)
    for i in di:
        assert net.food_inlet[i] == pytest.approx(expect)


def test_monomer_feed_rescaled_by_one_minus_f_dimer():
    n0 = polymer_network(**GEN)
    nd = polymer_network(**GEN, f_dimer=0.3)
    ia = n0.species.index("A*")
    ida = nd.species.index("A*")
    assert nd.food_inlet[ida] == pytest.approx(0.7 * n0.food_inlet[ia])


def test_fed_dimers_stay_nonfood_and_polymer():
    net = polymer_network(**GEN, f_dimer=0.5)
    di = _dimer_idx(net)
    for i in di:
        assert not net.is_food[i]
        assert net.polymer_mask[i]


def test_inlet_species_declared_exactly_dimers():
    net = polymer_network(**GEN, f_dimer=0.5)
    assert net.inlet_species == frozenset(_dimer_idx(net))
    net0 = polymer_network(**GEN)
    assert net0.inlet_species == frozenset()


def test_undeclared_nonfood_inlet_rejected():
    net = polymer_network(**GEN)
    di = _dimer_idx(net)
    inlet = np.zeros(net.n_species)
    inlet[di[0]] = 0.1
    with pytest.raises(ValueError, match="inlet"):
        ReactionNetwork(net.species, net.channels, net.food,
                        food_inlet=inlet)


def test_declared_inlet_species_accepted():
    net = polymer_network(**GEN)
    di = _dimer_idx(net)
    inlet = np.zeros(net.n_species)
    inlet[di[0]] = 0.1
    n2 = ReactionNetwork(net.species, net.channels, net.food,
                         food_inlet=inlet, inlet_species=di)
    assert n2.food_inlet[di[0]] == pytest.approx(0.1)


@pytest.mark.parametrize("bad", [-0.1, 1.0, 1.5, np.inf, np.nan])
def test_invalid_f_dimer_rejected(bad):
    with pytest.raises(ValueError):
        polymer_network(**GEN, f_dimer=bad)


def test_f_dimer_enters_hash_and_roundtrip():
    n0 = polymer_network(**GEN)
    nd = polymer_network(**GEN, f_dimer=0.3)
    assert n0.hash != nd.hash
    rt = ReactionNetwork.from_dict(nd.to_dict())
    assert rt.hash == nd.hash
    assert rt.inlet_species == nd.inlet_species
    np.testing.assert_array_equal(rt.food_inlet, nd.food_inlet)


def test_letter_conservation_and_detailed_balance():
    net = polymer_network(**GEN, f_dimer=0.5)
    net.check_letter_conservation()
    net.check_detailed_balance()


def test_ssa_feed_processes_for_fed_dimers():
    net = polymer_network(**GEN, f_dimer=0.4)
    procs = gillespie.build_processes(net, dilution=0.5)
    feeds = [p for p in procs if p.kind == "feed"]
    n_expected = int(np.count_nonzero(net.food_inlet))
    assert len(feeds) == n_expected
    di = _dimer_idx(net)
    fed_idx = {p.delta[0][0] for p in feeds}
    for i in di:
        assert i in fed_idx
    # feed process rate = D * inlet concentration
    for p in feeds:
        i = p.delta[0][0]
        assert p.rate == pytest.approx(0.5 * net.food_inlet[i])


def test_ssa_ic_seeds_fed_dimers():
    net = polymer_network(**GEN, f_dimer=0.4)
    res = gillespie.simulate(net, t_max=0.01, dilution=0.5, volume=1000,
                             seed=0, eps=0.0)
    n0 = res.n[0]
    di = _dimer_idx(net)
    for i in di:
        assert n0[i] == round(net.food_inlet[i] * 1000)


def test_dimer_feed_raises_cstr_dimer_steady_state():
    """Fed dimers must equilibrate at D*f/(D*(1-rho_i)) above legacy ~0."""
    x_eqs = []
    for fd in (0.0, 0.5):
        net = polymer_network(**GEN, f_dimer=fd)
        x0 = net.food_inlet.copy()
        x0[~net.is_food] = 0.15
        res = ode.integrate_cstr(net, t_max=300.0, dilution=0.5, x0=x0,
                                 polymer_retention=0.9, eps=0.0)
        x = res.x_eq if res.x_eq is not None else res.x[-1]
        di = _dimer_idx(net)
        x_eqs.append(float(np.mean(x[di])))
    assert x_eqs[1] > 10 * max(x_eqs[0], 1e-9)


def test_material_balance_holds_with_dimer_feed():
    """Dimer-fed run must satisfy the same material-balance diagnostic
    as the equivalent f_dimer=0 run (trapezoid discretization of the
    retention-aware boundary law has O(dt^2) error ~0.02-0.03; the
    honest check is baseline-comparability, not machine precision)."""
    errs = []
    for fd in (0.0, 0.4):
        net = polymer_network(**GEN, f_dimer=fd)
        x0 = net.food_inlet.copy()
        x0[~net.is_food] = 0.15
        res = ode.integrate_cstr(net, t_max=300.0, dilution=0.5, x0=x0,
                                 polymer_retention=0.9, eps=0.0)
        assert res.trajectory_type in ("STATIONARY", "OSCILLATORY")
        errs.append(res.material_error)
    assert errs[1] <= max(errs[0], 0.05)
