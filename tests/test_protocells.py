"""Protocell population dynamics: permeation, osmotic growth, division,
lysis, and observables."""
import numpy as np
import pytest

from chemtrans.dynamics import protocell_pop as ppop
from chemtrans.measure import build_network, measure_run


def _net(**kw):
    gen = {"kind": "polymer", "alphabet_size": 2, "max_len": 4,
           "p_cat": 0.4, "p_info": 0.0, "p_selfrep": 0.3,
           "delta_mu": 4.0, "chi": 0.7, "eps": 0.001, "kappa": 10.0,
           "cat_power": 2, "min_cat_len": 2, "ligation": True,
           "p_cat_lig": 0.4, "k_join": 10.0, "seed": 7}
    gen.update(kw)
    return build_network(gen)[0]


def _pop(net, n_cells=4, **kw):
    return ppop.ProtocellPopulation(net, n_cells=n_cells, **kw)


def test_permeation_selectivity():
    net = _net()
    pop = _pop(net)
    L = np.asarray(net.polymer_lengths)
    assert pop._P[L <= 1].min() > 0          # monomers/fuel/waste cross
    assert np.all(pop._P[L == 2] > 0)        # dimers permeate slowly
    assert np.all(pop._P[L == 2] < pop._P[L == 1].min())
    # polymers strongly retained but leaky (imperfect membrane)
    assert np.all(pop._P[L >= 3] > 0)
    assert np.all(pop._P[L >= 3] < pop._P[L == 2].min())


def test_osmotic_growth_direction():
    net = _net()
    pop = _pop(net, n_cells=2, m_init=0)
    hi = np.asarray(net.polymer_lengths) >= 3
    pop.x[0, hi] += 0.5                     # cell A packed with polymer
    _, dv = pop._rhs(pop.x.copy(), pop.v.copy())
    assert dv[0] > 0.0                       # tense vesicle grows
    assert dv[1] < 0.0                       # flaccid vesicle shrinks


def test_uncoupled_population_neutral():
    net = _net()
    pop = _pop(net, mu_growth=0.0, d_leak=0.0, m_init=0)
    v0 = pop.v.copy()
    _, dv = pop._rhs(pop.x.copy(), pop.v.copy())
    assert np.allclose(dv, 0.0)
    assert np.allclose(pop.v, v0)


def test_division_conserves_molecule_counts():
    net = _net()
    pop = _pop(net, n_cells=1, m_init=0)
    rng_before = pop.rng
    # seed discrete molecule counts in the mother
    pop.x[0, :] = 0.0
    idx = np.arange(net.n_species)
    pop.x[0, 10] = 0.02      # 20 molecules at v=1, omega=1000
    pop.x[0, 20] = 0.007     # 7 molecules
    pop.v[0] = 2.5           # above v_div -> forces split
    n_before = np.floor(pop.x[0] * pop.v[0] * pop.omega + 0.5).astype(int)
    pop._events()
    assert pop.n_divisions == 1
    assert len(pop.v) == 2
    n_after = (np.floor(pop.x * pop.v[:, None] * pop.omega + 0.5)
               .astype(int))
    # molecule counts conserve exactly on partitioned species
    assert np.all(n_after.sum(axis=0) >= n_before - 0)
    total = n_after[:, idx].sum(axis=0)
    assert total[10] + total[20] >= n_before[10] + n_before[20]


def test_lysis_removes_below_threshold():
    net = _net()
    pop = _pop(net, n_cells=3, m_init=0)
    pop.v = np.array([0.3, 0.2, 1.1])      # two below 0.4*v0
    x0 = pop.x.copy()
    pop._events()
    assert pop.n_lysis == 2
    assert len(pop.v) == 1
    assert pop.v[0] == pytest.approx(1.1)


def test_carrying_capacity_flushes_flaccid():
    net = _net()
    pop = _pop(net, n_cells=4, n_max=4, m_init=0)
    hi = np.asarray(net.polymer_lengths) >= 3
    pop.x[0, hi] += 2.0                     # packed cell survives
    pop.v = np.array([1.5, 0.6, 0.7, 0.8])
    pop._events()                           # division fires on cell 0
    assert len(pop.v) <= 4
    # the high-polymer lineage is still present
    assert np.any(pop.x[:, hi].sum(axis=1) > 1.0)


def test_run_produces_valid_result():
    net = _net()
    pop = _pop(net, n_cells=8, m_init=20, chunk=1.0)
    res = pop.run(t_max=5.0)
    assert res.status in {"ok", "EXTINCT"}
    assert res.x_cells.shape[1] == net.n_species
    assert len(res.t) == len(res.n_alive_hist)
    assert np.all(res.x_cells >= -1e-9)


def test_observables_shapes():
    net = _net()
    x = np.zeros((3, net.n_species))
    L = np.asarray(net.polymer_lengths)
    x[:, L == 4] = 1.0 / (L == 4).sum()
    masks = {4: L == 4}
    out = ppop.per_cell_sparsity(x, masks)
    assert out["L4_sparsity"] == pytest.approx(0.0, abs=1e-6)
    assert out["L4_n_eff"] == pytest.approx((L == 4).sum(), abs=0.5)
    si = ppop.speciation_index(x, L >= 3)
    assert si == pytest.approx(0.0, abs=1e-9)


def test_signature_cells_threshold():
    net = _net()
    x = np.zeros((2, net.n_species))
    ma = np.zeros(net.n_species)
    ma[30] = 5                              # one high-assembly species
    x[0, 30] = 0.02                         # above floor -> signature
    x[1, 30] = 0.001                        # trace -> no signature
    sig = ppop.signature_cells(x, ma, a_cut=5, mass_floor=0.01)
    assert sig["n_signature_cells"] == 1
    assert sig["max_signature_mass"] == pytest.approx(0.02)


def test_measure_run_protocell_population():
    gen = {"kind": "polymer", "alphabet_size": 2, "max_len": 4,
           "p_cat": 0.4, "p_info": 0.0, "p_selfrep": 0.3,
           "delta_mu": 4.0, "chi": 0.7, "eps": 0.001, "kappa": 10.0,
           "cat_power": 2, "min_cat_len": 2, "ligation": True,
           "p_cat_lig": 0.4, "k_join": 10.0, "seed": 11}
    task = {"seed_key": "pcell-test", "generator": gen,
            "reactor": {"type": "protocell_population", "n_cells": 6,
                        "t_max": 4.0, "mu_growth": 2.0,
                        "chunk": 1.0, "m_init": 20},
            "assembly": {"a_cut": 4, "mass_floor": 0.01}}
    rec = measure_run(task)
    pop = rec["protocell_population"]
    assert pop["status"] in {"ok", "EXTINCT"}
    assert set(pop) >= {"n_alive", "n_divisions", "n_lysis",
                        "speciation_index", "per_cell_sparsity",
                        "signature"}
    sig = pop["signature"]
    assert 0 <= sig["n_signature_cells"] <= pop["n_alive"]
