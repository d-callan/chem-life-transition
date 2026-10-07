"""Dynamics: falling-factorial hazards, CSTR mass balance, structured
convergence, exact SSA bookkeeping."""

import numpy as np

from chemtrans.dynamics import gillespie, ode
from chemtrans.networks.benchmarks import (
    immigration_death, transcritical_motif)
from chemtrans.networks.core import Channel, ReactionNetwork
from chemtrans.networks.generators import polymer_network


def _pick(procs, kind=None, channel=None):
    out = []
    for i, p in enumerate(procs):
        if kind is not None and p.kind != kind:
            continue
        if channel is not None and p.channel != channel:
            continue
        out.append((i, p))
    return out


def homodimer_net():
    """A + A <=> P catalyzed by P (forces merged multiplicity in SSA)."""
    ch = Channel(reactants=((0, 2),), products=((1, 1),),
                 catalysts=(1,), k_f=1.0, k_r=0.5, synth=True)
    return ReactionNetwork(
        ("A", "P"), [ch], (0,), comp=np.array([[1.0, 2.0]]),
        g=np.zeros(2), food_inlet=np.array([1.0, 0.0]),
        kappa=5.0, eps=0.0)


def test_falling_factorial_identical_reactants():
    net = homodimer_net()
    procs = gillespie.build_processes(net, dilution=0.0)
    # catalyzed fwd: participants A, A, P -> hazard k*K * nA(nA-1)nP / V^2
    cat_fwd = [p for _, p in _pick(procs, kind="channel")
               if p.participants == (0, 0, 1)]
    assert len(cat_fwd) == 1
    p = cat_fwd[0]
    n = np.array([3, 0])
    assert gillespie.hazard(p, n, 10.0) == 0.0   # no catalyst present
    n = np.array([3, 1])
    expect = 1.0 * 5.0 * (3 * 2) * 1 / 10.0 ** 2
    assert np.isclose(gillespie.hazard(p, n, 10.0), expect)


def test_impossible_event_cannot_fire():
    """Every process needing species 1 (P) must have zero hazard when
    n_P = 0 — impossible events are never clamped into existence."""
    net = homodimer_net()
    procs = gillespie.build_processes(net, dilution=0.0, eps=0.3)
    needs_p = [p for _, p in _pick(procs)
               if 1 in p.participants]
    assert needs_p
    for p in needs_p:
        assert gillespie.hazard(p, np.array([5, 0]), 10.0) == 0.0


def test_ssa_never_negative_counts():
    net = immigration_death(k_birth=3.0, k_death=1.0, f_inlet=0.5)
    res = gillespie.simulate(net, t_max=20.0, dilution=0.5, volume=50,
                             seed=0)
    assert (res.n >= 0).all()
    assert res.completed


def test_ssa_event_ordering_nonanticipatory():
    net = immigration_death(k_birth=2.0, k_death=1.0, f_inlet=0.0)
    res = gillespie.simulate(net, t_max=10.0, dilution=0.0, volume=10,
                             seed=1, record_dt=0.5)
    # non-anticipatory: X[0] is the pre-first-event state
    assert res.n[0, 1] == 0  # X count starts at 0 (no feed, D=0)


def test_cstr_converges_stationary_and_mass_balances():
    net, _ = transcritical_motif(f_inlet=1.0, k_lig=1.0, kappa=2.0,
                                 keq=np.inf, eps=0.01)
    r = ode.integrate_cstr(net, t_max=300.0, dilution=0.5)
    assert r.trajectory_type == "STATIONARY", r.message
    assert r.x_eq is not None
    assert r.material_error < 1e-3


def test_transcritical_lambda_c_numeric():
    """Numerical onset must bracket the analytic lambda_c = D/F^2."""
    D, F = 0.5, 1.0
    kappas = [0.2, 0.35, 0.5, 0.7, 1.0]
    cvals = []
    for kap in kappas:
        net, an = transcritical_motif(f_inlet=F, k_lig=1.0, kappa=kap,
                                      keq=np.inf, eps=1e-4)
        # near-critical points relax slowly; give them a long horizon
        r = ode.integrate_cstr(net, t_max=3000.0, dilution=D,
                               n_chunks=30)
        assert r.converged, (kap, r.message, r.trajectory_type)
        cvals.append(r.x_eq[2])
        assert np.isclose(an["kappa_c_of_D"](D), D / (1.0 * F ** 2))
    cvals = np.array(cvals)
    lam_c = D / F ** 2  # = 0.5; kappa_c same since k_lig=1
    below = cvals[np.array(kappas) < lam_c]
    above = cvals[np.array(kappas) > lam_c]
    assert below.max() < 0.05
    assert above.min() > 0.1
    # near-critical growth: c* ~ F - sqrt(D/(k_f kappa))
    assert np.isclose(cvals[-1], F - np.sqrt(D / (1.0 * kappas[-1])),
                      atol=0.02)


def test_schlogl_bistability():
    from chemtrans.networks.benchmarks import schlogl
    net = schlogl()
    # low seed vs high seed initial X
    x_lo = net.food_inlet.copy()
    x_hi = net.food_inlet.copy()
    x_hi[2] = 5.0
    r_lo = ode.integrate_cstr(net, t_max=200.0, dilution=1.0, x0=x_lo)
    r_hi = ode.integrate_cstr(net, t_max=200.0, dilution=1.0, x0=x_hi)
    assert r_lo.converged and r_hi.converged
    assert abs(r_lo.x_eq[2] - r_hi.x_eq[2]) > 0.5  # two distinct attractors


def test_closed_system_detailed_balance():
    """Closed (D=0) delta_mu=0 dynamics must approach zero net flux."""
    net = polymer_network(alphabet_size=2, max_len=2, seed=0,
                          delta_mu=0.0, eps=1.0, kappa=5.0)
    x0 = np.full(net.n_species, 0.1)
    r = ode.integrate_cstr(net, t_max=200.0, dilution=0.0, x0=x0,
                           eps=1.0)
    assert r.converged
    flux = np.abs(net.net_flux(r.x_eq, eps=1.0))
    assert flux.max() < 1e-5


def test_sparse_path_bitwise_equivalent():
    """Sparse mode must reproduce dense-path results: S, C, rates,
    d_rates, Jacobian, and the certified equilibrium."""
    from chemtrans.networks.core import ReactionNetwork
    kw = dict(max_len=4, alphabet_size=2, p_cat=0.6, p_info=0.7,
              p_selfrep=0.3, ligation=True, n_cat_max=3, p_multi=0.5,
              chi=0.7, delta_mu=4.0, kappa=500.0, cat_power=2, seed=5)
    nd = polymer_network(**kw)
    ns = ReactionNetwork(nd.species, nd.channels, nd.food,
                         structure=nd.structure, comp=nd.comp,
                         g=nd.g, food_inlet=nd.food_inlet,
                         kappa=nd.kappa, eps=nd.eps,
                         cat_power=nd.cat_power, sparse=True)
    assert ns.use_sparse and not nd.use_sparse
    rng = np.random.default_rng(0)
    x = np.abs(rng.normal(0, 0.3, nd.n_species))
    eps = 1e-3
    np.testing.assert_allclose(ns.stoichiometry().toarray(),
                               nd.stoichiometry())
    np.testing.assert_allclose(ns.catalysis_matrix().toarray(),
                               nd.catalysis_matrix())
    f1, b1 = nd.rates(x, eps)
    f2, b2 = ns.rates(x, eps)
    np.testing.assert_allclose(f2, f1)
    np.testing.assert_allclose(b2, b1)
    for xx in (x, np.zeros(nd.n_species)):
        Df1, Db1 = nd.d_rates(xx, eps)
        Df2, Db2 = ns.d_rates(xx, eps)
        np.testing.assert_allclose(Df2.toarray(), Df1)
        np.testing.assert_allclose(Db2.toarray(), Db1)
    np.testing.assert_allclose(ns.jacobian(x, 0.3, eps).toarray(),
                               nd.jacobian(x, 0.3, eps))


def test_sparse_cstr_equilibrium_matches_dense():
    """End-to-end: BDF+sparse-LU trajectory and lsmr polish must land
    on the same certified equilibrium as the dense path."""
    base = {"max_len": 5, "alphabet_size": 2,
            "p_cat": 0.6, "p_info": 0.7, "p_selfrep": 0.3,
            "ligation": True, "n_cat_max": 3, "p_multi": 0.5,
            "chi": 0.7, "delta_mu": 4.0, "kappa": 500.0,
            "cat_power": 2, "seed": 5}
    nd = polymer_network(**base)
    ns = polymer_network(**base, sparse=True)
    x0 = ns.food_inlet.copy()
    x0[~ns.is_food] = 0.15
    rd = ode.integrate_cstr(nd, t_max=200.0, dilution=0.3, x0=x0,
                            eps=0.0)
    rs = ode.integrate_cstr(ns, t_max=200.0, dilution=0.3, x0=x0,
                            eps=0.0)
    assert rs.trajectory_type == rd.trajectory_type == "STATIONARY"
    np.testing.assert_allclose(rs.x_eq, rd.x_eq, atol=1e-8, rtol=1e-5)
