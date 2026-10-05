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
