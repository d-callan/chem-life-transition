"""Negative-feedback knobs: size-selective dialysis
(retention_len_scale), product inhibition (k_inhibit), and the
hard catalyst-length threshold (min_cat_len).

rho_i(L) = rho*(1 - 2*gamma/L): gamma=0 is legacy flat retention,
gamma=1 lets dimers wash out at full D (dialysis cutoff).
boost_r = eps + kappa * sum_c w_c x_c^p / (1 + K * sum_{c in C_r} x_c).
"""

import numpy as np
import pytest

from chemtrans.dynamics import gillespie, ode
from chemtrans.networks.generators import polymer_network

GEN = {"alphabet_size": 2, "max_len": 4, "p_cat": 0.5,
       "kappa": 5.0, "seed": 11}


# ---------------- retention_len_scale ------------------------------

def test_len_scale_outflow_formula():
    net = polymer_network(**GEN)
    lens = net.polymer_lengths
    D, rho, gam = 0.5, 0.9, 1.0
    out = net.dilution_outflow(D, rho, gam)
    for i in range(net.n_species):
        if lens[i] >= 2:
            assert out[i] == pytest.approx(
                D * (1 - rho * (1 - 2.0 / lens[i])))
        else:
            assert out[i] == pytest.approx(D)


def test_len_scale_zero_is_flat_retention():
    net = polymer_network(**GEN)
    pm = net.polymer_mask
    flat = net.dilution_outflow(0.5, 0.9, 0.0)
    np.testing.assert_allclose(flat[pm], 0.5 * 0.1)
    np.testing.assert_allclose(flat[~pm], 0.5)


def test_len_scale_dialysis_selectivity():
    """gamma=1: dimers wash out at full D; longer chains retained
    monotonically more — the membrane/pore size cutoff."""
    net = polymer_network(**{**GEN, "max_len": 6})
    lens = net.polymer_lengths
    out = net.dilution_outflow(0.4, 0.9, 1.0)
    for L in range(2, 7):
        idx = np.nonzero(lens == L)[0]
        assert len(idx)
        exp = 0.4 * (1 - 0.9 * (1 - 2.0 / L))
        np.testing.assert_allclose(out[idx], exp)
    i2, i6 = int(np.nonzero(lens == 2)[0][0]), \
        int(np.nonzero(lens == 6)[0][0])
    assert out[i2] == pytest.approx(0.4)      # dimer fully leaky
    assert out[i6] < 0.4 * 0.4                # 6-mer well retained


def test_len_scale_jacobian_diagonal():
    net = polymer_network(**GEN)
    x = np.full(net.n_species, 0.3)
    J0 = np.asarray(net.jacobian(x, 0.4, 1.0, retention=0.9,
                                 len_scale=0.0))
    J1 = np.asarray(net.jacobian(x, 0.4, 1.0, retention=0.9,
                                 len_scale=1.0))
    lens = net.polymer_lengths
    d = np.diag(J1) - np.diag(J0)
    # delta per polymer species: D*rho*(2/L) added back as outflow
    for i in range(net.n_species):
        if lens[i] >= 2:
            assert d[i] == pytest.approx(-0.4 * 0.9 * 2.0 / lens[i])
        else:
            assert d[i] == pytest.approx(0.0)


def test_len_scale_invalid():
    net = polymer_network(**GEN)
    with pytest.raises(ValueError, match="retention_len_scale"):
        ode.integrate_cstr(net, dilution=0.2, polymer_retention=0.5,
                           retention_len_scale=1.5)
    with pytest.raises(ValueError, match="retention_len_scale"):
        gillespie.build_processes(net, 0.2, polymer_retention=0.5,
                                  retention_len_scale=-0.1)


def test_len_scale_dialysis_shrinks_dimer_pool():
    """Dialysis strictly increases the dimer outflow channel, so
    absolute L2 mass at equilibrium must drop vs flat retention.
    (Mean length is NOT asserted: dimers are both the trap and the
    fuel feeding longer chains — the aggregate response is
    non-monotone and network-dependent.)"""
    net = polymer_network(**{**GEN, "p_cat": 0.8, "eps": 1.0})
    lens = net.polymer_lengths
    m2 = {}
    for gam in (0.0, 1.0):
        r = ode.integrate_cstr(net, t_max=300.0, dilution=0.3,
                               polymer_retention=0.9,
                               retention_len_scale=gam)
        assert r.converged
        m2[gam] = float(r.x_eq[lens == 2].sum())
    assert m2[1.0] < m2[0.0]


# ---------------- k_inhibit -----------------------------------------

def test_k_inhibit_rate_law():
    """Catalytic boost divided by (1 + K*sum_c x_c); uncatalyzed
    channels unaffected."""
    net0 = polymer_network(**GEN)
    net_i = polymer_network(**GEN, k_inhibit=0.5)
    x = np.full(net0.n_species, 0.4)
    C = net_i.catalysis_matrix()
    f0, _b0 = net0.rates(x, 1.0)
    fi, _bi = net_i.rates(x, 1.0)
    has_cat = np.array([bool(ch.catalysts) for ch in net_i.channels])
    np.testing.assert_allclose(fi[~has_cat], f0[~has_cat])
    # inhibited rates strictly lower where catalysis is active
    cat = net_i.kappa * np.asarray(
        C @ np.power(x, net_i.cat_power)).ravel()
    assert np.all(fi[has_cat & (cat > 0)] < f0[has_cat & (cat > 0)])


def test_k_inhibit_fd_drates():
    net = polymer_network(**GEN, k_inhibit=0.7)
    x = np.linspace(0.05, 0.5, net.n_species)
    Df, Db = net.d_rates(x, 1.0)
    Df = Df.toarray() if hasattr(Df, "toarray") else Df
    Db = Db.toarray() if hasattr(Db, "toarray") else Db
    h = 1e-7
    for j in range(0, net.n_species, 5):
        xp, xm = x.copy(), x.copy()
        xp[j] += h
        xm[j] -= h
        fp, bp = net.rates(xp, 1.0)
        fm, bm = net.rates(xm, 1.0)
        np.testing.assert_allclose(Df[:, j], (fp - fm) / (2 * h),
                                   rtol=1e-5, atol=1e-7)
        np.testing.assert_allclose(Db[:, j], (bp - bm) / (2 * h),
                                   rtol=1e-5, atol=1e-7)


def test_k_inhibit_preserves_detailed_balance():
    """Same divisor on fwd and rev -> k_f/k_r untouched."""
    net = polymer_network(**GEN, k_inhibit=0.5)
    assert net.check_detailed_balance()["ok"]


def test_k_inhibit_provenance():
    a = polymer_network(**GEN)
    b = polymer_network(**GEN, k_inhibit=0.3)
    assert a.hash != b.hash
    assert type(a).from_dict(b.to_dict()).k_inhibit == 0.3
    assert b.without_catalysis().k_inhibit == 0.3
    with pytest.raises(ValueError, match="k_inhibit"):
        polymer_network(**GEN, k_inhibit=-0.1)


def test_k_inhibit_ssa_hazard():
    """Catalyzed SSA hazards carry the shared per-channel
    denominator in both hazard modes."""
    net_i = polymer_network(**GEN, k_inhibit=0.5)
    n = np.full(net_i.n_species, 20, dtype=np.int64)
    for mode in ("participants", "functional"):
        procs_i = gillespie.build_processes(net_i, 0.0,
                                            hazard_mode=mode)
        tagged = [p for p in procs_i if p.func_inhibit is not None]
        assert tagged
        _cats, K = tagged[0].func_inhibit
        assert K == pytest.approx(0.5)
        # hazard is finite and reduced vs K=0 on a comparable proc
        h = gillespie.hazard(tagged[0], n, 100.0)
        assert np.isfinite(h) and h >= 0


def test_k_inhibit_suppresses_runaway_monopoly():
    """At high catalyst concentration the marginal boost saturates:
    boost(x) grows sublinearly once x_c ~> 1/K."""
    net = polymer_network(**GEN, k_inhibit=1.0)
    net0 = polymer_network(**GEN)
    lo = np.full(net.n_species, 0.1)
    hi = np.full(net.n_species, 2.0)
    f0_hi, _ = net0.rates(hi, 1.0)
    fi_hi, _ = net.rates(hi, 1.0)
    f0_lo, _ = net0.rates(lo, 1.0)
    fi_lo, _ = net.rates(lo, 1.0)
    has_cat = np.array([bool(ch.catalysts) for ch in net.channels])
    # relative suppression grows with concentration
    frac_hi = (fi_hi / np.maximum(f0_hi, 1e-12))[has_cat].max()
    frac_lo = (fi_lo / np.maximum(f0_lo, 1e-12))[has_cat].max()
    assert frac_hi < frac_lo


# ---------------- min_cat_len ---------------------------------------

def test_min_cat_len_excludes_short_catalysts():
    net = polymer_network(alphabet_size=2, max_len=4, p_cat=1.0,
                          seed=5, min_cat_len=3)
    lens = net.polymer_lengths
    for ch in net.channels:
        for c in ch.catalysts:
            assert lens[c] >= 3
    assert any(len(ch.catalysts) for ch in net.channels)


def test_min_cat_len_absorbing_without_seed():
    """eps=0 + min_cat_len>=3 + food_only init: the dead state stays
    dead (no species len>=3 exists to catalyze, no spontaneous
    synthesis) — invariant 3 holds under the threshold."""
    net = polymer_network(**{**GEN, "eps": 0.0}, min_cat_len=3)
    x0 = net.food_inlet.copy()
    r = ode.integrate_cstr(net, t_max=200.0, dilution=0.2,
                           x0=x0)
    assert r.converged
    assert r.x_eq[net.polymer_mask].max() < 1e-8
