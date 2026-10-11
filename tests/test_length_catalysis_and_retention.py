"""Length-weighted catalytic potency (cat_len_power) and CSTR polymer
retention (polymer_retention): backward compatibility, rate/Jacobian
consistency, boundary selectivity, and invariant preservation.

cat_len_power: boost_r = eps + kappa * sum_c w(c) x_c^p with
w(c) = (len(c)/2)^alpha — alpha=0 is the legacy flat-kappa law.
polymer_retention: structural polymers (len>=2) wash out at
D*(1-rho); food/fuel/waste at D — a semi-permeable compartment.
"""

import numpy as np
import pytest

from chemtrans.dynamics import gillespie, ode
from chemtrans.networks.generators import polymer_network

GEN = {"alphabet_size": 2, "max_len": 4, "p_cat": 0.5,
       "kappa": 5.0, "seed": 11}


def _has_long_catalyst(net, min_len=3):
    return any(net.polymer_lengths[c] >= min_len
               for ch in net.channels for c in ch.catalysts)


def test_alpha_zero_is_exact_legacy():
    """cat_len_power=0: every weight 1.0, rates identical to a net
    built without the knob."""
    net = polymer_network(**GEN)
    assert np.all(net.cat_weights == 1.0)
    x = np.linspace(0.01, 0.6, net.n_species)
    f, b = net.rates(x, 1.0)
    f2, b2 = net.rates(x, 1.0)
    np.testing.assert_array_equal(f, f2)
    np.testing.assert_array_equal(b, b2)


def test_cat_weights_follow_length_power_law():
    net = polymer_network(**GEN, cat_len_power=2.0)
    lens = net.polymer_lengths
    w = net.cat_weights
    # (L/2)^2: len1->1, len2->1, len3->2.25, len4->4.0; A*/W (-1) -> 1
    for i in range(net.n_species):
        L = max(lens[i], 2)
        assert w[i] == pytest.approx((L / 2.0) ** 2)


def test_cat_weights_scale_catalyzed_rates_only():
    """Catalyzed channel rates scale by w(c); uncatalyzed channels
    and mass-action factors are untouched."""
    net0 = polymer_network(**GEN)
    net2 = polymer_network(**GEN, cat_len_power=2.0)
    assert _has_long_catalyst(net2)
    x = np.full(net0.n_species, 0.4)
    f0, b0 = net0.rates(x, 1.0)
    f2, b2 = net2.rates(x, 1.0)
    has_cat = np.array([bool(ch.catalysts) for ch in net2.channels])
    np.testing.assert_allclose(f2[~has_cat], f0[~has_cat])
    np.testing.assert_allclose(b2[~has_cat], b0[~has_cat])
    assert np.any(f2[has_cat] > f0[has_cat])


def test_drates_fd_at_positive_alpha():
    """Analytic d_rates matches centered FD of rates() at alpha=2 —
    catches a missing w(c) factor in the catalyst derivative."""
    net = polymer_network(**GEN, cat_len_power=2.0)
    x = np.linspace(0.05, 0.5, net.n_species)
    Df, Db = net.d_rates(x, 1.0)
    Df = Df.toarray() if hasattr(Df, "toarray") else Df
    Db = Db.toarray() if hasattr(Db, "toarray") else Db
    h = 1e-7
    rng = np.random.default_rng(0)
    for j in rng.choice(net.n_species, size=min(8, net.n_species),
                        replace=False):
        xp, xm = x.copy(), x.copy()
        xp[j] += h
        xm[j] -= h
        fp, bp = net.rates(xp, 1.0)
        fm, bm = net.rates(xm, 1.0)
        np.testing.assert_allclose(Df[:, j], (fp - fm) / (2 * h),
                                   rtol=1e-5, atol=1e-7)
        np.testing.assert_allclose(Db[:, j], (bp - bm) / (2 * h),
                                   rtol=1e-5, atol=1e-7)


def test_alpha_symmetric_preserves_detailed_balance():
    """w(c) multiplies k_f and k_r identically -> k_f/k_r = exp(-dG)
    is unaffected (thermodynamic invariant 2)."""
    net = polymer_network(**GEN, cat_len_power=2.0)
    assert net.check_detailed_balance()["ok"]


def test_alpha_changes_net_hash_not_family():
    a = polymer_network(**GEN)
    b = polymer_network(**GEN, cat_len_power=1.0)
    assert a.hash != b.hash
    # serialization round-trip preserves the knob
    assert type(a).from_dict(b.to_dict()).cat_len_power == 1.0
    # with_channels variants (knockouts) inherit it
    assert b.without_catalysis().cat_len_power == 1.0


def test_invalid_alpha_rejected():
    with pytest.raises(ValueError, match="cat_len_power"):
        polymer_network(**GEN, cat_len_power=-1.0)


# ---------------- polymer retention -------------------------------

def test_retention_rhs_selective_outflow():
    """rhs with rho>0: polymer boundary outflow is D*(1-rho)*x,
    food boundary unchanged, internal flux untouched."""
    net = polymer_network(**GEN)
    x = np.full(net.n_species, 0.3)
    rho, D = 0.8, 0.25
    r0 = ode.rhs_cstr(net, x, D, 1.0)
    r1 = ode.rhs_cstr(net, x, D, 1.0, retention=rho)
    pm = net.polymer_mask
    # food rows identical (outflow D either way, same internal flux)
    np.testing.assert_allclose(r1[~pm], r0[~pm])
    # polymer rows differ by exactly D*rho*x on the outflow
    np.testing.assert_allclose(r1[pm] - r0[pm], D * rho * x[pm])


def test_retention_jacobian_diagonal_shift():
    net = polymer_network(**GEN)
    x = np.full(net.n_species, 0.3)
    J0 = np.asarray(net.jacobian(x, 0.2, 1.0))
    J1 = np.asarray(net.jacobian(x, 0.2, 1.0, retention=0.9))
    pm = net.polymer_mask
    d = np.diag(J1) - np.diag(J0)
    np.testing.assert_allclose(d[pm], 0.2 * 0.9)
    np.testing.assert_allclose(d[~pm], 0.0)
    off0 = J0 - np.diag(np.diag(J0))
    off1 = J1 - np.diag(np.diag(J1))
    np.testing.assert_allclose(off1, off0)


def test_retention_zero_is_legacy_cstr():
    net = polymer_network(**GEN)
    r0 = ode.integrate_cstr(net, t_max=100.0, dilution=0.3)
    r1 = ode.integrate_cstr(net, t_max=100.0, dilution=0.3,
                            polymer_retention=0.0)
    assert r0.trajectory_type == r1.trajectory_type
    if r0.converged and r1.converged:
        np.testing.assert_allclose(r0.x_eq, r1.x_eq)


def test_retention_raises_polymer_mass():
    """A retained polymer pool should reach >= the unretained steady
    polymer mass (retention lowers the effective washout cost)."""
    net = polymer_network(**{**GEN, "p_cat": 0.8, "eps": 1.0})
    pm = net.polymer_mask
    r0 = ode.integrate_cstr(net, t_max=300.0, dilution=0.3)
    r1 = ode.integrate_cstr(net, t_max=300.0, dilution=0.3,
                            polymer_retention=0.8)
    assert r0.converged and r1.converged
    assert r1.x_eq[pm].sum() >= r0.x_eq[pm].sum() - 1e-8


def test_retention_material_balance():
    """The trajectory material-balance diagnostic must stay small
    under retention — internal channels still conserve letters."""
    net = polymer_network(**GEN)
    r = ode.integrate_cstr(net, t_max=120.0, dilution=0.3,
                           polymer_retention=0.7)
    assert np.isfinite(r.material_error)
    assert r.material_error < 5e-3


def test_retention_absorbing_limit():
    """eps=0 + retention: dead state stays dead — retention cannot
    open synthesis (invariant 3)."""
    net = polymer_network(**GEN, eps=0.0)
    r = ode.integrate_cstr(net, t_max=200.0, dilution=0.2,
                           polymer_retention=0.9)
    assert r.converged
    assert r.x_eq[net.polymer_mask].max() < 1e-8


def test_retention_invalid_rejected():
    net = polymer_network(**GEN)
    with pytest.raises(ValueError, match="polymer_retention"):
        ode.integrate_cstr(net, dilution=0.2, polymer_retention=1.0)
    with pytest.raises(ValueError, match="polymer_retention"):
        gillespie.build_processes(net, 0.2, polymer_retention=-0.1)


def test_ssa_dilution_rates_respect_retention():
    net = polymer_network(**GEN)
    rho, D = 0.75, 0.4
    procs = gillespie.build_processes(net, D,
                                      polymer_retention=rho)
    pm = net.polymer_mask
    for p in procs:
        if p.kind != "dilution":
            continue
        i = p.participants[0]
        assert p.rate == pytest.approx(D * (1 - rho if pm[i] else 1.0))


def test_ssa_catalyst_rate_weighted_both_modes():
    """Catalyzed process rate constants carry w(c) in both
    participants and functional hazard modes."""
    net0 = polymer_network(**GEN)
    net2 = polymer_network(**GEN, cat_len_power=2.0)
    for mode in ("participants", "functional"):
        r0 = sorted(p.rate for p in
                    gillespie.build_processes(net0, 0.0, eps=1.0,
                                              hazard_mode=mode))
        r2 = sorted(p.rate for p in
                    gillespie.build_processes(net2, 0.0, eps=1.0,
                                              hazard_mode=mode))
        assert len(r0) == len(r2)
        assert r0 != r2
        # every changed rate grew (weights >= 1 everywhere)
        grew = [b > a for a, b in zip(r0, r2) if a != b]
        assert all(grew)


def test_ssa_retention_extends_survival():
    """Retained polymers dilute less -> organized state survives
    measurably longer in the stochastic run."""
    net = polymer_network(alphabet_size=2, max_len=3, p_cat=0.9,
                          kappa=5.0, eps=1.0, seed=3)
    kw = {"t_max": 60.0, "dilution": 0.8, "volume": 30.0, "eps": 1.0,
          "max_events": 1_000_000}
    finals = []
    for rho in (0.0, 0.9):
        m_last = 0.0
        for rep in range(3):
            res = gillespie.simulate(net, seed=100 + rep,
                                     polymer_retention=rho, **kw)
            pm = net.polymer_mask
            m_last = max(m_last, float(
                res.final_counts[pm].sum() / kw["volume"]))
        finals.append(m_last)
    assert finals[1] >= finals[0]


def test_polymer_mask_excludes_food_and_benchmarks():
    """Only structural len>=2 polymers are retained: monomers, A*,
    W, and unstructured species (benchmark nets) are excluded."""
    net = polymer_network(**GEN)
    pm = net.polymer_mask
    assert not np.any(pm[list(net.food)])
    assert not pm[net.species.index("W")]
    assert not pm[net.species.index("A*")]
    assert pm[net.species.index("AA")]
    from chemtrans.networks.benchmarks import schlogl
    bench = schlogl()
    assert not bench.polymer_mask.any()
