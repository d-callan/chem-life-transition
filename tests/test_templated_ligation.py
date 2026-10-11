"""Templated ligation: Watson-Crick edge weights and sequence sparsity.

Physical model: a strand T catalyzes ligation P1+P2 -> P12 with weight
w = min(exp(m*delta_g_pair - u*delta_g_mismatch), w_max_template),
where (m, u) count matches/mismatches of T against the complement of
P12 (A<->B). Replaces the random catalyst lottery on ligation channels
with syntax-derived potency.

Invariants exercised:
- Legacy networks are bit-identical when templated_ligation is off.
- Weights are monotonic in pairing affinity and capped (no overflow).
- Catalysis stays symmetric fwd/rev; detailed balance untouched.
- k_inhibit denominator uses binary incidence, not edge weights.
- cat_w round-trips through serialization and enters the hash.
"""
import numpy as np
import pytest

from chemtrans.dynamics import gillespie
from chemtrans.networks.core import ReactionNetwork
from chemtrans.networks.generators import _complement, polymer_network
from chemtrans.observables.sparsity import sequence_entropy

GEN = {"seed": 1, "max_len": 4, "p_cat": 0.5, "p_selfrep": 0.3,
       "p_info": 0.0, "delta_mu": 4.0, "chi": 0.7, "eps": 0.001,
       "kappa": 100.0, "cat_power": 2}
LIG = {"ligation": True, "p_cat_lig": 0.5, "k_join": 5.0}
TEM = {"templated_ligation": True, "delta_g_pair": 2.0,
       "delta_g_mismatch": 2.0, "w_max_template": 100.0,
       "w_min_template": 1.5}


def _lig_channels(net):
    return [ch for ch in net.channels if ch.mechanism == "ligation"]


def test_complement_mapping():
    assert _complement("ABBA") == "BAAB"
    assert _complement("AAAA") == "BBBB"
    assert _complement(_complement("ABAB")) == "ABAB"


def test_weight_matches_mismatch_counts():
    """Perfect complement gets cap; one mismatch -> exp(3*g - 1*u)."""
    net = polymer_network(**GEN, **LIG, **TEM)
    for ch in _lig_channels(net):
        prod = net.structure[ch.products[0][0]]
        comp = _complement(prod)
        for k, c in enumerate(ch.catalysts):
            t = net.species[c]
            m = sum(1 for a, b in zip(t, comp) if a == b)
            u = len(t) - m
            expect = min(np.exp(m * 2.0 - u * 2.0), 100.0)
            assert ch.cat_w[k] == pytest.approx(expect)


def test_weights_monotonic_in_delta_g_pair():
    # dg=0 with a mismatch penalty: no strand reaches w_min — the
    # binding-recognition regime is simply absent (empty edge set)
    net0 = polymer_network(**GEN, **LIG, templated_ligation=True,
                           delta_g_pair=0.0, delta_g_mismatch=2.0)
    assert all(not ch.cat_w for ch in _lig_channels(net0))
    weights = []
    for dg in (1.0, 2.0):
        net = polymer_network(**GEN, **LIG, templated_ligation=True,
                              delta_g_pair=dg, delta_g_mismatch=2.0)
        w = [w_ for ch in _lig_channels(net) for w_ in ch.cat_w]
        weights.append(float(np.max(w)))
    assert weights[0] < weights[1]


def test_boost_capped_no_overflow():
    net = polymer_network(**GEN, **LIG, templated_ligation=True,
                          delta_g_pair=200.0, delta_g_mismatch=0.0,
                          w_max_template=50.0, w_min_template=1.0)
    for ch in _lig_channels(net):
        assert all(w <= 50.0 for w in ch.cat_w)
        assert all(np.isfinite(w) and w > 0 for w in ch.cat_w)


def test_threshold_prunes_weak_templates():
    """w_min_template excludes low-affinity strands entirely."""
    net = polymer_network(**GEN, **LIG, templated_ligation=True,
                          delta_g_pair=1.0, delta_g_mismatch=2.0,
                          w_min_template=10.0, w_max_template=100.0)
    for ch in _lig_channels(net):
        for k, c in enumerate(ch.catalysts):
            assert ch.cat_w[k] >= 10.0


def test_legacy_ligation_bit_identical():
    n0 = polymer_network(**GEN, **LIG)
    n1 = polymer_network(**GEN, **LIG, templated_ligation=False)
    assert n0.hash == n1.hash


def test_cat_w_default_empty_legacy_catalysts():
    net = polymer_network(**GEN, **LIG)
    assert all(not ch.cat_w for ch in net.channels)
    C = net.catalysis_matrix()
    assert set(np.unique(C)) <= {0.0, 1.0}


def test_template_params_enter_hash_and_roundtrip():
    n0 = polymer_network(**GEN, **LIG, **TEM)
    n1 = polymer_network(**GEN, **LIG, **{**TEM, "delta_g_pair": 3.0})
    n2 = polymer_network(**GEN, **LIG)
    assert n0.hash != n1.hash
    assert n0.hash != n2.hash
    rt = ReactionNetwork.from_dict(n0.to_dict())
    assert rt.hash == n0.hash
    assert all(a.cat_w == b.cat_w
               for a, b in zip(rt.channels, n0.channels))


def test_dense_sparse_rates_and_derivatives():
    rng = np.random.default_rng(0)
    net = polymer_network(**GEN, **LIG, **TEM)
    nets = polymer_network(**GEN, **LIG, **TEM, sparse=True)
    x = np.abs(rng.normal(0.5, 0.3, net.n_species)) + 0.01
    f, b = net.rates(x)
    fs, bs = nets.rates(x)
    np.testing.assert_allclose(f, fs, atol=1e-12)
    np.testing.assert_allclose(b, bs, atol=1e-12)
    Df, Db = nets.d_rates(x)
    Df = Df.toarray() if hasattr(Df, "toarray") else Df
    # central FD on a catalyst column
    j = [c for ch in _lig_channels(net) for c in ch.catalysts][0]
    h = 1e-7
    xp, xm = x.copy(), x.copy()
    xp[j] += h; xm[j] -= h
    fp, _ = nets.rates(xp); fm, _ = nets.rates(xm)
    fd = (fp - fm) / (2 * h)
    np.testing.assert_allclose(Df[:, j], fd,
                               rtol=1e-4, atol=1e-4)


def test_catalysis_symmetric_and_detailed_balance():
    net = polymer_network(**GEN, **LIG, **TEM)
    rng = np.random.default_rng(2)
    x = np.abs(rng.normal(0.5, 0.3, net.n_species)) + 0.01
    f, b = net.rates(x)
    net0 = polymer_network(**{**GEN, "kappa": 0.0}, **LIG, **TEM)
    f0, b0 = net0.rates(x)
    # catalytic boost identical on both directions: (f/b)/(f0/b0) == 1
    mask = (b > 0) & (f0 > 0) & (b0 > 0)
    ratio = (f[mask] / b[mask]) / (f0[mask] / b0[mask])
    np.testing.assert_allclose(ratio, 1.0, rtol=1e-12)
    assert net.check_detailed_balance()["ok"]


def test_k_inhibit_uses_binary_incidence():
    """Denominator counts catalyst concentrations once each, never
    weighted by edge affinity."""
    net = polymer_network(**GEN, **LIG, **TEM, k_inhibit=0.5)
    rng = np.random.default_rng(3)
    x = np.abs(rng.normal(0.5, 0.3, net.n_species)) + 0.01
    Cb = net._catalysis_binary()
    Cw = net.catalysis_matrix()
    assert set(np.unique(Cb)) <= {0.0, 1.0}
    assert np.max(Cw) > 1.0  # templated edges present
    # replicate the denominator computation manually
    denom = 1.0 + 0.5 * np.asarray(Cb @ x).ravel()
    r = next(i for i, ch in enumerate(net.channels)
             if ch.mechanism == "ligation" and len(ch.catalysts) > 1)
    expect = 1.0 + 0.5 * sum(x[c] for c in net.channels[r].catalysts)
    assert denom[r] == pytest.approx(expect)


def test_ssa_kcat_folds_edge_weight():
    net = polymer_network(**GEN, **LIG, **TEM)
    procs = gillespie.build_processes(net, dilution=0.5,
                                      hazard_mode="participants")
    for r, ch in enumerate(net.channels):
        if not ch.catalysts:
            continue
        for kc, c in enumerate(ch.catalysts):
            edge = float(ch.cat_w[kc]) if ch.cat_w else 1.0
            expect = ch.k_f * net.kappa * net.cat_weights[c] * edge
            hits = [p for p in procs
                    if p.channel == r and p.rate > ch.k_f
                    and set(p.delta) == set(ch.nu.items())]
            # at least one catalyzed fwd process with folded weight
            assert any(p.rate == pytest.approx(expect)
                       for p in hits), (r, ch)


def test_validation_rejects_bad_params():
    with pytest.raises(ValueError, match="ligation"):
        polymer_network(**GEN, templated_ligation=True)
    with pytest.raises(ValueError, match="alphabet"):
        polymer_network(**GEN, **LIG, templated_ligation=True,
                        alphabet_size=1)
    with pytest.raises(ValueError):
        polymer_network(**GEN, **LIG, templated_ligation=True,
                        delta_g_pair=-1.0)
    with pytest.raises(ValueError):
        polymer_network(**GEN, **LIG, templated_ligation=True,
                        w_min_template=10.0, w_max_template=1.0)


def test_cat_w_length_mismatch_rejected():
    net = polymer_network(**GEN)
    ch = net.channels[0]
    bad = type(ch)(reactants=ch.reactants, products=ch.products,
                   catalysts=(), k_f=ch.k_f, k_r=ch.k_r,
                   synth=ch.synth, mechanism=ch.mechanism,
                   cat_w=(1.0,))  # cat_w without catalysts
    with pytest.raises(ValueError):
        ReactionNetwork(net.species, [bad] + list(net.channels[1:]),
                        net.food, structure=net.structure,
                        comp=net.comp, g=net.g)


# ---------------- sequence sparsity observables --------------------

def _net_l4():
    return polymer_network(**GEN, **LIG, **TEM)


def test_entropy_uniform_is_max():
    net = _net_l4()
    idx = np.nonzero(net.polymer_lengths == 4)[0]
    x = np.zeros(net.n_species)
    x[idx] = 1.0
    out = sequence_entropy(net, x, 4)
    assert out["entropy"] == pytest.approx(4.0)
    assert out["n_eff"] == pytest.approx(16.0)
    assert out["sparsity"] == pytest.approx(0.0)


def test_entropy_single_sequence():
    net = _net_l4()
    idx = np.nonzero(net.polymer_lengths == 4)[0]
    x = np.zeros(net.n_species)
    x[idx[3]] = 2.5
    out = sequence_entropy(net, x, 4)
    assert out["entropy"] == pytest.approx(0.0)
    assert out["n_eff"] == pytest.approx(1.0)
    assert out["sparsity"] == pytest.approx(1.0)
    assert out["top_sequences"][0]["mass_fraction"] == pytest.approx(1.0)


def test_entropy_empty_pool_is_null_not_zero():
    net = _net_l4()
    out = sequence_entropy(net, np.zeros(net.n_species), 4)
    assert out["sparsity"] is None
    assert out["entropy"] == 0.0
    assert out["top_sequences"] == []


def test_entropy_two_peak_mixture():
    net = _net_l4()
    idx = np.nonzero(net.polymer_lengths == 4)[0]
    x = np.zeros(net.n_species)
    x[idx[0]] = 1.0
    x[idx[1]] = 1.0
    out = sequence_entropy(net, x, 4)
    assert out["entropy"] == pytest.approx(1.0)
    assert out["n_eff"] == pytest.approx(2.0)


def test_negative_solver_excursions_clipped_and_diagnosed():
    """ODE endpoints can carry small negative concentrations. They must
    not produce entropy < 0 or sparsity > 1 — clip and report via
    neg_mass_abs."""
    net = _net_l4()
    idx = np.nonzero(net.polymer_lengths == 4)[0]
    x = np.zeros(net.n_species)
    x[idx[0]] = 1.0
    x[idx[1]] = -0.05   # solver excursion, not a physical concentration
    out = sequence_entropy(net, x, 4)
    assert out["sparsity"] == pytest.approx(1.0)
    assert 0.0 <= out["entropy"] <= out["entropy_max"]
    assert out["neg_mass_abs"] == pytest.approx(0.05)


def test_numerically_absent_pool_is_null():
    net = _net_l4()
    idx = np.nonzero(net.polymer_lengths == 4)[0]
    x = np.zeros(net.n_species)
    x[idx[0]] = 1e-12   # below MASS_FLOOR — shell effectively absent
    x[idx[1]] = -1e-11
    out = sequence_entropy(net, x, 4)
    assert out["sparsity"] is None


def test_pathological_negative_pool_is_null():
    """Negative excursion comparable to pool mass: the clipped vector
    is not a distribution — report undefined, not sparsity=1."""
    net = _net_l4()
    idx = np.nonzero(net.polymer_lengths == 4)[0]
    x = np.zeros(net.n_species)
    x[idx[0]] = 0.01
    x[idx[1]] = -0.02   # neg_mass (0.02) > 0.5 * pool_mass (0.01)
    out = sequence_entropy(net, x, 4)
    assert out["sparsity"] is None
    assert out["neg_mass_abs"] == pytest.approx(0.02)
