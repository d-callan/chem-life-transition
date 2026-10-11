"""Non-pairing monomer C and RNA-like folding competence.

pair_score: only A.B / B.A are Watson-Crick matches — self-pairs and
all C-involving pairs are mismatches.
folded_fraction: two-state Boltzmann occupancy of an intramolecular
hairpin (s[k] paired with s[L-1-k], min central loop 2). Floppy
homopolymers fold poorly; alternating/palindromic strands fold.

Effective template weight w_eff = f_fold(T) * w_pair(T, P12):
a catalyst must BOTH fold AND match the product.
"""
import numpy as np
import pytest

from chemtrans.networks.generators import (
    folded_fraction, hairpin_fold_stability, pair_score,
    polymer_network)
from chemtrans.observables.sparsity import (
    _is_complement, _stem_pairs, sequence_entropy)

GEN2 = {"seed": 1, "max_len": 6, "alphabet_size": 2, "p_cat": 0.5,
        "p_selfrep": 0.3, "p_info": 0.0, "delta_mu": 4.0, "chi": 0.7,
        "eps": 0.001, "kappa": 100.0, "cat_power": 2}
GEN3 = {**GEN2, "max_len": 4, "alphabet_size": 3}
LIG = {"ligation": True, "p_cat_lig": 0.5, "k_join": 5.0}
TEM = {"templated_ligation": True, "delta_g_pair": 2.0,
       "delta_g_mismatch": 2.0, "w_max_template": 100.0,
       "w_min_template": 1.5}


def _lig_channels(net):
    return [ch for ch in net.channels if ch.mechanism == "ligation"]


# ------------------- pairing rules --------------------------------

def test_pair_score_truth_table():
    assert pair_score("A", "B") == 1
    assert pair_score("B", "A") == 1
    assert pair_score("A", "A") == 0
    assert pair_score("B", "B") == 0
    assert pair_score("A", "C") == 0
    assert pair_score("C", "A") == 0
    assert pair_score("B", "C") == 0
    assert pair_score("C", "C") == 0


def test_is_complement_requires_pure_AB():
    assert _is_complement("AABB", "BBAA")
    assert not _is_complement("AABB", "AABB")
    assert not _is_complement("AACB", "BBAB")   # C never pairs
    assert not _is_complement("AABB", "BBA")    # length mismatch


# ------------------- hairpin folding ------------------------------

def test_stem_pairs_count():
    # L4, loop 2 -> l_stem=1: s[0] pairs s[3]
    assert _stem_pairs("ABAB") == 1   # A-B
    assert _stem_pairs("AAAA") == 0   # A-A
    assert _stem_pairs("ACCB") == 1   # A-B (C in loop, not stem)
    # L6, loop 2 -> l_stem=2: s[0].s[5], s[1].s[4]
    assert _stem_pairs("ABBAAB") == 2   # A.B and B.A both pair
    assert _stem_pairs("ABBBAB") == 2   # A.B match, B.A match
    assert _stem_pairs("ABBBBB") == 1   # A.B match, B.B mismatch
    assert _stem_pairs("CBAABC") == 0   # C.C mismatch, B.B mismatch
    assert _stem_pairs("BBBBBB") == 0


def test_fold_free_energy_sign():
    dg_folded = hairpin_fold_stability("ABBAAB", 2.0, 2.0, 1.0)
    dg_coil = hairpin_fold_stability("AAAAAA", 2.0, 2.0, 1.0)
    assert dg_folded < 0 < dg_coil   # folder favored, homopolymer not


def test_folded_fraction_monotonic_and_bounds():
    f0 = folded_fraction("ABBAAB", 0.0)
    f2 = folded_fraction("ABBAAB", 2.0)
    f5 = folded_fraction("ABBAAB", 5.0)
    coil5 = folded_fraction("AAAAAA", 5.0)
    assert f0 == 1.0                      # folding disabled at dg=0
    assert 0.0 < f2 < f5 <= 1.0           # stronger pairing -> more folded
    assert coil5 < 0.01                   # floppy strand stays unfolded


def test_fold_gate_prunes_coils_and_scales_weights():
    """Higher dg_fold strengthens imperfect stems too — the gate
    removes true coils (m_stem=0) while boosting folders, so edge
    count saturates but homopolymer catalysts disappear."""
    net0 = polymer_network(**GEN2, **LIG, **TEM, delta_g_fold=0.0)
    net5 = polymer_network(**GEN2, **LIG, **TEM, delta_g_fold=5.0)
    cats0 = {net0.species[c] for ch in _lig_channels(net0)
             for c in ch.catalysts}
    cats5 = {net5.species[c] for ch in _lig_channels(net5)
             for c in ch.catalysts}
    # homopolymers can never fold: excluded at dg_fold=5
    assert "AAAAAA" in cats0 and "AAAAAA" not in cats5
    assert "BBBBBB" in cats0 and "BBBBBB" not in cats5
    # every remaining L6 catalyst has at least one stem pair —
    # shorter coils retain w_eff = 0.047*cap >= w_min at these
    # params (the gate attenuates, it does not hard-exclude)
    assert all(_stem_pairs(s) >= 1 for s in cats5 if len(s) == 6)
    # folding weakens a weak folder's edge weight vs the ungated case
    def edge_w(net, prod_s, cat_s):
        for ch in _lig_channels(net):
            if net.structure[ch.products[0][0]] == prod_s:
                for k, c in enumerate(ch.catalysts):
                    if net.species[c] == cat_s:
                        return ch.cat_w[k]
        return None
    # "BBBBAB" (m_stem=1) templates "AAAABB" with one end-mismatch;
    # folding gate attenuates its weight relative to the ungated case
    w0 = edge_w(net0, "AAAABB", "BBBBAB")
    w5 = edge_w(net5, "AAAABB", "BBBBAB")
    assert w0 is not None and w5 is not None and w5 < w0


def test_fold_params_enter_hash_and_default_legacy():
    n0 = polymer_network(**GEN2, **LIG, **TEM)
    n1 = polymer_network(**GEN2, **LIG, **TEM, delta_g_fold=2.0)
    n2 = polymer_network(**GEN2, **LIG, **TEM, delta_g_loop=5.0)
    assert n0.hash != n1.hash
    # dg_fold=0 disables folding: loop penalty alone cannot change net
    n3 = polymer_network(**GEN2, **LIG, **TEM,
                         delta_g_fold=0.0, delta_g_loop=99.0)
    assert n0.hash != n2.hash      # param recorded
    assert n0.hash != n3.hash
    # but the CHANNELS are physically identical (folding off)
    assert [ch.cat_w for ch in n0.channels] == \
           [ch.cat_w for ch in n3.channels]


# ------------------- K=3 non-pairing monomer ----------------------

def test_k3_templated_network_builds_and_validates():
    net = polymer_network(**GEN3, **LIG, **TEM)
    net.check_letter_conservation()
    assert net.check_detailed_balance()["ok"]
    assert any("C" in s for s in net.species)


def test_c_contaminated_products_get_weak_or_no_templates():
    """A product containing C can never full-match a template; pure
    AB products of equal length get strictly stronger coverage."""
    net = polymer_network(**GEN3, **LIG, **TEM)
    by_prod = {}
    for ch in _lig_channels(net):
        prod = net.structure[ch.products[0][0]]
        w = max(ch.cat_w, default=0.0)
        by_prod.setdefault(prod, []).append(w)
    ab_max = max(w for s, ws in by_prod.items()
                 if "C" not in s for w in ws)
    c_max = max((w for s, ws in by_prod.items()
                 if "C" in s for w in ws), default=0.0)
    assert ab_max == pytest.approx(100.0)     # capped perfect matches
    assert c_max < ab_max                     # C products capped below


def test_c_positions_always_mismatch():
    net = polymer_network(**GEN3, **LIG, **TEM)
    for ch in _lig_channels(net):
        prod = net.structure[ch.products[0][0]]
        n_c = prod.count("C")
        for k, c in enumerate(ch.catalysts):
            t = net.species[c]
            m = sum(pair_score(t[i], prod[i]) for i in range(len(prod)))
            u = len(prod) - m
            expect = min(np.exp(m * 2.0 - u * 2.0), 100.0)
            assert ch.cat_w[k] == pytest.approx(expect)
            if n_c:
                assert m <= len(prod) - n_c


# ------------------- observables ----------------------------------

def test_ab_purity_and_stem_observables():
    net = polymer_network(**GEN3, **LIG, **TEM)
    idx = np.nonzero(net.polymer_lengths == 4)[0]
    x = np.zeros(net.n_species)
    pure = [i for i in idx if "C" not in net.species[i]]
    junk = [i for i in idx if "C" in net.species[i]]
    x[pure[0]] = 3.0
    x[junk[0]] = 1.0
    out = sequence_entropy(net, x, 4)
    assert out["ab_purity"] == pytest.approx(0.75)
    assert out["mean_stem_pairs"] is not None
    # pure-AB pool only -> purity 1
    x2 = np.zeros(net.n_species)
    x2[pure[0]] = 1.0
    assert sequence_entropy(net, x2, 4)["ab_purity"] == pytest.approx(1.0)


def test_comp_pairs_detected_in_top():
    net = polymer_network(**GEN2, **LIG, **TEM)
    idx = np.nonzero(net.polymer_lengths == 4)[0]
    i_ab = [i for i in idx if net.species[i] == "AABB"][0]
    i_ba = [i for i in idx if net.species[i] == "BBAA"][0]
    x = np.zeros(net.n_species)
    x[i_ab] = 2.0
    x[i_ba] = 2.0
    out = sequence_entropy(net, x, 4)
    assert ["AABB", "BBAA"] in out["comp_pairs_in_top"]


def test_k3_invalid_params_rejected():
    with pytest.raises(ValueError, match="alphabet_size"):
        polymer_network(**{**GEN3, "alphabet_size": 1}, **LIG,
                        templated_ligation=True)
    with pytest.raises(ValueError):
        polymer_network(**GEN3, **LIG, **TEM, delta_g_fold=-1.0)
    with pytest.raises(ValueError):
        polymer_network(**GEN3, **LIG, **TEM, delta_g_loop=-0.5)
