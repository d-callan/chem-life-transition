"""Sequence-space sparsity / pruning observables.

Quantifies whether polymer mass at a given length is distributed
uniformly across sequence isomers (unguided abiotic soup) or
concentrated on a few sequences (selection / heritable lineages).

For length L over alphabet size K there are K^L sequences.
    p_i   = x_i / sum_{j: len(j)=L} x_j
    H_L   = -sum_i p_i log2 p_i          (bits; max = L*log2 K)
    N_eff = 2^{H_L}                      (effective populated count)
    O_sparsity = 1 - H_L / (L*log2 K)    (0 = flat, 1 = single seq)

Empty pools return sparsity=None (undefined, not zero) so callers
never mistake 'no polymers of this length' for 'perfect selection'.
"""
from __future__ import annotations

import numpy as np


def _len_pool(net, length: int) -> np.ndarray:
    lens = net.polymer_lengths
    return np.nonzero(lens == length)[0]


MASS_FLOOR = 1e-9


def sequence_entropy(net, x: np.ndarray, length: int,
                     top_k: int = 5) -> dict:
    """Shannon entropy, effective diversity, and sparsity for one
    polymer-length shell at concentration/state vector x.

    Negative solver excursions are clipped to 0 and reported via
    `neg_mass_abs` (a measurement diagnostic — this is an observable,
    not the vector field). A pool whose positive mass is numerically
    absent reports sparsity=None, never a fabricated 0 or 1."""
    x = np.asarray(x, dtype=float)
    idx = _len_pool(net, length)
    h_max = float(np.log2(max(len(idx), 1)))
    out = {"length": length, "n_sequences": int(len(idx)),
           "entropy": 0.0, "entropy_max": h_max,
           "n_eff": 1.0, "sparsity": None, "top_sequences": [],
           "pool_mass": 0.0, "neg_mass_abs": 0.0,
           "ab_purity": None, "mean_stem_pairs": None,
           "comp_pairs_in_top": []}
    if len(idx) == 0:
        return out
    xi = x[idx]
    out["neg_mass_abs"] = float(-xi[xi < 0].sum())
    xi = np.maximum(xi, 0.0)
    tot = float(xi.sum())
    out["pool_mass"] = tot
    # numerically absent shell, or solver excursions so large the
    # clipped vector is no longer a meaningful distribution — either
    # way the observable is undefined, not "perfect selection"
    if tot <= MASS_FLOOR or out["neg_mass_abs"] > 0.5 * tot:
        return out
    p = xi / tot
    nz = p[p > 0]
    h = float(-np.sum(nz * np.log2(nz)))
    pos_order = np.argsort(-xi)
    out["entropy"] = h
    out["n_eff"] = float(2.0 ** h)
    out["sparsity"] = (1.0 - h / h_max) if h_max > 0 else None
    out["top_sequences"] = [
        {"species": net.species[idx[pos]],
         "mass_fraction": float(xi[pos] / tot),
         "concentration": float(xi[pos])}
        for pos in pos_order[:top_k]
    ]
    # alphabet purification: fraction of pool mass in pure {A,B}
    # sequences (non-pairing letters like C are 'abiotic junk')
    ab_mask = np.array([
        all(ch in "AB" for ch in net.species[i]) for i in idx])
    if ab_mask.any():
        out["ab_purity"] = float(xi[ab_mask].sum() / tot)
        # folding enrichment: mass-weighted hairpin stem-pair count
        stems = np.array([_stem_pairs(net.species[i]) for i in idx])
        out["mean_stem_pairs"] = float((p * stems).sum())
        # reciprocal-pair co-enrichment: complementary pairs within
        # the top-k mass leaders (candidate replicating hypercycles)
        top_s = [t["species"] for t in out["top_sequences"]]
        out["comp_pairs_in_top"] = [
            [a, b] for a in top_s for b in top_s
            if a < b and _is_complement(a, b)]
    return out


def _is_complement(a: str, b: str) -> bool:
    """Pointwise Watson-Crick complement: every position pairs.
    Requires both strings pure {A,B} — any C is always mismatch."""
    if len(a) != len(b):
        return False
    return all((x == "A" and y == "B") or (x == "B" and y == "A")
               for x, y in zip(a, b))


def _stem_pairs(s: str, l_loop_min: int = 2) -> int:
    """Intramolecular hairpin stem matches: s[k] vs s[L-1-k]."""
    L = len(s)
    l_stem = (L - l_loop_min) // 2
    return sum(1 for k in range(l_stem)
               if (s[k] == "A" and s[L - 1 - k] == "B")
               or (s[k] == "B" and s[L - 1 - k] == "A"))


def sparsity_report(net, x: np.ndarray,
                    lengths: tuple = (2, 3, 4, 5, 6),
                    top_k: int = 5) -> dict:
    """Per-length sparsity summary plus pooled max-length metrics."""
    rep = {}
    for L in lengths:
        rep[f"L{L}"] = sequence_entropy(net, x, L, top_k=top_k)
    return rep
