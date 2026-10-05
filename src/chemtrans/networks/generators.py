"""Polymer network generator for the chemtrans 2.0 model.

Baseline chemistry: monomer-addition polymerization driven by activated
monomers (see plan §5.4). Species:

  letters  A, B, ...                (monomers, food)
  A*, B*, ...                       (activated monomers = fuel, food)
  W                                 (waste / leaving group, food)
  strings of length >= 2            (polymers; only these may catalyze)

Channels (all reversible, rates from species free energies g_i):

  elongation     P + A* <=> PA + W   driven; dG = g_bond - dmu; may be
                                     catalyzed (k_f/k_r multiplied by
                                     (eps + kappa*x_C) — same factor both
                                     directions, so K_eq is untouched).
  thermal        PA <=> P + A        spontaneous condensation / hydrolysis;
                                     dG = -g_bond; eps-scaled.
  deactivation   A* <=> A + W        fuel decay; dG = -dmu; NOT eps-scaled.

Detailed balance: dG_elong - dG_deact + dG_hyd = 0 identically, so rate
ratios derived from g satisfy Wegscheider cycle consistency by
construction. `eps` scales the UNcatalyzed rate of polymer-synthesis
channels (elongation + thermal); at eps=0 the polymer-free state is
absorbing (catalysts are non-food polymers only).

Food inlet (concentrations, strictly mass-conserving across sweeps):
  [A*]_inlet = chi * F_inlet / |Sigma|
  [A]_inlet  = (1-chi) * F_inlet / |Sigma|
  [W]_inlet = [P]_inlet = 0        (total letter influx = F_inlet always)

Catalyst assignment (per elongation channel, latent-score design so that
p_cat / p_selfrep / p_info sweeps nest rather than reshuffle):
  u < p_cat            -> channel is catalyzed
  v < p_selfrep        -> catalyst is the product PA itself
  else w < p_info      -> catalyst shares sequence with PA (sub/superstring)
  else                 -> random polymer of length >= min_cat_len
"""

from __future__ import annotations

import itertools

import numpy as np

from .core import Channel, ReactionNetwork

GENERATOR_VERSION = "polymer_v2"


def _all_strings(alphabet: list[str], max_len: int) -> list[str]:
    out: list[str] = []
    for L in range(1, max_len + 1):
        out.extend("".join(t) for t in itertools.product(alphabet, repeat=L))
    return out


def polymer_network(
    alphabet_size: int = 2,
    max_len: int = 4,
    *,
    p_cat: float = 0.3,
    p_selfrep: float = 0.3,
    p_info: float = 0.7,
    min_cat_len: int = 2,
    delta_mu: float = 4.0,
    chi: float = 1.0,
    f_inlet: float = 1.0,
    eps: float = 1.0,
    kappa: float = 1.0,
    g_bond: float = 1.0,
    k_lig: float = 1.0,
    k_hyd: float = 0.1,
    k_deact: float = 0.1,
    seed: int | None = None,
) -> ReactionNetwork:
    """Build a conserved, reversible, energy-driven polymer network.

    All species strings up to `max_len` over `alphabet_size` letters are in
    the pool; every elongation of a shorter string is a channel, so the
    pool is closed by construction (no silent truncation).

    Rates (dimensionless free energies in k_BT, standard concentration = 1):
      g(monomer) = 0, g(A*) = delta_mu, g(W) = 0,
      g(P) = g_bond * (len(P) - 1).
      k_r is derived: k_r = k_f * exp(dG)  (k_f/k_r = exp(-dG) = K_eq).
    """
    if alphabet_size < 1 or max_len < 2:
        raise ValueError("need alphabet_size>=1 and max_len>=2")
    if delta_mu < 0:
        raise ValueError("delta_mu must be >= 0")
    if not (0.0 <= chi <= 1.0):
        raise ValueError("chi must be in [0,1]")

    alphabet = [chr(65 + i) for i in range(alphabet_size)]
    letters = list(alphabet)
    act = [a + "*" for a in alphabet]
    strings = [s for s in _all_strings(alphabet, max_len)]

    species: list[str] = letters + act + ["W"] + strings[alphabet_size:]
    index = {s: i for i, s in enumerate(species)}
    n = len(species)

    structure = {i: s for s, i in index.items()
                 if "*" not in s and s != "W"}
    # structure for activated monomer counts its letter for `comp`
    comp = np.zeros((alphabet_size, n))
    for i, s in enumerate(species):
        if s == "W":
            continue
        base = s.replace("*", "")
        for ch in base:
            comp[alphabet.index(ch), i] += 1

    g = np.zeros(n)
    g[index["W"]] = 0.0
    for a in alphabet:
        g[index[a]] = 0.0
        g[index[a + "*"]] = delta_mu
    for s in strings[alphabet_size:]:
        g[index[s]] = g_bond * (len(s) - 1)

    food = {index[a] for a in alphabet} | {index[a + "*"] for a in alphabet}
    food.add(index["W"])
    food_inlet = np.zeros(n)
    for a in alphabet:
        food_inlet[index[a + "*"]] = chi * f_inlet / alphabet_size
        food_inlet[index[a]] = (1.0 - chi) * f_inlet / alphabet_size

    # deterministic per-channel latent scores: seed -> one rng; arrays drawn
    # in a fixed order so p_cat/p_selfrep/p_info sweeps nest.
    rng = np.random.default_rng(seed)
    n_elong = sum(
        alphabet_size for L in range(1, max_len)
        for _ in range(alphabet_size ** L)
    )
    u_cat = rng.random(n_elong)
    u_self = rng.random(n_elong)
    u_info = rng.random(n_elong)
    u_kin = rng.random(n_elong)

    candidates = [
        index[s]
        for s in strings
        if len(s) >= min_cat_len
    ]
    cand_strings = [s for s in strings if len(s) >= min_cat_len]

    channels: list[Channel] = []
    e = 0
    for L in range(1, max_len):
        for P in itertools.product(alphabet, repeat=L):
            P = "".join(P)
            iP = index[P]
            for A in alphabet:
                PA = P + A
                iPA, iA, iW = index[PA], index[A + "*"], index["W"]
                # --- driven elongation P + A* <=> PA + W -----------------
                dG = g[iPA] + g[iW] - g[iP] - g[iA]  # = g_bond - delta_mu
                k_f = k_lig
                k_r = k_lig * np.exp(dG)
                cats: tuple[int, ...] = ()
                if u_cat[e] < p_cat and candidates:
                    if u_self[e] < p_selfrep:
                        cats = (iPA,)
                    elif u_info[e] < p_info:
                        kin = [
                            index[s] for s in cand_strings
                            if s != PA and (PA in s or s in PA)
                        ]
                        if kin:
                            cats = (kin[int(u_kin[e] * len(kin))],)
                        else:
                            cats = (
                                candidates[int(u_kin[e] * len(candidates))],
                            )
                    else:
                        cats = (candidates[int(u_kin[e] * len(candidates))],)
                channels.append(Channel(
                    reactants=tuple(sorted(((iP, 1), (iA, 1)))),
                    products=tuple(sorted(((iPA, 1), (iW, 1)))),
                    catalysts=cats,
                    k_f=k_f, k_r=k_r,
                    synth=True, mechanism="elongation",
                ))
                e += 1

                # --- thermal condensation / hydrolysis PA <=> P + A -----
                # P may equal A (e.g. AA <=> A + A): merge stoichiometry.
                prods: dict[int, int] = {}
                for i in (iP, index[A]):
                    prods[i] = prods.get(i, 0) + 1
                dG3 = sum(v * g[i] for i, v in prods.items()) - g[iPA]
                channels.append(Channel(
                    reactants=((iPA, 1),),
                    products=tuple(sorted(prods.items())),
                    catalysts=(),
                    k_f=k_hyd, k_r=k_hyd * np.exp(dG3),
                    synth=True, mechanism="thermal",
                ))

    for a in alphabet:
        iA_, iA0, iW = index[a + "*"], index[a], index["W"]
        dG2 = g[iA0] + g[iW] - g[iA_]  # = -delta_mu
        channels.append(Channel(
            reactants=((iA_, 1),),
            products=tuple(sorted(((iA0, 1), (iW, 1)))),
            catalysts=(),
            k_f=k_deact, k_r=k_deact * np.exp(dG2),
            synth=False, mechanism="deactivation",
        ))

    return ReactionNetwork(
        species,
        channels,
        food,
        structure=structure,
        comp=comp,
        g=g,
        food_inlet=food_inlet,
        kappa=kappa,
        eps=eps,
        params=dict(
            generator=GENERATOR_VERSION,
            alphabet_size=alphabet_size, max_len=max_len,
            p_cat=p_cat, p_selfrep=p_selfrep, p_info=p_info,
            min_cat_len=min_cat_len,
            delta_mu=delta_mu, chi=chi, f_inlet=f_inlet, eps=eps,
            kappa=kappa, g_bond=g_bond,
            k_lig=k_lig, k_hyd=k_hyd, k_deact=k_deact,
            seed=seed,
        ),
    )
