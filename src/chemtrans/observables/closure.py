"""Formal food-generated RAF (maxRAF) — Hordijk-Steel fixed point.

Definition (per the standard RAF algorithm):

  R' <- set of catalyzed channels (channels with >= 1 catalyst assigned)
  repeat:
      cl <- closure of food F under R', ignoring catalysts
            (grow: a reaction fires if all reactants are in cl)
      R' <- { r in R' : reactants(r) subset cl
                        AND catalysts(r) intersect cl != {} }
  until fixpoint; the unique maximal RAF is the result (empty allowed).

Catalysts are alternatives: ANY one catalyst in the closure suffices.
Food species are inside the closure by construction, so a food-catalyzed
reaction is always eligible (the baseline generator never creates them).

This module also reports a *separate* kinetic-accessibility closure
(which channels could run at the given eps and catalyst availability) and
an explicitly-named internally-catalyzed subset — never conflated with
the RAF.
"""

from __future__ import annotations

import numpy as np

from ..networks.core import ReactionNetwork


def _closure(net: ReactionNetwork, chan_set: set[int]) -> set[int]:
    """Food-generated species closure under chan_set, ignoring catalysts."""
    avail = set(net.food)
    changed = True
    while changed:
        changed = False
        for r in chan_set:
            rxn = net.channels[r]
            if all(i in avail for i, _ in rxn.reactants):
                new = {i for i, _ in rxn.products} - avail
                if new:
                    avail |= new
                    changed = True
    return avail


def raf_reactions(net: ReactionNetwork) -> set[int]:
    """Indices of channels in the unique maximal RAF (may be empty)."""
    R = {r for r, ch in enumerate(net.channels) if ch.catalysts}
    while True:
        cl = _closure(net, R)
        pruned = {
            r
            for r in R
            if all(i in cl for i, _ in net.channels[r].reactants)
            and any(c in cl for c in net.channels[r].catalysts)
        }
        if pruned == R:
            return R
        R = pruned


def raf_species(net: ReactionNetwork) -> set[int]:
    """Non-food species participating in the maxRAF."""
    out: set[int] = set()
    for r in raf_reactions(net):
        ch = net.channels[r]
        out.update(i for i, _ in ch.products)
        out.update(i for i, _ in ch.reactants)
        out.update(ch.catalysts)
    return out - net.food


def kinetic_accessible(net: ReactionNetwork) -> set[int]:
    """Channels that can fire at all: reactants reachable from food under
    eps-scaled rates and catalyst availability. Uses a food-generated
    reachability closure; kinetic accessibility != RAF membership.

    At eps > 0 every channel fires once its reactants are reachable
    (the eps background term supplies uncatalyzed synthesis). At
    eps == 0 a channel needs a reachable catalyst only when its FORWARD
    direction carries the eps gate (synth_dir in {"fwd", "both"});
    non-synth channels and decay-forward channels (deactivation,
    thermal hydrolysis with synth_dir="rev" — bg = 1.0 on the forward
    direction) always fire once reactants are present.
    """
    avail = set(net.food)
    accessible: set[int] = set()
    changed = True
    while changed:
        changed = False
        for r, ch in enumerate(net.channels):
            if r in accessible:
                continue
            if ch.k_f <= 0 and ch.k_r <= 0:
                continue
            if not all(i in avail for i, _ in ch.reactants):
                continue
            if net.eps == 0 and ch.synth and ch.synth_dir != "rev" \
                    and not any(c in avail for c in ch.catalysts):
                continue
            accessible.add(r)
            new = {i for i, _ in ch.products} - avail
            if new:
                avail |= new
                changed = True
    return accessible


def phi_decomposition(
    net: ReactionNetwork, x: np.ndarray, eps: float | None = None
) -> dict:
    """Channel-derived flux decomposition at state x.

    Returns dict with forward/reverse/net flux inside vs outside the RAF:
      mass_raf   primary order parameter: monomer mass in maxRAF species /
                 total non-food monomer mass
      phi_gross  RAF-channel |flux| / total |flux| (structural diagnostic)
      phi_cat    catalyzed flux inside RAF channels / total |flux|
      phi_net    net forward RAF flux / total forward flux (diagnostic only)
    All four are returned with explicit denominators; phi_net can exceed
    [0,1] and is labeled diagnostic.
    """
    x = np.asarray(x, dtype=float)
    f, b = net.rates(x, eps)
    fwd, rev = np.abs(f), np.abs(b)
    net_flux = f - b
    core = raf_reactions(net)
    core_idx = np.array(sorted(core), dtype=int)

    # which flux is catalyzed vs background? decompose boost factors:
    # f_r = k_f*mass*(bg + kappa*C x^p) -> catalyzed part uses the same
    # cat_power nonlinearity as rates() (x^p, not x)
    Rf, Rr, _ = net._factors()
    kf, kr, _gate_f, _gate_r = net._rate_vectors()
    x_app = np.append(x, 1.0)
    mf = np.prod(x_app[Rf], axis=1)
    mr = np.prod(x_app[Rr], axis=1)
    cat_amp = net.kappa * np.asarray(
        net.catalysis_matrix() @ np.power(x, net.cat_power)).ravel()
    cat_f = kf * mf * cat_amp
    cat_b = kr * mr * cat_amp

    tot_f = fwd.sum()
    tot_abs = fwd.sum() + rev.sum()
    denom = {
        "tot_forward": float(tot_f),
        "tot_absflux": float(tot_abs),
        "tot_catflux": float(cat_f.sum() + cat_b.sum()),
    }
    out = {
        "n_raf_channels": int(len(core)),
        "phi_gross": (float((fwd + rev)[core_idx].sum()) / tot_abs
                      if tot_abs > 0 and len(core) else 0.0),
        "phi_cat": (float((cat_f + cat_b)[core_idx].sum()) / tot_abs
                    if tot_abs > 0 and len(core) else 0.0),
        "phi_net": (float(net_flux[core_idx].sum()) / tot_f
                    if tot_f > 0 and len(core) else 0.0),
        "mass_raf": _mass_raf(net, x),
        "denominators": denom,
    }
    return out


def _mass_raf(net: ReactionNetwork, x: np.ndarray) -> float:
    """Monomer mass inside maxRAF species / total non-food monomer mass."""
    if net.comp is None:
        raise ValueError("mass_raf requires a composition matrix")
    comp_mass = net.comp.sum(axis=0)         # monomer mass per species
    nonfood = ~net.is_food
    total = float(comp_mass[nonfood] @ x[nonfood])
    if total <= 0:
        return 0.0
    core_s = sorted(raf_species(net))
    if not core_s:
        return 0.0
    idx = np.array(core_s, dtype=int)
    return float(comp_mass[idx] @ x[idx] / total)
