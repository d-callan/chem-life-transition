"""Exact stochastic simulation (Gillespie SSA) for chemtrans 2.0.

Processes (expanded from channels, single source of truth):

  per channel r with participant multiset alpha (fwd) / beta (rev):
    background fwd:  k = k_f * (eps if fwd direction gated else 1)
    background rev:  k = k_r * (eps if rev direction gated else 1)
    (synth_dir selects which direction carries the eps gate — a thermal
    PA <=> P + A channel with synth_dir="rev" keeps hydrolysis running
    at eps=0 while condensation stays gated)
    catalyzed fwd by c in catalysts: k = k_f * kappa,  participants alpha+[c]
    catalyzed rev by c:              k = k_r * kappa,  participants beta+[c]

  Hazard convention (matches ODE k * prod x^{a}):
      h(n; V) = k * V^{1-|a|} * prod_i (n_i)_{a_i}
      with (n)_q the falling factorial. A hazard is exactly 0 whenever
      n_i < a_i for any participant i (impossible events cannot fire).

  Boundary processes (CSTR):
      feed on food species i: h = D * food_inlet_i * V   (delta +1)
      dilution on species i:  h = D * n_i                (delta -1)

Recording is non-anticipatory: record times strictly before an event's
time are filled with the PRE-event state; the post-event state applies
from the event time onward. Censoring (max_events) is reported as
`truncated` with the actual censor time — never padded to the horizon.

Observers: callables `obs(t, n, kind, proc_idx)` invoked AFTER each event
(and once at t=0). kind in {"init", "channel", "feed", "dilution"}; use
for online first-passage/dwell statistics without storing full paths.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..networks.core import ReactionNetwork


@dataclass(frozen=True)
class Process:
    """One elementary stochastic process."""
    rate: float                    # stochastic rate constant k
    participants: tuple[int, ...]  # flat species list (with repetition)
    offsets: tuple[int, ...]       # per-position count of same species earlier
    delta: tuple[tuple[int, int], ...]  # (species, dcount) net stoich
    kind: str                      # "channel" | "feed" | "dilution"
    channel: int                   # source channel index or -1
    # functional-boost catalysis (hazard_mode="functional"): h is
    # multiplied by (n_s / V)^power — the exact SSA analogue of the
    # ODE's phenomenological kappa * x_c^p rate law, rather than the
    # default literal p-molecule falling-factorial mechanism. None for
    # "participants" mode and all non-catalyzed processes.
    func_boost: tuple[int, int] | None = None
    # product inhibition / site saturation (net.k_inhibit > 0): h is
    # divided by (1 + K * sum_{c in channel catalysts} n_c / V) — the
    # SSA analogue of the ODE's shared per-channel denominator.
    # (catalyst_indices, K); None when uninhibited.
    func_inhibit: tuple[tuple[int, ...], float] | None = None


def _mk_offsets(parts: list[int]) -> tuple[int, ...]:
    seen: dict[int, int] = {}
    out = []
    for s in parts:
        out.append(seen.get(s, 0))
        seen[s] = seen.get(s, 0) + 1
    return tuple(out)


def build_processes(net: ReactionNetwork, dilution: float,
                    eps: float | None = None,
                    hazard_mode: str = "participants",
                    polymer_retention: float = 0.0,
                    retention_len_scale: float = 0.0) -> list[Process]:
    """Expand channels + boundary processes for the SSA.

    hazard_mode selects the cooperative-catalysis convention:
      "participants" (default) — cat_power catalyst COPIES are literal
        reactant participants; hazard uses falling factorials (n_c)_p,
        identically zero at n_c < p. This is the mechanistically exact
        interpretation of cooperative catalysis.
      "functional" — the ODE's phenomenological boost kappa * (n_c/V)^p
        multiplies the rate; a single catalyst copy produces nonzero
        catalytic hazard for any n_c > 0, matching the deterministic
        rate law pointwise. Use to isolate whether a stochastic claim
        (small-V extinction, nucleation barriers) is physics or an
        artifact of the (n_c)_p zero at n_c <= 1.

    polymer_retention is the semi-permeable-compartment fraction:
    structural polymer (len>=2) dilution events run at rate
    D*(1-rho); food/feed species at D. Catalyst hazards carry the
    network's cat_weights (length-scaled potency), matching the ODE
    rate law.
    """
    if hazard_mode not in ("participants", "functional"):
        raise ValueError(
            f"hazard_mode must be 'participants'|'functional', "
            f"got {hazard_mode!r}")
    if not (0.0 <= polymer_retention < 1.0):
        raise ValueError("polymer_retention must be in [0, 1)")
    if not (0.0 <= retention_len_scale <= 1.0):
        raise ValueError("retention_len_scale must be in [0, 1]")
    if eps is None:
        eps = net.eps
    procs: list[Process] = []
    nu_cache = {}
    for r, ch in enumerate(net.channels):
        nu = ch.nu
        fwd_delta = tuple(sorted(nu.items()))
        rev_delta = tuple(sorted((i, -a) for i, a in nu.items()))
        bg_f = eps if (ch.synth and ch.synth_dir != "rev") else 1.0
        bg_r = eps if (ch.synth and ch.synth_dir != "fwd") else 1.0
        for direction, parts, delta, k in (
                ("fwd", ch.factors("fwd"), fwd_delta, ch.k_f * bg_f),
                ("rev", ch.factors("rev"), rev_delta, ch.k_r * bg_r)):
            if k > 0:
                procs.append(Process(
                    rate=k, participants=parts, offsets=_mk_offsets(list(parts)),
                    delta=delta, kind="channel", channel=r))
        for c in ch.catalysts:
            if hazard_mode == "participants":
                # cooperative catalysis: cat_power copies of c participate
                cat_parts = (c,) * net.cat_power
                boost = None
            else:
                cat_parts = ()
                boost = (c, net.cat_power)
            # length-weighted potency folds into the rate constant:
            # identical multiplier on both directions, matching the
            # ODE's symmetric catalysis (detailed balance untouched)
            kcat = net.kappa * float(net.cat_weights[c])
            # shared per-channel inhibition denominator (product
            # inhibition / site saturation) — same divisor in both
            # hazard modes, matching the ODE boost/(1+K*sum x_c)
            finh = ((ch.catalysts, net.k_inhibit)
                    if net.k_inhibit > 0 else None)
            for direction, parts, delta, k in (
                    ("fwd", ch.factors("fwd") + cat_parts, fwd_delta,
                     ch.k_f * kcat),
                    ("rev", ch.factors("rev") + cat_parts, rev_delta,
                     ch.k_r * kcat)):
                if k > 0:
                    procs.append(Process(
                        rate=k, participants=parts,
                        offsets=_mk_offsets(list(parts)), delta=delta,
                        kind="channel", channel=r,
                        func_boost=boost, func_inhibit=finh))
    if dilution > 0:
        for i in sorted(net.food | net.inlet_species):
            if net.food_inlet[i] > 0:
                procs.append(Process(
                    rate=dilution * net.food_inlet[i], participants=(),
                    offsets=(), delta=((i, 1),), kind="feed", channel=-1))
        out_rates = net.dilution_outflow(
            dilution, polymer_retention, retention_len_scale)
        for i in range(net.n_species):
            procs.append(Process(
                rate=float(out_rates[i]), participants=(i,),
                offsets=(0,), delta=((i, -1),), kind="dilution",
                channel=-1))
    return procs


def hazard(p: Process, n: np.ndarray, volume: float) -> float:
    """h = k * V^{1-m} * prod_p (n_{i_p} - offset_p); zero if infeasible.

    func_boost processes additionally multiply by (n_c/V)^power —
    positive for any n_c > 0, zero at n_c = 0 (unlike the participants
    falling factorial, which also vanishes at 0 < n_c < p).
    """
    h = p.rate
    m = len(p.participants)
    for pos, s in enumerate(p.participants):
        f = n[s] - p.offsets[pos]
        if f <= 0:
            return 0.0
        h *= f
    if m > 0:
        h *= volume ** (1 - m)
    else:
        h *= volume  # zeroth-order influx: conc/time * V
    if p.func_boost is not None:
        c, pw = p.func_boost
        if n[c] <= 0:
            return 0.0
        h *= (n[c] / volume) ** pw
    if p.func_inhibit is not None:
        cats, k_inh = p.func_inhibit
        h /= 1.0 + k_inh * float(
            n[np.asarray(cats, dtype=int)].sum()) / volume
    return h


def counts_from_concentration(
    x: np.ndarray, volume: float, policy: str = "round",
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Explicit concentration -> integer count policy (recorded in runs)."""
    x = np.asarray(x, dtype=float)
    if policy == "round":
        return np.maximum(np.rint(x * volume), 0).astype(np.int64)
    if policy == "poisson":
        if rng is None:
            raise ValueError("poisson init requires rng")
        return rng.poisson(np.maximum(x, 0.0) * volume).astype(np.int64)
    raise ValueError(f"unknown init policy {policy!r}")


@dataclass
class SSAResult:
    t: np.ndarray                  # record times (non-anticipatory)
    n: np.ndarray                  # counts at record times
    events: int
    completed: bool                # reached t_max with time remaining
    truncated: bool                # hit max_events
    censor_time: float             # actual last simulated time
    final_counts: np.ndarray
    net_hash: str
    seed_key: str
    diagnostics: dict = field(default_factory=dict)


def simulate(
    net: ReactionNetwork,
    t_max: float,
    dilution: float = 0.1,
    volume: float = 100.0,
    n0: np.ndarray | None = None,
    eps: float | None = None,
    seed: int | None = None,
    max_events: int = 5_000_000,
    record_dt: float | None = None,
    observers: tuple = (),
    hazard_mode: str = "participants",
    polymer_retention: float = 0.0,
    retention_len_scale: float = 0.0,
) -> SSAResult:
    """Run the exact SSA. See module docstring for conventions and
    `build_processes` for the hazard_mode catalysis convention."""
    rng = np.random.default_rng(seed)
    if eps is None:
        eps = net.eps
    procs = build_processes(net, dilution, eps, hazard_mode=hazard_mode,
                            polymer_retention=polymer_retention,
                            retention_len_scale=retention_len_scale)
    n_species = net.n_species

    if n0 is None:
        n = np.zeros(n_species, dtype=np.int64)
        for i in net.food | net.inlet_species:
            n[i] = int(round(net.food_inlet[i] * volume))
    else:
        n = np.asarray(n0, dtype=np.int64).copy()
        if (n < 0).any() or n.shape != (n_species,):
            raise ValueError("n0 must be nonnegative integer counts")

    if record_dt is None:
        record_dt = t_max / 2000.0
    n_records = int(np.floor(t_max / record_dt)) + 1
    rec_times = np.arange(n_records) * record_dt
    X = np.empty((n_records, n_species), dtype=np.int64)
    filled = 0

    def fill_upto(t_event: float, state: np.ndarray):
        """Fill record times strictly < t_event with the pre-event state."""
        nonlocal filled
        while filled < n_records and rec_times[filled] < t_event:
            X[filled] = state
            filled += 1

    t = 0.0
    for obs in observers:
        obs(0.0, n, "init", -1)
    fill_upto(0.0 + np.finfo(float).eps, n)  # t=0 record

    events = 0
    truncated = False
    while t < t_max:
        h = np.array([hazard(p, n, volume) for p in procs])
        total = h.sum()
        if total <= 0.0:
            break  # absorbing state: nothing can fire
        dt = rng.exponential(1.0 / total)
        if t + dt > t_max:
            t = t_max
            break
        # fill record times < event time with PRE-event state
        fill_upto(t + dt, n)
        t += dt
        u = rng.random() * total
        idx = int(np.searchsorted(np.cumsum(h), u, side="right"))
        idx = min(idx, len(procs) - 1)
        p = procs[idx]
        for i, dcount in p.delta:
            n[i] += dcount
        if (n < 0).any():
            raise RuntimeError(
                f"impossible event produced negative counts "
                f"(process kind={p.kind} channel={p.channel})")
        events += 1
        for obs in observers:
            obs(t, n, p.kind, idx)
        if events >= max_events:
            truncated = True
            break

    # fill remaining records through horizon with final constant state
    while filled < n_records:
        X[filled] = n
        filled += 1

    return SSAResult(
        t=rec_times, n=X, events=events,
        completed=not truncated, truncated=truncated,
        censor_time=t, final_counts=n.copy(),
        net_hash=net.hash, seed_key=str(seed),
        diagnostics={"dilution": dilution, "volume": volume, "eps": eps,
                     "hazard_mode": hazard_mode,
                     "max_events": max_events, "record_dt": record_dt},
    )
