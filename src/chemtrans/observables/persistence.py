"""First-passage / survival analysis with honest censoring.

Conventions:

- Event times are extracted by `first_passage_times` (dwell-conditioned:
  a trajectory counts as having entered the organized region only if it
  stays at or above `theta_hi` for `dwell` consecutive recorded time).
- Observations are (time, event_flag): censored trajectories carry
  their censoring time, never imputed event times.
- Kaplan-Meier and Nelson-Aalen are standard right-censored estimators.
- RMST integrates the KM curve to a fixed horizon t_star (never beyond
  the data).
- Scaling fits (ln tau ~ V etc.) are the caller's responsibility; this
  module reports estimators, not claimed exponents.
"""

from __future__ import annotations

import numpy as np


def region_state(t: np.ndarray, o: np.ndarray,
                 theta_hi: float, theta_lo: float) -> np.ndarray:
    """Hysteretic in/out state of the organized region.

    in-region turns ON when o >= theta_hi and OFF when o < theta_lo;
    inside the hysteretic band the previous state persists.
    """
    t = np.asarray(t, float)
    o = np.asarray(o, float)
    state = np.zeros(len(t), dtype=bool)
    inside = False
    for k in range(len(t)):
        if o[k] >= theta_hi:
            inside = True
        elif o[k] < theta_lo:
            inside = False
        state[k] = inside
    return state


def first_passage_times(
    series: list[tuple[np.ndarray, np.ndarray]],
    theta_hi: float,
    dwell: float,
    theta_lo: float | None = None,
) -> list[tuple[float, bool]]:
    """Dwell-conditioned first-entry times into the organized region.

    Entry is credited when the hysteretic region state stays ON for
    `dwell` units of time continuously. Returns [(t_entry, True)] for
    observed entries, [(t_censor, False)] for censored trajectories.
    """
    if theta_lo is None:
        theta_lo = theta_hi
    out = []
    for t, o in series:
        t = np.asarray(t, float)
        state = region_state(t, o, theta_hi, theta_lo)
        entry = None
        on_time = None
        for k in range(len(t)):
            if state[k]:
                if on_time is None:
                    on_time = t[k]
                if t[k] - on_time >= dwell:
                    entry = t[k]
                    break
            else:
                on_time = None
        out.append((entry if entry is not None else float(t[-1]),
                    entry is not None))
    return out


def kaplan_meier(times: np.ndarray, events: np.ndarray):
    """KM survival curve. Returns (t_knots, S, se_greenwood)."""
    times = np.asarray(times, float)
    events = np.asarray(events, bool)
    order = np.argsort(times)
    times, events = times[order], events[order]
    t_knots, S, se = [0.0], [1.0], [0.0]
    surv, var = 1.0, 0.0
    uniq = np.unique(times)
    for tt in uniq:
        at_risk = np.sum(times >= tt)
        d = np.sum(events & (times == tt))
        if at_risk == 0:
            continue
        surv *= 1.0 - d / at_risk
        if at_risk - d > 0:
            var += d / (at_risk * (at_risk - d))
        t_knots.append(float(tt))
        S.append(surv)
        se.append(surv * np.sqrt(var))
    return np.array(t_knots), np.array(S), np.array(se)


def nelson_aalen(times: np.ndarray, events: np.ndarray):
    """Cumulative hazard estimate. Returns (t_knots, H)."""
    times = np.asarray(times, float)
    events = np.asarray(events, bool)
    t_knots, H = [0.0], [0.0]
    h = 0.0
    for tt in np.unique(times):
        at_risk = np.sum(times >= tt)
        d = np.sum(events & (times == tt))
        if at_risk == 0:
            continue
        h += d / at_risk
        t_knots.append(float(tt))
        H.append(h)
    return np.array(t_knots), np.array(H)


def rmst(times: np.ndarray, events: np.ndarray,
         t_star: float) -> float:
    """Restricted mean survival time: integral of KM S(t) on [0, t_star].

    t_star must not exceed the largest observed time; caller chooses it
    (typically min over trajectories of censor/event times).
    """
    if t_star <= 0:
        return 0.0
    t_knots, S, _ = kaplan_meier(times, events)
    t_ext = np.append(np.clip(t_knots, 0, t_star), t_star)
    t_ext = np.unique(t_ext)
    S_at = np.empty_like(t_ext)
    for i, tt in enumerate(t_ext):
        j = np.searchsorted(t_knots, tt, side="right") - 1
        S_at[i] = S[max(j, 0)]
    return float(np.trapezoid(S_at, t_ext))


def survival_table(obs: list[tuple[float, bool]]) -> dict:
    """Bundle times/events arrays + counts for run records."""
    times = np.array([o[0] for o in obs])
    events = np.array([o[1] for o in obs])
    return {
        "times": times,
        "events": events,
        "n": int(len(times)),
        "n_events": int(events.sum()),
        "n_censored": int((~events).sum()),
        "max_observed_time": float(times.max()) if len(times) else 0.0,
    }


def extinction_scale_fit(volumes: np.ndarray, tau: np.ndarray) -> dict:
    """Fit ln(tau_ext) = a + b*V (mean-field exponential survival scaling).

    For background-driven NUCLEATION times the expected trend is
    decreasing in V (more chances to nucleate); that direction must be
    fit with tau ~ 1/V or a decreasing function, NOT ln tau ~ V. This
    function is for extinction/persistence times only.
    """
    volumes = np.asarray(volumes, float)
    tau = np.asarray(tau, float)
    if len(volumes) < 3:
        return {"status": "insufficient_points", "n": int(len(volumes))}
    y = np.log(tau)
    A = np.polyfit(volumes, y, 1)
    pred = np.polyval(A, volumes)
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return {"status": "ok", "slope_V": float(A[0]),
            "intercept": float(A[1]),
            "r2": 1.0 - ss_res / max(ss_tot, 1e-300),
            "n": int(len(volumes))}


def nucleation_scale_fit(volumes: np.ndarray, tau: np.ndarray) -> dict:
    """Fit tau_nuc ~ C / V^g for background-driven nucleation."""
    volumes = np.asarray(volumes, float)
    tau = np.asarray(tau, float)
    if len(volumes) < 3:
        return {"status": "insufficient_points", "n": int(len(volumes))}
    y = np.log(tau)
    x = np.log(volumes)
    A = np.polyfit(x, y, 1)
    pred = np.polyval(A, x)
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return {"status": "ok", "gamma": float(-A[0]),
            "log_C": float(A[1]),
            "r2": 1.0 - ss_res / max(ss_tot, 1e-300),
            "n": int(len(volumes))}
