"""Co-criticality: do independently-defined life-observables share λc?

Pipeline:
  1. group records into network FAMILIES (family_id; falls back to
     net_hash). A family = one fixed network content swept over the
     control axis. Networks are never pooled into a single regression:
     different random networks have different intrinsic λc, and pooled
     fits smear family-specific transitions into false crossovers.
  2. per-family critical-point estimation:
       continuous -> fit_all_models model selection per family
       discrete   -> per-family changepoint: earliest λ where O > 0
  3. PAIRED differences D_jk(fam) = λc_j - λc_k over families where
     both are defined. Verdicts:

       "coincident_within_tol" : max |D| <= tol_pair   (small-N exact)
       "coincident"            : TOST — 90% CI of mean D inside
                                 ±tol_pair                            (N large)
       "separated"             : min |D| > tol_pair    (small-N) OR
                                 95% CI of mean D excludes 0 (N large)
       "mixed"                 : some families inside tol, some out
       "inconclusive"          : CI straddles 0 AND ±tol_pair
       "incomparable"          : no family has both λc defined

Observable classes: CONTINUOUS/DISCRETE are fit by default.
STRUCTURAL (O_raf) is a fixed-network-ensemble property — excluded
from the default set but fittable if explicitly requested.
CONDITIONAL (O_select) is only meaningful inside the multistable
regime and is NEVER fit as a global λc curve.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .critical_points import fit_all_models

CONTINUOUS = ("O_maint", "O_persist", "O_assembly")
DISCRETE = ("O_hered",)
CONDITIONAL = ("O_select",)
STRUCTURAL = ("O_raf",)
DEFAULT_OBSERVABLES = CONTINUOUS + DISCRETE

# map selected model -> (params key holding λc, transition class)
_LC_MAP = {"hinge": ("lc", "onset"),
           "itrans": ("lc", "onset_imperfect"),
           "fold": ("lf", "saddle_node"),
           "sigmoid": ("lc", "crossover")}


@dataclass
class LcEstimate:
    lc: float | None
    transition_class: str      # onset/onset_imperfect/saddle_node/
                               # crossover/none/discrete_changepoint
    detail: dict


def group_by_network(records: list[dict]) -> dict[str, list[dict]]:
    """Group records into families for PAIRED analysis.

    Key: rec["family_id"] if present (composite seed|content identity
    written by the campaign runner), else net_hash. A family = one
    fixed network swept over the control axis; sweeps that touch
    generator params produce different net_hash and therefore
    different families — they are never paired.
    """
    by_net: dict[str, list[dict]] = {}
    for r in records:
        key = r.get("family_id") or r["net_hash"]
        by_net.setdefault(key, []).append(r)
    return by_net


group_by_family = group_by_network


def extract_series(by_net: dict[str, list[dict]],
                   obs: str) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """{net_hash: (lam, y)} keeping only records where O_k is defined.

    Requires each record's life_vector (call observables.life_vector
    first and attach as rec["life_vector"]) and a scalar control value
    at rec["control"]["lam"].
    """
    out = {}
    for h, recs in by_net.items():
        lam, y = [], []
        for r in recs:
            lv = r.get("life_vector", {})
            comp = lv.get(obs, {})
            if comp.get("defined") and comp["value"] is not None:
                lam.append(float(r["control"]["lam"]))
                y.append(float(comp["value"]))
        if len(lam) >= 3:
            order = np.argsort(lam)
            out[h] = (np.asarray(lam)[order], np.asarray(y)[order])
    return out


def _clean_aicc(aicc: dict) -> dict:
    return {m: (float(v) if np.isfinite(v) else None)
            for m, v in aicc.items()}


def estimate_lc_continuous(lam: np.ndarray, y: np.ndarray,
                           seed: int = 0) -> LcEstimate:
    fit = fit_all_models(lam, y, seed=seed)
    sel, params = fit["selected"], fit["params"]
    detail = {"selected": sel, "aicc": _clean_aicc(fit["aicc"])}
    if sel in _LC_MAP:
        key, cls = _LC_MAP[sel]
        detail["delta_aicc"] = _clean_aicc(fit["delta_aicc"])
        lc = float(params[key])
        # flag fits that extrapolate a critical point outside the
        # sampled grid — on near-flat data the fitter can place lc far
        # below/above the sweep, which is not a measurement
        if lc < lam.min() or lc > lam.max():
            detail["extrapolated"] = True
        return LcEstimate(lc, cls, detail)
    return LcEstimate(None, "none" if sel == "null" else "crossover",
                      detail)


def estimate_lc_discrete(by_net_series: dict, active_pred,
                         threshold: float = 0.5) -> LcEstimate:
    """Changepoint for a discrete observable (e.g. O_hered > 0).

    λc = smallest λ where the fraction of networks active >=
    `threshold`. Per-point Clopper-Pearson-ish uncertainty reported as
    the λ-range over which the active fraction straddles `threshold`.
    """
    grid = sorted({v for lam, _ in by_net_series.values() for v in lam})
    fracs = []
    for lv in grid:
        act = sum(1 for lam, y in by_net_series.values()
                  for l0, y0 in zip(lam, y) if l0 == lv and active_pred(y0))
        tot = sum(1 for lam, _ in by_net_series.values()
                  for l0 in lam if l0 == lv)
        fracs.append(act / tot if tot else np.nan)
    fracs = np.asarray(fracs)
    frac_dict = {float(g): (float(f) if np.isfinite(f) else None)
                 for g, f in zip(grid, fracs)}
    hits = [grid[i] for i, f in enumerate(fracs)
            if np.isfinite(f) and f >= threshold]
    if not hits:
        return LcEstimate(None, "none",
                          {"active_fraction": frac_dict})
    lc = hits[0]
    # uncertainty: neighboring grid values where frac < threshold
    below = [g for g, f in zip(grid, fracs)
             if g < lc and np.isfinite(f)]
    above = [g for g, f in zip(grid, fracs)
             if g > lc and np.isfinite(f) and f < threshold]
    return LcEstimate(lc, "discrete_changepoint",
                      {"lc_lo": below[-1] if below else None,
                       "lc_hi": above[0] if above else None,
                       "active_fraction": frac_dict})


def _lc_discrete_family(lam: np.ndarray, y: np.ndarray) -> LcEstimate:
    """Per-family changepoint: earliest λ where the value differs from
    its baseline at the smallest swept λ.

    This catches both 0->onset transitions AND structure changes
    (e.g. n_attractors 2->3) that a simple y>0 threshold misses when the
    observable is already active at the grid edge. A family whose value
    never changes is 'always_active' (baseline > 0) or 'none'."""
    order = np.argsort(lam)
    lam, y = np.asarray(lam)[order], np.asarray(y)[order]
    y0 = y[0]
    chg = np.nonzero(y != y0)[0]
    if not len(chg):
        return LcEstimate(
            None, "always_active" if y0 > 0 else "none",
            {"baseline": float(y0)})
    i = int(chg[0])
    return LcEstimate(
        float(lam[i]), "discrete_changepoint",
        {"baseline": float(y0),
         "new_value": float(y[i]),
         "lc_lo": float(lam[i - 1]) if i else None,
         "lc_hi": float(lam[i + 1]) if i + 1 < len(lam) else None})


def lc_by_family(by_net: dict[str, list[dict]], obs: str,
                 seed: int = 0) -> dict[str, LcEstimate]:
    """{family: LcEstimate} — a fit per family, never pooled."""
    if obs in CONDITIONAL:
        return {}
    series = extract_series(by_net, obs)
    out = {}
    for fam, (lam, y) in series.items():
        if obs in DISCRETE:
            out[fam] = _lc_discrete_family(lam, y)
        else:
            out[fam] = estimate_lc_continuous(lam, y, seed=seed)
    return out


def _lc_for_obs(by_net: dict[str, list[dict]], obs: str,
                seed: int) -> LcEstimate:
    """Pooled estimate over families — kept for reference only; pooled
    curves mix heterogeneous networks and are NOT used for verdicts."""
    if obs in CONDITIONAL:
        return LcEstimate(None, "conditional",
                          {"note": "only defined inside the multistable "
                                   "regime; never fit as a global curve"})
    series = extract_series(by_net, obs)
    if not series:
        return LcEstimate(None, "undefined", {})
    if obs in DISCRETE:
        return estimate_lc_discrete(series, lambda v: v > 0)
    lam = np.concatenate([s[0] for s in series.values()])
    y = np.concatenate([s[1] for s in series.values()])
    return estimate_lc_continuous(lam, y, seed=seed)


def _grid_resolution(by_net: dict[str, list[dict]]) -> float | None:
    """Median spacing of the swept control grid across all records."""
    lams = sorted({float(r["control"]["lam"])
                   for recs in by_net.values() for r in recs
                   if "control" in r and r["control"].get("lam") is not None})
    if len(lams) < 2:
        return None
    d = np.diff(np.asarray(lams))
    return float(np.median(d))


def coincidence_report(by_net: dict[str, list[dict]],
                       observables: tuple[str, ...] = DEFAULT_OBSERVABLES,
                       n_boot: int = 400, seed: int = 0,
                       tol_pair: float | None = None,
                       exact_table_max: int = 12) -> dict:
    """Pairwise coincidence matrix from PER-FAMILY paired λc diffs.

    tol_pair: equivalence margin in λ units. None -> 2x median grid
    spacing (below one grid step, coincidence cannot be resolved).
    Small N (<= exact_table_max defined pairs): exact per-family
    difference table, no bootstrap. Large N: bootstrap the mean
    difference over families; TOST (90% CI inside ±tol_pair) =>
    "coincident"; 95% CI excluding 0 => "separated"; else inconclusive.
    """
    rng = np.random.default_rng(seed)
    fams = sorted(by_net)
    n = len(fams)

    grid_res = _grid_resolution(by_net)
    if tol_pair is None:
        tol_pair = 2 * grid_res if grid_res else 0.1

    per_fam = {o: lc_by_family(by_net, o, seed) for o in observables}
    point_ests = {o: _lc_for_obs(by_net, o, seed) for o in observables}

    matrix = {}
    for j, oj in enumerate(observables):
        for ok in observables[j + 1:]:
            diffs = {f: (per_fam[oj][f].lc - per_fam[ok][f].lc)
                     for f in fams
                     if f in per_fam[oj] and f in per_fam[ok]
                     and per_fam[oj][f].lc is not None
                     and per_fam[ok][f].lc is not None
                     and not per_fam[oj][f].detail.get("extrapolated")
                     and not per_fam[ok][f].detail.get("extrapolated")}
            n_pairs = len(diffs)
            entry = {"n_pairs": n_pairs,
                     "paired_diffs": {f: float(d)
                                      for f, d in diffs.items()}}
            if n_pairs == 0:
                entry["verdict"] = "incomparable"
            elif n_pairs <= exact_table_max:
                d = np.asarray(list(diffs.values()))
                mx, mn = float(np.abs(d).max()), float(np.abs(d).min())
                if mx <= tol_pair:
                    verdict = "coincident_within_tol"
                elif mn > tol_pair:
                    verdict = "separated"
                else:
                    verdict = "mixed"
                entry.update({"verdict": verdict,
                              "max_abs_diff": mx,
                              "mean_diff": float(d.mean()),
                              "estimator": "exact_table"})
            else:
                d = np.asarray(list(diffs.values()))
                boot = np.array([
                    d[rng.integers(0, n_pairs, n_pairs)].mean()
                    for _ in range(n_boot)])
                lo90, hi90 = np.quantile(boot, [0.05, 0.95])
                lo95, hi95 = np.quantile(boot, [0.025, 0.975])
                if lo90 >= -tol_pair and hi90 <= tol_pair:
                    verdict = "coincident"          # TOST equivalence
                elif lo95 > 0 or hi95 < 0:
                    verdict = "separated"
                else:
                    verdict = "inconclusive"
                entry.update({"verdict": verdict,
                              "mean_diff": float(d.mean()),
                              "ci90": [float(lo90), float(hi90)],
                              "ci95": [float(lo95), float(hi95)],
                              "estimator": "paired_bootstrap"})
            matrix[f"{oj}|{ok}"] = entry
    return {
        "family_estimates": {
            o: {f: {"lc": e.lc, "class": e.transition_class,
                    **({"extrapolated": True}
                       if e.detail.get("extrapolated") else {})}
                for f, e in per_fam[o].items()}
            for o in observables},
        "point_estimates": {
            o: {"lc": e.lc, "class": e.transition_class,
                "detail": e.detail, "pooled": True}
            for o, e in point_ests.items()},
        "pairwise": matrix,
        "n_families": n, "n_boot": n_boot,
        "tol_pair": tol_pair, "grid_resolution": grid_res,
        "note": ("verdicts use per-family paired differences, not "
                 "pooled fits; 'coincident*' is an equivalence verdict "
                 "inside tol_pair, not failure-to-reject"),
    }


def coincidence_2d(records: list[dict], secondary_param: str,
                   observables: tuple[str, ...] = DEFAULT_OBSERVABLES,
                   n_boot: int = 400, seed: int = 0) -> dict:
    """Perturbational co-criticality over a SECOND swept axis.

    The decisive gated-collinearity test: partition the records by the
    value of `secondary_param` in each record's sweep dict, run the
    per-family paired analysis at every level, and ask whether the
    paired differences D_jk(fam; level) stay within tolerance as the
    second axis varies. Coincidence at a single grid point can be
    faked by gating (one observable's transition merely enables the
    next); a coincidence that PERSISTS across a 2-D scan — transition
    lines that move together or stay pinned — is the signature of a
    shared underlying critical locus.

    Per-pair verdicts:
      "tracking"            : paired diffs within tol at >= 2 levels
                              (and no level violates tol)
      "drift"               : diffs within tol at some level(s) and
                              outside at others — lines separate as
                              the second axis varies
      "separated_at_all_levels": outside tol everywhere defined
      "underdetermined"     : paired diffs defined at < 2 levels

    Also emits `lc_curves` — per-observable, per-family λc(secondary)
    trajectories — so parallel vs diverging transition lines are
    directly inspectable.
    """
    levels: dict[float, list[dict]] = {}
    for r in records:
        v = r.get("sweep", {}).get(secondary_param)
        if v is not None:
            levels.setdefault(float(v), []).append(r)
    if len(levels) < 2:
        return {"status": "insufficient_levels",
                "secondary_param": secondary_param,
                "n_levels": len(levels)}

    per_level: dict[float, dict] = {}
    for k in sorted(levels):
        per_level[k] = coincidence_report(
            group_by_network(levels[k]), observables,
            n_boot=n_boot, seed=seed)

    # union of pair keys across levels
    pairs = sorted({p for rep in per_level.values()
                    for p in rep["pairwise"]})
    matrix2d = {}
    for pair in pairs:
        level_detail = {}
        diffs_at_levels = {}
        for k, rep in per_level.items():
            ent = rep["pairwise"].get(pair)
            if ent is None:
                continue
            diffs = ent.get("paired_diffs") or {}
            diffs_at_levels[k] = diffs
            level_detail[str(k)] = {
                "verdict": ent["verdict"],
                "n_pairs": ent["n_pairs"],
                "tol_pair": rep["tol_pair"],
                "paired_diffs": diffs}
        n_levels_with_pairs = sum(1 for d in diffs_at_levels.values() if d)
        flat = [(k, d, per_level[k]["tol_pair"])
                for k, dd in diffs_at_levels.items()
                for d in dd.values()]
        within = sum(1 for _, d, t in flat if abs(d) <= t)
        if n_levels_with_pairs < 2:
            verdict = "underdetermined"
        elif flat and within == len(flat):
            verdict = "tracking"
        elif within == 0:
            verdict = "separated_at_all_levels"
        else:
            verdict = "drift"
        matrix2d[pair] = {"verdict": verdict,
                          "n_levels_with_pairs": n_levels_with_pairs,
                          "n_diffs_total": len(flat),
                          "n_diffs_within_tol": within,
                          "per_level": level_detail}

    # per-family λc trajectories along the secondary axis
    lc_curves: dict[str, dict[str, list]] = {}
    for o in observables:
        fams: dict[str, list] = {}
        for k in sorted(per_level):
            ests = per_level[k]["family_estimates"].get(o, {})
            for f, e in ests.items():
                if e["lc"] is not None:
                    fams.setdefault(f, []).append(
                        [k, e["lc"], e["class"]])
        if fams:
            lc_curves[o] = fams

    return {"status": "ok",
            "secondary_param": secondary_param,
            "levels": sorted(levels),
            "n_levels": len(levels),
            "pairwise_2d": matrix2d,
            "lc_curves": lc_curves,
            "per_level": {str(k): v for k, v in per_level.items()},
            "note": ("'tracking' = coincidence persists across the "
                     "secondary axis (signature of a shared critical "
                     "locus); 'drift' = lines separate as the second "
                     "axis varies (gated collinearity at best)")}
