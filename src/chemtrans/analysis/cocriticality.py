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

from .critical_points import fit_all_models, lc_from_fit

CONTINUOUS = ("O_maint", "O_persist", "O_assembly")
DISCRETE = ("O_hered",)
CONDITIONAL = ("O_select",)
STRUCTURAL = ("O_raf",)
DEFAULT_OBSERVABLES = CONTINUOUS + DISCRETE

# map selected model -> (params key holding λc, transition class)
# "fold" is handled specially: the model only counts as a saddle-node
# when the fitted jump discontinuity c1-c0 is a real positive jump
# (>= 10% of the data span). c1 ~ c0 makes it a continuous sqrt onset
# (mean-field beta=1/2 second-order), c1 < c0 a downward jump — neither
# is a first-order saddle-node and neither may be labeled one.
_LC_MAP = {"hinge": ("lc", "onset"),
           "itrans": ("lc", "onset_imperfect"),
           "fold": ("lf", "saddle_node"),
           "sigmoid": ("lc", "crossover"),
           "hinge_rev": ("lc", "onset"),
           "itrans_rev": ("lc", "onset_imperfect"),
           "fold_rev": ("lf", "saddle_node")}


def _fold_class(params: dict, y: np.ndarray) -> str:
    """Transition class for a selected 'fold' fit from its jump size.

    The fold model is c0 below lf and c1 + A*sqrt(lam - lf) above, so
    c1 - c0 is the discontinuity at the bifurcation. Only a positive
    jump of >= 10% of the observed data span is a genuine first-order
    saddle-node; a vanished or negative jump is a continuous sqrt
    onset (beta = 1/2) or a downward discontinuity."""
    span = float(y.max() - y.min()) if len(y) else 0.0
    jump = float(params.get("c1", 0.0) - params.get("c0", 0.0))
    if jump >= 0.1 * max(span, 1e-12):
        return "saddle_node"
    return "continuous_half_power"


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


def dropped_lams(by_net: dict[str, list[dict]],
                 obs: str) -> dict[str, list[float]]:
    """{net_hash: [lam]} for records where O_k is undefined — the
    non-certified / UNRESOLVED points a certified-only fit cannot see.
    Their positions relative to the fitted lc matter: drops clustered
    near lc are exactly where survivorship bias bites."""
    out = {}
    for h, recs in by_net.items():
        drops = []
        for r in recs:
            lam = (r.get("control") or {}).get("lam")
            if lam is None:
                continue
            comp = r.get("life_vector", {}).get(obs, {})
            if not (comp.get("defined") and comp["value"] is not None):
                drops.append(float(lam))
        if drops:
            out[h] = sorted(drops)
    return out


def _clean_aicc(aicc: dict) -> dict:
    return {m: (float(v) if np.isfinite(v) else None)
            for m, v in aicc.items()}


def estimate_lc_continuous(lam: np.ndarray, y: np.ndarray,
                           seed: int = 0, n_boot: int = 0,
                           censor_lams=None) -> LcEstimate:
    fit = fit_all_models(lam, y, seed=seed)
    sel, params = fit["selected"], fit["params"]
    detail = {"selected": sel, "aicc": _clean_aicc(fit["aicc"])}
    # Akaike-weighted model average over lc-bearing models: when the
    # winner's margin is thin (delta_aicc < 2) the point lc is a
    # model-selection coin flip — the weighted estimate and its
    # weighted spread make that variance explicit instead of hiding
    # it inside a single model's breakpoint definition.
    delta = {m: v for m, v in fit["delta_aicc"].items()
             if np.isfinite(v)}
    if delta:
        w_raw = {m: np.exp(-0.5 * v) for m, v in delta.items()}
        w_tot = sum(w_raw.values())
        weights = {m: v / w_tot for m, v in w_raw.items()}
        detail["akaike_weights"] = {m: float(v)
                                    for m, v in weights.items()}
        detail["winner_weight"] = float(weights.get(sel, 0.0))
        deltas_sorted = sorted(v for m, v in delta.items() if v > 0)
        detail["runner_up_delta_aicc"] = (
            float(deltas_sorted[0]) if deltas_sorted else None)
        lcs = []
        for m, fm in fit["_fits"].items():
            key = _LC_MAP.get(m)
            if key and np.isfinite(fm.aicc) and np.isfinite(
                    fm.params.get(key[0], np.nan)):
                lcs.append((weights.get(m, 0.0),
                            float(fm.params[key[0]])))
        if lcs:
            wsum = sum(w for w, _ in lcs)
            if wsum > 0:
                mean = sum(w * v for w, v in lcs) / wsum
                var = sum(w * (v - mean) ** 2 for w, v in lcs) / wsum
                detail["lc_aicc_avg"] = float(mean)
                detail["lc_model_spread"] = float(np.sqrt(var))
    if censor_lams:
        # survivorship bracket: a dropped point hides either a
        # pre-transition or a post-transition state. Refit with ALL
        # censored points imputed at each extreme of the family's own
        # range — the two fitted lc's bound the displacement the
        # certified-only estimate conceals. Drops near lc widen the
        # bracket; drops far from it barely move it.
        cl = np.asarray(list(censor_lams), float)
        detail["n_censored"] = len(cl)
        detail["censor_lams"] = [float(v) for v in cl]
        bounds = []
        for fill in (float(y.min()), float(y.max())):
            lam_b = np.concatenate([lam, cl])
            y_b = np.concatenate([y, np.full(len(cl), fill)])
            lcb = lc_from_fit(fit_all_models(lam_b, y_b, seed=seed))
            if lcb is not None:
                bounds.append(lcb)
        if len(bounds) == 2:
            detail["lc_censor_bracket"] = [
                float(min(bounds)), float(max(bounds))]
            detail["censor_robust"] = True
        elif len(bounds) == 1:
            detail["lc_censor_bracket"] = None
            detail["censor_robust"] = False
            detail["censor_note"] = (
                "transition destroyed under extreme imputation")
        else:
            detail["lc_censor_bracket"] = None
            detail["censor_robust"] = False
            detail["censor_note"] = (
                "no critical point under either extreme imputation")
    if sel in _LC_MAP:
        key, cls = _LC_MAP[sel]
        if sel in ("fold", "fold_rev"):
            # the model name is not the physics: verify the jump
            cls = _fold_class(params, y)
            detail["fold_jump"] = float(
                params.get("c1", 0.0) - params.get("c0", 0.0))
        detail["delta_aicc"] = _clean_aicc(fit["delta_aicc"])
        lc = float(params[key])
        # transition direction: monotone models go up ("onset") or
        # down ("offset") across the grid — an emergence and an
        # extinction at the same λ are opposite events, never a
        # shared critical point
        order = np.argsort(lam)
        detail["direction"] = ("onset"
                               if y[order][-1] > y[order][0]
                               else "offset")
        # flag fits that extrapolate a critical point outside the
        # sampled grid — on near-flat data the fitter can place lc far
        # below/above the sweep, which is not a measurement
        if lc < lam.min() or lc > lam.max():
            detail["extrapolated"] = True
        if n_boot > 0:
            detail.update(
                _lc_bootstrap(lam, y, n_boot=n_boot, seed=seed))
        return LcEstimate(lc, cls, detail)
    detail_boot = detail
    if n_boot > 0:
        detail_boot = dict(detail)
        detail_boot.update(_lc_bootstrap(lam, y, n_boot=n_boot,
                                         seed=seed))
    return LcEstimate(None, "none" if sel == "null" else "crossover",
                      detail_boot)


def _lc_bootstrap(lam: np.ndarray, y: np.ndarray, n_boot: int,
                  seed: int) -> dict:
    """Point-resample bootstrap of the fitted critical point.

    Resamples (lam, y) pairs with replacement and refits the full
    candidate set — the resulting distribution mixes point-resample
    variance AND model-selection variance (the selected model can
    change across resamples), which is exactly the uncertainty a
    single-model lc conceals. `n_finite` reports resamples whose
    winner carried an lc — a small fraction means the "critical
    point" is itself unstable across resamples.
    """
    rng = np.random.default_rng(seed)
    n = len(lam)
    lcs, sels = [], []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        f = fit_all_models(lam[idx], y[idx],
                           seed=int(rng.integers(2**31)))
        sels.append(f["selected"])
        key = _LC_MAP.get(f["selected"])
        if key and np.isfinite(f["params"].get(key[0], np.nan)):
            lcs.append(float(f["params"][key[0]]))
    labels, counts = np.unique(sels, return_counts=True)
    out = {"lc_boot_n": n_boot,
           "lc_boot_n_finite": len(lcs),
           "selection_freq": {m: int(c)
                              for m, c in zip(labels, counts)}}
    if len(lcs) >= 10:
        out["lc_ci95"] = [float(v)
                          for v in np.quantile(lcs, [0.025, 0.975])]
        out["lc_boot_median"] = float(np.median(lcs))
    else:
        out["lc_ci95"] = None
        out["lc_boot_median"] = (float(np.median(lcs))
                                 if lcs else None)
    return out


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
    never changes is 'always_active' (baseline > 0) or 'none'.

    The transition direction ("onset" = value increased, "offset" =
    decreased) is recorded in detail["direction"] — an emergence at λ
    and an extinction at λ are NOT the same critical event and must
    never be declared coincident."""
    order = np.argsort(lam)
    lam, y = np.asarray(lam)[order], np.asarray(y)[order]
    y0 = y[0]
    chg = np.nonzero(y != y0)[0]
    if not len(chg):
        return LcEstimate(
            None, "always_active" if y0 > 0 else "none",
            {"baseline": float(y0)})
    i = int(chg[0])
    direction = "onset" if y[i] > y0 else "offset"
    return LcEstimate(
        float(lam[i]), f"discrete_{direction}",
        {"baseline": float(y0),
         "direction": direction,
         "new_value": float(y[i]),
         "lc_lo": float(lam[i - 1]) if i else None,
         "lc_hi": float(lam[i + 1]) if i + 1 < len(lam) else None})


def lc_by_family(by_net: dict[str, list[dict]], obs: str,
                 seed: int = 0, lc_boot: int = 0
                 ) -> dict[str, LcEstimate]:
    """{family: LcEstimate} — a fit per family, never pooled."""
    if obs in CONDITIONAL:
        return {}
    series = extract_series(by_net, obs)
    drops = dropped_lams(by_net, obs)
    out = {}
    for fam, (lam, y) in series.items():
        if obs in DISCRETE:
            out[fam] = _lc_discrete_family(lam, y)
        else:
            out[fam] = estimate_lc_continuous(
                lam, y, seed=seed, n_boot=lc_boot,
                censor_lams=drops.get(fam))
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
                       tol_frac: float = 0.10,
                       exact_table_max: int = 12,
                       min_pairs: int = 3,
                       lc_boot: int = 0) -> dict:
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
    # swept domain across all family records — the equivalence margin
    # is reported against it so "within tol" is auditable as a
    # fraction of the axis, not an invisible constant
    all_lams = [float(r["control"]["lam"]) for recs in
                by_net.values() for r in recs if r.get("control")]
    domain = (max(all_lams) - min(all_lams)) if all_lams else None
    tol_pair_auto = tol_pair is None
    if tol_pair_auto:
        tol_pair = 2 * grid_res if grid_res else 0.1
    if tol_pair_auto and domain and domain > 0:
        # bound the margin to a fraction of the swept domain: a bare
        # 2x-grid-spacing default can span 30-40% of the axis on
        # coarse grids, calling macroscopically separated transitions
        # "coincident_within_tol" (review §3.5). Never below one grid
        # step — finer coincidence is unresolvable on this grid.
        tol_pair = float(min(tol_pair, max(tol_frac * domain,
                                          grid_res or 0.0)))
    tol_frac_of_domain = (tol_pair / domain) if domain else None

    per_fam = {o: lc_by_family(by_net, o, seed, lc_boot=lc_boot)
               for o in observables}
    point_ests = {o: _lc_for_obs(by_net, o, seed) for o in observables}

    # survivorship-bias visibility: count how many records per family
    # never produced a defined observable (non-STATIONARY, truncated,
    # missing x_eq). A λc fit on the survivors hides exactly the
    # critical-region dropouts where the discriminating signal lives.
    coverage: dict = {}
    for o in observables:
        cov = {}
        for h, recs in by_net.items():
            n_def = 0
            for r in recs:
                comp = (r.get("life_vector") or {}).get(o) or {}
                if comp.get("defined") and comp.get("value") is not None:
                    n_def += 1
            cov[h] = {"n_defined": n_def,
                      "n_dropped": len(recs) - n_def}
        coverage[o] = cov

    matrix = {}
    for j, oj in enumerate(observables):
        for ok in observables[j + 1:]:
            paired = [(f, per_fam[oj][f], per_fam[ok][f])
                      for f in fams
                      if f in per_fam[oj] and f in per_fam[ok]
                      and per_fam[oj][f].lc is not None
                      and per_fam[ok][f].lc is not None
                      and not per_fam[oj][f].detail.get("extrapolated")
                      and not per_fam[ok][f].detail.get("extrapolated")]
            # an onset in one observable and an offset in the other at
            # the same λ are opposite events, not a shared critical
            # point — any such family disqualifies a coincidence claim
            dir_mismatch = {
                f for f, ej, ek in paired
                if ej.detail.get("direction") and ek.detail.get("direction")
                and ej.detail["direction"] != ek.detail["direction"]}
            diffs = {f: ej.lc - ek.lc for f, ej, ek in paired
                     if f not in dir_mismatch}
            n_pairs = len(diffs)
            entry = {"n_pairs": n_pairs,
                     "paired_diffs": {f: float(d)
                                      for f, d in diffs.items()},
                     "n_direction_mismatch": len(dir_mismatch)}
            if n_pairs == 0:
                entry["verdict"] = ("discordant_direction"
                                    if dir_mismatch else "incomparable")
            elif n_pairs < min_pairs:
                # below the evidence floor no equivalence or separation
                # claim is defensible — report the diffs but do not
                # pronounce on coincidence
                entry.update({"verdict": "insufficient_families",
                              "n_pairs": n_pairs,
                              "min_pairs": min_pairs,
                              "estimator": "exact_table"})
            elif n_pairs <= exact_table_max:
                d = np.asarray(list(diffs.values()))
                mx, mn = float(np.abs(d).max()), float(np.abs(d).min())
                if mx <= tol_pair:
                    verdict = ("discordant_direction" if dir_mismatch
                               else "coincident_within_tol")
                elif np.all(d > tol_pair) or np.all(d < -tol_pair):
                    # every family's diff clears tol in the SAME
                    # direction — a consistent ordering of transitions
                    verdict = "separated"
                elif mn > tol_pair:
                    # all diffs large but in OPPOSITE directions:
                    # bidirectional scatter, not an ordered cascade
                    verdict = "discordant"
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
                    verdict = ("discordant_direction" if dir_mismatch
                               else "coincident")   # TOST equivalence
                elif lo95 > 0 or hi95 < 0:
                    verdict = "separated"
                else:
                    verdict = "inconclusive"
                # bootstrap p-values for multiple-testing correction:
                # TOST equivalence p = max tail beyond ±tol; two-sided
                # separation p = 2 * min tail across 0
                p_tost = float(max(
                    np.mean(boot <= -tol_pair),
                    np.mean(boot >= tol_pair)))
                p_sep = float(min(1.0, 2.0 * min(
                    np.mean(boot <= 0.0), np.mean(boot >= 0.0))))
                entry.update({"verdict": verdict,
                              "mean_diff": float(d.mean()),
                              "ci90": [float(lo90), float(hi90)],
                              "ci95": [float(lo95), float(hi95)],
                              "p_tost": p_tost, "p_separation": p_sep,
                              "estimator": "paired_bootstrap"})
            matrix[f"{oj}|{ok}"] = entry
    # Benjamini–Hochberg across the cells that carry bootstrap
    # p-values. Exact-table cells (n <= exact_table_max) are
    # deterministic classifications, not hypothesis tests — they are
    # annotated rather than adjusted.
    for p_key, q_key in (("p_tost", "q_tost"),
                         ("p_separation", "q_separation")):
        tested = sorted(
            ((e[p_key], k) for k, e in matrix.items()
             if p_key in e and np.isfinite(e[p_key])),
            key=lambda t: t[0])
        m = len(tested)
        if m:
            # enforce monotonicity: q_(i) = min over j>=i of p_(j)*m/j
            adj = [min(1.0, p * m / rank)
                   for rank, (p, _) in enumerate(tested, start=1)]
            for i in range(m - 2, -1, -1):
                adj[i] = min(adj[i], adj[i + 1])
            for (p, k), q in zip(tested, adj):
                matrix[k][q_key] = float(q)
    n_tested = sum(1 for e in matrix.values() if "p_tost" in e)
    return {
        "family_estimates": {
            o: {f: {"lc": e.lc, "class": e.transition_class,
                    **({"extrapolated": True}
                       if e.detail.get("extrapolated") else {}),
                    **({"n_censored": e.detail["n_censored"],
                        "censor_lams": e.detail["censor_lams"],
                        "censor_robust": e.detail.get("censor_robust", True),
                        **({"lc_censor_bracket": e.detail["lc_censor_bracket"]}
                           if e.detail.get("lc_censor_bracket") else {}),
                        **({"censor_note": e.detail["censor_note"]}
                           if e.detail.get("censor_note") else {})}
                       if e.detail.get("n_censored") else {})}
                for f, e in per_fam[o].items()}
            for o in observables},
        "point_estimates": {
            o: {"lc": e.lc, "class": e.transition_class,
                "detail": e.detail, "pooled": True}
            for o, e in point_ests.items()},
        "pairwise": matrix,
        "n_pairwise_tests": len(matrix),
        "n_bootstrap_tested": n_tested,
        "multiple_testing": (
            "Benjamini-Hochberg FDR applied across the "
            f"{n_tested} bootstrap-tested cells (q_tost / "
            "q_separation); exact-table cells are deterministic "
            "classifications, not hypothesis tests, and are not "
            "adjusted"),
        "coverage": coverage,
        "n_families": n, "n_boot": n_boot, "min_pairs": min_pairs,
        "tol_pair": tol_pair, "grid_resolution": grid_res,
        "tol_pair_frac_of_domain": tol_frac_of_domain,
        "swept_domain": domain,
        "note": ("verdicts use per-family paired differences, not "
                 "pooled fits; 'coincident*' is an equivalence verdict "
                 "inside tol_pair, not failure-to-reject"),
    }


def coincidence_2d(records: list[dict], secondary_param: str,
                   observables: tuple[str, ...] = DEFAULT_OBSERVABLES,
                   n_boot: int = 400, seed: int = 0,
                   min_pairs: int = 3,
                   tol_frac: float = 0.10,
                   lc_boot: int = 0) -> dict:
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
            n_boot=n_boot, seed=seed, min_pairs=min_pairs,
            tol_frac=tol_frac, lc_boot=lc_boot)

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
        # "tracking" must be evidenced WITHIN a family: count families
        # whose paired diff is defined at >= 2 secondary levels.
        # Without one, within-tol agreement could be fabricated by
        # pooling different families at different levels (fam A at
        # level 1, fam B at level 2) — a shared locus no family
        # actually exhibited.
        fams_seen = sorted({f for dd in diffs_at_levels.values()
                            for f in dd})
        n_fam_multilevel = 0
        for fam in fams_seen:
            n_lv = sum(1 for k in diffs_at_levels if fam in
                       diffs_at_levels[k])
            if n_lv >= 2:
                n_fam_multilevel += 1
        # a tracking/drift verdict needs at least TWO families actually
        # followed across >= 2 levels — with one, within-tol agreement
        # is a single-network observation, not an ensemble result
        if n_fam_multilevel < 2:
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
                          "n_families_multilevel": n_fam_multilevel,
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
