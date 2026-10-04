"""Critical-point detection and cross-observable alignment.

THE core analysis for the Model A / B / C discrimination experiment.

For each observable O(lambda) measured on an ensemble of networks:

1. Estimate lambda_c per observable — peak of ensemble susceptibility
   chi = Var(O) across replicates, or peak of |d<O>/dlambda|.
2. Bootstrap over replicates to get a DISTRIBUTION of lambda_c estimates.

Discrimination:

- Model B (one transition): all observables' lambda_c CIs overlap.
- Model A (cascade): lambda_c distributions staggered and separable.
- Model C (nested): one primary lambda_c plus secondary thresholds
  detectable only within the organized sub-ensemble — "phases within a
  phase". Tested by conditioning: among runs above the primary lambda_c,
  do secondary observables have their own critical points?
"""

import numpy as np
import pandas as pd


def get_observable(df: pd.DataFrame, name: str) -> np.ndarray:
    """Extract an observable column, supporting dotted paths into dicts."""
    if name in df.columns:
        return df[name].astype(float).to_numpy()
    head, *rest = name.split(".")
    return df[head].map(lambda d: _dig(d, rest)).astype(float).to_numpy()


def _dig(d, keys):
    for k in keys:
        d = d.get(k) if isinstance(d, dict) else None
        if d is None:
            return np.nan
    return d


def susceptibility(matrix: np.ndarray) -> np.ndarray:
    """chi(lambda) = variance across ensemble replicates at each lambda.

    matrix: (n_lambda, n_reps). NaN-safe; returns nan where <3 samples.
    """
    return np.array([
        np.nanvar(row) if np.isfinite(row).sum() >= 3 else np.nan
        for row in matrix
    ])


def estimate_critical_point(
    lambdas: np.ndarray,
    matrix: np.ndarray,
    method: str = "susceptibility",
) -> float:
    """lambda_c for one observable.

    method="susceptibility": argmax of ensemble variance (peaks where the
    ensemble is split between phases — works for percolation-type and
    fluctuation-driven transitions alike).
    method="gradient": argmax of |d<O>/dlambda| on the ensemble mean.
    """
    if method == "susceptibility":
        crit = susceptibility(matrix)
    elif method == "gradient":
        mean = np.array([
            np.nanmean(row) if np.isfinite(row).any() else np.nan
            for row in matrix
        ])
        crit = np.abs(np.gradient(mean, lambdas))
    else:
        raise ValueError(method)
    if np.all(~np.isfinite(crit)):
        return np.nan
    return float(lambdas[int(np.nanargmax(crit))])


def align_critical_points(
    df: pd.DataFrame,
    observables: list[str],
    lambda_col: str = "p_cat",
    n_bootstrap: int = 200,
    seed: int = 0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Bootstrap lambda_c distributions per observable; compare alignment.

    df: one row per run, containing lambda_col and the observable columns
    (dotted paths allowed, e.g. "assembly.tail_mass").

    Returns (estimates, pairwise): estimates has lc/ci_lo/ci_hi per
    observable; pairwise[i,j]["aligned"] = fraction of bootstrap lambda_c
    pairs within one grid step (high = co-located = Model B signature;
    low = staggered = Model A signature).
    """
    rng = np.random.default_rng(seed)
    lambdas_all = np.sort(df[lambda_col].unique())
    boots: dict[str, list[float]] = {o: [] for o in observables}

    mats = {}
    for o in observables:
        vals = get_observable(df, o)
        mat = np.full((len(lambdas_all), 0), np.nan)
        cols = []
        for lam in lambdas_all:
            cols.append(vals[(df[lambda_col] == lam).to_numpy()])
        mats[o] = lambdas_all, cols

    def mat_from(cols):
        m = max(len(c) for c in cols)
        out = np.full((len(cols), m), np.nan)
        for i, c in enumerate(cols):
            out[i, : len(c)] = c
        return out

    for o in observables:
        lambdas, cols = mats[o]
        full = mat_from(cols)
        for _ in range(n_bootstrap):
            resampled = np.array([
                rng.choice(col[np.isfinite(col)], size=len(col), replace=True)
                if np.isfinite(col).any() else np.full(len(col), np.nan)
                for col in cols
            ], dtype=object)
            boot_mat = mat_from([np.asarray(c, float) for c in resampled])
            lc = estimate_critical_point(lambdas, boot_mat)
            if np.isfinite(lc):
                boots[o].append(lc)
        # point estimate on the real data
        boots[o].insert(0, estimate_critical_point(lambdas, full))

    est_rows = []
    for o in observables:
        b = np.asarray(boots[o])
        point, dist = b[0], b[1:]
        est_rows.append({
            "observable": o,
            "lambda_c": point,
            "ci_lo": float(np.nanpercentile(dist, 2.5)) if len(dist) else np.nan,
            "ci_hi": float(np.nanpercentile(dist, 97.5)) if len(dist) else np.nan,
            "lc_std": float(np.nanstd(dist)) if len(dist) else np.nan,
        })
    estimates = pd.DataFrame(est_rows)

    step = float(np.median(np.diff(lambdas_all))) if len(lambdas_all) > 1 else 1.0
    pairs = []
    for i, oi in enumerate(observables):
        for oj in observables[i + 1:]:
            bi, bj = np.asarray(boots[oi][1:]), np.asarray(boots[oj][1:])
            bi, bj = bi[np.isfinite(bi)], bj[np.isfinite(bj)]
            if len(bi) and len(bj):
                # aligned: bootstrap lambda_c draws within one grid step
                aligned = float(np.mean(
                    np.abs(bi[:, None] - bj[None, :]) <= 1.5 * step))
                pairs.append({"obs_i": oi, "obs_j": oj,
                              "aligned": aligned})
    pairwise = pd.DataFrame(pairs)
    return estimates, pairwise
