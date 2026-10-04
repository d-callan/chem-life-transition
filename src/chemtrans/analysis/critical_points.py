"""Critical-point detection and cross-observable alignment.

THE core analysis for the Model A / B / C discrimination experiment.

For each observable O(lambda) measured on an ensemble of networks:

1. Estimate lambda_c per observable — e.g. peak of susceptibility
   chi = N * Var(O), peak of |dO/dlambda|, or bifurcation fitting.
2. Get a DISTRIBUTION of lambda_c estimates across ensemble replicates
   and bootstrap resamples.

Discrimination:

- Model B (one transition): all observables' lambda_c distributions
  overlap; susceptibilities peak together; order parameters are
  cross-correlated in the critical region.
- Model A (cascade): lambda_c distributions are staggered and
  statistically separable.
- Model C (nested): one primary lambda_c (organizational closure) plus
  secondary sub-thresholds detectable ONLY within the organized
  sub-ensemble — the "phases within a phase" signature. Test by
  conditioning: among runs that crossed the primary threshold, do
  secondary observables show their own critical points?
"""

import numpy as np
import pandas as pd


def susceptibility(values: np.ndarray, n_samples: int) -> np.ndarray:
    """chi = N * Var(O) across ensemble at each lambda — peaks at lambda_c."""
    raise NotImplementedError


def estimate_critical_point(
    lambdas: np.ndarray, observable: np.ndarray, method: str = "susceptibility"
) -> tuple[float, float]:
    """Return (lambda_c, uncertainty) for one observable."""
    raise NotImplementedError


def align_critical_points(
    df: pd.DataFrame, observables: list[str], n_bootstrap: int = 200
) -> pd.DataFrame:
    """Bootstrap lambda_c distributions per observable; return alignment table.

    Output feeds the Model A/B/C decision: pairwise overlap of lambda_c
    distributions + conditional sub-threshold analysis for Model C.
    """
    raise NotImplementedError
