"""Experiment 4 — THE headline experiment: Model A vs B vs C.

Uses the sweep_lambda dataset (or runs a denser one) and applies
analysis.critical_points.align_critical_points to all five order
parameters: A, S, H, E, C.

Decision logic:
- lambda_c distributions for all observables overlap   -> Model B
  (single underlying transition; the five observables are coupled
  manifestations of one critical point)
- lambda_c distributions staggered and separable       -> Model A
  (cascade of independent thresholds)
- one primary lambda_c + sub-thresholds appearing only
  within the post-transition ensemble                   -> Model C
  (nested: secondary transitions inside the organized phase)

Model C test detail: split runs at the primary lambda_c (from Phi),
then rerun critical-point estimation for S, H, E, C on the
organized-phase subset only. Secondary critical points there are the
"phases within a phase" signature.

Also worth running across GENERATORS (random_crn vs polymer_ligation):
if the alignment pattern is the same in chemically different ensembles,
that's the universality evidence the conversation ends on.
"""

if __name__ == "__main__":
    raise SystemExit("not implemented")
