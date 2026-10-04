"""Experiment 1 — the base lambda sweep (phase diagram).

For a grid of control-parameter values, generate an ensemble of networks,
simulate each to steady state, and record all order parameters.

Output schema (one jsonl row per run):
    run_id, generator, generator_params, seed,
    lambda_components {p_cat, feed_rate, dilution_rate},
    rho_K, phi, A, S, H, E, C_mean, C_tail_mass, assembly_histogram,
    steady_state_reached, notes

This dataset is the input to analysis.critical_points and the
Model A/B/C discrimination test.
"""

if __name__ == "__main__":
    raise SystemExit("not implemented — first milestone")
