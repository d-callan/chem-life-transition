"""Experiment 1 — the base lambda sweep (phase diagram).

For a grid of control-parameter values, generate an ensemble of networks,
simulate each to steady state, and record all order parameters.

Output schema (one jsonl row per run):
    run_id, generator, generator_params, seed,
    lambda {p_cat, feed_rate, dilution_rate},
    rho_K, phi, A, assembly {mean, tail_mass, hist},
    n_core_species, steady_state_reached, runaway

This dataset is the input to analysis.critical_points and the
Model A/B/C discrimination test.

Usage:
    python experiments/sweep_lambda.py --out data/sweep.jsonl --plot data/phi.png
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from chemtrans.measure import measure_run
from chemtrans.networks.generators import random_crn


def run_one(p_cat, feed_rate, dilution_rate, args, seed):
    net = random_crn(
        n_species=args.n_species,
        n_reactions=args.n_reactions,
        p_cat=p_cat,
        n_food=args.n_food,
        max_reactants=args.max_reactants,
        p_cat_required=args.p_cat_required, k_sat=args.k_sat,
        seed=seed,
    )
    row = {
        "generator": "random_crn",
        "generator_params": {
            "n_species": args.n_species,
            "n_reactions": args.n_reactions,
            "n_food": args.n_food,
            "max_reactants": args.max_reactants,
            "p_cat_required": args.p_cat_required,
            "k_sat": args.k_sat,
        },
        "seed": seed,
        "lambda": {"p_cat": p_cat, "feed_rate": feed_rate,
                   "dilution_rate": dilution_rate},
    }
    row.update(measure_run(net, feed_rate, dilution_rate, t_max=args.t_max))
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/sweep.jsonl")
    ap.add_argument("--plot", default=None, help="optional PNG path for Phi(p_cat)")
    ap.add_argument("--n-species", type=int, default=60)
    ap.add_argument("--n-reactions", type=int, default=120)
    ap.add_argument("--n-food", type=int, default=5)
    ap.add_argument("--max-reactants", type=int, default=2)
    ap.add_argument("--n-reps", type=int, default=8)
    ap.add_argument("--p-cat", type=float, nargs="+",
                    default=np.linspace(0.0, 0.5, 11).tolist())
    ap.add_argument("--p-cat-required", type=float, default=0.5)
    ap.add_argument("--k-sat", type=float, default=None,
                    help="Michaelis-Menten saturation conc; None = unbounded mass-action")
    ap.add_argument("--feed-rate", type=float, default=1.0)
    ap.add_argument("--dilution-rate", type=float, default=0.1)
    ap.add_argument("--t-max", type=float, default=300.0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    n_done = 0
    total = len(args.p_cat) * args.n_reps
    with out.open("w") as f:
        for p_cat in args.p_cat:
            for rep in range(args.n_reps):
                seed = args.seed + rep * 1009 + hash(round(p_cat, 6)) % 997
                row = run_one(p_cat, args.feed_rate, args.dilution_rate,
                              args, seed)
                row["run_id"] = f"s{seed}_p{p_cat:.3f}"
                row.pop("x_ss", None)  # don't write full state to jsonl
                f.write(json.dumps(row) + "\n")
                n_done += 1
                print(f"\r{n_done}/{total}  p_cat={p_cat:.3f} rep={rep}  "
                      f"phi={row.get('phi', float('nan')):.3f} "
                      f"rho_K={row['rho_K']:.2f}", end="", flush=True)
    print(f"\nwrote {out}")

    if args.plot:
        plot_phi(args)


def plot_phi(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    df = pd.read_json(args.out, lines=True)
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    for axi, col, label in zip(
        ax, ["phi", "A", "rho_K"],
        ["organizational closure $\\Phi$", "autocatalytic activity $A$",
         "spectral radius $\\rho(K)$"],
    ):
        axi.scatter(df["lambda"].str["p_cat"], df[col], alpha=0.4, s=12)
        g = df.groupby(df["lambda"].str["p_cat"])[col].mean()
        axi.plot(g.index, g.values, "k-", lw=2)
        axi.set_xlabel("$p_{cat}$"); axi.set_ylabel(label)
    fig.tight_layout()
    fig.savefig(args.plot, dpi=150)
    print(f"wrote {args.plot}")


if __name__ == "__main__":
    main()
