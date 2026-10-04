"""Experiment 3 — finite-size scaling.

A real phase transition sharpens with system size; a crossover doesn't.
Sweep p_cat at several network sizes with FIXED reaction density
(n_reactions = density * n_species) and check:

- does the transition region narrow as N grows? (P(Phi>thresh) -> step)
- does Var(Phi) in the critical region grow with N? (susceptibility)
- do lambda_c estimates at different N converge to one value?

Negative result is also informative: if nothing sharpens, the "transition"
is a generator artifact, not physics.

Usage:
    python experiments/finite_size.py --out data/finite_size.jsonl
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from chemtrans.measure import measure_run
from chemtrans.networks.generators import random_crn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/finite_size.jsonl")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--sizes", type=int, nargs="+", default=[20, 40, 80])
    ap.add_argument("--density", type=float, default=5.0,
                    help="n_reactions / n_species, held fixed across sizes")
    ap.add_argument("--n-food", type=int, default=8)
    ap.add_argument("--max-reactants", type=int, default=1)
    ap.add_argument("--p-cat", type=float, nargs="+",
                    default=np.linspace(0.1, 0.4, 7).tolist())
    ap.add_argument("--p-cat-required", type=float, default=0.5)
    ap.add_argument("--n-reps", type=int, default=10)
    ap.add_argument("--feed-rate", type=float, default=1.0)
    ap.add_argument("--dilution-rate", type=float, default=0.1)
    ap.add_argument("--t-max", type=float, default=150.0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    total = len(args.sizes) * len(args.p_cat) * args.n_reps
    done = 0
    with out.open("w") as f:
        for n_species in args.sizes:
            n_reactions = int(round(args.density * n_species))
            for p_cat in args.p_cat:
                for rep in range(args.n_reps):
                    seed = (args.seed + rep * 1009
                            + n_species * 7919 + int(p_cat * 1000))
                    net = random_crn(
                        n_species=n_species, n_reactions=n_reactions,
                        p_cat=p_cat, n_food=min(args.n_food, n_species // 4),
                        max_reactants=args.max_reactants,
                        p_cat_required=args.p_cat_required, seed=seed,
                    )
                    row = {"n_species": n_species, "n_reactions": n_reactions,
                           "p_cat": p_cat, "seed": seed,
                           "run_id": f"N{n_species}_s{seed}_p{p_cat:.3f}"}
                    row.update(measure_run(net, args.feed_rate,
                                           args.dilution_rate,
                                           t_max=args.t_max))
                    row.pop("x_ss", None)
                    f.write(json.dumps(row) + "\n")
                    done += 1
                    print(f"\r{done}/{total} N={n_species} "
                          f"p_cat={p_cat:.2f} rep={rep} "
                          f"phi={row.get('phi', float('nan')):.3f}",
                          end="", flush=True)
    print(f"\nwrote {out}")

    if args.plot:
        plot(args)


def plot(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    df = pd.read_json(args.out, lines=True)
    df["organized"] = df["phi"].fillna(0) > 0.1
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for n, g in df.groupby("n_species"):
        p = g.groupby("p_cat")["organized"].mean()
        v = g.groupby("p_cat")["phi"].var()
        ax[0].plot(p.index, p.values, "o-", label=f"N={n}")
        ax[1].plot(v.index, v.values, "o-", label=f"N={n}")
    ax[0].set_ylabel("P(organized)"); ax[0].set_xlabel("$p_{cat}$")
    ax[1].set_ylabel("Var($\\Phi$)"); ax[1].set_xlabel("$p_{cat}$")
    ax[0].legend(); ax[1].legend()
    fig.tight_layout(); fig.savefig(args.plot, dpi=150)
    print(f"wrote {args.plot}")


if __name__ == "__main__":
    main()
