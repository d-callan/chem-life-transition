"""Experiment 4 — THE headline experiment: Model A vs B vs C.

Runs its own ensemble sweep (with persistence observables enabled) and
asks how many independent critical points the five order parameters have:

  A(λ) autocatalytic activity, S(λ) shock survival, R0(λ) reproductive
  number, E(λ) evolvability proxy, C(λ) assembly tail mass, Phi(λ) closure.

Decision logic (from analysis.critical_points):
- bootstrap lambda_c CIs all overlap            -> Model B (one transition)
- consistent staggered ordering in lambda_c     -> Model A (cascade)
- primary lambda_c + secondary critical points
  inside the organized-phase subset only        -> Model C (nested)

Caveat: H (heritable information) is NOT measurable in random_crn —
there is no information carrier. It needs the polymer_ligation generator
(milestone 4). E is a mutational-sensitivity proxy, not true selection.

Usage:
    python experiments/transition_discrimination.py --out data/discrim.jsonl
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from chemtrans.analysis.critical_points import align_critical_points
from chemtrans.measure import measure_run
from chemtrans.networks.generators import random_crn

OBSERVABLES = ["phi", "A", "S", "R0", "E",
               "assembly.tail_mass", "n_core_species"]


def run_ensemble(args):
    rows = []
    total = len(args.p_cat) * args.n_reps
    done = 0
    for p_cat in args.p_cat:
        for rep in range(args.n_reps):
            seed = args.seed + rep * 1009 + hash(round(p_cat, 6)) % 997
            net = random_crn(
                n_species=args.n_species, n_reactions=args.n_reactions,
                p_cat=p_cat, n_food=args.n_food,
                max_reactants=args.max_reactants,
                p_cat_required=args.p_cat_required, k_sat=args.k_sat,
                exclusive_catalysis=args.exclusive_catalysis,
                seed=seed,
            )
            row = {"p_cat": p_cat, "seed": seed,
                   "run_id": f"s{seed}_p{p_cat:.3f}",
                   "exclusive_catalysis": args.exclusive_catalysis}
            row.update(measure_run(
                net, args.feed_rate, args.dilution_rate,
                t_max=args.t_max, with_persistence=True,
                persistence_seed=seed,
            ))
            rows.append(row)
            done += 1
            print(f"\r{done}/{total} p_cat={p_cat:.2f} rep={rep} "
                  f"phi={row.get('phi', float('nan')):.3f}",
                  end="", flush=True)
    print()
    return rows


def discriminate(df: pd.DataFrame, args) -> dict:
    est, pairs = align_critical_points(df, OBSERVABLES,
                                       n_bootstrap=args.n_bootstrap)
    print("\n=== lambda_c estimates (bootstrap 95% CI) ===")
    print(est.to_string(index=False))
    print("\n=== pairwise alignment (≈1 = co-located critical point) ===")
    print(pairs.to_string(index=False) if len(pairs) else "n/a")

    # primary transition = lambda_c of Phi (the organizational order param)
    lc_phi = est.loc[est.observable == "phi", "lambda_c"].iloc[0]

    # Model C test: within organized-phase runs, do secondary order
    # parameters retain internal critical structure?
    organized = df[df["phi"].fillna(0) > args.phi_threshold]
    nested = {}
    if len(organized) >= args.n_reps * 3:
        sub_est, _ = align_critical_points(
            organized, [o for o in OBSERVABLES if o != "phi"],
            n_bootstrap=args.n_bootstrap)
        nested = sub_est.set_index("observable")["lambda_c"].to_dict()
        print("\n=== Model C check: critical pts within organized phase ===")
        print(sub_est.to_string(index=False))
    else:
        print("\nModel C check skipped: too few organized-phase runs")

    # verdict heuristic
    aligned = pairs["aligned"] if len(pairs) else pd.Series([1.0])
    frac_aligned = float(np.mean(aligned > 0.5))
    verdict = {
        "lambda_c_primary": lc_phi,
        "frac_pairs_aligned": frac_aligned,
        "nested_critical_points": nested,
    }
    if frac_aligned < 0.3:
        verdict["model"] = "A (staggered cascade)"
    elif nested:
        verdict["model"] = "C (nested) — coupled primary + internal thresholds"
    else:
        verdict["model"] = "B (single coupled transition)"
    return verdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/discrimination.jsonl")
    ap.add_argument("--n-species", type=int, default=40)
    ap.add_argument("--n-reactions", type=int, default=200)
    ap.add_argument("--n-food", type=int, default=8)
    ap.add_argument("--max-reactants", type=int, default=1)
    ap.add_argument("--p-cat", type=float, nargs="+",
                    default=np.linspace(0.1, 0.4, 7).tolist())
    ap.add_argument("--p-cat-required", type=float, default=0.5)
    ap.add_argument("--k-sat", type=float, default=None,
                    help="Michaelis-Menten saturation conc; None = unbounded mass-action")
    ap.add_argument("--exclusive-catalysis", action="store_true")
    ap.add_argument("--n-reps", type=int, default=12)
    ap.add_argument("--feed-rate", type=float, default=1.0)
    ap.add_argument("--dilution-rate", type=float, default=0.1)
    ap.add_argument("--t-max", type=float, default=150.0)
    ap.add_argument("--phi-threshold", type=float, default=0.1)
    ap.add_argument("--n-bootstrap", type=int, default=100)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    rows = run_ensemble(args)
    with out.open("w") as f:
        for r in rows:
            r2 = {k: v for k, v in r.items() if k != "x_ss"}
            f.write(json.dumps(r2) + "\n")
    print(f"wrote {out}")

    df = pd.DataFrame(rows)
    verdict = discriminate(df, args)
    print("\n=== VERDICT ===")
    print(json.dumps(verdict, indent=2, default=float))


if __name__ == "__main__":
    main()
