"""Experiment 6 — Question 0: is the chemistry->organization shift a
bifurcation (phase transition) or a smooth crossover?

Three diagnostics along a 1-D control-parameter slice:

1. SPECTRAL GAP CLOSING: alpha = max Re(eig(J)) at the measured steady
   state. A stable attractor has alpha < 0; alpha -> 0 as lambda -> lc is
   critical slowing down — the signature of an approaching bifurcation.
   tau ~ -1/alpha is the relaxation time.
2. ENSEMBLE SUSCEPTIBILITY + BIMODALITY: Var(phi) and Sarle's bimodality
   coefficient BC across reps. A transition shows a var peak and BC > 5/9
   in the critical region; a crossover stays smooth.
3. WITHIN-NETWORK BISTABILITY: for networks near the boundary, simulate
   unseeded (food only) vs seeded (x_ss of an organized run as the seed)
   — different outcomes = a real dynamical bifurcation in that network,
   not just ensemble percolation.

Usage:
    python experiments/criticality.py --generator polymer_ligation \
        --lambda-param p_info --p-cat-fixed 0.6
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from chemtrans.dynamics.ode import simulate
from chemtrans.measure import measure_run
from chemtrans.networks.generators import polymer_ligation, random_crn
from chemtrans.observables.closure import organizational_closure
from chemtrans.observables.spectral import (
    invasion_pressure,
    jacobian_spectral_abscissa,
)


def bimodality_coefficient(vals) -> float:
    """Sarle's BC = (skew^2 + 1) / kurtosis; > 5/9 suggests bimodality."""
    v = np.asarray(vals, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) < 4 or v.std() == 0:
        return 0.0
    m = v - v.mean()
    g1 = float(np.mean(m ** 3) / m.std() ** 3)
    g2 = float(np.mean(m ** 4) / m.std() ** 4)  # non-excess kurtosis
    return (g1 * g1 + 1.0) / g2 if g2 > 0 else 0.0


def make_net(args, lam, seed):
    if args.generator == "polymer_ligation":
        kw = {"alphabet_size": args.alphabet_size,
              "max_length": args.max_length,
              "n_ligations": args.n_reactions // 2,
              "p_cat": args.p_cat_fixed, "p_info": args.p_info,
              "p_selfrep": args.p_selfrep,
              "p_cat_required": args.p_cat_required,
              "k_sat": args.k_sat, "food_len": args.food_len,
              "seed": seed}
        kw[args.lambda_param] = lam
        return polymer_ligation(**kw)
    if args.lambda_param != "p_cat":
        raise ValueError("--lambda-param=p-info is polymer-only")
    return random_crn(
        n_species=args.n_species, n_reactions=args.n_reactions,
        p_cat=lam, n_food=args.n_food, max_reactants=args.max_reactants,
        p_cat_required=args.p_cat_required, k_sat=args.k_sat,
        exclusive_catalysis=args.exclusive_catalysis, seed=seed)


def bistability_test(net, args, x_seed):
    """Same network, organized seed vs food-only start. Returns
    (phi_seeded, phi_unseeded); a split = dynamical bistability."""
    def phi_of(x0):
        _, X = simulate(net, (0.0, args.t_max),
                        feed_rate=args.feed_rate,
                        dilution_rate=args.dilution_rate,
                        x0=np.asarray(x0), steady_check=False)
        finite = X[~np.isnan(X).any(axis=1)]
        if not len(finite):
            return float("nan")
        return organizational_closure(net, net.propensity(finite[-1]))

    food_x0 = np.zeros(net.n_species)
    food_x0[list(net.food)] = args.feed_rate / args.dilution_rate
    return phi_of(np.asarray(x_seed)), phi_of(food_x0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/criticality.jsonl")
    ap.add_argument("--generator", choices=["random_crn", "polymer_ligation"],
                    default="polymer_ligation")
    ap.add_argument("--lambda-param", choices=["p_cat", "p_info"],
                    default="p_info")
    ap.add_argument("--lam-grid", type=float, nargs="+",
                    default=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0])
    ap.add_argument("--n-reps", type=int, default=5)
    ap.add_argument("--n-species", type=int, default=40)
    ap.add_argument("--n-reactions", type=int, default=400)
    ap.add_argument("--n-food", type=int, default=8)
    ap.add_argument("--max-reactants", type=int, default=1)
    ap.add_argument("--alphabet-size", type=int, default=2)
    ap.add_argument("--max-length", type=int, default=10)
    ap.add_argument("--p-cat-fixed", type=float, default=0.6)
    ap.add_argument("--p-info", type=float, default=0.7)
    ap.add_argument("--p-selfrep", type=float, default=0.3)
    ap.add_argument("--p-cat-required", type=float, default=0.5)
    ap.add_argument("--food-len", type=int, default=2)
    ap.add_argument("--k-sat", type=float, default=1.0)
    ap.add_argument("--exclusive-catalysis", action="store_true")
    ap.add_argument("--feed-rate", type=float, default=5.0)
    ap.add_argument("--dilution-rate", type=float, default=1.0)
    ap.add_argument("--t-max", type=float, default=60.0)
    ap.add_argument("--n-bistability", type=int, default=3,
                    help="networks per lambda to bistability-test")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rows = []
    total = len(args.lam_grid) * args.n_reps
    done = 0
    for lam in args.lam_grid:
        for rep in range(args.n_reps):
            seed = args.seed + rep * 1009 + int(lam * 1000) * 13
            net = make_net(args, lam, seed)
            row = {"lam": lam, "lambda_param": args.lambda_param,
                   "seed": seed, "generator": args.generator}
            row.update(measure_run(
                net, args.feed_rate, args.dilution_rate,
                t_max=args.t_max, with_persistence=False))
            x_ss = row.get("x_ss")
            if x_ss is not None and not row["runaway"]:
                row["alpha"] = jacobian_spectral_abscissa(
                    net, np.asarray(x_ss), args.feed_rate,
                    args.dilution_rate)
                row["invasion"] = invasion_pressure(
                    net, np.asarray(x_ss), args.feed_rate,
                    args.dilution_rate)
                if rep < args.n_bistability:
                    ps, pu = bistability_test(net, args, x_ss)
                    row["phi_seeded"], row["phi_unseeded"] = ps, pu
                    row["bistable"] = (
                        np.isfinite(ps) and np.isfinite(pu)
                        and abs(ps - pu) > 0.1)
            rows.append(row)
            done += 1
            print(f"\r{done}/{total} {args.lambda_param}={lam:.2f} "
                  f"phi={row.get('phi', float('nan')):.3f} "
                  f"alpha={row.get('alpha', float('nan')):.3g}",
                  end="", flush=True)
    print()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        for r in rows:
            f.write(json.dumps(
                {k: v for k, v in r.items() if k != "x_ss"}) + "\n")
    print(f"wrote {out}")

    df = pd.DataFrame(rows)
    print(f"\n=== per-lambda summary ({args.lambda_param}) ===")
    print(f"{'lam':>5} {'<phi>':>7} {'var':>7} {'BC':>5} {'<alpha>':>9} "
          f"{'<invas>':>9} {'frac_bist':>9}")
    for lam in args.lam_grid:
        sub = df[df.lam == lam]
        al = sub["alpha"].dropna()
        inv = sub["invasion"].dropna()
        bi = sub["bistable"].dropna()
        print(f"{lam:5.2f} {sub['phi'].mean():7.3f} {sub['phi'].var():7.3f} "
              f"{bimodality_coefficient(sub['phi']):5.2f} "
              f"{al.mean() if len(al) else float('nan'):9.3g} "
              f"{inv.mean() if len(inv) else float('nan'):9.3g} "
              f"{bi.mean() if len(bi) else float('nan'):9.2f}")

    print("\n=== interpretation guide ===")
    print("alpha -> 0 near the boundary + var(BC) peak + bistable split")
    print("= dynamical bifurcation. Flat alpha/var + no bistability")
    print("= ensemble-level structural transition (percolation), not a")
    print("bifurcation of any single network's dynamics.")


if __name__ == "__main__":
    main()
