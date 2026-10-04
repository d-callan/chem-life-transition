"""Experiment 5 — 2-D critical-manifold test (polymer_ligation).

M4 showed the A/B/C verdict depends on the control axis: sweeping p_info
gives staggered lambda_c's (Model-A-like), sweeping p_cat gives coupled
behavior (Model-B-like). If organization lives on a 2-D critical manifold
in (p_cat, p_info) space, the right question is not "are the 1-D lambda_c
estimates aligned" but "do all order parameters switch on the SAME
boundary curve in the plane".

Method:
  - sweep a (p_cat x p_info) grid, n_reps networks per cell
  - binarize each observable per run (organized vs not, natural threshold)
  - fit logistic boundary P(on) = sigmoid(b0 + b1*p_cat + b2*p_info) per
    observable, with bootstrap CIs on the boundary parameters
  - compare boundaries: co-located -> one critical manifold (Model B);
    separated -> genuinely distinct thresholds (Model A)

Usage:
    python experiments/manifold_2d.py --out data/manifold_2d.jsonl
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from chemtrans.measure import measure_run
from chemtrans.networks.generators import polymer_ligation

# observable column -> threshold for "this property is ON"
THRESHOLDS = {
    "phi": 0.1,
    "A": 0.5,
    "S": 0.5,
    "R0": 1.0,
    "E": 0.3,
    "H": 0.5,
    "n_core_species": 0.5,   # >0 means a RAF core exists
    "assembly.tail_mass": 0.1,
}


def get_obs(row: dict, name: str):
    if name.startswith("assembly."):
        return row.get("assembly", {}).get(name.split(".")[1])
    return row.get(name)


def logistic_fit(df: pd.DataFrame, obs: str, thresh: float,
                 n_boot: int = 0, seed: int = 0):
    """Fit P(obs>thresh) = sigmoid(b0 + b1 p_cat + b2 p_info).

    Returns (params, boot_array): params = (b0,b1,b2) MLE on standardized
    axes; boot_array shape (n_boot,3) for CI. Boundary line in raw space:
    b0 + b1*zc + b2*zi = 0 -> solve for p_info given p_cat.
    """
    y = np.array([1.0 if (v := get_obs(r, obs)) is not None
                  and not pd.isna(v) and v > thresh else 0.0
                  for r in df.to_dict("records")])
    X = df[["p_cat", "p_info"]].to_numpy()
    mu, sd = X.mean(0), X.std(0) + 1e-12
    Z = np.column_stack([np.ones(len(X)), (X - mu) / sd])

    def nll(beta):
        z = np.clip(Z @ beta, -30, 30)
        return float(np.sum(np.logaddexp(0, z) - y * z))

    res = minimize(nll, np.zeros(3), method="Nelder-Mead")
    beta = res.x
    boots = []
    rng = np.random.default_rng(seed)
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        yb, Zb = y[idx], Z[idx]
        if yb.std() == 0:
            continue

        def nll_b(beta):
            z = np.clip(Zb @ beta, -30, 30)
            return float(np.sum(np.logaddexp(0, z) - yb * z))

        boots.append(minimize(nll_b, beta, method="Nelder-Mead").x)
    return beta, np.array(boots), (mu, sd)


def boundary_point(beta, mu, sd, p_cat):
    """p_info value where P=0.5 at given p_cat, or None if out of [0,1]."""
    zc = (p_cat - mu[0]) / sd[0]
    zi = -(beta[0] + beta[1] * zc) / max(beta[2], 1e-12)
    p_info = zi * sd[1] + mu[1]
    return p_info if 0.0 <= p_info <= 1.0 else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/manifold_2d.jsonl")
    ap.add_argument("--p-cat-grid", type=float, nargs="+",
                    default=[0.2, 0.4, 0.6, 0.8])
    ap.add_argument("--p-info-grid", type=float, nargs="+",
                    default=[0.0, 0.25, 0.5, 0.75, 1.0])
    ap.add_argument("--n-reps", type=int, default=4)
    ap.add_argument("--alphabet-size", type=int, default=2)
    ap.add_argument("--max-length", type=int, default=10)
    ap.add_argument("--n-ligations", type=int, default=200)
    ap.add_argument("--p-selfrep", type=float, default=0.3)
    ap.add_argument("--p-cat-required", type=float, default=0.5)
    ap.add_argument("--food-len", type=int, default=2)
    ap.add_argument("--k-sat", type=float, default=1.0)
    ap.add_argument("--feed-rate", type=float, default=5.0)
    ap.add_argument("--dilution-rate", type=float, default=1.0)
    ap.add_argument("--t-max", type=float, default=60.0)
    ap.add_argument("--n-bootstrap", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rows = []
    total = len(args.p_cat_grid) * len(args.p_info_grid) * args.n_reps
    done = 0
    for pc in args.p_cat_grid:
        for pi in args.p_info_grid:
            for rep in range(args.n_reps):
                seed = (args.seed + rep * 1009
                        + int(pc * 1000) * 7 + int(pi * 1000) * 13)
                net = polymer_ligation(
                    alphabet_size=args.alphabet_size,
                    max_length=args.max_length,
                    n_ligations=args.n_ligations,
                    p_cat=pc, p_info=pi, p_selfrep=args.p_selfrep,
                    p_cat_required=args.p_cat_required,
                    k_sat=args.k_sat, food_len=args.food_len, seed=seed)
                row = {"p_cat": pc, "p_info": pi, "seed": seed,
                       "run_id": f"mani_s{seed}_c{pc:.2f}_i{pi:.2f}"}
                row.update(measure_run(
                    net, args.feed_rate, args.dilution_rate,
                    t_max=args.t_max, with_persistence=True,
                    persistence_seed=seed))
                rows.append(row)
                done += 1
                print(f"\r{done}/{total} p_cat={pc} p_info={pi} "
                      f"phi={row.get('phi', float('nan')):.3f}",
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

    # raw hit-rate tables first — the eyeball test matters
    for obs, thresh in THRESHOLDS.items():
        if get_obs(rows[0], obs) is None:
            continue
        tab = np.full((len(args.p_cat_grid), len(args.p_info_grid)), np.nan)
        for i, pc in enumerate(args.p_cat_grid):
            for j, pi in enumerate(args.p_info_grid):
                sub = df[(df.p_cat == pc) & (df.p_info == pi)]
                vals = [get_obs(r, obs) for r in sub.to_dict("records")]
                vals = [v for v in vals if v is not None and not pd.isna(v)]
                tab[i, j] = np.mean([v > thresh for v in vals]) \
                    if vals else np.nan
        print(f"\n=== P({obs} > {thresh}) rows=p_cat, cols=p_info ===")
        print("        " + "  ".join(f"{v:5.2f}" for v in args.p_info_grid))
        for i, pc in enumerate(args.p_cat_grid):
            print(f"{pc:5.2f}  " + "  ".join(f"{v:5.2f}" for v in tab[i]))

    # logistic boundary fits
    print("\n=== boundary params beta (standardized) ===")
    print(f"{'obs':22s} {'b0':>7s} {'b_cat':>7s} {'b_info':>7s}  "
          f"boundary p_info(p_cat) [boot CI]")
    rng = np.random.default_rng(args.seed)
    for obs, thresh in THRESHOLDS.items():
        if get_obs(rows[0], obs) is None:
            continue
        beta, boots, (mu, sd) = logistic_fit(
            df, obs, thresh, n_boot=args.n_bootstrap,
            seed=int(rng.integers(1e9)))
        mids = []
        for pc in args.p_cat_grid:
            mid = boundary_point(beta, mu, sd, pc)
            ci = ""
            if mid is not None and len(boots):
                bpts = [boundary_point(b, mu, sd, pc) for b in boots]
                bpts = [b for b in bpts if b is not None]
                if bpts:
                    lo, hi = np.percentile(bpts, [5, 95])
                    ci = f"[{lo:.2f},{hi:.2f}]"
            mids.append(f"{mid:.2f}{ci}" if mid is not None else "  -  ")
        print(f"{obs:22s} {beta[0]:7.2f} {beta[1]:7.2f} {beta[2]:7.2f}  "
              + " ".join(f"{m:>12s}" for m in mids))


if __name__ == "__main__":
    main()
