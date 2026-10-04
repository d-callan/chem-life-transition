"""Experiment 2 — hysteresis and bistability.

p_cat is a STRUCTURAL parameter — it can't be ramped adiabatically on a
fixed network. So two complementary probes of whether the transition is
first-order-like (discontinuous, needs nucleation) vs continuous:

1. Seeded vs unseeded, same network, same lambda:
   - unseeded: food-only initial condition — can organization nucleate?
   - seeded: food + RAF-core species boosted — does organization persist?
   If seeded maintains Phi>0 where unseeded gives Phi=0, the organized
   state is a genuine attractor with a nucleation barrier (bistability).

2. Feed-rate ramp with continuation on a fixed network:
   up-ramp starts from the food-only state at low feed; down-ramp starts
   from the seeded/organized state at high feed; each step continues from
   the previous steady state. The loop area int|Phi_up - Phi_down| d(feed)
   is the hysteresis statistic.

Usage:
    python experiments/hysteresis.py --out data/hyst.jsonl --plot data/hyst.png
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from chemtrans.measure import measure_run
from chemtrans.networks.generators import random_crn
from chemtrans.observables.closure import raf_core


def seeded_x0(net, feed_rate, dilution_rate, seed_conc):
    """Food at chemostat level + RAF-core species boosted to seed_conc."""
    x0 = np.zeros(net.n_species)
    x0[list(net.food)] = feed_rate / max(dilution_rate, 1e-12)
    core = raf_core(net)
    for i in core:
        x0[i] = max(x0[i], seed_conc)
    return x0, len(core)


def ramp(net, feeds, feed0_state, dilution_rate, t_max):
    """Adiabatic continuation: each feed starts from previous steady state."""
    phis, x = [], feed0_state
    for f in feeds:
        row = measure_run(net, f, dilution_rate, t_max=t_max, x0=x)
        phis.append(row.get("phi"))
        if "x_ss" in row:
            x = np.array(row["x_ss"])
    return phis


def run_one(p_cat, args, seed):
    net = random_crn(
        n_species=args.n_species, n_reactions=args.n_reactions, p_cat=p_cat,
        n_food=args.n_food, max_reactants=args.max_reactants,
        p_cat_required=args.p_cat_required, seed=seed,
    )
    row = {"p_cat": p_cat, "seed": seed,
           "generator_params": vars(args)}

    unseeded = measure_run(net, args.feed_rate, args.dilution_rate,
                           t_max=args.t_max)
    x0_seed, n_core = seeded_x0(net, args.feed_rate, args.dilution_rate,
                              args.seed_conc)
    seeded = measure_run(net, args.feed_rate, args.dilution_rate,
                         t_max=args.t_max, x0=x0_seed)

    row.update({
        "raf_size": n_core,
        "phi_unseeded": unseeded.get("phi"),
        "phi_seeded": seeded.get("phi"),
        "bistable": (
            (unseeded.get("phi") or 0) < args.phi_threshold
            and (seeded.get("phi") or 0) >= args.phi_threshold
        ),
    })

    # feed-rate ramp, only worth it if the network CAN be organized
    if n_core > 0:
        feeds = np.linspace(args.feed_lo, args.feed_hi, args.ramp_steps)
        x_food = np.zeros(net.n_species)
        x_food[list(net.food)] = args.feed_lo / args.dilution_rate
        phi_up = ramp(net, feeds, x_food, args.dilution_rate, args.ramp_t)
        x_org = np.array(seeded.get("x_ss", x0_seed))
        phi_down = ramp(net, feeds[::-1], x_org, args.dilution_rate,
                        args.ramp_t)
        area = float(np.trapezoid(
            np.abs(np.array(phis_or_zero(phi_up))
                   - np.array(phis_or_zero(phi_down))),
            x=feeds,
        ))
        row["ramp"] = {
            "feeds": feeds.tolist(),
            "phi_up": phi_up,
            "phi_down": phi_down,
            "loop_area": area,
        }
    return row


def phis_or_zero(phis):
    return [p if p is not None else 0.0 for p in phis]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/hysteresis.jsonl")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--n-species", type=int, default=40)
    ap.add_argument("--n-reactions", type=int, default=200)
    ap.add_argument("--n-food", type=int, default=8)
    ap.add_argument("--max-reactants", type=int, default=1)
    ap.add_argument("--p-cat", type=float, nargs="+",
                    default=[0.15, 0.2, 0.25, 0.3])
    ap.add_argument("--p-cat-required", type=float, default=0.5)
    ap.add_argument("--n-reps", type=int, default=10)
    ap.add_argument("--feed-rate", type=float, default=1.0)
    ap.add_argument("--dilution-rate", type=float, default=0.1)
    ap.add_argument("--seed-conc", type=float, default=1.0)
    ap.add_argument("--phi-threshold", type=float, default=0.1)
    ap.add_argument("--feed-lo", type=float, default=0.2)
    ap.add_argument("--feed-hi", type=float, default=2.0)
    ap.add_argument("--ramp-steps", type=int, default=8)
    ap.add_argument("--ramp-t", type=float, default=100.0,
                    help="integration time per ramp step")
    ap.add_argument("--t-max", type=float, default=150.0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w") as f:
        for p_cat in args.p_cat:
            for rep in range(args.n_reps):
                seed = args.seed + rep * 1009 + hash(round(p_cat, 6)) % 997
                row = run_one(p_cat, args, seed)
                row["run_id"] = f"s{seed}_p{p_cat:.3f}"
                f.write(json.dumps(row) + "\n")
                print(f"p={p_cat:.2f} rep={rep} "
                      f"phi_uns={row['phi_unseeded']!r:>7} "
                      f"phi_seed={row['phi_seeded']!r:>7} "
                      f"bistable={row['bistable']} "
                      f"loop={row.get('ramp',{}).get('loop_area')}", flush=True)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
