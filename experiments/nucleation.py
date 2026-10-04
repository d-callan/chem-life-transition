"""Experiment 7 — P(nucleation): mapping the bistable window under SSA.

M5 showed finite copy number restores the nucleation barrier in
exclusive-regime polymer nets (p_cat_required=1): discreteness closes
the cleavage-fragment backdoor in a bounded volume window. This
experiment maps that window: for each network, sweep (volume,
cleavage_rate) and measure

    P(nucleation) = fraction of unseeded SSA trajectories that reach
                    an organized state (phi > phi_threshold)

plus a seeded check per cell (does the organized state even persist at
this V? — distinguishes the nucleation boundary from the extinction
boundary).

The bistable window is where P(nucleation)~0 AND seeded survival ~1.

Usage:
    python experiments/nucleation.py --out data/nucleation.jsonl
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from chemtrans.dynamics.gillespie import simulate as ssa
from chemtrans.dynamics.ode import simulate as ode
from chemtrans.networks.generators import polymer_ligation
from chemtrans.observables.closure import organizational_closure


def phi_at_end(net, Xs, volume):
    conc = Xs[-1] / volume
    return organizational_closure(net, net.propensity(conc))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/nucleation.jsonl")
    ap.add_argument("--seeds", type=int, nargs="+",
                    default=[11118, 9100, 4909, 12127, 6500, 9527])
    ap.add_argument("--volumes", type=float, nargs="+",
                    default=[0.5, 1.0, 2.0, 5.0])
    ap.add_argument("--cleavage-rates", type=float, nargs="+",
                    default=[0.002, 0.02, 0.2])
    ap.add_argument("--n-reps", type=int, default=5)
    ap.add_argument("--alphabet-size", type=int, default=2)
    ap.add_argument("--max-length", type=int, default=10)
    ap.add_argument("--n-ligations", type=int, default=200)
    ap.add_argument("--p-cat", type=float, default=0.6)
    ap.add_argument("--p-info", type=float, default=0.7)
    ap.add_argument("--p-selfrep", type=float, default=0.3)
    ap.add_argument("--food-len", type=int, default=2)
    ap.add_argument("--k-sat", type=float, default=1.0)
    ap.add_argument("--feed-rate", type=float, default=5.0)
    ap.add_argument("--dilution-rate", type=float, default=1.0)
    ap.add_argument("--t-max", type=float, default=150.0)
    ap.add_argument("--phi-threshold", type=float, default=0.05)
    ap.add_argument("--max-events", type=int, default=800_000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rows = []
    total = (len(args.seeds) * len(args.volumes)
             * len(args.cleavage_rates) * (args.n_reps + 1))
    done = 0
    for net_seed in args.seeds:
        for cleav in args.cleavage_rates:
            net = polymer_ligation(
                alphabet_size=args.alphabet_size,
                max_length=args.max_length,
                n_ligations=args.n_ligations,
                p_cat=args.p_cat, p_info=args.p_info,
                p_selfrep=args.p_selfrep, p_cat_required=1.0,
                cleavage_rate=cleav, k_sat=args.k_sat,
                food_len=args.food_len, seed=net_seed)
            _, Xo = ode(net, (0, 200), feed_rate=args.feed_rate,
                        dilution_rate=args.dilution_rate,
                        steady_check=False, n_points=500)
            x_ss = Xo[-1]
            for V in args.volumes:
                # unseeded: nucleation probability
                hits = 0
                for rep in range(args.n_reps):
                    t, Xs = ssa(
                        net, t_max=args.t_max,
                        feed_rate=args.feed_rate,
                        dilution_rate=args.dilution_rate,
                        volume=V, seed=args.seed + rep * 31 + net_seed,
                        max_events=args.max_events)
                    phi = phi_at_end(net, Xs, V)
                    hit = phi > args.phi_threshold
                    hits += hit
                    rows.append({"seed": net_seed, "V": V,
                                 "cleavage_rate": cleav, "rep": rep,
                                 "condition": "unseeded", "phi": phi,
                                 "nucleated": bool(hit)})
                    done += 1
                    print(f"\r{done}/{total} s{net_seed} V={V} "
                          f"clv={cleav} un phi={phi:.3f}",
                          end="", flush=True)
                # seeded control: does the organized state persist here?
                t, Xs = ssa(net, t_max=args.t_max,
                            feed_rate=args.feed_rate,
                            dilution_rate=args.dilution_rate,
                            volume=V, x0=np.asarray(x_ss) * V,
                            seed=args.seed + net_seed,
                            max_events=args.max_events)
                phi_seeded = phi_at_end(net, Xs, V)
                rows.append({"seed": net_seed, "V": V,
                             "cleavage_rate": cleav, "rep": -1,
                             "condition": "seeded", "phi": phi_seeded,
                             "nucleated": phi_seeded > args.phi_threshold})
                done += 1
                print(f"\r{done}/{total} s{net_seed} V={V} "
                      f"clv={cleav} seeded phi={phi_seeded:.3f}",
                      end="", flush=True)
    print()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    print(f"wrote {out}")

    df = pd.DataFrame(rows)
    un = df[df.condition == "unseeded"]
    se = df[df.condition == "seeded"]
    print("\n=== P(nucleation) unseeded | seeded phi — per net ===")
    for seed in args.seeds:
        print(f"\nnet seed={seed}   rows=V, cols=cleavage_rate "
              f"{' '.join(f'{c:>12g}' for c in args.cleavage_rates)}")
        for V in args.volumes:
            cells = []
            for c in args.cleavage_rates:
                pn = un[(un.seed == seed) & (un.V == V)
                        & (un.cleavage_rate == c)]["nucleated"].mean()
                ps = se[(se.seed == seed) & (se.V == V)
                        & (se.cleavage_rate == c)]["phi"].iloc[0]
                cells.append(f"{pn:4.2f}/{ps:4.2f}")
            print(f"  V={V:<5} " + "  ".join(f"{c:>12s}" for c in cells))
    print("\nformat: P(nucleated)/phi_seeded. Bistable window = cells "
          "with ~0 P(nucleation) and high phi_seeded.")


if __name__ == "__main__":
    main()
