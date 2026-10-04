"""Experiment 8 — mean first-passage times: the barrier as a rate.

The nucleation sweep showed the cleavage-door 'barrier' is not binary:
discreteness turns it into a Poisson process gamma(V, cleavage, net).
Bistability is a timescale statement:

    window = { 1/gamma_nuc  >>  t_obs  AND  1/gamma_ext  >>  t_obs }

This experiment measures both rates directly on the phi(t) trajectory:
  - unseeded SSA: tau_nuc = first t with phi > threshold (censored at
    t_max if never)
  - seeded SSA:   tau_ext = first t after which phi stays below
    threshold (censored if still organized at t_max)

Exponential-MLE rate estimate under right-censoring:
    gamma = n_events / sum(min(tau_i, t_max))

Usage:
    python experiments/first_passage.py --out data/first_passage.jsonl
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
from chemtrans.observables.closure import organizational_closure, raf_reactions


def phi_series(net, raf_idx, Xs, volume):
    """phi(t) along an SSA trajectory (counts -> concentrations)."""
    total_flux = []
    core_flux = []
    for row in Xs:
        a = net.propensity(row / volume)
        total_flux.append(a.sum())
        core_flux.append(a[raf_idx].sum() if len(raf_idx) else 0.0)
    total_flux = np.asarray(total_flux)
    core_flux = np.asarray(core_flux)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(total_flux > 0, core_flux / total_flux, 0.0)


def first_above(ts, phi, thresh):
    idx = np.argmax(phi > thresh)
    return float(ts[idx]) if phi[idx] > thresh else np.nan


def last_above_then_below(ts, phi, thresh):
    """First t where phi drops below thresh and never recovers."""
    above = phi > thresh
    if not above.any():
        return 0.0
    last = np.where(above)[0][-1]
    if last == len(phi) - 1:
        return np.nan  # still organized at t_max (censored)
    return float(ts[last + 1])


def gamma_mle(times, t_max):
    """Exponential rate MLE with right-censoring: events/exposure."""
    t = np.asarray(times, dtype=float)
    events = np.isfinite(t)
    exposure = np.where(events, t, t_max).sum()
    return float(events.sum() / exposure) if exposure > 0 else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/first_passage.jsonl")
    ap.add_argument("--seeds", type=int, nargs="+", default=[4909, 9100, 11118])
    ap.add_argument("--volumes", type=float, nargs="+", default=[0.5, 1.0, 2.0, 5.0])
    ap.add_argument("--cleavage-rates", type=float, nargs="+", default=[0.02, 0.2])
    ap.add_argument("--n-reps", type=int, default=4)
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
    ap.add_argument("--t-max", type=float, default=120.0)
    ap.add_argument("--phi-threshold", type=float, default=0.05)
    ap.add_argument("--max-events", type=int, default=800_000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rows = []
    total = (len(args.seeds) * len(args.volumes)
             * len(args.cleavage_rates) * 2 * args.n_reps)
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
            raf_idx = list(raf_reactions(net))
            _, Xo = ode(net, (0, 200), feed_rate=args.feed_rate,
                        dilution_rate=args.dilution_rate,
                        steady_check=False, n_points=500)
            x_ss = Xo[-1]
            for V in args.volumes:
                for cond in ["unseeded", "seeded"]:
                    taus = []
                    for rep in range(args.n_reps):
                        x0 = (np.asarray(x_ss) * V if cond == "seeded"
                              else None)
                        t, Xs = ssa(net, t_max=args.t_max,
                                    feed_rate=args.feed_rate,
                                    dilution_rate=args.dilution_rate,
                                    volume=V, x0=x0,
                                    seed=args.seed + rep * 31 + net_seed,
                                    max_events=args.max_events,
                                    record_dt=0.5)
                        phi = phi_series(net, raf_idx, Xs, V)
                        tau = (first_above(t, phi, args.phi_threshold)
                               if cond == "unseeded"
                               else last_above_then_below(
                                   t, phi, args.phi_threshold))
                        taus.append(tau)
                        done += 1
                        print(f"\r{done}/{total} s{net_seed} V={V} "
                              f"clv={cleav} {cond[:4]} tau={tau:.1f}",
                              end="", flush=True)
                    g = gamma_mle(taus, args.t_max)
                    rows.append({"seed": net_seed, "V": V,
                                 "cleavage_rate": cleav,
                                 "condition": cond,
                                 "taus": [None if np.isnan(x) else x
                                          for x in taus],
                                 "n_events": int(np.isfinite(taus).sum()),
                                 "gamma": g,
                                 "mean_lifetime": 1.0 / g if g > 0
                                 else float("inf")})
    print()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    print(f"wrote {out}")

    df = pd.DataFrame(rows)
    print("\n=== mean lifetime 1/gamma (time units; inf = no event in t_max) ===")
    for seed in args.seeds:
        print(f"\nnet seed={seed}   rows=V, cols=cleavage_rate   "
              f"{' '.join(f'{c:>14g}' for c in args.cleavage_rates)}")
        for V in args.volumes:
            cells = []
            for c in args.cleavage_rates:
                tn = df[(df.seed == seed) & (df.V == V)
                        & (df.cleavage_rate == c)
                        & (df.condition == "unseeded")]["mean_lifetime"].iloc[0]
                te = df[(df.seed == seed) & (df.V == V)
                        & (df.cleavage_rate == c)
                        & (df.condition == "seeded")]["mean_lifetime"].iloc[0]
                cells.append(f"{tn:6.0f}/{te:6.0f}")
            print(f"  V={V:<5} " + "  ".join(f"{c:>14s}" for c in cells))
    print("\nformat: tau_nuc/tau_ext. Bistable window = both >> t_obs "
          f"({args.t_max}).")


if __name__ == "__main__":
    main()
