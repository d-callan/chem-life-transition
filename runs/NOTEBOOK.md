# Experiment notebook

Lab notebook for simulation runs. Each entry records the exact command,
the result, and a tentative interpretation. Data files live in `data/`
(gitignored — local artifacts; regenerate by rerunning the command).
Row-level detail is in the jsonl; this file is for the *meaning*.

Conventions: runs listed newest-last. `PY` = chemtrans conda env python.

---

## 2026-10-04 — first sweep attempt (failed, instructive)

**Cmd:** `python experiments/sweep_lambda.py --out data/sweep_test.jsonl --n-reps 8 --t-max 300`

**Result:** hung — first run pegged CPU for 7+ min with zero output rows.

**Interpretation:** catalytic mass-action with `cat_boost=20` has no resource
saturation; supercritical networks diverge and LSODA crawls at tiny steps.
Fixed by adding a terminal blowup event (`x_cap`, default 1e5) in
`ode.simulate`; runaway is a *recorded outcome*, not an integration error.

## 2026-10-04 — sparse-regime probes

**Cmd:** `python experiments/sweep_lambda.py --out data/probe.jsonl --p-cat 0.1 0.3 0.5 --n-reps 3 --t-max 100` (defaults: N=60, R=120, max_reactants=2)
and full sweep `data/sweep_m1.jsonl` (110 runs, p_cat 0.0→0.5, +`--plot data/phi_m1.png`)

**Result:** Φ ≈ 0 almost everywhere; only 3/110 runs nonzero (max 0.129,
one runaway-ish 1.0). `rho_K` = 15–50 whenever any catalysis exists.

**Interpretation:** strict RAF cores are percolation-limited here — a core
reaction needs ALL reactants producible, and 2-reactant random reactions
rarely close. Also: `rho_K >> 1` even in dead networks means `cat_boost=20`
makes ρ(K)=1 a meaningless literal threshold; treat rho_K as a monotone
proxy only (calibration TODO).

## 2026-10-04 — single-reactant probe

**Cmd:** `--out data/probe_r1.jsonl --p-cat 0.05 0.15 0.25 0.35 0.5 --n-reps 6 --max-reactants 1 --t-max 150`

**Result:** still mostly Φ=0; occasional nonzero (0.063 @0.25, 0.840 @0.5).

**Interpretation:** food-fed chains help but closure still rare at N=60,
R=120. RAF needs directed catalytic cycles — density is the binding
constraint, not reactant arity alone.

## 2026-10-04 — dense regime: organized phase exists

**Cmd:** `--out data/probe_dense.jsonl --n-species 40 --n-reactions 200 --n-food 8 --max-reactants 1 --p-cat 0.5 0.7 0.9 --n-reps 6`

**Result:** Φ = 0.93–0.99 on EVERY run at all three p_cat values.

**Interpretation:** the organized regime is real and dominant at density
5 rxn/species. Transition must live at lower p_cat in this regime.

## 2026-10-04 — THE transition curve (milestone-1 headline)

**Cmd:** `--out data/sweep_dense.jsonl --plot data/phi_dense.png --n-species 40 --n-reactions 200 --n-food 8 --max-reactants 1 --p-cat 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.4 0.5 --n-reps 8 --t-max 150`

**Result:**
| p_cat | outcome |
|---|---|
| ≤0.10 | Φ=0, all 8 reps |
| 0.15 | 2/8 nonzero (0.046, 0.463) |
| 0.20 | bimodal: {0,0.489,0.051,0.066,0.003,0,0.002,0} |
| 0.25 | all nonzero, 0.06–0.92 |
| ≥0.30 | reliably 0.9–0.97 |

**Interpretation:** transition region p_cat ≈ 0.2–0.3 with strong bimodality —
ensemble members either contain a RAF motif or don't. Looks like a
threshold, but at this point we can't tell ensemble-percolation from a
dynamical transition (M2 answers this).

## 2026-10-04 — hysteresis probe (M2a)

**Cmd:** `python experiments/hysteresis.py --out data/hyst_probe.jsonl --p-cat 0.2 0.25 --n-reps 4` (dense regime defaults)

**Result:** `phi_seeded ≈ phi_unseeded` on every network (e.g. 0.489 vs
0.492, 0.051 vs 0.051); all feed-rate ramp loop areas ≈ 0. No bistability.

**Interpretation:** organization nucleates from food alone when the RAF
motif exists — `p_cat_required=0.5` means half the catalyzed reactions can
fire spontaneously, which seeds the catalysts. The Φ(p_cat) transition is
**structural percolation over the ensemble**, not dynamical bistability.

## 2026-10-04 — hysteresis with fully-required catalysis (M2b)

**Cmd:** `python experiments/hysteresis.py --out data/hyst_req1.jsonl --p-cat 0.3 0.4 --n-reps 4 --p-cat-required 1.0`

**Result:** STILL no bistability — unseeded runs reach Φ≈0.97; ramp loops
~0. `bistable=False` everywhere.

**Interpretation:** deeper cause: a catalyst *species* can be produced by a
spontaneous non-core reaction as well as its catalyzed core reaction —
topology requires "catalyst producible by core" but dynamics only needs
"catalyst producible". RAF existence ≠ nucleation barrier. A real barrier
needs a generator option where complex species have ONLY catalyzed
production routes (`exclusive_catalysis`, planned M3).

## 2026-10-04 — finite-size scaling (M2c)

**Cmd:** `python experiments/finite_size.py --out data/fs_probe.jsonl --plot data/fs_probe.png --sizes 20 40 80 --p-cat 0.1 0.2 0.25 0.3 0.4 --n-reps 8` (density held at 5 rxn/species)

**Result:** P(Φ>0.1) rises through p_cat ≈ 0.2–0.3 at all three N; maybe
slightly delayed at N=80. No clear sharpening visible at 8 reps.

**Interpretation:** transition at fixed mean degree regardless of N is
textbook percolation behavior — consistent with the hysteresis finding.
Sharpening needs more reps to resolve (bimodal variance dominates).

---

## Open threads going into M3

1. `exclusive_catalysis` knob: species producible ONLY via catalyzed
   reactions → real nucleation barriers → hysteresis test becomes live.
2. `P(organized)` is arguably the better order parameter than mean Φ for a
   percolation transition (probability vs magnitude) — critical_points.py
   should support both.
3. `rho_K` needs a dimensionless normalization (vs dilution rate, typical
   concentration scale) before it can predict λc.
4. Width/sharpening question is statistically unresolved — needs ≥20 reps.

---

## 2026-10-04 — exclusive_catalysis v1 fails to create barrier

**Cmd:** `python experiments/hysteresis.py --out data/hyst_excl.jsonl --p-cat 0.25 0.3 --n-reps 6 --p-cat-required 1.0 --exclusive-catalysis`

**Result:** still `bistable=False`, phi_unseeded up to 0.89.

**Interpretation:** v1 only removed spontaneous routes to exclusive
*products*; catalysts of required reactions were still random species
reachable spontaneously — one spontaneous catalyst ignites the whole core.
RAF topology ≠ nucleation barrier; the barrier lives in the *catalysts'*
production routes. Fixed by rewiring required-reaction catalysts into the
exclusive set.

## 2026-10-04 — BISTABILITY ACHIEVED (exclusive v2)

**Cmd:** `python experiments/hysteresis.py --out data/hyst_excl2.jsonl --p-cat 0.25 0.3 --n-reps 6 --p-cat-required 1.0 --exclusive-catalysis`

**Result:** 12/12 runs `bistable=True`. phi_unseeded = 0 everywhere;
phi_seeded = 0.13–0.96. Feed-rate ramp loop areas 0.2–1.7 (vs ~1e-4
without exclusivity).

**Interpretation:** with catalysts confined to the exclusive set, the
organized state is a true attractor requiring nucleation — the model now
exhibits BOTH regimes: ensemble percolation (default) and dynamical
first-order-like transition (exclusive). Two different kinds of "origin of
life" transitions, one knob apart.

## 2026-10-04 — runaway audit: everything was measured mid-explosion

**Discovery:** all 30 rows of `data/discrim_probe.jsonl` (and, on
inspection, effectively all earlier organized runs) were `runaway: true` —
the blowup event fires because mass-action autocatalysis is unbounded.
Φ was being measured at the blowup state, not a steady state.

**Two divergence mechanisms found:**
1. `A -> A+B` spontaneous replicators (reactant∈products) — removed:
   `random_crn` now forbids reactant∩product overlap; amplification lives
   only in the catalysis channel.
2. Directed cycles with net species gain (A->B+C, B->C+A...) — fixed by
   `k_sat`: Michaelis-Menten saturation applied to the WHOLE propensity
   (substrates + catalyst), bounding every rate → guaranteed steady states.

**Cmd:** debug — `random_crn(40,200,p_cat=0.4,...,k_sat=1.0,seed=5)` +
`simulate(t_span=(0,300))`

**Result:** clean steady state, max conc ~689, no warnings.

## 2026-10-04 — first discrimination run (M3 headline)

**Cmd:** `python experiments/transition_discrimination.py --out data/discrim_ksat.jsonl --p-cat 0.1 0.2 0.25 0.3 0.4 --n-reps 6 --n-bootstrap 50 --t-max 120 --k-sat 1.0`

**Result:** transition persists under saturating kinetics at the same
location (Φ: 0@0.1, bimodal@0.2, ~0.6-0.88@0.4). λc estimates:
phi 0.30 [0.2,0.3], A 0.30 [0.1,0.3], S 0.20 [0.1,0.2], R0 0.25 [0.2,0.4],
E 0.30 [0.2,0.3], tail_mass 0.20, n_core 0.30. All 21 pairwise alignments
≥0.55. **Verdict: Model B** (single coupled transition). Model C check
skipped — too few organized-phase runs in 30.

**Interpretation (cautious):** co-location is partly *by construction* —
S/R0/E are all functions of the same RAF core, so they inherit its
threshold. "No evidence of staggered thresholds in random_crn" is the
honest claim. The observables need decoupling (or a different generator)
before Model B is a strong claim. H is still unmeasurable (no information
carrier) — polymer_ligation territory.

## Open threads going into M4

1. The discriminating power of the A/B/C test currently rests on whether
   secondary observables can vary independently of the RAF core. In
   random_crn they largely can't — need mechanisms where e.g. assembly
   tail or evolvability decouple (polymer_ligation should do this).
2. `x_cap` semantics: for k_sat models, runaway should be near-impossible;
   keep cap high (1e7) as a sanity alarm only.
3. Grid resolution: 5-7 lambda points gives coarse λc; discrimination
   conclusions need ≥10 points and ≥12 reps.
4. Exclusive-catalysis regime changes the whole experiment: unseeded never
   organizes → for that regime the discrimination sweep should seed runs
   (or measure P(nucleation) instead of P(organized)).
