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

# ============================ MILESTONE 4 ============================
# Polymer/string chemistry: species as sequences, ligation/cleavage,
# template catalysis, true assembly index, heredity observable H.

## 2026-10-04 — polymer_ligation redesign (three iterations)

**Problem found:** first polymer implementation produced essentially no
RAF cores (n_core=0 in ~20 networks across densities), so all dynamics
probes were dead (phi=0, S=0) and the first discrimination probe's
"Model A" verdict was sparse-data noise. Root causes identified and fixed:

1. **One ligation route per product.** Each string had a single recorded
   (a,b) split; a self-rep catalyst closes a RAF only if THAT split's
   reactants are food/core-produced (prob ~ p_cat^depth — never).
   FIX: enumerate every split of every product (a+b->ab forall splits),
   catalyst assigned per PRODUCT (template matching is a property of
   what's assembled). Cleavage = all-splits hydrolysis. This is more
   honest chemistry anyway — a string can form by any split.

2. **No self-replicators.** Template catalysts were chosen among
   SUPERSTRINGS of the product only (longer strings containing the
   motif) — RAF impossible (catalyst needs even longer precursors).
   FIX: `p_selfrep` knob — fraction of informational catalysts that are
   exact self-templates (cat == product, the minimal replicator RAF);
   remainder are kin (super- OR substring — bidirectional sequence
   sharing, both count as informational for H).

3. **Pool too sparse at short lengths.** Uniform pair-sampling over a
   growing pool makes monomer+monomer ligations rare.
   FIX: length-biased ligation sampling (~1/len) + `food_len` knob
   (strings with len <= food_len are food; food_len=2 = "rich broth").

4. **p_info knob added** (the actual deliverable): p_cat = catalytic
   connectivity, p_info = fraction of catalysts that are
   template-matched vs random. H≈1 by construction before; now H varies
   0->1 with p_info and can decouple from phi.

**Result after fix (probe):**
`polymer_ligation(alphabet=2, max_len=10, n_lig=200, p_cat=0.5, k_sat=1.0)`
RAF cores: food_len=1 -> 2/8 networks (sizes 0,1,0,0,1,0,0,0);
           food_len=2 -> 8/8 networks (sizes 1,1,7,32,14,2,7,29).
Networks ~500 species / ~6300 reactions.

## 2026-10-04 — propensity vectorization + sparse matrices

Polymer nets (6300 rxns) made the Python-loop propensity the bottleneck:
52.9s per t_max=60 sim. Fixes: vectorized fast path in
ReactionNetwork.propensity (padded reactant-index array + catalyst CSR
matrix; falls back to loop for >2 reactants / nu!=1 / mixed k_sat —
verified bit-identical to slow path), stoichiometry S as csr in ode.simulate.
Result: ~0.4-1s per LSODA sim (~50-100x). Also fixed: evolvability_proxy
Reaction copies were dropping k_sat.

**Probe:** `measure_run` over p_info x p_cat (food_len=2):
```
p_info=0.0: phi=0.000 H~0.01  S=0     n_core=0   (both p_cat)
p_info=0.3: phi<=0.001 H~0.1-0.9 S~0-0.999 n_core 0-1
p_info=0.7: phi 0.23-0.33 H 0.72-0.81 S=0.999 R0~1.2-1.3 n_core 12-29
p_info=1.0: phi 0.05-0.36 H 0.09-0.75 S=0.999 R0~1.5-1.7 n_core 25
```
~5s per run incl. persistence. NOTE: random catalysis alone (p_info=0)
never organizes — template catalysis does the organizational work here.

## 2026-10-04 — polymer discrimination probes (the M4 headline)

**Cmd A:** `python experiments/transition_discrimination.py --generator polymer_ligation --lambda-param p_info --p-cat 0.0 0.1 0.2 0.3 0.5 0.7 1.0 --p-cat-fixed 0.5 --alphabet-size 2 --max-length 10 --n-reactions 400 --food-len 2 --k-sat 1.0 --feed-rate 5.0 --dilution-rate 1.0 --t-max 60 --n-reps 6 --n-bootstrap 50 --out data/discrim_polymer_pinfo.jsonl`

**Result A (sweeping p_info):** phi: all 0 @ p_info<=0.1; bimodal
@0.2-0.3 (0 / 0.1-0.33); organized @>=0.5 (0.03-0.56). lc estimates:
S 0.2, H 0.2 << phi 0.5, A 0.5, R0 0.5 << E 0.7, tail_mass 1.0, n_core 1.0.
frac_pairs_aligned = 0.286 -> VERDICT: Model A (staggered cascade).
Model C check also found internal critical pts in organized subset
(S 0.2; R0/E/H ~0.7; n_core 1.0) — consistent with nested structure.

**Cmd B:** `python experiments/transition_discrimination.py --generator polymer_ligation --lambda-param p_cat --p-cat 0.2 0.4 0.6 0.8 --p-info 0.7 --alphabet-size 2 --max-length 10 --n-reactions 400 --food-len 2 --k-sat 1.0 --feed-rate 5.0 --dilution-rate 1.0 --t-max 60 --n-reps 4 --n-bootstrap 30 --out data/discrim_polymer_pcat.jsonl`

**Result B (sweeping p_cat):** phi rises gradually 0.01->~0.3 over the
whole grid; ALL observables' lc at grid edge 0.8 except S=0.2.
frac_aligned = 0.75 -> VERDICT: Model B — but lc AT the boundary means
the transition may not have been crossed; inconclusive.

**Interpretation:** the striking thing is the AXIS DEPENDENCE.
- Sweeping p_info: apparent staged ordering — lone self-replicators
  (S,H early ~0.2) precede collective closure (phi,A,R0 ~0.5), which
  precedes complexity/evolvability (~0.7-1.0). Physically sensible:
  a single self-template replicator maintains itself and copies
  sequence info before collectively-autocatalytic cores exist.
- Sweeping p_cat: everything rides together (when it moves at all).
- S precedes closure in BOTH sweeps — S=0.999 whenever ANY self-rep
  exists, so its threshold is structurally lower, not noise.
- Honest caveat: if organization lives on a 2-D critical manifold in
  (p_cat, p_info), different slices produce different apparent
  orderings — "Model A" along one axis could be "Model B" on the
  manifold. True Model A requires staggered thresholds along ALL
  control directions. Next: fine p_info grid at several p_cat values
  to map the manifold, and seed the exclusive/nucleation regimes.

## Open threads going into M5

1. 2-D sweep (p_cat x p_info) to test the critical-manifold hypothesis —
   do the staggered lc's collapse onto one curve in 2-parameter space?
2. S is structurally "easy" (any self-rep suffices) — consider a
   collective-S variant (shock survival of the FULL organized state,
   not just one motif) to avoid it always preceding phi.
3. Polymer food_len=2 regime: which species actually carry the flux?
   Check whether organized states are dominated by short self-reps
   (trivial) vs extended autocatalytic sets (interesting).
4. Gillespie/SSA on polymer nets for stochastic nucleation + extinction
   (P(nucleation) as an order parameter for the exclusive regime).

# ============================ MILESTONE 5 ============================
# 2-D critical manifold test + Gillespie SSA.

## 2026-10-04 — stoichiometry bugfix

`A+A -> AA` ligations were recorded as {A:1} -> {AA:1} (dict literal
collapsed the duplicate key) — consumed only ONE A. Fixed: reactants and
cleavage products now get nu=2 for a==b splits. random_crn unaffected
(samples without replacement). Fast-propensity path extended to expand
nu copies into the padded index array (x[i]*x[i] == x[i]**2), so nu=2
reactions keep the vectorized path. Verified bit-identical vs slow loop.

## 2026-10-04 — gillespie.py implemented (direct SSA)

Counts-space direct method: reaction propensities = net.propensity(n/V)*V,
food influx = feed_rate*V/species, dilution = per-molecule death.
Validated: food-only mean 994 vs expected 1000; rate=0 catalyst-required
reaction NEVER nucleates (exact chicken-and-egg at finite copy number);
self-rep channel with small spontaneous rate nucleates and grows.
Tests in tests/test_milestone5.py (4 passing).

## 2026-10-04 — THE 2-D MANIFOLD RESULT (headline)

**Cmd:** `python experiments/manifold_2d.py --p-cat-grid 0.2 0.4 0.6 0.8 --p-info-grid 0.0 0.25 0.5 0.75 1.0 --n-reps 4 --n-ligations 200 --n-bootstrap 50 --out data/manifold_2d.jsonl`

80 runs (4x5 grid x 4 reps), polymer_ligation, alphabet=2 max_len=10
food_len=2 k_sat=1.0 feed=5 dilution=1 t_max=60.

**Result — hit-rate tables P(observable ON), rows=p_cat cols=p_info:**

```
phi>0.1:                S>0.5:                  H>0.5:
    .00 .25 .50 .75 1.0     .00 .25 .50 .75 1.0     .00 .25 .50 .75 1.0
0.2  0   0   0   0  .5  0.2  0   1   1   1   1  0.2  0  .5   1   1   1
0.4  0   0   1   1   1  0.4 .75  1   1   1   1  0.4  0  .5   1   1   1
0.6  0  .5   1   1   1  0.6 .25  1   1   1   1  0.6  0 .25   1   1   1
0.8 .5   1   1   1   1  0.8  1   1   1   1   1  0.8  0   0   1   1   1

R0>1 & n_core>0 ~ S-map (any motif suffices). A>0.5 & E>0.3 only in the
upper-right corner (need BOTH high). tail_mass always on (degenerate).
```

**Boundary fits (p_info needed at p_cat = 0.2/0.4/0.6/0.8):**
- phi:    0.88, 0.56, 0.24, —   -> CURVED boundary, substitutability
- H:      ~0.25 flat across p_cat -> pure information-axis threshold
- S,R0,core: ~0.0 (on everywhere except the no-info low-connectivity corner)
- A,E:    only p_cat>=0.6-0.8 AND p_info>=0.5

**Interpretation:**
1. The boundaries are clearly NOT one shared curve — Model B (single
   transition where everything co-emerges) is refuted in this regime.
2. The natural ordering along any path from origin to corner:
   {S,R0,core: any autocatalytic motif} -> {H: template flux, p_info~.25}
   -> {phi: closure flux, needs p_cat AND p_info} -> {A,E: both high}.
   This is the CASCADE/NESTED structure (Models A/C).
3. **The standout cell p_cat=0.8, p_info=0: phi~0.5, cores exist,
   H=0 — organized autocatalytic chemistry with ZERO heredity.**
   Dense random catalysis closes RAFs without templates. Organization
   without information is a real regime — the decoupling the
   experiment was designed to find.
4. phi's boundary slopes DOWN with p_cat: connectivity and information
   are SUBSTITUTABLE for closure — consistent with a manifold like
   p_cat*(a + b*p_info) or similar, but each observable sits on its
   own surface, not the shared one.
5. Honest caveat: these are activation boundaries of observable
   thresholds (smooth logistic surfaces), not proven singularities.
   S/R0 are structurally "easy" (any motif). The strongest claim:
   the properties have DIFFERENT activation surfaces in parameter
   space, so a 1-D sweep's A-vs-B verdict is slice-dependent —
   resolving A vs B requires the 2-D picture, which favors cascade/nested.

## Open threads

1. Exclusive/nucleation regime on polymer nets: p_cat_required=1 +
   no seeding — is there bistability along p_info too?
2. Which species carry organized flux — short self-reps or extended
   autocatalytic sets? (trivial vs interesting organization)
3. Tau-leaping for SSA on polymer-scale nets; P(nucleation) observable.
4. random_crn 2-D analog: p_cat x p_cat_required or p_cat x feed_rate —
   does its manifold also show substitutable axes?

## 2026-10-04 — QUESTION 0: is it a bifurcation or a crossover?

User pointed out the A/B/C question presupposes a genuine bifurcation.
Built experiments/criticality.py measuring, per network along a lambda
slice: phi + var(phi) + Sarle BC, jacobian_spectral_abscissa alpha at
the steady state, invasion_pressure, and within-network seeded-vs-
unseeded bistability.

**Diagnostic bug found+fixed:** naive full-space Jacobian reported
alpha > 0 at organized steady states — artifact of finite-differencing
at EXTINCT species under the max(x,0) clamp (one-sided derivatives at
x=0 create spurious unstable directions). Verified: phi sits rock-
steady 400+ time units; resident-restricted Jacobian gives alpha < 0
and zero dead-species invasion pressure. jacobian_spectral_abscissa now
restricts to the support (alive species); invasion_pressure added.

**Cmd A:** `python experiments/criticality.py --generator polymer_ligation --lambda-param p_info --p-cat-fixed 0.6 --lam-grid 0.0 0.1 0.2 0.3 0.5 0.7 1.0 --n-reps 5 --n-bistability 3 --out data/criticality_polymer2.jsonl`

**Result A (p_info sweep, corrected alpha):**
  lam   <phi>  var    <alpha>  frac_bist
  0.00  0.010  0.000  -0.62    0
  0.10  0.042  0.001  -0.33    0
  0.20  0.100  0.005  -0.27    0
  0.30  0.115  0.003  -0.17    0
  0.50  0.140  0.002  -0.14    0
  0.70  0.224  0.004  -0.17    0
  1.00  0.275  0.006  -0.11    0

- Organized states are GENUINE stable attractors (alpha<0, invasion-
  proof) — real dynamical regimes, not transients. Good.
- BUT alpha never approaches 0, var(phi) rises smoothly (no
  susceptibility peak), zero within-network bistability ->
  NO dynamical-bifurcation signature along p_info. The transition is
  ensemble-level/structural: which attractors exist and are reachable,
  not a bifurcation of a fixed system's dynamics.

**Cmd B:** `python experiments/criticality.py --generator polymer_ligation --lambda-param p_info --p-cat-fixed 0.6 --p-cat-required 1.0 --lam-grid 0.1 0.3 0.5 0.7 --n-reps 4 --n-bistability 4 --out data/criticality_polymer_excl.jsonl`

**Result B (exclusive regime):** STILL frac_bistable=0. Reason found:
the CLEAVAGE BACKDOOR — a catalyst-required product P can be produced
spontaneously as a cleavage fragment of any longer string that forms
without catalysis. All-splits hydrolysis diffuses the system through
sequence space, washing out nucleation barriers. To get bistability in
polymer chemistry we'd need exclusive sets closed under cleavage
(catalysts whose superstrings are also unreachable) — a new knob.

**Question-0 answer so far:**
- Within-system dynamical bifurcation: PROVEN possible in this model
  family (exclusive_catalysis random_crn, M3: seeded/unseeded split,
  feed-ramp hysteresis). NOT found along polymer structural axes —
  cleavage backdoor prevents barriers.
- Ensemble/structural transition: YES everywhere — the fraction of
  networks reaching organized attractors switches sharply with
  (p_cat, p_info). This is a percolation-like transition in the SPACE
  OF SYSTEMS (which attractors exist), not in one system's dynamics.
- Conceptual clarification worth keeping: 'phase transition' math
  applies at two levels. (a) Bifurcation theory applies to a FIXED
  network's dynamics — needs hysteresis/gap-closing; only exclusive
  random_crn shows it. (b) The chemistry->life 'transition' as we
  measure it is level-(b): a landscape/percolation transition where
  new attractor TYPES become reachable. Both are legitimate physics
  but they are different objects, and claims should name which.

## 2026-10-04 — SSA restores the nucleation barrier (finite copy number)

Question: does the cleavage backdoor need artificial blocking, or does
discreteness close it naturally? Ran exclusive-regime polymer nets
(p_cat_required=1.0, all organizing unseeded in ODE) under direct SSA,
seeded (ODE steady state x V) vs unseeded (food only), V in {0.5,1,2,5},
t_max=150.

**Result (phi at trajectory end, 6 nets x 4 V x seeded/unseeded):**
- 4909 @ V=1: CLEAN BISTABLE (unseeded 0.000 / seeded 0.181) — the
  fragment never materializes as a molecule; window closes by V=2.
- 9100: both dead at V<=1 (demographic extinction of ~1-molecule
  species), organizes unseeded at V>=2 — window below its reach.
- 12127: dead at EVERY V — its ODE phi=0.151 'organized' state lived
  on sub-molecular concentrations: a continuum phantom, never existed.
- 6500: organizes at all V — fragment door wide open in this net.
- 9527, 11118: marginal/noisy.

**Interpretation:**
1. Finite copy number DOES create bistability where the ODE showed
   none — the barrier is real physics, bounded in a VOLUME WINDOW:
   below it, demographic noise kills ~1-molecule organized species;
   above it, the fragment door reopens. No artificial blocking needed;
   blocking would model the wrong thing.
2. SSA is the arbiter of genuine organization: some ODE 'organized'
   states are phantoms living at sub-molecular concentrations.
   Caveat on ALL prior phi measurements at this scale — a fraction
   are probably phantom. Worth a phi_SSA audit of earlier regimes.
3. The measurable knob is door WIDTH, not a binary block:
   cleavage_rate scales the fragment-supply flux; P(nucleation)
   vs fragment flux maps where the barrier actually sits.

**Open thread promoted:** P(nucleation) as an order parameter under
SSA — sweep (V, cleavage_rate) for exclusive polymer nets; the
bistable window's boundaries ARE the answer to 'how much suppression
of the fragment door does the barrier require'.

## 2026-10-04 — P(nucleation) sweep: the barrier is a rate, not a wall

**Cmd:** `python experiments/nucleation.py --seeds 11118 9100 4909 6500 --n-reps 4 --out data/nucleation.jsonl`
240 SSA runs: 4 nets x 4 V x 3 cleavage_rate x (4 unseeded + 1 seeded),
exclusive polymer regime (p_cat_required=1.0), t_max=150.

**Result (P_nuc / phi_seeded per cell):**
- 4909: P_nuc 0.75-1.0 almost everywhere including V=1 clv=0.02 —
  the earlier 'clean bistable' single run was stochastic luck
  (P~0.25 that particular draw). NOT a hard barrier.
- 9100: seeded phi~0 at nearly every cell — below the EXTINCTION
  boundary, not bistable. Only nucleates at widest door (clv=0.2, V>=2).
- 11118, 6500: P_nuc ~0.75-1.0 everywhere — door open at all V.

**Corrected interpretation:** discreteness does NOT make the barrier
binary — it makes nucleation a POISSON PROCESS with rate
gamma(V, cleavage_rate, net). 'Bistability' is only meaningful as a
timescale statement: 1/gamma_nuc >> t_obs AND 1/gamma_ext >> t_obs.
For these nets/gammas the window is narrow and net-dependent; most
cells show P_nuc ~ O(1) at t_max=150.
- The single-run M5 'bistable' was real for that draw but not the
  cell's typical behavior — single-run probes overstate hard barriers.
- Extinction and nucleation boundaries both vary with V; some nets
  (9100) have NO window: organized state can't persist at small V
  even while nucleation is still possible at large V.
- cleavage_rate weakly modulates door width within this range; net
  topology (which motifs exist) dominates over door width.

**What this means for blocking:** neither a binary block NOR pure
discreteness gives a clean barrier — the meaningful experiment is
measuring gamma_nuc(V, cleavage) and gamma_ext(V) directly
(mean first-passage times), which is the next experiment.

## 2026-10-04 — mean first-passage: gamma_nuc and gamma_ext measured

**Cmd:** `python experiments/first_passage.py --out data/first_passage.jsonl`
192 SSA runs: 3 nets x 4 V x 2 cleavage x 2 conditions x 4 reps,
t_max=120, phi threshold 0.05, record_dt=0.5.

**Result (tau_nuc / tau_ext, time units; inf = censored):**
- 4909, 11118: tau_nuc = 0 EVERYWHERE — degenerate. These nets carry a
  FOOD-LEVEL self-rep (a food dimer catalyzing its own production
  satisfies the RAF criterion trivially); 'organized' at t=0. Their
  earlier 'bistability' was against a trivial organization.
- 9100 (no food self-rep — real structure):
      clv=0.02:  V=0.5: 32/75   V=1: 17/113   V=2: 8/147   V=5: 6/117
      clv=0.2:   V=0.5: 7/155   V=1: 4/116    V=2: 3/157   V=5: 2/477

**Interpretation:**
1. The rates behave like physics says they should:
   - tau_nuc DECREASES with V (more fragment molecules -> faster
     nucleation) and with cleavage_rate (wider door -> faster).
   - tau_ext INCREASES with V (bigger organized populations are more
     stable against demographic noise).
2. For 9100 the barrier is WEAK: tau_nuc ~ 2-30 vs t_obs=120 —
   organization nucleates quickly wherever it can persist. The
   window tau_nuc >> t_obs essentially doesn't exist at these
   conditions; bistability would need a sparser sequence space
   (bigger alphabet -> rarer fragments) or narrower door.
3. tau_ext is a measurable rate too (75 -> 477 with V) — extinction
   is not a boundary but a timescale, matching theory.
4. Diagnostic value: phi(t) trajectories exposed food-level self-reps
   that made earlier 'bistability' readings partly trivial. RAF-based
   phi counts them; a 'nontrivial-organization' variant of phi (excl.
   motifs entirely inside the food set) may be worth adding.

**Prediction for next run:** alphabet_size 3-4 should steepen tau_nuc
sharply (fragment space grows exponentially with alphabet) — the
bistable window should open where it couldn't at alphabet=2.

**Verification:** raf_reactions audit — 4909 has 8 and 11118 has 17
RAF reactions catalyzed by FOOD species (supplied by feed from t=0);
9100 has 0. The food-catalyst mechanism for instant organization is
confirmed, not inferred.
