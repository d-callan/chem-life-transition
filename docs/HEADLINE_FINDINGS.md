# Headline Scientific Findings — chemtrans 2.0 (Post-Audit Status)

**Generated**: Post-Remediation Audit (Round 7 / `r6` Validation Epoch)  
**Code Version**: `4616bc4` + Round 7 audit fixes  
**Data Artifacts**: `data/validation/r6/` (21 campaigns regenerated under unclipped vector field, exact event-time SSA, order-invariant clustering, and anti-pooling inference guards)

---

## Executive Summary: Is there a Phase Transition?

The central question of this repository is whether the boundary between non-living chemistry and organized biochemistry can be meaningfully characterized as a physical phase transition.

Based on the full suite of audited and regenerated `r6` campaigns, the answer is:

1. **Ignition is real, reproducible, and threshold-governed per network**:
   Individual reaction networks exhibit reproducible thresholds along dilution ($D$) and chemical driving ($\Delta\mu$). An organized, self-maintaining RAF state exists below a family-specific critical dilution $D_c$, collapsing to an absorbing dead state when dilution exceeds catalytic replenishment.
2. **There is NO universal critical point ($D_c$) and NO co-criticality across life criteria**:
   Across all 21 regenerated campaigns, **zero** observable pairs show coincident critical transitions under paired, family-preserving tests. Every pairwise cell in `summary.json` is `incomparable`, `underdetermined`, or `insufficient_families`. Different life-vector observables ($O_{\text{maint}}$, $O_{\text{hered}}$, $O_{\text{persist}}$, $O_{\text{assembly}}$) ignite at disparate dilutions or fail to co-occur within the same network.
3. **The ensemble boundary is a quenched-disorder crossover band, not a sharp line**:
   Rather than a shared critical temperature $T_c$, the chemistry-to-biochemistry boundary in random polymer ensembles is best described as a **percolation-like ignition probability** across catalytic efficacy $\kappa$:
   - At $\kappa = 100$: **0/16 (0%)** of networks ignite ($O_{\text{hered}} > 0$).
   - At $\kappa = 500$: **9/16 (56.3%)** of networks ignite.
   - At $\kappa = 1500$: **9/16 (56.3%)** of networks ignite.
   Ignition thresholds scatter widely across network architectures ($D_c \in [0.2, 0.9]$), reflecting frozen disorder in the randomly generated kinetic topologies.
4. **Transition classes are decoupled and heterogeneous**:
   Where onsets occur, they do not belong to a single universality class. Attractor diversity ($O_{\text{hered}}$) exhibits discrete, step-like onsets; maintenance mass ($O_{\text{maint}}$) exhibits smooth continuous crossovers; and persistence ($O_{\text{persist}}$) exhibits stochastic metastability.

---

## Detailed Comparison: Notebook Claims vs. Audited Evidence

| # | Scientific Claim | Notebook / Historical Status | Audited Status (`r6` Data) | Verdict / Reason |
|---|---|---|---|---|
| **1** | **Ordered Co-Critical Cascade** | Claimed in E9: ordered cascade $D_c(\text{maint}) \sim 1.0 > D_c(\text{persist}) \sim 0.26 > D_c(\text{hered}) \sim 0.1$. | **Refuted as general claim.** All 21 campaigns yield `incomparable` or `insufficient_families`. Observed diffs do not replicate across network seeds. | **UNSUPPORTED / OVERCLAIM**. Reflects network-specific idiosyncrasies, not an ensemble cascade law. |
| **2** | **2D Co-Criticality & Tracking Across $\kappa$** | Claimed in E9 / `n8_zoom`: $O_{\text{assembly}} \mid O_{\text{hered}}$ "tracking" verdict across $\kappa \in \{100, 500, 1500\}$. | **Refuted as artifact of pooling.** Old tracking drew on net5 at $\kappa=100,500$ and net6 at $\kappa=1500$. Under powered-up $N=16$ runs with multi-level guards, $n_{\text{fam\_multilevel}} = 0$. | **REFUTED**. The claim was fabricated by cross-family pooling. |
| **3** | **Continuous Transcritical Transition ($\beta=1$) at $\epsilon=0$** | Claimed in E2b: `hinge` model selected on $m_{\text{poly}}$ with linear turn-on ($\beta=1$) at $D_c \approx 1.13$. | **Partially supported as functional fit; overclaimed as statistical certification.** Turn-on is linear for specific RAF+ nets, but Gaussian AICc on noise-free ODE data cannot provide formal likelihood odds. | **QUALIFIED / APPROXIMATION ONLY**. Model ranking is informative, but p-values and odds ratios are unjustified. |
| **4** | **Finite-Volume Extinction Scaling** | Claimed in E3: extinction time grows strongly with volume ($V=4, 16, 64$), demonstrating first-order metastability. | **Supported on single seed.** Event-time Kaplan-Meier medians scale exponentially with volume ($R^2 = 0.9995$, $\tau_{50} = 44.5\text{s}$ at $V=4$, $54.6\text{s}$ at $V=16$, $112.6\text{s}$ at $V=64$). | **CONFIRMED (SINGLE SEED)**. Mechanism holds for seed `e2coop|s3`; multi-seed generality remains unmeasured. |
| **5** | **Critical Slowing Down (CSD)** | Claimed in E2b: leading Jacobian eigenvalue $\to 0^-$ monotonically approaching $D_c$. | **Recorded but unfitted.** Eigenvalues soften toward zero, but points are too sparse near bifurcation to fit a scaling exponent $\nu$. | **QUALIFIED / UNCALIBRATED**. Eyeballed trend; no formal critical exponent extractable. |
| **6** | **Compositional Multistability & Heredity ($p=2$)** | Claimed in E2: multiple attractors separated by absorbing dead state; lineages survive serial transfer. | **Confirmed.** Under order-invariant single-linkage clustering and mass floor ($m > 0.01$), 56% of networks at $\kappa \ge 500$ exhibit genuine multistability. | **SUPPORTED**. Cooperative catalysis ($p=2$) reliably breaks monostability. |
| **7** | **Assembly Theory Complexity Signature** | Claimed in E7/assembly_l6: zero high-assembly molecules observed at $L \le 6$ (null result). | **Confirmed null, but conceptually inverted.** At $L \le 6$, maximum assembly index requires zero motif reuse, which is mathematically anticorrelated with autocatalytic repetition. | **STRUCTURAL NULL BY DESIGN**. The observable is self-defeating at short polymer lengths ($L \le 6$). |

---

## Critical Methodological & Implementation Bugs Resolved

During the audit, several critical implementation defects and statistical traps were uncovered that directly impacted scientific conclusions:

1. **Clipped Vector Field Artifacts (§1.1)**:
   - *Defect*: `np.maximum(x, 0)` inside ODE and Jacobian evaluations concealed unphysical negative excursions and fabricated spurious limit cycles at low dilution ($e2d$, $pilot$).
   - *Resolution*: Removed all clips inside solver evaluations; excursions now trigger explicit `UNRESOLVED` status with `neg_excursion` diagnostics. Fake limit cycles reverted to certified `STATIONARY` equilibria.
2. **False 2D Tracking via Heterogeneous Pooling (§2.13, §2.16)**:
   - *Defect*: Campaign summaries flattened secondary sweep axes, combining paired differences from completely different network families at different $\kappa$ levels to fabricate "tracking" verdicts on $n \approx 1\text{--}2$ usable seeds.
   - *Resolution*: Enforced `min_pairs = 3` floor and `n_families_multilevel \ge 2` guard. 1D reports are suppressed when secondary axes are declared.
3. **Circular Assembly Theory Null Calibration (§1.2)**:
   - *Defect*: `derive_a_cut` fitted null background spectra to the same organized equilibrium state being scored.
   - *Resolution*: Defaulted to independent abiotic reference (`without_catalysis()`), isolating genuine catalytic enrichment.
4. **Seeded SSA Extinction Blindness (Round 7 Fix)**:
   - *Defect*: SSA event-time observers did not initialize `st["t_first"] = 0.0` for dense-seeded starts at $t=0$. Because extinction was guarded by `t_first is not None`, extinction events were never recorded, causing all seeded runs at $V \ge 16$ to be falsely marked as surviving until $t_{\text{max}}$.
   - *Resolution*: Initialized `st["t_first"] = 0.0` when `st["org"]` is active at $t=0$. Extinction scaling recovered with $R^2 = 0.9995$.
5. **Censoring Bracket Overconfidence (Round 7 Fix)**:
   - *Defect*: `lc_censor_bracket` collapsed to a single point `[lcb, lcb]` if one extreme imputation eliminated the critical point.
   - *Resolution*: Requires both bounding fits to detect transitions; otherwise flags `censor_robust = False` ("transition destroyed under extreme imputation").
6. **Attractor Clustering Order Dependence (§2.5)**:
   - *Defect*: First-center greedy clustering caused attractor counts and heredity indices to depend on initial-condition array ordering.
   - *Resolution*: Replaced with deterministic single-linkage connected-component clustering via union-find with medoid representatives.

---

## What Would Be Needed to Prove a True Phase Transition?

To advance beyond "quenched-disorder finite-size crossover" to a rigorous physical phase transition, future work must establish:

1. **Finite-Size Scaling Series ($L \to \infty$)**:
   - System size in polymer chemistry is governed by max chain length $L$ and alphabet size $|\Sigma|$ ($N_{\text{species}} \sim |\Sigma|^L$).
   - A phase transition requires demonstrating that the crossover band $\Delta D_c(L)$ narrows as $L \to \infty$, and that maximum response susceptibilities diverge according to finite-size scaling forms $L^{\gamma/\nu}$.
2. **Ensemble Extinction Scaling**:
   - The $R^2 = 0.9995$ exponential volume scaling currently confirmed on seed `e2coop|s3` must be measured across an ensemble of $\ge 30$ independent network topologies to determine whether the metastability barrier is generic or seed-dependent.
3. **Uncoupled Order Parameters**:
   - Future experimental protocols must retire `mass_raf` as a continuous order parameter and standardize on $m_{\text{poly}}$ (non-food monomer mass concentration) or catalytic flux $\Phi_{\text{cat}}$ to prevent artificial saturation at $\approx 1.0$.

---

## Round 7+ Update: Mechanism Deep-Dives — Can Recognition, Folding, or Compartments Prune the Soup?

**Epoch**: `r7` | **Commits**: `1192868`, `416184c`, `e8a91bb`, `d790d6d`, `e605bf2`
**Question**: after the audited campaigns established ignition/multistability, can *any* mechanism in this model family concentrate a winning sequence lineage — the step from high-assembly soup to life-like sequence selection?

### Summary of mechanism tests (all cheap scouts + one campaign, K2/L6 unless noted)

| Mechanism | Implementation | Result | Verdict |
|---|---|---|---|
| Base-pairing templated ligation | `templated_ligation`, `delta_g_pair` sweep (`template_lig_l6`, 80 cells) | Sparsity *dropped* 0.85→0.05; N_eff rose ~2→52 | **INVERTED**: complementarity is a symmetry, not a selector — every sequence has a complement, catalysis equalizes |
| Folding-competence gate | `f_fold(T)` stem-loop Boltzmann weight on template edges (`folding_l6`, 40 cells) | Flat across delta_g_fold 0→5 | **NULL**: gate restricts who templates, not what gets made; coverage stays complete |
| Non-pairing monomer C | K=3 alphabet, `purify_l5`/`K3/L4` scout | templated purity 0.46 < blind 0.88 | **NULL**: blind topological lottery out-purifies recognition |
| Recognition stringency | uncapped w_max, perfect-only w_min, kappa 50↔500 | identical flat distributions | **NULL**: not a cap/saturation artifact — complete symmetric coverage |
| Initial-condition asymmetry | 8 asymmetric ICs, 10-vesicle encapsulation scout | all converge to one attractor | **NULL**: single dominant basin; passive trapping is not a bottleneck |
| Osmotic protocell competition | `protocell_population` ODE interior (lipid-stealing, division, lysis) | real divisions/lysis, speciation ~1e-10 | **NULL for speciation**: reward is colligative/sequence-blind |
| SSA interiors (discrete counts) | `interior="ssa"` exact Gillespie per vesicle | **speciation 0.5–0.67**; extinctions real | **PARTIAL**: true founder-drift divergence — but winners are luck, not function |
| Identity-coupled growth | `mu_func` fold-weighted growth term (heuristic) | pure → extinction; balanced → first-mover collapse | **NULL**: no sustained identity differential exists to select on; amplifies founder luck |
| Poisson trap statistics | analytic + equilibrium concentrations | ≥1 MA≥5 molecule: ~3 vesicles; complementary pair seed: ~120 | **MIXED**: single molecules cheap to trap but decay; viable seeds 40x rarer; scaling is quadratic in dilution, exponential in length |

### The two durable positive findings

1. **Semi-permeable retention prunes deterministically.** Removing uniform CSTR dilution and trapping L≥3 species shifts interior sequence sparsity from ~0.04 to ~0.7+ (at long timescale). Selection *by boundary conditions* — convergent (same winner in every cell) and parameter-sensitive (growth dilution at mu=40 flattens it; purification is a transient that needs time).
2. **Discrete molecule counts create real founder effects.** SSA interiors produce speciation ~0.5 vs ODE's ~1e-10 — extinction and lineage lock-in are genuine absorbing-state physics unavailable to mean field.

### The distilled negative result

Every mechanism tested selects *something* — but never for sustained sequence function. Osmosis is colligative (blind), encapsulation is drift (lucky), identity-coupled growth amplifies first-mover luck into collapse. The missing ingredient is **heredity within lineages**: nothing makes a daughter's interior resemble its mother long enough for cumulative selection. This is an architectural gap (lineage tracking across generations), not a parameter.

### The Assembly-Theory answer

Protocells did **not** drive MA≥5 species into detectable abundance — they eliminated them (trapped molecules hydrolyze; sub-single-molecule species are never inside to begin with). Chance encapsulation of a viable complementary-pair seed costs ~10² vesicles at L=6 but scales as ~1/pool² — thermodynamically impossible at realistic replicator lengths. Compartments cannot bootstrap the high-assembly transition; they can only amplify a lineage that already regenerates its own tail.
