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
