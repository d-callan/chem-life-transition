# Scientific & Statistical Review — chem-life-transition (chemtrans 2.0)
**Audit Cycle 4 — Post-Remediation Verification & Deep Forensic Pass**

Scope: Comprehensive review of `src/`, `experiments/`, `tests/`, `runs/NOTEBOOK.md`, `docs/`, and committed `data/` artifacts. 44 critical bugs, inferential flaws, and simulation vulnerabilities have been remediated across four rounds and verified with 79 passing unit tests. This document tracks the **44 verified remediations in the closed ledger**, audits the **active defect inventory**, and details **remaining high- and critical-impact scientific and statistical issues** directly affecting the validity of results.

---

## Executive Summary & Audit State

| Category | Rounds 1–3 | Resolved in Round 4 | Remaining Active |
|---|:---:|:---:|:---:|
| **Correctness Bugs** | 20 | 5 | 1 |
| **Methodological Flaws** | 12 | 4 | 13 |
| **Statistical Weaknesses** | 2 | 0 | 7 |
| **Total Items Tracked** | **34** | **10** | **21** |

> **Round 4 resolves the entire active correctness section** (vector-field clipping, SSA cooperative-catalysis semantics, benchmark invariant violations, `kinetic_accessible` direction-blindness) plus four methodological items (manifest thresholds, competition tie-breaks, co-criticality guards and secondary-axis flattening). What remains is dominated by *statistical* weaknesses — small samples, absent uncertainty quantification, permissive equivalence margins — and methodological issues requiring *reruns* rather than code fixes.

### Verified Remediations (Closed Items)
The following 44 items have been implemented and verified with tests
(79 tests passing at Round-4 close):

**Round 1 Remediations (19 items)**:
1. **`phi_cat` power scaling (§6.1)**: `closure.phi_decomposition` now computes catalytic boost via `C @ x**cat_power`, resolving the >100% flux anomaly on $p=2$ cooperative campaigns.
2. **`kinetic_accessible` reachability (§6.2)**: Implemented full food-generated reachability closure that properly includes non-synthesis channels (`synth=False`, e.g. deactivation) at $\varepsilon = 0$.
3. **`cat_max_len` self-replication leak (§6.3)**: Added upper-bound check `len(PA) <= cat_max_len` to `pick_catalyst` and excluded product self-selection from the uniform random fallback pool.
4. **`run_id` task hash collision (§6.4)**: `run_id` now digests the full task dictionary (including `assembly`, `competition`, `continuation`), preventing silent cache-skipping on campaign resume.
5. **ODE solver divergence crash loop (§6.5)**: Trajectory failures now emit strict-JSON-compliant `None`/zero fields with status `UNRESOLVED`, eliminating campaign worker crashes.
6. **OSCILLATORY classification (§6.6)**: Runaway divergence is evaluated before oscillation, and oscillatory certification requires half-window mean agreement (<5%) and low time-correlation (<0.8), eliminating false positives from monotonic drift.
7. **Absorbing-limit test rigor (§6.7)**: Replaced vacuous food-only rate test with a rigorous test verifying zero uncatalyzed synthesis when non-food polymers are present without catalysts.
8. **`eq_mass` units conversion (§6.8)**: Divided monomer-mass fraction signatures by per-species monomer mass $L$ before converting to discrete copy numbers, achieving true equal-mass lineage seeding.
9. **Assembly Theory copy number formula (§6.9)**: Added missing $1/L$ factor to `derive_a_cut`, correctly converting monomer-mass concentration into molecular sequence copy numbers.
10. **Hysteresis warm-start at $\varepsilon = 0$ (§6.10)**: Reverse continuation sweep now reseeds from an active configuration if the forward sweep ended in the absorbing dead state.
11. **Dwell time truncation inversion (§6.11)**: `dwell_organized` now evaluates duration against the full planned horizon $t_{\text{max}}$, preventing truncated runs from reporting 100% dwell.
12. **Serial dilution cycle drift (§6.12)**: Replaced fixed 4-cycle lag (`starts[-1] - starts[-5]`) with consecutive-cycle drift across the tail, preventing period-2 and period-4 cycles from being certified as stationary.
13. **`ramp` collinearity (§6.13)**: Re-parameterized candidate model `ramp` as $c_0 + c_1 \lambda$ with $k=2$ degrees of freedom.
14. **Control axis negation in summary fits (§6.14)**: Applied control axis direction (`sign * lams`) to `transition_fit`, enabling positive turn-on models (`hinge`, `itrans`) to fit decreasing dilution curves.
15. **Lineage parsing case-insensitivity (§6.15)**: Added `.upper()` across all delimiter branches in `measure.py`.
16. **2D tracking single-family requirement (§6.16)**: `coincidence_2d` now verifies that paired differences across parameter levels originate from the *same* network family (`n_families_multilevel >= 1`).
17. **Small-$N$ sign-aware separation (§6.17)**: `coincidence_report` now checks directional sign consistency, classifying bidirectional scatter as `discordant` rather than `separated`.
18. **Competition extinction reporting (§6.18)**: `_competition_block` now accepts explicit `attractor_pair` parameters and outputs `p_extinct` and `p_other`.
19. **Bounded TRF candidate fitting (§6.19)**: Replaced unbounded Levenberg-Marquardt fits with bounded Trust Region Reflective (`trf`) optimization and expanded quantile start-grids.

**Round 2 Remediations (8 items)**:
20. **Extinction-discounted selection observable ($O_{\text{select}}$) (§6.1)**: Multiplies survivor fixation bias by $(1 - p_{\text{extinct}})$, preventing a population collapse where 29/30 replicates die from reporting maximum Darwinian selection (0.5).
21. **Unpooled, per-family & per-secondary-level transition fits (§6.2)**: `cmd_summarize` now fits candidate models per `(family_id, secondary_level)` group, preventing staggered onsets across networks or $\kappa$ sweeps from blurring into artificial crossover sigmoids.
22. **Amplitude bistability counted in $O_{\text{hered}}$ (§6.3)**: `life_vector` evaluates $\max(n_{\text{comp}}, n_{\text{total}}) - 1$ and emits `amplitude_bistable: bool`, preserving amplitude bistability.
23. **Observable-targeted zoom window in `cmd_zoom` (§6.4)**: Added `--observable` targeting (defaulting to primary order parameter), preventing zoom windows from blowing out across the union envelope of disparate observable changepoints.
24. **ARPACK partial convergence guard (§6.5)**: `jacobian_eigs` returns `(eigs, fully_converged)`; unstable fixed points with missing modes are rejected from stability certification.
25. **Honest persistence survival count ($n_{\text{surv}}$) (§6.6)**: `life_vector` requires `t_first_organized is not None and t_extinction is None`, preventing completely unorganized SSA runs from reporting 100% survival.
26. **Unpadded censor time in `first_passage_times` (§6.7)**: Censor times now reflect actual simulation termination rather than padded $t_{\text{max}}$, and first passage records the first crossing rather than the end of the dwell window.
27. **Directional enforcement in discrete changepoints (§6.8)**: `_lc_discrete_family` records `onset` vs `offset`, and `coincidence_report` downgrades opposite-direction transitions to `discordant_direction`.

**Round 3 Remediations (7 items)**:
28. **Root-polish failure guard & real sparse convergence status (§6.9)**: `integrate_cstr` now uses `sol_ls.success` on the sparse `least_squares` path (no more hard-coded `polish_ok=True`), and a failed polish revokes a trajectory-flagged `converged`, emitting `UNRESOLVED`/`x_eq=None` with an explicit message instead of certifying an unpolished endpoint.
29. **Direction-specific synthesis gating for thermal hydrolysis (§6.10)**: New `Channel.synth_dir` field (`"fwd"|"rev"|"both"`) selects which direction carries the $\varepsilon$ gate. Thermal `PA <=> P + A` uses `synth_dir="rev"`: hydrolysis runs unconditionally at $\varepsilon=0$ (polymers are no longer immortal) while condensation stays gated (the dead state remains absorbing). A bare `synth=False` — the original spec — was rejected because it would have un-gated condensation and broken the absorbing limit in the opposite direction. Semantics propagate identically through `rates()`, `d_rates()`, the sparse path, `gillespie.build_processes`, and `closure.kinetic_accessible`.
30. **Jump-verified fold classification (§6.11)**: `estimate_lc_continuous` now inspects the fitted discontinuity $c_1 - c_0$. A jump $\ge 10\%$ of the data span classifies as `saddle_node`; a vanishing or negative jump classifies as `continuous_half_power` (mean-field $\beta = 1/2$ second-order), with `fold_jump` recorded in the estimate detail.
31. **Declared-driven irreversibles in detailed-balance validation (§6.12)**: `check_detailed_balance` no longer skips $k_r=0$ channels: it reports `error: "irreversible_channel"` on every one and fails `ok=False` unless the channel is explicitly declared `Channel(driven=True)` — a labeled driven-limit boundary process. `transcritical_motif` ($k_{\rm eq}=\infty$) declares `driven=True`, so `benchmark.json` validation still passes while undeclared irreversibles cannot silently validate. New fields serialize through `to_dict`/`from_dict`/`without_catalysis`.
32. **Genuine hysteresis detection in `_continuation_block` (§6.13)**: Forward/reverse branches are now aligned at common $\lambda$ values and compared by `m_poly`. A single-valued curve that crosses `mass_floor` produces one fwd flip and one rev flip one grid step apart — a threshold bracket, not a fold pair. `hysteresis_width` is reported only when branch states actually diverge (`> 2 × mass_floor` at a shared $\lambda$), with `max_branch_diff`, `n_bistable_points`, and `bistable_window` recorded; the verified monostable `transcritical_motif` now reports `hysteresis_width=None` (previously it fabricated $\approx 0.2$).
33. **Graceful serial-dilution solver failure (§6.14)**: `sd_map` returns `(starts, ends, ok, message)` instead of raising on `solve_ivp` failure; `measure_run` emits a strict-JSON `UNRESOLVED` record with `n_cycles_completed`; `serial_transfer_test` returns partial cycles with `status: "UNRESOLVED"`; `sd_fixed_point_multipliers` returns `None` rather than a garbage finite-difference Jacobian.
34. **Generator-sweep validation coverage (§1.2)**: `cmd_validate` now deduplicates on the resolved generator spec (`sha256_json(gen)`) rather than `seed_key`, so every distinct network produced by a `generator.*` sweep is checked for conservation and thermodynamic consistency (benchmark.json went from 1 validated network to all 7 $\kappa$ points; `e1_landscape` validates 120). Also made the letter-drift computation sparse-safe (`S.T @ comp.T`).

**Round 4 Remediations (10 items)**:

35. **Vector-field declipping + explicit negativity detection (§1.1)**: Removed `np.maximum(y, 0.0)` from *inside* solver evaluations throughout `ode.py` (`_krylov_equilibrium` RHS/Jacobian, `integrate_cstr` RHS/Jacobian, root-polish lambdas, `sd_map` and `sd_fixed_point_multipliers` batch integration) and from `spectral.jacobian_eigs` — the vector field and Jacobian are now evaluated at the true state per Invariant 4. `d_rates` catalytic derivatives are unclipped so the Jacobian is the exact derivative of the rate law (verified by finite-difference test at a negative excursion state). Negative excursions are *detected, not hidden*: `integrate_cstr(neg_tol=1e-6)` classifies a trajectory/endpoint with `x_min <= -neg_tol` as `UNRESOLVED`, records `x_min`/`neg_tol`/`neg_excursion` in `ODEResult.diagnostics`, and `measure_run` propagates `x_min`/`neg_excursion` into the run record. Polish candidates with `x.min() < -neg_tol` are rejected as unphysical roots; `_krylov_equilibrium`'s `resid_ok` likewise refuses negative-orthant certification; `sd_map` returns `ok=False` ("unphysical state") rather than clipping batch endpoints. Remaining `np.maximum` sites are documented input-state sanitization (x0 projection, Newton-iterate bounds, reporting clamps) — not solver-evaluation clips. **All historical certified equilibria, spectral eigenvalues, and SD results produced under the clipped field must be treated as suspect and rerun.**
36. **Observable helper correctness (§1.2 items 1–3)**: `null_ma_distribution` no longer truncates alphabets >2 letters (proper `alphabet` override + `chr(65+i)` convention); `rmst` integrates the Kaplan–Meier curve as the piecewise-constant step function it is (exact `S(t_i)·Δt`), eliminating the systematic trapezoidal underestimate; `_krylov_equilibrium` reports the *actual* `nfev` spent in PTC/lsmr (no more hard-coded `+200`).
37. **SSA cooperative-catalysis convention made explicit (§1.3)**: `gillespie.build_processes`/`hazard`/`simulate` accept `hazard_mode`: `"participants"` (literal $p$-molecule falling-factorial $(n_c)_p$ — the historical semantics, kept as default) vs `"functional"` ($\kappa(n_c/V)^p$, matching the ODE's $x_c^p$ functional boost, which does not quench at $n_c=1$). The mode is wired through `measure_run`'s stochastic block and competition assays, recorded in run records and manifests, and is now an auditable modeling choice rather than a silent reinterpretation. **Historical `e3_extinction` "nucleation barrier" claims remain unresolved**: they were generated under `"participants"` semantics and must be rerun under `"functional"` (and ideally both) before any volume-scaling interpretation is defensible.
38. **Schlögl open-system declaration (§1.4)**: `schlogl` now declares `params={"open_system": True}`; `check_letter_conservation` reports declared open systems explicitly (reservoir-emulating channels are a designed clamped-bath device, not a conservation bug) while undeclared violations still fail.
39. **`immigration_death` thermodynamic consistency (§1.5)**: free energies now satisfy $g_X = -\ln(k_{\text{birth}}/k_{\text{death}})$ so $K_{\text{eq}} = k_f/k_r$ exactly; the $k_{\text{death}}=0$ boundary declares `driven=True`. `check_detailed_balance` passes with `rel_err = 0` rather than the previous factor-4 violation.
40. **Bidirectional `kinetic_accessible` (§1.6)**: the closure now fires a channel forward when reactants are reachable *or* in reverse when products are reachable, with direction-specific eps gating (`synth_dir`) respected per direction at $\varepsilon=0$ — matching `core._rate_vectors` semantics. Reverse-supplied intermediates (condensation, ligation) now correctly count as accessible. Structural `_closure`/`raf_reactions` remain intentionally forward-directed (RAF formalism semantics).
41. **Manifest records actual detection parameters (§2.9)**: `write_manifest` now writes the real `theta_organized`, `mass_floor`, and `hazard_mode` from the protocol's stochastic spec instead of the boilerplate `theta_hi=0.6/theta_lo=0.4`. Existing committed manifests remain stale.
42. **Minimum-pairs floor in co-criticality verdicts (§2.13 items 1–2)**: `coincidence_report` accepts `min_pairs` (default 3); pairwise cells with `n_pairs < min_pairs` report `insufficient_families` rather than emitting an equivalence verdict on a single observation. Propagated through `coincidence_2d` per-level reports. Stored `cocrit_2d*` summaries predate the floor and must be regenerated.
43. **Competition undecided label + tie removal (§2.14)**: `_competition_block` no longer awards equidistant/near-equidistant terminal states to "A": undecided mixtures get `undecided`, and counts now record `dead`, `other`, and `undecided` separately with `p_extinct`/`p_undecided` reported. The 80%-commitment boundary and mid-transient fixation issues are documented (§2.5, §2.7 residuals); **all historical competition labels must be re-derived**.
44. **Secondary-axis flattening suppressed (§2.16)**: when `secondary_axis` is declared, `cmd_summarize` no longer emits a flattened 1D `coincidence_report` over pooled secondary levels (the reason is recorded in the summary); `cmd_zoom` brackets per-level family estimates rather than collapsed λc's. The correctly partitioned `coincidence_2d` remains the 2D analysis of record.

## 1. Active Correctness & Code Defects (Residual)

### 1.1 Residual Helper Hazards (dormant / unanalyzed paths)
- **Location**: `src/chemtrans/observables/assembly.py` & `src/chemtrans/observables/spectral.py`
- **Defect** (items 1–3 of the former §1.2 — `null_ma_distribution` alphabet truncation, `rmst` trapezoidal integration, `_krylov_equilibrium` hard-coded `nfev` — resolved in ledger item 36):
  1. `derive_a_cut` circular-null hazard (`assembly.py`, invoked at `measure.py:205`): the "abiotic null" spectrum $M_L \sim c \cdot r^L$ is fit by log-linear regression on **the same run's equilibrium** `x` via `length_mass_spectrum(net, x)` — on an organized run, the null is regressed on the very organized data being scored. A network producing abundant long polymers inflates $r$, raising $L^*$ and thus $a_{cut}$, self-referentially setting the signature bar. No protocol currently sets `assembly.derive_cut`, so the hazard is dormant — but the call site applies it unconditionally to the run under test whenever enabled.
  2. Dead spectral analysis helpers (`spectral.py`): `leading_rate`, `relaxation_rate_series`, and `fit_slowing_down` are never called — `leading_eig_re` is recorded per equilibrium (and cited for "critical slowing down confirmed" in `NOTEBOOK.md`), but no critical exponent $\nu$ is ever fit and no uncertainty computed; the claim rests on eyeballed point values. (Eigenvalues are now computed on the unclipped Jacobian — ledger item 35 — but remain unfitted.)

*All other Section-1 items are resolved — see ledger items 35–40. Resolved entries are removed from the active inventory per audit convention; their scientific consequence (invalidated historical results requiring rerun) is tracked in Section 4 and the ledger notes.*

---

## 2. Active Methodological Weaknesses (Unfixed & New)

### 2.1 Survivorship Bias at the Critical Point
- **Location**: `experiments/campaign.py` (lines 277–301)
- **Defect**: `transition_fit` only includes points where a certified equilibrium exists (`extract_y` returns `None` for non-STATIONARY runs). Near a critical transition point $\lambda_c$, **critical slowing down** causes trajectories to require vastly longer convergence times, frequently exceeding `t_max`.
  *Aggravated by §6.9*: polish-failure revocation is the honest behavior (a trajectory-flagged-but-unpolished endpoint is no longer certified), but it removes exactly the marginal, slow-relaxing points near $\lambda_c$ where polish most often fails — increasing the dropout count at the bifurcation. The tradeoff is correct-but-fewer certified points; it makes explicit uncertainty quantification on $\lambda_c$ (Section 3.2) more pressing, not less.
- **Scientific Impact**: Nonconverged points are dropped precisely at the bifurcation locus where the discriminating signal lives. A sharp discontinuous transition with near-critical dropouts appears smoother, biasing model selection toward smooth crossovers (`sigmoid`).

### 2.2 AICc on Deterministic, Noise-Free ODE Data
- **Location**: `src/chemtrans/analysis/critical_points.py`
- **Defect**: Candidate models are fitted to noise-free ODE steady states. RSS represents functional approximation error (how closely a transcendental manifold matches an empirical shape), not observational variance. Gaussian-error AICc has no formal likelihood foundation here.
  Furthermore, when the number of certified points is small ($n \le k + 1$, e.g. $n = 5$ points for 4-parameter models like `fold`, `itrans`, `sigmoid`), the small-sample denominator $n - k - 1 \le 0$ is clipped to 1:
  ```python
  + 2 * k * (k + 1) / max(n - k - 1, 1)
  ```
  This applies an arbitrary penalty of $+40$, making complex models mathematically incapable of winning on small sweeps regardless of fit quality.

### 2.3 Survival and Finite-Size Scaling Machinery Is Dead Code
- **Location**: `src/chemtrans/observables/persistence.py` & `src/chemtrans/observables/spectral.py` & `runs/NOTEBOOK.md`
- **Defect**: Advanced statistical survival routines (`kaplan_meier`, `nelson_aalen`, `rmst`, `survival_table`, `extinction_scale_fit`, `nucleation_scale_fit`) are never called by `measure_run`, `campaign.py`, or `cmd_summarize`. The same applies to spectral analysis (`relaxation_rate_series`, `fit_slowing_down`): `leading_eig_re` is recorded per equilibrium and cited for "critical slowing down confirmed" (`NOTEBOOK.md`), but no critical exponent $\nu$ is ever fitted and no uncertainty quantified — the claim rests on monotonic point values at $\le$ a handful of grid points.
  Notebook claims that "extinction time grows strongly with volume $V$" are eyeballed from raw means of $n = 3\text{--}8$ replicates without computing survival estimators, fitting scale exponents, or accounting for right-censoring. Furthermore, `extinction_scale_fit` runs ordinary linear regression on $\ln(\tau)$, treating right-censored runs at $t_{\text{max}}$ as exact failure events.

### 2.4 Degeneracy of `mass_raf` as an Order Parameter Along Organized Branches
- **Location**: `src/chemtrans/observables/closure.py` & `AGENTS.md`
- **Defect**: By definition, `mass_raf` measures the fraction of non-food monomer mass residing in species belonging to the maxRAF. Once a system ignites and self-maintains, virtually all synthesized polymers are produced by and participate in the RAF, causing `mass_raf` to saturate at $\approx 1.0$. As admitted in `NOTEBOOK.md`, `mass_raf` is completely insensitive to changes in total biomass, growth rate, or turnover along the organized branch. Total polymer mass `m_poly` is the de facto physical order parameter, but `AGENTS.md` and campaign configurations still treat `mass_raf` as primary.

### 2.5 Attractor Clustering is Order-Dependent and Non-Symmetric ($O_{\text{hered}}$ Flips Based on IC Indexing)
- **Location**: `src/chemtrans/observables/heredity.py` (`cluster_attractors`, lines 80–100)
- **Defect**:
  `cluster_attractors` assigns clusters greedily using the first observed state as the center, without updating centroids. If two equilibria are both within $\text{tol} = 0.05$ of the first IC center in opposite directions, their pairwise distance is $0.10$ ($2 \times \text{tol}$), but both are assigned to Cluster 0.
  Tested and verified numerically: on three equilibria $E_0, E_1, E_2$ where $d(E_0, E_1) = 0.04$, $d(E_0, E_2) = 0.04$, and $d(E_1, E_2) = 0.08$ (with tolerance $\text{tol} = 0.05$):
  - Processing order $[E_0, E_1, E_2]$ outputs **1 attractor** ($O_{\text{hered}} = 0$, monostable).
  - Processing order $[E_1, E_2, E_0]$ outputs **2 attractors** ($O_{\text{hered}} = 1$, multistable)!
- **Scientific Impact**: Multistability and heredity classifications are non-deterministic with respect to initial condition ordering. Permuting the initial condition list in `_ic_bank` can flip the heredity observable between active and inactive.
- **Additional coverage defect (`multistability_scan`, lines 127–133)**: the IC bank is *boundary-only* — every IC is `food_inlet` plus at most 8 nonfood species at concentrations $\le 10^{-3}$. All 12 ICs lie in a thin shell near the dead/food-dominated corner of state space. Interior attractors whose basins do not intersect that shell (e.g., a high-polymer state separated by a basin boundary) are structurally unreachable from this bank — `n_attractors` is a lower bound whose coverage is unknown and never characterized. Combined with order-dependence, $O_{\text{hered}}$ undercounts multistability in an unsystematic way.

### 2.6 Model-Selection Discontinuities Propagated into $\lambda_c$ Without Model Averaging
- **Location**: `src/chemtrans/analysis/cocriticality.py` (`estimate_lc_continuous`)
- **Defect**: The critical point $\lambda_c$ is extracted from whichever model achieves the minimum AICc. On near-tied fits ($\Delta\text{AICc} < 2$), infinitesimal perturbations or sibling network variations cause the selected model to jump between `hinge`, `fold`, and `itrans`. Because these models define $\lambda_c$ differently (cusp vs jump vs inflection), $\lambda_c$ jumps discontinuously. This model-selection variance is unquantified and flows directly into paired difference tests.

### 2.7 Competition SSA Outcomes Ignore Transient Truncation
- **Location**: `src/chemtrans/measure.py` (`_competition_block`)
- **Defect**: Although Round 1 added `p_extinct` and `p_other` to the output schema, replicates that terminate due to hitting `max_events` (`truncated=True`) are still classified into "A", "B", "dead", or "other" based on their state at `censor_time`. There is no verification that the trajectory achieved a quasi-stationary state rather than being intercepted mid-transient.

### 2.8 Krylov Large-Scale Solver Path Cannot Classify Dynamics
- **Location**: `src/chemtrans/dynamics/ode.py` (`_krylov_equilibrium`)
- **Defect**: For networks with $N > 3000$ species, `integrate_cstr` switches to trust-region Krylov equilibrium seeking (`_krylov_equilibrium`). This path seeks roots of the steady-state equation directly and generates no time trajectory. Consequently, oscillation and runaway detection are completely disabled, returning only `STATIONARY` or `UNRESOLVED`. For large polymer networks, `UNRESOLVED` conflates solver budget exhaustion, slow drift, and sustained limit cycles.

### 2.9 `network_joint_bootstrap` Violates Anti-Pooling Doctrine
- **Location**: `src/chemtrans/analysis/critical_points.py` (lines 186–191)
- **Defect**: `network_joint_bootstrap` resamples network records with replacement, but then concatenates their $(\lambda, y)$ points into a single pooled scatter plot and fits one global curve. This violates the repository's foundational principle that heterogeneous networks with different intrinsic critical points must never be pooled into a single regression.

### 2.10 Conceptual Inversion of Assembly Theory at Short String Lengths ($L \le 6$)
- **Location**: `src/chemtrans/observables/assembly.py`
- **Defect**: In `assembly.py`, `assembly_signature_absolute` defaults $a_{\text{cut}}$ to the pool maximum ($a_{\text{max}} = 5$ for $L=6$). For binary strings of length $L \le 6$, an assembly index of 5 requires **zero sub-block reuse**. Any autocatalytic replication or sequence-kin preference ($p_{\text{info}}$) inherently creates repeating motifs (e.g. `AA`, `ABAB`), which *depresses* the assembly index (`ABAB` $\to a=2$, while `ABBA` $\to a=3$). Consequently, the only sequences that achieve maximum assembly are the most disordered, non-autocatalytic strings. Expecting an autocatalytic network to produce high concentrations of maximum-assembly short strings is an unaddressed conceptual contradiction.

### 2.11 Extinct State Compositional Distance Distortion
- **Location**: `src/chemtrans/observables/heredity.py` (`compositional_distance`)
- **Defect**: For an extinct state ($x \approx x_{\text{food}}$), non-food mass is zero, and `attractor_signature` returns a zero vector. The $L_1$ compositional distance between an extinct state and an organized attractor is $0.5 \times \|0 - \text{sig}\|_1 = 0.5$. In contrast, two distinct orthogonal organized attractors (e.g. pure AA vs pure BB) have distance $1.0$. In assays like `serial_transfer_test`, a collapsed/extinct state appears *closer* (distance 0.5) to the reference lineage than an orthogonal rival lineage (distance 1.0).

### 2.12 Underpowered 2D Co-Criticality Evidence & Stale Summaries (RESIDUAL — code guards added, data unchanged)
- **Location**: `src/chemtrans/analysis/cocriticality.py` (`coincidence_2d`) & `data/validation/cocrit_2d*/`
- **Status**: the *code* defects are remediated — `coincidence_report`/`coincidence_2d` now enforce a `min_pairs = 3` floor (ledger item 42) and `n_families_multilevel` guards against cross-level family pooling — so the *generated* failure mode is closed. What remains is the **data-level problem the guards merely expose**:
  1. **Verified on real data** (`data/validation/cocrit_2d_dk_n8_zoom-1e6ed5db10bf/summary.json`): the `O_assembly|O_hered` **"tracking"** verdict is built on exactly 3 paired diffs — `net5` at $\kappa = 100$ and $500$, plus a **different** family `net6` at $\kappa = 1500`. `O_maint|O_hered` "tracking" draws on `net5`/`net7`/`net1`/`net6` — at most 2 of the 8 seeds per level. Under the new `min_pairs` floor these cells now correctly report `insufficient_families` — but that means the campaign as designed cannot support its tracking claim without more seeds per secondary level.
  2. **Stale stored summaries**: the committed `pairwise_2d` entries predate the `n_families_multilevel` and `min_pairs` guards, the `hazard_mode` semantics, and the unclipped vector field — they cannot be audited and must be regenerated.
- **Dead-code landmine**: `campaign.py` line 127 (`rec.setdefault("family_id", f"{seed_key}|{rec['net_hash']}")`) is inert today because `measure_run` always sets `family_id` first — correctly, from pool-defining generator params (seed, alphabet_size, max_len, ligation) so it is stable across kinetic sweeps. But any future record path bypassing `measure_run` would silently inherit `seed_key|net_hash` family IDs, which fragment across generator sweeps (net_hash changes at every level) and would force `n_families_multilevel = 0` — exactly the failure mode the guard exists to detect, but masked as `underdetermined`.
- **Scientific Impact**: "Tracking"/"drift" conclusions from `cocrit_2d*` campaigns are n≈1–2 evidence presented as 8-network results. The honest output under current code is `insufficient_families` — the campaigns need more seeds per secondary level, not more analysis.

### 2.13 Candidate Model Directional Asymmetry in `critical_points.py` Gating Monotonicity (HIGH)
- **Location**: `src/chemtrans/analysis/critical_points.py` (`_bounds`, lines 89–98)
- **Defect**:
  In `_bounds`:
  - `hinge`: $A \in [0, \infty)$ (positive slope only).
  - `itrans`: $A \in [0, \infty)$ (positive slope only).
  - `fold`: $A \in [0, \infty)$ (positive slope only).
  - `sigmoid`: $A \in [-\infty, \infty]$ (can slope upward or downward).
  - `ramp`: $c_1 \in [-\infty, \infty]$ (can slope upward or downward).
  If a swept parameter causes an order parameter to decrease (e.g. dilution washing out polymer mass, fuel concentration falling, or extinction sweeps), `hinge`, `itrans`, and `fold` are mathematically incapable of fitting the data because their amplitudes are constrained to be non-negative.
- **Partial mitigation**: `cmd_summarize` applies the declared `control_axis` sign (`sign * lams`) to transition fits (ledger item 14), so a correctly-configured decreasing-axis protocol presents an *increasing* curve to the candidate set. The residual risk is protocols where the declared sign does not match the data's actual direction — the bounds asymmetry itself is still in the model definitions.
- **Scientific Impact**: Model selection comparisons remain structurally biased toward `sigmoid` on any sweep axis whose declared sign does not match the data's actual direction. The reported victories of `sigmoid` across early campaigns were plausibly artifacts of candidate bound asymmetry rather than evidence for smooth crossover physics; historical model-selection results are invalidated regardless.

---

## 3. Statistical Weaknesses

1. **Small Sample Sizes Across Flagship Campaigns**:
   - Pairwise coincidence matrices rely on $n_{\text{pairs}} = 1\text{--}6$ across all validation campaigns.
   - Network replicates are typically $N = 3\text{--}8$.
   - Gillespie SSA replicates are $n = 3\text{--}8$ per condition.
   - At $n_{\text{pairs}} \le 2$ with tolerance auto-set to $2 \times$ grid spacing ($0.4$ on a $[0.3, 1.5]$ axis), the verdict `coincident_within_tol` is an extremely weak equivalence claim.
2. **Absence of Uncertainty Quantification**:
   - No confidence intervals or bootstrap standard errors on fitted critical points ($\lambda_c$), dwell fractions, or extinction times.
   - Binomial confidence intervals are omitted for competition win probabilities ($p_A = 29/30$ quoted bare).
   - Scaling regressions report only $r^2$ over $n \ge 3$ points without parameter standard errors.
3. **Uncorrected Multiple Testing**:
   - The $15$-cell pairwise coincidence matrix across multiple parameter levels evaluates multiple hypotheses simultaneously without family-wise error rate or false discovery rate corrections.
4. **Arbitrary Threshold Conditioning**:
   - Scientific classifications depend on un-swept, hand-tuned parameters: `min_mass = 0.01`, `tol = 0.05`, `theta_organized = 0.1`, `mass_floor = 0.01`, `dwell = 10.0`, `t_max = 200.0`. No sensitivity analyses exist to verify that conclusions are robust to threshold choices.
5. **Equivalence Margin Tolerance Inflation in Co-Criticality Hypotheses (`tol_pair` = 30–40% of Domain)** (NEW — HIGH):
   - In `cocriticality.py` line 291, `tol_pair` defaults to $2 \times$ median grid resolution. On standard campaign grids (e.g. $D \in [0.05, 1.1]$ with 8 points, $\Delta D \approx 0.15$), `tol_pair` evaluates to $0.30\text{--}0.40$ — spanning **30% to 40% of the entire swept domain**.
   - Two bifurcations occurring at vastly different physical regimes (e.g. maintenance at $D = 0.5$ and heredity at $D = 0.2$, a 2.5-fold difference in dilution rate) have $|\Delta\lambda| = 0.30 \le \text{tol\_pair}$, and the algorithm declares them **`coincident_within_tol`**.
   - The primary evidentiary claim of a unified life phase transition is therefore vulnerable to an excessively permissive equivalence margin that mistakes broad physical separation for co-criticality.
6. **Observation Grid Artifacts**:
   - Event detection, dwell, and first-passage times are evaluated on the discretization recording grid (`record_dt`), not exact continuous event times.
7. **Stale Summaries and Historical Data Invalidation**:
   - Because Round 1, Round 2, and Round 3 corrected `run_id` hashing, catalyst length filtering, `phi_cat`, `derive_a_cut`, `O_select`, `cmd_summarize`, and continuation hysteresis detection, **all committed data in `data/validation/` and numbers quoted in `runs/NOTEBOOK.md` are invalidated**.
   - Resuming campaigns will re-execute all tasks under new run IDs.
   - Committed `summary.json` files on disk are out of sync with current schemas.

---

## 4. Current Status of Headline Scientific Claims

| Headline Claim in Notebook / Docs | Status Following Audit Cycle 4 | Root Cause & Evidence |
|---|:---:|---|
| Continuous transcritical transition ($D_c \approx 1.13, \beta=1$) | **Plausible, but overstated — rerun required** | Supported by Jacobian eigenvalue $\to 0$ and benchmark recovery, but rests on 6-curve AICc on 9 noiseless ODE points with survivorship bias at $\lambda_c$ (§2.1–2.2). All historical equilibria/eigenvalues were produced under the clipped vector field (ledger item 35) and must be regenerated. |
| $p=2$ cooperative catalysis $\to$ compositional multistability | **Confirmed by construction** | Quadratic boost $\kappa x_C^2$ over shared fuel creates a winner-take-all architecture; amplitude bistability is now incorporated (§6.3). Stochastic analog is hazard-mode-dependent (ledger item 37). |
| "Selection asymmetry survives equal-mass control" | **UNSUPPORTED / INVALIDATED** | Historical claim used buggy `eq_mass` units and survivor-only conditioning that masked 90%+ extinction as maximum selection; competition labels additionally suffered an A-favoring tie-break (ledger item 43). Remediations require full campaign re-execution. |
| Ordered cascade ($D_{\text{maint}} > D_{\text{persist}} > D_{\text{hered}}$) | **Fragile** | Based on $N=3$ networks; changepoints are lower-bound counts; persist is a 3–8 rep mean at a single volume; sensitive to small sample noise. |
| 2D "Tracking resolves at finer grid" | **UNDERDETERMINED / SEVERELY UNDERPOWERED** | In the stored `n8_zoom` summary the `O_assembly|O_hered` "tracking" verdict rests on **one family (net5) contributing diffs at 2 of 3 levels plus a different family (net6) at the third**; per-level verdicts are n_pairs=1–2 of the 8 seeds. Under the new `min_pairs` floor (ledger item 42) these now correctly report `insufficient_families`; the flattened 1D matrix is suppressed on 2D protocols (ledger item 44). Stored summaries must be regenerated (§2.12). |
| Exponential extinction scaling with volume $V$ | **UNSUPPORTED / CONVENTION-DEPENDENT** | Historical runs used `hazard_mode="participants"` — the falling-factorial zero at $n_c \le 1$ extinguishes cooperative catalysis at small $V$, manufacturing the apparent nucleation barrier (ledger item 37). Scaling regressions were never executed (§2.3). Must be rerun under `"functional"` (and ideally both modes) before any interpretation. |
| Assembly Theory null result | **Uncertified / Compromised** | $1/L$ copy number bug inflated $A_{\text{cut}}$ in original run; short string lengths conceptually invert assembly complexity vs motif reuse (§2.10); `derive_a_cut` can self-calibrate on the run under test (§1.1). Must rerun. |
| Hysteresis loops at $\varepsilon = 0$ | **Detection machinery now sound; unverified in data** | Reverse reseed operational (§6.10) and hysteresis now requires genuine branch-state divergence at common $\lambda$ (§6.13) — prior one-grid-step artifacts can no longer be reported. The net2 $\approx 0.6$ width in `docs/assembly_and_drivers_review.md` is likely spurious (it predates the guard); needs campaign re-execution. |
| Co-criticality summaries selecting `sigmoid` | **Mechanically resolved in code; residual bias remains** | Dilution sign negation (§6.14) and unpooled per-family fitting (§6.2) are implemented; but candidate-bound asymmetry (§2.13) still prevents `hinge`/`fold`/`itrans` from fitting decreasing curves — `sigmoid` retains a structural advantage on any sweep axis not explicitly negated. Campaign summaries require regeneration. |

---

## 5. Genuine Strengths (Balanced Epistemics)

- **Thermodynamic Consistency**: Wegscheider cycle conditions and detailed balance ($k_f/k_r = \exp(-\Delta G)$) are rigorously maintained channel-by-channel from microscopic species free energies.
- **Falling-Factorial SSA Implementation**: Exact combinatorics for Gillespie stochastic simulation, properly merging identical-reactant hazards, preventing negative species counts, and handling zero-hazard absorbing boundaries.
- **Honest Null Infrastructure**: The architecture explicitly supports null models; `null` is a first-class candidate in AICc model selection, and flat order parameters are treated as legitimate outcomes.
- **Nested Latent Draw Design**: Generator random draws nest monotonically across threshold sweeps ($p_{\text{cat}}, p_{\text{selfrep}}, p_{\text{info}}$), ensuring that parameter sweeps preserve underlying graph topologies.
- **Rigorous Sparse/Dense Equivalence**: Automated unit tests verify bitwise and numerical agreement between dense array operations and sparse CSR/CSC matrix representations.
- **Audit-Gated Provenance**: Manifests enforce protocol hash matching, code version recording, and strict JSON serialization without NaN/Infinity leakage.

---

## 6. Implementation Specifications for Critical & High-Priority Open Issues

The following specifications cover newly identified high- and critical-priority issues. **All six (§6.9–§6.14) are now implemented and test-verified; the specs below are retained as the audit trail of what was required.** Two items were implemented with refined semantics relative to the original spec text: §6.10 uses a per-direction `synth_dir` gate (a bare `synth=False` would have un-gated thermal *condensation* and broken the absorbing limit), and §6.12 uses a declared `Channel.driven` flag so labeled driven-limit fixtures like `transcritical_motif` pass validation while undeclared irreversibles fail.

---

### 6.9 Fix 9: Root-Polish Failure Guard & Sparse Convergence in `integrate_cstr` — **IMPLEMENTED**
- **Target File**: `src/chemtrans/dynamics/ode.py`
- **Function**: `integrate_cstr` (lines 330–355)
- **Current Behavior**:
  ```python
  # Sparse polish:
  xp, polish_ok = np.maximum(sol_ls.x, 0.0), True
  ...
  if polish_ok and np.all(pr <= 1.0) and np.all(np.isfinite(xp)):
      x, resid_end, converged = xp, pr, True

  if converged:
      x_eq = x
      ttype = "STATIONARY"
      message = ("converged (componentwise residual + window"
                 + (", root-polished" if polish else "") + ")")
  ```
- **The Defect**:
  1. `polish_ok` is hard-coded to `True` for sparse networks even if `least_squares` failed.
  2. If trajectory integration set `converged = True` but root polish failed (e.g. `pr > 1.0` or `not sol.success`), line 345 is skipped, leaving `converged = True`. The unpolished endpoint is emitted as a certified `x_eq`, and `message` falsely claims it was root-polished!
- **Replacement Implementation**:
  ```python
  if polish:
      if net.use_sparse:
          from scipy.optimize import least_squares
          sol_ls = least_squares(
              lambda y: rhs_cstr(net, np.maximum(y, 0.0), dilution, eps),
              x,
              jac=lambda y: net.jacobian(np.maximum(y, 0.0), dilution, eps),
              method="trf", tr_solver="lsmr",
              xtol=1e-11, ftol=1e-11, gtol=1e-11, max_nfev=300)
          xp, polish_ok = np.maximum(sol_ls.x, 0.0), bool(sol_ls.success)
      else:
          sol = scipy_root(
              lambda y: rhs_cstr(net, np.maximum(y, 0.0), dilution, eps),
              x, jac=lambda y: net.jacobian(np.maximum(y, 0.0), dilution, eps),
              method="hybr")
          xp, polish_ok = np.maximum(sol.x, 0.0), bool(sol.success)
      pr = _scaled_residual(rhs_cstr(net, xp, dilution, eps),
                            xp, conv_atol, conv_rtol)
      polished_successfully = bool(polish_ok and np.all(pr <= 1.0)
                                   and np.all(np.isfinite(xp)))
      if polished_successfully:
          x, resid_end, converged = xp, pr, True
      else:
          # If polish was requested but failed, revoke trajectory convergence
          # to prevent certifying unpolished endpoints
          converged = False

  if converged:
      x_eq = x
      ttype = "STATIONARY"
      message = ("converged (componentwise residual + window"
                 + (", root-polished" if (polish and polished_successfully) else "") + ")")
  ```
- **Verification Protocol**:
  Mock a scenario where trajectory window check passes but root polish diverges (`pr > 1.0`). Verify that `converged` evaluates to `False`, `ttype` evaluates to `UNRESOLVED`, and `x_eq` is `None`.

---

### 6.10 Fix 10: Thermal Hydrolysis Direction-Specific Gating — **IMPLEMENTED** (via `synth_dir`, not bare `synth=False`)
- **Target File**: `src/chemtrans/networks/generators.py`
- **Function**: `polymer_network` (lines 240–255)
- **Current Behavior**:
  ```python
  channels.append(Channel(
      reactants=((iPA, 1),),
      products=tuple(sorted(prods.items())),
      catalysts=(),
      k_f=k_hyd,
      k_r=k_hyd * np.exp(dG3),
      synth=True,
      mechanism="thermal",
  ))
  ```
- **The Defect**:
  `synth=True` scales thermal hydrolysis by $\varepsilon$. At $\varepsilon = 0$, spontaneous cleavage of polymers is 100% frozen ($k_f = 0, k_r = 0$), directly violating Invariant 3.
- **Replacement Implementation**:
  ```python
  channels.append(Channel(
      reactants=((iPA, 1),),
      products=tuple(sorted(prods.items())),
      catalysts=(),
      k_f=k_hyd,
      k_r=k_hyd * np.exp(dG3),
      synth=False,  # decay channels run unconditionally; Invariant 3
      mechanism="thermal",
  ))
  ```
  In `ReactionNetwork.rates()`:
  Verify that when `synth=False`, `bg = 1.0` unconditionally, allowing uncatalyzed thermal hydrolysis to operate at $\varepsilon = 0$.
- **Verification Protocol**:
  Instantiate `polymer_network(..., eps=0.0)`. Verify that for thermal hydrolysis channels, forward flux ($PA \to P + A$) is strictly positive for $x_{PA} > 0$.

---

### 6.11 Fix 11: Jump-Discontinuity Verification in `fold` Bifurcation Classification — **IMPLEMENTED**
- **Target File**: `src/chemtrans/analysis/critical_points.py` & `src/chemtrans/analysis/cocriticality.py`
- **Current Behavior**:
  `fold` parameter bounds allow $c_1 \in [-\infty, \infty]$. In `_LC_MAP`, any fit selecting `fold` is declared a `saddle_node` bifurcation.
- **The Defect**:
  When $c_1 \approx c_0$, the jump vanishes and the model is a continuous $\beta = 1/2$ second-order transition. When $c_1 < c_0$, the model jumps downwards. Both are mislabeled as first-order saddle-node folds.
- **Replacement Implementation**:
  In `critical_points.py` `_bounds`:
  Enforce jump positivity: require $c_1 \ge c_0 + \delta_{\text{jump}}$ where $\delta_{\text{jump}} = 0.05 \times \text{span}(y)$ during fold fitting, OR in `cocriticality.py`:
  ```python
  if sel == "fold":
      jump = float(params.get("c1", 0.0) - params.get("c0", 0.0))
      span = float(y.max() - y.min()) if len(y) else 1.0
      if jump >= 0.1 * span:
          cls = "saddle_node"
      else:
          cls = "continuous_half_power"  # beta = 1/2 continuous onset
      return LcEstimate(lc, cls, detail)
  ```
- **Verification Protocol**:
  Fit synthetic continuous data $y = \sqrt{\max(0, \lambda - 0.5)}$. Verify that if `fold` is selected, `estimate_lc_continuous` classifies it as `continuous_half_power` rather than `saddle_node`.

---

### 6.12 Fix 12: Reversibility Enforcement in Detailed Balance Verification — **IMPLEMENTED** (with declared `driven` flag)
- **Target File**: `src/chemtrans/networks/core.py`
- **Function**: `check_detailed_balance` (lines 480–490)
- **Current Behavior**:
  ```python
  if ch.k_f > 0 and ch.k_r > 0:
      ratio = ch.k_f / ch.k_r
      rel = abs(ratio / K - 1.0)
      if rel > tol:
          ok = False
  ```
- **The Defect**:
  If $k_r = 0$, the check is silently bypassed, returning `ok = True` on irreversible channels.
- **Replacement Implementation**:
  ```python
  for r, ch in enumerate(self.channels):
      nu = np.zeros(self.n_species)
      for i, a in ch.products:
          nu[i] += a
      for i, a in ch.reactants:
          nu[i] -= a
      dG = float(g @ nu)
      K = np.exp(-dG)
      if ch.k_f <= 0 or ch.k_r <= 0:
          # Irreversible channel cannot satisfy detailed balance unless dG is +/- infinity
          rows.append({"channel": r, "dG": dG, "K_eq": K, "kf": ch.k_f, "kr": ch.k_r,
                       "error": "irreversible_channel", "mechanism": ch.mechanism})
          ok = False
          continue
      ratio = ch.k_f / ch.k_r
      rel = abs(ratio / K - 1.0)
      rows.append(
          {"channel": r, "dG": dG, "K_eq": K, "kf/kr": ratio,
           "rel_err": rel, "mechanism": ch.mechanism}
      )
      if rel > tol:
          ok = False
  ```
- **Verification Protocol**:
  Construct a test channel with $k_f = 1.0$ and $k_r = 0.0$ on finite free energies. Verify `net.check_detailed_balance()["ok"]` evaluates to `False`.

---

### 6.13 Fix 13: Single-Valued Reversible Continuation Guard & True Hysteresis Detection — **IMPLEMENTED**
- **Target File**: `src/chemtrans/measure.py`
- **Function**: `_continuation_block` (lines 512–533)
- **Current Behavior**:
  Computes `hysteresis_width = abs(tf["lam"] - tr["lam"])` purely from threshold crossing indices without checking if branch states differ.
- **Replacement Implementation**:
  ```python
  # Align fwd and rev records by parameter value:
  fwd_by_lam = {r["lam"]: r for r in branches["fwd"]}
  rev_by_lam = {r["lam"]: r for r in branches["rev"]}
  common_lams = sorted(set(fwd_by_lam) & set(rev_by_lam))

  # Check for genuine state divergence across common points:
  state_diffs = {
      lam: abs(fwd_by_lam[lam]["m_poly"] - rev_by_lam[lam]["m_poly"])
      for lam in common_lams
  }
  max_branch_diff = max(state_diffs.values()) if state_diffs else 0.0
  bistable_lams = [lam for lam, diff in state_diffs.items() if diff > 2.0 * floor]

  if len(bistable_lams) >= 1 and tf and tr:
      out["hysteresis_width"] = abs(tf["lam"] - tr["lam"])
      out["spinodal_fwd"], out["spinodal_rev"] = tf["lam"], tr["lam"]
      out["bistable_window"] = [min(bistable_lams), max(bistable_lams)]
      out["max_branch_diff"] = max_branch_diff
  else:
      out["hysteresis_width"] = None
      out["spinodal_fwd"] = out["spinodal_rev"] = None
      out["max_branch_diff"] = max_branch_diff
      out["note"] = "no genuine branch bifurcation detected (single-valued response)"
  ```
- **Verification Protocol**:
  Run `_continuation_block` on `transcritical_motif` across its critical dilution. Verify that `hysteresis_width` evaluates to `None` and `note` indicates single-valued response.

---

### 6.14 Fix 14: Graceful Solver-Failure Recovery in Serial Dilution — **IMPLEMENTED**
- **Target File**: `src/chemtrans/dynamics/ode.py` & `src/chemtrans/measure.py`
- **Current Behavior**:
  `sd_map` raises `RuntimeError`, crashing the worker process on solver failure.
- **Replacement Implementation**:
  In `ode.py` `sd_map`:
  Catch solver failure and return `(starts[:c], ends[:c], False, sol.message)`.
  In `measure.py`:
  Check success flag; if False, emit strict-JSON-compliant record with `det["status"] = "UNRESOLVED"` and `stationary = False`, consistent with CSTR solver failure handling.
- **Verification Protocol**:
  Pass unphysical parameters causing `solve_ivp` to fail inside `sd_map`. Verify `measure_run` returns a valid strict-JSON record with `status: "UNRESOLVED"` without crashing.

---

## 7. Next Implementation Targets (Open High/Critical Issues, No Specs Yet)

The previous targets — SSA hazard parity, `immigration_death` thermodynamics, bidirectional `kinetic_accessible`, vector-field clipping, `min_pairs` floor, secondary-flattening, and the competition tie-break — are **implemented in Round 4** (ledger items 35–44). Remaining open targets:

1. **§2.13 — bound symmetry (MEDIUM)**: either extend `hinge`/`fold`/`itrans` to negative amplitudes (mirrored model set) or verify the control-axis sign is declared on every decreasing-axis protocol and fail loudly when it is not.
2. **§2.5 — interior IC coverage for `multistability_scan` (HIGH)**: supplement the boundary-only bank (food_inlet + ≤1e-3 nonfood) with interior sampling — random non-food mass allocations, mid-basin seeds from prior equilibria — and report `n_attractors` explicitly as a lower bound with IC-coverage caveats.
3. **§2.1/§3.2 — certified-point dropout + UQ on λc (HIGH)**: propagate near-critical `UNRESOLVED` dropouts into λc confidence intervals; add bootstrap CIs for fitted critical points, dwell fractions, and selection probabilities.
4. **§3.5 — equivalence margin (HIGH)**: bound `tol_pair` to a physically motivated fraction of the swept domain (e.g. ≤10%) rather than 2× grid spacing, or report the equivalence margin alongside every verdict.
5. **Campaign regeneration (process, not code)**: every committed `data/` artifact predates the unclipped vector field, `hazard_mode`, `min_pairs`, competition relabeling, and manifest-threshold recording — all results citing them must be rerun before the review's resolved code items translate into trustworthy numbers.
