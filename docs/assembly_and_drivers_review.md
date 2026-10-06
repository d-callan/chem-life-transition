# Review: Control Axes, Assembly Measurement, and the Nonequilibrium
# Physics of the Life Vector

Status: **tenth revision**. The assembly threshold is **derived, not
imported**: $A_{\text{cut}} \approx L^* - 2$ where $L^*$ is the length at
which the enumerable abiotic null's expected copy number of any specific
string falls below 1 (computable, V-dependent, falsifiable). Compute
requirements are now a **three-tier ladder** — laptop ($L \le 8$),
workstation+ ($L \sim 10\text{--}12$, sparse solvers), cluster
(contingency only: literal MA 15, million-trajectory scaling, large-N
confirmation). Corrects the alphabet misconception (length, not alphabet,
sets the MA ceiling — binary is the cheapest route) and reframes fragment
ligation (Knob E) as the mechanism that makes MA dynamically meaningful
under monomer-addition-only kinetics.
1. **Emergence vs. Maintenance**: $\varepsilon$ acts as an external symmetry-breaking
   field. $\varepsilon = 0$ measures maintenance/extinction of pre-existing states;
   $\varepsilon > 0$ measures de novo ignition/emergence from sterile prebiotic food.
2. **Causal Gating vs. Co-Criticality**: Nested causal pipeline
   ($\text{Topology} \to \text{Maintenance} \to \text{Heredity} \to \text{Selection}$,
   $\text{Maintenance} \to \text{Assembly}$). Gated collinearity must be distinguished
   from genuine co-criticality.
3. **Within-Network Pairing**: `_lc_for_obs` must not pool heterogeneous networks into
   a single regression scatter. Critical points must be fit per network family to avoid
   smearing sharp transitions into artificial smooth crossovers.
4. **Laptop vs. Cluster Compute Boundaries**: Testing literal $MA \ge 15$ in linear binary
   strings requires $L \ge 16$ ($> 130,000$ species; HPC cluster scope). On a laptop
   ($L \le 8$), the assay evaluates whether mass concentrates in the pool's top
   complexity class with an absolute mass floor.
5. **Structural & Hereditary Connectivity**: Strict `min_cat_len` collapses the graph
   unless accompanied by catalytic multiplicity ($\mu_{\text{cat}} > 1$). Heredity in
   longer polymers ($L \ge 4$) collapses into a single monostable web unless protected
   by sequence-selective catalysis ($p_{\text{info}}$).
6. **Hysteresis & First-Order Protocols**: Distinguishes continuous bifurcations from
   first-order fold bifurcations ($p=2$), requiring forward/reverse sweeps to measure
   hysteresis loops and nucleation barriers.
7. **Experimental Rigor & Methodology**: Enforces equivalence testing (TOST) against
   false co-criticality, exact paired reporting for $N < 10$, two-stage scout/zoom sweeps,
   lineage-biased ICs, and shared limiting resources for selection.
8. **Measurement-Level Safeguards**: Fixes six verified estimator biases —
   padded-trajectory dwell inflation, `O_persist` fraction-only threshold,
   undefined extinction at ε>0, `O_hered` blindness to amplitude bistability,
   residual-floor "equilibria" near criticality, and volume-rounding loss
   of seeded species.

---

## 1. Purpose, Observable Taxonomy, and the Causal Pipeline

This repository is a **controlled hypothesis-testing laboratory**: can the
chemistry $\to$ biochemistry boundary be modeled as a nonequilibrium phase transition —
one master transition, a cascade of bifurcations, orthogonal axes, or a conditional hierarchy?

### Emergence vs. Maintenance: The Role of $\varepsilon$

A fundamental distinction must govern all experimental designs:

- **The Absorbing Limit ($\varepsilon = 0$) — Maintenance & Extinction**:
  When uncatalyzed background ligation is zero ($\varepsilon = 0$), the polymer-free
  state ($x_{\text{poly}} \equiv 0$) is a strict absorbing boundary. Sterile food
  cannot synthesize any polymer. Initializing from dense seeds (`x0_policy="dense"`)
  tests **ecological survival / maintenance of an already existing system**. It measures
  washout and extinction, not de novo emergence.
- **The Driven Prebiotic Regime ($\varepsilon > 0$) — Spontaneous Ignition**:
  Uncatalyzed background synthesis $\varepsilon$ plays the exact role of an
  **external symmetry-breaking field $h$** (analogous to a magnetic field in the Ising model).
  At small $\varepsilon > 0$, thermal chemistry nucleates trace oligomers. If the drive
  ($\Delta\mu, \chi, \kappa$) exceeds a critical ignition threshold, these trace polymers
  autocatalytically explode into macroscopic organization.

| Parameter Regime | Initial Condition | Physical Question | Transition Type |
|---|---|---|---|
| $\varepsilon = 0$ | Dense polymer seed | Does living organization maintain itself against dilution? | Sharp absorbing-state extinction ($\text{DP}$) |
| $\varepsilon > 0$ | Sterile monomer food | Does living organization spontaneously ignite from non-life? | Imperfect bifurcation / ignition crossover |

### The Causal Hierarchy (The Gated Collinearity Trap)

The user’s central question is:
> *"If independently defined properties all become discontinuous/critical at the same $\lambda_c$, then you have something much closer to evidence for a genuine underlying chemistry $\to$ life phase transition."*

However, in an open reactor at $\varepsilon = 0$, these observables form a strict, nested causal DAG:

```
Network Topology (O_raf)
       │
       ▼
Self-Maintenance (O_maint, O_persist)
       ├───► Heredity (O_hered) ───► Selection (O_select)
       │
       └───► Assembly Complexity (O_assembly)
```

- **$\mathcal{O}_{\text{maint}}$ vs. $\mathcal{O}_{\text{persist}}$**: These are not two distinct biological properties. $\mathcal{O}_{\text{maint}}$ is the deterministic ($V \to \infty$) limit, while $\mathcal{O}_{\text{persist}}$ is the finite-$V$ stochastic dwell of the exact same self-maintenance state.
- **$\mathcal{O}_{\text{hered}}$ requires $\mathcal{O}_{\text{maint}}$**: At $\varepsilon = 0$, the dead state is unique. To have multiple attractors ($\mathcal{O}_{\text{hered}} \ge 1$), at least one attractor *must* bear nonzero polymer mass. Heredity cannot precede maintenance.
- **$\mathcal{O}_{\text{select}}$ requires $\mathcal{O}_{\text{hered}}$**: Competition bias is only defined between two distinct attractors.
- **$\mathcal{O}_{\text{assembly}}$ requires $\mathcal{O}_{\text{maint}}$**: High-assembly polymers cannot be maintained if total polymer mass is zero.

**The Null Trap:** If a single master bifurcation ignites polymer maintenance at $\lambda_c$, all downstream observables are unblocked simultaneously. Finding that $\lambda_c^{\text{maint}} = \lambda_c^{\text{hered}} = \lambda_c^{\text{assembly}}$ could simply reflect a **master gate opening** (turning the main valve feeds all downstream pipes), rather than independent biological properties converging co-critically.

**The Test:** To demonstrate genuine co-criticality, the framework must test whether heredity bifurcates *away* from the primary maintenance branch via a subsequent secondary bifurcation ($\lambda_c^{\text{maint}} \ne \lambda_c^{\text{hered}}$), or whether it nucleates simultaneously from the fold.

**The Stronger Test — Perturbational Co-Criticality:** Onset coincidence on a
single axis is weak evidence — a gating cascade can produce exact
coincidence. The discriminating claim is **surface-level**: if maintenance
and heredity share a genuine critical structure, then
$\Delta\lambda_c(\kappa) = \lambda_c^{\text{maint}}(\kappa) - \lambda_c^{\text{hered}}(\kappa) \approx 0$
should hold *as a second control axis varies* — across a $(D \times \kappa)$
or $(D \times \Delta\mu)$ scan, not just at one point. If the onsets coincide
on one axis but split under variation of a conditioning axis, the
coincidence was accidental, not co-critical.

### Observable Taxonomy

| Class | Observables | Role in Analysis |
|---|---|---|
| **Class A: Structural Topology** | `O_raf` | **Ensemble graph invariant**. Static on any fixed network. Invariant to $D, \Delta\mu, \chi, \kappa, p, V$. Evaluated only during generator sweeps (e.g. $p_{\text{cat}}$ percolation). **Removed from fixed-network co-criticality matrix**. |
| **Class B: Dynamical State** | `O_maint`<br>`O_hered`<br>`O_assembly` | **Primary co-criticality triad**. Fit with continuous/discrete onset models along physical driver sweeps ($\Delta\mu, \kappa, D$). |
| **Class B-Scale: Finite-Size** | `O_persist` | Finite-volume stochastic realization of `O_maint`. Evaluated across $V \in [16, 1024]$ for finite-size scaling; not pooled as an independent biological axis. |
| **Class C: Conditional State** | `O_select` | **Evaluated only where $O_{\text{hered}} > 0$**. Evaluates lineage asymmetry on the multistable island. Dropped from global $\lambda_c$ curve fitting. |

### Basin-Conditioned Observables (resolved)

*Was a verified defect — fixed.* `multistability.per_attractor_obs`
now evaluates `m_poly` and the absolute assembly signature on EVERY
certified attractor center (from `attractor_x0`), and `life_vector`
exposes it under `O_maint["basin"]` / `O_assembly["basin"]` with
declared reductions `max` and `basin_weighted_mean` (weights = IC
count per cluster — a basin-volume proxy). `value` still reads the
policy-seeded endpoint for backward compatibility. First use already
matters: on a 3-attractor net only the dominant-basin attractor
carried signature mass (0.025 vs 0.0) — the basin block makes that
dependence explicit instead of hidden.

*Original defect, for the record:*

`O_maint` and `O_assembly` are computed from the `equilibrium` block of a
**single** `integrate_cstr` run under one `x0_policy` — while `O_hered`
comes from a *different* code path that scans 16 ICs. In the multistable
regime, the reported `O_maint`/`O_assembly` therefore depend on **which
attractor the policy seed landed in** — an arbitrary basin sample, not a
representative value. Two identical networks can report different
`O_maint` purely from seeding, and the fitted $\lambda_c^{\text{maint}}$
can shift where the seed crosses a basin boundary mid-sweep.

Requirement: observables are **per-attractor**, not per-run. `O_maint`
reports `{attractor → m_poly}` over the certified stable equilibria from
the IC-bank scan, with a declared scalar reduction (e.g. `max` organized
mass); `O_assembly` is evaluated per attractor. The multistability block —
not the primary integration — is the source of truth in the multistable
regime (which is only known post hoc).

---

## 2. Analysis Architecture: Within-Network Pairing vs. Pooling

In the current `src/chemtrans/analysis/cocriticality.py`, `_lc_for_obs` concatenates all data points across all networks into a single regression array:
```python
lam = np.concatenate([s[0] for s in series.values()])
y = np.concatenate([s[1] for s in series.values()])
return estimate_lc_continuous(lam, y, seed=seed)
```

### The Multi-Network Pooling Defect
- Each random network realization has its own distinct topology, and therefore its own critical threshold: $\lambda_c(\text{net}_1) \ne \lambda_c(\text{net}_2)$.
- If $\text{net}_1$ has $\lambda_c^{\text{maint}} = 0.8$ and $\lambda_c^{\text{hered}} = 0.8$, while $\text{net}_2$ has $\lambda_c^{\text{maint}} = 1.2$ and $\lambda_c^{\text{hered}} = 1.2$:
  - On each individual network, the two properties are **100% coincident** ($\Delta\lambda_c = 0$).
  - But pooling their points together into one scatter plot smears the response into a broad, noisy curve, artificially inflating the residual variance and causing AICc to favor a smooth crossover (`sigmoid` or `itrans`) over a sharp transition!
- Furthermore, discrete observables (`O_hered`) snap to the parameter grid, while continuous regression yields a floating-point value. Comparing grid-snapped values to continuous fits introduces an artificial discretization bias bounded by $\Delta\lambda$.

### The Solution: Per-Network-Family Fitting
Co-criticality must be evaluated **per network family first**:
1. For each network family $k$, extract paired critical points $\lambda_c^{\text{maint}}(k)$, $\lambda_c^{\text{hered}}(k)$, $\lambda_c^{\text{assembly}}(k)$.
2. Compute the paired difference distribution:
   $$\Delta\lambda_c(k) = \lambda_c^{(i)}(k) - \lambda_c^{(j)}(k)$$
3. Test whether the paired difference bounds bracket zero.
4. This requires a stable `family_id` (derived from topology, alphabet, and random seed) while retaining `net_hash` for exact run provenance.

---

## 3. Candidate Control Axes (Physical Drivers vs. Structural Generators)

We treat the control axis as an explicit experimental knob. Response direction and transition class are measured outcomes, not assumptions.

| Axis | Parameter | Physical Meaning | Axis Class |
|---|---|---|---|
| **Washout Clearance** | `reactor.dilution` ($D$) | Clearance AND replenishment: $\dot{x} = Sv + D(x_{\text{in}} - x)$. Increasing $D$ raises fuel influx while shortening residence time. | Continuous dynamical |
| **Free Energy Drive** | `generator.delta_mu` ($\Delta\mu$) | Standard-state free energy gap of activated food ($\Delta\mu = \mu_{A^*} - \mu_A$). Note: $\Delta\mu = 0$ is NOT equilibrium in an open CSTR. | Continuous dynamical |
| **Inlet Fuel Fraction** | `generator.chi` ($\chi$) | Fraction of inlet food that is activated ($[A^*]_{\text{in}} / [A]_{\text{tot}}$). Pure boundary disequilibrium. | Continuous dynamical |
| **Catalytic Acceleration** | `generator.kappa` ($\kappa$) | Kinetic rate acceleration of catalyzed channels. Symmetrical fwd/rev; equilibrium constants $K_{\text{eq}}$ are unchanged. | Continuous dynamical |
| **Kinetic Cooperativity** | `generator.cat_power` ($p$) | Catalytic rate exponent $(x_c)^p$. Discrete switch enabling multistable fold bifurcations. | Model class (categorical) |
| **Prebiotic Background** | `task.eps` ($\varepsilon$, `net.eps` fallback) | External symmetry-breaking field — floor of `(eps + κ·x_C)` on BOTH directions of every channel. $\varepsilon=0 \implies$ absorbing limit; $\varepsilon>0 \implies$ de novo emergence. | Dynamical drive |
| **System Volume** | `stochastic.volume` ($V$) | Copy-number scale. Finite-size scaling axis to probe the Kurtz thermodynamic limit. | Finite-size axis |
| **Network Connectivity** | `p_cat`, `min_cat_len`, `p_info` | Random graph parameters controlling hypergraph percolation and RAF closure. | Structural ensemble axes |

---

## 4. Assembly Complexity: Compute Tiers, Not Cut Corners

We do not cut corners, lower thresholds arbitrarily, or pretend $MA=5$ is $MA=15$. With a **derived** threshold (below), the cluster stops being a scientific prerequisite and becomes a **scaling contingency**: the principle is testable at any pool size; what large pools add is dynamic range above the abiotic detection floor.

### The Combinatorial Reality of Assembly Index 15

In Assembly Theory (Marshall et al. 2021, Sharma et al. 2023), $MA \ge 15$ is
the threshold where the probability of abiotic production of an object at
non-trace copy number falls below the probability that selection/copying
produced it — partly empirical (the ~$10^8$-copy floor is an observational
constraint), partly theoretical (combinatorial growth of assembly space vs.
undirected sampling).

### The Toy-Model Threshold Is Derivable, Not Imported

The string analogue can be computed exactly rather than borrowed. Under our
sequence-blind null (composition-blind energetics; $\kappa=0$ or eps-only
dynamics), every length-$L$ string is produced with equal probability:

$$P(\text{specific } s \mid L\text{-mer}) = 2^{-L}$$

Expected copy number of the most abundant specific string at length $L$
under the abiotic null:

$$\mathbb{E}[n_s] = V \cdot M_L^{\text{null}} \cdot 2^{-L},
\qquad M_L^{\text{null}} \sim c\,r^{L}$$

where $r < 1$ is the effective abiotic elongation yield per length step
(measurable from a $\kappa=0$/eps-only null run). The threshold is where
this falls below the detection floor:

$$L^* = \min_L \{ L : V \cdot c\, r^L \cdot 2^{-L} < 1 \}, \qquad
A_{\text{cut}} \approx L^* - 2$$

Properties of the derived threshold:
- **Computable, not fitted**: measure $(c, r)$ from the null, solve for
  $L^*$, record the derivation in the run manifest.
- **Volume-dependent** (physically sensible): larger $V$ lowers the
  detection floor, so abiotic chemistry can occasionally populate longer
  strings in bigger reactors — the toy analogue of Walker's ~$10^8$-copy
  floor. $A_{\text{cut}}$ must be re-derived whenever $V$ or the null
  changes.
- **Falsifiable**: if the null is wrong (e.g. $r$ is length- or
  composition-dependent), update the derivation and restate the cut. The
  rule is "say why it is what it is; change it when the reasoning fails."

- In organic molecules, $MA \ge 15$ is reached at molecular weights $\sim 300\text{--}500\,\text{Da}$ because of **3D graph complexity** (branching, heteroatoms, rings, functional groups).
- In **linear binary polymers** (our toy string model), the exact-exhaustive max
  assembly index (recomputed over all strings through L=10) is:

  | L | max MA | top-class size / pool |
  |---:|---:|---|
  | 6 | 5 | 10 / 64 |
  | 7 | 6 | 4 / 128 |
  | 8 | 6 | 80 / 256 |
  | 9 | 7 | 56 / 512 |
  | 10 | 8 | 16 / 1024 |

  Growth is roughly $L-1$ to $L-2$ with plateaus (L8 is a plateau at max=6)
  — **not** a clean linear law and definitely not "logarithmic." The
  $L = 16 \implies \max(MA) \approx 15$ figure is an **unverified
  extrapolation** (likely ~13–14 given the trend); verify L=12–14
  (16k strings, trivially enumerable) before quoting it.

- **Length, not alphabet, sets the MA ceiling**: $a(s) \le L - 1$
  regardless of alphabet size. Binary is the *cheapest* pool for reaching
  $MA \ge 15$ ($L \approx 17 \implies 2^{18}-2 \approx 262\text{k}$
  species); a 4-letter alphabet at $L=16$ spans $4^{16} \approx 4.3$
  billion species. Larger alphabets are for pool richness/topology, not
  for reaching the threshold.

- **MA's pathway space exceeds our kinetics (deepest caveat)**:
  $a(s)$ counts steps in the shortest *block-joining* pathway
  ($P_1 + P_2 \to P_1P_2$), but the generator only allows monomer
  elongation ($P + A^* \to PA + W$). Under monomer-addition-only kinetics
  every length-$L$ string has identical kinetic distance $L-1$ — MA is a
  **compressibility measure decoupled from kinetic accessibility**, with
  the only real reuse being shared prefixes. Fragment ligation (Knob C)
  is what makes MA dynamically meaningful: with $P_1+P_2 \to P_1P_2$
  channels, MA literally counts the minimum catalytic steps to construct
  the object.
- For an exhaustive binary polymer network up to $L = 16$:
  $$\text{Species Count} = 2^{17} - 2 = \mathbf{131,070\text{ species}} \quad (> 260,000\text{ reactions}).$$
- Integrating a stiff, nonlinear chemical ODE with 130,000 species and dense Jacobians requires hundreds of gigabytes of RAM and massive parallel compute. It is fundamentally impossible on a local machine.

### The Estimator Also Can't Reach 15 — Not Just the Solver

`assembly_index` returns an explicit **upper bound** (estimator `"bound"`,
not the index) for strings longer than `exact_limit = 12`, and the
signature logic requires `est == "exact"` to qualify. So even if the
solver problem were solved, **every L>12 species currently fails the
signature test by construction** — the observable structurally cannot
detect $MA \ge 15$ in its current form. The cluster phase must include an
estimator upgrade (the concat-DAG search is exponential in length) or a
bound-aware criterion ($\text{bound} \ge A_{\text{cut}}$ with stated
caveats), plus substring caching across the pool.

### Three-Tier Scope Ladder

```
┌────────────────────────────────────────────────────────────────────────┐
│                        TIER 1: LAPTOP SCOPE                            │
│                           (max_len <= 8, N <= 510)                     │
│                                                                        │
│ • Full causal DAG co-criticality across physical axes (D, dmu, kappa) │
│ • Kurtz limit, multistability, Jacobian spectra, exact SSA             │
│ • Assembly metric: derived A_cut (top complexity class, e.g. a=5/L6)  │
│   with an ABSOLUTE mass concentration floor                            │
│ • All mechanism tests: ligation, min_cat_len, gamma_ma, p_info         │
│ • Emergence vs maintenance, hysteresis, gated vs co-critical           │
│ • Causal question: Can mass escape the dimer trap without cheats?      │
└────────────────────────────────────┬───────────────────────────────────┘
                                     │
                                     ▼
┌────────────────────────────────────────────────────────────────────────┐
│                     TIER 2: WORKSTATION+ SCOPE                         │
│                   (L ~ 10-12, N ~ 2k-8k species)                       │
│                                                                        │
│ • Sparse Jacobians (dense 8k x 8k ~ 512MB -> sparse solvers)           │
│ • Ligation channels (~10-100k) feasible                                │
│ • Several MA classes between trivial and maximal -> wider dynamic      │
│   range above the derived abiotic floor                                │
│ • Moderate SSA ensembles and V-scaling                                 │
│   CVODE / SUNDIALS with sparse LU — still a single machine             │
└────────────────────────────────────┬───────────────────────────────────┘
                                     │
                                     ▼
┌────────────────────────────────────────────────────────────────────────┐
│                 TIER 3: HPC CLUSTER (CONTINGENCY ONLY)                 │
│                      (max_len >= 16, N >= 130,000)                     │
│                                                                        │
│ • Literal MA >= 15 stress test — interpretive value only, the          │
│   derived-cut version is already tested in Tiers 1-2                   │
│ • Sparse stiff ODE solvers + rule-based/on-the-fly generation          │
│   (BioNetGen style): only populated species instantiated               │
│ • Million-trajectory extinction / finite-size scaling                  │
│ • N >= 100 confirmatory sweeps — THROUGHPUT, not feasibility           │
│   If Tier-1 results are flat, a bigger pool will not rescue the        │
│   physics; if they show structure, Tier 2 is the next escalation.      │
└────────────────────────────────────────────────────────────────────────┘
```

### The Corrected Laptop Observable (`O_assembly`)

1. **Abolish length-conditional quantiles**: The previous criterion ($a(s) > Q_{0.95}$) was mathematically impossible because random strings maximize Kolmogorov/assembly complexity ($Q_{0.95} = \max(MA)$). It also **samples the null per network family** (`seed_key + "|manull"`), so the qualifying threshold differs across families — injecting uncontrolled variance into cross-family `O_assembly` comparisons. An absolute cut removes both the unattainability and this variance source.
2. **Derived absolute complexity threshold**: In the laptop pool ($L \le 6$), define high assembly via the **derived cut** $A_{\text{cut}} \approx L^* - 2$ from the null computation above — at $L=6$, the pool's maximal class is $a=5$ (10 of 64 strings, e.g. `AAABBA`, `AABABB`), so $A_{\text{cut}} = 5$ is the laptop analogue, and its justification must be stated as "pool-maximal complexity class under the derived-null rule," not as an approximation of Walker 15.
3. **Absolute mass concentration floor**:
   $$\mathcal{O}_{\text{assembly}} = \sum_{s \in \text{High-MA}} m(s) x_s \quad \text{subject to } \mathcal{O}_{\text{assembly}} \ge M_{\text{floor}}$$
   Using relative mass fraction alone was flawed: if total polymer mass was $10^{-12}$, a relative fraction could still register $1.0$. The metric must measure **absolute concentration of complex material**.
4. **Report the full spectrum**: Output the complete $(MA, \text{mass})$ scatter and the chain-length distribution ($L \ge 4$ mass fraction) alongside the scalar order parameter.

---

## 5. Mechanism Knobs: Avoiding Graph Collapse & Loss of Heredity

### Stated model assumptions (context for all knobs)

- **Composition-blind energetics**: every elongation has identical
  $\Delta G = g_{\text{bond}} - \Delta\mu$; free energies depend on length
  only, never on sequence. The dimer trap is purely **kinetic/topological**,
  never thermodynamic. Any claim about why complex sequences lose must live
  in catalytic topology and stoichiometry.
- **A↔B permutation symmetry** holds at the generator level; symmetry is
  broken only by random catalyst assignment per network. The sign of
  $P(A \text{ wins})$ is meaningless across networks (A/B are arbitrary
  cluster indices set by IC order); $|P - 0.5|$ is the correct pooled
  statistic.


### Knob A: `min_cat_len` and Catalytic Multiplicity

- **The Leak (fixed)**: `p_selfrep` could assign a product shorter than
  `min_cat_len` as its own catalyst. Now self-rep requires
  `len(PA) >= min_cat_len`; short products fall through to the
  kin/random draw.
- **The Graph Collapse — corrected, stronger result**: `min_cat_len >= 3`
  is not merely unlikely to form RAFs — under elongation-only chemistry
  it is **structurally impossible**. RAF closure starts from food
  (monomers + activated monomers); the first producible polymer is a
  dimer, but the dimer channel's catalyst must have length $\ge 3$,
  which cannot be in the closure. Catch-22: 0/10 was an understatement
  — it is 0/$\infty$. Neither self-rep loops (product must already be
  in the closure) nor ligation (needs two polymers) can break it.
  **Consequence: `min_cat_len = 2` is the only viable catalyst-length
  floor in this chemistry; higher values are vacuous, not strict.**
- **Multiplicity (implemented)**: `n_cat_max` / `p_multi` add extra
  catalyst slots per catalyzed channel (same nested policy,
  deduplicated; `n_cat_max=1` reproduces the old ensemble
  bit-for-bit modulo the params record). Measured RAF frequency at
  `min_cat_len=2`, $p_{\text{cat}}=0.5$: 5/10 $\to$ 7/10 with
  `n_cat_max=3, p_multi=0.5` — a modest density boost, not a rescue;
  the deeper limit is the single-catalyst-per-channel × food-reachability
  structure, not slot count alone.

### Knob B: Sequence-Selective Catalysis ($p_{\text{info}}$) for Heredity

- Heredity requires multiple stable compositional attractors (e.g. Attractor A dominated by $A$-rich polymers, Attractor B dominated by $B$-rich polymers).
- In short polymers ($L \le 2$), there are only 4 dimers; mutual competition readily sustains multistability.
- In longer polymers ($L \ge 4$, 124 species), uniform random catalyst assignment causes widespread **catalytic cross-talk**: all species cross-catalyze each other into a single homogeneous web, collapsing multistability into a single monostable attractor ($\mathcal{O}_{\text{hered}} = 0$).
- **The Solution**: Maintaining heredity in large pools strictly requires **sequence-selective catalysis**:
  - Preferential catalysis of sequence kin ($p_{\text{info}}$ in `generators.py`), OR
  - Motif/complementarity-based catalytic assignment.
- **Exact `p_info` semantics (precisely weaker than the name suggests)**:
  the kin pool is `PA in s or s in PA` — substring/superstring kinship only,
  NOT complementarity or lineage identity. `AABB` is kin of `AB`, but so is
  `BBAA`. Whether substring kinship preserves A-vs-B lineage coherence is an
  open empirical question — do not assume it protects heredity until tested.

### Knob C: Complexity-Weighted Catalysis ($\gamma_{\text{ma}}$)

- Proposed law: $\kappa(s) = \kappa_0 [MA(s) / MA_{\max}(|s|)]^{\gamma_{\text{ma}}}$.
- **Warning**: Under relative scaling, dimers have $MA / MA_{\max} = 1.0$ (full weight), while complex hexamers have $5 / 5 = 1.0$ and intermediate hexamers have $4 / 5 = 0.8$. This law **does not disfavor dimers**; it only disfavors simple hexamers relative to complex hexamers.
- **Unbounded Alternative**: Scale by absolute assembly index: $\kappa(s) = \kappa_0 [MA(s)]^{\gamma}$.
- **Scientific Guardrail**: This is an *imposed intervention*, not spontaneous emergence.

### Knob D: Structural Protection ($\alpha_{\text{prot}}$)

- Must modify forward and reverse thermal rates symmetrically to preserve detailed balance ($\Delta G$ and $K_{\text{eq}}$ unchanged).
- **Caveat**: In the absorbing limit ($\varepsilon = 0$), thermal channels carry zero flux. Protection implemented on thermal channels has zero effect in the absorbing regime; it only functions when $\varepsilon > 0$.

### Knob E: Fragment Ligation & Templating (the MA-meaningfulness knob)

Two distinct extensions that must not be conflated:

- **Fragment ligation** ($P_1 + P_2 \rightleftharpoons P_1P_2$, no $W$)
  — **implemented** (`ligation=True`, `p_cat_lig`, `k_join` in
  `generators.py`; `GENERATOR_VERSION = "polymer_v3"`):
  the mechanism that makes MA *dynamically meaningful*. With only monomer
  addition, all length-$L$ strings share kinetic distance $L-1$ and MA
  measures compressibility, not assembly difficulty. Ligation channels
  align MA with the real pathway space — and give long strings a kinetic
  shortcut that can break the dimer trap (block assembly vs. serial
  addition).

  As built:
  - Substrates restricted to len≥2 polymers (monomer incorporation stays
    activation-gated via $A^*$; the thermal channel already covers
    polymer+monomer condensation).
  - Ordered pairs deduplicated on (reactant multiset, product):
    $(p_1,p_2)$ and $(p_2,p_1)$ collapse iff $p_1p_2 = p_2p_1$
    (e.g. $AA+AAA$ vs $AAA+AA \to AAAAA$); distinct products from the
    same reactants stay separate channels.
  - Catalyst assignment reuses the nested self/kin/random lottery on a
    dedicated latent-draw block, with `p_cat_lig` (default: `p_cat`)
    and the same `n_cat_max`/`p_multi` multiplicity policy.
  - Verified: letter conservation, detailed balance, absorbing limit at
    $\varepsilon=0$, no $W$, no monomer substrates; ODE and SSA run on
    ligated nets; live attractors reach $L=4$ including `ABAB` via the
    block join.
  - Channel counts (deduplicated): 16 at $L=4$, 78 at $L=5$, 266 at
    $L=6$ — steeper than species growth, as predicted.
  - Residual caveat: MA is a *lower bound* on kinetic distance — realized
    difficulty = pathway length × catalytic enhancement per join. A
    low-MA string on an uncatalyzed join can still lose to a high-MA
    string on a catalyzed path.
  - Interacts with Knob A: a single catalyzed ligation ($AA+BB \to AABB$)
    replaces three sequential catalyzed elongations, so strict
    `min_cat_len` may survive when it couldn't under elongation-only.
    Ligation and length-competence are NOT orthogonal — test jointly.
  - Scaling: directed ligations $\sum_{L=2}^{L_{\max}} 2^L (L-1) \approx
    516$ at $L=6$ but $\sim 2\times10^6$ at $L=17$ — channel count, not
    species count, is what pushes ligation-enabled pools toward cluster
    scale. Catalyst assignment: dedicated latent-draw block with
    `p_cat_lig` override (implemented).
- **Templating** (sequence-directed ligation, mismatch error probabilities,
  complementary duplex thermodynamics): a fidelity mechanism layered on top
  of ligation. Underspecified until its stoichiometry is defined; remains
  deferred.

---

## 6. Experimental Rigor, Statistics, and Methodological Safeguards

To prevent false discoveries, artifacts, and unrepeatable results, the experimental protocols must adhere to five methodological safeguards:

### 1. Equivalence Testing (TOST) vs. The "Indistinguishable" Fallacy
- **The Pitfall**: In `coincidence_report`, checking whether the 95% bootstrap CI of $\Delta\lambda_c$ covers zero and reporting `"indistinguishable"` is a classical Type II error. In small samples ($N=3$), CIs are wide purely from low power; covering zero is absence of evidence, not evidence of co-criticality.
- **The Protocol**: Enforce Two One-Sided Tests (TOST) with an explicit equivalence margin $\pm \delta_{\text{tol}}$ (e.g. grid spacing $\Delta\lambda = 0.1$):
  - **`coincident`**: 90% CI is completely contained within $[-\delta_{\text{tol}}, +\delta_{\text{tol}}]$.
  - **`separated`**: 95% CI completely excludes zero.
  - **`inconclusive`**: CI covers zero but extends beyond $\pm\delta_{\text{tol}}$ (insufficient statistical power).

### 2. Small-Sample Integrity: Exact Paired Tables for $N < 10$
- **The Pitfall**: For $N=3$ networks, only $\binom{3+3-1}{3} = 10$ unique bootstrap multisets exist. Running 400 bootstrap iterations generates an illusion of continuous tail coverage.
- **The Protocol**: For exploratory sweeps ($N < 10$), report the **exact empirical table of paired differences per network family** ($\Delta\lambda_c(k)$ for each $k \in \{1, \dots, N\}$). Reserve bootstrap resampling and CI estimation for confirmatory campaigns with $N \ge 20$.

### 3. Transition-Zone Sampling: Two-Stage Sweeps (Scout vs. Zoom)
- **The Pitfall**: Sweeping 10 points across a broad range ($D \in [0.05, 1.6]$) places only 1 point inside the steep transition zone. A 4-parameter nonlinear regression (`fold`, `sigmoid`, `itrans`) on 10 points has a data-to-parameter ratio of $2.5:1$, making model discrimination unreliable.
- **The Protocol**:
  - **Stage 1 (Scout)**: Coarse 10-point sweep across broad parameter space to bracket the active boundary.
  - **Stage 2 (Zoom)**: Dense 15–20 point sweep focused tightly on the transition window ($\Delta\lambda_{\text{zoom}} \approx 0.02$) to place at least 5 data points across the gradient, enabling reliable AICc model selection.

### 4. Lineage-Biased Initial Conditions in Multistability Scans
- **The Pitfall**: Starting multi-IC ODE scans solely from an unbiased dense mixture ($x_i = 0.15$ for all species) starts directly on the symmetry plane between lineage attractors, while random sparse seeds of scale $10^{-3}$ can fall below the nucleation threshold.
- **The Protocol**: The initial condition bank (`_ic_bank`) must include **lineage-biased seeds**:
  - $A$-biased seed: high concentration in $A$-homopolymers/rich oligomers, zero in $B$.
  - $B$-biased seed: high concentration in $B$-homopolymers/rich oligomers, zero in $A$.
  - This ensures that asymmetric composomes separated by a central separatrix are reliably discovered.

### 5. Shared Limiting Resources for Darwinian Selection ($\mathcal{O}_{\text{select}}$)
- **The Pitfall**: When $A^*$ and $B^*$ monomers are supplied via independent inlets, pure $A$-rich and $B$-rich lineages do not compete for food. They are uncoupled parallel systems; competition SSA measures demographic drift or kinetic asymmetry rather than true resource exclusion.
- **The Protocol**: To measure genuine Darwinian selection in $\mathcal{O}_{\text{select}}$, lineages must compete for a **shared limiting constraint** (e.g. shared non-specific energy activator $E^*$, shared catalytic binding sites, or mixed sequences competing for common ligation fragments).

### 6. Measurement-Level Safeguards (verified estimator biases)

Six defects verified directly in the implementation; each silently biases
results in exactly the regime of interest.

**6a. `O_persist` censoring bias (verified).**
`gillespie.simulate` pads truncated trajectories (runs stopped at
`max_events` before `t_max`) with the final state held constant across the
remaining record grid. `dwell_organized = above.mean()` then counts padded
rows as organized — systematically **inflating dwell for slow,
high-activity runs**, i.e. precisely the interesting regime. Requirement:
compute dwell only over `t ∈ [0, censor_time]`, report `n_truncated`, and
treat truncated runs as right-censored (or exclude with count stated).

**6b. `O_persist` fraction-floor bug (verified).**
`above = mass_raf > theta` uses a **fraction**, not absolute mass. A run
with near-zero polymer mass concentrated entirely in RAF species scores
`mass_raf ≈ 1` and registers "organized." This is the same relative-fraction
trap fixed for `O_assembly`. Requirement: joint criterion — absolute
nonfood mass above a floor **AND** `mass_raf > θ`.

**6c. Extinction is undefined at ε>0 (verified).**
`t_extinction` fires only when all nonfood counts hit zero. At ε>0 that
never happens — the observable's semantics silently flip between regimes.
Also: eps enters as the floor of the catalytic multiplier
`(eps + κ·x_C)` on **both directions of every channel** — it is a
symmetry-breaking field for formation AND decay, not a pure one-way
source. Requirement: at ε>0 redefine persistence as quasi-extinction
(time below an absolute mass floor) and describe dwell as a residency-time
measure. Also note the spec path is `task.eps` (with `net.eps` fallback),
not `reactor.eps`.

**6d. `O_hered` is blind to amplitude bistability (verified).**
`attractor_signature` is a **normalized** mass-fraction vector — amplitude
is discarded at `v / s`. A p=2 saddle-node produces high/low-mass branches
of the same composition → identical signatures → clustered as one
attractor → `O_hered = 0` **exactly where multistability is expected**.
Competition labeling inherits the same degeneracy (nearest-signature
within `4·tol`). Requirement: cluster on a joint feature
(composition + `log m`), or report `n_attractors_compositional` and
`n_attractors_total` separately; competition reports `P(dead)` and
`P(other)` alongside `P(A)`.

**6e. Residual floor admits slow-drift "equilibria" (verified).**
`STATIONARY` certification requires `|dx_i/dt| ≤ 1e-8 + 1e-5·|x_i|`. Near a
fold, critical slowing lets drift at |ẋ| ~ 1e-6 certify as equilibrium;
certified λ_c inherits bias bounded by this floor. Requirement: record the
scaled residual in each `equilibrium` block; run a doubled-`t_max`
sensitivity check at grid points adjacent to candidate critical points.

**6f. Integer rounding drops dilute seeded species (minor).**
`counts_from_concentration` uses `rint(x·V)` per species; at V=16 a seeded
species at x<0.031 rounds to zero and silently disappears from the dense
initial state. Requirement: report zeroed seeded species; require
`x_i·V ≥ 0.5` for species intended to be present.

**6g. Oscillatory regimes silently truncate the series (verified).**
`O_maint` and `O_assembly` are undefined when the deterministic status is
`OSCILLATORY`/`UNRESOLVED`/`RUNAWAY`. A Hopf onset in the swept range
drops every oscillatory point from the λ–O curve, so model selection can
read a Hopf transition as an extinction edge. Requirement: report the
fraction of sweep points per status class; a contiguous OSCILLATORY
region is itself a transition class (report it, don't drop it);
`osc_amplitude` is already recorded and can carry an `O_maint_osc`
extension.

**6h. Greedy attractor clustering is order-dependent (verified).**
`cluster_attractors` assigns each equilibrium to the first center within
`tol`; on a distance chain (A–B 0.04, B–C 0.04, A–C 0.08, tol=0.05) the
count depends on IC ordering, not physics. Requirement: use
agglomerative/complete-linkage clustering (order-invariant), or ship the
pairwise `center_distances` matrix so borderline merges are auditable.

**6i. Uncalibrated constants and headline choices (minor).**
- `theta_organized = 0.1` is a magic constant; justify it against the
  bimodal `mass_raf` distribution (organized states sit near ~0.7+),
  same requirement as `abundance_floor`.
- Competition headline currently reports the `raw` protocol
  (`out["p_a_wins"] = raw`), but `eq_mass` exists precisely because raw
  mixtures confound basin asymmetry with mass asymmetry. Report both
  side-by-side; make `eq_mass` the headline estimator.
- `itrans`'s `c_field` is a free imperfect-bifurcation handle: at fixed ε,
  fitted `c_field` should scale with ε — a cheap internal consistency
  check for emergence-regime fits.

**6j. Cross-axis normalization and protocol details.**
- `O_maint` is absolute nonfood mass; sweeps that change the mass budget
  (`f_inlet`, `chi`, and `D` via residence-flux product) change the
  achievable maximum. Report realized food influx `D·x_in` alongside
  `O_maint` for cross-axis amplitude comparisons.
- The `fold` model (`c1 + A·√(L−lf)`) describes a spinodal branch shape;
  model selection on single-direction data says nothing about the
  coexistence region. Branch-shape verdicts and hysteresis measurements
  answer different questions.
- `competition.volume` defaults to 16 independently of
  `stochastic.volume`; competition V must be reported, and matched to
  `O_persist`'s V where the two are compared.
- The serial-dilution mode (`sd_map`) exists but is unused; it is the
  natural protocol for propagation/inheritance questions (lineage
  survival under serial bottlenecks — a more literal life criterion than
  CSTR dwell) and should appear in the roadmap.

---

## 7. Statistical Mechanics: Bifurcations, Hysteresis, and Coexistence

1. **The Kurtz Thermodynamic Limit**: In a well-mixed CSTR (0D, no spatial gradients), Kurtz's theorem proves that as volume $V \to \infty$, the Gillespie jump Markov process converges to the deterministic ODE system on any finite time horizon $[0, T]$.
2. **Continuous Transitions (Transcritical, Pitchfork)**:
   - For continuous branch bifurcations (e.g. $p=1$ at $D_c \approx 1.13$), the deterministic bifurcation point **is** the true infinite-volume phase transition ($V \to \infty$).
   - Belong to the mean-field Directed Percolation universality class ($\beta = 1$).
3. **Discontinuous Transitions (Fold / Saddle-Node at $p=2$)**:
   - The deterministic fold point $D_{\text{fold}}$ is a **spinodal limit of metastability** (where $\text{Re}(\lambda_{\max}) = 0$).
   - In the true thermodynamic limit under noise, the first-order transition occurs at a **Maxwell coexistence point** $D^* \ne D_{\text{fold}}$, where the large-deviation quasi-potentials $\Phi(\text{active}) = \Phi(\text{absorbing})$ are equal.
4. **Hysteresis Protocols for First-Order Candidates**:
   - Single-point dense seeding only measures the upper extinction threshold $\lambda_{\text{extinct}}$.
   - To detect first-order transitions, the protocol must execute **forward and reverse continuation sweeps**:
     - *Forward*: start from sterile/low-activity state, increasing drive to find the ignition point $\lambda_{\text{ignite}}$.
     - *Reverse*: start from active state, decreasing drive to find the washout point $\lambda_{\text{extinct}}$.
     - Measure the hysteresis width $\Delta\lambda_{\text{hyst}} = |\lambda_{\text{ignite}} - \lambda_{\text{extinct}}|$.

---

## 8. Execution Roadmap (Structured by Compute Scope)

### Phase 1: Local / Laptop Engine (Foundation & Observable Fixes)
- [x] **Step 0: Observable Taxonomy & Assay Gates** *(implemented)*:
  - ~~Implement the **absolute $MA \ge 5$ + absolute mass floor**
    observable~~ → `assembly_signature_absolute(net, x, a_cut, mass_floor)`
    in `observables/assembly.py`; default `a_cut` = pool's maximal exact
    index, `mass_floor` is an absolute concentration. Full per-species
    `(a, estimator, mass_abs, mass_fraction)` spectrum emitted;
    `est="bound"` never qualifies.
  - ~~Positive-control validation~~ → `test_observables.py`:
    all-mass-in-max-MA-species signs; trace-mass fraction-1.0 species
    fails on the absolute floor; explicit `a_cut` respected.
  - ~~`O_raf` out of the co-criticality matrix~~ → `STRUCTURAL` class in
    `cocriticality.py`; excluded from `DEFAULT_OBSERVABLES`.
  - ~~`O_select` conditional~~ → `CONDITIONAL` class; never fit as a
    global λc curve (`lc_by_family` returns {} and pooled estimate
    reports `class="conditional"`).
  - ~~Measurement biases~~ → `measure.py`: `dwell_organized` masked to
    `t <= censor_time`; organized requires `mass_raf > θ` AND
    `nf_mass >= mass_floor`; `t_extinction` = quasi-extinction below the
    mass floor (valid at ε>0), `t_absorbing` = strict `nf==0` (ε=0 only);
    both masked post-organization; `scaled_resid_max` recorded on
    certified equilibria; `n_seed_dropped` reports species rounded to
    zero at volume V.
  - `heredity.cluster_attractors` gained `mass_ratio_tol` (joint
    composition+mass clustering); `multistability` block reports
    `n_attractors` (compositional), `n_attractors_total` (joint),
    `attractor_masses`, `n_kept`.
- [x] **Step 1: Network-Family ID, Within-Network Fitting & TOST** *(implemented)*:
  - `family_id = seed_key|net_hash` written by `campaign.py`;
    `group_by_network` keys on it (net_hash fallback). Generator sweeps
    change net_hash → distinct families, never paired.
  - `lc_by_family`: per-family fits (continuous → model selection;
    discrete → per-family changepoint). No pooled regression.
  - `coincidence_report` verdicts on paired differences
    $\Delta\lambda_c(f)$: `coincident_within_tol` (exact table,
    $N \le 12$ pairs), `coincident` (TOST: 90% CI inside ±`tol_pair`),
    `separated`, `mixed`, `inconclusive`, `incomparable`.
    `tol_pair` defaults to 2× median grid spacing; `grid_resolution`
    reported. Pooled estimates retained but flagged `"pooled": true`.
  - `control_axis` block: `{param, direction, label}` in the protocol
    drives `rec["control"]` (sign, axis, label); legacy
    `control_negate` still honored.
- [x] **Step 2: Catalytic Multiplicity & Lineage-Biased ICs** *(implemented)*:
  - ~~Self-rep length leak~~ → `pick_catalyst` requires
    `len(PA) >= min_cat_len` for self-rep; `min_cat_len >= 3` shown to be
    *structurally impossible* under elongation-only chemistry (see Knob A).
  - ~~Catalytic multiplicity~~ → `n_cat_max` / `p_multi` generator params;
    measured RAF frequency 5/10 → 7/10 at $(p_{\text{cat}}=0.5,
    \text{min\_cat\_len}=2)$ — modest boost, not a rescue.
  - ~~Lineage-biased ICs~~ → `_ic_bank` adds one seed per letter seeding
    species where that letter is $\ge 2/3$ of composition.
  - *Remaining (empirical)*: does $p_{\text{info}}$ preserve
    multistability in $L \ge 4$ pools? Needs a dedicated sweep — the
    biased ICs are now in place to detect it if present.
- [~] **Step 3: Small-Ensemble Driver & Emergence Scans ($N=3$)** *(in progress)*:
  - ~~Sweep $\Delta\mu$ at fixed $D$ under both regimes~~ *(done)*:
    `cocrit_dmu_maint` (ε=0, dense) + `cocrit_dmu_emerge` (ε=10⁻³,
    food-only), Δμ ∈ [0,8], N=3, ligation+multiplicity on.
    **First results**: maintenance onsets at Δμ ≈ 1.4–2.5; on net2
    $O_{\text{maint}}$ (2.51), $O_{\text{assembly}}$ (2.81) and the
    $O_{\text{hered}}$ 2→3-attractor jump (3.0) co-localize within
    `tol_pair` — consistent with gated collinearity, needs the 2-D test.
    net1 shows `always_active` heredity: genuine kinetic bistability at
    **Δμ = 0** (open CSTR, catalytic autocatalysis — not forbidden).
    Emergence: $m_{\text{poly}} \approx 0.002$ flat at all Δμ — no de
    novo organization at ε=10⁻³, consistent with the dilute null.
    Discrete estimator upgraded: changepoint = first *change* from
    baseline (catches 2→3 transitions), `always_active` reported;
    continuous fits flag `extrapolated` λc outside the grid.
  - ~~Forward/reverse continuation for hysteresis~~ *(implemented)*:
    `reactor.continuation = {"param": "dilution"|"eps", "values", ...}`
    warm-starts each point from the previous endpoint, rev from the fwd
    endpoint; records per-branch series + spinodal pair +
    `hysteresis_width`. **Verified**: ε=0 gives *infinite* width
    (absorbing — rev branch never reignites); ε=0.01 gives a finite
    loop on net2 (organized to D≈1.0 up, reignites at D≈0.5 down,
    width ≈0.6). Transitions located this way are spinodals, not the
    coexistence point.
  - ~~2-D perturbational scan ($D \times \kappa$)~~ *(done, first
    result)*: `coincidence_2d` partitions records by the secondary
    axis and aggregates per-level paired diffs — verdicts `tracking`
    (coincidence persists: shared critical locus) vs `drift` (gated
    collinearity at best) / `separated_at_all_levels` /
    `underdetermined`; `lc_curves` gives per-family λc(secondary)
    trajectories. Campaign `cocrit_2d_dk` (3 nets × 8 D × 3 κ, ε=0,
    deterministic only — SSA dropped for runtime): **all comparable
    pairs report `tracking`** — maint|assembly, maint|hered,
    assembly|hered coincide at κ=100 AND κ=500, consistent with a
    shared critical locus rather than gated collinearity. Caveats:
    only 1/3 families contributes pairs; κ=1500 degrades transitions
    to crossovers (higher catalytic rate is NOT monotonically
    favorable); O_persist unmeasured (needs the SSA campaign on the
    resolved transition line).

    **n=8 replication (decisive)**: the single-family `tracking` did
    NOT generalize. At N=8, `O_maint|O_assembly` and `O_maint|O_hered`
    report **`drift`**: coincidence at κ=100, mixed at κ=1500
    (one family +0.72), separated at κ=500 (one family −0.65).
    Reading: coincidence is *family-specific*, not universal — gated
    collinearity in some networks, genuine separation in others. The
    framework discriminated as designed. Also fixed: extrapolated λc
    fits are now excluded from paired diffs (one out-of-grid fit had
    poisoned a κ level). Remaining separations are ~1.5–1.8× grid
    spacing — the `zoom` subcommand exists to resolve exactly these.

    **n=8 zoom (D ~0.10 spacing over [0.14,1.05])**: the drift
    RESOLVED — all three comparable pairs `tracking` at a *tighter*
    tol (0.20 vs scout 0.40), every diff within ±0.1 λ at all κ
    levels. The coarse-grid separations were grid-snap artifacts (the
    known discrete-vs-continuous resolution bias) plus out-of-window
    fits. Caveats: the contributing family set shifts with κ (some
    λc leave the window or go undefined at κ=1500), so the claim is
    "coincidence wherever both observables are resolved", not "all 8
    families share one locus"; O_persist still unmeasured.

    **O_persist transition-line campaign** (`persist_transition_line`,
    κ=500, D∈{0.3..0.9}, 3 SSA reps × V=64 × t=200): O_persist joins
    the co-critical set — `O_maint|O_persist` −0.14,
    `O_persist|O_assembly` +0.13, `O_persist|O_hered` +0.33, all
    within tol; n_pairs=1 each (most families keep O_persist defined
    on or off throughout the window). Note: the pooled
    transition_fit selected `fold` on only n=4 points with
    rss≈1e-35 — a 4-param interpolation of 4 points, i.e. overfit,
    not evidence for a fold.
  - ~~Serial-dilution propagation protocol~~ *(implemented)*: the `sd`
    reactor now records per-cycle `m_poly_cycles`, `comp_dist_cycle`
    (cycle-to-cycle L1 signature drift — a stable lineage converges to
    ~0), `comp_dist_to_seed` (does the initial imprint survive?),
    `final_signature_top`, and Floquet stability of the periodic fixed
    point. `x0_policy="lineage:<letter>"` seeds a letter-dominated
    lineage. **First result**: on a monostable net, organized material
    survives 30 dilution cycles (m_poly≈0.23, stable fixed point) but
    lineage:a and lineage:b seeds converge to the *same* signature —
    **propagation without heredity**. The diagnostic is live; a
    positive heredity result requires a net multistable under the SD
    map.
  - *Remaining*: two-stage Scout→Zoom sweeps; κ axis sweeps.

### Phase 2: Workstation+ Engine (Tier 2 — dynamic range, not thresholds)
- [ ] **Step 4: Sparse-ODE upgrade**:
  - Sparse Jacobians for $L \sim 10\text{--}12$ ($N \sim 2\text{k--}8\text{k}$ species).
  - Ligation channels (~10–100k) + explicit catalyst-assignment rule for
    the new channel class (`u_cat` extension or `p_cat_lig`).
  - MA estimator upgrade or bound-aware criterion (L>12 returns bounds).
- [ ] **Step 5: Derived-threshold assembly tests at scale**:
  - Re-derive $A_{\text{cut}}$ at the larger pool and at the volumes used;
    test whether driven chemistry populates above-cut species.

### Phase 3: HPC Cluster (Tier 3 — contingency, not prerequisite)
- [ ] **Step 6: Confirmatory scaling** (only if Tier 1–2 warrant):
  - Literal $MA \ge 15$ stress test ($L \ge 16$, $N \ge 130\text{k}$
    species, ~2M ligation channels) — interpretive value only.
  - Million-trajectory extinction scaling ($V = 10^2 \to 10^6$) for
    critical exponents and quasi-potential coexistence.
  - $N \ge 100$ confirmatory sweeps — a throughput win, not feasibility.
