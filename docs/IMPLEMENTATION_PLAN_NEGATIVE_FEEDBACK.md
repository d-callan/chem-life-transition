# Implementation Specification: Negative Feedback & Length-Selective Constraints

> **STATUS: IMPLEMENTED** — see repo `runs/NOTEBOOK.md` and the pilot
> `neg_feedback_l6` for as-built results. Formula for Knob A was
> corrected during implementation (see §2.1 note).

**Target Branch/Repo**: `chemtrans 2.0` (`/home/dcallan-adm/chem-life-transition`)  
**Objective**: Introduce physically-grounded negative feedback and length-selective filter mechanisms into the model to break the "dimer sink" (where 2-mers consume all catalytic and monomer flux) and enable flux to reach longer, higher-assembly polymers ($L \ge 4$).

---

## 1. Problem Statement & Mathematical Diagnosis

In the pilot study (`catret_l6`), we demonstrated that:
1. **Positive feedback alone ($\kappa \gg 0$) accelerates the shortest loops first.** Because dimers ($L=2$) require only 1 bond, they ignite and deplete the monomer/fuel pool before higher-order polymers ($L \ge 4$) can form.
2. **Flat compartmental retention ($\rho$) amplifies whatever is already present.** It concentrated 2-mers and 3-mers by $15\times$, but left $L \ge 4$ at zero concentration because 2-mers are retained just as effectively as 6-mers.
3. **Smooth catalytic reweighting ($(L/2)^\alpha$) cannot bootstrap absent species.** If 6-mers are at concentration $0$, their $9\times$ potency multiplier yields $9.0 \times 0 = 0$ boost.

To route mass and catalytic activity into higher-assembly species, the system requires **negative constraints** that suppress, penalize, or filter out short trivial cycles.

---

## 2. Proposed Negative Feedback & Constraint Knobs

We specify three complementary, physically-grounded negative constraints. Each knob defaults to inactive/legacy behavior ($0.0$ or existing defaults) to maintain exact backward compatibility.

### 2.1 Knob A: Size-Dependent Dialysis / Length-Graded Retention (`retention_len_scale`, $\gamma \ge 0$)

#### Physical Rationale
In prebiotic environments (lipid vesicles with nanopores, mineral pore networks under thermal gradients, or coacervate droplets), small oligomers (2-mers, 3-mers) diffuse across membranes or partition weakly into droplets. Only large macromolecules ($L \ge 4\text{--}6$) are sterically trapped or phase-separated.

#### Mathematical Formulation
Currently, in `ReactionNetwork.jacobian` and `ode.rhs_cstr`, polymer outflow is governed by:
$$\text{outflow}_i = D \cdot (1 - \rho_i)$$
where $\rho_i = \text{polymer\_retention}$ for all species with $L_i \ge 2$.

Under length-graded retention with parameter $\gamma = \text{retention\_len\_scale} \in [0, 1]$:
$$\rho_i(L_i) = \rho \cdot \left(1 - \frac{2\gamma}{L_i}\right) \quad \text{for } L_i \ge 2$$
*(Note: an earlier draft of this section used $\rho(1 - (2/L)^\gamma)$, which degenerates to zero retention at $\gamma=0$ and breaks backward compatibility. The mixing form above is what was implemented — it interpolates between flat retention and the size cutoff.)*

- When $\gamma = 0$: $\rho_i = \rho$ for all polymers (exact legacy flat retention).
- When $\gamma = 1.0$:
  - $L=2$: $\rho_2 = 0$ (2-mers wash out at full dilution rate $D$).
  - $L=3$: $\rho_3 = \rho/3$.
  - $L=4$: $\rho_4 = \rho/2$.
  - $L=6$: $\rho_6 = 2\rho/3$.
- $\gamma = 0.5$: $\rho_i = \rho(1 - 1/L)$ — intermediate gradient.

**Effect**: 2-mers face the full open-CSTR wash-out rate $D$, preventing them from accumulating into an impenetrable sink, while 4-mers and 6-mers are selectively retained inside the compartment.

---

### 2.2 Knob B: Hard Catalyst Threshold Exposure (`min_cat_len`, $L_{\text{min}} \ge 2$)

#### Physical Rationale
In aqueous chemistry, 2-mers and 3-mers cannot fold into stable 3D catalytic pockets. Enzymatic competence requires a minimum polymer chain length.

#### Current Code & Necessary Fix
`src/chemtrans/networks/generators.py` already accepts `min_cat_len: int = 2`.
However, in campaign protocols and `measure.py`:
- `min_cat_len` is not swept or exposed in task defaults.
- When `min_cat_len >= 3`, short homopolymers (`AA`, `BB`, `AB`) are **strictly prohibited from acting as catalysts**.
- **Critical Architectural Requirement**: At $\varepsilon = 0$, if `min_cat_len >= 3` and the reactor is initialized with `x0_policy="food_only"`, uncatalyzed synthesis is zero and no species of length $\ge 3$ exists. The system would be trapped at dead absorbing state. Therefore:
  - If `min_cat_len >= 3`, the system requires either $\varepsilon > 0$ (spontaneous thermal background allows rare synthesis of 3-mers), OR a seeded initialization (`x0_policy="dense"` or lineage seed).

---

### 2.3 Knob C: Product Inhibition / Resource Saturation (`k_inhibit`, $K_{\text{inh}} \ge 0$)

#### Physical Rationale
In real enzyme kinetics (and prebiotic template replication), products bind to catalytic sites and inhibit further turnover unless displaced. In simple autocatalytic systems, when a single species (e.g. `AA`) reaches high concentration, it crowds the catalyst or clogs the active site.

#### Mathematical Formulation
The catalytic acceleration factor for channel $r$ with catalyst set $C_r$ becomes:
$$\text{boost}_r = \varepsilon + \frac{\kappa \sum_{c \in C_r} w(c) x_c^p}{1 + K_{\text{inh}} \sum_{c \in C_r} x_c}$$
- When $K_{\text{inh}} = 0.0$: standard uninhibited rate law (exact backward compatibility).
- When $K_{\text{inh}} > 0.0$: when catalyst $x_c$ reaches high concentration ($x_c \gg 1/K_{\text{inh}}$), its marginal acceleration saturates and drops, preventing runaway monomer monopolization by the fastest-replicating short species.

---

## 3. Concrete Code Changes Required

### 3.1 `src/chemtrans/networks/core.py`
1. Update `ReactionNetwork.__init__`:
   - Add `k_inhibit: float = 0.0` (validated $\ge 0$).
2. Update `_rate_vectors` / `rates`:
   - Compute denominator: `denom = 1.0 + self.k_inhibit * (C @ x)` if `self.k_inhibit > 0` else `1.0`.
   - Divide catalytic term by `denom`.
3. Update `d_rates` and `jacobian`:
   - Include quotient rule derivative for catalyst and participant species when `k_inhibit > 0`.
4. Update `to_dict` / `from_dict`:
   - Add `k_inhibit` to schema and content hash.

### 3.2 `src/chemtrans/dynamics/ode.py`
1. Update `rhs_cstr` and `integrate_cstr`:
   - Add `retention_len_scale: float = 0.0` parameter.
   - Compute per-species retention vector:
     $$\text{outflow}_i = D \cdot \left(1.0 - \rho \cdot \left[1.0 - \left(\frac{2}{\max(L_i, 2)}\right)^\gamma\right]\right)$$
     for $i \in \text{polymer\_idx}$, with $\rho_i = 0$ for food/fuel/waste.
   - Update `jacobian`: apply identical diagonal shift $- \text{outflow}_i$.
   - Update trajectory material-balance diagnostic integration to use the length-graded $\text{outflow}_i$.

### 3.3 `src/chemtrans/dynamics/gillespie.py`
1. Update `build_processes`:
   - Pass `retention_len_scale` into dilution rate calculation:
     $$\text{rate}_i = D \cdot (1 - \rho_i(L_i))$$
   - For `k_inhibit > 0`: In SSA "functional" mode, evaluate hazard with inhibition denominator.

### 3.4 `src/chemtrans/measure.py`
1. Unpack `retention_len_scale` from `reactor` dictionary (default `0.0`).
2. Pass `retention_len_scale` to `integrate_cstr`, `spectral.jacobian_eigs`, and SSA `build_processes` / `simulate`.

---

## 4. Verification & Testing Plan

1. **Regression Baseline**:
   - Ensure all 110 existing unit tests pass with defaults `k_inhibit=0.0` and `retention_len_scale=0.0`.
2. **Dedicated Unit Tests (`tests/test_negative_feedback.py`)**:
   - **Test 1: Length-graded retention selectivity**:
     Verify $L=2$ species have dilution rate $D$ (unretained), while $L=6$ species have dilution rate $D(1 - \frac{2}{3}\rho)$.
   - **Test 2: Jacobian & finite-difference consistency**:
     Verify analytic Jacobian matches centered numerical differentiation of `rhs_cstr` under both `retention_len_scale > 0` and `k_inhibit > 0`.
   - **Test 3: Detailed balance preservation**:
     Verify $k_f / k_r = \exp(-\Delta G)$ holds exactly with product inhibition.
   - **Test 4: Absorbing limit preservation**:
     Verify dead state remains absorbing at $\varepsilon = 0$.

---

## 5. Pilot Campaign Protocol: `neg_feedback_l6.json`

- **Pool**: $L_{\text{max}} = 6$, Alphabet $\{A, B\}$, $N=4$ network seeds.
- **Factorial Grid**:
  - `reactor.polymer_retention`: $[0.8, 0.95]$
  - `reactor.retention_len_scale`: $[0.0, 1.0, 2.0]$
  - `generator.min_cat_len`: $[2, 3]$
- **Success Criteria**:
  - Does mass in $L \ge 4$ turn strictly positive?
  - Does $O_{\text{assembly}}$ (`mass_signature_abs`) exceed the detection threshold?

---

## 6. As-Implemented Details & Pilot Results

### 6.1 Implementation notes (deviations from §2-§5)

- **Knob A formula corrected**: the draft's `rho_i = rho*(1 - (2/L)^gamma)` gives zero retention at `gamma=0`, breaking backward compatibility. Implemented instead `rho_i(L) = rho*(1 - 2*gamma/L)` with `gamma in [0,1]`: `gamma=0` is exact legacy flat retention, `gamma=1` gives the full dialysis cutoff (L2 unretained, L3 -> rho/3, L4 -> rho/2, L6 -> 2rho/3). Single source of truth: `ReactionNetwork.dilution_outflow(dilution, retention, len_scale)` used by `rhs_cstr`, `jacobian`, `gillespie.build_processes`, and the material-balance diagnostic.
- **Knob B (`min_cat_len`)**: already existed in the generator — no new code, just protocol exposure. It enters `net.params`/`net.hash`/`task_identity` automatically.
- **Knob C (`k_inhibit`)**: implemented as a network property (generator spec -> `ReactionNetwork`), boost_r = eps + kappa*sum w_c x_c^p / (1 + K*sum_{c in C_r} x_c). Quotient-rule derivative added to dense+sparse `d_rates`. SSA: `Process.func_inhibit = (channel_catalysts, K)` divides the hazard by `1 + K*sum n_c/V` in both hazard modes, exactly matching the per-channel ODE denominator (including multi-catalyst channels).
- Provenance: `k_inhibit` in `to_dict`/`params` -> `net.hash`; `retention_len_scale` recorded resolved in `rec["reactor"]`; `min_cat_len` was already in `params`.

### 6.2 Verification

- 110/110 pre-existing tests pass unchanged; 14 new tests in `tests/test_negative_feedback.py` (124 total): outflow formula/selectivity, Jacobian diagonal shifts, flat-retention equivalence at gamma=0, k_inhibit rate law + FD `d_rates`, detailed balance, hash/serialization, SSA `func_inhibit` hazards, invalid-value rejection, `min_cat_len` absorbing limit.

### 6.3 Pilot `neg_feedback_l6` (96 records; `data/pilots/r6/neg_feedback_l6-eb3e57eb04a1/`)

Factorial `rho in {0.8, 0.95}` x `gamma in {0, 0.5, 1.0}` x `min_cat_len in {2, 3}` x `k_inhibit in {0, 0.5}` on the same 4 seeds / L=6 networks as `catret_l6`.

**Result: complete null on high-assembly mass — L>=4 mass is 0.000 in every one of the 96 cells.** Beyond that, each constraint hurt rather than helped:

| Knob | Effect |
|---|---|
| `min_cat_len = 3` | **Kills the network dead everywhere** (m_poly = 0, 0 attractors, all seeds/all cells). With eps=0, producing any len-3 catalyst requires an existing len>=3 catalyst; the RAF exists structurally but the ignition bootstrap fails — catalytic boost at seed concentrations cannot sustain the pathway. The "enzymes must be long" constraint is so strong it prevents ignition entirely. |
| `retention_len_scale` (dialysis) | **Monotonically weakens the organized phase**: net7 at rho=0.95 goes m_poly 4.87 (gamma=0) -> 1.16 (gamma=0.5) -> 0.33 (gamma=1). Starving dimers starves everything downstream — every pathway to L>=4 runs *through* L=2. No L>=4 mass appears at any gamma. |
| `k_inhibit = 0.5` | **Nearly inert** (~0.5% m_poly suppression): at these concentrations the denominator `1 + 0.5*x_c ~ 1.15` barely trims a boost that is already saturated by kappa=500. |

**Interpretation (this is the scientifically important part)**: the dimer "trap" is not merely a pathology — it is also the *entry point and fuel* of the polymer economy. Negative constraints that starve short chains starve the entire organized phase before any long chain can benefit. What the model still lacks is a mechanism where long chains have a **direct kinetic pathway** (fragment ligation — `ligation=True` exists in the generator, joining two polymers in one step) or an **essential functional role** (e.g. a reaction only a long catalyst can run, or energetic stabilization/folding that lowers their free energy). Constraint *strength* cannot substitute for pathway *structure*.
