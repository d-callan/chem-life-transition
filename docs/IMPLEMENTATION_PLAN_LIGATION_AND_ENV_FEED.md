# Implementation Plan: Ligation Nucleation Sweeps & Environmental Dimer Influx

This document provides a self-contained, context-free technical specification
for two sequential investigations:
1. **Phase 1: Ligation × Nucleation / Multistability Sweep** using existing
   machinery in the codebase without code modifications.
2. **Phase 2: Environmental Oligomer / Dimer Inflow** introducing an explicit
   feed term for simple oligomers in the reactor inlet and CSTR boundary conditions.

---

## 1. Context & Motivation

### 1.1 The Short-Polymer Trap
Previous campaigns (`assembly_l6`, `catret_l6`, `neg_feedback_l6`) established
that in a pure sequential elongation regime ($A^* + P_L \leftrightarrow P_{L+1} + W$):
- Autocatalytic networks stall at short polymers ($L=2, 3$).
- Increasing catalytic potency on longer chains (`cat_len_power`) does not pull
  mass into $L \ge 4$ because 4-mers are not synthesized in the first place.
- Negative constraints (dialysis `retention_len_scale`, hard catalyst floor
  `min_cat_len=3`) starve the dimer pool and collapse the organized state entirely,
  because every path to length $\ge 4$ must pass sequentially through length 2.

### 1.2 The Ligation Hypothesis & Nucleation Barrier
In `ReactionNetwork` and `polymer_network()`, fragment ligation channels
($P_1 + P_2 \leftrightarrow P_{1+2}$) already exist when `ligation=True`.
Exploratory integration reveals a critical physical property:
- **From dimer-only or unseeded ICs**: The system settles into an $L=2$ state;
  higher lengths remain at zero.
- **When inoculated with $L \ge 4$ polymers**: The high-assembly organized state
  persists at steady state ($m_{L=4} > 0, m_{L=5} > 0$).

Thus, high assembly can be an **alternative stable attractor (bistable / hysteresis)**
separated from the dimer state by a **nucleation barrier**.

### 1.3 The Environmental Oligomer Influx Concept
Real prebiotic compartments (vesicles, mineral pores, coacervate droplets) do
not exist in empty void; they sit in an active geological background where
abiotic, unselective mineral chemistry continuously generates short oligomers.
In such an environment:
- Dimers and trimers are present in the **bulk fluid** and flow or permeate into
  the compartment.
- The compartment's job is not to invent dimers from scratch, but to **trap
  and assemble them into macromolecules**.

---

## 2. Phase 1: Ligation × Nucleation Sweep (Existing Machinery)

### 2.1 Objectives
1. Map the coexistence of low-assembly and high-assembly attractors in $L=6$
   networks with `ligation=True`.
2. Measure whether increasing ligation velocity (`k_join`), polymer retention
   (`polymer_retention`), or catalyst potency scaling (`cat_len_power`) expands
   the basin of attraction of high-assembly states.
3. Test multistability scanning across standard vs. interior (inoculated) initial
   conditions.

### 2.2 Protocol Specification (`.devin/campaigns/ligation_nuc_l6.json`)

Create campaign protocol file:
```json
{
  "name": "ligation_nuc_l6",
  "experiment": "E1_landscape",
  "task": {
    "generator": {
      "kind": "polymer",
      "alphabet_size": 2,
      "max_len": 6,
      "p_cat": 0.5,
      "kappa": 500.0,
      "cat_power": 2,
      "delta_mu": 4.0,
      "chi": 1.0,
      "eps": 0.0,
      "ligation": true,
      "p_cat_lig": 0.5,
      "cat_len_power": 2.0
    },
    "reactor": {
      "type": "cstr",
      "dilution": 0.1,
      "t_max": 300.0,
      "polymer_retention": 0.9,
      "retention_len_scale": 0.0,
      "multistability": {
        "n_ic": 16,
        "n_interior": 8,
        "ic_scale": 0.1,
        "tol": 0.05
      }
    }
  },
  "sweeps": [
    {
      "param": "generator.k_join",
      "values": [1.0, 5.0, 10.0]
    },
    {
      "param": "reactor.polymer_retention",
      "values": [0.8, 0.95]
    },
    {
      "param": "reactor.retention_len_scale",
      "values": [0.0, 0.5]
    }
  ],
  "seed_keys": [
    "e1land|net1",
    "e1land|net3",
    "e1land|net7",
    "catret|net0"
  ]
}
```

### 2.3 Execution & Verification Commands
```bash
# 1. Validate the protocol
.devin/env/venv/bin/python experiments/campaign.py validate .devin/campaigns/ligation_nuc_l6.json

# 2. Plan task count (3 * 2 * 2 * 4 = 48 tasks)
.devin/env/venv/bin/python experiments/campaign.py plan .devin/campaigns/ligation_nuc_l6.json

# 3. Run campaign
.devin/env/venv/bin/python experiments/campaign.py run .devin/campaigns/ligation_nuc_l6.json --workers 4
```

### 2.4 Observable Analysis
Inspect `n_attractors` and the per-attractor mass spectrum:
$$\{L: \sum_{i, \text{len}(i)=L} c_i x_{i, \text{eq}}\}$$
Assert whether any certified attractor in `rec["multistability"]["per_attractor_obs"]`
exhibits:
$$m_{L \ge 4} > 0.01$$

---

## 3. Phase 2: Environmental Oligomer / Dimer Inflow

### 3.1 Physical Rationale & Boundary Term Formulation

In a standard CSTR, `food_inlet` represents the feed concentration $f_i^{\text{inlet}}$.
Currently, the codebase enforces an invariant:
`food_inlet may only supply food species` (`core.py:177`).
In the generator, `food` is defined strictly as $\{A, B, A^*, B^*, W\}$.

To allow an environmental background of simple oligomers without breaking
invariants or thermodynamic consistency:
1. Introduce a configurable inlet feed parameter: `f_dimer` (or `dimer_inlet_fraction`).
2. If `f_dimer > 0`:
   - Dimers ($L=2$) are designated as permeable inlet species (or supplied in the feed).
   - In CSTR mass balance:
     $$\frac{dx_i}{dt} = (S v(x))_i + D \cdot \left(f_i^{\text{inlet}} - (1 - \rho_i) x_i\right)$$
     For dimers in the feed, $f_{\text{dimer}}^{\text{inlet}} > 0$.
   - **Letter conservation**: Total letters fed per unit time is:
     $$F_{\text{letters}} = \sum_i c_i f_i^{\text{inlet}}$$
     If total inlet flux is fixed to $f_{\text{inlet}}$, a fraction $\phi_{\text{dimer}}$ of letters
     is apportioned to activated monomers and $\phi_{\text{dimer}}$ to dimers.
   - **Absorbing limit & reversibility**: The inlet is an open boundary flux, NOT an internal
     stoichiometric channel. It does not alter internal thermodynamic detailed balance ($k_f/k_r = \exp(-\Delta G)$).
     However, when $f_{\text{dimer}} > 0$, the $x=0$ state is no longer an absorbing fixed point
     (dimers are continuously replenished from the environment, exactly like food monomers).
     When $f_{\text{dimer}} = 0$, the exact legacy absorbing limit is preserved.

### 3.2 File-by-File Changes Required

#### 1. `src/chemtrans/networks/core.py`
- Modify validation in `ReactionNetwork.__post_init__`:
  Allow `inlet` (or relaxed `food_inlet`) to supply species that have length $\le 2$
  if explicitly flagged, OR relax the constraint from `nonzero_inlet <= self.food`
  to checking that letter conservation accounts for all inlet species.
  *Preferred clean abstraction*:
  Allow `inlet_species: frozenset[int] | None = None` (defaults to `food`).
  Then validate:
  ```python
  allowed_inlet = self.food if self.inlet_species is None else (self.food | self.inlet_species)
  if not nonzero_inlet <= allowed_inlet:
      raise ValueError("inlet may only supply allowed feed species")
  ```

#### 2. `src/chemtrans/networks/generators.py`
- Add parameter to `polymer_network(...)`:
  `f_dimer: float = 0.0` (fraction or concentration of dimers in the inlet).
  When `f_dimer > 0`:
  - Dimers ($L=2$) are added to `allowed_inlet` (and to `food` if they act as external feed).
  - Inlet vector apportioning:
    Total food inlet letters conserved.
    ```python
    # For alphabet size K, there are K^2 dimers.
    # Each dimer contains 2 letters.
    # Let monomer inlet carry (1 - f_dimer) of inlet letters,
    # and dimer inlet carry f_dimer of inlet letters.
    ```
  - Pass `f_dimer` into network `params` and `to_dict` so network hash reflects the feed.

#### 3. `src/chemtrans/dynamics/ode.py`
- Ensure `material_balance` check correctly accounts for letter flux of dimer inlet:
  $$\text{letter\_flux\_in} = \sum_i \text{comp\_sum}_i \cdot D \cdot f_i^{\text{inlet}}$$
  Already written as:
  `inflow_letters = float(net.letter_mass(net.food_inlet)) * dilution`
  Since `net.letter_mass(x)` evaluates `comp @ x`, this automatically works for any
  species with non-zero composition!

#### 4. `src/chemtrans/dynamics/gillespie.py`
- Ensure SSA feed processes are created for every species with `food_inlet[i] > 0`:
  ```python
  if net.food_inlet[i] > 0:
      procs.append(Process(
          name=f"feed_{net.species[i]}",
          rate=dilution * net.food_inlet[i],
          participants=(),
          deltas=((i, 1),),
          ...
      ))
  ```
  This is ALREADY the implementation in `gillespie.py:153-156`! Any non-zero inlet
  automatically spawns an SSA feed event.

#### 5. `src/chemtrans/measure.py`
- Thread `f_dimer` from generator/reactor task specs into network creation.
- Ensure `f_dimer` is recorded in output records and invalid ranges ($< 0$ or $\ge 1$)
  are rejected.

---

## 4. Verification & Test Plan

### 4.1 Unit Tests (`tests/test_ligation_and_dimer_feed.py`)
1. **Ligation Presence**: Verify `ligation=True` produces exact expected number of
   ligation channels ($P_1 + P_2 \leftrightarrow P_{1+2}$) and matches thermodynamic
   $k_r = k_f \exp(\Delta G)$.
2. **Backward Compatibility**: `f_dimer=0.0` reproduces bit-identical baseline
   behavior and passes all 124 existing unit tests.
3. **Letter Mass Balance with Dimer Feed**:
   Verify CSTR integration with `f_dimer > 0` satisfies:
   $$\frac{d}{dt} \text{letters}_{\text{total}} = D \sum c_i (f_i^{\text{inlet}} - (1-\rho_i) x_i) = 0 \quad \text{at steady state}.$$
4. **SSA Process Invariance**:
   Verify Gillespie process list includes correct Poisson feed events for dimers
   when `f_dimer > 0`.

---

## 5. Escalation: From Single CSTR with Influx to Two-Compartment Model

If Phase 2 demonstrates that environmental dimer feed enables sustained $L \ge 4$
growth inside a retaining reactor, escalate to a formal 2-compartment model:
- **Compartment 1 (Bulk Environment)**: Large volume $V_1$, standard CSTR, $A^* \to$ Dimers.
- **Compartment 2 (Protocell / Mineral Pore)**: Small volume $V_2$, diffusive exchange
  with Compartment 1:
  $$\frac{dx_{i,2}}{dt} = (S v)_2 + k_{\text{perm}, i} \frac{A_{\text{pore}}}{V_2} (x_{i,1} - x_{i,2}) - D_2 (1-\rho_i) x_{i,2}$$
  where $k_{\text{perm}, i}$ is large for monomers and dimers, but zero for $L \ge 4$.
EOF
cat >> runs/NOTEBOOK.md <<'EOF'

---

## Ligation Nucleation & Environmental Inflow Plan

Per `docs/IMPLEMENTATION_PLAN_LIGATION_AND_ENV_FEED.md`:
- Documented Phase 1 (Ligation x Nucleation/Multistability sweep using existing
  ligation=True, k_join, retention, and multistability machinery).
- Documented Phase 2 (Environmental oligomer feed f_dimer in inlet stream,
  extending CSTR boundary term to reflect prebiotic background chemistry).
EOF
echo "Plan written successfully"
---

## 6. As-Implemented Details & Results

### 6.1 Implementation notes

- **`inlet_species` abstraction** (per §3.2 preferred option):
  `ReactionNetwork` accepts `inlet_species` (frozenset of species
  indices) — non-food species the CSTR inlet may legally supply.
  Validation becomes `nonzero_inlet <= food | inlet_species`.
  Serialized in `to_dict`/`from_dict`/`with_channels` -> enters
  `net.hash`. Fed dimers remain non-food: `is_food`, `polymer_mask`,
  and all order parameters are unaffected.
- **`generator.f_dimer` in [0,1)**: fraction of inlet letter flux
  arriving as pre-formed dimers, uniform over all K^2 dimers.
  Monomer feed scaled by (1 - f_dimer) so total inlet letter flux is
  exactly f_inlet. f_dimer=0 is bit-identical legacy.
- **SSA fix**: `build_processes` feed loop and IC seeding now iterate
  `net.food | net.inlet_species` (previously `net.food` only — fed
  dimers would have been silently unsupplied).
- `measure.py` needed no changes: `build_network` signature-filters
  generator specs, `f_dimer` enters `net.params` -> hash -> run_id.
- Semantic note: f_dimer>0 means x=0 is not an absorbing fixed point
  (dimers flow in like food). Intended physics, documented in
  AGENTS.md invariant 3 exception.

### 6.2 Verification

19 new tests in `tests/test_dimer_feed.py`: legacy bit-identity at
f_dimer=0, inlet letter-flux conservation, uniform dimer feed,
monomer rescaling, non-food status, inlet_species declaration/
enforcement, invalid-value rejection, hash/roundtrip, letter
conservation + detailed balance, SSA feed processes + IC seeding,
dimer steady-state elevation, material-balance baseline comparability
(the trapezoid retention diagnostic has O(dt^2) ~0.02-0.03 error on
ALL retention runs including f_dimer=0 — not a feed artifact).

### 6.3 Smoke result (pre-campaign)

seed=7, rho=0.9, alpha=2, kappa=500, eps=0:

| ligation | k_join | f_dimer | L2 | L3 | L4 | L5 | L6 |
|---|---|---|---|---|---|---|---|
| False | - | 0.0 | 0 | 0 | 0 | 0 | 0 (dead basin) |
| False | - | 0.3 | 0.61 | 1.74 | 0 | 0 | 0 |
| False | - | 0.6 | 1.54 | 2.07 | 0 | 0 | 0 |
| True | 10 | 0.0 | 0 | 0 | 0 | 0 | 0 (dead basin) |
| True | 10 | 0.3 | 0.51 | 1.74 | **0.10** | **0.09** | **0.18** |
| True | 10 | 0.6 | 1.29 | 1.79 | **0.38** | **0.30** | **0.20** |

**First nonzero L>=4 mass in any experiment.** Mechanism: fed dimers
are simultaneously the elongation substrate (2 -> 3) and the ligation
monomer pool (2+2 -> 4, 2+3 -> 5, 2+4 -> 6) — the "outside chemistry"
supplies what the interior cannot bootstrap. Feed alone reaches only
L=3; ligation alone dies in the dead basin; the combination ignites
the full length spectrum.

### 6.4 Campaigns

- `ligation_nuc_l6` (48 records): ligation x retention x dialysis x
  k_join WITHOUT environmental feed — nucleation/attractor mapping.
- `dimerfeed_l6` (48 records): f_dimer {0,0.2,0.4} x ligation
  {false,true} x rho {0.8,0.95} x 4 seeds — the decisive interaction.

### 6.5 ligation_nuc_l6 results (48 records, no environmental feed)

`data/pilots/r6/ligation_nuc_l6-5f858aac09d9/`:
k_join {1,5,10} x rho {0.8,0.95} x gamma {0,0.5} x 4 seeds.

- **L>=4 mass is real with ligation alone** (dense/inoculated IC):
  net7 rho=0.95 reaches mL4+ = 2.62 (kj=5), 2.14 (kj=10); catret|net1
  ~0.21. Without ligation every prior campaign had exactly zero.
- **Six certified attractors on net7** (m_poly 2.23-6.43) — all
  reached only from organized/inoculated ICs; the dense-IC
  equilibrium matches the top attractor (m_poly 6.43).
- **Absolute nucleation barrier**: the food-only IC (index 0) is in
  `organized_ic` in ZERO of 48 cells — bare monomer feed never
  reaches any organized attractor regardless of knobs.
- **Dialysis still suppresses**: gamma=0.5 cuts attractor count
  (net7: 6 -> 2) and zeroes mL4+ everywhere — consistent with
  neg_feedback_l6 (starving dimers starves everything downstream).
- **k_join is nearly inert** for attractor counts and mL4+ (1 vs 5 vs
  10 changes little) — ligation's effect is topological (channels
  exist), not rate-limited at these concentrations.

Interpretation: the compartment can HOST high-assembly attractors but
cannot bootstrap them from bare monomers — exactly the "two regimes"
picture. The missing rung is an environmental supply of short
oligomers, which Phase 2 (f_dimer) now tests directly.
