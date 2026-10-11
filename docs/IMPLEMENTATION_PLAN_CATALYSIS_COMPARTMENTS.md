# Implementation Specification: Length-Weighted Catalysis & Compartmental Retention

> **STATUS: IMPLEMENTED** (see §6 below for the as-built parameter
> locations and the first pilot results).

**Target Branch/Repo**: `chemtrans 2.0` (`/home/dcallan-adm/chem-life-transition`)  
**Objective**: Introduce two physically-grounded prebiotic knobs to overcome the short-polymer trap (where 2-mers dominate all catalytic flux) and the open-dilution barrier (which penalizes long polymers exponentially):
1. **Length-Weighted Catalytic Power** (`cat_len_power`, $\alpha \ge 0$): Catalytic acceleration scales with polymer chain length, modeling the physical requirement of chain length for stable catalytic folding.
2. **Compartmental Polymer Retention** (`polymer_retention`, $\rho \in [0, 1)$): Semi-permeable compartment boundary where monomers/fuel wash through at standard dilution rate $D$, but polymers are retained with efficiency $\rho$, modeling protocells, coacervate droplets, or mineral pores.

Both knobs must be strictly backward-compatible: defaulting to `0.0` reproduces the existing codebase and all 92 unit tests bit-for-bit.

---

## 1. Mathematical & Physical Formulation

### 1.1 Knob 1: Length-Weighted Catalytic Power

#### Current Rate Law
In `ReactionNetwork._rate_vectors`, the catalytic multiplier for channel $r$ with catalyst species set $C_r$ is:
$$\text{boost}_r = \varepsilon + \kappa \sum_{c \in C_r} x_c^p$$
where $p = \text{cat\_power} \in \{1, 2\}$, and $\kappa$ is a constant scalar. Every catalyst $c$ has identical weight $1.0$, regardless of whether $c = \text{'AA'}$ (length 2) or $c = \text{'ABABAB'}$ (length 6).

#### New Length-Weighted Rate Law
Let $L(c) = \text{len}(c)$ be the monomer length of species string $c$. Define the length-dependent weight relative to a reference length $L_{\text{ref}} = 2$ (the minimum possible polymer catalyst):
$$w(c) = \left(\frac{L(c)}{L_{\text{ref}}}\right)^\alpha$$
where $\alpha = \text{cat\_len\_power} \ge 0.0$.

The catalytic boost on channel $r$ becomes:
$$\text{boost}_r = \varepsilon + \kappa \sum_{c \in C_r} w(c) \cdot x_c^p$$

- When $\alpha = 0.0$: $w(c) = 1.0$ for all $c$ (exact backward compatibility).
- When $\alpha = 1.0$: A 2-mer has weight $1.0$; a 4-mer has weight $2.0$; a 6-mer has weight $3.0$ (linear folding reward).
- When $\alpha = 2.0$: A 2-mer has weight $1.0$; a 4-mer has weight $4.0$; a 6-mer has weight $9.0$ (superlinear active-site scaling).

#### Jacobian Consistency
The derivative with respect to catalyst concentration $x_k$ must include $w(k)$:
$$\frac{\partial (\text{boost}_r)}{\partial x_k} = \kappa \cdot p \cdot w(k) \cdot x_k^{p-1} \quad (\text{for } k \in C_r)$$

#### Gillespie SSA Hazard Analogue
In `src/chemtrans/dynamics/gillespie.py`, for a catalyzed process with catalyst $c$:
- **`hazard_mode="participants"`**: The reaction rate constant $k_{\text{eff}}$ receives the multiplier $w(c)$.
- **`hazard_mode="functional"`**: The functional multiplier is $w(c) \cdot \kappa \cdot (n_c / V)^p$.

---

### 1.2 Knob 2: Compartmentalization / Semi-Permeable Retention

#### Current CSTR Dilution
In `src/chemtrans/dynamics/ode.py` (`rhs_cstr`):
$$\frac{dx_i}{dt} = \left(S \cdot v(x)\right)_i + D \cdot \left(f_i^{\text{inlet}} - x_i\right)$$
where $f_i^{\text{inlet}}$ is non-zero only for food monomers ($A, B$). All polymers dilute at the exact same rate $-D \cdot x_i$.

#### New Retention Formulation
Let $\rho \in [0, 1.0)$ be `polymer_retention`. For each species $i$:
$$\rho_i = \begin{cases} 0.0 & \text{if } i \text{ is food / activated monomer / waste } (A, B, A^*, B^*, W) \\ \rho & \text{if } i \text{ is a polymer } (\text{length} \ge 2) \end{cases}$$

The CSTR differential equation becomes:
$$\frac{dx_i}{dt} = \left(S \cdot v(x)\right)_i + D \cdot \left(f_i^{\text{inlet}} - (1 - \rho_i) x_i\right)$$

- When $\rho = 0.0$: All species dilute at rate $D$ (standard CSTR).
- When $\rho \in (0, 1)$: Small food/fuel molecules enter and leave at rate $D$, but synthesized polymers wash out at a reduced effective rate $D_{\text{eff}} = D \cdot (1 - \rho)$.
- When $\rho \to 1.0$: Approaching a closed compartment (polymers only decay through internal thermal hydrolysis).

#### Letter Mass Conservation (Invariant 1 Compliance)
Letter mass balance across the compartment boundary:
$$\frac{dM_{\text{total}}}{dt} = D \left( M_{\text{inlet}} - M_{\text{eff\_out}} \right)$$
where $M_{\text{eff\_out}} = \sum_i (1 - \rho_i) \cdot \text{mass}_i \cdot x_i$.
Internal reactions strictly conserve letters ($S^T \cdot \text{comp} = 0$). The retention knob alters only the physical exchange rates at the boundary, preserving Invariant 1.

#### SSA Gillespie Dilution Processes
In `gillespie.build_processes`:
- Monomer dilution processes remain: rate $= D$, delta $= ((i, -1),)$.
- Polymer dilution processes become: rate $= D \cdot (1 - \rho)$, delta $= ((i, -1),)$.

---

## 2. File-by-File Code Modification Specifications

### 2.1 `src/chemtrans/networks/core.py`

1. **`ReactionNetwork` Data Structure**:
   - Add attribute `cat_len_power: float = 0.0` to `ReactionNetwork.__init__`.
   - Store or precompute species string lengths or catalyst weights vector $w \in \mathbb{R}^{N_{\text{species}}}$:
     ```python
     self.cat_len_power = float(cat_len_power)
     # Lengths of species strings: monomers=1, polymers=len(str)
     lens = np.array([len(s) if not s.endswith("*") and s != "W" else 1 for s in self.species], dtype=float)
     # Reference length L_ref = 2.0 (smallest polymer catalyst)
     self.cat_weights = np.where(lens >= 2, (lens / 2.0) ** self.cat_len_power, 1.0)
     ```
   - Update `_rate_vectors`:
     In dense `_rate_vectors`:
     ```python
     if self.cat_power == 1:
         # Include catalyst length weight
         boost = self.eps + self.kappa * (self.cat_mask * self.cat_weights) @ x
     elif self.cat_power == 2:
         boost = self.eps + self.kappa * (self.cat_mask * self.cat_weights) @ (x ** 2)
     ```
   - Update `jacobian` / `d_rates`:
     Ensure catalytic derivative terms multiply by `self.cat_weights[c]`.
   - Update `to_dict` and `from_dict` to serialize `cat_len_power`.

---

### 2.2 `src/chemtrans/networks/generators.py`

1. **`polymer_network` signature**:
   - Add parameter: `cat_len_power: float = 0.0`.
   - Forward `cat_len_power=cat_len_power` to `ReactionNetwork(...)`.
   - Pass `cat_len_power` into `ReactionNetwork.params` metadata dictionary.

---

### 2.3 `src/chemtrans/dynamics/ode.py`

1. **`rhs_cstr` signature & implementation**:
   - Signature:
     ```python
     def rhs_cstr(net: ReactionNetwork, x: np.ndarray, dilution: float,
                  eps: float | None = None,
                  polymer_retention: float = 0.0) -> np.ndarray:
     ```
   - Boundary loss term:
     ```python
     # Construct retention mask: 0 for food/feed/fuel, polymer_retention for polymers
     if polymer_retention > 0.0:
         ret = np.where(net.is_food | (net.food_inlet > 0), 0.0, polymer_retention)
         washout = (1.0 - ret) * x
     else:
         washout = x
     return net.stoichiometry() @ net.net_flux(x, eps) + dilution * (net.food_inlet - washout)
     ```
2. **`integrate_cstr`**:
   - Accept `polymer_retention: float = 0.0`.
   - Propagate to `rhs_cstr` closures, Jacobian closures, Krylov routines, and root polish.
   - Record `polymer_retention` in `ODEResult.diagnostics`.

---

### 2.4 `src/chemtrans/dynamics/gillespie.py`

1. **`build_processes`**:
   - Accept `polymer_retention: float = 0.0` and read `cat_len_power` from `net.cat_len_power`.
   - For catalytic channels:
     Apply `cat_weight = (len(net.species[c]) / 2.0) ** net.cat_len_power` to the process rate constant (or `func_boost`).
   - For dilution processes:
     ```python
     for i in range(net.n_species):
         is_polymer = (not net.is_food[i]) and (net.food_inlet[i] == 0)
         rho_i = polymer_retention if is_polymer else 0.0
         rate_dil = dilution * (1.0 - rho_i)
         if rate_dil > 0.0:
             procs.append(Process(rate=rate_dil, participants=(i,), offsets=(0,),
                                  delta=((i, -1),), kind="dilution", channel=-1))
     ```
2. **`simulate`**:
   - Accept `polymer_retention: float = 0.0` and propagate to `build_processes`.

---

### 2.5 `src/chemtrans/measure.py`

1. **Task Spec Unpacking**:
   - In `measure_run(task)`:
     - Read `cat_len_power` from `task.get("generator", {}).get("cat_len_power", 0.0)`.
     - Read `polymer_retention` from `task.get("reactor", {}).get("polymer_retention", 0.0)`.
   - Propagate `polymer_retention` to `integrate_cstr`, `_multistability_block`, `_competition_block`, and `gillespie.simulate`.
   - Record both knobs in the output record under `rec["reactor"]["polymer_retention"]` and `rec["network"]["cat_len_power"]`.

---

## 3. Backward Compatibility & Invariant Checklist

1. **Invariant 1: Exact Letter Mass Conservation**:
   - Verify `comp @ S == 0` for all internal channels.
   - Verify that CSTR dilution mass accumulation matches $\frac{dM}{dt} = D(M_{\text{in}} - M_{\text{eff\_out}})$.
2. **Invariant 2: Thermodynamic Consistency**:
   - Length-dependent catalytic weights multiply forward and reverse catalytic rates **identically** ($w(c)$ multiplies the channel multiplier $\kappa$).
   - The ratio $k_f / k_r = \exp(-\Delta G)$ is unaffected. Detailed balance holds exactly.
3. **Invariant 3: Absorbing Limit at $\varepsilon = 0$**:
   - If polymer concentration is 0, all $w(c) x_c^p = 0$. Uncatalyzed synthesis remains impossible at $\varepsilon = 0$. Dead state remains absorbing.
4. **Invariant 4: Structured Dynamic Statuses**:
   - No negative clipping inside the vector field. Negativity excursions continue to trigger `UNRESOLVED`.
5. **Invariant 6: Strict JSON & Provenance**:
   - New parameters must be included in `task_identity` so that altering `cat_len_power` or `polymer_retention` automatically hashes to a fresh `run_id`.

---

## 4. Verification and Test Suite Plan

### Step 1: Regression Baseline
Run existing test suite:
```bash
.devin/env/venv/bin/python -m pytest tests -q
```
Must pass 92/92 tests without modification.

### Step 2: Unit Tests (`tests/test_length_catalysis_and_retention.py`)
1. **`test_cat_len_power_rate_scaling()`**:
   - Construct a network with a 2-mer and a 4-mer catalyst.
   - Compare catalytic flux at $\alpha = 0$ vs $\alpha = 2.0$.
   - Verify the 4-mer channel rate is multiplied by $(4/2)^2 = 4.0$.
   - Verify detailed balance passes.
2. **`test_jacobian_fd_consistency_with_length_power()`**:
   - Finite-difference check of `net.jacobian(x)` with $\alpha = 2.0$ against numerical perturbation of `rhs_cstr`.
3. **`test_polymer_retention_dilution_mass_balance()`**:
   - Run CSTR integration with $\rho = 0.8$.
   - Verify steady-state polymer mass is strictly higher than $\rho = 0.0$.
   - Verify letter mass balance matches the analytical inlet/outlet flux.
4. **`test_ssa_retention_and_length_power()`**:
   - Verify SSA process table carries $(1 - \rho) D$ dilution rates for polymers.
   - Verify exact SSA simulation runs and preserves integer non-negativity.

---

## 5. Pilot Campaign: Testing for High Assembly Index

Create a test protocol: `.devin/campaigns/pilot_assembly_knobs.json`:
- Network: $N=8$ seeds, `max_len=6`, $\kappa = 500$, `cat_power=2`.
- Sweeps:
  - Sweep 1: `generator.cat_len_power` $\in \{0.0, 1.0, 2.0, 3.0\}$.
  - Sweep 2: `reactor.polymer_retention` $\in \{0.0, 0.5, 0.8\}$.
- Observable to monitor: $O_{\text{assembly}}$ (`mass_signature_abs` with $a_{\text{cut}} \ge 4$).

**Expected Hypothesis**:
At $\alpha = 0, \rho = 0$, $O_{\text{assembly}} = 0$ (the existing null).  
As $\alpha \ge 2.0$ and $\rho \ge 0.8$, the combination of high catalytic payoff for long chains and reduced dilution washout will drive non-zero concentration into length-5 and length-6 species, triggering the first non-zero $O_{\text{assembly}}$ signatures in the history of the repo.

---

## 6. As-Implemented Status & First Pilot Results

### 6.1 Where the knobs landed

| Knob | Type | Default | Plumbed through |
|---|---|---|---|
| `cat_len_power` (alpha >= 0) | generator spec | `0.0` | `polymer_network` -> `ReactionNetwork.cat_len_power` -> `cat_weights` -> `rates`, `d_rates` (dense+sparse), `jacobian`, `gillespie.build_processes` (rate constant in both hazard modes), `to_dict`/`from_dict`/`with_channels`, `params` record |
| `polymer_retention` (rho in [0,1)) | reactor spec | `0.0` | `measure_run` -> `ode.integrate_cstr` (`rhs_cstr`, PTC + polish Jacobians, Krylov equilibrium seek), `spectral.jacobian_eigs`, `gillespie.build_processes`/`simulate` (dilution-process rate `D*(1-rho)`), `_multistability_block`, `_continuation_block`, `_competition_block`, abiotic-null calibration |

- Retention mask: `net.polymer_mask` = species with `structure` length >= 2. Food monomers, `A*`, `W`, and unstructured (benchmark) species always dilute at `D`. Serial-dilution (`sd_map`/`sd_step`) is a different reactor protocol and is intentionally unchanged.
- Material-balance diagnostic: at `rho > 0` the scalar closed form no longer applies (internal flux moves letters between food and polymer pools), so the check integrates the predicted boundary rate `D*(feed - sum letter_i*(1-rho_i)*x_i)` along the trajectory (trapezoid; duplicate chunk timestamps contribute zero-width intervals).
- Provenance: `cat_len_power` enters `to_dict` -> `net.hash` and `params`; `polymer_retention` is recorded resolved in `rec["reactor"]`; both are in `task_identity` -> fresh `run_id`s. `AUDIT_ROUND` stays `r6` — defaults are bit-identical to prior code, so unchanged tasks produce unchanged numbers.
- Sweeps: campaign `_set_path` accepts `generator.cat_len_power` and `reactor.polymer_retention` dotted paths with no whitelist changes.

### 6.2 Verification performed

- 92/92 pre-existing tests pass unchanged (bit-for-bit compatibility at defaults).
- 18 new tests in `tests/test_length_catalysis_and_retention.py` (110 total): weight formula, catalyzed-only rate scaling, centered-FD `d_rates` at alpha=2, detailed balance, hash/serialization/knockout provenance, selective rhs/Jacobian outflow, retention material balance, absorbing limit at eps=0+rho, SSA rate tables in both hazard modes, invalid-value rejection.
- End-to-end `measure_run` with both knobs: STATIONARY + certified `x_eq`, stable leading eig, event-time SSA replicates, strict JSON.

### 6.3 Pilot `catret_l6` (64 records; `data/pilots/r6/catret_l6-b619f2902de2/`)

Factorial `cat_len_power in {0,1,2,3}` x `polymer_retention in {0,0.5,0.8,0.95}`, L=6 polymer nets (`max_len=6`, `kappa=500`, `cat_power=2`, `eps=0`, `D=0.5`, dense init), 4 seeds (two shared with `assemblyl6`) + 8-IC multistability per cell.

**What moved:**
- **Retention strongly increases organized polymer mass**: m_poly ~0.28-0.33 at rho=0 -> 2.2-4.9 at rho=0.95 (~10-15x) on ignitable networks.
- **Retention increases attractor multiplicity**: `catret|net1` 2->4 certified attractors, `e1land|net7` 1->3 as rho -> 0.95. First compartment-enables-heredity signal in the repo.
- **Null topology stays null**: `catret|net0` is dead (m_poly=0, 0 attractors) across the whole grid — neither knob manufactures RAF structure.

**What did NOT move (the important null):**
- **High-assembly mass did not turn on.** Length-resolved x_eq on `net7` at rho=0.95: L=2 mass 2.68, L=3 mass 2.19, L=4-6 ~0. Retention accumulates mass but piles it at the *shortest retained lengths* — flat rho removes the washout penalty without creating a length gradient; `g_bond` thermodynamics still penalize each additional bond. `mass_signature_abs` = 0 everywhere.
- **`cat_len_power` has almost no systematic effect** on m_poly or length distribution, and scrambles attractor counts non-monotonically. Mechanistic reading: the *abundant* catalysts are themselves short, so a 6-mer's 9x potency weight is moot when no 6-mer is ever synthesized — the knob reweights catalysts that exist, it cannot bootstrap ones that don't.

**Interpretation (cautious):** the two knobs as designed do not produce high-assembly polymers. Retention alone amplifies whatever the RAF topology already produces — mostly short chains. Candidate next levers suggested by the data: (a) length-*graded* retention rho_i = rho(1 - 1/L_i) or a molecular-weight cutoff so long chains are retained preferentially; (b) `min_cat_len >= 3` hard threshold forcing catalysis through longer species; (c) sequence-specific catalytic competence rather than length-only. Each is a one-line/one-knob variant on the machinery now in place.
