# Implementation Plan: Symmetry Breaking via Non-Pairing Monomers ($C$) & RNA Folding Competence

This specification provides a complete, context-free technical blueprint for
breaking the combinatorial symmetry of template-directed ligation in `chemtrans` 2.0.
It enables an autonomous agent to implement, verify, and run experimental campaigns
testing whether the combination of a **non-pairing abiotic monomer ($C$)** and
**intramolecular RNA-like folding competence ($f_{\text{fold}}$)** produces genuine
sequence-space pruning, alphabet purification, and heritable sequence selection.

---

## 1. Problem Statement & Theoretical Background

### 1.1 The Null Finding of Campaign `template_lig_l6`
In campaign `template_lig_l6` (audit round `r7`), introducing base-pairing-derived
templated ligation under an alphabet $\Sigma = \{A, B\}$ ($A \leftrightarrow B$)
ignited high polymer mass ($m_{\text{poly}} \approx 3.3$, a $+26\%$ boost over blind ligation).
However, sequence sparsity **collapsed toward zero** ($O_{\text{sparsity}}(L=6) \approx 0.05$,
$N_{\text{eff}} \approx 52$ out of $64$ sequences), inverting the pruning hypothesis.

### 1.2 The Root Cause: Isotropic Graph Symmetry
In a 2-letter alphabet where every $A$ pairs with $B$ and vice versa:
1. **Compositional symmetry**: Every sequence $S$ of length $L$ has an exact, unique
   complement $\bar{S}$ of identical length in the exact same pool.
2. **Topological regularity**: Every ligation product $P_{12}$ has approximately
   the same number of complementary template candidates ($\sim L+1$ when accounting
   for single-mismatch tolerance).
3. **No selection pressure**: Complementarity acts as an isotropic kinetic amplifier
   rather than a filter. Everyone gets boosted equally, flattening the distribution.

### 1.3 The Two Physical Symmetry-Breaking Mechanisms
In real prebiotic chemistry, nature breaks this isotropic symmetry through two universal laws:

1. **Non-Pairing / Prebiotic "Dirty Soup" Monomers ($C$)**:
   Abiotic chemistry does not deliver pure pairing nucleotides. For every Watson-Crick
   pairing base ($A, B$), there are non-pairing organics, amino acids, or damaged bases ($C$).
   - $A \cdot B$ and $B \cdot A$ are pairing matches ($-\Delta g_{\text{pair}}$).
   - Any pairing involving $C$ ($A \cdot C, B \cdot C, C \cdot C$) cannot form hydrogen bonds
     and incurs a severe mismatch penalty ($+\Delta g_{\text{mismatch}}$).
   - **Hypothesis A (Alphabet Purification)**: Templated ligation selectively accelerates
     pure $\{A, B\}$ information polymers, purging $C$-contaminated junk into low concentrations.

2. **Intramolecular Folding Competence (RNA Hairpin Secondary Structure)**:
   In physical chemistry (RNA ribozymes, aptamers), a floppy, disordered linear polymer
   cannot act as an active-site catalyst. To hold two substrate fragments positioned $\sim 2\text{ \AA}$
   apart for ligation, the template strand must collapse into a **rigid 3D fold**.
   - In single-stranded nucleic acids, folding is driven by **intramolecular self-complementarity**
     (stem-loop / hairpin formation).
   - Floppy homopolymers (`AAAAAA`, `BBBBBB`) have zero self-complementarity and cannot fold
     ($f_{\text{fold}} \approx 0$).
   - Palindromic/alternating sequences fold into stable double-helical stems ($f_{\text{fold}} \to 1$).
   - **Hypothesis B (Function-Gated Pruning)**: Gating catalytic potency by folding competence
     ($w_{\text{eff}} = f_{\text{fold}} \cdot w_{\text{pair}}$) eliminates the equal-opportunity
     lottery, concentrating mass onto the rare, self-folding sequence lineages.

---

## 2. Mathematical Formulation

### 2.1 Alphabet & Pairing Rules
Alphabet $\Sigma = \{A, B, C\}$ ($K=3$):
- Pairing relation: $\mathcal{M} = \{(A, B), (B, A)\}$.
- Non-pairing base: $C$.
- For any two bases $x, y \in \Sigma$:
  $$\text{pair\_score}(x, y) = \begin{cases}
  1 & \text{if } (x, y) \in \mathcal{M} \quad (\text{Watson-Crick match}) \\
  0 & \text{otherwise} \quad (\text{mismatch, including all pairs with } C)
  \end{cases}$$

### 2.2 Template-Directed Intermolecular Pairing ($w_{\text{pair}}$)
Let template strand $T$ and ligation product $P_{12} = P_1 + P_2$ both have length $L$:
- Number of matches:
  $$m(T, P_{12}) = \sum_{k=1}^L \text{pair\_score}(T_k, (P_{12})_k)$$
- Number of mismatches:
  $$u(T, P_{12}) = L - m(T, P_{12})$$
- Free energy of duplex hybridization:
  $$\Delta G_{\text{bind}}(T, P_{12}) = - m \cdot \Delta g_{\text{pair}} + u \cdot \Delta g_{\text{mismatch}}$$
- Raw binding weight:
  $$z = m \cdot \Delta g_{\text{pair}} - u \cdot \Delta g_{\text{mismatch}}$$
  $$w_{\text{pair}}(T, P_{12}) = \begin{cases}
  W_{\max} & \text{if } z \ge \ln(W_{\max}) \\
  \exp(z) & \text{otherwise}
  \end{cases}$$

*(When $T$ or $P_{12}$ contains $C$, $u \ge 1$, exponentially penalizing the binding affinity).*

### 2.3 Intramolecular Folding Competence ($f_{\text{fold}}$)
A single strand $T = s_1 s_2 \dots s_L$ folds back on itself into a stem-loop hairpin.
For a hairpin with a minimal central turn loop of size $l_{\text{loop}} = 2$:
- Maximum possible stem pairs: $L_{\text{stem}} = \lfloor (L - l_{\text{loop}}) / 2 \rfloor$.
- For $k = 0, \dots, L_{\text{stem}} - 1$, the $k$-th base from the 5' end pairs with the
  $k$-th base from the 3' end:
  $$m_{\text{stem}}(T) = \sum_{k=0}^{L_{\text{stem}}-1} \text{pair\_score}(T[k], T[L - 1 - k])$$
  $$u_{\text{stem}}(T) = L_{\text{stem}} - m_{\text{stem}}(T)$$
- Folding free energy difference relative to the unfolded coil state:
  $$\Delta G_{\text{fold}}(T) = - m_{\text{stem}}(T) \cdot \Delta g_{\text{fold}} + u_{\text{stem}}(T) \cdot \Delta g_{\text{fold,mm}} + \Delta G_{\text{loop}}$$
  where:
  - $\Delta g_{\text{fold}} \ge 0$ is the stabilization free energy per stem pair (in $k_B T$).
  - $\Delta g_{\text{fold,mm}} \ge 0$ is the destabilization penalty per mismatched stem position (in $k_B T$).
  - $\Delta G_{\text{loop}} > 0$ is the entropic penalty for closing a loop (in $k_B T$, default $1.0$).

- **Folded fraction ($f_{\text{fold}} \in (0, 1)$)**:
  Derived from the two-state Boltzmann partition function:
  $$f_{\text{fold}}(T) = \frac{1}{1 + \exp\left(\Delta G_{\text{fold}}(T)\right)} = \frac{\exp\left(-\Delta G_{\text{fold}}(T)\right)}{1 + \exp\left(-\Delta G_{\text{fold}}(T)\right)}$$

*(If $\Delta g_{\text{fold}} = 0.0$ and $\Delta G_{\text{loop}} = 0.0$, $f_{\text{fold}} \equiv 1.0$, recovering the legacy unfolding-independent behavior).*

### 2.4 Effective Catalytic Weight
The net catalytic potency of strand $T$ on channel $P_1 + P_2 \leftrightarrow P_{12}$ is:
$$w_{\text{eff}}(T, P_{12}) = f_{\text{fold}}(T) \cdot w_{\text{pair}}(T, P_{12})$$
Candidate templates are assigned as catalysts if:
$$w_{\text{eff}}(T, P_{12}) \ge w_{\text{min\_template}}$$
(where $w_{\text{min\_template}} \ge 0.5$, default $1.0$).

### 2.5 Thermodynamic Consistency & Conservation Invariants
1. **Detailed Balance**:
   The catalytic factor $1 + \kappa \cdot w_{\text{eff}} \cdot x_T^p$ multiplies **both** forward
   and reverse rate constants of channel $P_1 + P_2 \leftrightarrow P_{12}$ identically.
   $$k_f / k_r = \exp(-\Delta G_{\text{rxn}})$$
   is strictly preserved across all channels.
2. **Letter Conservation**:
   Internal channel stoichiometry $\mathbf{S}$ satisfies $\mathbf{A} \cdot \mathbf{S} = 0$,
   where $\mathbf{A}_{j, i}$ is the composition count of letter $j \in \{A, B, C\}$ in species $i$.
3. **Product Inhibition Denominator**:
   The inhibition denominator $1 + K_{\text{inh}} \sum_{c \in \mathcal{C}_r} x_c$ uses the
   binary catalyst incidence matrix (`_catalysis_binary()`), never $w_{\text{eff}}$.

---

## 3. Computational Scale & Reactor Regime

### 3.1 Network Size with $K=3$ ($\{A, B, C\}$)
- At $L_{\max} = 6$: Total species $= 3^1 + 3^2 + 3^3 + 3^4 + 3^5 + 3^6 + 4 = 1{,}096$ species;
  $\approx 2{,}754$ ligation channels. This size is computationally heavy for multi-IC scans.
- **Optimal Choice: $L_{\max} = 5$ with $K=3$**:
  - Species: $3 + 9 + 27 + 81 + 243 + 4 = 367$ species.
  - Ligation channels: 567 channels.
  - At $L=4$ and $L=5$, stem-loop hairpins with $l_{\text{loop}} = 2$ have $L_{\text{stem}} = 1$.
- **Alternatively: $K=2$ ($A, B$) at $L_{\max} = 6$ with pure folding**:
  - $L_{\max} = 6$ has $L_{\text{stem}} = 2$, allowing 2-pair stems.
  - 126 species, 266 ligation channels (runs in seconds per cell).

**Strategy**: The code should support **both**:
- Any `alphabet_size` ($K=2$ or $K=3$).
- Hairpin folding score calculated dynamically for any length $L \ge 3$.

---

## 4. File-by-File Implementation Details

### 4.1 `src/chemtrans/networks/generators.py`

#### A. Complementarity Function
Replace the hardcoded 2-letter `_complement` with a general Watson-Crick pair scoring function:
```python
def pair_score(a: str, b: str) -> int:
    """1 if (a, b) is (A, B) or (B, A); 0 otherwise (including C or mismatches)."""
    return 1 if (a == "A" and b == "B") or (a == "B" and b == "A") else 0
```

#### B. Intramolecular Hairpin Folding Function
```python
def hairpin_fold_stability(
    s: str,
    delta_g_fold: float = 0.0,
    delta_g_fold_mm: float = 2.0,
    delta_g_loop: float = 1.0,
    l_loop_min: int = 2,
) -> float:
    """Compute Delta G_fold (in kT) for an intramolecular hairpin.

    Pairs the outermost bases inward: s[k] with s[L - 1 - k]
    for k in 0 .. L_stem - 1, leaving at least l_loop_min in the center.
    Returns Delta G_fold (positive = unfolded coil favored, negative = fold favored).
    """
    L = len(s)
    if L < l_loop_min + 2:
        # Too short to form a stem with a loop
        return 10.0  # Unfolded

    l_stem = (L - l_loop_min) // 2
    if l_stem <= 0:
        return 10.0

    m_stem = sum(pair_score(s[k], s[L - 1 - k]) for k in range(l_stem))
    u_stem = l_stem - m_stem
    return -m_stem * delta_g_fold + u_stem * delta_g_fold_mm + delta_g_loop


def folded_fraction(
    s: str,
    delta_g_fold: float = 0.0,
    delta_g_fold_mm: float = 2.0,
    delta_g_loop: float = 1.0,
) -> float:
    """Boltzmann equilibrium fraction of strand in the folded catalytic state."""
    if delta_g_fold == 0.0 and delta_g_loop == 0.0:
        return 1.0  # Legacy / folding disabled
    dg = hairpin_fold_stability(s, delta_g_fold, delta_g_fold_mm, delta_g_loop)
    # Clip dg to prevent overflow in exp
    if dg > 30.0:
        return 0.0
    if dg < -30.0:
        return 1.0
    return float(1.0 / (1.0 + np.exp(dg)))
```

#### C. `polymer_network(...)` Parameters & Ligation Logic
Add parameters:
- `delta_g_fold: float = 0.0` (folding stabilization per stem pair)
- `delta_g_fold_mm: float = 2.0` (mismatch penalty in stem)
- `delta_g_loop: float = 1.0` (loop closing penalty)

In the ligation channel loop:
```python
if templated_ligation:
    t_idx, t_w = [], []
    for s in strings:
        if len(s) != len(p12):
            continue
        # 1. Folding competence of template T
        f_fold = folded_fraction(s, delta_g_fold, delta_g_fold_mm, delta_g_loop)
        if f_fold < 1e-4:
            continue

        # 2. Intermolecular Watson-Crick pairing against product P12
        m = sum(pair_score(s[k], p12[k]) for k in range(len(s)))
        u = len(s) - m
        z = m * delta_g_pair - u * delta_g_mismatch
        w_pair = w_max_template if z >= np.log(w_max_template) else float(np.exp(z))

        # 3. Combined effective weight
        w_eff = f_fold * w_pair
        if w_eff >= w_min_template:
            t_idx.append(index[s])
            t_w.append(w_eff)
```

Add `delta_g_fold`, `delta_g_fold_mm`, `delta_g_loop` to `generator.params` dictionary so that
`net.hash` deterministically tracks folding parameters.

---

### 4.2 `src/chemtrans/observables/sparsity.py`

Extend the report returned by `sequence_entropy` and `sparsity_report` to include:

1. **Alphabet Purity (for $K=3$)**:
   Fraction of length-$L$ polymer mass residing in pure $\{A, B\}$ sequences (containing zero $C$'s):
   ```python
   def ab_purity(net, x: np.ndarray, length: int) -> float | None:
       """Fraction of length-L mass that contains only A and B (no C)."""
       idx = _len_pool(net, length)
       if len(idx) == 0:
           return None
       tot_mass = float(np.maximum(x[idx], 0.0).sum())
       if tot_mass <= 1e-9:
           return None
       ab_idx = [i for i in idx if "C" not in net.species[i]]
       ab_mass = float(np.maximum(x[ab_idx], 0.0).sum())
       return float(ab_mass / tot_mass)
   ```

2. **Mean Fold Competency ($\langle f_{\text{fold}} \rangle_L$)**:
   Mass-weighted average folding fraction:
   $$\langle f_{\text{fold}} \rangle_L = \sum_{i \in \mathcal{P}_L} p_i \cdot f_{\text{fold}}(\text{species}_i)$$

3. **Reciprocal Pair Co-Enrichment**:
   Check if the top 5 sequences at length $L$ contain reciprocal Watson-Crick pairs ($S_1 \cdot S_2$).

---

## 5. Verification & Test Plan (`tests/test_folding_and_purification.py`)

1. **Legacy Equivalence**:
   - `delta_g_fold=0.0` and `delta_g_loop=0.0` reproduces existing `template_lig_l6` behavior bit-for-bit.
2. **Hairpin Stem Scoring**:
   - At $L=4$: `ABAB` has stem match ($A \cdot B$, $m=1$); `AAAA` has stem mismatch ($A \cdot A$, $m=0$).
   - At $L=6$: `ABCCBA` has $T[0]=A, T[5]=A$ (mismatch); `ABBA` palindrome has maximal stem pairs.
3. **Folding Fraction Monotonicity**:
   - As `delta_g_fold` increases, $f_{\text{fold}}$ for `ABAB` rises from $0.27 \to 0.99$, while $f_{\text{fold}}$ for `AAAA` stays $\approx 0$.
4. **Alphabet $C$ Exclusion**:
   - When $K=3$, sequences containing $C$ have $m=0$ for all $C$ positions.
   - Channel $AA + BB \to AABB$ gets high catalytic weight from $BBAA$, but near-zero weight from $BBCA$.
5. **Detailed Balance & Conservation**:
   - All generated networks with $K=3$ and folding competence pass `net.check_letter_conservation()`
     and `net.check_detailed_balance()`.

---

## 6. Experimental Campaign Protocols

### Campaign 1: Folding Competence in Binary Alphabet (`folding_l6`)
Protocol file: `.devin/campaigns/folding_l6.json`
- Alphabet: $\{A, B\}$ ($K=2$, $L_{\max}=6$).
- Fix: `f_dimer=0.4`, `polymer_retention=0.9`, `k_join=10.0`, `kappa=500.0`, `delta_g_pair=2.0`.
- Sweep:
  - `generator.delta_g_fold`: `[0.0, 1.0, 2.0, 3.5, 5.0]`
- Replicates: 8 network seeds.
- **Target**: Does sequence sparsity $O_{\text{sparsity}}(L=6)$ rise from $0.05$ (at $\Delta g_{\text{fold}}=0$)
  to $\ge 0.60$ as homopolymers are extinguished and mass condenses onto self-folding hairpins?

### Campaign 2: Alphabet Purification with Non-Pairing Monomer $C$ (`purify_l5`)
Protocol file: `.devin/campaigns/purify_l5.json`
- Alphabet: $\{A, B, C\}$ ($K=3$, $L_{\max}=5$).
- Fix: `f_dimer=0.3`, `polymer_retention=0.9`, `k_join=10.0`, `kappa=500.0`, `delta_g_pair=2.0`, `delta_g_fold=2.5`.
- Sweep:
  - `generator.templated_ligation`: `[false, true]`
  - `generator.delta_g_mismatch`: `[1.0, 2.0, 3.0]`
- Replicates: 4 network seeds.
- **Target**: Does templated ligation selectively concentrate mass into pure $\{A, B\}$ sequences,
  driving $\text{Purity}(L=4, 5) \to 1.0$ while $C$-containing polymers are washed out?

---

## 7. Execution Checklist for the Implementer

1. [ ] Implement `pair_score`, `hairpin_fold_stability`, and `folded_fraction` in `src/chemtrans/networks/generators.py`.
2. [ ] Wire `delta_g_fold`, `delta_g_fold_mm`, and `delta_g_loop` into `polymer_network` and verify `net.hash` sensitivity.
3. [ ] Extend `src/chemtrans/observables/sparsity.py` with `ab_purity` and fold score diagnostics.
4. [ ] Write `tests/test_folding_and_purification.py` covering stem scoring, $C$ exclusion, and conservation invariants.
5. [ ] Run `python -m compileall -q src experiments` and full test suite (`pytest -q`).
6. [ ] Author `.devin/campaigns/folding_l6.json` and validate with `experiments/campaign.py validate`.
7. [ ] Execute `folding_l6` campaign across 8 worker processes.
8. [ ] Summarize and report results in `runs/NOTEBOOK.md`.

---

## 8. As-Built Divergences (implemented)

- **Folding-disable semantics**: `folded_fraction` returns 1.0 whenever
  `delta_g_fold <= 0`, *regardless* of `delta_g_loop`. The plan's
  version required both `delta_g_fold=0` AND `delta_g_loop=0`, which
  would have silently enabled the fold gate (loop penalty with no stem
  benefit) under default parameters and broken `template_lig_l6`
  reproducibility. `delta_g_fold` is the master switch: zero stem
  stabilization means there is no folded state to occupy.
- **alphabet_size**: relaxed from `== 2` to `>= 2`; `pair_score` treats
  every letter beyond `B` as non-pairing, so `C` (and any further
  letters) are abiotic junk by construction. Products containing `C`
  can never full-match a template (u >= count of C).
- **Gate is soft, not hard**: `w_eff = f_fold * w_pair` attenuates
  rather than excludes — at `delta_g_fold=5`, L4/L5 coils
  (f_fold ~ 0.047) still pass `w_min=1.5` with w_eff ~ 4.7. Only L6
  true coils (m_stem=0, f_fold ~ 0.007) are fully gated at default
  parameters. Selection strength is continuous in `delta_g_fold`.
- **Observables** (`sequence_entropy` per-length dict):
  `ab_purity` (fraction of pool mass in pure {A,B} sequences),
  `mean_stem_pairs` (mass-weighted hairpin stem matches — syntax-only,
  parameter-free), and `comp_pairs_in_top` (reciprocal Watson-Crick
  pairs within top-k mass leaders). All default to null/[] for
  absent pools.
- **purify_l5 protocol**: uses `f_dimer=0.3` (per plan), which feeds
  ALL nine K=3 dimers uniformly — including C-containing junk dimers,
  i.e. the feed itself is a dirty soup; purification must act against
  a sustained junk influx.
