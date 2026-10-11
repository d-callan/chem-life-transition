# Implementation Plan: Physical Base-Pairing, Emergent Templated Ligation & Sequence Pruning

This document provides a self-contained, context-free technical specification
for modeling **base-pairing thermodynamics**, **emergent templated ligation**,
and **sequence-space pruning (entropy collapse)** in the `chemtrans` 2.0 framework.

---

## 1. Problem Statement & Conceptual Rationale

### 1.1 The Shortcoming of Random Catalytic Assignment
In the current network generator (`generators.py`), catalytic competence on
ligation channels ($P_1 + P_2 \leftrightarrow P_{12}$) is assigned via a
random Bernoulli coin flip (`p_cat_lig`). Every sequence has an identical,
sequence-independent probability of acting as a catalyst.
Under this law:
- All sequences of equal length are kinetically interchangeable.
- The network cannot focus mass on specific sequence motifs.
- Sequence diversity remains high-entropy (flat distribution over all $K^L$ isomers).

### 1.2 The Physics of Molecular Recognition (Watson-Crick Pairing)
In physical chemistry (RNA, DNA, peptide $\beta$-sheets), templating is not
an arbitrary rule; it is an **automatic thermodynamic consequence of base-pairing**:
1. Two complementary strands form hydrogen bonds: $A \cdot B$ (and $B \cdot A$).
2. A potential template $T$ binds two fragments $P_1$ and $P_2$ if $T$ is
   complementary to their joined sequence $P_{12} = P_1 + P_2$.
3. When bound, the reactive termini of $P_1$ and $P_2$ are held $\sim 2\text{ \AA}$ apart,
   producing a massive effective local concentration boost (proximity effect).
4. **Length threshold**: At $L \le 3$, binding energy is too weak to overcome
   thermal agitation ($|\Delta G_{\text{bind}}| < k_B T$), so ligation remains unguided.
   At $L \ge 4\text{--}6$, cumulative binding energy crosses the thermal threshold
   ($|\Delta G_{\text{bind}}| \gg k_B T$), causing templated alignment to
   **emerge spontaneously**.

### 1.3 Hypotheses to Test
1. **Hypothesis 1 (Emergent Templating)**: As duplex affinity $\Delta g_{\text{pair}}$
   increases, the network transitions from unguided, random ligation to
   template-directed replication.
2. **Hypothesis 2 (Sequence Pruning & Entropy Collapse)**: The onset of templated
   copying breaks the combinatorial symmetry, pruning away random sequence junk
   and concentrating mass into a narrow, self-maintaining lineage
   (e.g., reciprocal complementary hypercycles: $S_1 \leftrightarrow S_2$).

---

## 2. Mathematical & Physical Formulation

### 2.1 Sequence Complementarity Definition
For alphabet size $K=2$ ($\{A, B\}$):
$$\text{comp}(A) = B, \quad \text{comp}(B) = A$$
For a sequence $S = s_1 s_2 \dots s_L$, its Watson-Crick complement is:
$$\overline{S} = \text{comp}(s_1)\text{comp}(s_2)\dots\text{comp}(s_L)$$
*(Note: in antiparallel nucleic acids strands run in reverse, but in an abstract 2-letter
polymer alphabet direct point-to-point complementarity captures the identical physical
cooperativity).*

### 2.2 Duplex Hybridization Free Energy ($\Delta G_{\text{bind}}$)
For a potential template strand $T$ of length $L_T$ aligned against a target product
$P_{12} = P_1 + P_2$ of length $L_{12}$:
- If $L_T \neq L_{12}$, or $L_{12} < L_{\min,\text{pair}}$:
  $$\Delta G_{\text{bind}}(T, P_{12}) = 0$$
- If $L_T = L_{12} = L$:
  Let $m = \sum_{k=1}^L \mathbb{I}[T_k = \overline{(P_{12})}_k]$ be the number of matches,
  and $u = L - m$ be the number of mismatches.
  $$\Delta G_{\text{bind}}(T, P_{12}) = - m \cdot \Delta g_{\text{pair}} + u \cdot \Delta g_{\text{mismatch}}$$
  *(in thermal units $k_B T$, where $\Delta g_{\text{pair}} \ge 0$ is the per-match binding energy,
  and $\Delta g_{\text{mismatch}} \ge 0$ is the mismatch penalty).*

### 2.3 Proximity-Boost Rate Law
The template $T$ acts as a catalyst on the channel $P_1 + P_2 \leftrightarrow P_{12}$
with catalytic weight:
$$w_{\text{template}}(T, P_{12}) = \exp\left(-\frac{\Delta G_{\text{bind}}(T, P_{12})}{k_B T}\right) = \exp\left(m \cdot \Delta g_{\text{pair}} - u \cdot \Delta g_{\text{mismatch}}\right)$$

To maintain a finite maximum catalytic acceleration and prevent unbounded exponential divergence:
$$w(T, P_{12}) = \min\left(w_{\text{template}}(T, P_{12}), \ W_{\max}\right)$$
where $W_{\max} = 100.0$ (or configurable).

- **Full Complement ($m=L, u=0$)**: Max boost $\exp(L \cdot \Delta g_{\text{pair}})$.
- **Single Mismatch ($u=1$)**: Boost suppressed by factor $\exp(-(\Delta g_{\text{pair}} + \Delta g_{\text{mismatch}}))$.
- **No Match ($m=0$)**: Boost $\le 1.0$ (negligible).

### 2.4 Thermodynamic Consistency & Invariants
- The template acceleration multiplies forward and reverse rates **identically**:
  $$k_f = k_{\text{join}} \cdot (1 + \kappa \cdot w \cdot x_T^p), \quad k_r = k_{\text{join}} e^{\Delta G_{\text{rxn}}} \cdot (1 + \kappa \cdot w \cdot x_T^p)$$
- **Detailed balance** ($k_f / k_r = \exp(-\Delta G_{\text{rxn}})$) is preserved exactly.
- **Letter mass conservation** is unaffected ($T$ is a catalyst, net stoichiometry is $P_1 + P_2 \to P_{12}$).

---

## 3. Observable: Sequence Sparsity & Pruning Metrics

Add a new observable module: `src/chemtrans/observables/sparsity.py`.

### 3.1 Polymer Sequence Entropy ($H_{\text{seq}}$)
Let $\mathcal{P}_L$ be the set of polymer species of length $L$.
The relative concentration distribution across length $L$ is:
$$p_i = \frac{x_i}{\sum_{j \in \mathcal{P}_L} x_j}$$
The Shannon entropy of length $L$ is:
$$H_L = -\sum_{i \in \mathcal{P}_L, p_i > 0} p_i \log_2 p_i$$
The maximum possible entropy is:
$$H_{L,\max} = \log_2(K^L) = L \cdot \log_2(K) = L \text{ bits (for } K=2\text{)}$$

### 3.2 Normalized Sequence Sparsity Index ($O_{\text{sparsity}}$)
$$O_{\text{sparsity}}(L) = 1.0 - \frac{H_L}{H_{L,\max}} \in [0.0, 1.0]$$
- **$O_{\text{sparsity}} \approx 0.0$ (Unguided Soup)**: All sequences at length $L$ have equal
  concentrations ($H_L \approx H_{L,\max}$).
- **$O_{\text{sparsity}} \to 1.0$ (Extreme Pruning / Heredity)**: Mass is concentrated into
  a single sequence or reciprocal pair, while the other $K^L - 2$ sequences are zero.

### 3.3 Effective Diversity ($N_{\text{eff}}$)
$$N_{\text{eff}}(L) = 2^{H_L}$$
Measures the effective number of actively populated sequence variants.

---

## 4. File-by-File Changes Required

### 4.1 `src/chemtrans/observables/sparsity.py` (New File)
Implement:
```python
def sequence_entropy(net: ReactionNetwork, x: np.ndarray, length: int) -> dict:
    """Compute Shannon entropy, effective diversity, and sparsity index for length L."""
```
Returns:
- `entropy`: float ($H_L$)
- `entropy_max`: float ($H_{L,\max}$)
- `sparsity`: float ($O_{\text{sparsity}} \in [0, 1]$)
- `n_eff`: float ($2^{H_L}$)
- `top_sequences`: list of (sequence_string, mass_fraction)

### 4.2 `src/chemtrans/networks/generators.py`
Add parameters to `polymer_network(...)`:
- `templated_ligation: bool = False`
- `delta_g_pair: float = 0.0` (binding energy per base-pair match in $k_B T$)
- `delta_g_mismatch: float = 2.0` (penalty per mismatch)
- `w_max_template: float = 100.0`

When `templated_ligation=True` and `ligation=True`:
- Instead of drawing random catalysts via `v_cat < p_cl`:
- For each product $P_{12}$, compute its Watson-Crick complement string $\overline{P_{12}}$.
- Search the species pool for strands $T$ of length $\text{len}(P_{12})$.
- For each candidate $T$, calculate match count $m$ and mismatch count $u$.
- If $w(T, P_{12}) \ge 1.5$:
  Designate $T$ as an eligible catalyst for channel $P_1 + P_2 \leftrightarrow P_{12}$!
- Add all templated parameters to generator `params` so network hash derives deterministically.

### 4.3 `src/chemtrans/measure.py`
- In `measure_run`:
  Compute sequence sparsity for lengths $L=2, 3, 4, 5, 6$ at equilibrium:
  ```python
  if rec.get("equilibrium", {}).get("x_eq"):
      rec["sparsity"] = {
          f"L{L}": sequence_entropy(net, rec["equilibrium"]["x_eq"], length=L)
          for L in (2, 3, 4, 5, 6)
      }
  ```

---

## 5. Verification & Test Plan

Create `tests/test_templated_ligation.py`:
1. **Legacy Equivalence**: `templated_ligation=False` or `delta_g_pair=0.0` reproduces
   exact existing behavior bit-for-bit.
2. **Complementarity Calculation**: Verify `AABB` matches complement `BBAA` with $m=4, u=0$.
3. **Channel Catalyst Assignment**:
   Verify channel $AA + BB \to AABB$ is catalyzed by `BBAA` when `templated_ligation=True`.
4. **Detailed Balance & Letter Conservation**:
   Ensure all networks with templated ligation pass `check_letter_conservation()`
   and `check_detailed_balance()`.
5. **Entropy & Sparsity Bounds**:
   Verify an equiprobable distribution has $O_{\text{sparsity}} = 0.0$, and a single-sequence
   delta distribution has $O_{\text{sparsity}} = 1.0$.

---

## 6. Experimental Protocol: The Pruning Campaign (`.devin/campaigns/templated_pruning_l6.json`)

Cross:
- `generator.templated_ligation`: `[true]`
- `generator.delta_g_pair`: `[0.0, 1.0, 2.0, 3.0]` (sweeping duplex binding affinity)
- `generator.f_dimer`: `[0.3]` (ensuring environmental feed breaks the short-polymer trap)
- `reactor.polymer_retention`: `[0.9]` (semi-permeable protocell)
- 4 seeds (`e1land|net3`, `e1land|net7`, `catret|net0`, `catret|net1`)

### Success Criteria:
- Does $N_{\text{eff}}(L=4)$ drop from $\approx 16$ at $\Delta g_{\text{pair}}=0$ to $\le 4$ at $\Delta g_{\text{pair}}=3.0$?
- Does $O_{\text{sparsity}}(L=4)$ rise above $0.70$?
- Does the system spontaneously condense onto complementary self-replicating pairs?

---

## 7. Implementation Status (as-built, divergences from §§3–6)

**Implemented.** The mechanism landed with per-(channel, catalyst) edge
weights rather than a separate channel type:

- `Channel.cat_w`: optional tuple aligned with `catalysts`; empty =
  legacy unit weight. Positive finite weights required; length must
  match `catalysts`. Positional constructors unaffected (appended last).
- `net.catalysis_matrix()` now returns the *weighted* incidence
  `C[r,c] = cat_w` (or 1.0). A new `_catalysis_binary()` returns the
  unweighted 0/1 incidence used by the `k_inhibit` denominator
  (`1 + K*sum_{c in C_r} x_c` counts catalyst concentration once,
  never affinity-weighted) — and by `func_inhibit` in the SSA.
- `rates()`, `_d_rates_dense()`, `_d_rates_sparse()`, `_cat_idx_padded()`
  all thread `cat_w`; dense and sparse agree to 1e-10 relative on
  central finite differences.
- SSA `kcat = kappa * cat_weights[c] * edge_w` in both hazard modes —
  edge weight folds into the rate constant, symmetric fwd/rev.
- `to_dict`/`from_dict` serialize `cat_w`; `params` carry all five
  template knobs so `net.hash` changes under any of them.

Generator (`templated_ligation=True`, requires `ligation=True` and
`alphabet_size=2`):
- Ligation fragments are len >= 2, so every templated product
  satisfies len(P12) >= 4 — the plan's length threshold is structural.
- Candidates restricted to len(T) == len(P12) (full-length alignment,
  no sliding); `w = min(exp(m*dg - u*dm), w_max)`, edges kept iff
  `w >= w_min_template` (new knob, default 1.5, plan hardcoded 1.5).
- `w` computed via `z = m*dg - u*dm` with the cap applied on `z` —
  no `exp` overflow for extreme `delta_g_pair`.
- `delta_g_pair=0` (with default mismatch penalty) yields **no**
  template edges at all — this is the "recognition absent" control and
  differs from blind ligation (random catalysts, `templated_ligation=False`).
- The random catalyst lottery (`p_cat_lig`) is bypassed entirely in
  templated mode; `min_cat_len`/`cat_max_len` do not apply to template
  candidates (binding affinity, not the catalyst-pool lottery, is the
  physics).

Observables:
- `src/chemtrans/observables/sparsity.py`: `sequence_entropy(net, x, L)`
  and `sparsity_report(net, x)`. Empty length-pools return
  `sparsity=None` (undefined, not 0/1) — an absent shell is not
  "perfect selection".
- `measure_run` records `equilibrium.sequence_sparsity` (and on
  oscillatory-endpoint equilibria) — plan said `rec["sparsity"]`;
  placed inside the equilibrium block so nonconverged runs never
  report sparsity as if it were a certified equilibrium observable.

Campaign: `.devin/campaigns/template_lig_l6.json` (name differs from
§6's `templated_pruning_l6`):
- Fixes the productive `dimerfeed_l6` regime: `f_dimer=0.4`,
  `polymer_retention=0.9`, `k_join=10`, `kappa=500`, `cat_power=2`,
  `cat_len_power=2`, `min_cat_len=2`, `dilution=0.5`, `multistability`
  12 ICs.
- Sweeps `templated_ligation` {false, true} x `delta_g_pair`
  {0.0, 0.5, 1.0, 2.0, 3.0} x 8 network seeds (`n_networks: 8`,
  `template_lig_l6|net0..7`) — 80 cells. Blind-ligation cells under
  the dg sweep are deterministic replicates (dg inert when blind);
  `templated + dg=0` is the zero-recognition control.
- `delta_g_mismatch`, `w_max_template`, `w_min_template` fixed at
  defaults for the primary sweep; a mismatch sweep is the natural
  follow-up (error-threshold scan).