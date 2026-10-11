# Implementation Plan: Multilevel Protocell Selection (Compartment-Metabolism Coupling)

This document provides a complete, context-free technical specification for
implementing **multilevel protocell population dynamics** in `chemtrans` 2.0.
It enables an autonomous agent to model a discrete population of semi-permeable
vesicles that grow, compete for a shared lipid pool via osmotic swelling, lyse,
and divide, testing whether compartment-level selection breaks the combinatorial
commons and drives Darwinian sequence pruning.

---

## 1. Problem Statement & Conceptual Rationale

### 1.1 The "Tragedy of the Commons" and Assembly Dilution
Previous campaigns (`template_lig_l6`, `folding_l6`, `purify_l5`) tested sequence
selection in a single, homogeneous Continuous Stirred-Tank Reactor (CSTR).
All parameter sweeps resulted in **uniform, high-entropy soups** ($N_{\text{eff}} \approx 52\text{--}60$
out of 64 sequences at $L=6$).

Crucially, this produced an **Assembly Theory paradox**:
- The bulk templated reactor successfully synthesized complex macromolecules with
  Molecular Assembly Index $MA = 5$ (e.g. `ABBBAA`, `AAABBB`, `BAABBB`).
- However, because molecular recognition is symmetric in an unconstrained commons,
  that polymer mass was divided democratically across all 64 six-mer isomers.
- As a result, every individual high-assembly molecule remained at trace abundance
  ($x \sim 0.0017$), well below the detection floor ($x_{\text{floor}} = 0.01$).
- In Assembly Theory terms, the signature count collapsed to **$n_{\text{signature}} = 0$**.
  The system produced abiotic assembly complexity, but failed to produce the
  **high abundance** that characterizes living biochemistry.

### 1.2 The Biophysical Hypothesis: Protocell Amplification of High-Assembly Lineages
In experimental prebiotic biophysics (Irene Chen & Jack Szostak, *Science* 2004;
Budin & Szostak, *PNAS* 2011), protocells couple internal chemistry to container survival:
1. **Semi-Permeable Retention**: Small abiotic precursors (activated monomers, dimers)
   permeate through the lipid bilayer; longer polymers ($L \ge 3$) cannot escape.
2. **Osmotic Membrane Tension & Lipid Stealing**:
   - As internal reactions turn permeating food into retained polymers, the internal
     osmotic pressure $\Pi$ rises.
   - Mechanically tense vesicles draw lipids out of the membranes of flaccid, inactive vesicles.
   - Tense vesicles grow; inactive vesicles shrink, lyse (burst), and die.
3. **Growth-Triggered Division & High-Assembly Abundance**:
   - When a vesicle's volume doubles ($V \ge 2 V_0$), hydrodynamic shear forces it
     to pinch off into two daughter vesicles.
   - **Hypothesis**: By coupling replication rate to container survival, protocell
     competition acts as a non-linear amplifier. It takes trace high-assembly
     molecules ($MA=5$) and concentrates a winning lineage into high intra-vesicle
     abundance ($x \ge 0.01, n_{\text{signature}} > 0$), converting an ambiguous
     abiotic soup into a definitive molecular biosignature.

---

## 2. Mathematical & Physical Formulation

### 2.1 The Two-Scale System Architecture
The system consists of an **External Bath** and a discrete **Population of $N(t)$ Vesicles**:
- **External Bath**: A well-mixed reservoir with fixed volume $V_{\text{bath}}$ and
  constant nutrient concentrations $x_{\text{food}}$ and $x_{\text{dimer}}$.
- **Vesicle Population**: A set of $N(t)$ compartments indexed by $j \in \{1, \dots, N(t)\}$.
  Each vesicle $j$ has:
  - Current volume $v_j(t) \in [v_{\min}, v_{\text{div}}]$.
  - Internal concentration vector $\mathbf{x}_j(t) \in \mathbb{R}_{\ge 0}^S$ (where $S$ is the number of species).
  - Absolute molecule count vector $\mathbf{n}_j(t) = v_j(t) \cdot \mathbf{x}_j(t) \cdot N_A$
    (or dimensionless mole count $m_{i, j} = v_j x_{i, j}$).

### 2.2 Membrane Permeation Kinetics
Substrates exchange between the external bath and the vesicle interior across the membrane
surface area $A_j(t) \propto v_j(t)^{2/3}$:
- Species $i$ has a length-dependent membrane permeability coefficient $P_i \ge 0$:
  - Food & monomers ($L=1$): $P_i = P_{\text{fast}}$ (high permeability, e.g. $1.0$).
  - Dimers ($L=2$): $P_i = P_{\text{dimer}}$ (moderate permeability, e.g. $0.2$).
  - Polymers ($L \ge 3$): $P_i = 0.0$ (strictly retained inside the protocell).
- Permeation flux vector $\mathbf{J}_{\text{perm}, j}$:
  $$J_{i, j}^{\text{perm}} = \frac{A_j}{v_j} P_i \left(x_{i, \text{bath}} - x_{i, j}\right)$$

### 2.3 Internal Reaction Kinetics
Inside vesicle $j$, chemistry follows the exact detailed-balance, mass-action rate laws
of the underlying `ReactionNetwork`:
$$\left(\frac{d\mathbf{x}_j}{dt}\right)_{\text{rxn}} = \mathbf{S} \cdot \left(\mathbf{v}_f(\mathbf{x}_j) - \mathbf{v}_r(\mathbf{x}_j)\right)$$
where $\mathbf{S}$ is the stoichiometric matrix and $\mathbf{v}_f, \mathbf{v}_r$ include
the templated ligation and folding boosts already implemented in `core.py`.

### 2.4 Osmotic Swelling and Competitive Volume Growth ($dv_j/dt$)
1. **Internal Osmotic Pressure**:
   $$\Pi_j(t) = \sum_{i \in \text{retained}} x_{i, j}(t) = \sum_{i: L_i \ge 3} x_{i, j}(t)$$
2. **Population Mean Osmotic Pressure**:
   $$\bar{\Pi}(t) = \frac{1}{N(t)} \sum_{k=1}^{N(t)} \Pi_k(t)$$
3. **Membrane Lipid Dynamics (Zero-Sum Competition)**:
   The total lipid pool across all vesicles is finite (carrying capacity $V_{\text{tot}} = N_0 \cdot v_0$).
   A vesicle grows if its internal osmotic pressure exceeds the population average:
   $$\frac{d v_j}{dt} = \mu \cdot v_j \left(\Pi_j - \bar{\Pi}\right) - D_{\text{leak}} v_j$$
   where:
   - $\mu$ is the osmotic growth coupling constant.
   - $D_{\text{leak}}$ is the basal lipid loss / dilution rate.
4. **Volume-Concentration Coupling (Dilution by Growth)**:
   By the product rule on total mass $m_{i, j} = v_j x_{i, j}$:
   $$\frac{d x_{i, j}}{dt} = \left(\frac{d x_{i, j}}{dt}\right)_{\text{rxn}} + J_{i, j}^{\text{perm}} - \frac{1}{v_j}\frac{d v_j}{dt} x_{i, j}$$
   *(Note: when a vesicle grows, $dv_j/dt > 0$, diluting internal concentrations;
   when it shrinks, internal concentrations concentrate until lysis).*

### 2.5 Protocell Lifecycle: Lysis & Division

#### A. Lysis (Death Threshold)
If a vesicle loses lipids such that its volume drops below a critical collapse threshold:
$$v_j(t) \le v_{\text{lysis}} = 0.4 \cdot v_0$$
The vesicle **lyses**:
- Its membrane disintegrates.
- Its internal contents spill into the external waste sink (or bath).
- It is removed from the population ($N \gets N - 1$).

#### B. Division (Reproduction Threshold)
When a vesicle synthesizes enough polymers that its volume doubles:
$$v_j(t) \ge v_{\text{div}} = 2.0 \cdot v_0$$
The vesicle **divides**:
1. It splits into two daughter vesicles $j_a$ and $j_b$, each with initial volume $v_0 = v_j / 2$.
2. **Stochastic Molecular Partitioning**:
   Let $n_i = \text{round}(v_j x_{i, j} \cdot \Omega)$ be the discrete copy count of species $i$
   (where $\Omega$ is a simulation volume scale, e.g. 1000).
   The copies partitioned into daughter $j_a$ follow a Hypergeometric / Binomial distribution:
   $$n_{i, j_a} \sim \text{Binomial}\left(n_i, \ p = 0.5\right), \quad n_{i, j_b} = n_i - n_{i, j_a}$$
   The daughter concentrations become:
   $$x_{i, j_a} = \frac{n_{i, j_a}}{v_0 \cdot \Omega}, \quad x_{i, j_b} = \frac{n_{i, j_b}}{v_0 \cdot \Omega}$$
3. The population gains one vesicle ($N \gets N + 1$).
4. If $N(t) > N_{\max}$ (finite rock-pore carrying capacity), a randomly chosen flaccid
   vesicle is flushed out to keep $N(t) \le N_{\max}$.

---

## 3. Observables: Quantifying Darwinian Protocell Selection

A new observable module: `src/chemtrans/observables/protocells.py`.

### 3.1 Inter-Vesicle Sequence Divergence (Speciation Index)
For two vesicles $j$ and $k$, compute the Bray-Curtis or cosine distance between their
length-6 sequence vectors $\mathbf{p}_j(L=6)$ and $\mathbf{p}_k(L=6)$:
$$d(j, k) = 1.0 - \frac{\mathbf{p}_j \cdot \mathbf{p}_k}{\|\mathbf{p}_j\| \|\mathbf{p}_k\|}$$
The **Population Speciation Index** is the average pairwise distance across all surviving vesicles:
$$\bar{D}_{\text{pop}} = \frac{2}{N(N-1)} \sum_{j < k} d(j, k) \in [0.0, 1.0]$$
- $\bar{D}_{\text{pop}} \approx 0.0$: All vesicles are identical (converged to the flat commons).
- $\bar{D}_{\text{pop}} \to 1.0$: Different vesicles have specialized into distinct, non-overlapping lineages!

### 3.2 Intra-Vesicle Sequence Sparsity ($\bar{O}_{\text{sparsity}}$)
The average sequence sparsity across all active vesicles:
$$\bar{O}_{\text{sparsity}}(L=6) = \frac{1}{N} \sum_{j=1}^N O_{\text{sparsity}, j}(L=6)$$
- **Test**: Does protocell competition drive $\bar{O}_{\text{sparsity}}$ from $\approx 0.04$
  (the CSTR null) up to $\ge 0.70$ (clean sequence pruning)?

### 3.3 High-Assembly Abundance Signature ($n_{\text{sig,pop}}$)
Integrate with the existing `chemtrans.observables.assembly` signature definitions:
- For each surviving protocell $j$, evaluate whether any molecule with
  $MA(s) \ge a_{\text{cut}} = 5$ achieves concentration $x_{s, j} \ge x_{\text{floor}} = 0.01$.
- Record:
  - `n_signature_cells`: Number of protocells possessing a certified high-assembly signature.
  - `max_signature_mass`: Maximum concentration of any $MA \ge 5$ molecule observed across the population.
  - `signature_concordance`: Does high protocell division count correlate with the presence
    of high-assembly signature molecules?

### 3.4 Lineage Survival & Extinction Curves
Track the genealogy of each protocell:
- Number of divisions over time ($N_{\text{div}}(t)$).
- Shannon diversity of surviving maternal lineages (clonal dominance vs neutral drift).

---

## 4. Implementation Structure

### 4.1 New Simulation Module: `src/chemtrans/dynamics/protocell_pop.py`
Create a high-performance vectorized population simulator:
```python
class ProtocellPopulation:
    def __init__(
        self,
        net: ReactionNetwork,
        n_cells: int = 50,
        v0: float = 1.0,
        omega: float = 1000.0,
        mu_growth: float = 1.0,
        p_fast: float = 1.0,
        p_dimer: float = 0.2,
        seed: int = 0,
    ):
        ...
    def step(self, dt: float):
        """Advance the population by dt:
        1. Vectorized ODE step for internal reactions in all N cells.
        2. Substrate permeation from bath.
        3. Osmotic pressure calculation & volume growth dv/dt.
        4. Concentration update with dilution-by-growth.
        5. Check lysis (v < 0.4*v0) -> remove dead cells.
        6. Check division (v >= 2.0*v0) -> stochastic binomial split.
        """
```

### 4.2 Integration into `measure_run` (`src/chemtrans/measure.py`)
Add support for `reactor.type == "protocell_population"`:
- Reads `n_cells`, `t_max`, `mu_growth`, `omega`, `p_dimer`.
- Initializes vesicles by encapsulating random small subsets from the abiotic feed.
- Simulates population for $t \in [0, t_{\max}]$.
- Returns strict JSON record containing:
  - `surviving_cells`: count of live protocells.
  - `total_divisions`: cumulative division events.
  - `mean_sparsity_L6`: average intra-vesicle sequence sparsity.
  - `speciation_index`: inter-vesicle sequence divergence $\bar{D}_{\text{pop}}$.
  - `dominant_lineages`: top 3 sequence motifs across the population.

---

## 5. Verification & Test Suite (`tests/test_protocells.py`)

1. **Permeation Selectivity**:
   - Verify monomers and dimers permeate across the membrane; verify $L \ge 3$ polymers have zero flux.
2. **Osmotic Growth Directionality**:
   - Construct two vesicles: Cell A filled with $L=4$ polymers, Cell B empty.
   - Verify Cell A grows ($dv_A/dt > 0$) and Cell B shrinks ($dv_B/dt < 0$).
3. **Mass Conservation During Division**:
   - When a cell divides, verify $n_{i, j_a} + n_{i, j_b} == n_{i, \text{mother}}$ exactly for all species $i$.
4. **Lysis Boundary**:
   - Verify a cell with negative growth shrinks and is deleted upon crossing $v \le 0.4 v_0$.

---

## 6. Experimental Campaign Protocol: `protocell_selection_l6.json`

File: `.devin/campaigns/protocell_selection_l6.json`
- **Network**: Alphabet $\{A, B\}$, $L_{\max} = 6$, `templated_ligation = True`, `delta_g_pair = 2.0`.
- **Reactor**: `type = "protocell_population"`, `n_cells = 40`, `t_max = 500.0`.
- **Sweeps**:
  - `reactor.mu_growth`: `[0.0, 0.5, 2.0, 5.0]` (sweeping osmotic coupling strength; $\mu = 0$ is the uncoupled control).
  - `reactor.omega`: `[200.0, 1000.0]` (testing small vs large vesicle discrete sampling).
- **Seeds**: 4 network seeds.

### Success Criteria:
1. **The Control ($\mu = 0.0$)**:
   - Zero growth/division; vesicles drift or relax.
   - High assembly ($MA \ge 5$) remains diluted across all 64 sequences.
   - Population signature count: $n_{\text{sig,pop}} = 0$.
2. **The Living Biosignature ($\mu \ge 2.0$)**:
   - Osmotic competition breaks the commons: fast-growing vesicles out-compete and lyse inactive ones.
   - Intra-vesicle sequence sparsity $\bar{O}_{\text{sparsity}}(L=6)$ rises from $0.04 \to \ge 0.50$.
   - **High Assembly in High Abundance**: Specific $MA=5$ sequences concentrate to $x \ge 0.01$
     inside the winning vesicles, producing $n_{\text{sig,pop}} \ge 1$ (the definitive biosignature threshold).
   - Inter-vesicle divergence $\bar{D}_{\text{pop}} \ge 0.70$ (different surviving protocells
     specialize into distinct, non-overlapping genetic lineages).

---

## 7. As-Built Implementation Record (deviations from sections 2-4)

Implemented in `src/chemtrans/dynamics/protocell_pop.py`,
wired into `measure_run` as `reactor.type == "protocell_population"`.

1. **Basal polymer leakage `p_polymer`** (not in §2.2): retained L>=3
   polymers permeate at a small nonzero rate (default 0.01, 100x slower
   than monomers). Rationale: with P=0 the internal polymer mass has no
   outflow — concentration grows without bound even under frequent
   division (division halves molecule counts but preserves
   concentration). Real fatty-acid membranes are imperfect; the leak
   bounds the internal steady state.
2. **Integration**: per-chunk adaptive LSODA per cell with the
   mean-field osmotic coupling pi_bar frozen within a chunk
   (operator splitting), instead of fixed-step RK4. Fixed-step RK4
   overflowed on the stiff templated rates; the per-cell split keeps
   each LSODA problem (S+1)-dimensional instead of one coupled
   (N*S+N)-dimensional system whose numerical Jacobian was
   prohibitive.
3. **Division partitioning**: molecule counts n_i =
   floor(v_j x_i omega + 0.5); daughter A gets Binomial(n_i, 0.5),
   daughter B gets the remainder plus sub-molecule residual mass.
4. **Carrying capacity**: when N > n_max, the lowest-osmotic-pressure
   cells are flushed (implemented, matches plan).
5. **Observables**: `speciation_index` is cosine distance over the
   retained L>=3 pool (plan said L=6 only); `per_cell_sparsity`
   reports L_eff/N_eff per length >= 3; `signature_cells` uses the
   task's assembly a_cut/mass_floor.
6. **Encapsulation init**: each vesicle traps m_init molecules drawn
   multinomial-uniform over the non-food polymer pool plus fresh bath
   concentrations — a random thin slice of abiotic soup, not a sample
   from a converged bulk equilibrium.

### Scout result (single network, t_max=80, n_cells=20)
- mu=0: 20/20 alive, 0 divisions, speciation 0.035 — control clean.
- mu=2.0: 20 divisions, 10 lyses, population grew to n_max=30, but
  speciation ~ 2e-7: osmotic competition selected productive cells
  yet every survivor converged to the same internal lineage.
  Reading: osmotic reward is sequence-blind (it senses total retained
  mass, not identity); internal single-attractor chemistry
  homogenizes cells faster than division-time partition noise
  differentiates them. Sequence speciation may require smaller omega
  (stronger sampling noise — swept) or identity-coupled growth.
