# chem-life-transition

Toy computational models for a hypothesis from [docs/conversation.md](docs/conversation.md):

> **There exists a nonequilibrium transition in sufficiently connected chemical
> reaction networks at which self-maintaining, heritable organization becomes
> dynamically stable; increasing assembly complexity is an emergent signature of
> crossing this transition rather than the defining property of life itself.**

## The three competing models

The repo exists to discriminate between three hypotheses about what happens
between "chemistry" and "biochemistry":

| Model | Claim | Signature in the data |
|---|---|---|
| **A — cascade** | Several independent thresholds (autocatalysis, self-maintenance, heredity, evolution), each with its own critical point | `λc` estimates for different order parameters are **staggered** |
| **B — single transition** | One deep organizational transition; autocatalysis, heredity, evolvability, high assembly are coupled manifestations of one critical point | `λc` distributions **overlap**; susceptibilities peak together |
| **C — nested** | One primary organizational transition, then secondary transitions *within* the organized phase | One primary `λc` + sub-critical points detectable only in the post-transition ensemble |

## Repo layout

```
src/chemtrans/
  networks/    core.py        ReactionNetwork: stoich matrix S, catalysis matrix K
               generators.py  random_crn / raf_biased / polymer_ligation
  dynamics/    ode.py         mass-action chemostat (workhorse for sweeps)
               gillespie.py   SSA for extinction/persistence at low copy number
  observables/ spectral.py    rho(K) static threshold
               closure.py     organizational closure Phi, RAF core, A(lambda)
               assembly.py    assembly-index distribution P(A)
               persistence.py self-maintenance S, reproductive number R0
  analysis/    critical_points.py  per-observable lambda_c estimation + alignment
experiments/   sweep_lambda.py          phase diagram base dataset
               hysteresis.py            first-order signature (up/down ramps)
               finite_size.py           does lambda_c sharpen with N?
               transition_discrimination.py   Model A vs B vs C — headline test
runs/          NOTEBOOK.md — lab notebook: every run, command, result, reading
docs/          conversation.md — origin conversation (the spec)
               vocabulary.md — glossary: physics/chem terms vs. biological
               meanings ("nucleus" here ≠ the organelle). Read first.
data/          run outputs (jsonl, gitignored)
notebooks/     plots / exploration
```

## The order parameters (per run, per λ)

| Symbol | Meaning | Module |
|---|---|---|
| `Φ` | organizational closure: flux fraction inside self-maintaining (RAF) core | `closure.py` |
| `A(λ)` | autocatalytic activity: catalyzed fraction of production flux | `closure.py` |
| `S(λ)` | self-maintenance: deepest dilution shock survived | `persistence.py` |
| `R₀(λ)` | spectral radius of next-generation matrix (offspring per catalyst-lifetime) | `persistence.py` |
| `E(λ)` | evolvability: fraction of catalyst-edge mutations altering the RAF core | `persistence.py` |
| `H(λ)` | heritable information: template-matching fraction of catalyzed flux (polymer only — `random_crn` has no information carrier) | `heredity.py` |
| `C(λ)` | assembly-complexity distribution P(A); watch the **high-assembly tail** | `assembly.py` |
| `ρ(K)` | spectral radius of catalysis matrix — static threshold predictor | `spectral.py` |

## Generator knobs (the things we fiddle with)

`n_species`, `n_reactions` (finite-size scaling), `p_cat` (catalytic
connectivity → pushes `ρ(K)` through 1), `feed_rate` (energy/matter flux),
`dilution_rate` (mortality pressure), `n_food`, `max_reactants`,
`p_cat_required` (catalysis as precondition vs speedup),
`exclusive_catalysis` (catalysts get NO spontaneous production route →
true nucleation barriers → dynamical bistability), `k_sat`
(Michaelis–Menten saturation; None = unbounded mass-action → runaway
growth), and the generator itself (`random_crn` vs `raf_biased` vs
`polymer_ligation`) — comparing across generators is the universality test.

`polymer_ligation`-specific knobs: `alphabet_size`, `max_length`,
`n_ligations`, `p_info` (fraction of catalysts that are
template-matched — the information channel, decoupled from `p_cat`),
`p_selfrep` (fraction of informational catalysts that are exact
self-templates `cat == product`), `food_len` (strings with
`len ≤ food_len` are food), `cleavage_rate`, `ligation_rate`.
Chemistry: every split `a+b → ab` of every pool string is a ligation
route (reverse = cleavage); catalysis is assigned per product since
template matching is a property of what is assembled.

## Caveats we bake in

- Generator assumptions can manufacture or hide a transition. **Every run
  records full generator params in its output metadata.** A transition that
  only appears under one generator is a modeling artifact until proven
  otherwise.
- ODE steady states can report "alive" at 1e-9 concentration. Persistence
  claims need the Gillespie runs to confirm survival at finite copy number.
- `ρ(K) = 1` is a linearization guess; the nonlinear transition may be
  displaced — that gap is itself a measurement.

## Roadmap

1. **Milestone 1** ✅ — `random_crn` + `ode.simulate` + `closure`/`spectral`/`assembly`; `Φ(λ)` shows a sharp, bimodal transition region at `p_cat ≈ 0.2–0.3` in the dense regime.
2. **Milestone 2** ✅ — `hysteresis.py` + `finite_size.py`. Finding: NO bistability in default regime — the transition is structural percolation across the ensemble, not dynamical; consistent with fixed-degree percolation (λc stable across N).
3. **Milestone 3** ✅ — `exclusive_catalysis` (true nucleation barriers → genuine bistability), `k_sat` MM saturation (bounded kinetics → real steady states), S/R₀/E observables, `critical_points.py` alignment analysis, `transition_discrimination.py` runs end-to-end (first verdict: Model B, weakly supported — observables not yet decoupled).
4. **Milestone 4** ✅ — `polymer_ligation` (string species, all-splits ligation/cleavage chemistry, template catalysis); `molecular_assembly_index` + `H(λ)` measurable; `p_info` decouples information channel from catalytic connectivity. Finding: organization turns on along `p_info` not `p_cat`; sweeping `p_info` gives staggered λc (S,H ≈ 0.2 → Φ,A,R₀ ≈ 0.5 → E,assembly ≈ 0.7–1.0, Model-A-like) while `p_cat` sweeps show coupled behavior — verdict is **axis-dependent**, motivating 2-D critical-manifold mapping (see `runs/NOTEBOOK.md`).
5. **Milestone 5** ✅ — `manifold_2d.py` (the A/B/C question done right: activation boundaries in the 2-D `(p_cat, p_info)` plane, not 1-D slices) + `gillespie.py` (direct SSA, validated). Finding: observables sit on **different activation surfaces** — cascade/nested, not one transition; organized chemistry with zero heredity exists at `(p_cat=0.8, p_info=0)`.

## Question 0 — before A vs B vs C: is there a bifurcation at all?

The A/B/C question only makes sense if the chemistry→organization shift
is a *phase transition / bifurcation* in the dynamical-systems sense —
not just a smooth crossover. The evidence required:

- **Spectral gap closing**: the dominant Jacobian eigenvalue at the
  organized steady state approaches 0 as λ→λc (critical slowing down).
  A crossover has no such signature.
- **Within-network bistability**: a FIXED network shows two attractors
  (seeded organized / unseeded dead) — a saddle-node-like structure.
  NB: this differs from *ensemble percolation* (fraction of random
  networks containing a motif), which is a structural transition of the
  ensemble, not a bifurcation of any one system's dynamics.
- **Susceptibility peak**: ensemble Var(Φ) peaks near the boundary and
  sharpens with system size; a crossover stays flat/smooth.
- **Hysteresis**: adiabatic continuation in a dynamical knob
  (feed_rate, dilution) traces different up/down branches.

`experiments/criticality.py` measures these along a control slice.
Current status (see `runs/NOTEBOOK.md`):

- **Organized states are genuine stable attractors** — resident-Jacobian
  α < 0, no dead-species invasion (the naive full-space Jacobian is
  misleading at extinct species; use the support-restricted one).
- **No dynamical bifurcation signature along `p_info`** (polymer):
  α stays clearly negative, `Var(Φ)` smooth, `frac_bistable = 0` — even
  at `p_cat_required=1`. Reason found: **the cleavage backdoor** —
  all-splits hydrolysis produces any string spontaneously as a fragment
  of a longer one, washing out nucleation barriers in polymer chemistry.
- **Within-network bistability confirmed** only for `exclusive_catalysis`
  `random_crn` (M3): the model family *can* bifurcate, just not along
  the polymer structural axes probed so far.
- **Two-level picture** (worth keeping): (a) bifurcation of a FIXED
  network's dynamics — hysteresis/gap-closing required, seen only in
  exclusive regime; (b) a landscape/percolation transition in the space
  of networks — which attractor *types* become reachable — seen
  everywhere. Claims must name which level they mean.

## Where we are — findings so far (for future selves)

1. **Both transition types exist in the same model**, one knob apart: ensemble percolation of RAF motifs (default `random_crn`) vs dynamical first-order-like bistability (`exclusive_catalysis` — a real nucleation barrier).
2. **Watch out for runaway artifacts**: unbounded mass-action made every early "organized" measurement a mid-explosion snapshot. `k_sat` (MM saturation on the whole propensity) is required for meaningful steady states. Also: `A→A+B` self-replicators were a generator bug, now banned.
3. **Polymer chemistry organizes along the information axis** (`p_info`), not catalytic connectivity (`p_cat`) alone — random catalysis at `p_info=0` never organized except at very high `p_cat`.
4. **The 2-D manifold picture (headline)**: in `(p_cat, p_info)` space the order parameters have *different* activation surfaces —
   - `S`, `R₀`, `n_core`: on almost everywhere (any self-rep motif suffices)
   - `H`: pure `p_info ≈ 0.25` threshold, independent of `p_cat`
   - `Φ`: curved boundary — `p_cat` and `p_info` *substitute* for each other
   - `A`, `E`: only when *both* are high
   - **`(p_cat=0.8, p_info=0)`: Φ≈0.5, RAF cores, but H=0 — organized autocatalytic chemistry with zero heredity.**
   → Evidence favors cascade/nested (Models A/C) over single-transition (B) *in this regime*; caveat: these are smooth activation boundaries, not proven singularities.
5. **Verdicts are slice-dependent**: any single 1-D sweep can look Model-A or Model-B depending on which parameter it rides. A-vs-B discrimination needs the multi-parameter picture.

## Open threads / next steps

1. **Mean-first-passage measurement**: the barrier is a *rate*, not a wall — measure `γ_nucleation(V, cleavage)` and `γ_extinction(V)` directly (Poisson waiting times), which defines the bistable window as `1/γ ≫ t_obs` on both sides. Single-run seeded/unseeded tests overstate hard barriers.
2. **What carries the organized flux?** Are "organized" states dominated by trivial short self-reps (`AA → AA`) or extended collectively-autocatalytic sets? Needs a flux-decomposition diagnostic on polymer nets.
3. **`random_crn` 2-D analog**: `p_cat × p_cat_required` or `p_cat × feed_rate` — does the simple generator also show substitutable axes, or is this polymer-specific?
4. **Tau-leaping SSA** for polymer-scale nets (direct SSA is O(n_rxn)/event — fine at 200 rxns, impractical at 6000); `P(nucleation)` as an order parameter.
5. **Collective-S variant**: `S` currently turns on if *any* motif survives a shock — a version requiring the whole organized state to recover would compare more fairly with `Φ`.
6. **Sharper transition criteria**: current boundaries are logistic fits to observable thresholds; finite-size scaling of the 2-D boundary width would test whether these are true critical surfaces vs smooth crossovers.

## Setup

```bash
pip install -e ".[dev]"
pytest tests/
```
