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
| `H(λ)` | heritable information retention — NOT measurable in `random_crn`, needs `polymer_ligation` | — |
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
4. **Milestone 4** — `polymer_ligation` generator; real assembly index; H(λ) becomes measurable; cross-generator comparison (universality test).
5. **Milestone 5** — `gillespie.py`; noise, extinction, nucleation times, and the fate of the organized state near λc at finite copy number.

## Setup

```bash
pip install -e ".[dev]"
pytest tests/
```
