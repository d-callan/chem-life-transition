# Vocabulary

This project borrows terms from statistical physics, chemistry, and
biology. Several of them mean **different things in different fields** —
"the nucleus" being the canonical collision. This document fixes what
each term means *here*, and flags where a biologist's default reading
would be wrong.

Terms are grouped by origin field. Within each, "⚠ collision" marks
words with a competing technical meaning the reader likely knows first.

---

## 1. Statistical physics — the transition vocabulary

**Order parameter.** A scalar that is ~0 on one side of a transition
and ~1 on the other, used to locate the transition and characterize
its nature. Our order parameters are `Φ`, `A`, `S`, `R₀`, `E`, `H`,
and the assembly tail (see §3). *Not* a biological "parameter" — it is
an observable output.

**Critical point / threshold `λc`.** The value of a control parameter
at which an order parameter switches on. Estimated here as the λ where
an observable crosses a cutoff, e.g. `Φ > 0.1`.

**Control parameter `λ`.** The knob we sweep (`p_cat`, `p_info`,
`feed_rate`, ...). Analogous to temperature or pressure in a physical
phase diagram — the axis *we* move along, not something the system
computes.

**First-order vs. continuous transition.** First-order: two stable
phases coexist at the transition; the system jumps discontinuously and
shows hysteresis. Continuous: the order parameter grows smoothly from
zero; no coexistence. Which kind the chemistry→life transition is
(if either) is a core question of this repo.

**Bistability.** Two stable steady states exist at the *same* parameter
values for the *same* network; which one you land in depends on initial
conditions. In our runs: unseeded stays dead, seeded stays organized.
This is the deterministic signature of a first-order transition.

**Hysteresis.** The state depends on history: ramping λ up crosses the
transition at a different point than ramping λ down. Signature of
first-order transitions; measured in `experiments/hysteresis.py`.

**Metastable state.** A stable state that is not the "true" lowest-
free-energy state: it persists indefinitely under small perturbations
but a rare large fluctuation (nucleation) can tip it into the other
phase. Our unseeded alphabet-3 soups are metastable.

**Nucleation.** ⚠ collision. *Here*: the rare stochastic event in
which a metastable dead (unorganized) state produces the **first
molecule of a self-amplifying set** — the seed from which the
organized phase grows. The "nucleus" is the first RAF-closing
catalyst molecule, *not* the cellular organelle. Once one copy
exists, autocatalysis amplifies it; before that copy exists, the
catalyzed reaction has propensity exactly zero. Classical analogy:
the first microscopic crystal that lets a supercooled liquid freeze.

**Nucleation rate `γ_nuc` / time `τ_nuc`.** The Poisson rate of the
seed event, and its inverse mean waiting time. Measured as
first-passage time from a food-only state to `Φ > threshold`
(`experiments/first_passage.py`).

**Extinction rate `γ_ext` / time `τ_ext`.** ⚠ collision. *Here*: the
inverse process — a seeded organized state losing its core species to
demographic noise and falling below `Φ` threshold. Not species
extinction in the ecological sense; it is the death of one stochastic
realization of a chemical network.

**Demographic noise / extinction.** Fluctuations from integer copy
numbers: a species at ~2 molecules can hit 0 by chance. Present in
SSA, absent in ODE (concentrations never truly reach 0).

**Susceptibility.** Ensemble variance of an order parameter,
`Var(Φ)` across networks at fixed λ. Peaks at a true critical point;
rises smoothly for a crossover.

**Finite-size scaling.** How the transition's sharpness changes with
system size `N`. A true phase transition sharpens with N; a crossover
does not. (`experiments/finite_size.py`)

**Critical manifold.** When the transition boundary lives in a
multi-dimensional control-parameter space — here a curve in
`(p_cat, p_info)` — rather than at a single λc on one axis.
(`experiments/manifold_2d.py`)

**Crossover vs. true transition.** The distinction Question 0 exists
to answer. A crossover is a smooth, finite-size-dependent change with
no singularity; a true transition has a bifurcation or a singularity
in the N→∞ limit. Our "Φ turns on" observations can be either — the
order parameter alone cannot tell them apart.

**Fixed-network bifurcation vs. ensemble/structural transition.**
The repo's central subtlety. A dynamical bifurcation happens *within*
one system (its own equations develop a new stable state — detectable
via spectral gap closing, hysteresis, bistability). An ensemble
transition happens *across* systems (the fraction of random networks
whose dynamics ever reach an organized state rises — percolation-
like). Both are legitimate "transitions" but they are different
objects; conflating them is the main interpretive hazard.

**Phase / phase diagram.** ⚠ collision. *Here*: "phase" = a regime of
qualitatively different dynamical behavior (dead vs. organized), not a
cell-cycle phase or a physical phase like liquid/gas — though we use
it deliberately in that physics sense.

## 2. Chemistry / dynamical systems — the machinery

**CRN — chemical reaction network.** A set of species + reactions with
stoichiometry and rate constants. `ReactionNetwork` in
`networks/core.py`.

**Mass-action kinetics.** Reaction rate ∝ product of reactant
concentrations. The default kinetics in `dynamics/ode.py`.

**Chemostat.** The open-reactor idealization our ODE implements:
food species flow in at `feed_rate`, everything dilutes out at
`dilution_rate`. An *open* nonequilibrium system — required because
nothing interesting self-maintains at equilibrium.

**Propensity.** The stochastic event rate of a reaction given current
molecule counts; the SSA analog of a deterministic rate.

**SSA / Gillespie.** `dynamics/gillespie.py` — exact stochastic
simulation in integer molecule counts; rates are the ODE propensities
evaluated at `x = counts/volume`. Where discreteness lives.

**ODE.** Deterministic continuum approximation — real-valued
concentrations. Fast, but silently leaks through barriers that are
real at integer copy number (the fragment-door bug).

**Saturation `k_sat`.** Michaelis–Menten cap on catalyzed rates,
`rate ∝ x_cat·x_sub/(k_sat + x_sub)`. Prevents unbounded mass-action
blowup; itself a sweepable knob.

**RAF set.** *Reflexively autocatalytic and food-generated*: a subset
of reactions where every reaction is catalyzed by a molecule the set
produces, and every reactant is food or produced by the set. The
formal core of "self-maintaining organization" — see `raf_core` in
`observables/closure.py`. The theoretical backbone borrowed from
Kauffman/Hordijk–Steel autocatalytic-set theory.

**Stoichiometry / `S` matrix.** Column per reaction, entries = net
molecular change. `nu=2` matters: `A+A→AA` consumes two copies —
getting this right was a real bug fix.

**Catalysis matrix `K`, spectral radius `ρ(K)`.** `K[i,j] = 1` if
species j catalyzes a reaction producing i. `ρ(K) > 1` is the static
predictor for autocatalytic blowup (`observables/spectral.py`).

**Jacobian / spectral abscissa `α`.** Largest real part of the
linearized dynamics' eigenvalues at a steady state. `α < 0` = stable;
`α → 0` at a bifurcation ("spectral gap closing"). We use the
**resident-restricted** Jacobian (only alive species) plus **invasion
pressure** (can extinct species invade?) because the full Jacobian is
poisoned by extinct species pinned at the `max(x,0)` boundary.

**Seeded vs. unseeded.** Unseeded: start food-only, ask whether
organization emerges on its own. Seeded: start at the ODE-organized
steady state (scaled to counts by volume), ask whether it persists.
The pair diagnoses bistability.

**Exclusive catalysis `p_cat_required`.** Fraction of reactions that
*cannot* run uncatalyzed (rate = 0 without a catalyst). Creates the
chicken-and-egg structure that produces a real barrier: the product
needs a catalyst, the catalyst needs the product.

**Fragment door / backdoor.** Model-specific jargon. Because every
string cleaves at every split, a catalyst-required molecule `P` can
appear as a *fragment* of a longer string `XPY` made by spontaneous
ligation — a backdoor around the exclusivity barrier. In ODE it is a
continuous trickle (door open); in SSA it is a rare Poisson event
(door narrow); at larger alphabet it effectively closes. Not a
deliberate switch — an emergent property of all-splits chemistry.

**Phantom organized state.** An ODE "organized" steady state whose
core species sit at sub-molecular concentrations (x < 1 molecule at
the tested volume). SSA exposes it: the state never actually existed.
A caveat on all ODE-based Φ measurements at small volume.

**Food / `food_len`.** Species supplied by the feed (strings of length
≤ `food_len`, or the food subset in `random_crn`). The nonequilibrium
boundary condition — where matter/energy enter.

## 3. Biology-mapped observables — what we claim corresponds to life

| Symbol | Our definition | Biological concept we're modeling |
|---|---|---|
| `Φ` | Fraction of reaction flux inside the RAF core | Metabolic closure: the fraction of chemistry that is self-produced/self-sustaining |
| `A` | Catalyzed fraction of production flux | Degree to which metabolism is enzyme-driven vs. spontaneous |
| `S` | Deepest dilution shock the organized state survives | Homeostasis / robustness to dilution |
| `R₀` | Offspring per catalyst-lifetime (next-generation spectral radius) | Reproductive number — ⚠ borrowed from epidemiology, analogous role |
| `E` | Fraction of catalyst-edge mutations that alter the RAF core | Evolvability: mutational access to new organization |
| `H` | Template-matching fraction of catalyzed flux | Heredity: catalysis that copies a specific sequence |
| Assembly index | Min. recursive construction steps of a string; P(A) tail | Assembly-theory complexity measure (Cronin/Walker) — molecules that "require a history" |

**Heredity `H`.** ⚠ collision. *Here*: a much weaker property than
biological heredity. `H` measures whether catalysts are sequence-
matched to products — "the chemistry prefersentially copies certain
strings." It has no mutation, no descent, no lineages. A high `H`
means informational catalysis dominates, not that anything "reproduces"
in the biological sense.

**Self-replication / `p_selfrep`.** ⚠ collision. *Here*: a single
reaction `a+b→P` catalyzed by `P` itself — a *molecular* autocatalyst.
Not organismal reproduction, not a self-copying program. It is the
minimal RAF motif: one reaction, one loop.

**Evolution / evolvability `E`.** ⚠ collision. *Here*: static
sensitivity of the RAF core to rewiring one catalytic edge. A proxy
for "how much novel organization is one mutation away" — not
selection, not population dynamics, not adaptation.

**Organization / the "organized phase".** Our umbrella term for the
regime where `Φ` is high: most flux runs through a self-maintaining
core. Deliberately weaker than "life" — it is the physics-level
property that the hypothesis says is the transition's product.

**"Alive" / "dead" (informal).** Shorthand we use for the organized
vs. food-only dynamical phases. *Not* a biological claim — the repo
tests whether the organized phase has life-like properties, not
whether it "is" life.

## 4. Statistical/experimental terms

**Mean first-passage time (MFPT).** Expected time for a stochastic
trajectory to first reach a state. `τ_nuc`, `τ_ext` are MFPTs.

**Censored observation.** A trajectory where the event never occurred
before `t_max` — it tells us `τ > t_max`, not `τ = t_max`. Rates are
estimated via exposure MLE `γ = events/Σ exposure`, which correctly
handles censoring.

**Poisson process / rate.** Events that occur independently at
constant hazard: waiting times exponential, variance = mean. The
model we fit to both nucleation and extinction. "The barrier is a
rate, not a wall" = nucleation is a Poisson waiting time, not a
deterministic impossibility.

**MLE (maximum likelihood estimate).** The parameter value maximizing
the probability of the observed data. Our `γ` estimator is the MLE
under an exponential-hazard model.

**`inf` in rate tables.** No event observed in total exposure —
`γ` is bounded above (≲ 1/exposure), not proven zero. A bound, not
a measurement of impossibility.

**Replicate / rep.** Independent SSA runs of the same network and
condition with different RNG seeds — needed because single stochastic
trajectories are anecdotes. The seed-`4909` lesson.

**Seed.** Two distinct meanings, both in use: RNG `seed` (reproduces
one network instance — e.g., "net seed=4909") and a *seeded*
trajectory (initial-condition sense, §2). Context distinguishes them.

**Sweep / grid / cell.** Systematic scan over control-parameter
values; a "cell" is one (V, cleavage, ...) point of the grid.

---

## The three competing models (project-level vocabulary)

| Model | Claim | Signature |
|---|---|---|
| **A — cascade** | Autocatalysis, self-maintenance, heredity, evolvability switch on at *separate* thresholds | Staggered `λc` per observable |
| **B — single coupled transition** | One organizational transition produces all of them together | Overlapping `λc`s |
| **C — nested transitions** | One deep organizational transition creates a phase in which secondary transitions occur | One primary `λc` + sub-transitions inside the organized phase |

The evidence so far: axis-dependent ordering (staggered along `p_info`,
coupled along `p_cat`), distinct activation surfaces in the 2-D
manifold, and — at sufficient alphabet sparsity — genuine first-order
bistability. Leaning nested/staggered, on a real nucleation barrier.
