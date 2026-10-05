# chem-life-transition (chemtrans 2.0)

Toy computational model for one question from
[docs/conversation.md](docs/conversation.md):

> **Can the chemistry -> biochemistry boundary be modeled as a
> nonequilibrium phase transition at all — and if so, what kind?**

This is "theoretical chemi-bioinfo-physics": a deliberately toy
chemistry that keeps only the mathematical ingredients needed to ask the
question — mass conservation, thermodynamic consistency, catalytic
closure, a tunable energy drive, stochasticity. It is *not* a claim
about real prebiotic chemistry.

The A/B/C model-discrimination framing from the 1.0 codebase is retired.
What matters now is *which phenomena even occur*:

- continuous (transcritical-like, absorbing-state) transitions
- first-order transitions / bistability / hysteresis
- smooth crossovers
- oscillatory organized regimes

An honest null result — "the corrected model shows only crossovers" — is
a valid answer.

## The chemistry (src/chemtrans/networks/generators.py)

Strings over a finite alphabet, explicit `max_len`, closed pool:

    elongation:    P + A*  <=>  PA + W     (driven; dG = g_bond - delta_mu)
    hydrolysis:    PA      <=>  P + A      (thermal)
    deactivation:  A*      <=>  A + W      (fuel decay)

Rates derive from species free energies `g_i` so `k_f/k_r = exp(-dG)` on
every channel — local detailed balance, Wegscheider-consistent by
construction. Catalysis multiplies *both* directions by the same factor
`eps + kappa*x_C`, so it never alters K_eq. Catalysts are non-food
polymers only, so `eps -> 0` is a genuine absorbing limit for polymer
synthesis.

Controls: `p_cat` (catalytic connectivity), `p_selfrep` (product
catalyzes its own channel), `p_info` (sequence-kin catalysts),
`delta_mu`, `chi` (activated fraction of the feed — total incoming
monomer mass is fixed), `eps` (background synthesis rate), `D`
(dilution), `V` (volume / copy-number scale), reactor `CSTR | SD`.

## Observables

Primary order parameter: `mass_raf` = monomer mass in maxRAF species /
total non-food monomer mass. The maxRAF is the formal Hordijk–Steel
food-generated fixed point (`observables/closure.py`).

Secondary/diagnostic: `phi_gross`, `phi_cat`, `phi_net` (flux
decomposition with explicit denominators), assembly index (exact
concat-DAG join count for short strings, null-calibrated, requiring
joint high-assembly *and* high-abundance), Jacobian relaxation rate
(critical slowing down), Kaplan–Meier / Nelson–Aalen / RMST persistence,
compositional heredity via multistability + serial-transfer similarity.

## Repo layout

```
src/chemtrans/
  networks/   core.py        Channel/ReactionNetwork, thermo checks, hashing
              generators.py  polymer_network(...)  (the baseline chemistry)
              benchmarks.py  transcritical motif (exact lambda_c), Schlogl fold
  dynamics/   ode.py         CSTR integrator + SD map; structured ODEResult
                             with STATIONARY/OSCILLATORY/RUNAWAY/UNRESOLVED
              gillespie.py   exact SSA, falling-factorial hazards, observers
  observables/closure.py  persistence.py  assembly.py  spectral.py
              heredity.py
  analysis/   critical_points.py  candidate-curve fitting (null/ramp/
              sigmoid/hinge/fold/imperfect-transcritical) + AICc,
              network-joint bootstrap
  measure.py  task spec -> strict-JSON run record
  provenance.py
experiments/campaign.py     validate | plan | run | summarize
.devin/campaigns/*.json     protocols (validation, pilot, benchmark)
.devin/env/                 pinned environment
runs/NOTEBOOK.md            fresh lab notebook (2.0)
data/{validation,pilots,confirmation}/   manifest-backed outputs
```

## Running

```
python -m pytest tests -q                                     # gate A
python experiments/campaign.py validate --protocol .devin/campaigns/validation.json
python experiments/campaign.py plan      --protocol .devin/campaigns/pilot.json
python experiments/campaign.py run       --protocol .devin/campaigns/pilot.json --workers 2
python experiments/campaign.py summarize --manifest data/pilots/<id>/manifest.json
```

Campaigns are content-addressed, resumable, and refuse to mix protocols.
Run records are strict JSON (no NaN/Inf) and carry `run_id`, `net_hash`,
seed keys, statuses, and diagnostics. Legacy 1.0 results are archived in
git commit `116ce1d` for audit only — they were produced by a model with
broken stoichiometry and are never evidence.

## What "success" looks like

A defensible claim about the *type* of transition (or absence thereof)
in this toy model: continuous vs first-order vs crossover vs
oscillatory, with critical-slowing-down and finite-size evidence, plus
whether compositional heredity and assembly-complexity signatures ride
on the transition. Not finding a transition is a finding.
