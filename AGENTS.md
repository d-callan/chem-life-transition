# AGENTS.md — chemtrans 2.0

## What this repo is

A toy nonequilibrium-physics model asking whether the
chemistry->biochemistry boundary can be meaningfully modeled as a phase
transition. Toy chemistry, real math: conserved letters, reversible
channels, `k_f/k_r = exp(-dG)`, catalysis symmetric across directions,
exact falling-factorial SSA.

## Environment & commands

- Env: `.devin/env/venv` (uv venv, Python 3.12; pin list in
  `.devin/env/requirements.txt`). Interpreter: `.devin/env/venv/bin/python`.
- Tests: `.devin/env/venv/bin/python -m pytest tests -q`
- Compile check: `python -m compileall -q src experiments`
- Campaigns: `experiments/campaign.py {validate|plan|run|summarize}`
  with protocols in `.devin/campaigns/*.json`. Outputs land in
  `data/{validation,pilots,confirmation}/<audit_round>/<name>-<protocol_hash>/` —
  gitignored, resumable, manifest-gated. The `audit_round` epoch
  (`provenance.AUDIT_ROUND`) is injected into every task identity, so
  a bumped round produces fresh run_ids and its own directory —
  pre-remediation artifacts at the flat path are never overwritten or
  pooled with new data.

## Invariants to preserve (do not regress)

1. **Letter mass is conserved exactly.** `comp @ S == 0` for every
   internal channel. No species enters/leaves except via the CSTR feed
   and dilution boundary terms.
2. **Thermodynamic consistency.** `k_r` is derived from species free
   energies; `check_detailed_balance()` must pass. Catalysis is a
   common multiplier `(eps + kappa*w_C*x_C^p)` on both directions;
   `w_C = (len(C)/2)^cat_len_power` is the length-scaled potency
   weight (`cat_len_power=0` is the legacy flat-kappa law) and must
   stay symmetric across directions so `k_f/k_r` is unaffected.
2b. **Compartment retention & dialysis.** `reactor.polymer_retention`
   (rho in [0,1)) makes structural polymers (len>=2) wash out at
   `D*(1-rho)` while food/fuel/waste dilute at `D` — a semi-permeable
   compartment, not a change to internal stoichiometry.
   `reactor.retention_len_scale` (gamma in [0,1]) adds size-selective
   dialysis: `rho_i = rho*(1 - 2*gamma/L)` so gamma=1 lets dimers
   wash out at full `D` while long chains are trapped. Both must be
   threaded identically through `rhs_cstr`, `net.jacobian`,
   `spectral.jacobian_eigs`, and `gillespie.build_processes`/`simulate`
   (dilution-process rates) — `net.dilution_outflow` is the single
   source of truth for the outflow vector.
2c. **Product inhibition.** `generator.k_inhibit` (K >= 0) divides
   each channel's catalytic boost by `1 + K*sum_{c in C_r} x_c` —
   symmetric across directions (detailed balance preserved), carried
   into SSA hazards via `Process.func_inhibit` in both hazard modes.
2d. **Environmental oligomer feed.** `generator.f_dimer` (fraction of
   inlet letter flux in [0,1)) delivers pre-formed dimers in the CSTR
   feed — abiotic background chemistry outside the compartment. Fed
   dimers stay non-food (order parameters count them as polymer mass)
   but are declared via `net.inlet_species`, which relaxes the
   food-only inlet check; SSA feed processes and IC seeding must use
   `net.food | net.inlet_species`. Total inlet letter flux is
   conserved (monomer feed scaled by 1-f_dimer). f_dimer>0 means
   x=0 is no longer absorbing (dimers flow in like food) — that is
   the intended physical semantics.
3. **Absorbing limit.** `eps=0` + non-food catalysts => zero polymer
   synthesis is reachable with no catalysts: dead state absorbing for
   synthesis channels (decay channels still run; they're `synth=False`).
   Exception: `f_dimer>0` supplies dimers at the boundary, so the
   all-zero state is not a fixed point (by design).
4. **Structured statuses.** Dynamics return STATIONARY / OSCILLATORY /
   RUNAWAY / UNRESOLVED; only certified equilibria carry `x_eq`.
   Never use `max(x,0)` inside the vector field to hide negativity.
5. **mass_raf is the ignition order parameter; m_poly is the branch
   order parameter.** `mass_raf` (monomer mass in maxRAF species /
   non-food mass) saturates at ~1 once a system self-maintains — it
   detects ignition but cannot discriminate along organized branches.
   Use `order_parameter: "m_poly"` (nonfood monomer-mass
   concentration, emitted per equilibrium record) for branch-growth
   questions. `phi_*` are diagnostics with explicit denominators.
6. **Strict JSON run records**: no NaN/Infinity — use null + reason.
   run_id/net_hash are content-derived; campaigns resume by skipping
   existing ids and refuse protocol-hash mismatches.
7. **Honest nulls.** A flat-line order parameter is a legitimate
   outcome (model selection includes "null" and it can win). Never
   tune parameters to force a transition.
8. **Edge-weighted catalysis stays symmetric.** `Channel.cat_w`
   (template affinity etc.) is a per-edge multiplier inside the
   common catalytic boost — identical on fwd/rev, never touching
   `k_f/k_r`. The `k_inhibit` denominator uses the *binary* catalyst
   incidence (`_catalysis_binary`), not weighted affinity.
9. **Sequence observables report undefined honestly.** Empty
   length-pools yield `sparsity=null`, never 0/1; a single-sequence
   pool is sparsity=1, a uniform pool is 0.

## Known traps (learned during the rebuild)

- Do not put the sweep seed outside the task spec — seeds must derive
  from `seed_key` strings inside the spec so reordering tasks preserves
  run_ids.
- Reservoir species that are *consumed* by channels cannot emulate a
  clamped bath: use catalyst-role reservoirs (see `benchmarks.schlogl`).
- `np.gradient` on chunk-boundary-duplicated timestamps divides by
  zero; compare integrated letter mass to the closed-form CSTR law
  instead.
- Legacy 1.0 data is archived in commit `116ce1d` for audit only; it
  was generated by a model with broken stoichiometry.
- RAF structure != dynamics: `closure.raf_reactions` is a structural
  fixed point; self-maintenance is a separate question.
