# Vocabulary (chemtrans 2.0)

This project borrows terms from statistical physics, chemistry, and
biology. Several mean **different things in different fields**. This
document fixes what each term means *here*.

---

## 1. Statistical physics — the transition vocabulary

**Order parameter.** A scalar that is ~0 on one side of a transition
and ~1 on the other. Here the **primary** order parameter is
`mass_raf`: monomer mass in maxRAF species / total non-food monomer
mass. `phi_*` flux shares are diagnostics, not the order parameter.

**Control parameter `λ`.** The swept knob — `p_cat`, `p_selfrep`,
`p_info`, `delta_mu`, `chi`, `eps`, `D`, `V`, reactor geometry.

**Continuous transition.** The order parameter grows smoothly from
zero at `λc`; no coexistence. In a well-mixed reactor, the relevant
mean-field form is the transcritical/imperfect-bifurcation family

    O(λ, ε) = O0 + (A/2)[(λ − λc) + sqrt((λ − λc)^2 + 4 c ε)]

with mean-field exponents β = 1, δ = 2. We do **not** claim
directed-percolation exponents — there is no spatial structure.

**First-order transition / bistability / hysteresis.** Two stable
states at the same parameters; the order parameter jumps; up/down
sweeps disagree. In the code: a `fold` model in
`analysis/critical_points.py` and multistability measured by
attractor clustering in `observables/heredity.py`.

**Crossover.** A smooth, non-singular change — order parameter rises
without any bifurcation or multistability. A legitimate *finding*,
not a failed experiment.

**Absorbing state.** At `eps = 0` with non-food catalysts, the
zero-polymer state cannot be left: no synthesis channel can fire
without a catalyst and no catalyst exists. `eps` plays the role of a
symmetry-breaking field; the `eps -> 0` limit is where a genuine
absorbing-state transition could live.

**Critical slowing down.** The stable fixed point's leading Jacobian
eigenvalue satisfies Re λ₁ → 0⁻ at `λc`; relaxation time diverges.
Measured in `observables/spectral.py` at *certified* equilibria only.

**Nucleation / extinction.** ⚠ collision. *Here*: nucleation = the
rare stochastic event producing the first self-amplifying molecule of
an organized set (a seed), measured as dwell-conditioned first
passage into the organized region. Extinction = the reverse event.
Background-driven nucleation times *decrease* with volume (more
chances); extinction times may scale `ln τ_ext ∝ V`. The two must be
fit asymmetrically (`observables/persistence.py`).

**Thermodynamic limit / finite-size scaling.** `V -> ∞` sharpens a
true transition and leaves a crossover smooth. Volume here is a
copy-number scale, not a spatial extent — there are no spatial
correlations to claim universality over.

## 2. Chemistry — the machinery

**Channel.** One reversible elementary reaction `α ⇌ β` with optional
catalysts (alternatives — any one suffices, never consumed). Catalysis
multiplies both directions by `eps + κ·x_C`, leaving `K_eq` untouched.

**CSTR.** Continuous stirred-tank reactor: `dx/dt = S·v(x) +
D·(inlet − x)`. Open, driven, nonequilibrium — including at
`delta_mu = 0` (which is a *baseline-drive* condition, not equilibrium).

**SD (serial dilution / growth-division).** Closed batch reaction for
a cycle of length τ, then an exact mass-conserving dilution+refeed
map. Fixed points and Floquet multipliers of the map replace the CSTR
fixed-point/eigenvalue story.

**Detailed balance / local detailed balance.** `k_f/k_r = exp(−ΔG)`
per channel, ΔG from species free energies. A *closed* Δμ = 0 system
must reach zero net flux — used as a correctness test. Open-flow Δμ=0
is *not* equilibrium.

**Wegscheider consistency.** Products of K_eq around any reaction
cycle equal 1 — automatic here because all K_eq come from one `g`
vector.

**Falling-factorial hazards.** Exact stochastic kinetics:
`h = k·V^{1−Σα}·Π (n_i)_{α_i}` where `(n)_q` is the falling factorial.
`A + A → …` needs two molecules (`n_A(n_A−1)/V`); if the catalyst is
also a reactant, multiplicities merge (`n_A(n_A−1)(n_A−2)`). Zero
hazard = the event is impossible; impossible events are never
clamped into existence.

**RAF / maxRAF.** A set of catalyzed channels that is (a) food-
generated — every reactant reachable from the food set — and (b)
catalytically closed — every channel has a catalyst produced inside.
`observables/closure.py` computes the unique maximal RAF by iterated
pruning. Structural RAF membership is **not** the same as dynamic
self-maintenance; the code keeps both as separate observables.

**Assembly index.** Number of join steps in the minimum concat-DAG
producing a string, with built blocks reusable. Computed exactly for
short strings (`estimator="exact"`), an explicit upper bound
otherwise. High assembly alone is not a signature — the plan requires
joint high-assembly *and* high-abundance against null-calibrated
thresholds.

**Compositional heredity.** Transmission of the *composition* (mass-
fraction vector) of the organized state across serial-dilution
shocks — NOT survival of individual molecules and NOT sequence
copying (the model has no copying mechanism). Presupposes
multistability (Gate B).

## 3. Measurement hygiene

**Certified stationary.** A fixed point only exists when the
integrator's componentwise residual *and* window-drift checks pass,
optionally followed by a Newton polish that re-certifies the residual.
UNRESOLVED and OSCILLATORY are recorded statuses, never silently
promoted to equilibria.

**Censored.** A stochastic trajectory that reached `t_max` or
`max_events` without the event — recorded with its actual censoring
time, analyzed by right-censoring estimators (KM, Nelson–Aalen, RMST),
never padded or dropped.

**run_id / net_hash.** Content-derived identities; reruns skip
existing run_ids, and a manifest rejects an incompatible protocol
hash rather than mixing outputs.
