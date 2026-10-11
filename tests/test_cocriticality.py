"""Co-criticality apparatus: life_vector schema, multistability and
competition task blocks, ghost-signature regression, and the paired
coincidence analysis on synthetic data."""

import numpy as np
import pytest

from chemtrans.analysis import cocriticality as cc
from chemtrans.measure import measure_run
from chemtrans.networks.generators import polymer_network
from chemtrans.observables import heredity
from chemtrans.observables.life_vector import life_vector

P2_GEN = {"kind": "polymer", "max_len": 4, "alphabet_size": 2,
          "p_cat": 0.5, "p_info": 0.0, "p_selfrep": 0.3,
          "delta_mu": 4.0, "chi": 0.7, "eps": 0.001,
          "kappa": 500.0, "cat_power": 2}


def test_ghost_signature_regression():
    """Dead endpoints (~1e-35 mass) must not normalize into
    attractors when min_mass is enforced."""
    net = polymer_network(max_len=4, alphabet_size=2, p_cat=0.5,
                          seed=7, cat_power=2, kappa=500.0)
    dead = [net.food_inlet + 1e-35 * np.random.default_rng(i).random(
        net.n_species) for i in range(4)]
    cl = heredity.cluster_attractors(net, dead, tol=0.05,
                                     min_mass=0.01)
    assert cl["n_attractors"] == 0 and cl["n_kept"] == 0


def test_life_vector_undefined_semantics():
    # structural-only record: dynamical observables undefined, not 0
    rec = {"structural": {"raf_fraction_channels": 0.2,
                          "n_raf_species": 5}}
    lv = life_vector(rec)
    assert lv["O_raf"]["defined"] and lv["O_raf"]["value"] == 0.2
    for k in ("O_maint", "O_persist", "O_hered", "O_select",
              "O_assembly"):
        assert not lv[k]["defined"] and lv[k]["value"] is None


def test_life_vector_full_record():
    task = {"generator": dict(P2_GEN), "eps": 0.0,
            "seed_key": "e1land|net3",
            "reactor": {"type": "cstr", "dilution": 0.5,
                        "t_max": 400.0,
                        "multistability": {"n_ic": 16,
                                           "min_mass": 0.01,
                                           "t_max": 600.0}}}
    rec = measure_run(task)
    lv = life_vector(rec)
    assert lv["O_raf"]["defined"]
    assert lv["O_maint"]["defined"]          # dead eq is still defined
    assert lv["O_hered"]["kind"] == "discrete"
    assert isinstance(lv["O_hered"]["value"], int)
    assert not lv["O_persist"]["defined"]    # no stochastic block
    assert not lv["O_select"]["defined"]     # competition not requested


def test_multistability_deterministic():
    """Same seed_key -> identical attractor count and centers."""
    task = {"generator": dict(P2_GEN), "eps": 0.0,
            "seed_key": "e1land|net28",
            "reactor": {"type": "cstr", "dilution": 0.5,
                        "multistability": {"n_ic": 16,
                                           "min_mass": 0.01,
                                           "t_max": 600.0}}}
    a = measure_run(task)["multistability"]
    b = measure_run(task)["multistability"]
    assert a["n_attractors"] == b["n_attractors"] >= 2
    assert a["attractor_top_species"] == b["attractor_top_species"]
    assert a["min_mass"] == 0.01


def test_competition_undefined_when_monostable():
    """p=1 net has a single attractor -> competition must report
    incomplete with an explicit reason, never a silent zero."""
    gen = dict(P2_GEN, cat_power=1, kappa=20.0)
    task = {"generator": gen, "eps": 0.0, "seed_key": "e1land|net3",
            "reactor": {"type": "cstr", "dilution": 0.5,
                        "multistability": {"n_ic": 8,
                                           "min_mass": 0.01,
                                           "t_max": 400.0},
                        "competition": {"replicates": 2, "volume": 4,
                                        "t_max": 50.0}}}
    comp = measure_run(task)["competition"]
    assert comp["completed"] is False
    assert comp["reason"] == "fewer_than_two_organized_attractors"


def test_competition_mixture_init_equal_mass():
    """eq_mass protocol seeds equal monomer mass per lineage:
    mass-fraction signatures are converted to molecule concentrations
    by dividing by per-species monomer mass."""
    from chemtrans.dynamics import gillespie
    from chemtrans.measure import _multistability_block, _nonfood_mass, derive_seed
    net = polymer_network(**{k: v for k, v in P2_GEN.items()
                             if k != "kind"},
                          seed=derive_seed("e1land|net28|gen"))
    ms = _multistability_block(net, 0.5, 0.0, "e1land|net28",
                               n_ic=16, min_mass=0.01, t_max=600.0)
    assert ms["n_attractors"] >= 2
    xA = np.array(ms["attractor_x0"]["0"])
    xB = np.array(ms["attractor_x0"]["1"])
    sigA = heredity.attractor_signature(net, xA)
    sigB = heredity.attractor_signature(net, xB)
    m_t = min(_nonfood_mass(net, xA), _nonfood_mass(net, xB))
    comp_mass = net.comp.sum(axis=0)
    vol = 1024.0   # large V keeps integer-rounding error negligible
    n0 = gillespie.counts_from_concentration(
        0.5 * m_t * (sigA + sigB) / np.maximum(comp_mass, 1.0)
        + net.food_inlet, vol)
    xf = n0 / vol
    # realized monomer-mass concentration vector, projected onto the
    # two lineage signature axes: each lineage should carry m_t/2.
    # (The pre-fix code skipped the /comp_mass conversion and realized
    # ~L times too much mass, skewed toward longer species.)
    comp_mass_nf = np.where(net.is_food, 0.0, comp_mass)
    m_real = comp_mass_nf * xf
    assert float(m_real.sum()) == pytest.approx(m_t, rel=0.1)
    G = np.array([[sigA @ sigA, sigA @ sigB],
                  [sigA @ sigB, sigB @ sigB]])
    b = np.array([sigA @ m_real, sigB @ m_real])
    cA, cB = np.linalg.solve(G, b)
    assert float(cA) == pytest.approx(0.5 * m_t, rel=0.1)
    assert float(cB) == pytest.approx(0.5 * m_t, rel=0.1)


def _synthetic_net_series(lc, n_nets=4, n_pts=11, span=(0.5, 1.6)):
    """Fake records: hinge data with critical point lc."""
    recs = []
    rng = np.random.default_rng(0)
    for j in range(n_nets):
        for lam in np.linspace(*span, n_pts):
            recs.append({
                "net_hash": f"net{j}",
                "control": {"lam": float(lam)},
                "life_vector": {
                    "O_maint": {"defined": True,
                                "value": float(
                                    max(0, lam - lc)
                                    + 0.01 * rng.normal())},
                    "O_hered": {"defined": True,
                                "value": int(lam >= lc)},
                }})
    return recs


def test_coincidence_matrix_synthetic():
    """Coincident observables -> coincident_within_tol; separated ones
    -> separated. Verdicts come from per-family paired diffs."""
    recs = _synthetic_net_series(lc=1.0)
    for r in recs:  # add a second observable shifted by +0.4
        lam = r["control"]["lam"]
        r["life_vector"]["O_persist"] = {
            "defined": True,
            "value": float(max(0, lam - 1.4))}
    by_net = cc.group_by_network(recs)
    rep = cc.coincidence_report(
        by_net, observables=("O_maint", "O_persist", "O_hered"),
        n_boot=60, seed=1)
    assert rep["pairwise"]["O_maint|O_persist"]["verdict"] == \
        "separated"
    assert rep["pairwise"]["O_maint|O_hered"]["verdict"] == \
        "coincident_within_tol"
    # per-family paired diffs are reported exactly
    pd = rep["pairwise"]["O_maint|O_hered"]["paired_diffs"]
    assert pd and all(abs(d) <= rep["tol_pair"] for d in pd.values())
    pe = rep["point_estimates"]
    assert pe["O_maint"]["lc"] == pytest.approx(1.0, abs=0.15)
    assert pe["O_hered"]["lc"] == pytest.approx(1.0, abs=0.15)
    # per-family estimates exist for every family
    assert set(rep["family_estimates"]["O_maint"]) == set(by_net)


def test_conditional_observable_never_fit():
    """O_select is conditional on multistability: even when defined it
    is never fit as a global lambda_c curve."""
    recs = _synthetic_net_series(lc=1.0)
    for r in recs:
        lam = r["control"]["lam"]
        r["life_vector"]["O_select"] = {
            "defined": True, "value": float(abs(lam - 1.0) > 0.2)}
    rep = cc.coincidence_report(
        cc.group_by_network(recs),
        observables=("O_maint", "O_select"), n_boot=20, seed=0)
    assert rep["point_estimates"]["O_select"]["lc"] is None
    assert rep["family_estimates"]["O_select"] == {}
    assert rep["pairwise"]["O_maint|O_select"]["verdict"] == \
        "incomparable"


def test_family_id_latent_topology():
    """family_id = latent topology: stable across threshold/energy
    sweeps (p_cat, delta_mu), splits on pool changes (max_len)."""
    base = {"kind": "polymer", "max_len": 4, "alphabet_size": 2,
            "p_cat": 0.5, "delta_mu": 4.0}
    t = lambda **kw: measure_run({
        "generator": dict(base, **kw), "seed_key": "fam|n0",
        "reactor": {"type": "structural"}})
    a, b = t(p_cat=0.3), t(p_cat=0.9)
    c = t(delta_mu=2.0)
    d = t(max_len=5)
    assert a["family_id"] == b["family_id"] == c["family_id"]
    assert a["family_id"] != d["family_id"]
    assert a["net_hash"] != b["net_hash"]   # content differs, family same


def _synthetic_2d(lc_by_level, n_nets=3):
    """Fake 2-D records: hinge at lc that may shift per secondary level."""
    recs = []
    for j in range(n_nets):
        for k in (0.0, 1.0):
            for lam in np.arange(0.0, 4.5, 0.5):
                lc_m = lc_by_level["O_maint"](k)
                lc_a = lc_by_level["O_assembly"](k)
                recs.append({
                    "net_hash": f"net{j}",
                    "family_id": f"fam{j}",
                    "sweep": {"kappa": k},
                    "control": {"lam": float(lam)},
                    "life_vector": {
                        "O_maint": {"defined": True,
                                    "value": float(max(0, lam - lc_m))},
                        "O_assembly": {"defined": True,
                                       "value": float(max(0, lam - lc_a))},
                        "O_hered": {"defined": True,
                                    "value": int(lam >= lc_m)},
                    }})
    return recs


def test_coincidence_2d_tracking_and_drift():
    """Tracking: pinned diffs across levels. Drift: coincidence at
    level 0 that separates at level 1 — the gated-collinearity tell."""
    recs = _synthetic_2d({"O_maint": lambda k: 2.0,
                          "O_assembly": lambda k: 2.0 + 1.5 * k})
    rep = cc.coincidence_2d(recs, "kappa", n_boot=20, seed=0)
    assert rep["status"] == "ok" and rep["n_levels"] == 2
    # O_maint|O_hered pinned at 0 at both levels -> tracking
    assert rep["pairwise_2d"]["O_maint|O_hered"]["verdict"] == "tracking"
    # O_maint|O_assembly coincident at k=0, separated by 2 at k=1
    ent = rep["pairwise_2d"]["O_maint|O_assembly"]
    assert ent["verdict"] == "drift"
    assert ent["n_levels_with_pairs"] == 2
    # lc_curves record the per-family trajectory
    assert "O_maint" in rep["lc_curves"]


def test_coincidence_2d_needs_two_levels():
    recs = _synthetic_2d({"O_maint": lambda k: 2.0,
                          "O_assembly": lambda k: 2.0})
    for r in recs:
        r["sweep"] = {"kappa": 0.0}
    rep = cc.coincidence_2d(recs, "kappa")
    assert rep["status"] == "insufficient_levels"


def test_fold_continuous_sqrt_not_saddle_node():
    """A fold fit with c1 ~ c0 is a CONTINUOUS beta=1/2 onset, not a
    first-order saddle-node — the model name is not the physics."""
    rng = np.random.default_rng(5)
    lam = np.linspace(0, 1, 40)
    y = np.where(lam < 0.55, 0.2,
                 0.2 + 0.5 * np.sqrt(np.maximum(lam - 0.55, 0)))
    y += 0.002 * rng.standard_normal(40)
    est = cc.estimate_lc_continuous(lam, y, seed=0)
    assert est.detail["selected"] == "fold"
    assert est.transition_class == "continuous_half_power"
    assert abs(est.detail["fold_jump"]) < 0.1 * (y.max() - y.min())


def test_fold_real_jump_is_saddle_node():
    """A genuine positive jump discontinuity keeps the saddle_node
    class."""
    rng = np.random.default_rng(2)
    lam = np.linspace(0, 1, 40)
    y = np.where(lam < 0.55, 0.02,
                 0.7 + 0.4 * np.sqrt(np.maximum(lam - 0.55, 0)))
    y += 0.002 * rng.standard_normal(40)
    est = cc.estimate_lc_continuous(lam, y, seed=0)
    assert est.detail["selected"] == "fold"
    assert est.transition_class == "saddle_node"
    assert est.detail["fold_jump"] > 0.1 * (y.max() - y.min())


def test_continuation_monostable_no_hysteresis():
    """A single-valued response curve crossing mass_floor produces one
    fwd flip and one rev flip one grid step apart — that is a bracket,
    NOT a fold pair. The strictly reversible transcritical motif must
    report hysteresis_width=None."""
    from chemtrans.measure import _continuation_block
    from chemtrans.networks.benchmarks import transcritical_motif
    net, _ = transcritical_motif(f_inlet=1.0, k_lig=1.0, kappa=2.0,
                                 keq=np.inf, eps=0.01)
    # organized state as warm start: C present
    x0 = net.food_inlet.copy()
    x0[2] = 0.5
    # dilution axis: organized at low D, dead at high D; D_c = k_f*kappa
    values = [0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2]
    reactor = {"type": "cstr", "dilution": 0.5, "t_max": 400.0,
               "continuation": {"param": "dilution", "values": values,
                                "mass_floor": 0.05}}
    out = _continuation_block(net, reactor, dilution=0.5, eps=0.01,
                              x0_det=x0)
    assert out["hysteresis_width"] is None
    assert out["max_branch_diff"] < 1e-6
    assert "single-valued" in out["note"] or \
        "no fold pair" in out["note"]


def test_sd_solver_failure_is_unresolved_record(monkeypatch):
    """A solve_ivp failure inside sd_map must produce a strict-JSON
    UNRESOLVED record with partial cycle count, not crash."""
    from chemtrans.dynamics import ode as ode_mod
    from chemtrans.provenance import strict_jsonable

    class _Bad:
        success = False
        message = "mocked stiff failure"
        y = None

    monkeypatch.setattr(ode_mod, "solve_ivp", lambda *a, **k: _Bad())
    task = {"generator": dict(P2_GEN), "eps": 0.0,
            "seed_key": "sdfail|net0",
            "reactor": {"type": "sd", "tau": 2.0, "d_factor": 2.0,
                        "n_cycles": 10}}
    rec = measure_run(task)
    assert rec["det"]["status"] == "UNRESOLVED"
    assert rec["det"]["n_cycles_completed"] == 0
    assert "failed" in rec["det"]["message"]
    strict_jsonable(rec)      # must not raise: no NaN/inf in record


def test_serial_transfer_failure_is_unresolved(monkeypatch):
    """serial_transfer_test must return partial cycles + UNRESOLVED
    rather than raising when a batch integration fails."""
    from chemtrans.dynamics import ode as ode_mod
    net = polymer_network(**{k: v for k, v in P2_GEN.items()
                             if k != "kind"}, seed=3)
    x0 = net.food_inlet.copy()
    x0[~net.is_food] = 0.1
    calls = {"n": 0}
    real_ivp = ode_mod.solve_ivp

    class _Bad:
        success = False
        message = "mocked"
        y = None

    def flaky(*a, **k):
        calls["n"] += 1
        if calls["n"] >= 2:
            return _Bad()
        return real_ivp(*a, **k)

    monkeypatch.setattr(ode_mod, "solve_ivp", flaky)
    out = heredity.serial_transfer_test(
        net, x0, tau=1.0, d_factor=2.0, n_cycles=5, eps=0.0)
    assert out["status"] == "UNRESOLVED"
    assert out["n_cycles_completed"] == 1
    assert len(out["dist_post_dilution"]) == 1


def test_tol_pair_bounded_by_swept_domain():
    """Auto equivalence margin is capped at tol_frac of the swept
    domain — a bare 2x-grid default can span 30-40% of the axis and
    call macroscopically separated transitions coincident."""
    recs = _synthetic_net_series(lc=1.0)
    rep = cc.coincidence_report(
        cc.group_by_network(recs),
        observables=("O_maint",), n_boot=20, seed=0)
    # span (0.5,1.6): domain 1.1, grid_res 0.11 -> cap = max(0.11,0.11)
    assert rep["tol_pair"] == pytest.approx(0.11, abs=1e-9)
    assert rep["tol_pair_frac_of_domain"] == pytest.approx(0.1, abs=1e-9)
    assert rep["swept_domain"] == pytest.approx(1.1, abs=1e-9)


def test_coincidence_report_records_dropout_coverage():
    """Every family must report how many records produced a defined
    observable vs dropped out — survivorship at lambda_c made
    visible."""
    recs = _synthetic_net_series(lc=1.0)
    # knock out one family's observable near the critical point
    for r in recs:
        if r["net_hash"] == "net0" and r["control"]["lam"] > 0.95:
            r["life_vector"]["O_maint"] = {"defined": False,
                                           "value": None}
    rep = cc.coincidence_report(
        cc.group_by_network(recs),
        observables=("O_maint",), n_boot=20, seed=0)
    cov = rep["coverage"]["O_maint"]["net0"]
    assert cov["n_dropped"] == 6
    assert cov["n_defined"] == 5


def test_derive_cut_requires_independent_null():
    """assembly.derive_cut must not silently regress its 'abiotic
    null' on the organized run under test — null_source='abiotic'
    integrates the catalysis-removed reference; 'run' is allowed only
    as an explicitly flagged circular calibration."""
    base = {"generator": dict(P2_GEN), "eps": 0.05,
            "seed_key": "acut|net0",
            "reactor": {"type": "cstr", "dilution": 0.5,
                        "t_max": 400.0},
            "assembly": {"derive_cut": True}}
    rec_ab = measure_run({**base,
                          "assembly": {"derive_cut": True,
                                       "null_source": "abiotic"}})
    cal = rec_ab["assembly_calibration"]
    assert cal["null_source"] == "abiotic"
    # self-calibration requires an explicit opt-in and is flagged
    rec_run = measure_run({**base,
                           "assembly": {"derive_cut": True,
                                        "null_source": "run"}})
    cal2 = rec_run["assembly_calibration"]
    assert cal2["null_source"] == "run" and cal2["circular"] is True


def test_pairwise_bootstrap_cells_carry_fdr_correction():
    """Bootstrap-branch cells must report TOST/separation p-values AND
    Benjamini-Hochberg q-values across the pairwise matrix."""
    recs = _synthetic_net_series(lc=1.0, n_nets=4)
    for r in recs:
        lam = r["control"]["lam"]
        r["life_vector"]["O_persist"] = {
            "defined": True,
            "value": float(max(0, lam - 1.4))}
    rep = cc.coincidence_report(
        cc.group_by_network(recs),
        observables=("O_maint", "O_persist", "O_hered"),
        n_boot=50, seed=1, exact_table_max=2, min_pairs=3)
    # n_pairs=4 > exact_table_max=2 -> bootstrap branch -> p-values
    entry = rep["pairwise"]["O_maint|O_persist"]
    assert entry["estimator"] == "paired_bootstrap"
    assert 0.0 <= entry["p_tost"] <= 1.0
    assert 0.0 <= entry["p_separation"] <= 1.0
    # BH-adjusted q values annotated on every tested cell
    tested = [e for e in rep["pairwise"].values() if "p_tost" in e]
    assert tested and rep["n_bootstrap_tested"] == len(tested)
    for e in tested:
        assert 0.0 <= e["q_tost"] <= 1.0
        assert e["q_tost"] >= e["p_tost"] - 1e-12
    assert "Benjamini" in rep["multiple_testing"]


def test_ssa_event_time_observables_and_sensitivity():
    """Stochastic block must report event-time crossings (not grid
    times), exact dwell, and threshold-sensitivity cells from the
    same simulated trajectories."""
    task = {"generator": dict(P2_GEN), "eps": 0.05,
            "seed_key": "evt|net0",
            "reactor": {"type": "cstr", "dilution": 0.5,
                        "t_max": 300.0},
            "stochastic": {"enabled": True, "volume": 30.0,
                           "replicates": 1, "t_max": 80.0,
                           "max_events": 300_000,
                           "record_dt": 0.5},
            "sensitivity": {"stochastic":
                            {"theta_organized": [0.2]}}}
    rec = measure_run(task)
    st = rec["stochastic"]
    assert st["time_resolution"].startswith("event")
    rep = st["replicates"][0]
    assert "n_organized_toggles" in rep
    assert 0.0 <= rep["dwell_organized"] <= 1.0
    if rep["t_first_organized"] is not None:
        assert rep["t_first_organized"] <= rep["censor_time"] + 1e-9
    # sensitivity cell recorded from the same trajectory
    sens = rec.get("sensitivity", {}).get("stochastic")
    assert sens is not None
    assert "0.2|0.01" in sens["cells"]
    assert sens["primary"] == "0.1|0.01"
    assert len(sens["cells"]["0.2|0.01"]["replicates"]) == 1


def test_censor_bracket_bounds_certified_only_lc():
    """Dropped points near lc must widen the reported bracket: imputing
    all censored points low vs high bounds the displacement the
    certified-only estimate hides."""
    rng = np.random.default_rng(0)
    lam = np.linspace(0, 1, 40)
    y = np.where(lam < 0.5, 0.05, 0.05 + 0.9 * (lam - 0.5))
    y += 0.003 * rng.standard_normal(40)
    base = cc.estimate_lc_continuous(lam, y, seed=0)
    # 5 records dropped right at the transition
    drops = [0.48, 0.50, 0.52, 0.54, 0.56]
    est = cc.estimate_lc_continuous(lam, y, seed=0, censor_lams=drops)
    assert est.detail["n_censored"] == 5
    br = est.detail["lc_censor_bracket"]
    assert br[0] <= br[1]
    # the certified-only lc sits inside the censoring bracket
    assert br[0] <= base.lc <= br[1]
    # and a real drop near lc produces a non-vanishing bracket
    assert br[1] - br[0] > 0.05


def test_audit_round_labels_records_and_run_id():
    """Records carry the audit epoch and identical tasks hash to fresh
    run_ids per round — pre-remediation shards can't satisfy resume."""
    from chemtrans.provenance import AUDIT_ROUND, run_id
    task = {
        "generator": {"kind": "polymer", "max_len": 3,
                      "alphabet_size": 2, "p_cat": 0.4, "kappa": 5.0,
                      "seed": 11},
        "reactor": {"type": "cstr", "dilution": 1.0, "t_max": 50.0},
        "seed_key": "audit|t1",
    }
    rec = measure_run(task)
    assert rec["audit_round"] == AUDIT_ROUND
    # the epoch participates in the id: hashing the same task spec
    # with a different round must produce a different run_id
    gen = dict(task["generator"])
    from chemtrans.measure import derive_seed
    gen.setdefault("seed", derive_seed(task["seed_key"] + "|gen"))
    ident = {k: v for k, v in task.items() if k not in ("seed_key",
                                                       "name")}
    ident["generator"] = gen
    ident["audit_round"] = "r999"
    assert run_id(rec["net_hash"], ident,
                  task["seed_key"]) != rec["run_id"]
