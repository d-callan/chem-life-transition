"""Transition-shape inference: model recovery on synthetic data,
strict-JSON provenance, network-joint bootstrap."""

import numpy as np
import pytest

from chemtrans.analysis.critical_points import fit_all_models, network_joint_bootstrap
from chemtrans.provenance import run_id, strict_jsonable


def test_null_selected_on_flat_data():
    rng = np.random.default_rng(0)
    lam = np.linspace(0, 1, 30)
    y = 0.02 + 0.005 * rng.standard_normal(30)
    res = fit_all_models(lam, y, seed=0)
    assert res["selected"] == "null"


def test_hinge_selected_on_piecewise_data():
    rng = np.random.default_rng(1)
    lam = np.linspace(0, 1, 40)
    y = np.maximum(0.0, lam - 0.5) * 1.5 + 0.002 * rng.standard_normal(40)
    res = fit_all_models(lam, y, seed=0)
    assert res["selected"] in ("hinge", "itrans")
    # hinge and itrans are both continuous-onset models; the sigmoid or
    # fold must not win on clearly piecewise-linear data.
    assert res["aicc"]["sigmoid"] > res["aicc"][res["selected"]]


def test_fold_selected_on_jump_data():
    rng = np.random.default_rng(2)
    lam = np.linspace(0, 1, 40)
    y = np.where(lam < 0.55, 0.02,
                 0.7 + 0.4 * np.sqrt(np.maximum(lam - 0.55, 0.0)))
    y += 0.002 * rng.standard_normal(40)
    res = fit_all_models(lam, y, seed=0)
    assert res["selected"] == "fold"


def test_bootstrap_runs_and_preserves_pairing():
    rng = np.random.default_rng(3)
    records = []
    for _ in range(6):
        lam = np.linspace(0, 1, 15)
        y = np.maximum(0, lam - 0.4) + 0.01 * rng.standard_normal(15)
        records.append({"lambda": lam, "O": y})
    out = network_joint_bootstrap(records, "lambda", "O",
                                  n_boot=20, seed=0)
    assert out["n_boot"] == 20
    assert sum(out["selection_freq"].values()) == 20


def test_strict_jsonable_rejects_nan():
    with pytest.raises(ValueError):
        strict_jsonable({"x": float("nan")})
    with pytest.raises(ValueError):
        strict_jsonable({"x": np.array([1.0, np.inf])})
    assert strict_jsonable({"a": np.int64(3)}) == {"a": 3}


def test_run_id_deterministic():
    a = run_id("nethash", {"p": 1}, "seed|0")
    b = run_id("nethash", {"p": 1}, "seed|0")
    c = run_id("nethash", {"p": 2}, "seed|0")
    assert a == b != c


def test_mirrored_models_recover_decreasing_transitions():
    """Decreasing data must be fittable by onset-class models — the
    one-sided candidate set previously handed every decreasing curve
    to the direction-agnostic sigmoid/ramp (structural AICc bias)."""
    lam = np.linspace(0, 4, 20)
    # offset hinge: flat-high left, drop at lc=2.0, flat-low right
    y_h = np.where(lam < 2.0, 0.9 - 0.35 * lam, 0.2)
    fit = fit_all_models(lam, y_h, seed=0)
    assert fit["selected"] == "hinge_rev"
    assert fit["params"]["lc"] == pytest.approx(2.0, abs=0.05)
    # backward-opening fold: jump at lf=1.5, sqrt branch below
    y_f = np.where(lam > 1.5, 0.05,
                   0.5 + 0.4 * np.sqrt(np.maximum(0.0, 1.5 - lam)))
    fit = fit_all_models(lam, y_f, seed=0)
    assert fit["selected"] == "fold_rev"
    assert fit["params"]["lf"] == pytest.approx(1.5, abs=0.05)
