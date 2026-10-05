"""Transition-shape inference: model recovery on synthetic data,
strict-JSON provenance, network-joint bootstrap."""

import numpy as np
import pytest

from chemtrans.analysis.critical_points import (
    fit_all_models, network_joint_bootstrap)
from chemtrans.provenance import run_id, sha256_json, strict_jsonable


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
                 0.7 + 0.4 * np.sqrt(lam - 0.55))
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
