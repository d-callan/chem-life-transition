#!/usr/bin/env python3
"""Campaign runner: validate / plan / run / summarize.

A protocol JSON specifies the fixed generator + reactor spec, a sweep
(dotted param path -> values), network replicates (seeds), and an
output directory. Everything downstream is content-addressed:

- campaign dir: <out_dir>/<name>-<protocol_hash[:12]>
- manifest.json: protocol, hashes, environment, decisions — immutable;
  a rerun with an incompatible protocol hash is rejected.
- shard-*.jsonl: append-only run records; a run whose run_id already
  exists is skipped (resume-safe, reorder-safe — seeds derive from
  seed_key strings, not task order).
- summary.json: aggregated by sweep value; strict JSON (no NaN/Inf).

Dotted sweep paths index into the task spec, e.g. "generator.p_cat",
"reactor.dilution", "stochastic.volume".
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from itertools import product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import numpy as np

from chemtrans.measure import measure_run, derive_seed
from chemtrans.provenance import (
    code_version, environment, sha256_json, strict_jsonable,
)


def _set_path(spec: dict, path: str, value) -> dict:
    spec = json.loads(json.dumps(spec))  # deep copy
    keys = path.split(".")
    d = spec
    for k in keys[:-1]:
        d = d.setdefault(k, {})
    d[keys[-1]] = value
    return spec


def load_protocol(path: Path) -> dict:
    p = json.loads(path.read_text())
    for k in ("name", "experiment", "task", "sweeps", "out_dir"):
        if k not in p:
            raise ValueError(f"protocol missing {k!r}")
    return p


def expand_tasks(protocol: dict) -> list[dict]:
    """Cartesian product of sweep paths/values x network seeds."""
    sweeps = protocol["sweeps"]
    paths = [s["param"] for s in sweeps]
    values = [s["values"] for s in sweeps]
    # network replicates: explicit seed_keys (carry-forward from a prior
    # campaign) or n_networks generated from the protocol name
    seed_keys = protocol.get("seed_keys")
    if seed_keys is None:
        seed_keys = [f"{protocol['name']}|net{rep}"
                     for rep in range(int(protocol.get("n_networks", 1)))]
    tasks = []
    for rep, seed_key in enumerate(seed_keys):
        for combo in product(*values):
            spec = json.loads(json.dumps(protocol["task"]))
            lam = {}
            for path, v in zip(paths, combo):
                spec = _set_path(spec, path, v)
                lam[path] = v
            spec["seed_key"] = seed_key   # inside spec: drives all seeds
            tasks.append({"task": spec, "sweep": lam,
                          "seed_key": seed_key, "net_replicate": rep})
    return tasks


def campaign_dir(protocol: dict) -> Path:
    ph = sha256_json(protocol)[:12]
    return Path(protocol["out_dir"]) / f"{protocol['name']}-{ph}"


def write_manifest(protocol: dict, cdir: Path):
    manifest = {
        "protocol_hash": sha256_json(protocol),
        "protocol": protocol,
        "code_version": code_version(
            Path(__file__).resolve().parent.parent),
        "environment": environment(),
        "created": time.time(),
        "decisions": {
            "order_parameter": "mass_raf",
            "theta_hi": 0.6, "theta_lo": 0.4,
            "init_policy": "round",
            "primary_experiment": protocol["experiment"],
        },
    }
    (cdir / "manifest.json").write_text(
        json.dumps(strict_jsonable(manifest), indent=2, sort_keys=True))
    return manifest


def check_manifest(protocol: dict, cdir: Path):
    mpath = cdir / "manifest.json"
    if not mpath.exists():
        return None
    manifest = json.loads(mpath.read_text())
    if manifest["protocol_hash"] != sha256_json(protocol):
        raise SystemExit(
            f"manifest protocol_hash mismatch in {cdir}: refusing to "
            "mix outputs from different protocols")
    return manifest


def _task_record(task_entry: dict) -> dict:
    rec = measure_run(task_entry["task"])
    rec["sweep"] = task_entry["sweep"]
    rec["net_replicate"] = task_entry["net_replicate"]
    # measure_run emits family_id = latent topology id (seed + pool-
    # defining params): dynamical AND generator-threshold sweeps stay
    # paired; pool changes (alphabet_size/max_len/ligation) split it.
    rec.setdefault("family_id",
                   f"{task_entry['seed_key']}|{rec['net_hash']}")
    return rec


def cmd_validate(args):
    """E0-style validation: protocol sanity + thermodynamic/null checks."""
    protocol = load_protocol(Path(args.protocol))
    from chemtrans.networks import core as core_mod
    report = {"protocol": protocol["name"], "checks": []}
    tasks = expand_tasks(protocol)
    nets = {}
    for e in tasks:
        sk = e["seed_key"]
        if sk in nets:
            continue
        from chemtrans.measure import build_network
        gen = dict(e["task"]["generator"])
        gen.setdefault("seed", derive_seed(sk + "|gen"))
        net, _ = build_network(gen)
        nets[sk] = net
    ok = True
    for sk, net in nets.items():
        check = {"seed_key": sk, "net_hash": net.hash}
        try:
            db = net.check_detailed_balance()
            check["thermo_consistent"] = db["ok"]
            if not db["ok"]:
                ok = False
                check["bad_channels"] = [
                    r for r in db["channels"] if r["rel_err"] > 1e-6]
        except ValueError:
            check["thermo_consistent"] = None  # no g: not checkable
        if net.comp is not None:
            lc = net.check_letter_conservation()
            check["letter_conservation"] = lc["ok"]
            check["max_letter_drift"] = float(
                np.abs(net.comp @ net.stoichiometry()).max())
            if not lc["ok"]:
                ok = False
        report["checks"].append(check)
    report["status"] = "PASS" if ok else "FAIL"
    report["n_networks"] = len(nets)
    print(json.dumps(strict_jsonable(report), indent=2))
    sys.exit(0 if ok else 1)


def cmd_plan(args):
    protocol = load_protocol(Path(args.protocol))
    tasks = expand_tasks(protocol)
    cdir = campaign_dir(protocol)
    plan = {
        "campaign_dir": str(cdir),
        "n_tasks": len(tasks),
        "sweeps": protocol["sweeps"],
        "n_networks": protocol.get("n_networks", 1),
        "task_hashes": [sha256_json(t)[:16] for t in tasks],
    }
    print(json.dumps(strict_jsonable(plan), indent=2))


def cmd_run(args):
    protocol = load_protocol(Path(args.protocol))
    cdir = campaign_dir(protocol)
    cdir.mkdir(parents=True, exist_ok=True)
    manifest = check_manifest(protocol, cdir)
    if manifest is None:
        manifest = write_manifest(protocol, cdir)

    tasks = expand_tasks(protocol)
    done = set()
    shard = cdir / "shard-0.jsonl"
    if shard.exists():
        for line in shard.read_text().splitlines():
            if line.strip():
                done.add(json.loads(line)["run_id"])
    todo = []
    for t in tasks:
        rid = _preview_id(t)
        if rid not in done:
            todo.append(t)
            done.add(rid)      # also dedupe within the plan itself

    from concurrent.futures import ProcessPoolExecutor
    workers = int(args.workers)
    n_written = 0
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for rec in ex.map(_task_record, todo, chunksize=1):
                _append(shard, rec)
                n_written += 1
    else:
        for t in todo:
            _append(shard, _task_record(t))
            n_written += 1
    print(json.dumps({"campaign": str(cdir), "planned": len(tasks),
                      "skipped_existing": len(tasks) - len(todo),
                      "written": n_written}))


def _preview_id(task_entry: dict) -> str:
    """Compute run_id without running (same inputs as measure_run)."""
    from chemtrans.measure import build_network
    from chemtrans.provenance import run_id
    gen = dict(task_entry["task"]["generator"])
    gen.setdefault("seed", derive_seed(task_entry["seed_key"] + "|gen"))
    net, _ = build_network(gen)
    params = {"gen": gen, "reactor": task_entry["task"].get("reactor", {}),
              "eps": task_entry["task"].get("eps", net.eps),
              "stochastic": task_entry["task"].get("stochastic")}
    return run_id(net.hash, params, task_entry["seed_key"])


def _append(shard: Path, rec: dict):
    with shard.open("a") as fh:
        fh.write(json.dumps(strict_jsonable(rec), sort_keys=True) + "\n")


def cmd_summarize(args):
    mpath = Path(args.manifest)
    manifest = json.loads(mpath.read_text())
    cdir = mpath.parent
    records = []
    for sh in sorted(cdir.glob("shard-*.jsonl")):
        for line in sh.read_text().splitlines():
            if line.strip():
                records.append(json.loads(line))
    if not records:
        raise SystemExit("no run records found")

    sweep_paths = sorted({k for r in records for k in
                          r.get("sweep", {})})
    # explicit control axis: {"param", "direction", "label"}; falls back
    # to first sorted sweep path + legacy control_negate flag
    ax = manifest["protocol"].get("control_axis") or {}
    control = ax.get("param") or (sweep_paths[0] if sweep_paths else None)
    direction = ax.get("direction") or (
        "negate" if manifest["protocol"].get("control_negate")
        else "pos")
    sign = -1.0 if direction == "negate" else 1.0
    lam_label = ax.get("label") or (
        control.split(".")[-1] if control else "lambda")

    # status counts
    statuses = {}
    for r in records:
        s = r["det"]["status"]
        statuses[s] = statuses.get(s, 0) + 1

    # order parameter vs control (certified equilibria only)
    op = manifest["protocol"].get("order_parameter", "mass_raf")

    def extract_y(rec):
        """Resolve order parameter: 'mass_raf' (from equilibrium),
        'x_eq:IDX', or a dotted path into the record."""
        eq = rec.get("equilibrium")
        if op == "mass_raf":
            return eq["mass_raf"] if eq else None
        if op.startswith("x_eq:"):
            return eq["x_eq"][int(op.split(":")[1])] if eq else None
        node = rec
        for part in op.split("."):
            node = node.get(part) if isinstance(node, dict) else None
            if node is None:
                return None
        return node

    pts = {}
    for r in records:
        if control is None:
            continue
        lam = r["sweep"][control]
        pts.setdefault(float(lam), []).append(
            {"run_id": r["run_id"], "net": r["net_hash"],
             "lam": lam,
             "y": extract_y(r),
             "status": r["det"]["status"]})

    summary = {
        "campaign_dir": str(cdir),
        "protocol": manifest["protocol"]["name"],
        "experiment": manifest["protocol"]["experiment"],
        "n_records": len(records),
        "status_counts": statuses,
        "control": control,
        "control_label": lam_label,
        "control_axis": {"param": control, "direction": direction,
                         "label": lam_label},
        "order_parameter": op,
        "points": {str(k): v for k, v in sorted(pts.items())},
    }

    # transition-shape model comparison when we have a usable curve
    if control:
        lams, ys = [], []
        for lam, rows in pts.items():
            vals = [p["y"] for p in rows if p["y"] is not None]
            if vals:
                lams.append(lam)
                ys.append(float(np.mean(vals)))
        if len(lams) >= 4:
            from chemtrans.analysis.critical_points import (
                fit_all_models)
            fit = fit_all_models(np.array(lams), np.array(ys))
            summary["transition_fit"] = {
                k: v for k, v in fit.items() if k != "_fits"}
        else:
            summary["transition_fit"] = {
                "status": "insufficient_points",
                "n_lambda": len(lams)}

    # co-criticality summary: per-observable λc + pairwise matrix
    if manifest["protocol"].get("summarize_mode") == "cocriticality":
        if control is None:
            raise SystemExit(
                "cocriticality mode requires a sweep control")
        from chemtrans.analysis import cocriticality as cc
        from chemtrans.observables.life_vector import life_vector
        for r in records:
            r["life_vector"] = life_vector(r)
            r["control"] = {"axis": control,
                            "lam": sign * float(r["sweep"][control]),
                            "label": lam_label,
                            "direction": direction}
            r.setdefault(
                "family_id", f"{r.get('seed_key')}|{r['net_hash']}")
        summary["cocriticality"] = cc.coincidence_report(
            cc.group_by_network(records),
            n_boot=int(args.boot))

    (cdir / "summary.json").write_text(
        json.dumps(strict_jsonable(summary), indent=2, sort_keys=True))
    print(json.dumps({k: summary[k] for k in
                      ("campaign_dir", "n_records", "status_counts",
                       "transition_fit")}, indent=2))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("validate")
    p.add_argument("--protocol", required=True)
    p.set_defaults(fn=cmd_validate)
    p = sub.add_parser("plan")
    p.add_argument("--protocol", required=True)
    p.set_defaults(fn=cmd_plan)
    p = sub.add_parser("run")
    p.add_argument("--protocol", required=True)
    p.add_argument("--workers", type=int, default=1)
    p.set_defaults(fn=cmd_run)
    p = sub.add_parser("summarize")
    p.add_argument("--manifest", required=True)
    p.add_argument("--boot", type=int, default=400)
    p.set_defaults(fn=cmd_summarize)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
