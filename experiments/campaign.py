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

from chemtrans.measure import derive_seed, measure_run
from chemtrans.provenance import (
    code_version,
    environment,
    sha256_json,
    strict_jsonable,
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
            # the ACTUAL stochastic detection thresholds measure_run
            # applies (task.stochastic defaults at runtime), not
            # boilerplate — sweeps may override per task
            "theta_organized": float(
                (protocol.get("task", {}).get("stochastic") or {})
                .get("theta_organized", 0.1)),
            "mass_floor": float(
                (protocol.get("task", {}).get("stochastic") or {})
                .get("mass_floor", 0.01)),
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
    report = {"protocol": protocol["name"], "checks": []}
    tasks = expand_tasks(protocol)
    # dedupe on the RESOLVED generator spec, not seed_key: a
    # generator-param sweep (generator.delta_mu, generator.p_cat, ...)
    # produces different networks under one seed_key, and every one of
    # them needs conservation + thermodynamic validation
    nets = {}
    for e in tasks:
        sk = e["seed_key"]
        from chemtrans.measure import build_network
        gen = dict(e["task"]["generator"])
        gen.setdefault("seed", derive_seed(sk + "|gen"))
        gkey = sha256_json(gen)
        if gkey in nets:
            continue
        net, _ = build_network(gen)
        nets[gkey] = (sk, net, gen)
    ok = True
    for gkey, (sk, net, gen) in nets.items():
        check = {"seed_key": sk, "net_hash": net.hash,
                 "generator": {k: v for k, v in gen.items()
                               if k != "seed"}}
        try:
            db = net.check_detailed_balance()
            check["thermo_consistent"] = db["ok"]
            if not db["ok"]:
                ok = False
                check["bad_channels"] = [
                    r for r in db["channels"]
                    if r.get("rel_err", 0.0) > 1e-6
                    or (r.get("error") == "irreversible_channel"
                        and not r.get("declared_driven"))]
        except ValueError:
            check["thermo_consistent"] = None  # no g: not checkable
        if net.comp is not None:
            lc = net.check_letter_conservation()
            check["letter_conservation"] = lc["ok"]
            if lc.get("declared_open"):
                check["declared_open_system"] = True
            # sparse-safe: (S.T @ comp.T).T avoids a dense@sparse product
            S = net.stoichiometry()
            res = (np.asarray(S.T @ net.comp.T).T
                   if net.use_sparse else net.comp @ S)
            check["max_letter_drift"] = float(np.abs(res).max())
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
    # must hash the same task identity as measure_run: every
    # analysis-affecting key, not just gen/reactor/eps/stochastic
    task_identity = {k: v for k, v in task_entry["task"].items()
                     if k not in ("seed_key", "name")}
    task_identity["generator"] = gen
    return run_id(net.hash, task_identity, task_entry["seed_key"])


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
        if op == "m_poly":
            # nonfood monomer-mass concentration — the physical order
            # parameter that does not saturate along organized branches
            return eq["m_poly"] if eq else None
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

    # transition-shape model comparison: fit per (network family,
    # secondary level) — pooling staggered onsets across networks or
    # drive levels fabricates a smooth sigmoid out of sharp transitions
    if control:
        from chemtrans.analysis.critical_points import fit_all_models
        sec_ax = manifest["protocol"].get("secondary_axis")
        sec_param = (sec_ax if isinstance(sec_ax, str)
                     else (sec_ax or {}).get("param"))
        groups: dict = {}
        for r in records:
            lam_v = r.get("sweep", {}).get(control)
            y_v = extract_y(r)
            if lam_v is None or y_v is None:
                continue
            fam = f"{r.get('seed_key')}|{r['net_hash']}"
            sec = (float(r["sweep"][sec_param])
                   if sec_param and sec_param in r.get("sweep", {})
                   else None)
            groups.setdefault((fam, sec), {}).setdefault(
                float(lam_v), []).append(float(y_v))
        fits_by_group = {}
        for (fam, sec), lam_map in sorted(groups.items()):
            gname = fam if sec is None else f"{fam}|{sec_param}={sec}"
            xs = sorted(lam_map)
            ys_mean = [float(np.mean(lam_map[g])) for g in xs]
            if len(xs) < 4:
                fits_by_group[gname] = {
                    "status": "insufficient_points",
                    "n_lambda": len(xs)}
                continue
            # fit in CONTROL units: the control axis may negate the
            # raw sweep (e.g. dilution washes out polymer, so order
            # grows with -dilution). Onset models (hinge/fold/itrans)
            # are parameterized positive-on; fitting raw decreasing
            # data hands a spurious win to the direction-agnostic
            # sigmoid.
            fit = fit_all_models(sign * np.array(xs), np.array(ys_mean))
            tfit = {k: v for k, v in fit.items() if k != "_fits"}
            params = dict(tfit.get("params") or {})
            for key in ("lc", "lf"):
                if key in params:
                    params[key + "_raw"] = sign * params[key]
            tfit["params"] = params
            tfit["control_sign"] = sign
            # failed fits carry aicc=inf -> not strict-JSON; null them
            for kk in ("aicc", "rss"):
                if isinstance(tfit.get(kk), dict):
                    tfit[kk] = {m: (float(v) if np.isfinite(v) else None)
                                for m, v in tfit[kk].items()}
            fits_by_group[gname] = tfit
        sel_counts: dict = {}
        for g in fits_by_group.values():
            sel = g.get("selected")
            if sel:
                sel_counts[sel] = sel_counts.get(sel, 0) + 1
        summary["transition_fit"] = {
            "mode": "per_family_per_secondary",
            "grouping": ("seed_key|net_hash"
                         + (f" x {sec_param}" if sec_param else "")),
            "consensus_model": (max(sel_counts, key=sel_counts.get)
                                if sel_counts else None),
            "selection_counts": sel_counts,
            "fits": fits_by_group}

        # critical slowing down: |Re λ1| ~ |λ-λc|^ν per family on the
        # stable branch — wires the recorded leading_eig_re into an
        # actual exponent estimate instead of eyeballed point values.
        # One-sided (lam > lam_c): the post-critical stable branch.
        from chemtrans.observables import spectral as _spectral
        csd_groups: dict = {}
        for r in records:
            eq = r.get("equilibrium") or {}
            re1 = eq.get("leading_eig_re")
            lam_v = r.get("sweep", {}).get(control)
            if re1 is None or lam_v is None:
                continue
            fam = f"{r.get('seed_key')}|{r['net_hash']}"
            sec = (float(r["sweep"][sec_param])
                   if sec_param and sec_param in r.get("sweep", {})
                   else None)
            csd_groups.setdefault((fam, sec), []).append(
                (float(lam_v), float(re1)))
        csd: dict = {}
        for (fam, sec), pairs in sorted(csd_groups.items()):
            gname = (fam if sec is None
                     else f"{fam}|{sec_param}={sec}")
            gfit = fits_by_group.get(gname) or {}
            lam_c = (gfit.get("params") or {}).get("lc_raw")
            lams = np.array([p[0] for p in pairs])
            re1 = np.array([p[1] for p in pairs])
            if lam_c is None:
                csd[gname] = {"status": "no_lambda_c",
                              "n": len(pairs)}
                continue
            csd[gname] = _spectral.fit_slowing_down(
                lams, re1, float(lam_c))
            csd[gname]["lambda_c"] = float(lam_c)
        if csd:
            summary["slowing_down"] = csd

    # persistence scaling: Kaplan-Meier per volume + extinction-scale
    # fit across volumes — replaces eyeballing raw means of n=3-8
    # replicates and treats truncated runs as right-censored, not
    # events. Only emitted when stochastic records carry extinction
    # data.
    vol_groups: dict = {}
    for r in records:
        st = r.get("stochastic") or {}
        vol = st.get("volume")
        if vol is None:
            continue
        for rep in st.get("replicates") or []:
            t_ext = rep.get("t_extinction")
            t_cen = rep.get("censor_time")
            if t_ext is None and t_cen is None:
                continue
            vol_groups.setdefault(float(vol), []).append(
                (float(t_ext) if t_ext is not None else float(t_cen),
                 t_ext is not None))
    if vol_groups:
        from chemtrans.observables import persistence as _pers
        surv: dict = {}
        tau_v, tau_t = [], []
        for v, obs in sorted(vol_groups.items()):
            tab = _pers.survival_table(obs)
            t_k, S, _se = _pers.kaplan_meier(tab["times"], tab["events"])
            # KM median survival — the first crossing S < 0.5; None
            # when censoring keeps the curve above 0.5 throughout
            med = (float(t_k[int(np.argmax(S < 0.5))])
                   if np.any(S < 0.5) else None)
            surv[str(v)] = {
                "n": tab["n"], "n_events": tab["n_events"],
                "n_censored": tab["n_censored"],
                "median_survival": med,
                "rmst": _pers.rmst(tab["times"], tab["events"],
                                   float(tab["max_observed_time"]))}
            if med is not None:
                tau_v.append(v)
                tau_t.append(med)
        pers = {"per_volume": surv}
        if len(tau_v) >= 3:
            pers["extinction_scale"] = _pers.extinction_scale_fit(
                np.array(tau_v), np.array(tau_t))
        else:
            pers["extinction_scale"] = {
                "status": "insufficient_volumes",
                "n_volumes_with_median": len(tau_v)}
        summary["persistence"] = pers

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
        secondary = manifest["protocol"].get("secondary_axis")
        if secondary:
            if isinstance(secondary, str):
                secondary = {"param": secondary}
            lc_boot = int(manifest["protocol"].get("lc_boot", 0))
            summary["cocriticality_2d"] = cc.coincidence_2d(
                records, secondary["param"],
                n_boot=int(args.boot), lc_boot=lc_boot)
            # a pooled 1D report would flatten every secondary level
            # into each family's lambda-series — mixing heterogeneous
            # kinetic regimes. Per-level 1D reports live inside
            # cocriticality_2d["per_level"].
            summary["cocriticality"] = {
                "status": "suppressed",
                "reason": "secondary_axis declared; per-level reports "
                          "in cocriticality_2d.per_level"}
        else:
            summary["cocriticality"] = cc.coincidence_report(
                cc.group_by_network(records),
                n_boot=int(args.boot),
                lc_boot=int(manifest["protocol"].get("lc_boot", 0)))

    (cdir / "summary.json").write_text(
        json.dumps(strict_jsonable(summary), indent=2, sort_keys=True))
    print(json.dumps({k: summary[k] for k in
                      ("campaign_dir", "n_records", "status_counts",
                       "transition_fit")}, indent=2))


def cmd_zoom(args):
    """Scout->Zoom: read a campaign summary, bracket the detected
    critical region, and emit a refined-grid protocol.

    Uses per-family λc estimates in summary['cocriticality'] for ONE
    target observable (--observable; default: the protocol's declared
    order parameter / primary observable). Pooling λc's across
    observables takes the union envelope of unrelated changepoints,
    which widens the window across the whole original domain instead
    of zooming on a single transition. Per-family heterogeneity still
    widens the window — correct for a paired re-measurement.
    """
    mdir = Path(args.manifest).parent
    manifest = json.loads(Path(args.manifest).read_text())
    summary = json.loads((mdir / "summary.json").read_text())
    proto = manifest["protocol"]
    ax = proto.get("control_axis") or {}
    param = ax.get("param") or proto["sweeps"][0]["param"]
    direction = ax.get("direction") or "identity"

    fam_ests = (summary.get("cocriticality") or {}).get(
        "family_estimates") or {}
    if not fam_ests:
        # 2D campaign: the pooled 1D block is suppressed — union the
        # per-level family estimates instead. Each level's estimates
        # were computed on unflattened data; unioning only widens the
        # refinement envelope (conservative for bracketing).
        for lvl, rep in ((summary.get("cocriticality_2d") or {})
                         .get("per_level") or {}).items():
            for obs, fams in (rep.get("family_estimates") or {}).items():
                for fam, e in fams.items():
                    fam_ests.setdefault(obs, {})[
                        f"{fam}@k={lvl}"] = e
    # map record-path order parameters to the life-vector observable
    # that brackets the same transition
    op_map = {"equilibrium.m_poly": "O_maint", "mass_raf": "O_maint",
              "equilibrium.assembly_signature.mass_signature":
              "O_assembly",
              "equilibrium.assembly_signature.mass_signature_abs":
              "O_assembly"}
    target_obs = (args.observable
                  or proto.get("primary_observable")
                  or op_map.get(proto.get("order_parameter")))
    if target_obs and target_obs.upper() == "ALL":
        target_obs = None
    if target_obs:
        src = {target_obs: fam_ests[target_obs]} \
            if target_obs in fam_ests else {}
        if not src:
            raise SystemExit(
                f"zoom: observable '{target_obs}' has no family λc "
                f"estimates (available: {sorted(fam_ests)})")
    else:
        # no declared primary observable — fall back to the envelope
        # of ALL observables, flagged as such in the output metadata
        src = fam_ests
    lcs = [e["lc"] for obs in src.values()
           for e in obs.values()
           if e.get("lc") is not None and not e.get("extrapolated")]
    if not lcs:
        raise SystemExit("zoom: no in-range family λc estimates — "
                         "nothing to refine")
    pad = float(args.pad)
    lo, hi = min(lcs) - pad, max(lcs) + pad
    n = int(args.n_points)
    lam_zoom = np.linspace(lo, hi, n)
    raw = (-lam_zoom if direction == "negate" else lam_zoom)
    raw = np.sort(raw)
    # clip to the original sweep's domain — sampling outside it is
    # unconstrained extrapolation, and for physical axes like
    # dilution can be outright invalid (D < 0)
    orig = [float(v) for s in proto["sweeps"] if s["param"] == param
            for v in s["values"]]
    if orig:
        raw = raw[(raw >= min(orig) - 1e-12)
                  & (raw <= max(orig) + 1e-12)]
    if len(raw) < 3:
        raise SystemExit("zoom: refined grid <3 points inside the "
                         "original sweep domain")

    newp = json.loads(json.dumps(proto))
    newp["name"] = proto["name"] + "_zoom"
    for s in newp["sweeps"]:
        if s["param"] == param:
            s["values"] = [float(v) for v in raw]
    newp["zoom_of"] = {"campaign": str(mdir),
                       "lam_window": [float(lo), float(hi)],
                       "observable": target_obs or "ALL_POOLED",
                       "lc_inputs": [float(v) for v in sorted(lcs)]}
    out = Path(args.out or f"{proto['name']}_zoom.json")
    out.write_text(json.dumps(newp, indent=2, sort_keys=True))
    print(json.dumps({"protocol": str(out),
                      "control": param,
                      "observable": target_obs or "ALL_POOLED",
                      "lam_window": [float(lo), float(hi)],
                      "n_lc_inputs": len(lcs),
                      "values": [float(v) for v in raw]}, indent=2))


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
    p = sub.add_parser("zoom")
    p.add_argument("--manifest", required=True)
    p.add_argument("--n-points", type=int, default=9)
    p.add_argument("--pad", type=float, default=0.25)
    p.add_argument("--observable",
                   help="life-vector observable whose family λc's "
                        "bracket the zoom window (default: the "
                        "protocol's order_parameter / "
                        "primary_observable; 'ALL' pools every "
                        "observable — the old envelope behavior)")
    p.add_argument("--out")
    p.set_defaults(fn=cmd_zoom)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
