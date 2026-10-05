"""Provenance: protocol/code hashing, run IDs, manifest helpers.

Rules:
- run records are strict JSON: no NaN/Infinity; use null + reason field.
- run_id and campaign manifest hash are content-derived; incompatible
  hashes are rejected by the campaign runner, not silently overwritten.
- code_version = git commit SHA if available, else a source-tree hash.
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np


def strict_jsonable(obj):
    """Convert numpy scalars/arrays; reject non-finite floats."""
    if isinstance(obj, dict):
        return {str(k): strict_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [strict_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return strict_jsonable(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        f = float(obj)
        if not math.isfinite(f):
            raise ValueError(f"non-finite value {f!r} cannot be "
                             "serialized to strict JSON")
        return f
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if obj is None or isinstance(obj, (str, bool, int)):
        return obj
    raise TypeError(f"not JSON-serializable: {type(obj)!r}")


def sha256_json(obj) -> str:
    return hashlib.sha256(
        json.dumps(strict_jsonable(obj), sort_keys=True,
                   separators=(",", ":")).encode()
    ).hexdigest()


def code_version(root: Path | None = None) -> str:
    """git commit SHA (+ '-dirty' if worktree differs); else dir hash."""
    try:
        sha = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=root, check=True).stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain"], capture_output=True,
            text=True, cwd=root, check=True).stdout.strip()
        return sha + ("-dirty" if dirty else "")
    except Exception:
        h = hashlib.sha256()
        src = Path(root or ".") / "src"
        if src.exists():
            for f in sorted(src.rglob("*.py")):
                h.update(f.read_bytes())
        return "nogit-" + h.hexdigest()[:16]


def environment() -> dict:
    import scipy
    import networkx
    import pandas
    return {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "networkx": networkx.__version__,
        "pandas": pandas.__version__,
    }


def run_id(net_hash: str, params: dict, seed_key: str) -> str:
    """Deterministic run identity: hash of (network, params, seed key)."""
    return sha256_json({"net": net_hash, "params": params,
                        "seed": seed_key})[:24]
