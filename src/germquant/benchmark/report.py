"""Writing the benchmark outputs: per-gonad CSV, per-nucleus CSV and the pooled JSON."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

PER_GONAD = "per_gonad.csv"
PER_NUCLEUS = "per_nucleus.csv"
POOLED = "pooled.json"


def jsonable(obj):
    """Plain JSON types only: numpy scalars and arrays become Python values, NaN and infinities
    become null (so the file is standard JSON, readable outside Python), Paths become strings."""
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [jsonable(v) for v in obj.tolist()]
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        f = float(obj)
        return f if math.isfinite(f) else None
    if isinstance(obj, Path):
        return str(obj)
    return obj


def write_outputs(out_dir: str | Path, per_gonad: pd.DataFrame, per_nucleus: pd.DataFrame,
                  pooled: dict, extra_tables: dict | None = None) -> dict:
    """Write the three contract files (and any extra named tables) into ``out_dir``; returns the paths."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {"per_gonad": out_dir / PER_GONAD, "per_nucleus": out_dir / PER_NUCLEUS,
             "pooled": out_dir / POOLED}
    per_gonad.to_csv(paths["per_gonad"], index=False)
    per_nucleus.to_csv(paths["per_nucleus"], index=False)
    with open(paths["pooled"], "w", encoding="utf-8") as fh:
        json.dump(jsonable(pooled), fh, indent=1)
    for name, table in (extra_tables or {}).items():
        p = out_dir / f"{name}.csv"
        table.to_csv(p, index=False)
        paths[name] = p
    return paths
