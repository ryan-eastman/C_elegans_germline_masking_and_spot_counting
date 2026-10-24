"""Digests behind the sc_straighten v2 compatibility test (tests/test_sc_straighten_v2_compat.py) and its golden
(tests/golden/sc_straighten/v1_compat.json, made once with the v1 package code of commit 5cb3cb6).

Tables are digested over the v1 columns only (the declared columns of v1, plus the run's shared metadata for the
files the runner writes, without sc_straighten_version and sc_straighten_hash), formatted as the writer formats
them (``%.6g``), so columns added in v2 and the version and hash columns are not compared. Strips are digested
as the uint16 arrays the TIFF writer stores. Scene pins hash the image, the labels and the truth of a registered
synthetic family (the generator of protocol 2.3 stays identical in v2).
"""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd

TABLE_KEYS = ("strands", "samples", "gaps", "foci", "nuclei")
FAMILIES = ("iso", "dense", "gaps", "foci", "pairs", "noise", "snr")
SYNTH_CASES = tuple((f, "entry2", s) for f in FAMILIES for s in range(2000, 2005)) + \
    tuple((f, "registered", s) for f in ("dense", "snr") for s in range(2000, 2005))
SCENE_PINS = tuple((f, "entry2", s) for f in FAMILIES for s in (2000, 2001)) + \
    tuple((f, "registered", s) for f in ("dense", "snr") for s in (2000, 2001))
DROP_META = ("sc_straighten_version", "sc_straighten_hash")
# shared metadata that changes from run to run (where and when the germquant run of the fixture was made)
VOLATILE_META = ("file_path", "pipeline_version", "git_sha", "config_hash", "run_timestamp")
GONAD_IMAGE = "20251105_N2_nohs_HERM _.nd2"
_V1_COLUMNS: dict | None = None


def v1_columns() -> dict:
    """The declared columns of the tables (at v1 this is schema.TABLES; the golden stores them)."""
    from germquant.sc.straighten import schema as S

    return {k: list(v) for k, v in S.TABLES.items()}


def golden_path() -> Path:
    return Path(__file__).resolve().parent / "golden" / "sc_straighten" / "v1_compat.json"


def _columns(cols: dict | None) -> dict:
    if cols is not None:
        return cols
    p = golden_path()
    if p.is_file():
        return json.loads(p.read_text(encoding="ascii"))["columns"]
    return v1_columns()


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def table_digest(df: pd.DataFrame, cols: list[str]) -> str:
    d = df.copy()
    for c in cols:
        if c not in d.columns:
            d[c] = pd.NA
    buf = io.StringIO()
    d[list(cols)].to_csv(buf, index=False, float_format="%.6g", lineterminator="\n")
    return sha(buf.getvalue().encode())


def strips_digest(arr: np.ndarray) -> str:
    a = np.nan_to_num(np.asarray(arr, np.float64), nan=0.0)
    a = np.clip(np.rint(a), 0, 65535).astype(np.uint16)
    return sha(a.tobytes() + str(a.shape).encode())


def run_scene(sc, sets=()):
    from germquant.sc.straighten import params as PRM
    from germquant.sc.straighten.pipeline import straighten_image

    P = PRM.resolve(None, list(sets))
    roles = {k: v for k, v in sc.channels.items() if k != "central_element"}
    spots = {"spots": sc.spots} if len(sc.spots) else None
    return straighten_image(sc.channels["central_element"], sc.labels, sc.spacing, sc.nuclei, params=P, roles=roles,
                            spot_tables=spots, image_id=f"{sc.family}_{sc.seed}")


def result_digests(res, cols: dict | None = None) -> dict:
    cols = _columns(cols)
    out = {k: table_digest(getattr(res, k), cols[k]) for k in TABLE_KEYS}
    out["strips_flat"] = strips_digest(res.strips_flat)
    out["strips_axial"] = strips_digest(res.strips_axial)
    out["strip_index"] = sha(json.dumps([[int(x) for x in r] for r in res.strip_index]).encode())
    out["n_strands"] = int(len(res.strands))
    return out


def synthetic_digests(fam: str, gen: str, seed: int, sets=(), cols: dict | None = None) -> dict:
    from germquant.sc.straighten import synth

    sc = synth.make_scene(fam, seed, generator=gen)
    return result_digests(run_scene(sc, sets), cols)


def scene_digests(fam: str, gen: str, seed: int) -> dict:
    from germquant.sc.straighten import synth

    sc = synth.make_scene(fam, seed, generator=gen)
    truth = []
    for r in sc.ribbons:
        truth.append({"id": int(r.ribbon_id), "nucleus": int(r.nucleus_id),
                      "points": sha(np.round(np.asarray(r.points, float), 9).tobytes()),
                      "gaps": [[round(float(a), 9), round(float(b), 9)] for a, b in r.gaps],
                      "foci_s": [round(float(s), 9) for s in r.foci_s], "peak": round(float(r.peak), 6)})
    buf = io.StringIO()
    sc.spots.to_csv(buf, index=False, float_format="%.9g", lineterminator="\n")
    return {"channels": {k: sha(np.ascontiguousarray(v).tobytes()) for k, v in sorted(sc.channels.items())},
            "labels": sha(np.ascontiguousarray(sc.labels).tobytes()), "shape": list(sc.labels.shape),
            "truth": sha(json.dumps(truth, sort_keys=True).encode()), "spots": sha(buf.getvalue().encode())}


def finished_gonad(root: Path, image: str = GONAD_IMAGE):
    """process_image on the end-to-end synthetic gonad (tests/test_end_to_end.py) under ``image``, as the CLI
    tests make it: (results root, image id, stack)."""
    import pytest
    import test_end_to_end as e2e

    from germquant import pipeline
    from germquant.config import load_config

    mp = pytest.MonkeyPatch()
    stack = e2e._synthetic_stack()
    meta = {"sizes": {"C": 3, "Z": 20, "Y": 80, "X": 320}, "spacing": e2e.SPACING, "spacing_ok": True,
            "channel_names": e2e.CHANNELS, "dtype": "float32", "is_2d": False}
    mp.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: meta)
    mp.setattr(pipeline, "read_stack", lambda *a, **k: stack)
    cfg = load_config(str(Path(__file__).resolve().parents[1] / "config" / "config.yaml"))
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    res = pipeline.process_image(image, cfg, root / "img", prov=None)
    mp.undo()
    return root, res["image_id"], stack


def gonad_digests(tmp: Path, sets=(), cols: dict | None = None) -> dict:
    from germquant.sc.straighten import params as PRM
    from germquant.sc.straighten import runner
    from germquant.sc.straighten import schema as S

    cols = _columns(cols)
    root, iid, stack = finished_gonad(Path(tmp) / "results")

    def read(path, roles, xy_stride=1, z_range=None):
        return {r: stack.data[int(i)].copy() for r, i in roles.items()}

    oc = runner.run_results(root, params=PRM.resolve(None, list(sets)), reader=read, echo=lambda *_: None)
    assert [o.status for o in oc] == ["done"], [o.line for o in oc]
    d = root / "img"
    out = {"image_id": iid}
    for k in TABLE_KEYS:
        df = pd.read_csv(d / f"{iid}{S.SUFFIXES[k]}", keep_default_na=False, dtype=str)
        meta = [c for c in df.columns if c not in S.TABLES[k] and c not in DROP_META]
        want = [c for c in cols[k] if c in df.columns] + [c for c in meta if c not in cols[k] and c in _v1_meta(df, k)]
        out[k] = sha(df[want].to_csv(index=False, lineterminator="\n").encode())
        out[k + "_columns"] = want
    return out


def _v1_meta(df: pd.DataFrame, key: str) -> set:
    """The shared metadata columns of a written table (the run's SHARED_META block)."""
    from germquant.schema import SHARED_META

    return (set(SHARED_META) - set(VOLATILE_META)) if key != "samples" else set()
