"""The single-plane decision and the spacing of a single-plane run (tranche 5 package P2b; design 4.2, 4.3, 4.9 and
4.10; docs/PLANAR.md; decisions D6 to D8 and D11 of the Train B plan).

* planar.decide(cfg, meta, stack) is None unless io.planar is auto or true. auto is planar only when the header
  has one plane (no Z axis, or Z of 1) and the stack read has one plane, so a 3D file cut to one plane with
  --z-range stays 3D; true refuses more than one plane before any stage runs, naming --z-range and
  segmentation.nuclei.project.
* The spacing chain of a single-plane run: io.voxel_size_um, then Stack.spacing_yx, then stack.spacing[1:] when
  spacing_ok, else (1, 1) with the flag voxel:UNREADABLE_xy_assumed_1um; the three io.voxel_size_um forms with
  xy_stride applied after them; voxel:override and voxel:override_differs_from_file only when the key is set, in
  3D runs too (D6).
* recipe_keys refuses a malformed io.planar or io.voxel_size_um before any stage.
* batch.is_done's safe default (D8): a planar marker without segment_hash_planar is compared with segment_hash and
  reprocesses with the reason "planar flag"; without any segment hash is_done behaves as today.
* io.planar together with segmentation.nuclei.project: planar wins, with the flag segment:project_ignored_planar
  (D7).
* Stack.spacing_yx (the last field, None by default) and the bioio reader's _spacing_yx, with a fake BioImage
  (no bioio needed); `germquant info` prints a single-plane line for a one-plane header and nothing new for a
  3D one.
"""
from __future__ import annotations

import dataclasses
import json
import re
from collections import namedtuple
from pathlib import Path

import numpy as np
import pandas as pd
import planar_helpers as H
import pytest

from germquant import batch, cli, pipeline, planar, provenance, recipe_keys, stages
from germquant.config import load_config
from germquant.io import bioio_reader, read_image, read_image_metadata
from germquant.io.nd2_reader import Stack

ROOT = Path(__file__).resolve().parents[1]
ND2 = f"{H.STEM}.nd2"


def _stack(nz=1, *, spacing=(H.DZ, H.DY, H.DX), spacing_ok=True, spacing_yx=None, stride=1):
    data = np.zeros((2, nz, 8, 8), np.uint16)
    return Stack(data=data, spacing=spacing, channel_names=list(H.CHANNELS), path=Path(ND2), downsample_xy=stride,
                 spacing_ok=spacing_ok, spacing_yx=spacing_yx)


def _meta(z=None):
    sizes = {"C": 2, "Y": 8, "X": 8} if z is None else {"C": 2, "Z": z, "Y": 8, "X": 8}
    return {"sizes": sizes, "spacing": (H.DZ, H.DY, H.DX), "spacing_ok": True, "channel_names": list(H.CHANNELS)}


def _cfg(**overrides):
    cfg = load_config(ROOT / "config" / "profiles" / "generic.yaml")
    for k, v in overrides.items():
        cfg.set(k, v)
    return cfg


# ---- the decision ------------------------------------------------------------------------------------------

def test_decide_false_auto_true():
    assert planar.decide(_cfg(), _meta(), _stack()) is None                     # the key absent: every existing run
    assert planar.decide(_cfg(**{"io.planar": False}), _meta(), _stack()) is None
    auto = planar.decide(_cfg(**{"io.planar": "auto"}), _meta(), _stack())
    assert auto == {"planar": True, "dy": H.DY, "dx": H.DX, "spacing_source": "file", "flags": ["read:planar"]}
    assert planar.decide(_cfg(**{"io.planar": True}), _meta(z=1), _stack())["planar"] is True
    assert planar.PLANAR_DZ == 1.0
    # auto: one plane in the header AND in the stack read
    assert planar.decide(_cfg(**{"io.planar": "auto"}), _meta(z=1), _stack())["planar"] is True
    three_d = planar.decide(_cfg(**{"io.planar": "auto"}), _meta(z=5), _stack(nz=5))
    assert three_d["planar"] is False and three_d["flags"] == []
    assert planar.decide(_cfg(**{"io.planar": "auto"}), {"sizes": {}}, _stack())["planar"] is False


def test_auto_keeps_a_cut_3d_file_3d(tmp_path, monkeypatch):
    """Under auto, a 3D file read as one plane (--z-range) stays a 3D run: no planar record, today's flags."""
    cfg = _cfg(**{"io.planar": "auto"})
    assert planar.decide(cfg, _meta(z=5), _stack(nz=1))["planar"] is False
    H.patch_reader(monkeypatch, H.stack_3d(), two_d=False)
    res = pipeline.process_image(ND2, H.generic_2d(**{"spots.enabled": False}), tmp_path / "cut", prov=None,
                                 z_range=(2, 3))
    out = Path(res["out_dir"])
    manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert "planar" not in manifest and "spacing_source" not in manifest and manifest["z_range"] == [2, 3]
    assert "read:planar" not in res["qc_flags"]
    nuclei = res["tables"]["nuclei"]
    assert "area_um2" not in nuclei.columns and nuclei["volume_um3"].notna().all()
    assert "planar" not in json.loads((out / f"{H.STEM}__done.json").read_text(encoding="utf-8"))


def test_true_refuses_more_than_one_plane(tmp_path, monkeypatch):
    with pytest.raises(planar.PlanarConfigError) as err:
        planar.decide(_cfg(**{"io.planar": True}), _meta(z=5), _stack(nz=5))
    assert "--z-range" in str(err.value) and "segmentation.nuclei.project" in str(err.value)
    # in process_image: before any stage, before any file is written (a stale marker is still cleared)
    H.patch_reader(monkeypatch, H.stack_3d(), two_d=False)
    out = tmp_path / "r"
    out.mkdir()
    (out / f"{H.STEM}__done.json").write_text("{}", encoding="utf-8")
    seen = []
    monkeypatch.setattr(planar, "segment_planar", lambda *a, **k: seen.append("segment_planar"))
    monkeypatch.setattr(pipeline, "segment_nuclei", lambda *a, **k: seen.append("segment_nuclei"))
    with pytest.raises(planar.PlanarConfigError, match="--z-range"):
        pipeline.process_image(ND2, H.generic_2d(**{"io.planar": True}), out, prov=None)
    assert not seen and not list(out.iterdir())
    # cut to one plane, true runs it as a single plane
    assert planar.decide(_cfg(**{"io.planar": True}), _meta(z=5), _stack(nz=1))["planar"] is True


# ---- the spacing chain -------------------------------------------------------------------------------------

def test_the_spacing_chain():
    auto = {"io.planar": "auto"}
    yx = planar.decide(_cfg(**auto), _meta(z=1), _stack(spacing=(1.0, 1.0, 1.0), spacing_ok=False,
                                                         spacing_yx=(0.4, 0.5)))
    assert (yx["dy"], yx["dx"], yx["spacing_source"], yx["flags"]) == (0.4, 0.5, "file_yx", ["read:planar"])
    f = planar.decide(_cfg(**auto), _meta(), _stack(spacing=(2.0, 0.3, 0.35)))
    assert (f["dy"], f["dx"], f["spacing_source"]) == (0.3, 0.35, "file")
    guess = planar.decide(_cfg(**auto), _meta(), _stack(spacing=(1.0, 1.0, 1.0), spacing_ok=False))
    assert (guess["dy"], guess["dx"], guess["spacing_source"]) == (1.0, 1.0, "assumed")
    assert guess["flags"] == ["read:planar", "voxel:UNREADABLE_xy_assumed_1um"]
    strided = planar.decide(_cfg(**auto), _meta(), _stack(spacing=(1.0, 2.0, 2.0), spacing_ok=False, stride=2))
    assert (strided["dy"], strided["dx"]) == (2.0, 2.0)               # 1 um per file pixel, then the stride
    over = planar.decide(_cfg(**auto, **{"io.voxel_size_um": 0.25}), _meta(), _stack(spacing_yx=(0.4, 0.5)))
    assert (over["dy"], over["dx"], over["spacing_source"]) == (0.25, 0.25, "override")
    assert over["flags"] == ["read:planar", "voxel:override", "voxel:override_differs_from_file"]


def test_voxel_size_override_with_stride_and_flags(tmp_path, monkeypatch):
    s = 2
    file_sp = (0.5, 0.8 * s, 0.8 * s)                 # the reader's spacing after a stride of 2 (file 0.5, 0.8, 0.8)
    stack = _stack(nz=5, spacing=file_sp, stride=s)
    # the three forms, stride applied after; 3D runs keep the file's dz unless the override gives one (D6)
    assert planar.voxel_override(_cfg(**{"io.voxel_size_um": 0.8}), stack)["spacing"] == (0.5, 1.6, 1.6)
    assert planar.voxel_override(_cfg(**{"io.voxel_size_um": [0.8, 0.7]}), stack)["spacing"] == (0.5, 1.6, 1.4)
    assert planar.voxel_override(_cfg(**{"io.voxel_size_um": [0.3, 0.8, 0.8]}), stack)["spacing"] == (0.3, 1.6, 1.6)
    # voxel:override always, override_differs_from_file only beyond 1 percent of the file's own calibration
    assert planar.voxel_override(_cfg(**{"io.voxel_size_um": 0.805}), stack)["flags"] == ["voxel:override"]
    assert planar.voxel_override(_cfg(**{"io.voxel_size_um": 0.81}), stack)["flags"] == [
        "voxel:override", "voxel:override_differs_from_file"]
    assert planar.voxel_override(_cfg(**{"io.voxel_size_um": [0.6, 0.8, 0.8]}), stack)["flags"] == [
        "voxel:override", "voxel:override_differs_from_file"]
    # a file without calibration: nothing to differ from; a 3D dz that stays the 1 um guess keeps today's flag
    guess = _stack(nz=5, spacing=(1.0, 1.0, 1.0), spacing_ok=False)
    assert planar.voxel_override(_cfg(**{"io.voxel_size_um": 0.8}), guess) == {
        "spacing": (1.0, 0.8, 0.8), "flags": ["voxel:override", "voxel:UNREADABLE_assumed_isotropic_1um"],
        "spacing_source": "override"}
    assert planar.voxel_override(_cfg(**{"io.voxel_size_um": [0.4, 0.8, 0.8]}), guess)["flags"] == ["voxel:override"]
    assert planar.voxel_override(_cfg(), stack) is None                      # the key absent: every existing run
    # a planar run with stride: (dy, dx) of the file before the stride, times the stride
    d = planar.decide(_cfg(**{"io.planar": "auto", "io.voxel_size_um": [0.4, 0.5, 0.6]}), _meta(), _stack(stride=4))
    assert (d["dy"], d["dx"], d["spacing_source"]) == (2.0, 2.4, "override")
    # through process_image: a 3D run with the override (D6), and the same run without it records nothing new
    H.patch_reader(monkeypatch, H.stack_3d(), two_d=False)
    base = {"spots.enabled": False}
    res = pipeline.process_image(ND2, H.generic_2d(**base, **{"io.voxel_size_um": [0.25, 0.4, 0.4]}),
                                 tmp_path / "o", prov=None, xy_stride=2)
    flags = res["qc_flags"]
    assert "voxel:override" in flags and "voxel:override_differs_from_file" in flags and "read:planar" not in flags
    nuclei = pd.read_csv(tmp_path / "o" / f"{H.STEM}__nuclei.csv")
    assert nuclei[["voxel_dz_um", "voxel_dy_um", "voxel_dx_um"]].iloc[0].tolist() == [0.25, 0.8, 0.8]
    manifest = json.loads((tmp_path / "o" / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["spacing_source"] == "override" and "planar" not in manifest
    plain = pipeline.process_image(ND2, H.generic_2d(**base), tmp_path / "p", prov=None, xy_stride=2)
    assert not [f for f in plain["qc_flags"] if f.startswith("voxel:")]
    assert "spacing_source" not in json.loads((tmp_path / "p" / "run_manifest.json").read_text(encoding="utf-8"))
    assert (pd.read_csv(tmp_path / "p" / f"{H.STEM}__nuclei.csv")["voxel_dy_um"] == 2 * H.DY).all()


# ---- recipe keys -------------------------------------------------------------------------------------------

class _Cfg:
    def __init__(self, data):
        self._data = data

    def get(self, dotted, default=None):
        node = self._data
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node


@pytest.mark.parametrize("key,value", [
    ("planar", None), ("planar", False), ("planar", True), ("planar", "auto"),
    ("voxel_size_um", None), ("voxel_size_um", 0.5), ("voxel_size_um", 2), ("voxel_size_um", [0.5, 0.25]),
    ("voxel_size_um", [1, 0.5, 0.5]), ("voxel_size_um", (0.3, 0.2, 0.2)),
])
def test_io_keys_accept(key, value):
    cfg = _Cfg({"io": {key: value}})
    assert recipe_keys.problems(cfg) == []
    recipe_keys.check(cfg)


@pytest.mark.parametrize("key,value", [
    ("planar", "true"), ("planar", "false"), ("planar", "yes"), ("planar", "Auto"), ("planar", 1), ("planar", 0),
    ("voxel_size_um", 0), ("voxel_size_um", -0.5), ("voxel_size_um", "0.5"), ("voxel_size_um", True),
    ("voxel_size_um", [0.5]), ("voxel_size_um", [0.5, 0.5, 0.5, 0.5]), ("voxel_size_um", [0.5, -1]),
    ("voxel_size_um", [True, 0.5]), ("voxel_size_um", {"dy": 0.5}), ("voxel_size_um", float("nan")),
    ("voxel_size_um", [0.5, "0.5"]),
])
def test_io_keys_refuse(key, value):
    cfg = _Cfg({"io": {key: value}})
    (problem,) = recipe_keys.problems(cfg)
    assert problem.startswith(f"io.{key} must be"), problem
    with pytest.raises(ValueError, match=re.escape(f"io.{key}")):
        recipe_keys.check(cfg)


def test_a_malformed_io_key_stops_the_image_before_any_stage(tmp_path, monkeypatch):
    H.patch_reader(monkeypatch)
    for key, value in (("io.planar", "yes"), ("io.voxel_size_um", [0.5])):
        out = tmp_path / key
        with pytest.raises(ValueError, match=re.escape(key)):
            pipeline.process_image(ND2, H.generic_2d(**{key: value}), out, prov=None)
        assert out.is_dir() and not list(out.iterdir())


# ---- resume (D8) -------------------------------------------------------------------------------------------

def test_is_done_planar_marker_safe_default(tmp_path):
    cfg = load_config(ROOT / "config" / "config.yaml")
    tools = {t: "PINNED" for ts in stages._STAGE_TOOLS.values() for t in ts}
    h3d = stages.stage_hashes(cfg, {}, tools, "unused")["segment"]
    h2d = stages.stage_hashes(cfg, {}, tools, "unused", planar=True)["segment"]
    assert h3d != h2d

    def marker(iid, seg, **extra):
        provenance.write_done_marker(tmp_path, iid, config_hash=cfg.hash, enabled_stages=batch.enabled_stage_names(cfg),
                                     stage_hashes={"segment": seg}, extra=extra or None)

    marker("p", h2d, planar=True)
    done, why = batch.is_done(tmp_path, "p", cfg, segment_hash=h3d)
    assert done is False and "planar flag" in why
    assert batch.is_done(tmp_path, "p", cfg, segment_hash=h3d, segment_hash_planar=h2d) == (True, "done")
    assert batch.is_done(tmp_path, "p", cfg, segment_hash_planar=h2d) == (True, "done")
    assert batch.is_done(tmp_path, "p", cfg)[0] is True                         # no segment hash at all: today
    assert batch.is_done(tmp_path, "p", cfg, segment_hash=h3d, segment_hash_planar="other")[0] is False
    marker("d", h3d)                                                            # a 3D marker: today's rule
    assert batch.is_done(tmp_path, "d", cfg, segment_hash=h3d, segment_hash_planar=h2d) == (True, "done")
    assert batch.is_done(tmp_path, "d", cfg, segment_hash_planar=h2d) == (True, "done")
    assert batch.is_done(tmp_path, "d", cfg, segment_hash="zzz", segment_hash_planar=h2d)[0] is False


# ---- planar plus a projection (D7) -------------------------------------------------------------------------

def test_planar_wins_over_a_projection(tmp_path, monkeypatch):
    H.patch_reader(monkeypatch)
    res = pipeline.process_image(ND2, H.generic_2d(**{"segmentation.nuclei.project": "max", "spots.enabled": False}),
                                 tmp_path / "p", prov=None)
    out = Path(res["out_dir"])
    assert "segment:project_ignored_planar" in res["qc_flags"] and "segment:projected_max" not in res["qc_flags"]
    record = json.loads((out / f"{H.STEM}__stages.json").read_text(encoding="utf-8"))
    assert "segment:project_ignored_planar" in record["stages"]["segment"]["flags"]
    assert "segment_projection" not in record and record["planar"] is True
    assert res["n_nuclei"] == 9
    # an invalid projection mode is still refused before any pixel is processed
    with pytest.raises(planar.PlanarConfigError, match="segmentation.nuclei.project must be one of"):
        pipeline.process_image(ND2, H.generic_2d(**{"segmentation.nuclei.project": "median"}), tmp_path / "bad",
                               prov=None)


# ---- the readers (D11) -------------------------------------------------------------------------------------

def test_stack_spacing_yx_is_the_last_field_with_a_default():
    fields = dataclasses.fields(Stack)
    assert fields[-1].name == "spacing_yx" and fields[-1].default is None
    assert [f.name for f in fields[:-1]] == ["data", "spacing", "channel_names", "path", "downsample_xy",
                                             "spacing_ok", "t_index"]
    assert Stack(data=np.zeros((1, 1, 2, 2)), spacing=(1.0, 1.0, 1.0), channel_names=[], path=Path("x")).spacing_yx is None


class _Dims:
    def __init__(self, sizes):
        self.order = "".join(sizes)
        self.__dict__.update(sizes)


class _Fake:
    def __init__(self, z, pps):
        import dask.array as da

        data = np.zeros((2, z, 16, 16), np.uint16)
        self._dask = da.from_array(data, chunks=data.shape)
        self.scenes = ("Image:0",)
        self.channel_names = list(H.CHANNELS)
        self.dims = _Dims({"T": 1, "C": 2, "Z": z, "Y": 16, "X": 16})
        self.physical_pixel_sizes = namedtuple("PPS", "Z Y X")(*pps)
        self.dtype = data.dtype
        self.reader = self

    def get_image_dask_data(self, order, **kw):
        return self._dask


def test_the_bioio_reader_fills_spacing_yx_when_z_is_1(tmp_path, monkeypatch):
    p = tmp_path / "g.ome.tif"
    p.write_bytes(b"x")
    for z, pps, want in ((1, (None, 0.4, 0.5), (0.4, 0.5)), (1, (0.3, 0.4, 0.5), (0.4, 0.5)),
                         (1, (None, None, None), None), (1, (None, 0.0, 0.5), None), (1, (None, float("nan"), 0.5), None)):
        monkeypatch.setattr(bioio_reader, "_open", lambda path, z=z, pps=pps: _Fake(z, pps))
        st, meta = read_image(p), read_image_metadata(p)
        assert st.spacing_yx == want and meta["spacing_yx_um"] == want, (z, pps)
        assert st.spacing_ok is (pps[0] is not None and want is not None)       # _spacing's rule, untouched
        if want:
            assert read_image(p, xy_stride=2).spacing_yx == (2 * want[0], 2 * want[1])
    monkeypatch.setattr(bioio_reader, "_open", lambda path: _Fake(3, (0.3, 0.4, 0.5)))
    st, meta = read_image(p), read_image_metadata(p)
    assert st.spacing_yx is None and "spacing_yx_um" not in meta                # a 3D file: nothing new
    assert read_image(p, z_range=(1, 2)).spacing_yx is None


# ---- `germquant info` --------------------------------------------------------------------------------------

def test_info_prints_single_plane(monkeypatch, capsys):
    one = {"sizes": {"T": 1, "C": 2, "Z": 1, "Y": 256, "X": 256}, "spacing": (1.0, 1.0, 1.0), "spacing_ok": False,
           "channel_names": list(H.CHANNELS), "dtype": "uint16", "is_2d": True, "n_timepoints": 1,
           "spacing_yx_um": (0.8, 0.8)}
    monkeypatch.setattr(cli, "read_image_metadata", lambda path: dict(one))
    assert cli.main(["info", "x.ome.tif"]) == 0
    lines = capsys.readouterr().out.splitlines()
    plane = [ln for ln in lines if "single plane" in ln]
    assert len(plane) == 1 and "dy=0.8000 dx=0.8000" in plane[0]
    three = {**one, "sizes": {"T": 1, "C": 2, "Z": 5, "Y": 256, "X": 256}, "spacing": (0.5, 0.8, 0.8),
             "spacing_ok": True, "is_2d": False}
    three.pop("spacing_yx_um")
    monkeypatch.setattr(cli, "read_image_metadata", lambda path: dict(three))
    assert cli.main(["info", "x.ome.tif"]) == 0
    assert "single plane" not in capsys.readouterr().out
    nd2_2d = {**three, "sizes": {"C": 2, "Y": 256, "X": 256}, "is_2d": True}     # an .nd2 plane has no Z axis
    monkeypatch.setattr(cli, "read_image_metadata", lambda path: dict(nd2_2d))
    assert cli.main(["info", "x.ome.tif"]) == 0
    out = capsys.readouterr().out
    assert "single plane" in out and "dy=0.8000 dx=0.8000" in out
