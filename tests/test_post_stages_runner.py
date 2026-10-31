"""The post-stage runner on a hand-made result folder (design 5.6; plan 4.5, phase CORE; decisions D28 and D32 to
D34; P5b acceptance 4).

* The fixture post-stage writes only <image_id>__post_<name>_<table>.csv and <image_id>__post.json; every file
  that was there before is byte-identical afterwards, and no post-stage means no file at all.
* The sub-hash follows the parameters, the post_stages.yaml digest and the input digests (the main run's read
  sub-hash and run_timestamp, the tables and the label image the post-stage reads, the spacing), never the raw
  image (D32). The spacing comes from the main run's image_summary table on every route.
* Resume skips a stage whose sub-hash matches and whose tables are on disk; a failure is recorded in
  __post.json and the image stays done (acceptance 4). A post-stage that breaks the contract fails alone.
* The image is read through the run's reader at the run's stride and z range, once per image (D34).
* In a fresh interpreter, run_post_stages imports germquant.pipeline before the post-stage or a reader (D28).
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import re

import numpy as np
import pandas as pd
import pytest
from post_helpers import (  # noqa: F401
    GENERIC,
    SPACING,
    edit_marker,
    fixture_post_stage,
    make_result_folder,
    make_stack,
    post_record,
    run_child,
    snapshot,
)

from germquant import post_stages as PS
from germquant.batch import is_done
from germquant.post_stages import testing as T

N = T.FIXTURE_NAME
HEX12 = re.compile(r"[0-9a-f]{12}")


def _table_file(image_id: str, name: str = N, table: str = T.FIXTURE_TABLE) -> str:
    return f"{image_id}__post_{name}_{table}.csv"


def _run(rf, entries, **kw):
    kw.setdefault("stack", rf.stack)
    return PS.run_post_stages(rf.out_dir, rf.image_id, entries, cfg=rf.cfg, image_path=rf.image_path, **kw)


@pytest.mark.usefixtures("fixture_post_stage")
def test_the_fixture_writes_only_its_declared_table_and_the_record(tmp_path):
    rf = make_result_folder(tmp_path)
    before = snapshot(rf.out_dir)
    empty = _run(rf, [])
    assert empty["ran"] == empty["skipped"] == empty["failed"] == []
    assert snapshot(rf.out_dir) == before, "a run without post-stages wrote a file"
    res = _run(rf, [{"name": N}])
    after = snapshot(rf.out_dir)
    assert set(after) - set(before) == {_table_file(rf.image_id), f"{rf.image_id}__post.json"}
    assert {k: after[k] for k in before} == before, "a file of the main run changed"
    assert res["ran"] == [N] and res["failed"] == [] and res["skipped"] == []
    table = pd.read_csv(rf.out_dir / _table_file(rf.image_id))
    assert list(table.columns) == list(T.FIXTURE_COLUMNS)
    foci = rf.stack.data[1].astype(np.float64)
    voxel = SPACING[0] * SPACING[1] * SPACING[2]
    assert sorted(table["nucleus_id"]) == [1, 2, 3]
    for row in table.itertuples():
        inside = rf.labels == row.nucleus_id
        assert row.n_voxels == int(inside.sum())
        assert row.volume_um3 == pytest.approx(inside.sum() * voxel)
        assert row.foci_value == pytest.approx(float(foci[inside].mean()))
        assert row.foci_value_scaled == pytest.approx(row.foci_value)
    assert table.set_index("nucleus_id")["n_spots"].to_dict() == {1: 2, 2: 0, 3: 1}
    rec = post_record(rf.out_dir, rf.image_id)
    assert rec["format"] == 1 and rec["image_id"] == rf.image_id
    assert rec["post_stages_file_sha256"] == PS.post_stages_sha256([{"name": N}])
    assert set(rec["stages"]) == {N}
    st = rec["stages"][N]
    assert set(st) == {"version", "params", "sub_hash", "done", "error", "tables", "finished", "inputs"}
    assert st["version"] == T.FIXTURE_VERSION and st["done"] is True and st["error"] is None
    assert st["params"] == {"factor": 1.0, "min_voxels": 0, "statistic": "mean", "fail": False}
    assert HEX12.fullmatch(st["sub_hash"]) and st["tables"] == [_table_file(rf.image_id)] and st["finished"]
    assert res["post"] == rec and res["record"] == str(rf.out_dir / f"{rf.image_id}__post.json")
    # the parameters reach the post-stage
    _run(rf, [{"name": N, "params": {"factor": 2.0, "statistic": "max", "min_voxels": 200}}])
    table2 = pd.read_csv(rf.out_dir / _table_file(rf.image_id)).set_index("nucleus_id")
    kept = [i for i in (1, 2, 3) if (rf.labels == i).sum() >= 200]
    assert sorted(table2.index) == kept and len(kept) < 3
    for i in kept:
        assert table2.loc[i, "foci_value"] == pytest.approx(float(foci[rf.labels == i].max()))
        assert table2.loc[i, "foci_value_scaled"] == pytest.approx(2.0 * table2.loc[i, "foci_value"])


@pytest.mark.usefixtures("fixture_post_stage")
def test_the_sub_hash_follows_params_the_file_digest_and_the_inputs(tmp_path):
    """D32: the parameters, the post_stages.yaml digest, the main run's read sub-hash and run_timestamp, the
    digest of each table and of the label image the post-stage reads, and the spacing. Tables are digested on
    LF-normalised bytes; a table it does not read and the raw image itself never enter."""
    rf = make_result_folder(tmp_path)
    out, iid = rf.out_dir, rf.image_id

    def sub(entries=None):
        res = _run(rf, [{"name": N}] if entries is None else entries, resume=False)
        assert res["ran"] == [N], res["post"]
        return post_record(out, iid)["stages"][N]["sub_hash"]

    base = sub()
    assert sub() == base, "the sub-hash is not reproducible"
    assert sub([{"name": N, "params": {"factor": 2.0}}]) != base
    f = tmp_path / "post_stages.yaml"
    PS.write_post_stages_file(f, [{"name": N}])
    assert sub(f) == base, "a file written from the entries digests like the entries"
    f.write_bytes(b"# switched on for the screen\n" + f.read_bytes())
    assert sub(f) != base
    assert post_record(out, iid)["post_stages_file_sha256"] == hashlib.sha256(f.read_bytes()).hexdigest()

    def changed(path, mutate) -> bool:
        old = path.read_bytes()
        mutate(path)
        try:
            return sub() != base
        finally:
            path.write_bytes(old)

    def edit(column: str, value):
        def _do(p):
            df = pd.read_csv(p)
            assert column in df.columns
            df[column] = value
            df.to_csv(p, index=False)
        return _do

    def relabel(p):
        import tifffile

        lab = tifffile.imread(p)
        lab[lab == 2] = 0
        tifffile.imwrite(p, lab, compression="zlib", photometric="minisblack")

    def line_endings(p):
        data = p.read_bytes().replace(b"\r\n", b"\n")
        p.write_bytes(data if data != p.read_bytes() else data.replace(b"\n", b"\r\n"))

    marker = out / f"{iid}__done.json"
    assert changed(marker, lambda p: edit_marker(out, iid, stage_hashes={"read": "ffffffffffff"}))
    assert changed(out / f"{iid}__spots.csv", edit("y_um", 1.5))
    assert changed(out / f"{iid}__nuclei.csv", edit("run_timestamp", "2026-10-24T00:00:00+00:00"))
    assert changed(out / f"{iid}__nuclei_labels.tif", relabel)
    assert changed(out / f"{iid}__image_summary.csv", edit("voxel_dz_um", 0.6)), "the spacing is an input"
    assert changed(out / f"{iid}__image_summary.csv", edit("channel_map_json", '{"dna": 1, "foci": 0}')), \
        "the roles are an input"
    assert changed(out / f"{iid}__image_summary.csv", edit("run_timestamp", "2026-10-25T00:00:00+00:00")), \
        "the main run is an input"
    assert not changed(out / f"{iid}__image_summary.csv", edit("qc_pass", False)), \
        "only the spacing, the roles and the run_timestamp"
    assert not changed(out / f"{iid}__nuclei.csv", line_endings), "a table digest is not LF-normalised"
    assert not changed(out / f"{iid}__coloc.csv", edit("overlap", 0.75)), "a table it does not read"
    assert not changed(rf.image_path, lambda p: p.write_bytes(b"other raw pixels\n")), "the raw image is hashed"
    assert not changed(out / f"{iid}__montage.png", lambda p: p.write_bytes(b"x")), "a file it does not read"
    assert sub() == base


_RING = {"version": "0.1.0", "reads": {"roles": ["foci"], "tables": [], "labels": True, "image": True},
         "writes": {"rings": ["nucleus_id", "inside_mean"]}}


def _inside_mean(inputs, params):
    lab, foci = inputs.labels(), inputs.image("foci")
    ids = [int(i) for i in np.unique(lab) if i]
    return {"rings": pd.DataFrame({"nucleus_id": ids, "inside_mean": [float(foci[lab == i].mean()) for i in ids]})}


def test_a_new_main_run_restales_a_stage_that_reads_no_table(tmp_path):
    """Review S3: the read sub-hash never covers pixels, so a post-stage that reads no table (the shape of the
    doc's ring_intensity: labels and image only) kept its old numbers when the main run was redone on a
    re-exported image whose labels came out byte-identical. The main run's run_timestamp (image_summary) is an
    input, so a new main run restales every post-stage, while a resumed call on the same main run still skips."""
    from germquant.io.nd2_reader import Stack

    rf = make_result_folder(tmp_path)
    PS.register_post_stage("ring_intensity", _inside_mean, _RING, replace=True)
    e = [{"name": "ring_intensity"}]
    rings = rf.out_dir / _table_file(rf.image_id, "ring_intensity", "rings")
    try:
        assert _run(rf, e)["ran"] == ["ring_intensity"]
        old = pd.read_csv(rings)["inside_mean"].to_numpy()
        inputs = post_record(rf.out_dir, rf.image_id)["stages"]["ring_intensity"]["inputs"]
        assert inputs["main_run"] == "2026-10-23T00:00:00+00:00" and inputs["tables"] == {}
        # the main run redone on a re-exported image: the foci channel doubled, the DNA channel (so the labels)
        # unchanged, and every table written again with the new run's run_timestamp
        data = rf.stack.data.copy()
        data[1] *= 2.0
        redone = Stack(data=data, spacing=rf.stack.spacing, channel_names=list(rf.stack.channel_names),
                       path=rf.stack.path, spacing_ok=True)
        for table in ("nuclei", "spots", "image_summary", "coloc"):
            p = rf.out_dir / f"{rf.image_id}__{table}.csv"
            df = pd.read_csv(p)
            df["run_timestamp"] = "2026-10-25T09:00:00+00:00"
            df.to_csv(p, index=False)
        res = _run(rf, e, stack=redone)
        assert res["ran"] == ["ring_intensity"] and res["skipped"] == [], res
        assert pd.read_csv(rings)["inside_mean"].to_numpy() == pytest.approx(2.0 * old)
        assert _run(rf, e, stack=redone)["skipped"] == ["ring_intensity"], "the same main run must still resume"
    finally:
        PS.unregister_post_stage("ring_intensity")


def test_resume_skips_a_matching_sub_hash(tmp_path):
    rf = make_result_folder(tmp_path)
    calls = []
    fixture = T.fixture_function()

    def counting(inputs, params):
        calls.append(dict(params))
        return fixture(inputs, params)

    def second(inputs, params):
        calls.append({"second": True})
        return {"counts": pd.DataFrame({"n_nuclei": [len(inputs.tables("nuclei"))]})}

    PS.register_post_stage(N, counting, T.fixture_declaration(), replace=True)
    PS.register_post_stage("second_stage", second, {"version": "0.1.0", "reads": {"tables": ["nuclei"]},
                                                    "writes": {"counts": ["n_nuclei"]}, "params": []}, replace=True)
    record = rf.out_dir / f"{rf.image_id}__post.json"
    two = [{"name": N, "params": {"factor": 2.0}}]
    try:
        assert _run(rf, [{"name": N}])["ran"] == [N] and len(calls) == 1
        data = record.read_bytes()
        res = _run(rf, [{"name": N}])
        assert res["skipped"] == [N] and res["ran"] == [] and len(calls) == 1
        assert record.read_bytes() == data, "a resumed run rewrote the record"
        assert _run(rf, two)["ran"] == [N] and len(calls) == 2
        assert _run(rf, two)["skipped"] == [N] and len(calls) == 2
        (rf.out_dir / _table_file(rf.image_id)).unlink()
        assert _run(rf, two)["ran"] == [N] and len(calls) == 3, "a missing post table must rerun its stage"
        assert _run(rf, two, resume=False)["ran"] == [N] and len(calls) == 4
        fail = [{"name": N, "params": {"fail": True}}]
        assert _run(rf, fail)["failed"] == [N] and len(calls) == 5
        assert _run(rf, fail)["failed"] == [N] and len(calls) == 6, "a failed stage must be retried"
        # a stage taken out of post_stages.yaml leaves the record; its table stays on disk
        both = [{"name": N}, {"name": "second_stage"}]
        assert _run(rf, both)["ran"] == [N, "second_stage"]
        assert set(post_record(rf.out_dir, rf.image_id)["stages"]) == {N, "second_stage"}
        assert _run(rf, [{"name": N}])["ran"] == [N], "the file digest changed, so the stage reruns"
        assert set(post_record(rf.out_dir, rf.image_id)["stages"]) == {N}
        assert (rf.out_dir / _table_file(rf.image_id, "second_stage", "counts")).is_file()
    finally:
        PS.unregister_post_stage(N)
        PS.unregister_post_stage("second_stage")


@pytest.mark.usefixtures("fixture_post_stage")
def test_an_empty_list_leaves_a_record_with_no_stages(tmp_path):
    """Review S4: the record lists the post-stages of the list in force, also when that list is empty (every
    post-stage switched off): a record an earlier list left is rewritten once with no stages and the digest of
    the list or file in force, and the old tables stay on disk. A folder that never had a record still gets
    none, so a run without post-stages writes nothing."""
    rf = make_result_folder(tmp_path)
    record = rf.out_dir / f"{rf.image_id}__post.json"
    assert _run(rf, [{"name": N}])["ran"] == [N]
    res = _run(rf, [])
    rec = post_record(rf.out_dir, rf.image_id)
    assert rec == {"format": 1, "image_id": rf.image_id, "post_stages_file_sha256": PS.post_stages_sha256([]),
                   "stages": {}}, rec
    assert res["post"] == rec and res["ran"] == res["skipped"] == res["failed"] == [] and "error" not in res
    assert (rf.out_dir / _table_file(rf.image_id)).is_file(), "an old table must stay on disk"
    stamp = record.stat().st_mtime_ns
    assert _run(rf, [])["post"] == rec and record.stat().st_mtime_ns == stamp, "the record was written again"
    written = tmp_path / "post_stages.yaml"
    PS.write_post_stages_file(written, [])
    assert _run(rf, written)["post"] == rec and record.stat().st_mtime_ns == stamp, \
        "an empty file written by germquant digests like the empty list"
    by_hand = tmp_path / "by_hand.yaml"
    by_hand.write_bytes(b"# every post-stage off for now\n")
    assert _run(rf, by_hand)["post"] == {**rec, "post_stages_file_sha256": PS.post_stages_sha256(by_hand)}
    assert post_record(rf.out_dir, rf.image_id)["stages"] == {}
    fresh = make_result_folder(tmp_path / "fresh")
    before = snapshot(fresh.out_dir)
    for empty in ([], None, by_hand):
        assert _run(fresh, empty)["post"] is None
    assert snapshot(fresh.out_dir) == before, "a folder without a record gained a file"


@pytest.mark.usefixtures("fixture_post_stage")
def test_failure_recorded_image_stays_done(tmp_path):
    """Acceptance 4: a failing post-stage is recorded in __post.json and the image stays done (the completion
    marker is byte-identical and is_done still says done). Entries that cannot run are recorded, never raised."""
    rf = make_result_folder(tmp_path)
    marker = rf.out_dir / f"{rf.image_id}__done.json"
    marker_bytes = marker.read_bytes()
    assert is_done(rf.out_dir, rf.image_id, rf.cfg) == (True, "done")
    before = snapshot(rf.out_dir)
    res = _run(rf, [{"name": N, "params": {"fail": True}}])
    assert res["failed"] == [N] and res["ran"] == [] and res["skipped"] == []
    st = post_record(rf.out_dir, rf.image_id)["stages"][N]
    assert st["done"] is False and st["tables"] == [] and HEX12.fullmatch(st["sub_hash"])
    assert "RuntimeError" in st["error"] and "asked to fail" in st["error"]
    assert set(snapshot(rf.out_dir)) - set(before) == {f"{rf.image_id}__post.json"}
    assert marker.read_bytes() == marker_bytes
    assert is_done(rf.out_dir, rf.image_id, rf.cfg) == (True, "done")
    res = _run(rf, [{"name": "not_installed"}, {"name": N, "params": {"factor": 99.0}}, {"name": N}])
    assert res["failed"] == ["not_installed", N] and res["ran"] == []
    stages = post_record(rf.out_dir, rf.image_id)["stages"]
    assert stages["not_installed"]["done"] is False and "not installed" in stages["not_installed"]["error"]
    assert "factor" in stages[N]["error"] and "more than once" in stages[N]["error"]
    assert marker.read_bytes() == marker_bytes
    assert is_done(rf.out_dir, rf.image_id, rf.cfg) == (True, "done")


def test_two_post_stages_that_would_write_one_file_are_refused(tmp_path):
    """Review S2: the post table name <image_id>__post_<name>_<table>.csv is frozen and names may hold single
    underscores, so post-stage a with table b_c and post-stage a_b with table c would both write
    <image_id>__post_a_b_c.csv (the last one winning while both were recorded done, and resume then skipping
    both). The pair is refused wherever a list is checked: validation, the writer and the runner, which fails
    both and writes no table. Each one alone stays valid."""
    rf = make_result_folder(tmp_path)
    reads = {"tables": ["nuclei"]}
    PS.register_post_stage("a", lambda i, p: {"b_c": pd.DataFrame({"who": ["a"]})},
                           {"version": "1.0.0", "reads": reads, "writes": {"b_c": ["who"]}}, replace=True)
    PS.register_post_stage("a_b", lambda i, p: {"c": pd.DataFrame({"who": ["a_b"]})},
                           {"version": "1.0.0", "reads": reads, "writes": {"c": ["who"]}}, replace=True)
    assert PS.post_table_name(rf.image_id, "a", "b_c") == PS.post_table_name(rf.image_id, "a_b", "c")
    pair = [{"name": "a"}, {"name": "a_b"}]
    try:
        probs = PS.validate_post_stages(pair)
        assert len(probs) == 2 and all("__post_a_b_c.csv" in p for p in probs), probs
        assert probs[0].startswith("entry 0:") and probs[1].startswith("entry 1:"), probs
        with pytest.raises(ValueError, match="__post_a_b_c.csv"):
            PS.write_post_stages_file(tmp_path / "post_stages.yaml", pair)
        assert not (tmp_path / "post_stages.yaml").exists()
        before = snapshot(rf.out_dir)
        res = _run(rf, pair)
        assert res["failed"] == ["a", "a_b"] and res["ran"] == [], res
        stages = post_record(rf.out_dir, rf.image_id)["stages"]
        for name in ("a", "a_b"):
            assert "would both write" in stages[name]["error"] and stages[name]["tables"] == [], stages[name]
        assert set(snapshot(rf.out_dir)) - set(before) == {f"{rf.image_id}__post.json"}
        assert len(PS.validate_post_stages([*pair, {"name": "a_b"}])) == 3, "the pair once, plus the repeat"
        for alone in ("a", "a_b"):
            assert PS.validate_post_stages([{"name": alone}]) == []
            assert _run(rf, [{"name": alone}])["ran"] == [alone]
    finally:
        PS.unregister_post_stage("a")
        PS.unregister_post_stage("a_b")


def _frame(columns=T.FIXTURE_COLUMNS):
    return pd.DataFrame({c: [1] for c in columns})


def _writes_itself(inputs, params):
    (inputs.out_dir / "stray.txt").write_text("a post-stage must not write files itself\n", encoding="ascii")
    return {T.FIXTURE_TABLE: _frame()}


def _mutates_labels(inputs, params):
    inputs.labels()[0, 0, 0] = 5
    return {T.FIXTURE_TABLE: _frame()}


def _exits(inputs, params):
    raise SystemExit(3)                              # a post-stage must not end the caller's image loop


_BREACHES = {
    "reads_undeclared": (lambda i, p: i.tables("coloc"), "not declared"),
    "role_undeclared": (lambda i, p: i.image("dna"), "not declared"),
    "not_a_mapping": (lambda i, p: [1], "mapping"),
    "missing_table": (lambda i, p: {}, T.FIXTURE_TABLE),
    "extra_table": (lambda i, p: {T.FIXTURE_TABLE: _frame(), "other": _frame()}, "other"),
    "extra_column": (lambda i, p: {T.FIXTURE_TABLE: _frame().assign(extra=1)}, "extra"),
    "missing_column": (lambda i, p: {T.FIXTURE_TABLE: _frame()[list(T.FIXTURE_COLUMNS[:-1])]}, "n_spots"),
    "not_a_frame": (lambda i, p: {T.FIXTURE_TABLE: [1, 2]}, "DataFrame"),
    "writes_itself": (_writes_itself, "changed the results folder"),
    "mutates_labels": (_mutates_labels, "read-only"),
    "exits": (_exits, "SystemExit"),
}


def test_a_post_stage_that_breaks_the_contract_fails_alone(tmp_path):
    rf = make_result_folder(tmp_path)
    decl = T.fixture_declaration()
    names = list(_BREACHES) + ["reordered"]
    for name, (fn, _fragment) in _BREACHES.items():
        PS.register_post_stage(name, fn, decl, replace=True)
    PS.register_post_stage("reordered", lambda i, p: {T.FIXTURE_TABLE: _frame(T.FIXTURE_COLUMNS[::-1])}, decl,
                           replace=True)
    try:
        before = snapshot(rf.out_dir)
        res = _run(rf, [{"name": n} for n in names])
        assert res["failed"] == list(_BREACHES) and res["ran"] == ["reordered"]
        stages = post_record(rf.out_dir, rf.image_id)["stages"]
        for name, (_fn, fragment) in _BREACHES.items():
            assert stages[name]["done"] is False and stages[name]["tables"] == [], stages[name]
            assert fragment in stages[name]["error"], (name, stages[name]["error"])
        written = pd.read_csv(rf.out_dir / _table_file(rf.image_id, "reordered"))
        assert list(written.columns) == list(T.FIXTURE_COLUMNS), "declared column order"
        after = snapshot(rf.out_dir)
        assert set(after) - set(before) == {f"{rf.image_id}__post.json", "stray.txt",
                                            _table_file(rf.image_id, "reordered")}
        assert {k: after[k] for k in before} == before
    finally:
        for name in names:
            PS.unregister_post_stage(name)


def test_a_post_stage_that_writes_into_a_subfolder_fails(tmp_path):
    """Review N2: the folder guard compares the whole tree of the results folder, so a post-stage that writes a
    subfolder of its own fails, and so does the retry that writes into the subfolder the first call left."""
    rf = make_result_folder(tmp_path)
    sizes = iter(range(1, 10))

    def debug_mask(inputs, params):
        sub = inputs.out_dir / "debug"
        sub.mkdir(exist_ok=True)
        (sub / "mask.tif").write_bytes(b"x" * next(sizes))
        return {T.FIXTURE_TABLE: _frame()}

    PS.register_post_stage("debug_mask", debug_mask, T.fixture_declaration(), replace=True)
    try:
        for call in ("first", "retry"):
            res = _run(rf, [{"name": "debug_mask"}])
            assert res["failed"] == ["debug_mask"] and res["ran"] == [], (call, res)
            error = post_record(rf.out_dir, rf.image_id)["stages"]["debug_mask"]["error"]
            assert "changed the results folder" in error and "debug/mask.tif" in error, (call, error)
            assert not (rf.out_dir / _table_file(rf.image_id, "debug_mask")).exists(), call
    finally:
        PS.unregister_post_stage("debug_mask")


@pytest.mark.usefixtures("fixture_post_stage")
def test_every_stage_record_has_the_same_keys(tmp_path):
    """Review N4: a stage refused before it could run (here: not installed, a parameter out of range) carries
    the same keys as one that ran, with inputs null, so a reader of __post.json never meets a missing key."""
    rf = make_result_folder(tmp_path)
    PS.register_post_stage("out_of_range", T.fixture_function(), T.fixture_declaration(), replace=True)
    try:
        res = _run(rf, [{"name": N}, {"name": "not_installed"},
                        {"name": "out_of_range", "params": {"factor": 99.0}}])
    finally:
        PS.unregister_post_stage("out_of_range")
    assert res["ran"] == [N] and res["failed"] == ["not_installed", "out_of_range"], res
    stages = res["post"]["stages"]
    keys = {"version", "params", "sub_hash", "done", "error", "tables", "finished", "inputs"}
    assert all(set(rec) == keys for rec in stages.values()), {k: sorted(v) for k, v in stages.items()}
    assert stages[N]["inputs"]["read_sub_hash"] == "0123456789ab"
    assert stages["not_installed"]["inputs"] is None and stages["out_of_range"]["inputs"] is None
    assert post_record(rf.out_dir, rf.image_id) == res["post"]


@pytest.mark.usefixtures("fixture_post_stage")
def test_the_image_is_read_through_the_run_reader_at_its_stride_and_z_range(tmp_path, monkeypatch):
    """Without a stack handed over, the image is read with the run's reader (process_image's seam), at the
    xy_stride and z_range of the completion marker, once per image however many post-stages read it."""
    from germquant import pipeline
    from germquant.io.nd2_reader import Stack

    full, labels_full = make_stack(nz=6)
    run_stack = Stack(data=np.ascontiguousarray(full.data[:, 1:4, ::2, ::2]), spacing=(0.5, 0.4, 0.4),
                      channel_names=list(full.channel_names), path=full.path, spacing_ok=True)
    run_labels = np.ascontiguousarray(labels_full[1:4, ::2, ::2])
    rf = make_result_folder(tmp_path, xy_stride=2, z_range=[1, 4], stack=run_stack, labels=run_labels)
    calls = []

    def read_meta(path, **kw):
        calls.append(("meta", str(path), kw))
        return {"channel_names": list(full.channel_names), "spacing": full.spacing, "spacing_ok": True}

    def read_stack(path, **kw):
        calls.append(("stack", str(path), kw))
        return run_stack

    monkeypatch.setattr(pipeline, "read_nd2_metadata", read_meta)
    monkeypatch.setattr(pipeline, "read_stack", read_stack)
    PS.register_post_stage("second_reader", T.fixture_function(), T.fixture_declaration(), replace=True)
    try:
        res = PS.run_post_stages(rf.out_dir, rf.image_id, [{"name": N}, {"name": "second_reader"}], cfg=rf.cfg,
                                 image_path=rf.image_path)
    finally:
        PS.unregister_post_stage("second_reader")
    assert res["ran"] == [N, "second_reader"], res["post"]
    assert [c for c in calls if c[0] == "stack"] == [("stack", str(rf.image_path), {"xy_stride": 2,
                                                                                     "z_range": (1, 4)})]
    table = pd.read_csv(rf.out_dir / _table_file(rf.image_id)).set_index("nucleus_id")
    foci = run_stack.data[1].astype(np.float64)
    for i in table.index:
        assert table.loc[i, "foci_value"] == pytest.approx(float(foci[run_labels == i].mean()))


@pytest.mark.usefixtures("fixture_post_stage")
def test_the_spacing_comes_from_the_table_on_every_route(tmp_path, monkeypatch):
    """Review S1: the spacing a post-stage sees, so its numbers and its sub-hash, is the image_summary table's
    with a stack handed over (a batch pass that keeps it, D34) and without one (a post-only rerun reads the image
    again), even when the stack says otherwise. A single-plane run records voxel_dz_um as NaN, and dz is then
    PLANAR_DZ (1.0), what that run's own array code received; a table without voxel sizes gives None on both."""
    from germquant import pipeline
    from germquant.io.nd2_reader import Stack

    rf = make_result_folder(tmp_path)
    other = Stack(data=rf.stack.data, spacing=(1.5, 0.4, 0.4), channel_names=list(rf.stack.channel_names),
                  path=rf.stack.path, spacing_ok=True)
    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda p, **kw: {"channel_names": list(other.channel_names)})
    monkeypatch.setattr(pipeline, "read_stack", lambda p, **kw: other)
    summary = rf.out_dir / f"{rf.image_id}__image_summary.csv"

    def routes(dz, dy, dx):
        df = pd.read_csv(summary)
        df["voxel_dz_um"], df["voxel_dy_um"], df["voxel_dx_um"] = dz, dy, dx
        df.to_csv(summary, index=False)
        seen = []
        for stack in (other, None):
            res = _run(rf, [{"name": N}], stack=stack, resume=False)
            assert res["ran"] == [N], res["post"]
            seen.append((res["post"]["stages"][N], pd.read_csv(rf.out_dir / _table_file(rf.image_id))))
        (a, table), (b, again) = seen
        assert a["sub_hash"] == b["sub_hash"] and a["inputs"] == b["inputs"], (a, b)
        pd.testing.assert_frame_equal(table, again)
        return a["inputs"]["spacing"], table

    spacing, _table = routes(*SPACING)
    assert spacing == list(SPACING), "the table's spacing, never the stack's"
    nan = float("nan")
    spacing, table = routes(nan, SPACING[1], SPACING[2])
    assert spacing == [1.0, SPACING[1], SPACING[2]], "a single-plane run: dz is PLANAR_DZ"
    assert table["volume_um3"].to_numpy() == pytest.approx(table["n_voxels"].to_numpy() * SPACING[1] * SPACING[2])
    spacing, table = routes(nan, nan, nan)
    assert spacing is None and table["volume_um3"].isna().all(), "no recorded spacing on either route"


@pytest.mark.usefixtures("fixture_post_stage")
def test_unfinished_series_and_foreign_results_are_refused(tmp_path):
    """The runner's own guards, behind the entry points' refusals: no completion marker (the main run did not
    finish), a time-series timepoint (refused in T5, D33) and a folder made with another config_hash are each
    recorded as a failure of every stage, and nothing else is written."""
    cases = {}
    rf = make_result_folder(tmp_path / "series")
    edit_marker(rf.out_dir, rf.image_id, t_index=3, t_seconds=None)
    cases["series"] = (rf, "time series")
    rf = make_result_folder(tmp_path / "foreign")
    edit_marker(rf.out_dir, rf.image_id, config_hash="000000000000")
    cases["foreign"] = (rf, "config_hash 000000000000")
    rf = make_result_folder(tmp_path / "unfinished")
    (rf.out_dir / f"{rf.image_id}__done.json").unlink()
    cases["unfinished"] = (rf, "no completion marker")
    for what, (rf, fragment) in cases.items():
        before = snapshot(rf.out_dir)
        res = _run(rf, [{"name": N}])
        assert res["failed"] == [N] and res["ran"] == [], what
        st = post_record(rf.out_dir, rf.image_id)["stages"][N]
        assert st["done"] is False and fragment in st["error"], (what, st["error"])
        assert set(snapshot(rf.out_dir)) - set(before) == {f"{rf.image_id}__post.json"}, what
    assert PS.SERIES_REFUSAL in post_record(cases["series"][0].out_dir, "img_001")["stages"][N]["error"]


def test_run_post_stages_imports_the_pipeline_first(tmp_path):
    """D28: in a fresh interpreter, the first module run_post_stages imports is germquant.pipeline (torch before
    CuPy), before the post-stage's own module and before any reader (tifffile for the labels, bioio for an
    OME-TIFF). The import starts are recorded by a meta-path finder, because sys.modules is ordered by the end of
    each import."""
    import tifffile

    rf = make_result_folder(tmp_path)
    site = T.write_fixture_distribution(tmp_path / "site")
    image = tmp_path / "img_001.ome.tif"
    um = chr(0xB5) + "m"
    tifffile.imwrite(image, rf.stack.data, ome=True, photometric="minisblack",
                     metadata={"axes": "CZYX", "Channel": {"Name": ["DAPI", "foci"]}, "PhysicalSizeZ": SPACING[0],
                               "PhysicalSizeY": SPACING[1], "PhysicalSizeX": SPACING[2], "PhysicalSizeZUnit": um,
                               "PhysicalSizeYUnit": um, "PhysicalSizeXUnit": um})
    code = ("import importlib.abc, json, sys\n"
            "started = []\n"
            "class Spy(importlib.abc.MetaPathFinder):\n"
            "    def find_spec(self, name, path=None, target=None):\n"
            "        started.append(name)\n"
            "        return None\n"
            "sys.meta_path.insert(0, Spy())\n"
            "import germquant.post_stages as PS\n"
            "from germquant.config import load_config\n"
            f"cfg = load_config({str(GENERIC)!r})\n"
            "mark = len(started)\n"
            f"res = PS.run_post_stages({str(rf.out_dir)!r}, {rf.image_id!r}, [{{'name': {N!r}}}], cfg=cfg,"
            f" image_path={str(image)!r})\n"
            f"print(json.dumps({{'after': started[mark:], 'ran': res['ran'], 'post': res['post']}}))\n")
    r = run_child(["-c", code], pythonpath=[site])
    assert r.returncode == 0, r.stderr[-3000:]
    out = json.loads(r.stdout.strip().splitlines()[-1])
    after = out["after"]
    assert after and after[0] == "germquant.pipeline", after[:10]
    first = {m: after.index(m) for m in (T.FIXTURE_PACKAGE, "tifffile", "bioio", "nd2", "dask") if m in after}
    assert T.FIXTURE_PACKAGE in first and "tifffile" in first, first
    assert all(i > 0 for i in first.values()), first
    if importlib.util.find_spec("bioio") is not None and importlib.util.find_spec("bioio_ome_tiff") is not None:
        assert out["ran"] == [N], out["post"]
        assert "bioio" in first, first
