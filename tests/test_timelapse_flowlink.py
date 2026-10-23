"""The flow-prior linking backend (germquant.timelapse.flowlink, `germquant track-link --backend laptrack_flow`;
docs/TIMELAPSE.md "Flow-prior linking").

Three kinds of checks. The pieces (the settings, the mutual-nearest-neighbour links, the median field, the local
consensus, a step without enough links, predictions across a missing timepoint) on hand-made points. The linking on
synthetic nuclei with known identities: a quasi-random sheet of about 170 nuclei (dart throwing, at least 4 um
apart) under a smooth shear (the distal third still, the proximal third sliding 4 um per step, more than the
3 um laptrack cutoff and nearly the spacing) with a random walk of 0.12 um per axis per step and 4 percent
single-frame misses: plain laptrack joins hundreds of wrong pairs, the flow prior none; then the two series
folders of test_timelapse_piecewise (the 4.6 um proximal shear) and test_timelapse_tracking (drift, a dropout,
off-gonad nuclei) through link_series and the command line. And a regression check that every existing backend
gives, on those two folders, the tables it gave at origin/main 8eb73de, before this backend existed.
"""
from __future__ import annotations

import hashlib
import json

import numpy as np
import pandas as pd
import pytest
from test_timelapse_piecewise import PERIOD_S as PW_PERIOD_S
from test_timelapse_piecewise import STEM as PW_STEM
from test_timelapse_piecewise import _identity_errors, flow_per_step, write_piecewise_series
from test_timelapse_tracking import DROPOUT, _check_identities, _identity_of_points, write_synthetic_series

from germquant import cli
from germquant.timelapse import link_series, register_series, register_series_piecewise
from germquant.timelapse import flowlink as FL
from germquant.timelapse import track as TR

# ---- synthetic nuclei with known identities ------------------------------------------------------------------------


def _sheet(rng, n_target=170, size=(30.0, 120.0), min_d=4.0) -> np.ndarray:
    """(z, y, x) um of quasi-random nuclei on a thin sheet: dart throwing with a minimum spacing (no lattice)."""
    pts: list[np.ndarray] = []
    tries = 0
    while len(pts) < n_target and tries < 40000:
        tries += 1
        p = np.array([rng.uniform(2.6, 3.4), rng.uniform(0.0, size[0]), rng.uniform(0.0, size[1])])
        if all(np.hypot(p[1] - q[1], p[2] - q[2]) >= min_d for q in pts):
            pts.append(p)
    return np.asarray(pts)


def shear_objects(seed: int = 3, *, n_t: int = 7, shear_um: float = 4.0, walk_um: float = 0.12, p_miss: float = 0.04,
                  drop_timepoint: int | None = None):
    """An object table (t_index, nucleus_id, zr_um, yr_um, xr_um) of the sheared sheet and ``ids[t] = {nucleus_id:
    identity}``. Each step moves a nucleus by shear_um x smoothstep((x - 40) / 40) along x plus a random walk; a
    nucleus misses an inner frame with probability p_miss; nucleus ids are shuffled per frame. ``drop_timepoint``
    removes a whole timepoint (a failed frame of the series)."""
    rng = np.random.default_rng(seed)
    base = _sheet(rng)
    n = len(base)
    pos = base.copy()
    rows, ids = [], {}
    for t in range(n_t):
        if t:
            s = np.clip((pos[:, 2] - 40.0) / 40.0, 0.0, 1.0)
            pos = pos + np.column_stack([np.zeros(n), np.zeros(n), shear_um * s * s * (3 - 2 * s)])
            pos = pos + rng.normal(0.0, walk_um, (n, 3))
        present = np.ones(n, bool)
        if 0 < t < n_t - 1:
            present &= rng.uniform(size=n) >= p_miss
        if t == drop_timepoint:
            continue
        idx = np.flatnonzero(present)
        nid = rng.permutation(len(idx)) + 1
        ids[t] = {int(a): int(b) for a, b in zip(nid, idx)}
        rows += [{"t_index": t, "nucleus_id": int(a), "zr_um": pos[b, 0], "yr_um": pos[b, 1], "xr_um": pos[b, 2]}
                 for a, b in zip(nid, idx)]
    return pd.DataFrame(rows), ids


def link_errors(assign: dict, ids: dict) -> dict:
    """wrong: consecutive points of one track that are different identities; mixed: tracks holding more than one
    identity; split: identities spread over more than one track."""
    by_track: dict[int, list] = {}
    by_ident: dict[int, set] = {}
    for (t, nid), tid in assign.items():
        by_track.setdefault(tid, []).append((t, ids[t][nid]))
        by_ident.setdefault(ids[t][nid], set()).add(tid)
    wrong = mixed = 0
    for lst in by_track.values():
        idents = [i for _t, i in sorted(lst)]
        wrong += sum(1 for a, b in zip(idents[:-1], idents[1:]) if a != b)
        mixed += int(len(set(idents)) > 1)
    return {"wrong": wrong, "mixed": mixed, "split": sum(1 for v in by_ident.values() if len(v) > 1),
            "tracks": len(by_track)}


@pytest.fixture(scope="module")
def sheared():
    obj, ids = shear_objects()
    return {"objects": obj, "ids": ids}


# ---- the settings -----------------------------------------------------------------------------------------------

def test_flow_settings_defaults_keys_and_validation():
    s = FL.flow_settings(None)
    assert s == {"search_um": 6.0, "link_radius_um": 3.0, "ratio": 0.5, "k": 8, "iterations": 2, "cutoff_um": None,
                 "fb": False}
    assert set(FL.KEYS.values()) == {f"timelapse.track.flow_{n}" for n in FL.DEFAULTS}
    cfg = {"timelapse": {"track": {"flow_k": 6, "flow_cutoff_um": 2.5, "flow_fb": True}}}
    s = FL.flow_settings(cfg)
    assert s["k"] == 6 and s["cutoff_um"] == 2.5 and s["fb"] is True and s["ratio"] == 0.5
    s = FL.flow_settings(cfg, k=10, fb=False)                     # an explicit value wins over the key
    assert s["k"] == 10 and s["fb"] is False and s["cutoff_um"] == 2.5
    assert FL.flow_settings(None, search_um=0)["search_um"] == 0.0          # 0: start from the registration
    bad = [{"ratio": 0}, {"ratio": 1.5}, {"ratio": "x"}, {"k": 0}, {"k": 2.5}, {"k": True}, {"iterations": 0},
           {"search_um": -1}, {"search_um": float("nan")}, {"link_radius_um": 0}, {"cutoff_um": 0}, {"cutoff_um": -2},
           {"fb": "maybe"}, {"fb": 2}]
    for kw in bad:
        with pytest.raises(ValueError, match="timelapse.track.flow_"):
            FL.flow_settings(None, **kw)
    with pytest.raises(ValueError, match="unknown flow setting"):
        FL.flow_settings(None, radius=2.0)
    with pytest.raises(ValueError, match="timelapse.track.flow_ratio"):
        FL.flow_settings({"timelapse": {"track": {"flow_ratio": 2}}})


# ---- the pieces -------------------------------------------------------------------------------------------------

def test_mutual_links_radius_and_ratio_test():
    P = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 10.0], [0.0, 0.0, 30.0]])
    B = np.array([[0.0, 0.0, 0.5], [0.0, 0.0, 10.4], [0.0, 0.0, 11.0], [0.0, 0.0, 31.2], [0.0, 0.0, 32.0]])
    i, j = FL.mutual_links(P, B, radius_um=3.0, ratio=0.5)
    assert list(zip(i, j)) == [(0, 0), (1, 1)]          # 30 -> 31.2 has a rival at 32.0 (ratio 1.2 / 2.0 = 0.6)
    i, j = FL.mutual_links(P, B, radius_um=3.0, ratio=0.7)
    assert list(zip(i, j)) == [(0, 0), (1, 1), (2, 3)]
    i, j = FL.mutual_links(P, B, radius_um=0.45, ratio=0.5)
    assert list(zip(i, j)) == [(1, 1)]                  # 0.5 um is beyond the radius
    i, j = FL.mutual_links(P, B, radius_um=3.0, ratio=0.3)
    assert list(zip(i, j)) == [(0, 0)]                  # 0.4 / 1.0 fails a ratio of 0.3
    # not mutual: B[0] is nearest to both P points, P[1] is nearer to it
    i, j = FL.mutual_links(np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.3]]), np.array([[0.0, 0.0, 0.4]]), radius_um=3, ratio=1)
    assert list(zip(i, j)) == [(1, 0)]
    # a lone pair passes the ratio test; a non-finite point is never linked
    i, j = FL.mutual_links(np.array([[0.0, 0, 0], [np.nan, 0, 0]]), np.array([[0.0, 0, 1.0]]), radius_um=3, ratio=0.5)
    assert list(zip(i, j)) == [(0, 0)]
    assert len(FL.mutual_links(np.zeros((0, 3)), B, radius_um=3, ratio=0.5)[0]) == 0


def test_flow_field_is_the_median_of_the_nearest_links():
    S = np.array([[0.0, 0, 0], [0.0, 0, 1], [0.0, 0, 2], [0.0, 0, 50]])
    W = np.array([[0.0, 0, 1.0], [0.0, 0, 1.2], [0.0, 0, 9.0], [0.0, 3.0, 0]])
    f = FL.FlowField(S, W, k=3)
    np.testing.assert_allclose(f([[0.0, 0, 0.9]]), [[0.0, 0.0, 1.2]])       # median of 1.0, 1.2, 9.0: the outlier is out
    np.testing.assert_allclose(f([[0.0, 0, 49.0]]), [[0.0, 0.0, 1.2]])      # nearest three: 50, 2, 1
    out = f([[np.nan, 0, 0], [0.0, 0, 0.0]])
    assert np.isnan(out[0]).all() and np.isfinite(out[1]).all()
    assert FL.FlowField(S[:1], W[:1], k=8)([[5.0, 5, 5]]).tolist() == [[0.0, 0.0, 1.0]]     # fewer samples than k
    assert FL.FlowField.zero()([[1.0, 2, 3]]).tolist() == [[0.0, 0.0, 0.0]] and FL.FlowField.zero().n_samples == 0


def test_local_consensus_finds_a_shift_larger_than_the_spacing():
    rng = np.random.default_rng(11)
    A = _sheet(rng, n_target=120, size=(30.0, 60.0))
    shift = np.array([0.0, 1.5, 4.5])                    # 4.7 um, more than the 4 um spacing
    B = A + shift + rng.normal(0, 0.1, A.shape)
    u, support = FL.consensus_displacement(A, B, search_um=6.0, k=8)
    found = support > 0
    assert found.mean() > 0.9
    assert np.abs(u[found] - shift).max() < 0.4
    # out of reach (the true partners lie beyond the search radius): the shift cannot be found, and what some wrong
    # partners happen to agree on has a much weaker support (the search radius must cover the motion the
    # registration leaves; docs/TIMELAPSE.md)
    u_far, support_far = FL.consensus_displacement(A, B, search_um=3.0, k=8)
    assert not (np.abs(u_far[support_far > 0] - shift).max(axis=1) < 1.0).any()
    assert support_far.mean() < 0.5 * support.mean()
    assert FL.consensus_displacement(A, B[:0], search_um=6.0, k=8)[1].sum() == 0
    assert FL.consensus_displacement(A, B, search_um=0.0, k=8)[1].sum() == 0


def test_a_step_without_enough_links_keeps_the_registration():
    A = np.array([[0.0, 0, 0], [0.0, 0, 5], [0.0, 5, 0]])
    st = FL.estimate_step(A, A + 0.1, FL.flow_settings(None), t_from=3, t_to=4)
    assert st.info["fallback"].startswith("3 confident link(s), fewer than 8") and st.forward.n_samples == 0
    assert np.allclose(st.forward(A), 0.0) and (st.t_from, st.t_to) == (3, 4)
    st = FL.estimate_step(A, np.zeros((0, 3)), FL.flow_settings(None))
    assert "without nuclei" in st.info["fallback"]


def test_predictions_follow_the_steps_and_span_a_missing_timepoint():
    # 2.5 um per step: the 5 um across the missing timepoint stays inside the 6 um consensus search
    obj, _ids = shear_objects(drop_timepoint=3, p_miss=0.0, shear_um=2.5)
    assert sorted(obj["t_index"].unique()) == [0, 1, 2, 4, 5, 6]
    flows = FL.estimate_flow(obj, FL.flow_settings(None))
    assert sorted(flows) == [0, 1, 2, 4, 5] and flows[2].t_to == 4               # the step 2 -> 4 spans the gap
    P = FL.predict(obj, flows, 2)
    t = obj["t_index"].to_numpy()
    assert np.isnan(P[t == 2, 0]).all() and np.isfinite(P[t == 2, 1]).all()      # nothing at 3; 4 is reached
    assert np.isfinite(P[t == 1, 0]).all() and np.isnan(P[t == 1, 1]).all()      # 2 is reached, 3 is not on disk
    assert np.isfinite(P[t == 4]).all() and np.isnan(P[t == 6]).all()
    # two steps are the composition of the two step fields
    x0 = obj.loc[t == 0, ["zr_um", "yr_um", "xr_um"]].to_numpy(float)
    x1 = x0 + flows[0].forward(x0)
    np.testing.assert_allclose(P[t == 0, 0], x1)
    np.testing.assert_allclose(P[t == 0, 1], x1 + flows[1].forward(x1))
    Q = FL.predict(obj, flows, 2, backward=True)
    assert np.isnan(Q[t == 0]).all() and np.isfinite(Q[t == 4, 1]).all() and np.isnan(Q[t == 4, 0]).all()


# ---- linking: the shear laptrack cannot follow ------------------------------------------------------------------

def test_flow_prior_recovers_identities_that_laptrack_swaps(sheared):
    pytest.importorskip("laptrack")    # the timelapse extra; CI installs it
    obj, ids = sheared["objects"], sheared["ids"]
    plain, _ = TR._link_laptrack(obj, max_dist_um=3.0, gap_frames=2)
    e_plain = link_errors(plain, ids)
    assert e_plain["wrong"] >= 100 and e_plain["mixed"] >= 50, e_plain        # the sliding nuclei take their neighbours
    flow, info, extra = FL.link_flow(obj, max_dist_um=3.0, gap_frames=2, settings=FL.flow_settings(None))
    e_flow = link_errors(flow, ids)
    assert e_flow["wrong"] == 0 and e_flow["mixed"] == 0, e_flow
    assert e_flow["split"] <= 2                        # a miss in two consecutive frames cannot be bridged with gap 2
    assert len(flow) == len(obj) and info["n_steps"] == 6 and info["n_steps_registration_only"] == 0
    assert info["cutoff_um"] == 3.0 and info["cutoff_source"] == "max_dist_um" and info["residual_median_um"] < 0.5
    # the predicted step of a nucleus is the planted shear at its position
    s = np.clip((obj["xr_um"].to_numpy() - 40.0) / 40.0, 0, 1)
    planted = 4.0 * s * s * (3 - 2 * s)
    last = obj["t_index"].to_numpy() == obj["t_index"].max()
    assert np.isnan(extra.loc[last, "flow_step_x_um"]).all()
    err = np.abs(extra.loc[~last, "flow_step_x_um"].to_numpy() - planted[~last])
    assert np.median(err) < 0.25 and np.percentile(err, 95) < 0.8
    # without the consensus (the conservative first pass alone) the 4 um slide is out of reach, as for laptrack
    nv, _info, _ = FL.link_flow(obj, max_dist_um=3.0, gap_frames=2, settings=FL.flow_settings(None, search_um=0.0))
    assert link_errors(nv, ids)["wrong"] >= 100
    # deterministic: the same objects give the same tracks, and the rows in another order the same partition of the
    # nuclei into tracks (the track numbers may differ)
    again, _, _ = FL.link_flow(obj, max_dist_um=3.0, gap_frames=2, settings=FL.flow_settings(None))
    assert again == flow
    shuffled = obj.sample(frac=1.0, random_state=5).reset_index(drop=True)
    other, _, _ = FL.link_flow(shuffled, max_dist_um=3.0, gap_frames=2, settings=FL.flow_settings(None))

    def partition(assign):
        groups: dict[int, set] = {}
        for key, tid in assign.items():
            groups.setdefault(tid, set()).add(key)
        return {frozenset(g) for g in groups.values()}

    assert partition(other) == partition(flow)


def test_forward_backward_cost_and_an_explicit_cutoff(sheared):
    pytest.importorskip("laptrack")
    obj, ids = sheared["objects"], sheared["ids"]
    fb, info, _ = FL.link_flow(obj, max_dist_um=3.0, gap_frames=2, settings=FL.flow_settings(None, fb=True))
    assert link_errors(fb, ids)["wrong"] == 0 and info["flow"]["fb"] is True
    tight, info, _ = FL.link_flow(obj, max_dist_um=3.0, gap_frames=2, settings=FL.flow_settings(None, cutoff_um=0.3))
    e = link_errors(tight, ids)
    assert info["cutoff_um"] == 0.3 and info["cutoff_source"] == "flow_cutoff_um"
    assert e["wrong"] == 0 and e["tracks"] >= len(ids[0]) + 100      # a cutoff near the random walk cuts the tracks
    no_gap, _, _ = FL.link_flow(obj, max_dist_um=3.0, gap_frames=0, settings=FL.flow_settings(None))
    assert link_errors(no_gap, ids)["split"] > link_errors(
        FL.link_flow(obj, max_dist_um=3.0, gap_frames=2, settings=FL.flow_settings(None))[0], ids)["split"]


def test_a_missing_timepoint_is_bridged_by_gap_closing():
    pytest.importorskip("laptrack")
    obj, ids = shear_objects(drop_timepoint=3, p_miss=0.0, shear_um=2.5)
    flow, info, _ = FL.link_flow(obj, max_dist_um=3.0, gap_frames=2, settings=FL.flow_settings(None))
    e = link_errors(flow, ids)
    assert e["wrong"] == 0 and e["split"] == 0 and e["tracks"] == len(ids[0])
    assert [s["t_to"] - s["t_from"] for s in info["steps"]] == [1, 1, 2, 1, 1]


# ---- link_series and the command line on the series folders --------------------------------------------------------

@pytest.fixture(scope="module")
def pw_series(tmp_path_factory):
    return write_piecewise_series(tmp_path_factory.mktemp("flow_pw"))


def test_link_series_follows_the_shear_without_the_piecewise_registration(pw_series):
    pytest.importorskip("laptrack")
    out, ids = pw_series["out"], pw_series["ids"]
    plain = link_series(out, backend="laptrack", write=False)            # no shift table: the raw frame
    assert plain["record"]["shift_source"] == "none" and _identity_errors(plain["points"], ids)[1] >= 8
    res = link_series(out, backend="laptrack_flow")
    rec = res["record"]
    assert rec["backend"] == "laptrack_flow" and rec["link_source"] == "laptrack_flow" and rec["note"] == ""
    assert _identity_errors(res["points"], ids) == (0, 0)
    assert rec["n_tracks"] == 60 and rec["n_full_length"] == 60
    bi = rec["backend_info"]
    assert bi["cutoff_um"] == 3.0 and bi["gap_frames"] == 2 and bi["flow"]["k"] == 8 and len(bi["steps"]) == 4
    assert bi["n_steps_registration_only"] == 0 and all(s["n_confident"] >= 50 for s in bi["steps"])
    pts = res["points"]
    assert (pts["link_source"] == "laptrack_flow").all()
    # the predicted step is the planted flow; the table's registered frame is the registration's (zero here)
    first = pts[pts["t_index"] < pts["t_index"].max()]
    ident = _identity_of_points(first, ids)
    planted = np.array([flow_per_step(int(i)) for i in ident])
    err = np.abs(first["flow_step_x_um"].to_numpy() - planted)
    assert np.median(err) < 0.1 and err.max() < 0.5             # the ends of a row see neighbours on one side only
    np.testing.assert_allclose(pts["xr_um"], pts["x_um"])
    # the files: the flow columns ride along as extras, the record says what ran
    tp = pd.read_csv(out / f"{PW_STEM}__track_points.csv")
    assert {"flow_step_z_um", "flow_step_y_um", "flow_step_x_um"} <= set(tp.columns)
    on_disk = json.loads((out / f"{PW_STEM}__track_link.json").read_text())
    assert on_disk["backend"] == "laptrack_flow" and on_disk["backend_info"]["flow"]["search_um"] == 6.0
    tracks = res["tracks"].set_index("time_track_id")
    v = tracks["axis_velocity_raw_um_per_min"].to_numpy(float)
    t0 = pts[pts["t_index"] == 0].set_index("time_track_id")
    v_flow = np.array([flow_per_step(ids[0][int(t0.loc[t, "nucleus_id"])]) / (PW_PERIOD_S / 60.0) for t in tracks.index])
    np.testing.assert_allclose(v, v_flow, atol=0.01)
    # on top of the piecewise registration it is as clean (the residual flow is small there)
    register_series_piecewise(out, max_step_um=6.0)
    res_pw = link_series(out, backend="laptrack_flow", registration="piecewise", write=False)
    assert _identity_errors(res_pw["points"], ids) == (0, 0) and res_pw["record"]["n_full_length"] == 60
    assert res_pw["record"]["registration"] == "piecewise"


def test_link_series_flow_on_the_tracking_series(tmp_path):
    pytest.importorskip("laptrack")
    s = write_synthetic_series(tmp_path)
    res = link_series(s["out"], backend="laptrack_flow", write=False)
    linked = _check_identities(res["points"], s["ids"])               # one track per identity, the dropout bridged
    assert res["record"]["n_tracks"] == 18 and res["record"]["n_full_length"] == 17
    linked["identity"] = _identity_of_points(linked, s["ids"])
    drop = linked[linked["identity"] == DROPOUT[0]].sort_values("t_index")
    assert list(drop["t_index"]) == [0, 1, 3, 4] and list(drop["is_gap"]) == [False, False, True, False]
    # 18 germline nuclei: the steps have enough links; with every label (the two off-gonad nuclei) too
    assert res["record"]["backend_info"]["n_steps_registration_only"] == 0
    every = link_series(s["out"], backend="laptrack_flow", germline_only=False, write=False)
    assert every["record"]["n_tracks"] == 20
    # flow settings are refused with another backend, and bad values before anything is read
    with pytest.raises(ValueError, match="laptrack_flow backend only"):
        link_series(s["out"], backend="laptrack", flow={"k": 4}, write=False)
    with pytest.raises(ValueError, match="timelapse.track.flow_k"):
        link_series(s["out"], backend="laptrack_flow", flow={"k": 0}, write=False)
    assert link_series(s["out"], backend="laptrack", flow={"k": None}, write=False)["record"]["backend"] == "laptrack"
    # the config's keys reach the backend (here from a dict, as the series' run_manifest.json supplies one)
    cfg = {"timelapse": {"track": {"flow_k": 5, "flow_fb": True, "max_dist_um": 2.5}}}
    bi = link_series(s["out"], backend="laptrack_flow", cfg=cfg, write=False)["record"]["backend_info"]
    assert bi["flow"]["k"] == 5 and bi["flow"]["fb"] is True and bi["cutoff_um"] == 2.5


def test_cli_track_link_flow_options(tmp_path, capsys):
    pytest.importorskip("laptrack")
    s = write_synthetic_series(tmp_path)
    sd = str(s["out"])
    assert cli.main(["track-link", sd, "--backend", "laptrack_flow", "--flow-k", "6", "--flow-fb",
                     "--flow-search-um", "5", "--flow-iterations", "3", "--flow-ratio", "0.6",
                     "--flow-link-radius-um", "2.5", "--flow-cutoff-um", "2.8"]) == 0
    out = capsys.readouterr().out
    assert "backend laptrack_flow" in out and "18 track(s)" in out and "flow prior: 4 step(s)" in out
    assert "cutoff 2.80 um" in out
    rec = json.loads((s["out"] / f"{s['out'].name}__track_link.json").read_text())
    want = {"k": 6, "fb": True, "search_um": 5.0, "iterations": 3, "ratio": 0.6, "link_radius_um": 2.5, "cutoff_um": 2.8}
    assert {k: rec["backend_info"]["flow"][k] for k in want} == want
    # the defaults when no flag is given (the run_manifest.json of the folder sets max_dist_um 2.5: the cutoff)
    assert cli.main(["track-link", sd, "--backend", "laptrack_flow"]) == 0
    rec = json.loads((s["out"] / f"{s['out'].name}__track_link.json").read_text())
    assert rec["backend_info"]["flow"]["k"] == 8 and rec["backend_info"]["cutoff_um"] == 2.5
    assert rec["backend_info"]["flow"]["fb"] is False
    # a key in the series' config is used; --no-flow-fb overrides it
    (s["out"] / "run_manifest.json").write_text(json.dumps(
        {"config": {"timelapse": {"track": {"max_dist_um": 2.5, "flow_fb": True, "flow_k": 7}}}}))
    assert cli.main(["track-link", sd, "--backend", "laptrack_flow"]) == 0
    rec = json.loads((s["out"] / f"{s['out'].name}__track_link.json").read_text())
    assert rec["backend_info"]["flow"]["fb"] is True and rec["backend_info"]["flow"]["k"] == 7
    assert cli.main(["track-link", sd, "--backend", "laptrack_flow", "--no-flow-fb"]) == 0
    rec = json.loads((s["out"] / f"{s['out'].name}__track_link.json").read_text())
    assert rec["backend_info"]["flow"]["fb"] is False and rec["backend_info"]["flow"]["k"] == 7
    capsys.readouterr()
    # a flow option with another backend is refused, naming the option
    assert cli.main(["track-link", sd, "--backend", "laptrack", "--flow-k", "6"]) == 1
    assert "--flow-k apply to --backend laptrack_flow only" in capsys.readouterr().err
    assert cli.main(["track-link", sd, "--backend", "overlap_inhouse", "--flow-fb"]) == 1
    capsys.readouterr()
    # a value out of range names the option and its key
    assert cli.main(["track-link", sd, "--backend", "laptrack_flow", "--flow-ratio", "1.5"]) == 1
    assert "--flow-ratio (timelapse.track.flow_ratio) must be in (0, 1]" in capsys.readouterr().err
    assert cli.main(["track-link", sd, "--backend", "laptrack_flow", "--flow-k", "0"]) == 1
    with pytest.raises(SystemExit):
        cli.main(["track-link", sd, "--backend", "laptrack_flow", "--flow-k", "two"])     # argparse: not an integer


def test_the_backend_list_is_appended():
    assert TR.BACKENDS[:4] == ("overlap", "overlap_inhouse", "laptrack", "none") and TR.BACKENDS[-1] == "laptrack_flow"
    assert not any(k.startswith("flow") for k in TR.DEFAULTS)      # the backend's keys live in flowlink.DEFAULTS


# ---- the existing backends are unchanged -----------------------------------------------------------------------------

def _canonical_digest(df: pd.DataFrame, sort_cols) -> str:
    """sha256 of a table as text: rows sorted by ``sort_cols``, the columns in their order, floats rounded to 6
    decimals (NaN as 'nan', -0 as 0), so the same tables give the same digest on every CI platform."""
    d = df.sort_values(list(sort_cols), kind="stable").reset_index(drop=True)
    cols = []
    for c in d.columns:
        s = d[c]
        if pd.api.types.is_float_dtype(s):
            v = np.round(s.to_numpy(float), 6) + 0.0
            cols.append(["nan" if not np.isfinite(x) else f"{x:.6f}" for x in v])
        else:
            cols.append(["" if (x is None or (not isinstance(x, str) and pd.isna(x))) else str(x) for x in s.tolist()])
    lines = [",".join(str(c) for c in d.columns)] + [",".join(row) for row in zip(*cols)]
    return hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()


# Computed at origin/main 8eb73de (before the laptrack_flow backend) with this file's helpers on the two fixtures:
# the digests of the points table, the tracks table and the (t_index, nucleus_id, time_track_id) assignment.
EXPECTED_8EB73DE = {
    "tracking/overlap_inhouse_stride1": ("ccf0f0ce72537de0", "ee64e1fc22ecba48", "d5ef5c7e6420cf14", 19, 17),
    "tracking/overlap_inhouse": ("ccf0f0ce72537de0", "ee64e1fc22ecba48", "d5ef5c7e6420cf14", 19, 17),
    "tracking/laptrack": ("f0f7798c27313fad", "78c43396c3691284", "2fdb9e8273561a67", 18, 17),
    "tracking/laptrack_gap0": ("5cabb2f1bac50bad", "8edbce0ebbeb31bb", "109a008f736b3618", 19, 17),
    "tracking/laptrack_all_nuclei": ("28e23ac01055993a", "7707f93c8d373b32", "24afd3ce0ef66545", 20, 19),
    "tracking/none": ("2f454b3c264bb177", "ef6f8bc5f5bb13ae", "a49e64a9c3f58830", 0, 0),
    "tracking/overlap_stride1": ("ccf0f0ce72537de0", "2ddb7388bd4db754", "03fcd0336ea47f87", 19, 17),
    "piecewise/laptrack_rigid": ("d885bee631071110", "7e3d1d694e213c60", "efcfdf4c67467caa", 68, 52),
    "piecewise/laptrack_piecewise": ("ba886ffbcb850f03", "94a7035a88b19153", "eb6506a980a6909e", 60, 60),
    "piecewise/overlap_inhouse_piecewise": ("ba886ffbcb850f03", "cd38eadfff17d6a1", "f2dc4a064caeb46e", 60, 60),
    "piecewise/overlap_inhouse_rigid": ("16b751f42ff45c58", "b95eb0221c39159f", "be3c252ea1bdb602", 174, 26),
}


def test_existing_backends_give_the_tables_they_gave_before(tmp_path):
    pytest.importorskip("laptrack")
    s = write_synthetic_series(tmp_path / "tracking")
    runs = {"tracking/overlap_inhouse_stride1": (s["out"], dict(backend="overlap_inhouse", stride=1)),
            "tracking/overlap_inhouse": (s["out"], dict(backend="overlap_inhouse")),
            "tracking/laptrack": (s["out"], dict(backend="laptrack")),
            "tracking/laptrack_gap0": (s["out"], dict(backend="laptrack", gap=0)),
            "tracking/laptrack_all_nuclei": (s["out"], dict(backend="laptrack", germline_only=False)),
            "tracking/none": (s["out"], dict(backend="none"))}
    if TR.cellacdc_overlap_available()[0]:
        runs["tracking/overlap_stride1"] = (s["out"], dict(backend="overlap", stride=1))
    p = write_piecewise_series(tmp_path / "pw")
    register_series_piecewise(p["out"], max_step_um=6.0)
    register_series(p["out"], max_step_um=6.0)
    runs.update({"piecewise/laptrack_rigid": (p["out"], dict(backend="laptrack", gap=2, max_dist_um=3.0)),
                 "piecewise/laptrack_piecewise": (p["out"], dict(backend="laptrack", gap=2, max_dist_um=3.0,
                                                                  registration="piecewise")),
                 "piecewise/overlap_inhouse_piecewise": (p["out"], dict(backend="overlap_inhouse", stride=1,
                                                                         registration="piecewise")),
                 "piecewise/overlap_inhouse_rigid": (p["out"], dict(backend="overlap_inhouse", stride=1))})
    for name, (folder, kw) in runs.items():
        res = link_series(folder, write=False, **kw)
        pts, trk = res["points"], res["tracks"]
        got = (_canonical_digest(pts[["t_index", "nucleus_id", "time_track_id"]], ["t_index", "nucleus_id"])[:16],
               _canonical_digest(pts, ["t_index", "nucleus_id"])[:16], _canonical_digest(trk, ["time_track_id"])[:16],
               int(res["record"]["n_tracks"]), int(res["record"]["n_full_length"]))
        assert got == EXPECTED_8EB73DE[name], (name, got)
        assert "flow_step_x_um" not in pts.columns
    assert set(runs) <= set(EXPECTED_8EB73DE) and len(runs) >= len(EXPECTED_8EB73DE) - 1     # overlap: Cell-ACDC
