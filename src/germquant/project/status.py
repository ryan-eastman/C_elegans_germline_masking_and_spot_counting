"""Project status (docs/PROJECTS.md, "Status"; design 7.7; PLAN.md 9.3, decisions D12, D18, D44, D46).

``project_status(project)`` returns the record `germquant project status --json` prints (PLAN.md section 9.3):

    {"format": 1, "root", "name", "config_hash", "expected_config_hash", "edited", "validated",
     "pack": {"id", "pinned", "installed", "project_copy"},
     "models": [{"key", "ref", "sha256", "present", "verified", "mismatch", "needed"}],
     "images": [{"id", "path", "status", "reason_code", "reason", "out_dir", "qc_pass", "last_run_id",
                 "device_used", "error"}],
     "warnings": [str], "lock": {"held", "host", "pid", "since", "stale"},
     "runs": [{"run_id", "via", "started", "ended", "counts", "interrupted", "finished"}]}

An image is ``done`` when batch.is_done accepts its results folder under the project's recipe and run settings,
with the segment sub-hashes a run compares (project/run.py, ResumeHashes); ``stale`` with the reason is_done
gives (reason_code config_hash, stages, geometry, planar_flag or segment_hash); ``failed`` when its marker lists
failed stages (reason_code failed_stages) or, without results, when the last run that processed it failed (its
error); ``unreachable`` when it is not done and its file cannot be reached (a NAS or VPN not connected);
``excluded``; else ``missing`` (for a time series that the recipe cannot run as one, the reason adds the refusal
every run fails it with, run.series_refusal). ``path`` is the file a run reads (a relative manifest path joined
to the project root, an absolute or UNC path as it is). qc_pass comes from its image summary, device_used from
its stage record, last_run_id from the last run that processed it (done or failed).

``config_hash`` is computed now (open_project); ``pack.installed`` lists the versions of the pinned pack installed
in this environment outside the project, ``project_copy`` whether packs/<id>/<version>/ holds the project's copy.
``models`` is project/run.py's model check, the one a run refuses on. The warnings are open_project's (a hash
that changed and why, the project's copy of the pack missing or replaced, an unrecorded recipe.yaml, a held or
stale lock), a model the run needs that is missing or not the expected file, a segmentation model the engine
cannot resolve on this computer (said once, in the run's own words where a run would refuse, the classical
method included; status reports such a model and never fails on it, and where the segment hash needs the model
the images' segment stage is not compared), the labels_from rule of D12, and an environment mix: runs that name
more than one device, interpreter or engine git sha. ``runs`` merges every log under runs/ by time; a run whose
log has no run_end (it was stopped without one, or is still running in another process) is ``interrupted`` with
``finished`` false.
"""
from __future__ import annotations

import os
from pathlib import Path

from ..fsutil import long_path
from . import runlog
from .manifest import ProjectError

FORMAT = 1
# batch.is_done's reasons, by the words they start with, and the reason codes of PLAN.md 9.3
REASON_CODES = (("config_hash", "config_hash"), ("enabled stages", "stages"), ("xy_stride", "geometry"),
                ("z_range", "geometry"), ("t_index", "geometry"), ("stage(s) failed", "failed_stages"),
                ("planar flag", "planar_flag"), ("segment stage hash", "segment_hash"))


def reason_code(why: str) -> str | None:
    """The 9.3 reason code of a reason batch.is_done gives, None for done and for no marker."""
    text = str(why or "")
    return next((code for prefix, code in REASON_CODES if text.startswith(prefix)), None)


def _version_sort(versions) -> list[str]:
    try:
        from packaging.version import Version

        return sorted(versions, key=Version)
    except Exception:  # noqa: BLE001 - a version that is not PEP 440 sorts as text
        return sorted(versions)


def _last_events(events: list[dict]) -> tuple[dict, dict]:
    """``(the last image event per image id, the run id that last processed each image)``."""
    last, processed = {}, {}
    for e in events:
        kind = e.get("event")
        if kind in runlog.IMAGE_EVENTS and e.get("image_id") is not None:
            iid = str(e["image_id"])
            last[iid] = e
            if kind in ("image_done", "image_failed"):
                processed[iid] = e.get("run_id")
    return last, processed


def _runs(events: list[dict]) -> list[dict]:
    order: list[str] = []
    by_run: dict[str, list[dict]] = {}
    for e in events:
        rid = str(e.get("run_id"))
        if rid not in by_run:
            by_run[rid] = []
            order.append(rid)
        by_run[rid].append(e)
    out = []
    for rid in order:
        ev = by_run[rid]
        start = next((e for e in ev if e.get("event") == "run_start"), {})
        end = next((e for e in reversed(ev) if e.get("event") == "run_end"), None)
        tally = {"done": 0, "failed": 0, "skipped": 0}
        for e in ev:
            if e.get("event") in runlog.IMAGE_EVENTS:
                key = str(e["event"]).removeprefix("image_")
                tally[key] = tally.get(key, 0) + 1
        out.append({"run_id": rid, "via": start.get("via"), "started": start.get("time") or ev[0].get("time"),
                    "ended": end.get("time") if end else None,
                    "counts": dict(end.get("counts") or tally) if end else tally,
                    "interrupted": bool(end.get("interrupted")) if end else True, "finished": end is not None})
    return out


def _environment_mix(events: list[dict]) -> str | None:
    starts = [e for e in events if e.get("event") == "run_start"]
    seen = {"devices": sorted({str(s["device_selected"]) for s in starts if s.get("device_selected")}),
            "interpreters": sorted({str(s["interpreter"]) for s in starts if s.get("interpreter")}),
            "engine git shas": sorted({str((s.get("engine") or {}).get("git_sha")) for s in starts
                                       if (s.get("engine") or {}).get("git_sha")})}
    parts = [f"{what} {', '.join(vals)}" for what, vals in seen.items() if len(vals) > 1]
    if not parts:
        return None
    return ("environment mix: the project's runs used " + "; ".join(parts) + ". Results made on another device or "
            "with another engine can differ; rerun the images in one environment to compare them")


def _image(entry: dict, ctx: dict) -> dict:
    from .run import RESULTS_DIR, finished_fields, image_path, is_series, marker_id, reachable, unreachable_text

    root = ctx["root"]
    iid = str(entry["id"])
    path = image_path(root, entry)
    out = Path(root) / RESULTS_DIR / iid
    rec = {"id": iid, "path": str(path), "status": None, "reason_code": None, "reason": None, "out_dir": str(out),
           "qc_pass": None, "last_run_id": ctx["processed"].get(iid), "device_used": None, "error": None}
    if entry.get("excluded"):
        rec.update(status="excluded", reason=f"excluded: {entry.get('reason') or 'no reason given'}")
        return rec
    mid = marker_id(path)
    geo = {"xy_stride": ctx["xy_stride"], "z_range": ctx["z_range"]}
    if is_series(entry):
        whys = []
        for k in range(int(entry["timepoints"])):
            done, why = ctx["hashes"].is_done(out / f"t{k:03d}", f"{mid}_t{k:03d}", path, t=k, **geo)
            if not done:
                whys.append((k, why))
        if not whys:
            from .run import summary_qc_pass

            passes = [summary_qc_pass(out / f"t{k:03d}", f"{mid}_t{k:03d}") for k in range(int(entry["timepoints"]))]
            rec.update(status="done", reason="done", qc_pass=all(bool(p) for p in passes))
            return rec
        stale = [(k, w) for k, w in whys if w != "no completion marker"]
        done, why = False, (f"t{stale[0][0]:03d}: {stale[0][1]}" if stale else "no completion marker")
        code = reason_code(stale[0][1]) if stale else None
        fields = {}
    else:
        done, why = ctx["hashes"].is_done(out, mid, path, **geo)
        code = reason_code(why)
        fields = finished_fields(out, mid) if why != "no completion marker" else {}
        if done:
            rec.update(status="done", reason="done", qc_pass=fields.get("qc_pass"),
                       device_used=fields.get("segment_device_used"))
            return rec
    if not reachable(path):
        rec.update(status="unreachable", reason=unreachable_text(path))
        return rec
    if why != "no completion marker":
        rec.update(status="failed" if code == "failed_stages" else "stale", reason_code=code, reason=why,
                   qc_pass=fields.get("qc_pass"), device_used=fields.get("segment_device_used"))
        return rec
    last = ctx["last"].get(iid)
    if last is not None and last.get("event") == "image_failed":
        rec.update(status="failed", reason=f"the last run failed: {last.get('error')}", error=last.get("error"))
        return rec
    rec.update(status="missing", reason="no results yet")
    if is_series(entry):
        refused = _series_refused(entry, path, ctx)
        if refused:
            rec["reason"] = f"no results yet, and a run fails it: {refused}"
    return rec


def _series_refused(entry: dict, path, ctx: dict) -> str | None:
    """The refusal every run fails this series with (run.series_refusal), in the run's words, or None."""
    from .run import series_refusal, series_refusal_text

    try:
        why = series_refusal(ctx["cfg"], path)
    except Exception:  # noqa: BLE001 - io.planar auto reads the header; a run names that error itself
        return None
    return series_refusal_text(entry, why) if why else None


def _resolution_warning(cfg, manifest: dict, models: list[dict], problem: str | None) -> str | None:
    """The warning for a segmentation model the engine cannot resolve on this computer, said once (None when a
    needed model is already reported missing or mismatched). Where a run makes the check (no labels_from, and a
    fixed image to run: process_image resolves the model under every method, classical included), it is the
    run's own refusal; ``problem`` is ResumeHashes', whose segment stage comparison is then left out."""
    from .run import LABELS_KEY, _model_resolution_refusal, is_series

    if any(r["needed"] and (not r["present"] or r["mismatch"]) for r in models):
        return None
    fixed = any(not is_series(e) and not e.get("excluded") for e in manifest.get("images") or [])
    refusal = _model_resolution_refusal(cfg) if fixed and cfg.get(LABELS_KEY) is None else None
    if refusal is None and problem is None:
        return None
    text = refusal or (f"the recipe's segmentation model cannot be resolved here ({problem}): its images fail in "
                       f"a run until the model is back")
    if problem is not None:
        text += "; until then the images' segment stage is not compared"
    return text


def project_status(project) -> dict:
    """The status record of the module docstring (PLAN.md 9.3). ProjectError only when the project's recipe does
    not load (open_project names why); a model that is not on this computer is reported, never raised."""
    import germquant.pipeline  # noqa: F401  (torch first: the model lookup and the backends' provenance follow)

    from ..packs import discover, use_project
    from . import _root, open_project, project_config
    from .lock import lock_report
    from .recipe import resolve_pack
    from .run import ResumeHashes, mismatch_text, model_records

    root = _root(project)
    proj = open_project(root)
    m = proj.manifest
    if proj.config_hash is None:
        why = "; ".join(proj.warnings) or "the recipe does not load"
        raise ProjectError(f"the status of {root} cannot be computed: {why}")
    warnings = list(proj.warnings)
    events = runlog.read_logs(root)
    with use_project(root):
        cfg = project_config(proj)
        pack, recipe = resolve_pack(root, m)
        models, model_warnings = model_records(root, m, cfg, pack, recipe)
    for r in models:
        if r["needed"] and not r["present"]:
            warnings.append(f"model {r['ref']} ({r['key']}) is missing: a run is refused until it is back")
        elif r["needed"] and r["mismatch"]:
            warnings.append(f"model {r['ref']} ({r['key']}) {mismatch_text(r, pack)}; a run is refused")
    warnings += model_warnings
    mix = _environment_mix(events)
    if mix:
        warnings.append(mix)
    last, processed = _last_events(events)
    zr = m["run"].get("z_range")
    ctx = {"root": root, "cfg": cfg, "hashes": ResumeHashes(cfg, unresolved_ok=True),
           "xy_stride": int(m["run"].get("xy_stride") or 1),
           "z_range": (int(zr[0]), int(zr[1])) if zr else None, "last": last, "processed": processed}
    with use_project(root):
        images = [_image(e, ctx) for e in m.get("images") or []]
        unresolved = _resolution_warning(cfg, m, models, ctx["hashes"].problem)
    if unresolved:
        warnings.append(unresolved)
    p = m["pack"]
    installed = _version_sort({pk.version for pk in discover().packs if pk.id == p["id"]})
    return {"format": FORMAT, "root": str(root), "name": m.get("name"), "config_hash": proj.config_hash,
            "expected_config_hash": m["recipe"].get("expected_config_hash"), "edited": bool(m["recipe"].get("file")),
            "validated": m.get("validated"),
            "pack": {"id": p["id"], "pinned": f"{p['id']}@{p['version']}", "installed": installed,
                     "project_copy": os.path.isfile(long_path(Path(root) / "packs" / p["id"] / p["version"] /
                                                             "pack.yaml"))},
            "models": models, "images": images, "warnings": warnings, "lock": lock_report(root),
            "runs": _runs(events)}
