"""``germquant sc-straighten``: run the stage on finished results (inputs -> straighten_image -> new files)."""
from __future__ import annotations

import hashlib
import logging
import os
import time
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from ... import provenance
from ...fsutil import long_path
from . import VERSION, inputs, outputs
from . import schema as S
from .params import Params
from .pipeline import straighten_image

log = logging.getLogger(__name__)
HEAD_BYTES = 64 * 1024 * 1024


@dataclass
class Outcome:
    image_id: str
    status: str                       # done | skipped | refused | failed | dry_run
    code: str = ""
    message: str = ""
    out_dir: str = ""
    line: str = ""
    record: dict = field(default_factory=dict)


def _mirror(found_folder: Path, root: Path, out: Path | None) -> Path:
    if out is None:
        return found_folder
    try:
        rel = found_folder.resolve().relative_to(root.resolve())
    except ValueError:
        rel = Path(found_folder.name)
    return out / rel


def _imported_for(inp, labels, spacing, paths_from: list[str] | None, verified_default=False,
                  allow_sealed: bool = False) -> tuple[dict, list]:
    """{nucleus_id: [{"points", "verified", "id"}]} from --paths-from specs, and the file records."""
    from . import imports

    out: dict[int, list] = {}
    recs = []
    for spec in paths_from or []:
        kind, _, path = spec.partition(":")
        paths = imports.read_spec(kind, path, image_id=inp.image_id, spacing=spacing, allow_sealed=allow_sealed)
        recs.append({"spec": spec, "sha256": inputs.file_digest(path, HEAD_BYTES if kind == "ims" else None),
                     "n_paths": len(paths)})
        for ip in paths:
            nid = imports.assign_nucleus(ip.points, labels, spacing)
            if nid is None:
                continue
            out.setdefault(int(nid), []).append({"points": ip.points, "verified": bool(ip.verified or verified_default),
                                                 "id": ip.path_id})
    return out, recs


def code_digest() -> str:
    """sha256 (12 hex) over the sources of this package (line endings normalised), so the sub-hash and every row
    change with the code, not only with VERSION and the parameters."""
    h = hashlib.sha256()
    for f in sorted(Path(__file__).resolve().parent.glob("*.py")):
        h.update(f.name.encode())
        h.update(f.read_bytes().replace(b"\r\n", b"\n"))
    return h.hexdigest()[:12]


def _audit(path, rec: dict) -> None:
    """Append one line to the sealed-run audit log (``--allow-sealed``)."""
    import json

    from ...benchmark import sc_split

    line = {"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "git_sha": provenance.git_sha(),
            "split_sha256": sc_split.REGISTERED_SHA256, **rec}
    p = Path(path)
    os.makedirs(long_path(p.parent), exist_ok=True)
    with open(long_path(p), "a", encoding="utf-8") as fh:
        fh.write(json.dumps(line, default=str) + "\n")


def run_results(results_root: str | Path, image_ids: list[str] | None = None, *, params: Params,
                out: str | Path | None = None, geometry: str = "auto", paths_from: list[str] | None = None,
                refine_imported: bool = True, ends: str = "fit", nuclei: str = "germline", nucleus_ids=None,
                workers: int | None = None, force: bool = False, strips: bool | None = None,
                overview: bool | None = None, source_map=None, allow_sealed: bool = False,
                sealed_reason: str | None = None, audit_log: str | Path | None = None, dry_run: bool = False,
                reader=None, plane_reader=None, echo=print, progress=None) -> list[Outcome]:
    """Run sc-straighten on every finished image under ``results_root`` (or the named ones). ``reader(path, roles,
    xy_stride, z_range)`` overrides the pixel reader and ``plane_reader`` the per-plane reader of the other roles
    (tests). ``allow_sealed`` runs images of the sealed SC validation split; it needs ``sealed_reason`` and writes
    a start and an end line per sealed image to ``audit_log`` (default: ``sc_straighten_sealed_runs.jsonl`` in the
    output root). Returns one Outcome per image."""
    root = Path(results_root)
    out_root = Path(out) if out is not None else None
    if allow_sealed and not (sealed_reason and str(sealed_reason).strip()):
        raise ValueError("allow_sealed needs a sealed_reason (protocol section 10: every sealed run is logged)")
    audit_log = Path(audit_log) if audit_log is not None else (out_root or root) / "sc_straighten_sealed_runs.jsonl"
    found = inputs.select(inputs.discover(root), image_ids)
    if image_ids and not found:
        echo(f"no finished image under {root} matches {image_ids}")
    workers = int(workers if workers is not None else params["workers"])
    write_strips = bool(params["strips.write"] if strips is None else strips)
    write_overview = bool(params["overview"] if overview is None else overview)
    custom_reader = reader is not None          # a test reader brings its own pixels: no source file needed
    if not custom_reader and plane_reader is None:
        plane_reader = inputs.role_planes
    reader = reader or inputs.read_roles
    code = code_digest()
    results = []
    for f in found:
        t0 = time.perf_counter()
        oc = Outcome(image_id=f.image_id, status="failed")
        ctx = {"sealed": None}
        try:
            oc = _one(f, root, out_root, params, geometry, paths_from, refine_imported, ends, nuclei, nucleus_ids,
                      workers, force, write_strips, write_overview, source_map, allow_sealed, dry_run, reader,
                      progress, t0, need_source=not (dry_run or custom_reader), plane_reader=plane_reader,
                      code=code, sealed_reason=sealed_reason, audit_log=audit_log, ctx=ctx)
        except inputs.InputRefused as e:
            oc = Outcome(image_id=f.image_id, status="refused", code=e.code, message=e.message,
                         line=f"{f.image_id}: REFUSED ({e.code}) {e.message}")
        except Exception as e:  # noqa: BLE001 - one image must never stop the others
            log.exception("sc-straighten failed on %s", f.image_id)
            oc = Outcome(image_id=f.image_id, status="failed", code=type(e).__name__, message=str(e),
                         line=f"{f.image_id}: FAILED ({type(e).__name__}: {e})")
            try:
                od = _mirror(f.folder, root, out_root)
                if inputs.golden_root(od) is None:
                    os.makedirs(long_path(od), exist_ok=True)
                    outputs.write_record(outputs.paths_for(od, f.image_id)["record"], {
                        "format": "germquant.sc_straighten", "version": VERSION, "image_id": f.image_id,
                        "status": "failed", "done": False, "error": f"{type(e).__name__}: {e}",
                        "finished": provenance.run_timestamp()})
            except Exception:  # noqa: BLE001
                pass
        if ctx.get("sealed"):
            _audit(audit_log, {"event": "end", "image_id": f.image_id, "role": ctx["sealed"], "status": oc.status,
                               "code": oc.code, "reason": sealed_reason, "params_hash": params.hash})
        echo(oc.line)
        results.append(oc)
    return results


def _one(f, root, out_root, params, geometry, paths_from, refine_imported, ends, nuclei_mode, nucleus_ids,
         workers, force, write_strips, write_overview, source_map, allow_sealed, dry_run, reader, progress, t0,
         need_source=True, plane_reader=None, code="", sealed_reason=None, audit_log=None, ctx=None):
    ctx = ctx if ctx is not None else {}
    inp = inputs.resolve(f, source_map=source_map, need_source=need_source)
    role = inputs.sealed_role(inp.image_id, inp.file_path)
    if role and not allow_sealed:
        raise inputs.InputRefused("sealed", f"{inp.image_id} is a {role} image of the sealed SC validation split "
                                            "(docs/SC_VALIDATION_PROTOCOL.md); it runs only under --allow-sealed "
                                            "with --reason, which the scoring (scripts/benchmark_sc.py zones "
                                            "--final) passes")
    sealed = role if (role and allow_sealed) else None
    out_dir = _mirror(inp.folder, root, out_root)
    g = inputs.golden_root(out_dir)
    if g is not None:
        raise inputs.InputRefused("golden_folder", f"{out_dir} is inside the golden tree {g}; write elsewhere with --out")
    ids = inputs.select_nuclei(inp.nuclei, nuclei_mode, nucleus_ids)
    paths = outputs.paths_for(out_dir, inp.image_id)
    digests = {"nuclei_csv": inputs.file_digest(f.nuclei_path), "labels_tif": inputs.file_digest(inp.labels_path),
               "spots": {k: inputs.file_digest(v) for k, v in sorted(inp.spot_paths.items())}}
    src_stat = None
    if Path(long_path(inp.source_path)).exists():
        st = os.stat(long_path(inp.source_path))
        src_stat = {"path": str(inp.source_path), "size": st.st_size, "mtime": st.st_mtime,
                    "head_sha256": inputs.file_digest(inp.source_path, HEAD_BYTES)}
    payload = {"version": VERSION, "params_hash": params.hash, "tools": outputs.tool_versions(),
               "inputs": digests, "source": {k: src_stat[k] for k in ("size", "head_sha256")} if src_stat else None,
               "geometry": {"xy_stride": inp.xy_stride, "z_range": inp.z_range}, "mode": geometry,
               "paths_from": list(paths_from or []), "refine_imported": bool(refine_imported), "ends": ends,
               "nuclei": nuclei_mode, "nucleus_ids": ids if nucleus_ids else None,
               "strips": bool(write_strips), "code": code}
    sub = outputs.sub_hash(payload)
    prev = outputs.read_record(paths["record"])
    if prev and prev.get("done") and prev.get("sub_hash") == sub and not force:
        return Outcome(image_id=inp.image_id, status="skipped", code="resume", out_dir=str(out_dir),
                       line=f"{inp.image_id}: up to date (sub_hash {sub}), skipped; --force reruns", record=prev)
    if dry_run:
        return Outcome(image_id=inp.image_id, status="dry_run", out_dir=str(out_dir),
                       line=f"{inp.image_id}: would trace {len(ids)} nuclei -> {out_dir} (sub_hash {sub})")
    os.makedirs(long_path(out_dir), exist_ok=True)
    ctx["sealed"] = sealed
    if sealed:
        _audit(audit_log, {"event": "start", "image_id": inp.image_id, "role": sealed, "reason": sealed_reason,
                           "params_hash": params.hash, "sub_hash": sub, "out_dir": str(out_dir)})
    # pixels: the central_element channel whole (tracing, the noise model, the overview); the other roles only as
    # per-nucleus crops, read plane by plane (the whole channel is never in memory)
    cmap = {k: v for k, v in inp.channel_map.items() if v is not None}
    want = ["central_element"] + [r for r in list(params["profile.roles"]) + list(params["strips.roles"])
                                  if r != "central_element"]
    roles_idx = {r: cmap[r] for r in dict.fromkeys(want) if r in cmap}
    arrays = reader(inp.source_path, {"central_element": roles_idx["central_element"]}, xy_stride=inp.xy_stride,
                    z_range=inp.z_range)
    ce = arrays.pop("central_element")
    max_label = int(pd.to_numeric(inp.nuclei["nucleus_id"], errors="coerce").max())
    labels = inputs.read_labels(inp.labels_path, max_label)
    if tuple(labels.shape) != tuple(ce.shape):
        raise inputs.InputRefused("labels_mismatch", f"{inp.image_id}: labels {labels.shape} != image {ce.shape}")
    from .crops import boxes as _boxes
    from .crops import crop_slices

    bx = _boxes(labels)
    others = {r: i for r, i in roles_idx.items() if r != "central_element"}
    role_crops = {}
    if others:
        if plane_reader is not None:
            plane_fn = plane_reader(inp.source_path, others, xy_stride=inp.xy_stride, z_range=inp.z_range)
        else:
            full = reader(inp.source_path, others, xy_stride=inp.xy_stride, z_range=inp.z_range)

            def plane_fn(r, z, full=full):
                return full[r][int(z)]
        slices = [crop_slices(bx[i], inp.spacing, params["crop_pad_um"], labels.shape) for i in ids if i in bx]
        role_crops = inputs.crop_banks(plane_fn, list(others), labels.shape, slices)
    spots = inputs.read_spots(inp.spot_paths)
    imported, imp_recs = ({}, [])
    if geometry in ("import", "both"):
        imported, imp_recs = _imported_for(inp, labels, inp.spacing, paths_from, allow_sealed=bool(sealed))
    det_max = inputs.detector_max(inp, ce)
    wavelengths = inputs.emission_wavelengths(inp.source_path, roles_idx) if need_source else None
    res = straighten_image(ce, labels, inp.spacing, inp.nuclei, params=params, roles=role_crops,
                           spot_tables=spots, nucleus_ids=ids, imported=imported, geometry=geometry,
                           refine_imported=refine_imported, ends_mode=ends, workers=workers, detector_max=det_max,
                           image_id=inp.image_id, progress=progress, label_boxes=bx, wavelengths=wavelengths)
    # outputs
    from .qc import image_flags

    flags = image_flags({**res.record["flag_info"], "geometry_source": inp.geometry_source}, params)
    meta = {k: v for k, v in inp.meta.items()}
    meta.update({"sc_straighten_version": VERSION, "sc_straighten_hash": sub})
    strands = res.strands.copy()
    if len(strands):
        strands["marker"] = _marker_of(inp)
    outputs.write_table(strands, paths["strands"], S.STRANDS, meta)
    outputs.write_table(res.samples, paths["samples"], S.SAMPLES)
    outputs.write_table(res.gaps, paths["gaps"], S.GAPS, meta)
    outputs.write_table(res.foci, paths["foci"], S.FOCI, meta)
    outputs.write_table(res.nuclei, paths["nuclei"], S.NUCLEI, meta)
    if write_strips:
        desc = {"format": "germquant.sc_straighten.strips", "version": VERSION, "roles": res.strip_roles,
                "du_um": float(params["strips.du_um"]), "dv_um": float(params["strips.dv_um"]),
                "ds_um": float(params["strips.ds_um"]), "s_offset_um": -float(params["strips.extend_um"]),
                "half_um": list(params["strips.half_um"]), "image_id": inp.image_id, "sub_hash": sub,
                "index": res.strip_index}
        outputs.write_strips(paths["strips_flat"], res.strips_flat, dict(desc, kind="flat", rows="u",
                                                                         flat_v_um=float(params["strips.flat_v_um"])))
        outputs.write_strips(paths["strips_axial"], res.strips_axial, dict(desc, kind="axial", rows="v",
                                                                           axial_u_um=float(params["strips.axial_u_um"])))
    if write_overview:
        try:
            outputs.overview_png(paths["overview"], ce, strands, res.samples, res.nuclei, inp.nuclei, inp.spacing,
                                 inp.image_id)
        except Exception as e:  # noqa: BLE001 - a picture never fails the image
            log.warning("overview PNG failed for %s: %s", inp.image_id, e)
    elapsed = time.perf_counter() - t0
    c = res.record["counts"]
    record = {
        "format": "germquant.sc_straighten", "version": VERSION, "image_id": inp.image_id, "status": "done",
        "done": True, "sub_hash": sub, "params_hash": params.hash, "params": params.as_dict(),
        "inputs": {"results_dir": str(inp.folder), "nuclei_csv": {"path": str(f.nuclei_path), **{"sha256": digests["nuclei_csv"]}},
                   "labels_tif": {"path": str(inp.labels_path), "sha256": digests["labels_tif"]},
                   "spot_tables": {k: {"path": str(v), "sha256": digests["spots"][k]} for k, v in inp.spot_paths.items()},
                   "source": src_stat, "file_path_recorded": inp.file_path, "xy_stride": inp.xy_stride,
                   "z_range": inp.z_range, "geometry_source": inp.geometry_source, "channel_map": inp.channel_map,
                   "roles_read": sorted(roles_idx), "imported": imp_recs, "nuclei_mode": nuclei_mode,
                   "nucleus_ids": ids if nucleus_ids else None, "geometry": geometry,
                   "refine_imported": refine_imported, "ends": ends, "detector_max": det_max,
                   "emission_nm": wavelengths},
        "tool_versions": outputs.tool_versions(), "git_sha": provenance.git_sha(), "code_sha256": code,
        **{k: res.record[k] for k in ("noise_model", "null_q", "sigma_vox", "psf", "complete_median_integrated",
                                      "complete_median_fwhm", "n_cut_across_nuclei", "n_contacts_freed", "counts",
                                      "nucleus_errors", "timing")},
        "role_crops_mb": round(sum(b.nbytes for b in role_crops.values()) / 2 ** 20, 1),
        "flags": flags, "elapsed_s": round(elapsed, 1), "workers": int(workers),
        "sealed_override": bool(sealed), "sealed_role": sealed, "sealed_reason": sealed_reason if sealed else None,
        "outputs": {k: v.name for k, v in paths.items() if v.exists() or k == "record"},
        "finished": provenance.run_timestamp(),
    }
    outputs.write_record(paths["record"], record)
    n_nuc = c.get("nuclei", 0)
    st_txt = ", ".join(f"{k[7:]} {v}" for k, v in c.items() if k.startswith("status_"))
    cl_txt = ", ".join(f"{k[6:]} {v}" for k, v in c.items() if k.startswith("class_"))
    flag_txt = ", ".join(x for x in flags if x != "sc_straighten:uncalibrated") or "none"
    line = (f"{inp.image_id}: {n_nuc} nuclei ({st_txt}), {c.get('strands', 0)} strands ({cl_txt}), "
            f"flags: {flag_txt}, {elapsed / 60:.1f} min")
    return Outcome(image_id=inp.image_id, status="done", out_dir=str(out_dir), line=line, record=record)


def _marker_of(inp) -> str:
    """The central_element marker as the run's channel map names it (io.channel_map of the run's manifest,
    resolved against the repository); the role name when it cannot be found."""
    import json

    try:
        man = json.loads((inp.folder / "run_manifest.json").read_text(encoding="utf-8"))
        rel = man["config"]["io"]["channel_map"]
        p = Path(rel)
        if not p.is_absolute():
            root = provenance._repo_root()
            p = (root / rel) if root is not None else p
        from ...io.channel_map import ChannelMap

        return ChannelMap.from_yaml(p).marker("central_element")
    except Exception:  # noqa: BLE001 - a label only
        return "central_element"
