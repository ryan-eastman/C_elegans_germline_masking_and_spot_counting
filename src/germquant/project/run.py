"""Running a project (docs/PROJECTS.md, "Running a project"; design 7.5 to 7.8; PLAN.md 4.4 P5-RUN, decisions D9,
D12, D15, D17, D21, D28, D33, D40, D44, D46 to D49).

``run_project(project, images=None, resume=True, device=None, shard=None, progress=None, *, via="cli", run_id=None,
should_stop=None)`` runs the project's images into results/<image_id>/, one api.run_loaded call per image (the call
`germquant run` makes, so a project run writes what `germquant run --pack` writes) and an image with timepoints
through timelapse.process_series, and writes the run's log, runs/<run_id>.jsonl (project/runlog.py). It returns
``{"run_id", "log", "counts": {"done", "failed", "skipped"}, "exit_code", "interrupted", "warnings", "message"}``.

1. germquant.pipeline is imported first, so torch claims cuBLAS before any reader or backend can load CuPy (D28).
2. Refused before anything is written (ProjectError; the command's exit 2): a malformed run id, via, device or
   shard; an image named in ``images`` that the project does not hold, or excludes; a folder that is no project;
   a recipe that does not load; a pack that the pack check finds unusable (its model blockers aside: the model
   check below knows which models this project's recipe uses); a model the run needs that is missing, or whose
   resolved file has another sha256 than the manifest records (design 7.5: a file at the same relative path under
   $GERMQUANT_MODELS wins over the project's own), checked before any pixel is read; the sample sheet the manifest
   names (``samples``) when it cannot be opened or read (read_project_sheet); the project lock held by an edit in
   progress (a lock older than 12 h is a warning: it is never broken silently); a log of that run id that exists
   already. In a labels_from run the segmentation model is not loaded, so its check is skipped, with a warning
   (D12). The sheet is read once, and its warnings (germquant.io.sample_sheet: a cp1252 file, a non-ISO date, a
   germ_cell no count can be made from) join the run's warnings, each said once.
3. run_start, then each image in manifest order (``images`` picks some; ``shard`` "K/N" keeps every N-th of them
   from the K-th on, so the N shards cover the images exactly once): should_stop is asked before each image; an
   image already done is skipped when ``resume`` (batch.is_done, with the segment sub-hashes cli._batch computes:
   the planar one only for a recipe that sets io.planar, D44, and the per-image labels hash in labels_from runs,
   D46; in a project with a sheet, the image's row, as cli._batch compares it with --samples); an image that cannot
   be reached, or fails (in its resume check or, for a series, the header read that io.planar auto decides from, as
   well as in its run), is an image_failed event and the run goes on. In a project with a sheet each fixed image
   runs with its row (SampleSheet.row of its image id and the path the run reads), as `batch --samples` runs it:
   its design columns replace the file name's, and an image the sheet has no row for gets the empty row, which
   flags it metadata:not_in_sample_sheet. A project without a sheet makes the call it always made (D40).
4. run_end, with interrupted true when should_stop stopped the run between images, or when the run stopped part
   way (Ctrl-C, or an error that is no image's own): the image in progress is then an image_failed event, and
   an error raises RunError, which names the run id and the log.

An image with timepoints that this release cannot run as a series is refused before process_series is called,
which would resolve and hash the segmentation model before its loop (D33): a recipe with labels_from (D15), one
that runs single-plane images (io.planar true, or auto on a one-plane header; D9), a series that has a row in the
project's sample sheet (D21: the time-lapse loop gives its timepoints no row), or any series while post-stages are
switched on (design 5.6). A series without a row runs as before, its design columns from the file name.

Post-stages (docs/POST_STAGES.md, design 5.6): the project's own post_stages.yaml, next to recipe.yaml, replaces
the post_stages/<recipe>.yaml its pack ships, and without either none run. The list in force is read and checked
against this environment's post-stages once, before the run log exists: a list that cannot be read, or that holds
an entry that cannot run here, refuses the run (ProjectError).
Each fixed image then runs its post-stages through api.run_loaded(post_stages=), right after its main run; an image
that resume skips runs the due ones only (a changed parameter, file or input; a post-only rerun reads the image
again when a post-stage reads it). A post-stage's failure is recorded in the image's __post.json, never in the run
log: the image stays done.

GERMQUANT_DEVICE is set from ``device`` (the caller's choice, auto included, D47), else from the manifest's
run.device, for the run, and restored afterwards; with neither it is left as it is. The run settings xy_stride
and z_range come from the manifest. ``progress(event)`` is called with every log event as written. A run never
writes the manifest and takes no lock, so shards and Snakemake jobs run side by side.
"""
from __future__ import annotations

import contextlib
import logging
import os
import platform
import re
import socket
import sys
import time
from contextvars import ContextVar
from pathlib import Path

from ..fsutil import long_path
from .manifest import ProjectError, read_manifest
from .runlog import IMAGE_FIELDS, VIA, RunLog, check_run_id, done_sha256, new_run_id

log = logging.getLogger(__name__)

DEVICES = ("auto", "cuda", "mps", "cpu")
MODEL_KEY = "segmentation.nuclei.cellpose_model"
LABELS_KEY = "segmentation.nuclei.labels_from"
RESULTS_DIR = "results"
_SHARD_RE = re.compile(r"^\s*(\d+)\s*/\s*(\d+)\s*$")
# the command line of `germquant project run` (cli_projects sets it), recorded as run_start's argv; None from the API
RUN_ARGV: ContextVar[list | None] = ContextVar("germquant_project_run_argv", default=None)


class RunError(RuntimeError):
    """A run that stopped part way after its start (not a refusal: its log holds run_start, and run_end says
    interrupted when it could still be written). ``run_id`` and ``log`` name the run and its log."""

    def __init__(self, message: str, *, run_id: str | None = None, log: str | None = None) -> None:
        super().__init__(message)
        self.run_id, self.log = run_id, log


# ---- what to run ----------------------------------------------------------------------------------------------

def parse_shard(shard) -> tuple[int, int] | None:
    """``(K, N)`` from "K/N" (or a pair), 1 <= K <= N; None for None; ProjectError otherwise."""
    if shard is None:
        return None
    if isinstance(shard, (tuple, list)) and len(shard) == 2:
        try:
            k, n = (int(v) for v in shard)
        except (TypeError, ValueError):
            k = n = 0
        ok = all(not isinstance(v, bool) for v in shard)
    else:
        m = _SHARD_RE.fullmatch(str(shard))
        k, n = (int(m.group(1)), int(m.group(2))) if m else (0, 0)
        ok = m is not None
    if not ok or n < 1 or not 1 <= k <= n:
        raise ProjectError(f"shard {shard!r} is not K/N with 1 <= K <= N, such as 2/3 (the second of three shards)")
    return k, n


def select_images(manifest: dict, images=None, shard=None) -> list[dict]:
    """The manifest's image entries a run takes, in manifest order: every image not excluded, or those ``images``
    names (ProjectError for an id the project does not hold or excludes); then, with ``shard`` "K/N", every N-th
    of them from the K-th on."""
    entries = list(manifest.get("images") or [])
    if images is None:
        chosen = [e for e in entries if not e.get("excluded")]
    else:
        want = [str(images)] if isinstance(images, (str, os.PathLike)) else [str(i) for i in images]
        known = {str(e.get("id")): e for e in entries}
        unknown = [i for i in want if i not in known]
        if unknown:
            raise ProjectError(f"image {', '.join(unknown)} is not in the project (its images: "
                               f"{', '.join(known) or 'none'})")
        out = [i for i in want if known[i].get("excluded")]
        if out:
            why = "; ".join(f"{i}: {known[i].get('reason') or 'no reason given'}" for i in out)
            raise ProjectError(f"image {', '.join(out)} is excluded from the project ({why})")
        wanted = set(want)
        chosen = [e for e in entries if str(e.get("id")) in wanted]
    k_n = parse_shard(shard)
    if k_n is not None:
        k, n = k_n
        chosen = chosen[k - 1::n]
    return chosen


def is_series(entry: dict) -> bool:
    """Whether a manifest image has timepoints (project add records them for a header with T > 1)."""
    try:
        return int(entry.get("timepoints") or 1) > 1
    except (TypeError, ValueError):
        return False


def image_path(root, entry: dict) -> Path:
    from .images import resolve_path

    return resolve_path(root, entry["path"])


def reachable(path) -> bool:
    return os.path.exists(long_path(path))


def unreachable_text(path) -> str:
    return f"not reachable: {path} (is the NAS or VPN connected?)"


def series_refusal(cfg, path, *, row=None, post_stages=None) -> str | None:
    """Why an image with timepoints cannot run as a series under this recipe in this release, or None: a recipe
    that reads labels from files (D15), one that runs single-plane images (io.planar true, or auto on a header with
    one plane; D9), ``row``, the series' row of the project's sample sheet, when the sheet has one for it (D21:
    germquant.io.sample_sheet.SERIES_REFUSAL, the message `run --t ... --samples` and the API give), or
    ``post_stages``, the entries of the run's post-stage list, when it switches one on (D33:
    germquant.post_stages.SERIES_REFUSAL). The time-lapse loop is not changed in this release (design 4.13, PLAN.md
    Q1), so such a series must never reach it. `germquant run --t` refuses with these labels and single-plane
    checks too (cli._series_refusal; reviews/final_partA.md N2), so the two routes cannot drift apart."""
    if cfg.get(LABELS_KEY) is not None:
        from ..io.labels import SERIES_REFUSAL

        return SERIES_REFUSAL
    raw = cfg.get("io.planar")
    if raw is not None and raw is not False:          # planar_setting's own reading: a 0 is malformed, refused below
        from .. import planar

        try:
            mode = planar.planar_setting(cfg)
        except planar.PlanarConfigError as e:          # a malformed io.planar: refused here too
            return str(e)
        if mode == "true":
            return planar.PLANAR_SERIES_REFUSED
        if mode == "auto":
            from ..pipeline import read_input_metadata

            if planar.one_plane_header(read_input_metadata(path)):
                return planar.PLANAR_SERIES_REFUSED
    if row is not None and row.row_sha256 is not None:
        from ..io.sample_sheet import SERIES_REFUSAL as SHEET_SERIES_REFUSAL

        return SHEET_SERIES_REFUSAL
    if post_stages:
        from ..post_stages import SERIES_REFUSAL as POST_SERIES_REFUSAL

        return POST_SERIES_REFUSAL
    return None


def series_refusal_text(entry: dict, why: str, root=None) -> str:
    """The error a run logs for a series that series_refusal refuses (project_status shows it before any run). The
    sheet's refusal names --samples, which `project run` does not take, so here it gets the project's own remedy
    (its folder ``root``, else DIR)."""
    text = f"REFUSED: a time series of {int(entry['timepoints'])} timepoints: {why}"
    from ..io.sample_sheet import SERIES_REFUSAL as SHEET_SERIES_REFUSAL

    if why == SHEET_SERIES_REFUSAL:
        where = "DIR" if root is None else f'"{root}"'
        text += (f". In a project: take the row of {entry['id']} out of the project's sample sheet, or stop using the "
                 f"sheet with `germquant project edit {where} --samples none`; the series then runs, its design "
                 f"columns from the file name")
    elif why.startswith("post-stages are refused"):
        where = "DIR" if root is None else f'"{root}"'
        text += (f". In a project the post-stage list (its post_stages.yaml, else the pack's post_stages/<recipe>.yaml) "
                 f"applies to every image: switch the post-stages off with `python -m germquant.post_stages write "
                 f"{where} --entries []`, or run the series outside the project with `germquant run --t`")
    return text


def project_post_stages(root, pack, recipe) -> tuple[str | None, list]:
    """``(path, entries)`` of the post-stage list a project run uses (docs/POST_STAGES.md): the project's own
    post_stages.yaml, else the post_stages/<recipe>.yaml its pack (the project's copy) ships for the recipe
    (germquant.post_stages.post_stages_in_force); ``(None, [])`` without either. The list is read and checked
    against the post-stages of this environment here: a list that cannot be read, or holds an entry that cannot
    run here, is a ProjectError (the run's refusal, exit 2)."""
    from ..post_stages import entries_of, post_stages_in_force, validate_post_stages

    source, _origin = post_stages_in_force(project_root=root, pack_folder=pack.folder, recipe=recipe.id)
    if source is None:
        return None, []
    try:
        entries = entries_of(source)
    except ValueError as e:
        raise ProjectError(str(e)) from None
    except OSError as e:
        raise ProjectError(f"post-stages {Path(source).absolute()}: {e.strerror or e}") from None
    problems = validate_post_stages(entries)
    if problems:
        raise ProjectError(f"post-stages {Path(source).absolute()}: " + "; ".join(problems))
    return str(source), entries


def _post_only(out_dir: Path, mid: str, path, ctx: dict) -> None:
    """The due post-stages of an image whose main run resume skips (germquant.post_stages.run_post_stages: a
    changed parameter, file or input reruns that post-stage only). A failure is the image's __post.json's, and
    nothing here fails the image or the run."""
    from .. import post_stages as PS

    try:
        res = PS.run_post_stages(out_dir, mid, ctx["post"], cfg=ctx["cfg"], image_path=path)
    except Exception:  # noqa: BLE001 - a post-stage never fails the image (design 5.6)
        log.exception("%s: the post-stages could not be run", mid)
        return
    for name in res.get("failed") or []:
        log.warning("%s: post-stage %s failed (see %s)", mid, name, res.get("record"))


def read_project_sheet(root, manifest: dict):
    """The project's sample sheet (germquant.io.sample_sheet.SampleSheet), read once from the file the manifest
    names (``samples``, a path inside the project), or None when it names none. Its warnings stay in ``.warnings``:
    the SampleSheetWarning each is also issued as is held back, so a caller says each once. A sheet that cannot be
    read is a ProjectError with SampleSheetError's message; a file that cannot be opened one that names the sheet
    and the error's strerror, as `germquant batch --samples` does (never the long-path form the file is opened by)."""
    samples = manifest.get("samples")
    if not samples:
        return None
    import warnings

    from ..io.sample_sheet import SampleSheetError, SampleSheetWarning, read_sample_sheet
    from .images import resolve_path

    path = resolve_path(root, str(samples))
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SampleSheetWarning)      # each is in sheet.warnings, said once there
            return read_sample_sheet(path)
    except SampleSheetError as e:
        raise ProjectError(str(e)) from None
    except OSError as e:
        raise ProjectError(f"sample sheet {Path(path).absolute()}: {e.strerror or e}") from None


# ---- the digests a resume compares --------------------------------------------------------------------------

class ResumeHashes:
    """The segment sub-hashes batch.is_done compares, as cli._batch computes them: one pair for the run (the model
    file's sha256 and the segmentation backend's provenance; the planar hash only when the recipe sets io.planar,
    D44), computed on first use, or in a labels_from run the pair of each image's labels file (D46).

    A segmentation model that cannot be resolved here raises, as it would in a run; with ``unresolved_ok`` (project
    status, which reports such a model instead of failing) the reason is kept in ``problem`` and the pair is
    (None, None), so is_done leaves the segment stage out of its comparison (its default without a segment hash)."""

    def __init__(self, cfg, *, unresolved_ok: bool = False) -> None:
        self.cfg = cfg
        self.unresolved_ok = unresolved_ok
        self.problem: str | None = None
        self._pair: tuple[str | None, str | None] | None = None

    def for_image(self, path) -> tuple[str | None, str | None]:
        cfg = self.cfg
        if cfg.get(LABELS_KEY) is not None:
            from ..io.labels import segment_hashes

            seg, seg_planar, _problem = segment_hashes(cfg, path, model_sha="unused")
            return seg, seg_planar
        if self._pair is None:
            from .. import pipeline, provenance
            from ..models_registry import ModelPinError
            from ..stages import stage_hashes

            method = str(cfg.get("segmentation.nuclei.method", "auto"))
            try:
                model_sha = (provenance.file_sha256(pipeline.resolve_model_path(cfg)) if method != "classical"
                             else "unused")
            except (OSError, ValueError, ModelPinError) as e:
                if not self.unresolved_ok:
                    raise
                self.problem = f"{type(e).__name__}: {e}"
                self._pair = (None, None)
                return self._pair
            seg_prov = pipeline.segment_backend_provenance(cfg)
            tools = provenance.tool_versions()
            seg = stage_hashes(cfg, {}, tools, model_sha, segment_backend=seg_prov)["segment"]
            planar = (None if cfg.get("io.planar") in (None, False) else
                      stage_hashes(cfg, {}, tools, model_sha, segment_backend=seg_prov, planar=True)["segment"])
            self._pair = (seg, planar)
        return self._pair

    def is_done(self, out_dir, marker_id: str, path, *, xy_stride: int, z_range, t=None,
                row=None) -> tuple[bool, str]:
        """batch.is_done with this run's segment sub-hashes and, for ``row`` (the image's row of the project's
        sample sheet) when the sheet has one for it, its row_sha256: passed only then, as cli._batch passes it
        (D40, D44's call shape), and None standing for no row on the marker's side (PLAN.md 9.4 #13, reading 2)."""
        from .. import batch

        seg, seg_planar = self.for_image(path)
        kw = {} if seg_planar is None else {"segment_hash_planar": seg_planar}
        if t is not None:
            kw["t"] = t
        if row is not None and row.row_sha256 is not None:
            kw["sample_row_sha256"] = row.row_sha256
        return batch.is_done(out_dir, marker_id, self.cfg, xy_stride=xy_stride, z_range=z_range, segment_hash=seg,
                             **kw)


def marker_id(path) -> str:
    """The id process_image gives an image (and its done marker): the file name without its format suffix."""
    from ..io.formats import image_stem

    return image_stem(Path(path)).strip()


# ---- models (design 7.5, D12) ----------------------------------------------------------------------------------

class _Named:
    """A Config that names ``value`` under ``key``: resolve_model_path's view of a model the manifest records."""

    def __init__(self, cfg, key: str, value) -> None:
        self._cfg, self._key, self._value = cfg, key, value
        self.base_dir = cfg.base_dir

    def get(self, key, default=None):
        return self._value if key == self._key else self._cfg.get(key, default)


def _flat(data: dict, prefix: str = "") -> dict:
    out = {}
    for k, v in (data or {}).items():
        key = f"{prefix}.{k}" if prefix else str(k)
        if isinstance(v, dict) and v:
            out.update(_flat(v, key))
        else:
            out[key] = v
    return out


def _shadow_hint(path) -> str:
    root = os.environ.get("GERMQUANT_MODELS")
    if root:
        try:
            Path(os.path.abspath(path)).relative_to(Path(os.path.abspath(root)))
            return (f" The file lies under $GERMQUANT_MODELS ({root}), which the engine searches first for a relative "
                    f"model path, so it shadows the project's own file: move it away, or share the model through a "
                    f"pack (germquant pack derive)")
        except ValueError:
            pass
    return ""


def model_records(root, manifest: dict, cfg, pack, recipe) -> tuple[list[dict], list[str]]:
    """``([{"key", "ref", "sha256", "present", "verified", "mismatch", "needed"}], warnings)``: the models the
    manifest records (use_model: the file the run resolves, as the run resolves it, with its sha256 compared) and
    the pinned models the pack recipe needs (their blobs, as pack check finds them). ``needed`` is whether this
    run loads the model: its key still names it, and for the segmentation model the recipe is not classical and
    reads no labels from files (D12; such a model's digest is not checked, with a warning)."""
    from .. import pipeline, provenance
    from ..packs.blobs import check_blob

    flat = _flat(cfg.as_dict())
    labels_run = cfg.get(LABELS_KEY) is not None
    seg_loaded = str(cfg.get("segmentation.nuclei.method", "auto")) != "classical" and not labels_run
    out: list[dict] = []
    warnings: list[str] = []
    for e in manifest.get("models") or []:
        key, ref, want = str(e.get("key")), str(e.get("ref")), str(e.get("sha256") or "")
        named = flat.get(key) == ref
        needed = named and (seg_loaded if key == MODEL_KEY else True)
        try:
            path = pipeline.resolve_model_path(_Named(cfg, MODEL_KEY, ref)) if key == MODEL_KEY else Path(root) / ref
        except (OSError, ValueError):
            path = None
        # a plain path on purpose: Cellpose loads the model by it, and use-model refuses a longer one
        present = path is not None and Path(str(path)).is_file()
        got = provenance.file_sha256(path) if present else None
        verified = present and got == want
        mismatch = None
        if present and not verified:
            mismatch = (f"{path} has sha256 {str(got)[:12]}..., the manifest records {want[:12]}..."
                        + _shadow_hint(path))
        out.append({"key": key, "ref": ref, "sha256": want or None, "present": present, "verified": verified,
                    "mismatch": mismatch, "needed": needed})
        if named and key == MODEL_KEY and labels_run:
            warnings.append(f"the recipe reads its nuclei from labels files ({LABELS_KEY}), so the segmentation model "
                            f"{ref} is not needed: its digest is not checked (D12)")
        elif not named:
            warnings.append(f"the manifest records the model {ref} for {key}, but the recipe names "
                            f"{flat.get(key)!r} there, so that model's digest is not checked")
    for name in recipe.need("models"):
        spec = pack.model(name)
        if spec is None:
            continue                                         # a manifest problem, which pack check reports
        keys = sorted(k for k, v in flat.items() if v == name)
        key = keys[0] if keys else None
        needed = bool(keys) and (seg_loaded if key == MODEL_KEY else True)
        st = check_blob(spec)
        out.append({"key": key, "ref": name, "sha256": spec.sha256, "present": st["state"] != "missing",
                    "verified": st["state"] in ("verified", "trusted"),
                    "mismatch": st["detail"] if st["state"] in ("mismatch", "wrong_size") else None, "needed": needed})
        if key == MODEL_KEY and labels_run:
            warnings.append(f"the recipe reads its nuclei from labels files ({LABELS_KEY}), so the segmentation model "
                            f"{name} is not needed: its digest is not checked (D12)")
    return out, warnings


def mismatch_text(record: dict, pack) -> str:
    """What is wrong with a model whose file is there but is not the one expected: a pinned pack model (a ref
    with "@") whose blob is not the pinned file, or a project model whose file is not the one the manifest
    records."""
    if "@" in record["ref"]:
        return (f"is not the pinned blob ({record['mismatch']}): fetch it again with `germquant pack fetch-models "
                f"{pack.id}`")
    return f"is not the file the manifest records: {record['mismatch']}"


def model_refusal(records: list[dict], pack) -> str | None:
    """Why the run is refused for its models (a needed model missing, or not the one the manifest records)."""
    problems = []
    for r in records:
        if not r["needed"]:
            continue
        what = f"the model {r['ref']} ({r['key']})"
        if not r["present"]:
            fix = (f"germquant pack fetch-models {pack.id}" if "@" in r["ref"] else
                   "restore the file, or choose another model with `germquant project use-model`")
            problems.append(f"{what} is missing: {fix}")
        elif r["mismatch"]:
            problems.append(f"{what} {mismatch_text(r, pack)}")
    if not problems:
        return None
    return "; ".join(problems) + ". The run is refused before any pixel is read"


def _model_resolution_refusal(cfg) -> str | None:
    """Why the recipe's segmentation model cannot be resolved here (a registered or pinned model that is missing, a
    pin that does not match), or None: process_image resolves it before any stage of every image that reads no
    labels from files, so a run would fail on each such image."""
    from .. import pipeline
    from ..models_registry import ModelPinError

    try:
        pipeline.resolve_model_path(cfg)
    except (FileNotFoundError, ModelPinError) as e:
        return (f"the recipe's segmentation model ({MODEL_KEY}: {cfg.get(MODEL_KEY)!r}) cannot be resolved here: {e}. "
                f"The run is refused before any pixel is read")
    return None


def _pack_refusal(pack, recipe) -> str | None:
    """The pack check's blockers for this pack recipe (pack_status, shallow), its model blockers left out."""
    from ..packs import pack_status

    st = pack_status(pack, recipe=recipe.id)
    models = [f"model {name} is not usable" for name in recipe.need("models")]
    blockers = [b for b in st["blockers"] if not any(str(b["what"]).startswith(m) for m in models)]
    if not blockers:
        return None
    return (f"pack {pack.ref} recipe {recipe.id} is not usable here: "
            + "; ".join(f"{b['what']} (fix: {b['fix']})" for b in blockers))


# ---- the device -----------------------------------------------------------------------------------------------

def check_device(device) -> str | None:
    """``device`` lower case when it is auto, cuda, mps or cpu; None for None or ""; ProjectError otherwise."""
    if device is None or str(device).strip() == "":
        return None
    d = str(device).strip().lower()
    if d not in DEVICES:
        raise ProjectError(f"device {device!r} is not one of {', '.join(DEVICES)}")
    return d


@contextlib.contextmanager
def device_env(device: str | None):
    """GERMQUANT_DEVICE set to ``device`` for the block and restored afterwards (left alone for None)."""
    if device is None:
        yield
        return
    old = os.environ.get("GERMQUANT_DEVICE")
    os.environ["GERMQUANT_DEVICE"] = device
    try:
        yield
    finally:
        if old is None:
            os.environ.pop("GERMQUANT_DEVICE", None)
        else:
            os.environ["GERMQUANT_DEVICE"] = old


# ---- one image ------------------------------------------------------------------------------------------------

def _read_json(path) -> dict | None:
    import json

    try:
        with open(long_path(path), "rb") as fh:
            data = json.loads(fh.read().decode("utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def finished_fields(out_dir: Path, mid: str) -> dict:
    """What an image's results folder says: done_sha256 (LF-normalised, D48), config_hash, stage_hashes and
    failed_stages from its done marker, qc_pass from its image summary, segment_device_used from its stage
    record (None for what is not there)."""
    marker = _read_json(out_dir / f"{mid}__done.json") or {}
    stages_doc = _read_json(out_dir / f"{mid}__stages.json") or {}
    post = stages_doc.get("post_run") if isinstance(stages_doc.get("post_run"), dict) else {}
    return {"done_sha256": done_sha256(out_dir / f"{mid}__done.json") if marker else None,
            "config_hash": marker.get("config_hash"), "stage_hashes": marker.get("stage_hashes"),
            "failed_stages": marker.get("failed_stages"), "qc_pass": summary_qc_pass(out_dir, mid),
            "segment_device_used": post.get("segment_device_used")}


def summary_qc_pass(out_dir: Path, mid: str) -> bool | None:
    import csv

    try:
        with open(long_path(out_dir / f"{mid}__image_summary.csv"), newline="", encoding="utf-8") as fh:
            row = next(csv.DictReader(fh), None)
    except OSError:
        return None
    if not row or "qc_pass" not in row:
        return None
    return str(row["qc_pass"]).strip().lower() == "true"


def _series_fields(res: dict) -> tuple[str, dict]:
    rec = res.get("record") or {}
    entries = list(rec.get("timepoints") or [])
    failed_stages = sorted({s for e in entries for s in (e.get("failed_stages") or [])})
    done = [e for e in entries if e.get("status") == "done"]
    fields = {"done_sha256": None, "config_hash": rec.get("config_hash"), "stage_hashes": None,
              "failed_stages": failed_stages, "qc_pass": bool(done) and all(bool(e.get("qc_pass")) for e in done),
              "segment_device_used": None}
    if res.get("failed"):
        bad = ", ".join(f"t{int(f['t_index']):03d}" for f in res["failed"])
        fields["error"] = f"timepoint {bad} failed: " + "; ".join("; ".join(f.get("flags") or []) for f in res["failed"])
        return "image_failed", fields
    if entries and all(e.get("resumed") for e in entries):
        return "image_skipped", fields
    return "image_done", fields


def _run_one(entry: dict, ctx: dict) -> tuple[str, dict]:
    """``(event, fields)`` for one image (module docstring, step 3)."""
    from .. import api, timelapse

    root, cfg = ctx["root"], ctx["cfg"]
    iid = str(entry["id"])
    path = image_path(root, entry)
    out_dir = Path(root) / RESULTS_DIR / iid
    fields = {k: None for k in IMAGE_FIELDS}
    fields.update(image_id=iid, out_dir=str(out_dir))
    t0 = time.monotonic()

    def finish(kind: str, **more) -> tuple[str, dict]:
        fields.update(more)
        fields["status"] = {"image_done": "done", "image_failed": "failed", "image_skipped": "skipped"}[kind]
        fields["elapsed_s"] = round(time.monotonic() - t0, 3)
        return kind, fields

    mid = marker_id(path)
    sheet = ctx.get("sheet")
    # the image's row of the project's sample sheet, by its image id and the path the run reads (a row keyed by a
    # file path names that one file); None in a project without a sheet
    row = None if sheet is None else sheet.row(mid, path)
    if is_series(entry):
        try:
            why = series_refusal(cfg, path, row=row, post_stages=ctx.get("post_entries") or None)
        except Exception as e:  # noqa: BLE001 - io.planar auto reads the header: one image failing must not stop the run
            if not reachable(path):
                return finish("image_failed", error=unreachable_text(path))
            log.exception("FAILED %s", iid)
            return finish("image_failed", error=(f"the header of this time series cannot be read (io.planar auto "
                                                 f"decides from it): {type(e).__name__}: {e}")[:4000])
        if why:
            return finish("image_failed", error=series_refusal_text(entry, why, root))
        if not reachable(path):
            return finish("image_failed", error=unreachable_text(path))
        try:
            res = timelapse.process_series(path, cfg, out_dir, timepoints="all", xy_stride=ctx["xy_stride"],
                                           z_range=ctx["z_range"], resume=ctx["resume"])
        except Exception as e:  # noqa: BLE001 - one image failing must not stop the run
            log.exception("FAILED %s", iid)
            return finish("image_failed", error=f"{type(e).__name__}: {e}"[:4000])
        kind, more = _series_fields(res)
        return finish(kind, **more)
    if ctx["resume"]:
        try:
            done, why = ctx["hashes"].is_done(out_dir, mid, path, xy_stride=ctx["xy_stride"], z_range=ctx["z_range"],
                                              row=row)
        except Exception as e:  # noqa: BLE001 - the resume check of one image failing must not stop the run
            log.exception("FAILED %s", iid)
            return finish("image_failed", error=f"the resume check failed: {type(e).__name__}: {e}"[:4000])
        if done:
            if ctx.get("post") is not None:
                _post_only(out_dir, mid, path, ctx)        # its due post-stages only
            return finish("image_skipped", **finished_fields(out_dir, mid))
        if why != "no completion marker":
            log.info("%s: reprocessing: %s", iid, why)
    if not reachable(path):
        return finish("image_failed", error=unreachable_text(path))
    kw = {} if row is None else {"sample": row}          # a project with a sheet only (D40)
    if ctx.get("post") is not None:
        kw["post_stages"] = ctx["post"]                  # its post-stage list, only when there is one (D40)
    try:
        res = api.run_loaded(path, cfg, out_dir, xy_stride=ctx["xy_stride"], z_range=ctx["z_range"], **kw)
    except Exception as e:  # noqa: BLE001 - one image failing must not stop the run
        log.exception("FAILED %s", iid)
        return finish("image_failed", error=f"{type(e).__name__}: {e}"[:4000])
    more = finished_fields(out_dir, str(res.get("image_id") or mid))
    more["qc_pass"] = bool(res["qc_pass"]) if res.get("qc_pass") is not None else more["qc_pass"]
    return finish("image_done", **more)


# ---- the run --------------------------------------------------------------------------------------------------

def _disabled_stages(manifest: dict) -> list[str]:
    from ..stages import cli_switches

    names = {st.config_key: st.name for st in cli_switches()}
    over = manifest["recipe"].get("overrides") or {}
    return [names.get(k, k) for k in sorted(over) if over[k] is False]


def _start_fields(root, m: dict, cfg, entries: list[dict], *, via: str, shard: str | None, xy_stride: int) -> dict:
    from .. import __version__
    from ..device import requested_device, select_device
    from ..provenance import git_sha

    prefer = cfg.get("segmentation.nuclei.device")
    p = m["pack"]
    return {
        "via": via, "host": socket.gethostname(), "os": platform.platform(), "interpreter": sys.executable,
        "interpreter_source": "engine", "engine": {"germquant": __version__, "git_sha": git_sha()},
        "python": platform.python_version(), "constraints": None,
        "device_requested": requested_device(prefer), "device_selected": select_device(prefer, quiet=True),
        "config_hash": cfg.hash, "config_hash_source": "the engine: the project's own recipe, as this run loaded it",
        "recipe": {"kind": "project", "pack": f"{p['id']}@{p['version']}", "recipe": p["recipe"],
                   "edited": bool(m["recipe"].get("file")), "file": m["recipe"].get("file")},
        "xy_stride": int(xy_stride), "disabled_stages": _disabled_stages(m), "argv": RUN_ARGV.get(),
        "images": [{"path": str(image_path(root, e)), "image_id": str(e["id"]),
                    "out_dir": str(Path(root) / RESULTS_DIR / str(e["id"]))} for e in entries],
        "plugin": None, "log_location": "the project's runs/ folder",
        "project": {"root": str(root), "name": m.get("name"), "id": m.get("id")}, "shard": shard,
    }


def run_project(project, images=None, resume=True, device=None, shard=None, progress=None, *, via="cli",
                run_id=None, should_stop=None) -> dict:
    """Run the project's images (module docstring): ``{"run_id", "log", "counts": {"done", "failed", "skipped"},
    "exit_code", "interrupted", "warnings", "message"}``. exit_code is 0 when every image was done or skipped and
    1 when one failed; a refusal before any image raises ProjectError (exit 2 of `germquant project run`)."""
    import germquant.pipeline  # noqa: F401  (D28: torch claims cuBLAS before any reader or backend loads CuPy)

    from ..packs import use_project
    from . import Project, _root, project_config
    from .lock import ProjectLockedError, describe, is_stale, read_lock, unlock_hint
    from .recipe import resolve_pack

    if via not in VIA:
        raise ProjectError(f"via {via!r} is not one of {', '.join(VIA)}")
    run_id = new_run_id() if run_id is None else check_run_id(run_id)
    k_n = parse_shard(shard)
    asked_device = check_device(device)
    root = _root(project)
    m = read_manifest(root)                                  # fresh: what is on disk now
    warnings: list[str] = []
    holder = read_lock(root)
    if holder is not None:
        if not is_stale(root, holder):
            raise ProjectLockedError(f"the project {root} is locked by {describe(holder)}: a project edit is being "
                                     f"written; run again when that command has finished", holder)
        warnings.append(f"the project is locked by {describe(holder)}, and the lock is older than 12 h: "
                        f"{unlock_hint(root, holder)}")
    entries = select_images(m, images, k_n)
    run_dev = asked_device if asked_device is not None else check_device(m["run"].get("device"))
    xy_stride = int(m["run"].get("xy_stride") or 1)
    zr = m["run"].get("z_range")
    z_range = (int(zr[0]), int(zr[1])) if zr else None
    with use_project(root):
        cfg = project_config(Project(root=root, manifest=m))
        pack, recipe = resolve_pack(root, m)
        why = _pack_refusal(pack, recipe)
        if why:
            raise ProjectError(why)
        records, model_warnings = model_records(root, m, cfg, pack, recipe)
        why = model_refusal(records, pack)
        if why:
            raise ProjectError(why)
        if cfg.get(LABELS_KEY) is None and any(not is_series(e) for e in entries):
            why = _model_resolution_refusal(cfg)          # process_image resolves the model on every such image
            if why:
                raise ProjectError(why)
        sheet = read_project_sheet(root, m)               # once, before anything is written (D21; 9.4 #13)
        post, post_entries = project_post_stages(root, pack, recipe)    # read and checked once, before the log
        warnings += model_warnings
        if sheet is not None:
            warnings += sheet.warnings
        for w in warnings:
            log.warning("%s", w)
        rlog = RunLog(root, run_id)
        events: list[dict] = []

        def emit(rec: dict) -> None:
            events.append(rec)
            if progress is not None:
                try:
                    progress(rec)
                except Exception:  # noqa: BLE001 - a caller's progress display must not stop the run
                    log.exception("the progress callback failed")

        counts = {"done": 0, "failed": 0, "skipped": 0}
        interrupted = False
        failures: list[str] = []
        ctx = {"root": root, "cfg": cfg, "xy_stride": xy_stride, "z_range": z_range, "resume": bool(resume),
               "hashes": ResumeHashes(cfg), "sheet": sheet, "post": post, "post_entries": post_entries}
        with device_env(run_dev):
            emit(rlog.start(**_start_fields(root, m, cfg, entries, via=via,
                                            shard=None if k_n is None else f"{k_n[0]}/{k_n[1]}", xy_stride=xy_stride)))
            current = None
            try:
                for entry in entries:
                    if should_stop is not None and should_stop():
                        interrupted = True
                        break
                    current = entry
                    kind, fields = _run_one(entry, ctx)
                    current = None
                    counts[fields["status"]] += 1
                    if kind == "image_failed":
                        failures.append(f"{fields['image_id']} ({fields['error']})")
                    emit(rlog.event(kind, **fields))
            except BaseException as e:   # Ctrl-C, or an error outside the images: the run says it stopped part way
                with contextlib.suppress(Exception):
                    if current is not None:
                        emit(rlog.event("image_failed", **{**{k: None for k in IMAGE_FIELDS},
                                                            "image_id": str(current["id"]),
                                                            "out_dir": str(Path(root) / RESULTS_DIR /
                                                                           str(current["id"])),
                                                            "status": "failed",
                                                            "error": f"{type(e).__name__}: the run was stopped"}))
                        counts["failed"] += 1
                    emit(rlog.event("run_end", counts=dict(counts), interrupted=True))
                if isinstance(e, Exception):
                    raise RunError(f"the run {run_id} stopped part way ({type(e).__name__}: {e}); its log: "
                                   f"{rlog.path}", run_id=run_id, log=str(rlog.path)) from e
                raise
            emit(rlog.event("run_end", counts=dict(counts), interrupted=interrupted))
    exit_code = 1 if counts["failed"] else 0
    message = ""
    if failures:
        message = f"{len(failures)} of {len(entries)} images failed: " + "; ".join(failures)
    return {"run_id": run_id, "log": str(rlog.path), "counts": counts, "exit_code": exit_code,
            "interrupted": interrupted, "warnings": warnings, "message": message[:4000]}
