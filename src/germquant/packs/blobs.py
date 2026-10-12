"""The model blob store and `germquant pack fetch-models` (docs/DESIGN_system_packs.md, section 3.5).

Weights live once per digest at ``models_dir()/blobs/<sha256>/<file_name>``, with the verification record
``<sha256>.json`` beside them (source, size, mtime, when it was fetched and verified, the packs that use it).
models_dir() honours ``$GERMQUANT_MODELS``. A release model (kind release, such as germline_nuclei_combined)
keeps its place at ``models_dir()/<name>``, where `germquant get-model` puts it, and gets the same record.

fetch-models tries, in order: the file already in the store; Cellpose's own cache (CELLPOSE_LOCAL_MODELS_PATH
or ~/.cellpose/models), after a size and sha256 check, hard-linked when possible and copied otherwise; a file
given with --from-file, after the same check; a download of the pinned URL (datafetch: resumable, size and
sha256 checked). A mismatch is discarded and never registered, and a rerun does nothing.

pack status trusts a file whose size and mtime match its record and re-hashes it only when either changed,
so a card refresh in a new process does not re-hash 1.2 GB per model. The run path still hashes the file
once per process (provenance.file_sha256) and applies the pin check (models_registry.check_model_pin).
"""
from __future__ import annotations

import json
import logging
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

from .format import ModelSpec, Pack

log = logging.getLogger(__name__)


def blobs_dir() -> Path:
    from ..models_registry import models_dir

    return models_dir() / "blobs"


def model_location(spec: ModelSpec) -> Path:
    """Where a model's weights live: ``models_dir()/<name>`` for a release model, else its blob path."""
    from ..models_registry import models_dir

    if spec.kind == "release":
        return models_dir() / spec.name
    return blobs_dir() / spec.sha256 / spec.file_name


def record_path(sha256: str) -> Path:
    return blobs_dir() / sha256 / f"{sha256}.json"


def read_record(sha256: str) -> dict | None:
    try:
        with open(record_path(sha256), encoding="utf-8") as fh:
            rec = json.load(fh)
        return rec if isinstance(rec, dict) else None
    except (OSError, ValueError):
        return None


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_record(spec: ModelSpec, path: Path, *, source: str | None = None, packs: list[str] | None = None,
                 verified: bool = True) -> dict:
    """Write (or refresh) the verification record of ``spec`` for the file at ``path``."""
    st = Path(path).stat()
    old = read_record(spec.sha256) or {}
    rec = {
        "format": 1, "name": spec.name, "sha256": spec.sha256, "size": int(st.st_size),
        "mtime_ns": int(st.st_mtime_ns), "path": str(Path(path).resolve()), "file_name": spec.file_name,
        "source": source or old.get("source", "unknown"), "fetched_at": old.get("fetched_at") or _now(),
        "verified_at": _now() if verified else old.get("verified_at"),
        "packs": sorted(set(old.get("packs", [])) | set(packs or [])),
    }
    p = record_path(spec.sha256)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".json.part")
    try:
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(rec, fh, indent=1)
        os.replace(tmp, p)
    finally:
        try:
            tmp.unlink(missing_ok=True)                # left only when the write or the rename failed
        except OSError:
            pass
    return rec


def _save_record(spec: ModelSpec, path: Path, **kw) -> None:
    """write_record for a file already verified: a model store that cannot be written (a read-only
    $GERMQUANT_MODELS on a cluster or in a container) only costs a re-hash next time, never the check."""
    try:
        write_record(spec, path, **kw)
    except OSError as e:
        log.warning("could not save the verification record of %s (%s); it is re-hashed next time", path, e)


def check_blob(spec: ModelSpec, *, path: Path | None = None, rehash: bool = True) -> dict:
    """The state of a model's weights: ``missing``, ``wrong_size``, ``trusted`` (size and mtime match the
    verification record, not re-hashed), ``verified`` (re-hashed and equal; the record is refreshed when the
    store can be written), ``mismatch`` (re-hashed and different) or, with ``rehash=False`` and no matching
    record, ``unverified`` (present with the right size). Returns ``{"state", "path", "detail"}``."""
    from ..datafetch import sha256_file

    p = Path(path) if path is not None else model_location(spec)
    if not p.is_file():
        return {"state": "missing", "path": str(p), "detail": "not in the model store"}
    st = p.stat()
    if spec.size and st.st_size != spec.size:
        return {"state": "wrong_size", "path": str(p), "detail": f"{st.st_size} bytes, expected {spec.size}"}
    rec = read_record(spec.sha256)
    if (rec and rec.get("sha256") == spec.sha256 and rec.get("size") == st.st_size
            and rec.get("mtime_ns") == st.st_mtime_ns and rec.get("path") == str(p.resolve())):
        return {"state": "trusted", "path": str(p), "detail": f"verified {rec.get('verified_at')}"}
    if not rehash:
        return {"state": "unverified", "path": str(p), "detail": "right size, not verified yet"}
    got = sha256_file(p)
    if got != spec.sha256:
        return {"state": "mismatch", "path": str(p), "detail": f"sha256 {got[:12]}..., expected {spec.sha256[:12]}..."}
    _save_record(spec, p)
    return {"state": "verified", "path": str(p), "detail": "sha256 checked now"}


def find_model_spec(name: str) -> tuple[Pack, ModelSpec] | None:
    """The first installed or built-in pack (in search order) that declares a model ``name``."""
    from .discovery import discover

    for pack in discover().packs:
        spec = pack.model(name)
        if spec is not None:
            return pack, spec
    return None


def resolve_pack_model(name: str) -> Path:
    """The weights a pack model name ``<name>@<8 hex>`` names: the blob of the declaring pack's model spec,
    present with the right size (the run path hashes it and applies the pin check). FileNotFoundError naming
    `germquant pack fetch-models <pack>` when the blob is missing, and naming the packs available when no
    installed pack declares the name."""
    from ..models_registry import ModelPinError

    found = find_model_spec(name)
    if found is None:
        from .discovery import available_refs

        raise FileNotFoundError(f"model {name!r} is a pack model name, but no installed or built-in pack declares it; "
                                f"packs available: {', '.join(available_refs()) or 'none'} (install the pack that "
                                f"ships it with `germquant pack install`)")
    pack, spec = found
    st = check_blob(spec, rehash=False)
    if st["state"] in ("missing", "wrong_size"):
        raise FileNotFoundError(f"model {name!r} of pack {pack.ref} is not in the model store ({st['path']}: "
                                f"{st['detail']}); fetch it with `germquant pack fetch-models {pack.id}`")
    if st["state"] == "mismatch":
        raise ModelPinError(f"model {name!r}: {st['path']} {st['detail']}")
    return Path(st["path"])


def cellpose_cache_dirs() -> list[Path]:
    """Cellpose's model cache: CELLPOSE_LOCAL_MODELS_PATH when set, else ~/.cellpose/models."""
    env = os.environ.get("CELLPOSE_LOCAL_MODELS_PATH")
    return [Path(env)] if env else [Path.home() / ".cellpose" / "models"]


def _matches(spec: ModelSpec, path: Path) -> tuple[bool, str]:
    from ..datafetch import sha256_file

    if not path.is_file():
        return False, f"{path}: not a file"
    size = path.stat().st_size
    if spec.size and size != spec.size:
        return False, f"{path}: {size} bytes, expected {spec.size}"
    got = sha256_file(path)
    if got != spec.sha256:
        return False, f"{path}: sha256 {got[:12]}..., expected {spec.sha256[:12]}..."
    return True, "ok"


def _place(src: Path, dest: Path, *, link: bool) -> str:
    """Put a verified file at ``dest``: a hard link when ``link`` and the filesystem allows, else a copy."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".tmp")
    if tmp.exists():
        tmp.unlink()
    how = "copied"
    if link:
        try:
            os.link(src, tmp)
            how = "hard-linked"
        except OSError:
            how = "copied"
    if how == "copied":
        shutil.copyfile(src, tmp)
    os.replace(tmp, dest)
    return how


def needed_specs(pack: Pack, recipe: str | None = None) -> list[ModelSpec]:
    if recipe is None:
        return list(pack.models)
    names = pack.recipe(recipe).need("models")
    return [m for m in pack.models if m.name in names]


def fetch_models(pack: Pack, recipe: str | None = None, from_files: dict[str, str | Path] | None = None, *,
                 download: bool = True, say=None) -> list[dict]:
    """Put every model the recipe (or, without one, the pack) needs into the store. Returns one record per
    model: ``{"name", "status": present|adopted|from_file|downloaded|missing, "path", "detail", "tried"}``.
    PackFormatError for an unusable pack (its manifest has problems): store paths are built from its model
    names, digests and file names, which only a valid manifest vouches for."""
    from .format import PackFormatError

    if pack.problems:
        raise PackFormatError(f"pack {pack.ref} is unusable, so no model is fetched: " + "; ".join(pack.problems))
    say = say or (lambda *a: None)
    out = []
    files = {str(k): Path(v) for k, v in (from_files or {}).items()}
    unknown = sorted(set(files) - {m.name for m in pack.models})
    if unknown:
        raise ValueError(f"--from-file names models pack {pack.ref} does not declare: {unknown}")
    for spec in needed_specs(pack, recipe):
        out.append(_fetch_one(pack, spec, files.get(spec.name), download=download, say=say))
    return out


def _fetch_one(pack: Pack, spec: ModelSpec, from_file: Path | None, *, download: bool, say) -> dict:
    loc = model_location(spec)
    tried: list[str] = []
    st = check_blob(spec)
    if st["state"] in ("trusted", "verified"):
        _save_record(spec, loc, packs=[pack.ref], verified=st["state"] == "verified")
        say(f"{spec.name}: already in the store ({loc})")
        return {"name": spec.name, "status": "present", "path": str(loc), "detail": st["detail"], "tried": tried}
    if st["state"] in ("wrong_size", "mismatch"):
        tried.append(f"store: {st['detail']}")
        if spec.kind != "release":            # our own blob, wrong: discarded, never kept
            Path(st["path"]).unlink(missing_ok=True)
    if spec.kind != "release":
        for src in spec.sources:
            if src.get("type") != "cellpose_cache":
                continue
            for d in cellpose_cache_dirs():
                cand = d / str(src.get("name") or spec.file_name)
                if not cand.exists():
                    tried.append(f"cellpose cache: {cand} absent")
                    continue
                ok, why = _matches(spec, cand)
                if not ok:
                    tried.append(f"cellpose cache: {why}")
                    continue
                how = _place(cand, loc, link=True)
                write_record(spec, loc, source=f"cellpose_cache {cand}", packs=[pack.ref])
                say(f"{spec.name}: adopted from the Cellpose cache {cand} ({how})")
                return {"name": spec.name, "status": "adopted", "path": str(loc), "detail": how, "tried": tried}
    if from_file is not None:
        ok, why = _matches(spec, from_file)
        if ok:
            how = _place(from_file, loc, link=False)
            write_record(spec, loc, source=f"file {from_file.resolve()}", packs=[pack.ref])
            say(f"{spec.name}: copied from {from_file}")
            return {"name": spec.name, "status": "from_file", "path": str(loc), "detail": how, "tried": tried}
        tried.append(f"--from-file: {why}")
    if download:
        for src in spec.sources:
            kind = src.get("type")
            if kind == "url" and src.get("url"):
                from ..datafetch import DownloadError
                from ..datafetch import download as fetch

                url = str(src["url"])
                say(f"{spec.name}: downloading {spec.size} bytes from {url}")
                try:
                    fetch(url, loc, size=spec.size or None, sha256=spec.sha256)
                except DownloadError as e:
                    tried.append(f"url: {e}")
                    continue
                except Exception as e:  # noqa: BLE001 - a network failure: the .part file stays for a resume
                    tried.append(f"url: {type(e).__name__}: {e}")
                    continue
                write_record(spec, loc, source=f"url {url}", packs=[pack.ref])
                return {"name": spec.name, "status": "downloaded", "path": str(loc), "detail": url, "tried": tried}
            if kind == "release" and spec.kind == "release":
                from ..models_registry import get_model

                try:
                    path = get_model(spec.name, quiet=True)
                except Exception as e:  # noqa: BLE001 - reported with the others
                    tried.append(f"release: {type(e).__name__}: {e}")
                    continue
                write_record(spec, path, source="release asset", packs=[pack.ref])
                return {"name": spec.name, "status": "downloaded", "path": str(path), "detail": "release asset",
                        "tried": tried}
    say(f"{spec.name}: not fetched ({'; '.join(tried) or 'no source'})")
    return {"name": spec.name, "status": "missing", "path": str(loc), "detail": "no source matched", "tried": tried}
