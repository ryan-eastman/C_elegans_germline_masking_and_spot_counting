"""Where the trained models live and how to fetch them: `germquant get-model` and `germquant register-model`.

The trained Cellpose nucleus model is 1.2 GB and is not in the git repository; it is published as a
GitHub release asset of this repository. `get_model` downloads it to ``models/models/<name>`` under
the repository (or under ``$GERMQUANT_MODELS`` when that is set, the HPC / container binding),
verifies the sha256 and refuses a corrupt or truncated download, so the path config.yaml already
points at (``segmentation.nuclei.cellpose_model: models/models/germline_nuclei_combined``) just
works on a new computer. A new model version is a new release tag and a new entry here.

Local models (docs/TIMELAPSE.md, "Segmenting without DAPI: the training loop"): a model trained on this
computer with ``germquant train-nuclei`` is added to the registry with ``germquant register-model``, which
records its name, path, sha256, size and model card in ``local_models.json`` next to the release assets
(``models_dir() / "local_models.json"``). ``get_model(name)`` resolves a local name by verifying that file
instead of downloading, and ``segmentation.nuclei.cellpose_model: <name>`` in a profile runs it
(pipeline.resolve_model_path asks ``resolve_registered``).
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath, PureWindowsPath

log = logging.getLogger(__name__)

# A pack model name: <name>@<first 8 hex of its sha256> (docs/DESIGN_system_packs.md, section 3.5). The form is
# reserved for packs: register_local_model refuses it, and a value of this form resolves only to weights whose
# sha256 starts with the suffix and equals the pack spec's digest (check_model_pin).
PACK_MODEL_RE = re.compile(r"^([A-Za-z0-9][A-Za-z0-9_.-]{0,63})@([0-9a-f]{8})$")


class ModelPinError(RuntimeError):
    """A ``<name>@<8 hex>`` model value resolved to weights with another digest than its pin."""


def is_pack_model_name(name) -> bool:
    return isinstance(name, str) and PACK_MODEL_RE.fullmatch(name.strip()) is not None

REPO = "ryan-eastman/C_elegans_germline_masking_and_spot_counting"

MODELS = {
    "germline_nuclei_combined": {
        "tag": "model-germline_nuclei_combined-v1",
        "asset": "germline_nuclei_combined",
        "sha256": "e4e32320cfc321a145ccbc0c52138eded50e37b87130f34bae2296c88b148f57",
        "size": 1218646375,
        "what": "Cellpose-SAM fine-tuned on real hand-labelled germline nuclei plus synthetic nuclei "
                "(June 2026): the recommended nucleus model, used by every validated run.",
    },
}

LOCAL_REGISTRY_FILE = "local_models.json"
# names Cellpose resolves itself; a local model may not take one of them (it would shadow the stock model)
CELLPOSE_BUILTIN_NAMES = ("cpsam", "cpsam_v2", "cpdino", "cpdino-vitb", "sam")


def asset_url(name: str) -> str:
    m = MODELS[name]
    return f"https://github.com/{REPO}/releases/download/{m['tag']}/{m['asset']}"


def models_dir(base_dir: str | Path | None = None) -> Path:
    """``$GERMQUANT_MODELS/models`` when the variable is set, else ``<repo>/models/models``."""
    root = os.environ.get("GERMQUANT_MODELS")
    if root:
        return Path(root) / "models"
    if base_dir is None:
        from .provenance import _repo_root

        base_dir = _repo_root() or Path.cwd()
    return Path(base_dir) / "models" / "models"


def sha256_of(path: Path, progress: bool = False) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


# ---- local models (germquant register-model) ---------------------------------------------------------

def local_registry_path(base_dir: str | Path | None = None) -> Path:
    """``models_dir(base_dir) / local_models.json``: the file that lists the models registered on this
    computer (under ``$GERMQUANT_MODELS`` when that is set, as the release assets are)."""
    return models_dir(base_dir) / LOCAL_REGISTRY_FILE


def local_models(base_dir: str | Path | None = None) -> dict[str, dict]:
    """``{name: entry}`` of the locally registered models (``path``, ``sha256``, ``size``, ``card``,
    ``registered``, ``what``); empty when nothing has been registered or the file is unreadable."""
    p = local_registry_path(base_dir)
    try:
        with open(p, encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return {}
    models = data.get("models", {}) if isinstance(data, dict) else {}
    return {str(k): dict(v) for k, v in models.items() if isinstance(v, dict)}


def known_models(base_dir: str | Path | None = None) -> list[str]:
    """Every name ``get_model`` accepts: the release assets, then the local models, sorted."""
    return sorted(set(MODELS) | set(local_models(base_dir)))


def register_local_model(model: str | Path, card: str | Path | None = None, name: str | None = None,
                         base_dir: str | Path | None = None) -> dict:
    """Add a trained model file to the local registry (``germquant register-model``) and return its entry.

    ``model`` is the model file ``germquant train-nuclei`` wrote; ``card`` its model card JSON (optional but
    recommended: it is what says which slices the model was trained on and how it scored on the hold-out
    set). ``name`` defaults to the card's ``name``, else the file's name. The sha256 of the file is
    computed here and, when the card carries ``model_sha256``, must agree with it, so a card cannot be
    attached to the wrong weights. A name that is a release asset or a Cellpose built-in is refused;
    registering an existing local name again replaces its entry (a retrained model under the same name).
    The entry's ``path`` and ``card`` are absolute, so the registry file can be read from any working
    directory."""
    model_path = Path(model).resolve()
    if not model_path.is_file():
        raise FileNotFoundError(f"model file not found: {model_path}")
    card_path = Path(card).resolve() if card is not None else None
    card_data: dict = {}
    if card_path is not None:
        if not card_path.is_file():
            raise FileNotFoundError(f"model card not found: {card_path}")
        with open(card_path, encoding="utf-8") as fh:
            card_data = json.load(fh)
        if not isinstance(card_data, dict):
            raise ValueError(f"model card {card_path} is not a JSON object")
    name = str(name or card_data.get("name") or model_path.name).strip()
    if not name or any(ch in name for ch in "/\\") or name in (".", ".."):
        raise ValueError(f"model name {name!r} must be a bare name (no path separators)")
    if "@" in name:
        raise ValueError(f"{name!r}: names of the form <name>@<8 hex> are reserved for pack models "
                         f"(docs/PACKS.md); pick a name without '@'")
    if name in MODELS:
        raise ValueError(f"{name!r} is a release model; pick another name for a local model")
    if name.lower() in CELLPOSE_BUILTIN_NAMES:
        raise ValueError(f"{name!r} is a Cellpose built-in model name; pick another name")
    sha = sha256_of(model_path)
    card_sha = card_data.get("model_sha256")
    if card_sha and str(card_sha).lower() != sha:
        raise ValueError(f"the card's model_sha256 ({str(card_sha)[:12]}...) does not match the file "
                         f"({sha[:12]}...): this card belongs to other weights")
    entry = {
        "path": str(model_path), "sha256": sha, "size": model_path.stat().st_size,
        "card": str(card_path) if card_path is not None else None,
        "registered": datetime.now(timezone.utc).isoformat(),
        "what": str(card_data.get("what") or card_data.get("description") or ""),
        "base_model": card_data.get("base_model"),
        "trained": card_data.get("date"),
    }
    reg_path = local_registry_path(base_dir)
    reg_path.parent.mkdir(parents=True, exist_ok=True)
    models = local_models(base_dir)
    models[name] = entry
    tmp = reg_path.with_suffix(".json.part")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump({"format": 1, "models": models}, fh, indent=1, default=str)
    os.replace(tmp, reg_path)
    log.info("registered local model %r -> %s (%s)", name, model_path, reg_path)
    return {"name": name, **entry}


def resolve_registered(name: str, base_dir: str | Path | None = None) -> Path | None:
    """The file a registered model NAME points at, when that file is on disk: a local entry's ``path``,
    or ``models_dir(base_dir) / name`` for a release asset that ``get-model`` has downloaded. None for an
    unregistered name or a registered one whose file is missing (the caller then treats the value as a
    Cellpose built-in name or a path, as before). A local file whose size differs from the registered one
    is logged as a warning and still returned: the run's manifest records the digest that actually ran."""
    key = str(name).strip()
    if not key:
        return None
    entry = local_models(base_dir).get(key)
    if entry is not None:
        p = Path(str(entry.get("path", "")))
        if p.is_file():
            if int(entry.get("size", -1)) not in (-1, p.stat().st_size):
                log.warning("registered model %r at %s has size %d, registered as %s: re-register it",
                            key, p, p.stat().st_size, entry.get("size"))
            return p
        log.warning("registered model %r points at a missing file %s", key, p)
        return None
    if key in MODELS:
        p = models_dir(base_dir) / key
        if p.is_file():
            return p
    return None


# ---- pack models, the pin check and model cards (docs/DESIGN_system_packs.md, sections 3.5 and 3.6) --------------

def check_model_pin(name: str, sha256: str, *, spec_sha256: str | None = None) -> None:
    """Refuse weights that a ``<name>@<8 hex>`` value does not pin: their sha256 must start with the suffix
    and equal the digest of the installed pack's model spec (``spec_sha256``, else looked up in the installed
    packs; when no installed pack declares the name, only the suffix is checked). Any other name passes. A
    local registry entry is checked like any other file, so a wrong local entry is caught here (pack status
    never looks at the local registry). Raises ModelPinError."""
    m = PACK_MODEL_RE.fullmatch(str(name).strip())
    if not m:
        return
    got = str(sha256 or "").lower()
    if not got.startswith(m.group(2)):
        raise ModelPinError(f"model {name!r} resolved to weights with sha256 {got[:12] or 'absent'}..., which does "
                            f"not start with its pin {m.group(2)}: the file is not the pinned model (check the "
                            f"local model registry and $GERMQUANT_MODELS; `germquant pack fetch-models` fetches "
                            f"the pinned weights)")
    if spec_sha256 is None:
        from .packs.blobs import find_model_spec       # lazy: packs imports this module

        found = find_model_spec(name)
        spec_sha256 = found[1].sha256 if found else None
    if spec_sha256 and got != str(spec_sha256).lower():
        raise ModelPinError(f"model {name!r} resolved to weights with sha256 {got[:12]}..., but the pack spec pins "
                            f"{str(spec_sha256)[:12]}...: the file is not the pinned model")


def resolve_model_name(name: str, expect_sha256: str | None = None, base_dir: str | Path | None = None) -> Path:
    """The file a model NAME resolves to: the public helper a third-party backend calls for its own model key
    (docs/BACKENDS.md). A registered local or release model with its file on disk (resolve_registered), else a
    pack model ``<name>@<8 hex>``, which resolves to the verified blob of an installed pack's model spec. A pack
    name is always pin-checked (check_model_pin), and with ``expect_sha256`` the file's sha256 must equal it.
    FileNotFoundError when nothing resolves (a declared pack model whose blob is missing names
    ``germquant pack fetch-models <pack>``); ModelPinError or ValueError on a digest mismatch."""
    key = str(name).strip()
    path = resolve_registered(key, base_dir)
    if path is None and is_pack_model_name(key):
        from .packs.blobs import resolve_pack_model

        path = resolve_pack_model(key)
    if path is None:
        raise FileNotFoundError(f"model {key!r} is not a registered model with its file on disk "
                                f"({local_registry_path(base_dir)}) and not a pack model <name>@<8 hex>")
    sha = None
    if is_pack_model_name(key):
        from .provenance import file_sha256

        sha = file_sha256(path)
        check_model_pin(key, sha)
    if expect_sha256:
        if sha is None:
            from .provenance import file_sha256

            sha = file_sha256(path)
        if sha != str(expect_sha256).lower():
            raise ValueError(f"model {key!r} at {path} has sha256 {sha[:12]}..., expected {str(expect_sha256)[:12]}...")
    return Path(path)


CARD_NON_PORTABLE = ("model_file", "base_model.path")


def _is_machine_path(value) -> bool:
    """An absolute path of any platform (C:/..., \\\\server\\share, /home/...), which a published card must not hold."""
    if not isinstance(value, str) or not value.strip():
        return False
    v = value.strip()
    return PureWindowsPath(v).is_absolute() or v.startswith(("\\\\", "//")) or (
        PurePosixPath(v).is_absolute() and "://" not in v)


def validate_card(card, spec: dict | None = None, *, pack_written: bool = False) -> dict:
    """Check a model card against its spec, using the keys ``germquant train-nuclei`` actually writes
    (segment/training.py: name, model_file, model_sha256, model_size, base_model, ...), so every card written
    today passes. ``size`` is accepted for ``model_size`` only in a pack-written card (``pack_written``).
    ``spec`` is a pack model spec (name, sha256, size): the card must carry the same digest and size, and a
    name that is the spec's name or its part before ``@``. ``model_file`` and ``base_model.path`` hold
    absolute machine paths in a train_nuclei card: they are reported in ``non_portable`` (pack derive rewrites
    them), not as problems. Returns ``{"ok", "problems", "non_portable"}``."""
    if not isinstance(card, dict):
        return {"ok": False, "problems": ["the card is not a JSON object"], "non_portable": []}
    problems: list[str] = []
    name = str(card.get("name") or "").strip()
    if not name:
        problems.append("no name")
    sha = str(card.get("model_sha256") or "").lower()
    if not re.fullmatch(r"[0-9a-f]{64}", sha):
        problems.append("model_sha256 is not a 64-hex sha256")
    size = card.get("model_size")
    if size is None and pack_written:
        size = card.get("size")
    if size is None:
        problems.append("no model_size" + (" (or size)" if pack_written else ""))
    elif isinstance(size, bool) or not isinstance(size, int) or size <= 0:
        problems.append(f"model_size {size!r} is not a positive integer")
    if spec:
        spec_name = str(spec.get("name", ""))
        if name and name not in (spec_name, spec_name.split("@")[0]):
            problems.append(f"the card names {name!r}, the spec {spec_name!r}")
        if spec.get("sha256") and sha != str(spec["sha256"]).lower():
            problems.append(f"the card's model_sha256 {sha[:12]}... differs from the spec's {str(spec['sha256'])[:12]}...")
        if spec.get("size") is not None and size is not None and size != spec.get("size"):
            problems.append(f"the card's model_size {size} differs from the spec's {spec.get('size')}")
    non_portable = []
    if _is_machine_path(card.get("model_file")):
        non_portable.append("model_file")
    base = card.get("base_model")
    if isinstance(base, dict) and _is_machine_path(base.get("path")):
        non_portable.append("base_model.path")
    return {"ok": not problems, "problems": problems, "non_portable": non_portable}


def get_model(name: str = "germline_nuclei_combined", base_dir: str | Path | None = None, *, force: bool = False,
              quiet: bool = False) -> Path:
    """Download the model if it is not already present and correct. Returns the model path. A locally
    registered name (``germquant register-model``) is verified in place instead: its file must exist with
    the registered sha256, else a RuntimeError says to retrain or re-register (nothing is downloaded)."""
    say = (lambda *a: None) if quiet else (lambda *a: print(*a, flush=True))
    local = local_models(base_dir)
    if name in local and name not in MODELS:
        entry = local[name]
        p = Path(str(entry.get("path", "")))
        if not p.is_file():
            raise RuntimeError(f"{name}: registered at {p} but the file is missing; retrain it or register the "
                               f"moved file again with `germquant register-model`")
        if sha256_of(p) != str(entry.get("sha256", "")).lower():
            raise RuntimeError(f"{name}: the file at {p} does not match its registered sha256; the weights "
                               f"changed after registration, register the file again")
        say(f"{name}: local model verified at {p}")
        return p
    if name not in MODELS:
        raise KeyError(f"unknown model {name!r}; known: {known_models(base_dir)}")
    m = MODELS[name]
    dest = models_dir(base_dir) / name
    say = (lambda *a: None) if quiet else (lambda *a: print(*a, flush=True))
    if dest.is_file() and not force:
        if dest.stat().st_size == m["size"] and sha256_of(dest) == m["sha256"]:
            say(f"{name}: already present and verified at {dest}")
            return dest
        say(f"{name}: present but wrong size or checksum; downloading again")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".part")
    url = asset_url(name)
    say(f"{name}: downloading {m['size'] // 1_000_000} MB from {url}")
    _download(url, tmp, m["size"], say)
    got = sha256_of(tmp)
    if got != m["sha256"]:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"{name}: checksum mismatch after download (got {got[:12]}..., expected "
                           f"{m['sha256'][:12]}...); the file was discarded, try again")
    os.replace(tmp, dest)
    say(f"{name}: verified and installed at {dest}")
    return dest


def _download(url: str, tmp: Path, size: int, say) -> None:
    import urllib.request

    req = urllib.request.Request(url, headers={"User-Agent": "germquant get-model"})
    t0 = time.time()
    done = 0
    last = 0
    with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as out:
        while True:
            chunk = r.read(1 << 22)
            if not chunk:
                break
            out.write(chunk)
            done += len(chunk)
            if done - last >= 100_000_000:
                last = done
                rate = done / max(time.time() - t0, 1e-6) / 1e6
                say(f"  {done // 1_000_000} MB of {size // 1_000_000} MB ({rate:.0f} MB/s)")
    if size and done != size:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"download ended early ({done} of {size} bytes); check the connection and try again")


def main(argv=None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description="download and verify a trained model (a release asset), or verify a "
                                             "locally registered one in place")
    ap.add_argument("name", nargs="?", default="germline_nuclei_combined",
                    help=f"a release model {sorted(MODELS)} or a name added with `germquant register-model`")
    ap.add_argument("--force", action="store_true", help="download even if a verified copy exists")
    a = ap.parse_args(argv)
    try:
        get_model(a.name, force=a.force)
    except Exception as e:  # noqa: BLE001 - the message is the point for a lab user
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    return 0


def register_main(argv=None) -> int:
    """``germquant register-model --model FILE --card CARD.json [--name NAME]``."""
    import argparse

    ap = argparse.ArgumentParser(description="add a locally trained nucleus model to the model registry")
    ap.add_argument("--model", required=True, help="the model file `germquant train-nuclei` wrote")
    ap.add_argument("--card", help="its model card JSON (<name>.card.json next to the model file)")
    ap.add_argument("--name", help="the registry name (default: the card's name, else the file name); the value "
                                   "segmentation.nuclei.cellpose_model takes in a profile")
    a = ap.parse_args(argv)
    try:
        entry = register_local_model(a.model, a.card, a.name)
    except Exception as e:  # noqa: BLE001 - the message is the point for a lab user
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    print(f"registered {entry['name']}: {entry['path']} (sha256 {entry['sha256'][:12]}..., "
          f"{entry['size'] // 1_000_000} MB) in {local_registry_path()}")
    if entry.get("card"):
        print(f"  card: {entry['card']}")
    print(f"  use it with  segmentation.nuclei.cellpose_model: {entry['name']}  in a profile, "
          f"or check it with  germquant get-model {entry['name']}")
    return 0
