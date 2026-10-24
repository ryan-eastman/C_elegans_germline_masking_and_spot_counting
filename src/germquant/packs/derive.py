"""`germquant pack derive`: promote a model trained in a project into a derived pack (docs/PACKS.md, "Derived
packs"; docs/DESIGN_system_packs.md, sections 3.5, 3.6, 3.8 and 10, step 8).

    germquant pack derive --from ID@VER --recipe R --model PATH --card CARD --out DIR [--id NEW_ID] [--json]

writes ``DIR/<new id>/0.1.0/``, a pack with one recipe and one new model, in the installed layout that
discovery reads (``<project>/packs/``, a $GERMQUANT_PACKS folder or the per-user folder):

* pack.yaml: the id ``<parent>_derived_<sha8>`` (the first 8 hex digits of the model's sha256) unless --id
  names another, version 0.1.0, status provisional, the parent's system, object kind, engine range and expected
  sizes, ``derived_from`` {pack, recipe, pack_sha256}, and the model spec ``<name>@<sha8>`` with its sha256,
  size, card and ``sources: []``: the weights come from the model store, or from ``pack fetch-models
  --from-file`` on another computer;
* ``config/profiles/<recipe>.yaml``: ``extends: pack:<parent>@<version>/<recipe>`` and
  ``segmentation.nuclei.cellpose_model: <name>@<sha8>``, nothing else; its config_hash is computed here and
  recorded as expected_config_hash; the recipe keeps the parent's needs (the new model in place of the one it
  replaces), tunables and dims;
* ``models/<name>.card.json``: the card, made portable: model_file becomes the file's name, base_model.path and
  every other absolute path a bare name (a sha256 beside them stays), a path quoted inside a text the same, and
  the JSON is written with ASCII escapes for any non-ASCII text; a ``derived`` record names the parent, the
  source card's sha256 and the keys rewritten;
* README.md (a short card) and SHA256SUMS.

The parent version is copied beside the derived pack (packs.vendor, as a new project copies it), unless DIR
already holds the parent (the copy `germquant project new` made, or one an earlier derive wrote): the parent is
then used where it is, nothing is written there, and the record's ``parent.copied`` is false. Either way the
derived recipe keeps resolving after an engine update moves the parent's built-in version (section 3.8). The
model goes into the model store at ``models_dir()/blobs/<sha256>/<file_name>``, with its verification record,
where the pack-model lookup and the pin check find it; the store's copy is a hard link to the training file
when the disk allows, so that file is to be replaced, never overwritten in place. The pack is built in a
staging folder, put through every check of packs.install_pack (install.check_installable), and installed
through it; the same inputs give the same bytes, so deriving again is a no-op. Nothing is written before every
check has passed. If a write itself fails (a full disk, a path too long for Windows), the model in the store
and the parent's copy may stay; the next derive reuses them.

Refused: a parent recipe that never runs a Cellpose model (segmentation.nuclei.method classical, labels from
files, or another segmentation backend); a parent recipe that does not load or no longer hashes to its expected
value; a card whose model_sha256 or model_size is not the file's; an id that breaks the pack id pattern or is
the parent's; a DIR that already holds the id and version, or the parent's, with other content or damaged; a
DIR inside the engine's own packs folder or inside a pack.

While deriving, a DIR named ``packs`` that discovery does not search anyway (not a $GERMQUANT_PACKS folder or
the per-user folder) is searched first, as a project's packs folder is while the project is open
(packs.use_project on its parent folder), so a parent that only the project holds is found there.

Torch-free (gate G6): this module reads YAML and JSON, hashes and copies files, and loads configs, like the rest
of germquant.packs.
"""
from __future__ import annotations

import contextlib
import hashlib
import json
import os
import re
import shutil
import tempfile
from pathlib import Path

import yaml

from ..fsutil import long_path
from .format import MANIFEST, PACK_ID_RE, SUMS_FILE, ModelSpec, Pack, load_pack_folder, sums_sha256, verify_sums

DERIVED_VERSION = "0.1.0"
DERIVED_STATUS = "provisional"
MODEL_KEY = "segmentation.nuclei.cellpose_model"
# the part of a pack model name before "@<8 hex>" (models_registry.PACK_MODEL_RE)
MODEL_BASE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")


class PackDeriveError(ValueError):
    """Inputs that pack derive refuses (the message says why and what to do)."""


# an absolute path quoted inside a text: a drive path (C:\x, C:/x), a UNC path (\\server\share\x), a POSIX path of
# two or more components (/home/x) or a home path (~/x); never a part of a URL, a word or a number
_TEXT_PATH_RE = re.compile(
    r"(?<![A-Za-z0-9])[A-Za-z]:[\\/][^\s\"'<>|,;()]*"
    r"|(?<![\w\\])\\\\[^\s\\\"'<>|,;()]+\\[^\s\"'<>|,;()]*"
    r"|(?<![\w.~:/\\-])/(?:[\w.-]+/)+[\w.-]*"
    r"|(?<![\w.~:/\\-])~[\\/][^\s\"'<>|,;()]*")


def _ascii(text) -> str:
    """Text as ASCII for a Markdown file: any other character written as a \\uXXXX escape, as JSON writes it."""
    return "".join(ch if ord(ch) < 128 else f"\\u{ord(ch):04x}" for ch in str(text))


def _strip_paths(text: str) -> str:
    """``text`` with every absolute path quoted in it replaced by its last component (a final '.' or ':' kept)."""
    def repl(m: re.Match) -> str:
        s = m.group(0)
        core = s.rstrip(".:")
        return _bare_name(core) + s[len(core):] if core else s

    return _TEXT_PATH_RE.sub(repl, text)


def _sha256(path: Path, *, lf: bool = False) -> str:
    """sha256 of a file; ``lf``: of a text file with CRLF turned into LF (the rule for every text digest here)."""
    h = hashlib.sha256()
    with open(long_path(path), "rb") as fh:
        if lf:
            h.update(fh.read().replace(b"\r\n", b"\n"))
        else:
            for block in iter(lambda: fh.read(1 << 24), b""):
                h.update(block)
    return h.hexdigest()


def _write(path: Path, text: str) -> None:
    """An ASCII text file with LF line endings (raises when the text is not ASCII: a bug guard)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(long_path(path), "wb") as fh:
        fh.write(text.encode("ascii"))


def _bare_name(value: str) -> str:
    """The last component of a path of any platform (C:/a/b, \\\\server\\share\\b, /home/a/b, ~/b -> b)."""
    parts = [p for p in re.split(r"[\\/]+", value.strip()) if p]
    return parts[-1] if parts else value.strip()


def _is_machine_path(value: str) -> bool:
    """An absolute path of any platform, or a home path (~/...): what pack lint refuses in a published card."""
    from .lint import _machine_path

    return _machine_path(value)


_PATH_END = r"(?=[.:]?(?:[\s\"'<>|,;()]|$))"          # a path ends at a delimiter (a final '.' or ':' allowed)


def _shorten(text: str, subs: dict[str, str]) -> str:
    """``text`` with the absolute paths in it shortened to their last component: first the paths that are whole
    values somewhere in the card (``subs``, longest first, and only as whole paths, so a path with spaces goes
    as one and a path never loses part of a longer one), then any absolute path left (_TEXT_PATH_RE)."""
    for old in sorted(subs, key=len, reverse=True):
        text = re.sub(r"(?<![\w.~:/\\-])" + re.escape(old) + _PATH_END, lambda m, new=subs[old]: new, text)
    return _strip_paths(text)


def _rewrite(obj, path: str, subs: dict[str, str], changed: list[str]):
    """``obj`` with every text shortened (_shorten); the dotted paths of the texts changed are appended to
    ``changed``. A text that is a whole absolute path becomes its bare name."""
    if isinstance(obj, dict):
        return {k: _rewrite(v, f"{path}.{k}" if path else str(k), subs, changed) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_rewrite(v, f"{path}[{i}]", subs, changed) for i, v in enumerate(obj)]
    if not isinstance(obj, str):
        return obj
    new = _shorten(obj, subs)
    if new != obj:
        changed.append(path)
    return new


def _machine_leaves(obj, path: str = ""):
    """``(dotted path, value)`` of every text that is, or quotes, an absolute path."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from _machine_leaves(v, f"{path}.{k}" if path else str(k))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _machine_leaves(v, f"{path}[{i}]")
    elif isinstance(obj, str) and (_is_machine_path(obj) or _TEXT_PATH_RE.search(obj)):
        yield path, obj


def portable_card(card: dict, *, name: str, file_name: str, sha256: str, size: int) -> tuple[dict, list[str]]:
    """``(card, rewritten)``: the card as a derived pack publishes it (module docstring) and the sorted keys
    whose values changed. ``name``, ``file_name``, ``sha256`` and ``size`` are the model's in the pack."""
    subs = {v: _bare_name(v) for _, v in _machine_leaves(card) if _is_machine_path(v)}
    changed: list[str] = []
    out = _rewrite(json.loads(json.dumps(card)), "", subs, changed)
    for key, value in (("name", name), ("model_file", file_name), ("model_sha256", sha256), ("model_size", size)):
        if out.get(key) != value:
            out[key] = value
            if key not in changed:
                changed.append(key)
    return out, sorted(set(changed))


def _read_card(path: Path) -> dict:
    try:
        with open(long_path(path), encoding="utf-8-sig") as fh:
            data = json.load(fh)
    except (OSError, ValueError) as e:
        raise PackDeriveError(f"--card {path}: not a readable JSON file ({e})") from e
    if not isinstance(data, dict):
        raise PackDeriveError(f"--card {path}: not a JSON object (a model card is one)")
    return data


def _model_base(card: dict, model: Path) -> str:
    """The name part of the pack model name: the card's name, else the model file's name."""
    for cand in (str(card.get("name") or "").strip(), model.name):
        if MODEL_BASE_RE.fullmatch(cand):
            return cand
    raise PackDeriveError(f"neither the card's name {card.get('name')!r} nor the file name {model.name!r} is a model "
                          f"name ({MODEL_BASE_RE.pattern}: letters, digits, '.', '_' and '-'); rename the model file")


def _check_card(card: dict, sha: str, size: int, card_path: Path, model: Path) -> None:
    got = str(card.get("model_sha256") or "").lower()
    if got and got != sha:
        raise PackDeriveError(f"the card's model_sha256 ({got[:12]}...) is not the sha256 of --model {model} "
                              f"({sha[:12]}...): {card_path} belongs to other weights")
    msize = card.get("model_size")
    if msize is not None and msize != size:
        raise PackDeriveError(f"the card's model_size ({msize!r}) is not the size of --model {model} ({size} bytes): "
                              f"{card_path} belongs to other weights")


def _guard_out(out: Path) -> None:
    """Refuse a DIR inside the engine's own packs folder (built-in packs only, tracked by git) or inside a pack."""
    from ..config import builtin_roots

    real = Path(os.path.realpath(out))
    for root in builtin_roots():
        own = Path(os.path.realpath(root / "packs"))
        if real == own or own in real.parents:
            raise PackDeriveError(f"--out {out}: the engine's own packs folder holds the built-in packs only; derive "
                                  f"into a project's packs folder (<project>/packs), a $GERMQUANT_PACKS folder or "
                                  f"the per-user folder")
    for folder in (real, *real.parents):
        if (folder / MANIFEST).is_file():
            raise PackDeriveError(f"--out {out}: inside the pack folder {folder}; name the folder that holds packs "
                                  f"(<project>/packs, for example)")


def _named_packs(out: Path) -> bool:
    """Whether DIR is named ``packs``, as a project's packs folder is (case-blind on Windows, as os.path.normcase
    is there)."""
    return os.path.normcase(out.name) == "packs"


@contextlib.contextmanager
def _searching(out: Path):
    """Discovery searches DIR first while deriving when DIR is a project's packs folder: named ``packs`` and not
    a folder discovery searches anyway ($GERMQUANT_PACKS or the per-user folder, whose packs under development,
    ``DIR/<id>/pack.yaml``, a project's folder would not list)."""
    if _named_packs(out) and _where(out) != "searched":
        from .discovery import use_project

        with use_project(out.parent):
            yield
    else:
        yield


def _find_parent(ref: str, out: Path) -> Pack:
    """The parent: ``DIR/<id>/<version>/`` when DIR already holds a copy of that version, else discovery's."""
    from .discovery import find_pack, parse_ref

    pid, version, recipe = parse_ref(ref)
    if recipe is not None:
        raise PackDeriveError(f"--from {ref!r}: name the pack as ID@VER and the recipe with --recipe")
    if version:
        here = out / pid / version
        if (here / MANIFEST).is_file():
            pack = load_pack_folder(here, source="path")
            bad = verify_sums(here)
            if bad:
                raise PackDeriveError(f"the copy of {pack.ref} in {here} does not match its SHA256SUMS ({bad[0]})")
            return pack
    return find_pack(f"{pid}@{version}" if version else pid)


def _parent_recipe(parent: Pack, recipe: str):
    from ..config import load_config
    from .discovery import recipe_profile_path

    if parent.problems:
        raise PackDeriveError(f"pack {parent.ref} is unusable: " + "; ".join(parent.problems))
    r = parent.recipe(recipe)
    try:
        cfg = load_config(recipe_profile_path(parent, r))
    except Exception as e:  # noqa: BLE001 - reported as the refusal
        raise PackDeriveError(f"pack {parent.ref} recipe {r.id} does not load ({type(e).__name__}: {e})") from e
    if cfg.hash != r.expected_config_hash:
        raise PackDeriveError(f"pack {parent.ref} recipe {r.id} hashes to {cfg.hash}, its pack expects "
                              f"{r.expected_config_hash}: reinstall the parent pack version")
    method = str(cfg.get("segmentation.nuclei.method", "auto"))
    backend = str(cfg.get("segmentation.nuclei.backend", "cellpose"))
    why = None
    if cfg.get("segmentation.nuclei.labels_from"):
        why = "takes its labels from files (segmentation.nuclei.labels_from)"
    elif method == "classical":
        why = "segments with the classical watershed (segmentation.nuclei.method: classical)"
    elif backend != "cellpose":
        why = f"segments with the {backend!r} backend, not Cellpose (segmentation.nuclei.backend)"
    if why:
        raise PackDeriveError(f"pack {parent.ref} recipe {r.id} {why}, so it never runs a Cellpose model; derive "
                              f"from a recipe that segments with Cellpose (method cellpose or auto)")
    return r, cfg


def _derived_id(parent: Pack, sha: str, pack_id: str | None) -> str:
    new = pack_id if pack_id is not None else f"{parent.id}_derived_{sha[:8]}"
    if not PACK_ID_RE.fullmatch(new):
        hint = " (the default id is too long: name one with --id)" if pack_id is None else ""
        raise PackDeriveError(f"--id {new!r} does not match {PACK_ID_RE.pattern} (a pack id: a lower-case letter, then "
                              f"lower-case letters, digits and underscores){hint}")
    if new == parent.id:
        raise PackDeriveError(f"--id {new!r} is the parent's id; a derived pack needs an id of its own")
    return new


def _replaced(name: str, spec: ModelSpec | None, value) -> bool:
    """Whether the parent's needed model ``name`` is the one its recipe's cellpose_model ``value`` names."""
    v = str(value or "").strip()
    if not v:
        return False
    return name == v or name == _bare_name(v) or (spec is not None and spec.file_name == _bare_name(v))


def _readme(new_ref: str, parent: Pack, recipe_id: str, model_name: str, sha: str, size: int, card_rel: str,
            config_hash: str, licence: str, others: list[str]) -> str:
    lines = [
        f"# {new_ref}: a derived pack",
        "",
        f"Made with `germquant pack derive` (docs/PACKS.md, \"Derived packs\"): the recipe `{recipe_id}` of the pack",
        f"`{parent.ref}`, with one change, the segmentation model `{model_name}`.",
        "",
        f"- Recipe `{recipe_id}`: `extends: pack:{parent.ref}/{recipe_id}` plus the model; config_hash `{config_hash}`.",
        f"- Model `{model_name}`: {size} bytes, sha256 `{sha}`; its card is `{card_rel}`.",
        f"- Parent `{parent.ref}`: a copy of it sits beside this pack, so the recipe keeps resolving when the engine",
        "  moves on to a newer version of the parent.",
        "- Status: provisional. This pack has no recipe golden and no validation record of its own: the parent's",
        "  recipe golden runs the classical segmenter, which never loads this model.",
    ]
    if others:
        lines.append(f"- Other models the recipe needs, as the parent's: {', '.join(others)}.")
    lines += [
        "",
        "The weights are not in the pack. On this computer they sit in the model store; on another computer, copy",
        "the model file and run",
        "",
        f"    germquant pack fetch-models {new_ref.split('@')[0]} --from-file {model_name}=PATH_TO_THE_MODEL_FILE",
        "",
        "(add `--project DIR` when the pack lives in a project's packs folder).",
    ]
    if licence:
        lines += ["", f"Licence of the model, from its card: {_ascii(licence)}"]
    return "\n".join(lines) + "\n"


def _store_model(model: Path, spec: ModelSpec, pack_ref: str) -> dict:
    """Put the model into the store at its blob path (a hard link when the disk allows, else a copy), check it,
    and record it; a blob already there and verified is kept. A hard link shares its bytes with the training
    file, which saves the model's size on disk (about 1.2 GB for a Cellpose-SAM fine-tune): train-nuclei
    replaces its file (os.replace), and docs/PACKS.md says never to overwrite that file in place, since the pin
    check would then refuse the store's copy too. The verification record only saves a re-hash, so
    a store that cannot keep it (read-only, or a path past Windows' length limit) costs that and nothing else,
    as for fetch-models (blobs._save_record)."""
    from .blobs import _place, _save_record, check_blob, model_location

    loc = model_location(spec)
    st = check_blob(spec)
    placed = "already in the store"
    if st["state"] not in ("trusted", "verified"):
        if loc.exists():
            loc.unlink()                              # a wrong blob under this digest is discarded, never kept
        placed = _place(Path(long_path(model)), loc, link=True)    # a project's training folder may be deep
        st = check_blob(spec)
        if st["state"] != "verified":
            loc.unlink(missing_ok=True)
            raise PackDeriveError(f"the model copied into the store at {loc} does not check out ({st['detail']})")
    source = None if placed == "already in the store" else f"pack derive {os.path.abspath(model)}"
    _save_record(spec, loc, source=source, packs=[pack_ref], verified=st["state"] == "verified")
    return {"path": str(loc), "state": st["state"], "placed": placed}


def _where(out: Path) -> str:
    """``searched`` (a root discovery always searches), ``project`` (a project's packs folder, searched while
    the project is open, or with --project) or ``not searched``."""
    from .discovery import pack_roots

    real = os.path.realpath(out)
    if any(os.path.realpath(r) == real for kind, r in pack_roots() if kind != "project"):
        return "searched"
    return "project" if _named_packs(out) else "not searched"


def derive_pack(parent: str, recipe: str | None, model: str | Path, card: str | Path, out: str | Path, *,
                pack_id: str | None = None) -> dict:
    """`germquant pack derive` (module docstring). ``parent`` is ``ID@VER`` (or ``ID``: the newest compatible
    version, which the derived recipe then pins), ``recipe`` one of its recipes (None: its default), ``model``
    the model file and ``card`` its card, ``out`` the packs folder to write into. Returns the record that
    ``--json`` prints: ``{"ok", "status": "derived" | "already present", "pack", "id", "version", "folder",
    "recipe", "extends", "config_hash", "parent": {"pack", "folder", "copied"}, "model": {"name", "file_name",
    "sha256", "size", "blob", "placed"}, "card", "rewritten", "where"}``. PackDeriveError (a ValueError) or
    FileNotFoundError when refused, before anything is written."""
    from ..config import load_config
    from ..models_registry import validate_card
    from .discovery import ROOT_BUILTIN, discover, refresh
    from .format import build_sums_text
    from .install import PackInstallError, check_installable, install_pack, vendor
    from .lint import _inside
    from .status import engine_files_problems

    out = Path(os.path.abspath(out))
    model, card_path = Path(model), Path(card)
    if not os.path.isfile(long_path(model)):                   # long_path: a training folder on a deep share
        raise FileNotFoundError(f"--model {model}: no such file")
    if not os.path.isfile(long_path(card_path)):
        raise FileNotFoundError(f"--card {card_path}: no such file")
    _guard_out(out)
    card_data = _read_card(card_path)
    sha, size = _sha256(model), os.stat(long_path(model)).st_size
    _check_card(card_data, sha, size, card_path, model)
    base = _model_base(card_data, model)
    model_name = f"{base}@{sha[:8]}"
    card_rel = f"models/{base}.card.json"
    with _searching(out):
        refresh()
        par = _find_parent(parent, out)
        r, pcfg = _parent_recipe(par, recipe)
        inside = _inside(par.folder, out)              # the parent already sits in DIR: nothing to copy
        bad = verify_sums(par.folder) + engine_files_problems(par)
        if bad:
            raise PackDeriveError(f"pack {par.ref} cannot be used as the parent: " + "; ".join(bad))
        new_id = _derived_id(par, sha, pack_id)
        new_ref = f"{new_id}@{DERIVED_VERSION}"
        if any(p.source == ROOT_BUILTIN and p.id == new_id for p in discover().packs):
            raise PackDeriveError(f"--id {new_id!r} is a built-in pack's id; a derived pack needs an id of its own")
        old_value = pcfg.get(MODEL_KEY)
        kept = [m for m in r.need("models") if not _replaced(m, par.model(m), old_value)]
        kept_specs = [par.model(m) for m in kept]      # declared: a pack whose needs name others is refused above
        if any(s.card == card_rel for s in kept_specs):
            raise PackDeriveError(f"the parent's model card {card_rel} has the new model's card name; rename the model")
        pcard, rewritten = portable_card(card_data, name=base, file_name=model.name, sha256=sha, size=size)
        licence = str(pcard.get("licence") or pcard.get("licence_note") or "see the model card")
        spec = ModelSpec(name=model_name, kind="cellpose", file_name=model.name, sha256=sha, size=size, sources=[],
                         card=card_rel, licence=licence)            # safe_dump writes YAML escapes for non-ASCII
        pcard["derived"] = {"pack": new_ref, "from": par.ref, "recipe": r.id,
                            "card_sha256": _sha256(card_path, lf=True), "rewritten": rewritten}
        check = validate_card(pcard, spec.as_dict(), pack_written=True)
        left = list(_machine_leaves(pcard))
        if check["problems"] or check["non_portable"] or left:
            raise PackDeriveError(f"the rewritten card is not portable: {check['problems'] + check['non_portable']} "
                                  f"{[k for k, _ in left]}")
        with tempfile.TemporaryDirectory(prefix="germquant-derive-") as tmp:
            stage = Path(tmp) / new_id
            recipe_rel = f"config/profiles/{r.id}.yaml"
            _write(stage / recipe_rel, yaml.safe_dump(
                {"extends": f"pack:{par.ref}/{r.id}", "segmentation": {"nuclei": {"cellpose_model": model_name}}},
                sort_keys=False, allow_unicode=False))
            _write(stage / card_rel, json.dumps(pcard, indent=1, ensure_ascii=True) + "\n")
            for s in kept_specs:                       # the parent's other models keep their cards
                if s.card and (par.folder / s.card).is_file():
                    (stage / s.card).parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(long_path(par.folder / s.card), long_path(stage / s.card))
            post = par.folder / "post_stages" / f"{r.id}.yaml"
            if post.is_file():                          # the post-stages the parent recipe switches on (design 5.6)
                (stage / "post_stages").mkdir(parents=True, exist_ok=True)
                shutil.copyfile(long_path(post), long_path(stage / "post_stages" / post.name))
            try:
                dcfg = load_config(stage / recipe_rel)
            except FileNotFoundError as e:
                raise PackDeriveError(f"the derived recipe does not resolve its parent ({e}); derive into a folder this "
                                      f"engine searches: a project's packs folder, a $GERMQUANT_PACKS folder or the "
                                      f"per-user folder") from e
            except Exception as e:  # noqa: BLE001 - reported as the refusal
                raise PackDeriveError(f"the derived recipe does not load ({type(e).__name__}: {e})") from e
            if str(dcfg.get(MODEL_KEY)) != model_name:
                raise PackDeriveError(f"the derived recipe resolves {MODEL_KEY} to {dcfg.get(MODEL_KEY)!r}")
            needs = dict(r.needs or {})
            needs["models"] = [model_name, *kept]
            manifest = {
                "format": 1, "id": new_id, "version": DERIVED_VERSION,
                "title": f"{par.title}: derived, model {model_name}",
                "summary": f"The recipe {r.id} of {par.ref} with the segmentation model {model_name}, made with "
                           f"germquant pack derive.",
                "system": par.data.get("system") or {}, "object_kind": par.object_kind,
                "maintainer": "local: made with germquant pack derive", "licence": par.data.get("licence", ""),
                "citation": {"text": (par.data.get("citation") or {}).get("text", ""),
                             "doi": (par.data.get("citation") or {}).get("doi"), "cff": None},
                "requires": par.data.get("requires") or {}, "root": "pack",
            }
            if par.data.get("expects"):
                manifest["expects"] = par.data["expects"]
            if par.data.get("filename_regex"):
                manifest["filename_regex"] = par.data["filename_regex"]
            manifest["derived_from"] = {"pack": par.ref, "recipe": r.id, "pack_sha256": sums_sha256(par.folder)}
            manifest["recipes"] = [{
                "id": r.id, "title": f"{r.title or r.id} (model {model_name})", "profile": recipe_rel, "dims": r.dims,
                "default": True, "status": DERIVED_STATUS, "expected_config_hash": dcfg.hash, "needs": needs,
                "tunables": list(r.tunables or []), "validation": None}]
            manifest["models"] = [spec.as_dict(), *(s.as_dict() for s in kept_specs)]
            head = (f"# A derived pack, made with germquant pack derive: the recipe {r.id} of {par.ref} with the "
                    f"model {model_name}.\n")
            _write(stage / MANIFEST, head + yaml.safe_dump(manifest, sort_keys=False, allow_unicode=False, width=110))
            _write(stage / "README.md", _readme(new_ref, par, r.id, model_name, sha, size, card_rel, dcfg.hash,
                                                licence, kept))
            _write(stage / SUMS_FILE, build_sums_text(stage))
            try:
                check_installable(stage)          # every install check (validate_pack too), before any write
            except PackInstallError as e:
                raise PackDeriveError(f"the derived pack would not install: {e}") from e
            dest = out / new_id / DERIVED_VERSION
            present = dest.exists()
            if present and sums_sha256(dest) != sums_sha256(stage):
                raise PackDeriveError(f"{out} already holds {new_ref} with other content (another model, card or "
                                      f"parent): derive with --id NEW_ID, or remove {dest}")
            if present and verify_sums(dest):
                raise PackDeriveError(f"the copy of {new_ref} at {dest} no longer matches its SHA256SUMS "
                                      f"({verify_sums(dest)[0]}): remove it and derive again")
            here = out / par.id / par.version          # where vendor copies the parent, unless DIR holds it
            if not inside and here.exists():           # vendor's install checks, before any write
                if sums_sha256(here) != sums_sha256(par.folder):
                    raise PackDeriveError(f"{here} holds something other than a copy of {par.ref}: move it away, or "
                                          f"derive into another folder")
                if verify_sums(here):
                    raise PackDeriveError(f"the copy of {par.ref} at {here} no longer matches its SHA256SUMS "
                                          f"({verify_sums(here)[0]}): remove it and derive again")
            # every check passed: the model, the parent's copy, then the derived pack
            stored = _store_model(model, spec, new_ref)
            parent_folder = par.folder if inside else vendor(par.id, par.version, out)
            installed = install_pack(stage, target=out)
    where = _where(out)                 # outside _searching: DIR is a project root only while it is searched first
    return {"ok": True, "status": "already present" if present else "derived", "pack": installed.ref,
            "id": installed.id, "version": installed.version, "folder": str(installed.folder), "recipe": r.id,
            "extends": f"pack:{par.ref}/{r.id}", "config_hash": dcfg.hash,
            "parent": {"pack": par.ref, "folder": str(parent_folder), "copied": not inside},
            "model": {"name": model_name, "file_name": model.name, "sha256": sha, "size": size,
                      "blob": stored["path"], "placed": stored["placed"]},
            "card": card_rel, "rewritten": rewritten, "where": where}
