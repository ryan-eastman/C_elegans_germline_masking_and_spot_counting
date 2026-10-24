"""post_stages.yaml: the post-stages a run switches on (docs/POST_STAGES.md; design 5.6; plugin need I13).

The file is a list of {name, params}:

    - name: ring_intensity
      params: {radius_um: 0.8}

It sits in a project's root next to recipe.yaml, ships in a pack as post_stages/<recipe>.yaml, or is given with
--post-stages FILE. load_config never reads it, so it changes neither config_hash, nor a stage hash, nor the
done marker; its digest enters each post-stage's sub-hash instead. An empty list switches every post-stage off
(and so replaces a pack's file).

Validation is against the post-stages installed in the environment that checks: an entry must name a usable
post-stage, give only declared parameters, and give each a value of its declared type and range. The plugin
therefore writes the file through ``python -m germquant.post_stages write`` in the engine environment (D39).
The writer is strict and canonical: ASCII with YAML escapes, LF line endings, the parameter keys sorted, float
parameters as floats, an atomic replace through fsutil.long_path. The digest of a file is taken on its
LF-normalised bytes, and a list of entries digests like the file this module would write from it.
"""
from __future__ import annotations

import hashlib
import os
from pathlib import Path

import yaml

from ..fsutil import long_path
from .declaration import DECLARATION_KEYS, declaration_problems, normalise_declaration, normalise_value, param_problem

POST_STAGES_FILE = "post_stages.yaml"
ENTRY_KEYS = ("name", "params")
_HEADER = ("# post_stages.yaml: the post-stages this run switches on (docs/POST_STAGES.md).\n"
           "# Not part of config_hash. Written by germquant; check entries with python -m germquant.post_stages.\n")


def write_text_atomic(path, text: str) -> Path:
    """Write ASCII ``text`` with LF line endings to ``path`` through a temporary file and os.replace (no handle
    is open during the replace), so a reader never sees half a file."""
    path = Path(path)
    tmp = path.with_name(path.name + ".part")
    try:
        with open(long_path(tmp), "w", encoding="ascii", newline="\n") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(long_path(tmp), long_path(path))
    finally:
        try:
            os.remove(long_path(tmp))
        except FileNotFoundError:
            pass
    return path


def entries_problems(entries) -> list[str]:
    """Problems with the shape of a post-stage list (a list of mappings with a name and optional params), before
    any declaration is consulted."""
    if isinstance(entries, (str, bytes)) or not isinstance(entries, (list, tuple)):
        return [f"post-stages must be a list of {{name, params}} entries, got {type(entries).__name__}"]
    out = []
    for i, e in enumerate(entries):
        if not isinstance(e, dict):
            out.append(f"entry {i} must be a mapping {{name, params}}, got {type(e).__name__}")
            continue
        unknown = sorted(str(k) for k in e if k not in ENTRY_KEYS)
        if unknown:
            out.append(f"entry {i}: unknown key(s) {unknown} (an entry takes name and params)")
        if "name" not in e:
            out.append(f"entry {i} has no name")
        elif not isinstance(e["name"], str) or not e["name"].strip():
            out.append(f"entry {i}: name must be the post-stage's name as text, got {e['name']!r}")
        params = e.get("params")
        if params is not None and not isinstance(params, dict):
            out.append(f"entry {i}: params must be a mapping of parameter keys to values, got "
                       f"{type(params).__name__}")
    return out


def _entry(e: dict) -> dict:
    return {"name": e["name"], "params": dict(e.get("params") or {})}


def load_post_stages_file(path) -> list[dict]:
    """The entries of a post_stages.yaml, as [{"name", "params"}] (an empty file is an empty list). ValueError
    when the file is not UTF-8 YAML or not a list of {name, params}; the names and values are checked by
    validate_post_stages."""
    p = Path(path)
    raw = Path(long_path(p)).read_bytes()
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as e:
        raise ValueError(f"{p} is not UTF-8 text (byte {e.start})") from None
    try:
        data = yaml.safe_load(text)
    except yaml.YAMLError as e:
        mark = getattr(e, "problem_mark", None)
        where = f" (line {mark.line + 1})" if mark is not None else ""
        raise ValueError(f"{p} cannot be read as YAML{where}") from None
    if data is None:
        return []
    probs = entries_problems(data)
    if probs:
        raise ValueError(f"{p}: " + "; ".join(probs))
    return [_entry(e) for e in data]


def as_declarations(declarations=None) -> tuple[dict[str, dict], dict[str, list[str]]]:
    """({name: declaration in listing form}, {name: why it is unusable}) for validation. ``declarations``: None
    for the post-stages of this environment; a mapping {name: declaration} (raw, as under post_stages:, or a
    listing row); or a list of listing rows, such as list_post_stages returns. ValueError for a given
    declaration with problems."""
    if declarations is None:
        from .discovery import declarations as _installed
        from .discovery import unusable_reasons

        return _installed(), unusable_reasons()
    if isinstance(declarations, (list, tuple)):
        declarations = {str(r["name"]): r for r in declarations}
    out = {}
    for name, decl in dict(declarations).items():
        raw = {k: v for k, v in dict(decl).items() if k in DECLARATION_KEYS}
        probs = declaration_problems(name, raw)
        if probs:
            raise ValueError("; ".join(probs))
        out[name] = normalise_declaration(name, raw)
    return out, {}


def problems_by_name(entries, decls: dict[str, dict], reasons: dict[str, list[str]] | None = None
                     ) -> list[tuple[str | None, str]]:
    """(post-stage name or None, problem) for every problem of a post-stage list against ``decls``: the shape,
    a name listed twice, a post-stage that is not installed or not usable, an undeclared parameter, a value
    outside its declared type or range."""
    shape = entries_problems(entries)
    if shape:
        return [(None, p) for p in shape]
    reasons = reasons or {}
    out: list[tuple[str | None, str]] = []
    seen: set[str] = set()
    for i, e in enumerate(entries):
        name = e["name"]
        if name in seen:
            out.append((name, f"entry {i}: post-stage {name} is listed more than once (one entry per post-stage)"))
            continue
        seen.add(name)
        decl = decls.get(name)
        if decl is None:
            if name in reasons:
                out.append((name, f"entry {i}: post-stage {name} is installed but not usable: "
                                  f"{'; '.join(reasons[name])}"))
            else:
                known = ", ".join(sorted(decls)) or "none"
                out.append((name, f"entry {i}: post-stage {name} is not installed in this environment "
                                  f"(installed: {known})"))
            continue
        specs = {s["key"]: s for s in decl.get("params") or []}
        for key, value in (e.get("params") or {}).items():
            if key not in specs:
                out.append((name, f"entry {i}: post-stage {name} has no parameter {key} (it declares "
                                  f"{sorted(specs) or 'none'})"))
                continue
            prob = param_problem(specs[key], value)
            if prob:
                out.append((name, f"entry {i}: post-stage {name}: {prob}"))
    return out


def validate_post_stages(entries, declarations=None) -> list[str]:
    """Every problem of a post-stage list against ``declarations`` (None: the post-stages installed in this
    environment; see as_declarations). Empty when the list can run here."""
    decls, reasons = as_declarations(declarations)
    return [p for _name, p in problems_by_name(entries, decls, reasons)]


def canonical_text(entries, decls: dict[str, dict]) -> str:
    """The canonical text of a well-formed post-stage list against resolved declarations (as_declarations):
    the parameter keys sorted, float parameters as floats, everything ASCII (yaml.safe_dump with allow_unicode
    False), after a two-line comment header."""
    norm = []
    for e in entries:
        specs = {s["key"]: s for s in (decls.get(e["name"]) or {}).get("params") or []}
        params = dict(e.get("params") or {})
        norm.append({"name": e["name"],
                     "params": {k: normalise_value(specs[k], params[k]) if k in specs else params[k]
                                for k in sorted(params)}})
    body = yaml.safe_dump(norm, allow_unicode=False, sort_keys=False, default_flow_style=False) if norm else "[]\n"
    return _HEADER + body


def dump_post_stages(entries, declarations=None) -> str:
    """What write_post_stages_file writes for a well-formed post-stage list (see canonical_text)."""
    return canonical_text([_entry(e) for e in entries], as_declarations(declarations)[0])


def write_post_stages_file(path, entries, *, declarations=None) -> Path:
    """Validate ``entries`` (against this environment's post-stages unless ``declarations`` is given) and write
    them to ``path`` in the canonical form. ValueError naming every problem, and nothing written, when the list
    cannot run here."""
    decls, reasons = as_declarations(declarations)
    probs = [p for _name, p in problems_by_name(entries, decls, reasons)]
    if probs:
        raise ValueError(f"{POST_STAGES_FILE} not written: " + "; ".join(probs))
    return write_text_atomic(path, canonical_text([_entry(e) for e in entries], decls))


def text_sha256(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def file_sha256(path) -> str:
    """The sha256 of a file's LF-normalised bytes (the digest of a post_stages.yaml)."""
    return hashlib.sha256(Path(long_path(Path(path))).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def post_stages_sha256(source, declarations=None) -> str:
    """The digest that enters every post-stage's sub-hash: of a file (a path), its LF-normalised bytes; of a list
    of entries, the canonical text write_post_stages_file would write from it."""
    if isinstance(source, (str, os.PathLike)):
        return file_sha256(source)
    return text_sha256(dump_post_stages(source, declarations))
