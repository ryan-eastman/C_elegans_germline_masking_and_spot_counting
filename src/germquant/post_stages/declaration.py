"""The post-stage declaration, germquant_post_stage.yaml, and its parameter values (docs/POST_STAGES.md).

A post-stage declares itself in a file inside the top-level package its entry point names, so the engine reads
what it reads, writes and takes without importing it (decision D31):

    format: 1
    post_stages:
      <entry point name>:
        version: "1.0.0"
        title: "one line for a listing"                                       # optional
        reads: {roles: [foci], tables: [nuclei, spots], labels: true, image: true}
        writes: {<table>: [<column>, ...]}
        params: [{key: radius_um, type: float, unit: um, min: 0.1, max: 5.0, default: 1.0, label: "Radius"}]

Parameters are declared like recipe tunables (key, type, unit, range), so the engine gains no pydantic. The
types are float, int, bool, choice (with choices) and str; every parameter has a default. Values are checked
strictly: a bool is never a number, a choice keeps its YAML type, a float accepts an int and is stored as a
float. Names of post-stages, tables and parameters are lower-case snake case (no double or trailing
underscore), so the file names <image_id>__post_<name>_<table>.csv never look like a main table's.

Standard library only: listing post-stages never loads numpy, pandas or a backend.
"""
from __future__ import annotations

import math
import re

DECLARATION_FILE = "germquant_post_stage.yaml"
DECLARATION_FORMAT = 1
# lower-case snake case: no double underscore (file names split on "__"), no leading or trailing one
NAME_RE = re.compile(r"^[a-z][a-z0-9]*(?:_[a-z0-9]+)*$")
MAX_NAME = 40
COLUMN_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,63}$")
ROLE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,63}$")
# the semantic-version pattern of germquant.packs.format.SEMVER_RE (copied: this module imports no package)
SEMVER_RE = re.compile(r"^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(?:-([0-9A-Za-z.-]+))?(?:\+([0-9A-Za-z.-]+))?$")
PARAM_TYPES = ("float", "int", "bool", "choice", "str")
DECLARATION_KEYS = ("version", "title", "reads", "writes", "params")
READS_KEYS = ("roles", "tables", "labels", "image")
PARAM_KEYS = ("key", "type", "default", "unit", "min", "max", "choices", "label")


def _is_number(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _is_scalar(v) -> bool:
    return isinstance(v, (str, int, float, bool))


def name_problem(name, what: str = "post-stage name") -> str | None:
    """Why ``name`` cannot name a post-stage, a table or a parameter (None when it can)."""
    if not isinstance(name, str) or not NAME_RE.fullmatch(name) or len(name) > MAX_NAME:
        return (f"{what} {name!r} must be lower-case snake case of at most {MAX_NAME} characters "
                f"(letters, digits and single underscores, starting with a letter)")
    return None


def param_problem(spec: dict, value) -> str | None:
    """Why ``value`` is not a valid value of the declared parameter ``spec`` (None when it is)."""
    kind = spec.get("type")
    key = spec.get("key")
    if kind == "float":
        if not _is_number(value) or not math.isfinite(float(value)):
            return f"parameter {key} must be a number, got {value!r}"
    elif kind == "int":
        if not isinstance(value, int) or isinstance(value, bool):
            return f"parameter {key} must be an integer, got {value!r}"
    elif kind == "bool":
        if not isinstance(value, bool):
            return f"parameter {key} must be true or false, got {value!r}"
    elif kind == "choice":
        if not any(type(c) is type(value) and c == value for c in spec.get("choices") or []):
            return f"parameter {key} must be one of {list(spec.get('choices') or [])}, got {value!r}"
    elif kind == "str":
        if not isinstance(value, str):
            return f"parameter {key} must be text, got {value!r}"
    else:
        return f"parameter {key} has the unknown type {kind!r}"
    lo, hi = spec.get("min"), spec.get("max")
    if kind in ("float", "int"):
        if lo is not None and value < lo:
            return f"parameter {key} must be at least {lo}, got {value!r}"
        if hi is not None and value > hi:
            return f"parameter {key} must be at most {hi}, got {value!r}"
    return None


def normalise_value(spec: dict, value):
    """A valid value as it is recorded and hashed: a float parameter given as an int becomes a float."""
    if spec.get("type") == "float" and _is_number(value):
        return float(value)
    return value


def _param_spec_problems(i: int, spec) -> list[str]:
    if not isinstance(spec, dict):
        return [f"params[{i}] must be a mapping with key, type and default"]
    out = []
    key = spec.get("key")
    where = f"params[{i}] ({key})" if isinstance(key, str) else f"params[{i}]"
    unknown = sorted(str(k) for k in spec if k not in PARAM_KEYS)
    if unknown:
        out.append(f"{where}: unknown key(s) {unknown} (a parameter takes {list(PARAM_KEYS)})")
    prob = name_problem(key, "parameter key")
    if prob:
        out.append(f"{where}: {prob}")
    kind = spec.get("type")
    if kind not in PARAM_TYPES:
        out.append(f"{where}: type {kind!r} is not one of {list(PARAM_TYPES)}")
        return out
    for bound in ("min", "max"):
        if bound not in spec:
            continue
        v = spec[bound]
        if kind not in ("float", "int"):
            out.append(f"{where}: {bound} applies to float and int parameters only")
        elif kind == "float" and not (_is_number(v) and math.isfinite(float(v))):
            out.append(f"{where}: {bound} must be a number, got {v!r}")
        elif kind == "int" and (not isinstance(v, int) or isinstance(v, bool)):
            out.append(f"{where}: {bound} must be an integer, got {v!r}")
    lo, hi = spec.get("min"), spec.get("max")
    if _is_number(lo) and _is_number(hi) and lo > hi:
        out.append(f"{where}: min {lo} is above max {hi}")
    if kind == "choice":
        ch = spec.get("choices")
        if not isinstance(ch, list) or not ch or not all(_is_scalar(c) for c in ch):
            out.append(f"{where}: a choice parameter needs choices, a non-empty list of values")
        elif len({(type(c).__name__, c) for c in ch}) != len(ch):
            out.append(f"{where}: choices lists a value twice")
    elif "choices" in spec:
        out.append(f"{where}: choices applies to choice parameters only")
    for text in ("unit", "label"):
        if text in spec and not isinstance(spec[text], str):
            out.append(f"{where}: {text} must be text, got {spec[text]!r}")
    if "default" not in spec:
        out.append(f"{where}: no default (every parameter needs one)")
    elif not out:
        prob = param_problem(spec, spec["default"])
        if prob:
            out.append(f"{where}: the default is not valid ({prob})")
    return out


def _name_list(value, what: str, pattern, *, allow_post: bool = True) -> list[str]:
    if not isinstance(value, list):
        return [f"{what} must be a list"]
    out = []
    for v in value:
        if pattern is NAME_RE:
            prob = name_problem(v, f"{what} table name")
            if prob:
                out.append(prob)
            elif not allow_post and v.startswith("post_"):
                out.append(f"{what} names {v!r}: a post-stage reads the main run's tables, never post_ tables")
        elif not isinstance(v, str) or not pattern.fullmatch(v):
            out.append(f"{what} names {v!r}, which is not a role name")
    if len(set(map(str, value))) != len(value):
        out.append(f"{what} lists a name twice")
    return out


def declaration_problems(name: str, decl) -> list[str]:
    """Every problem with the declaration ``decl`` of post-stage ``name`` (one entry under post_stages:), each
    naming the post-stage. Empty when it is sound."""
    out = []
    prob = name_problem(name)
    if prob:
        out.append(prob)
    if not isinstance(decl, dict):
        return [f"post-stage {name}: the declaration must be a mapping"] + out
    unknown = sorted(str(k) for k in decl if k not in DECLARATION_KEYS)
    if unknown:
        out.append(f"unknown key(s) {unknown} (a declaration takes {list(DECLARATION_KEYS)})")
    version = decl.get("version")
    if not isinstance(version, str) or not SEMVER_RE.fullmatch(version):
        out.append(f"version {version!r} must be a quoted semantic version such as \"1.0.0\"")
    if "title" in decl and not isinstance(decl["title"], str):
        out.append("title must be text")
    reads = decl.get("reads", {})
    if not isinstance(reads, dict):
        out.append("reads must be a mapping {roles, tables, labels, image}")
        reads = {}
    unknown = sorted(str(k) for k in reads if k not in READS_KEYS)
    if unknown:
        out.append(f"reads: unknown key(s) {unknown} (reads takes {list(READS_KEYS)})")
    out += _name_list(reads.get("roles", []), "reads.roles", ROLE_RE)
    out += _name_list(reads.get("tables", []), "reads.tables", NAME_RE, allow_post=False)
    for flag in ("labels", "image"):
        if not isinstance(reads.get(flag, False), bool):
            out.append(f"reads.{flag} must be true or false")
    roles = reads.get("roles", [])
    if reads.get("image", False) is True and not roles:
        out.append("reads.image is true but reads.roles names no channel role to read")
    if roles and reads.get("image", False) is not True:
        out.append("reads.roles names channel roles but reads.image is not true")
    writes = decl.get("writes")
    if not isinstance(writes, dict) or not writes:
        out.append("writes must map each table the post-stage writes to its columns")
        writes = {}
    for table, cols in writes.items():
        prob = name_problem(table, "writes: table name")
        if prob:
            out.append(prob)
        if not isinstance(cols, list) or not cols:
            out.append(f"writes.{table}: a table needs a non-empty list of column names")
            continue
        bad = [c for c in cols if not isinstance(c, str) or not COLUMN_RE.fullmatch(c)]
        if bad:
            out.append(f"writes.{table}: column name(s) {bad} must start with a letter and hold letters, digits "
                       f"and underscores only")
        if len(set(map(str, cols))) != len(cols):
            out.append(f"writes.{table} lists a column twice")
    params = decl.get("params", [])
    if not isinstance(params, list):
        out.append("params must be a list of parameter mappings")
        params = []
    seen = set()
    for i, spec in enumerate(params):
        out += _param_spec_problems(i, spec)
        key = spec.get("key") if isinstance(spec, dict) else None
        if isinstance(key, str):
            if key in seen:
                out.append(f"params: key {key} is declared twice")
            seen.add(key)
    return [p if p.startswith(f"post-stage {name}") else f"post-stage {name}: {p}" for p in out]


def normalise_declaration(name: str, decl: dict) -> dict:
    """The listing form of a sound declaration: every reads key present, params as ordered mappings."""
    reads = dict(decl.get("reads") or {})
    params = []
    for spec in decl.get("params") or []:
        row = {"key": spec["key"], "type": spec["type"], "default": normalise_value(spec, spec["default"])}
        for k in ("unit", "min", "max", "choices", "label"):
            if k in spec:
                row[k] = list(spec[k]) if k == "choices" else spec[k]
        params.append(row)
    return {"name": name, "version": decl["version"], "title": str(decl.get("title", "")),
            "reads": {"roles": [str(r) for r in reads.get("roles", [])],
                      "tables": [str(t) for t in reads.get("tables", [])],
                      "labels": bool(reads.get("labels", False)), "image": bool(reads.get("image", False))},
            "writes": {str(t): [str(c) for c in cols] for t, cols in (decl.get("writes") or {}).items()},
            "params": params}


def declaration_from_document(doc, name: str) -> tuple[dict | None, list[str]]:
    """The declaration of post-stage ``name`` in a parsed germquant_post_stage.yaml, and the problems (format
    1, a post_stages mapping, an entry for ``name``, and that entry's own problems)."""
    if not isinstance(doc, dict):
        return None, [f"{DECLARATION_FILE} must be a mapping with format and post_stages"]
    probs = []
    unknown = sorted(str(k) for k in doc if k not in ("format", "post_stages"))
    if unknown:
        probs.append(f"{DECLARATION_FILE}: unknown key(s) {unknown} (it takes format and post_stages)")
    fmt = doc.get("format")
    if isinstance(fmt, bool) or fmt != DECLARATION_FORMAT:
        probs.append(f"{DECLARATION_FILE}: format {fmt!r} is not format {DECLARATION_FORMAT}")
    stages = doc.get("post_stages")
    if not isinstance(stages, dict):
        probs.append(f"{DECLARATION_FILE}: post_stages must map each entry point name to its declaration")
        return None, probs
    if name not in stages:
        declared = sorted(str(k) for k in stages)
        probs.append(f"{DECLARATION_FILE} declares no post-stage named {name!r} (it declares {declared})")
        return None, probs
    decl = stages[name]
    probs += declaration_problems(name, decl)
    return (None, probs) if probs else (decl, [])


def resolve_params(declaration: dict, given: dict | None) -> dict:
    """The parameters a post-stage receives: every declared default, overridden by the valid ``given`` values
    (checked beforehand with param_problem), floats as floats. ``declaration`` is the listing form."""
    out = {}
    given = dict(given or {})
    for spec in declaration.get("params") or []:
        key = spec["key"]
        out[key] = normalise_value(spec, given[key] if key in given else spec["default"])
    return out
