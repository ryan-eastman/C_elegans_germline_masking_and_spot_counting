"""Pack format 1: the manifest (pack.yaml), its validation, and SHA256SUMS (docs/DESIGN_system_packs.md, 3.2 to 3.7).

A pack is a folder of static files. ``pack.yaml`` names the recipes (each a profile file that load_config
reads like any other), the models by name and digest (weights never live in a pack), the examples and the
card. ``SHA256SUMS`` lists ``<sha256>  <relative/path>`` for every other file, sorted, with forward slashes
and LF; text files are hashed after CRLF to LF normalisation, so a Windows checkout (core.autocrlf=true) and
a Linux checkout agree. Nothing here imports more than yaml, hashlib and the standard library.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import stat
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath

import yaml

FORMAT_MAJOR = 1
MANIFEST = "pack.yaml"
SUMS_FILE = "SHA256SUMS"
PACK_ID_RE = re.compile(r"^[a-z][a-z0-9_]{1,40}$")
RECIPE_ID_RE = re.compile(r"^[a-z][a-z0-9_]{1,40}$")
SEMVER_RE = re.compile(r"^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(?:-([0-9A-Za-z.-]+))?(?:\+([0-9A-Za-z.-]+))?$")
HASH12_RE = re.compile(r"^[0-9a-f]{12}$")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
TEXT_SUFFIXES = (".yaml", ".yml", ".json", ".md", ".csv", ".txt", ".cff")
DIMS = ("3d", "2d", "projected", "3d+t")
STATUSES = ("provisional", "validated")
OBJECT_KINDS = ("nucleus", "cell")
ROOT_KINDS = ("pack", "checkout")
ACCELERATOR = ("optional", "recommended")
SOURCE_TYPES = ("cellpose_cache", "url", "release")
# junk: never hashed, never listed, never copied (build-sums, verify and install skip it; .gitignore keeps it
# out of git, design 3.7). SHA256SUMS itself is not hashed either, but it is part of the pack.
SKIP_DIRS = frozenset({"__pycache__"})
SKIP_NAMES = frozenset({".DS_Store"})
SKIP_SUFFIXES = (".pyc", ".pyo", ".pyd")
# where an outside recipe lives (design 3.2): config/profiles/ anchors its relative paths to the pack folder
PACK_PROFILE_RE = re.compile(r"^config/profiles/[^/]+\.ya?ml$")


class PackFormatError(ValueError):
    """A pack folder whose manifest cannot be read at all."""


# ---- digests ---------------------------------------------------------------------------------------------

def is_text(rel: str) -> bool:
    return PurePosixPath(rel).suffix.lower() in TEXT_SUFFIXES


def normalised_bytes(path: Path, rel: str | None = None) -> bytes:
    """The bytes a digest is taken on: text files (TEXT_SUFFIXES) with CRLF turned into LF, others as they are."""
    data = Path(path).read_bytes()
    return data.replace(b"\r\n", b"\n") if is_text(rel or Path(path).name) else data


def file_digest(path: Path, rel: str | None = None) -> str:
    return hashlib.sha256(normalised_bytes(path, rel)).hexdigest()


def lf_sha256(path: Path) -> str:
    """sha256 of a text file with CRLF turned into LF (the digest a root: checkout pack records per engine file)."""
    return hashlib.sha256(Path(path).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def skipped(rel_parts: tuple[str, ...]) -> bool:
    """Whether a path inside a pack is junk: a cache folder, compiled Python, .DS_Store (never hashed, listed
    or copied)."""
    if any(p in SKIP_DIRS for p in rel_parts[:-1]):
        return True
    name = rel_parts[-1]
    return name in SKIP_NAMES or name.endswith(SKIP_SUFFIXES)


def is_link(path: Path) -> bool:
    """Whether ``path`` is a symlink or a Windows directory junction. Python 3.11 reports a junction as neither
    a symlink nor a link (Path.is_symlink, os.path.islink) and walks into it, and creating one needs no
    privilege, so a junction inside a pack would pull files from outside it into the sums and an install."""
    p = Path(path)
    try:
        if p.is_symlink():
            return True
        isjunction = getattr(os.path, "isjunction", None)            # Python 3.12 and later
        if isjunction is not None:
            return bool(isjunction(p))
        if os.name == "nt":
            return getattr(os.lstat(p), "st_reparse_tag", 0) == stat.IO_REPARSE_TAG_MOUNT_POINT
    except OSError:
        return False
    return False


def _walk(folder: Path) -> tuple[list[str], list[str]]:
    """(files, links): sorted relative paths (forward slashes) of the plain files under ``folder`` (junk and
    SHA256SUMS left out) and of the symlinks and junctions, which are never followed."""
    files, links = [], []
    for dirpath, dirnames, filenames in os.walk(folder, followlinks=False):
        d = Path(dirpath)
        for name in list(dirnames):
            if is_link(d / name):
                dirnames.remove(name)                  # never walked into: its files are not the pack's
                links.append((d / name).relative_to(folder).as_posix())
        for name in filenames:
            rel = (d / name).relative_to(folder)
            if skipped(rel.parts) or rel.as_posix() == SUMS_FILE:
                continue
            (links if is_link(d / name) else files).append(rel.as_posix())
    return sorted(files), sorted(links)


def pack_files(folder: Path) -> list[str]:
    """Sorted relative paths (forward slashes) of every file of the pack that SHA256SUMS covers: every plain
    file except SHA256SUMS itself and junk. Symlinks and junctions are neither listed nor followed
    (pack_links names them, and verify_sums reports them)."""
    return _walk(Path(folder))[0]


def pack_links(folder: Path) -> list[str]:
    """Relative paths of the symlinks and Windows junctions inside a pack folder (a pack holds plain files only)."""
    return _walk(Path(folder))[1]


def build_sums_text(folder: Path) -> str:
    folder = Path(folder)
    lines = [f"{file_digest(folder / rel, rel)}  {rel}" for rel in pack_files(folder)]
    return "\n".join(lines) + ("\n" if lines else "")


def write_sums(folder: Path) -> Path:
    """(Re)write ``<folder>/SHA256SUMS``; returns its path."""
    p = Path(folder) / SUMS_FILE
    p.write_bytes(build_sums_text(folder).encode("ascii"))
    return p


def read_sums(folder: Path) -> dict[str, str]:
    """``{relative path: sha256}`` from SHA256SUMS; PackFormatError when it is missing or malformed."""
    p = Path(folder) / SUMS_FILE
    if not p.is_file():
        raise PackFormatError(f"{folder}: no {SUMS_FILE}")
    out: dict[str, str] = {}
    for i, line in enumerate(p.read_bytes().replace(b"\r\n", b"\n").decode("ascii", "replace").split("\n"), 1):
        if not line.strip():
            continue
        parts = line.split("  ", 1)
        if len(parts) != 2 or not SHA256_RE.fullmatch(parts[0]) or not parts[1]:
            raise PackFormatError(f"{p}: line {i} is not '<sha256>  <path>': {line[:80]!r}")
        out[parts[1]] = parts[0]
    return out


def sums_sha256(folder: Path) -> str | None:
    """pack_sha256: the sha256 of the LF-normalised SHA256SUMS text, or None without one."""
    p = Path(folder) / SUMS_FILE
    if not p.is_file():
        return None
    return hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def verify_sums(folder: Path) -> list[str]:
    """Problems with SHA256SUMS: missing or malformed, a listed file missing or changed, an unlisted file
    present, a path that is not a plain relative path, a symlink or junction inside the pack. Empty when the
    pack matches its sums."""
    folder = Path(folder)
    try:
        sums = read_sums(folder)
    except PackFormatError as e:
        return [str(e)]
    files, links = _walk(folder)
    problems = [f"{rel}: a symlink or junction; a pack holds plain files only" for rel in links]
    for rel, want in sorted(sums.items()):
        pp = PurePosixPath(rel)
        if pp.is_absolute() or ".." in pp.parts or "\\" in rel or re.match(r"^[A-Za-z]:", rel):
            problems.append(f"{SUMS_FILE} lists a path that is not a plain relative path: {rel!r}")
            continue
        f = folder / rel
        if not f.is_file():
            problems.append(f"{rel}: listed in {SUMS_FILE} but missing")
        elif any(rel == ln or rel.startswith(ln + "/") for ln in links):
            continue                                   # reached through a link: reported above
        elif file_digest(f, rel) != want:
            problems.append(f"{rel}: content differs from {SUMS_FILE}")
    for rel in files:
        if rel not in sums:
            problems.append(f"{rel}: not listed in {SUMS_FILE} (run `germquant pack build-sums`)")
    return problems


# ---- the manifest ----------------------------------------------------------------------------------------

@dataclass
class ModelSpec:
    name: str
    kind: str
    file_name: str
    sha256: str
    size: int
    sources: list = field(default_factory=list)
    card: str | None = None
    licence: str = ""

    @property
    def base_name(self) -> str:
        return self.name.split("@")[0]

    def as_dict(self) -> dict:
        return {"name": self.name, "kind": self.kind, "file_name": self.file_name, "sha256": self.sha256,
                "size": self.size, "sources": list(self.sources), "card": self.card, "licence": self.licence}


@dataclass
class Recipe:
    id: str
    profile: str
    title: str = ""
    dims: str = "3d"
    default: bool = False
    status: str = "provisional"
    expected_config_hash: str = ""
    needs: dict = field(default_factory=dict)
    tunables: list = field(default_factory=list)
    validation: str | None = None
    raw: dict = field(default_factory=dict)

    def need(self, key: str) -> list:
        v = (self.needs or {}).get(key) or []
        return [str(x) for x in v] if isinstance(v, list) else [str(v)]

    def as_dict(self) -> dict:
        return {"id": self.id, "title": self.title, "profile": self.profile, "dims": self.dims,
                "default": self.default, "status": self.status, "expected_config_hash": self.expected_config_hash,
                "needs": dict(self.needs or {}), "tunables": list(self.tunables or []), "validation": self.validation}


@dataclass
class Pack:
    id: str
    version: str
    folder: Path
    data: dict
    source: str = "path"
    recipes: list[Recipe] = field(default_factory=list)
    models: list[ModelSpec] = field(default_factory=list)
    problems: list[str] = field(default_factory=list)

    @property
    def ref(self) -> str:
        return f"{self.id}@{self.version}"

    @property
    def root(self) -> str:
        return str(self.data.get("root", "pack"))

    @property
    def title(self) -> str:
        return str(self.data.get("title", self.id))

    @property
    def object_kind(self) -> str:
        return str(self.data.get("object_kind", "nucleus"))

    @property
    def usable_manifest(self) -> bool:
        return not self.problems

    def default_recipe(self) -> Recipe:
        flagged = [r for r in self.recipes if r.default]
        if flagged:
            return flagged[0]
        if not self.recipes:
            raise FileNotFoundError(f"pack {self.ref} declares no recipe")
        return self.recipes[0]

    def recipe(self, recipe_id: str | None = None) -> Recipe:
        """The recipe ``recipe_id``, or the default recipe for None. FileNotFoundError naming the recipes."""
        if recipe_id is None or str(recipe_id) == "":
            return self.default_recipe()
        for r in self.recipes:
            if r.id == recipe_id:
                return r
        raise FileNotFoundError(f"pack {self.ref} has no recipe {recipe_id!r}; its recipes: "
                                f"{', '.join(r.id for r in self.recipes) or 'none'}")

    def model(self, name: str) -> ModelSpec | None:
        return next((m for m in self.models if m.name == name), None)

    def summary(self) -> dict:
        """What `pack list --json` prints for the pack."""
        return {"id": self.id, "version": self.version, "title": self.title,
                "summary": self.data.get("summary", ""), "system": self.data.get("system", {}),
                "object_kind": self.object_kind, "root": self.root, "source": self.source,
                "folder": str(self.folder), "requires": self.data.get("requires", {}),
                "licence": self.data.get("licence", ""), "citation": self.data.get("citation", {}),
                "recipes": [r.as_dict() for r in self.recipes], "models": [m.as_dict() for m in self.models],
                "problems": list(self.problems)}


def load_manifest(folder: Path) -> dict:
    """The parsed pack.yaml (yaml.safe_load; ASCII required). PackFormatError when unreadable."""
    p = Path(folder) / MANIFEST
    try:
        raw = p.read_bytes()
    except OSError as e:
        raise PackFormatError(f"{p}: cannot be read ({e})") from e
    try:
        text = raw.decode("ascii")
    except UnicodeDecodeError as e:
        raise PackFormatError(f"{p}: byte {e.start} is not ASCII (pack files are ASCII only)") from e
    try:
        data = yaml.safe_load(text)
    except yaml.YAMLError as e:
        raise PackFormatError(f"{p}: not valid YAML ({str(e).splitlines()[0]})") from e
    if not isinstance(data, dict):
        raise PackFormatError(f"{p}: not a mapping")
    return data


def _format_major(value) -> int | None:
    try:
        return int(str(value).split(".")[0])
    except (TypeError, ValueError):
        return None


def _is_plain_relative(path: str) -> bool:
    pp = PurePosixPath(str(path))
    return (bool(str(path)) and not pp.is_absolute() and ".." not in pp.parts and "\\" not in str(path)
            and not re.match(r"^[A-Za-z]:", str(path)))


def validation_record_passes(folder: Path, rel: str | None) -> tuple[bool, str]:
    """Whether a recipe's validation record exists and records a passing gate (verdict pass)."""
    if not rel:
        return False, "no validation record named"
    p = Path(folder) / rel
    if not _is_plain_relative(rel) or not p.is_file():
        return False, f"validation record {rel!r} not found"
    try:
        rec = json.loads(p.read_text(encoding="ascii"))
    except (OSError, ValueError) as e:
        return False, f"validation record {rel!r} unreadable ({e})"
    verdict = rec.get("verdict") if isinstance(rec, dict) else None
    if isinstance(verdict, dict):
        verdict = "pass" if verdict.get("pass") is True else verdict.get("verdict")
    return (str(verdict).lower() == "pass"), f"validation record {rel!r} verdict {verdict!r}"


def validate_manifest(data: dict, folder: Path | None = None) -> list[str]:
    """Problems with a pack.yaml (format 1). Refused: an unknown format major; a pack id, recipe id or model
    name that breaks its pattern; a model without sha256 and size, or whose 8-hex name suffix differs from the
    start of its sha256; a release model that is not the registry's; an outside recipe (root: pack) whose
    profile is not config/profiles/<recipe>.yaml, the layout that anchors its relative paths to the pack folder
    (config._base_dir_for); status validated (or validated: true) without a passing validation record (checked
    when ``folder`` is given). The tunables need the recipes loaded: packs.validate_pack checks them."""
    from ..models_registry import MODELS, PACK_MODEL_RE

    probs: list[str] = []
    major = _format_major(data.get("format"))
    if major != FORMAT_MAJOR:
        probs.append(f"format {data.get('format')!r} is not format {FORMAT_MAJOR} (this engine reads format "
                     f"{FORMAT_MAJOR} packs)")
        return probs
    pid = str(data.get("id", ""))
    if not PACK_ID_RE.fullmatch(pid):
        probs.append(f"pack id {pid!r} does not match {PACK_ID_RE.pattern}")
    version = str(data.get("version", ""))
    if not SEMVER_RE.fullmatch(version):
        probs.append(f"version {version!r} is not a semantic version (MAJOR.MINOR.PATCH)")
    if str(data.get("object_kind", "nucleus")) not in OBJECT_KINDS:
        probs.append(f"object_kind {data.get('object_kind')!r} is not one of {OBJECT_KINDS}")
    root = str(data.get("root", "pack"))
    if root not in ROOT_KINDS:
        probs.append(f"root {root!r} is not one of {ROOT_KINDS}")
    req = data.get("requires") or {}
    if not isinstance(req, dict):
        probs.append("requires is not a mapping")
        req = {}
    if req.get("germquant") is not None:
        try:
            from packaging.specifiers import SpecifierSet

            SpecifierSet(str(req["germquant"]))
        except Exception as e:  # noqa: BLE001 - any parse failure is a manifest problem
            probs.append(f"requires.germquant {req['germquant']!r} is not a version range ({e})")
    expects = data.get("expects") or {}
    for key in ("pixel_um", "z_step_um"):
        rng = expects.get(key) if isinstance(expects, dict) else None
        if rng is not None and not (isinstance(rng, list) and len(rng) == 2
                                    and all(isinstance(v, (int, float)) for v in rng) and rng[0] <= rng[1]):
            probs.append(f"expects.{key} must be [low, high] in micrometres, got {rng!r}")

    models = data.get("models") or []
    names = set()
    if not isinstance(models, list):
        probs.append("models is not a list")
        models = []
    for i, m in enumerate(models):
        if not isinstance(m, dict):
            probs.append(f"models[{i}] is not a mapping")
            continue
        name = str(m.get("name", ""))
        kind = str(m.get("kind", ""))
        sha = str(m.get("sha256", "")).lower()
        size = m.get("size")
        if name in names:
            probs.append(f"model {name!r} is declared twice")
        names.add(name)
        if not SHA256_RE.fullmatch(sha):
            probs.append(f"model {name!r}: sha256 missing or not 64 hex")
        if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
            probs.append(f"model {name!r}: size missing or not a positive integer")
        if kind == "release":
            reg = MODELS.get(name)
            if reg is None:
                probs.append(f"model {name!r}: kind release but not a release model {sorted(MODELS)}")
            elif reg["sha256"] != sha or reg["size"] != size:
                probs.append(f"model {name!r}: sha256 or size differs from the release record")
        else:
            mm = PACK_MODEL_RE.fullmatch(name)
            if not mm:
                probs.append(f"model {name!r} breaks the pack model pattern <name>@<first 8 hex of sha256> (release "
                             f"and Cellpose built-in names cannot be pack model names)")
            elif SHA256_RE.fullmatch(sha) and mm.group(2) != sha[:8]:
                probs.append(f"model {name!r}: the name suffix {mm.group(2)} is not the start of its sha256 {sha[:8]}")
        fname = str(m.get("file_name", "") or name.split("@")[0])
        if not fname or "/" in fname or "\\" in fname or fname in (".", ".."):
            probs.append(f"model {name!r}: file_name {fname!r} must be a bare file name")
        srcs = m.get("sources") or []
        if not isinstance(srcs, list) or any(not isinstance(s, dict) or s.get("type") not in SOURCE_TYPES
                                             for s in srcs):
            probs.append(f"model {name!r}: sources must be a list of {{type: {'|'.join(SOURCE_TYPES)}, ...}}")
        card = m.get("card")
        if card is not None and not _is_plain_relative(str(card)):
            probs.append(f"model {name!r}: card {card!r} must be a relative path inside the pack")

    recipes = data.get("recipes")
    if not isinstance(recipes, list) or not recipes:
        probs.append("recipes is missing or empty")
        recipes = []
    ids = set()
    n_default = 0
    for i, r in enumerate(recipes):
        if not isinstance(r, dict):
            probs.append(f"recipes[{i}] is not a mapping")
            continue
        rid = str(r.get("id", ""))
        if not RECIPE_ID_RE.fullmatch(rid):
            probs.append(f"recipe id {rid!r} does not match {RECIPE_ID_RE.pattern}")
        if rid in ids:
            probs.append(f"recipe {rid!r} is declared twice")
        ids.add(rid)
        prof = str(r.get("profile", ""))
        if not prof:
            probs.append(f"recipe {rid!r}: no profile")
        elif not _is_plain_relative(prof):
            probs.append(f"recipe {rid!r}: profile {prof!r} must be a relative path (inside the pack, or under the "
                         f"engine config root for root: checkout)")
        elif root == "pack" and not PACK_PROFILE_RE.fullmatch(prof):
            probs.append(f"recipe {rid!r}: profile {prof!r} must be config/profiles/<recipe>.yaml inside the pack "
                         f"(that layout anchors the recipe's relative paths to the pack folder)")
        if str(r.get("dims", "3d")) not in DIMS:
            probs.append(f"recipe {rid!r}: dims {r.get('dims')!r} is not one of {DIMS}")
        status = str(r.get("status", "provisional"))
        if status not in STATUSES:
            probs.append(f"recipe {rid!r}: status {status!r} is not one of {STATUSES}")
        if not HASH12_RE.fullmatch(str(r.get("expected_config_hash", ""))):
            probs.append(f"recipe {rid!r}: expected_config_hash must be 12 hex digits")
        if r.get("default"):
            n_default += 1
        needs = r.get("needs") or {}
        if not isinstance(needs, dict):
            probs.append(f"recipe {rid!r}: needs is not a mapping")
            needs = {}
        for key in ("segment_backends", "spot_backends", "models", "modules", "post_stages"):
            v = needs.get(key)
            if v is not None and not (isinstance(v, list) and all(isinstance(x, str) for x in v)):
                probs.append(f"recipe {rid!r}: needs.{key} must be a list of names")
        for mname in needs.get("models") or []:
            if mname not in names:
                probs.append(f"recipe {rid!r}: needs model {mname!r}, which the pack does not declare")
        if needs.get("accelerator") is not None and str(needs.get("accelerator")) not in ACCELERATOR:
            probs.append(f"recipe {rid!r}: needs.accelerator must be one of {ACCELERATOR}")
        tun = r.get("tunables") or []
        if not isinstance(tun, list) or any(not isinstance(t, dict) or not t.get("key") for t in tun):
            probs.append(f"recipe {rid!r}: tunables must be a list of mappings with a key")
        claims = status == "validated" or r.get("validated") is True
        if claims:
            if folder is None:
                probs.append(f"recipe {rid!r}: claims validated, which needs a passing validation record")
            else:
                ok, why = validation_record_passes(folder, r.get("validation"))
                if not ok:
                    probs.append(f"recipe {rid!r}: claims validated without a passing validation record ({why})")
    if n_default > 1:
        probs.append("more than one recipe is marked default")
    files = data.get("engine_files") or {}
    if root == "checkout":
        if not isinstance(files, dict) or any(not _is_plain_relative(k) or not SHA256_RE.fullmatch(str(v))
                                              for k, v in files.items()):
            probs.append("engine_files must map relative engine paths to LF-normalised sha256 digests")
    elif files:
        probs.append("engine_files is only for root: checkout packs")
    return probs


def parse_pack(folder: Path, data: dict, source: str = "path") -> Pack:
    """A Pack from a parsed manifest; ``problems`` holds what validate_manifest found (an unusable pack is
    still listed, never raised)."""
    folder = Path(folder)
    problems = validate_manifest(data, folder)
    recipes = []
    for r in data.get("recipes") or []:
        if not isinstance(r, dict):
            continue
        recipes.append(Recipe(id=str(r.get("id", "")), profile=str(r.get("profile", "")),
                              title=str(r.get("title", "")), dims=str(r.get("dims", "3d")),
                              default=bool(r.get("default", False)), status=str(r.get("status", "provisional")),
                              expected_config_hash=str(r.get("expected_config_hash", "")),
                              needs=dict(r.get("needs") or {}) if isinstance(r.get("needs"), dict) else {},
                              tunables=list(r.get("tunables") or []) if isinstance(r.get("tunables"), list) else [],
                              validation=r.get("validation"), raw=dict(r)))
    models = []
    for m in data.get("models") or []:
        if not isinstance(m, dict):
            continue
        name = str(m.get("name", ""))
        size = m.get("size")
        models.append(ModelSpec(name=name, kind=str(m.get("kind", "")),
                                file_name=str(m.get("file_name", "") or name.split("@")[0]),
                                sha256=str(m.get("sha256", "")).lower(),
                                size=int(size) if isinstance(size, int) and not isinstance(size, bool) else 0,
                                sources=list(m.get("sources") or []) if isinstance(m.get("sources"), list) else [],
                                card=m.get("card"), licence=str(m.get("licence", ""))))
    return Pack(id=str(data.get("id", "")), version=str(data.get("version", "")), folder=folder, data=data,
                source=source, recipes=recipes, models=models, problems=problems)


def load_pack_folder(folder: Path, source: str = "path") -> Pack:
    """Read and validate the pack in ``folder``. PackFormatError only when pack.yaml cannot be parsed."""
    return parse_pack(Path(folder), load_manifest(folder), source)
