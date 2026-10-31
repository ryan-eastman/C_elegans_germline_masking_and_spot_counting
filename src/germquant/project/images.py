"""Adding images to a project (docs/PROJECTS.md, "Images"; design 7.8; PLAN.md decision D27; P5 acceptance 7).

add_images takes files, folders and glob patterns. A folder is searched (recursively) for the recipe's
io.input_suffixes when the recipe sets them, and for every suffix germquant reads (io.formats.SUPPORTED_SUFFIXES)
otherwise, so a lab can add a folder of OME-TIFFs to a project on a pack whose base names only .nd2 files. A
recipe whose io.input_suffixes recipe_keys refuses (as every run of it does) is refused before anything is
added. The search always skips germquant's own outputs (files named <id>__nuclei_labels.tif and the other
OUTPUT_SUFFIXES) and the project's own results/, runs/, training/, tuning/, annotations/ and packs/ folders. An
OME-Zarr folder is one image and is not searched.

Every candidate's header is read with project.headers (light for .nd2 and TIFF files). Refused, each with its
reason, while the other files are still added: a missing file, an unsupported suffix, a germquant output, a file
inside the project's own folders, an unreadable header, a multi-scene file (more than one OME image or TIFF
series), an .nd2 with more than one position (P > 1), and a duplicate image id. The id is the file name without
its format suffix (io.formats.image_stem), and it names the results folder, so two ids may not differ only in
letter case either (Windows and macOS would put both into one folder).

A path inside the project is stored relative, with forward slashes, so the project can move; any other path is
stored absolute, a UNC path verbatim.
"""
from __future__ import annotations

import glob
import logging
import os
from pathlib import Path, PurePosixPath, PureWindowsPath

from ..fsutil import long_path
from ..io.formats import SUPPORTED_SUFFIXES, image_format, image_stem
from ..recipe_keys import input_suffixes_problem
from .manifest import PROJECT_DIRS

log = logging.getLogger(__name__)

# the image files germquant itself writes (pipeline.py, validate/imaris_ims.py): never an input
OUTPUT_SUFFIXES = ("__nuclei_labels.tif", "__spots.tif", "__granules_labels.tif", "__envelope_labels.tif",
                   "__sc_ribbon.tif", "__syp_aggregate.tif", "__imaris_labels.tif")
_DIR_SUFFIXES = (".zarr", ".ome.zarr")
_GLOB_CHARS = set("*?[")


def is_output_name(name: str) -> str | None:
    """The output suffix ``name`` ends with (a .tif or .tiff germquant writes), or None."""
    n = str(name).lower()
    return next((s for s in OUTPUT_SUFFIXES if n.endswith(s) or n.endswith(s + "f")), None)


def stored_path(root, path) -> str:
    """How the manifest records ``path``: relative with forward slashes inside ``root``, else absolute."""
    p = Path(os.path.abspath(path))
    try:
        return p.relative_to(Path(os.path.abspath(root))).as_posix()
    except ValueError:
        return str(p)


def resolve_path(root, stored: str) -> Path:
    """The file a manifest path names: an absolute path (of any platform, UNC included) as it is, a relative one
    under the project root."""
    s = str(stored)
    if PureWindowsPath(s).is_absolute() or PurePosixPath(s).is_absolute():
        return Path(s)
    return Path(root) / PurePosixPath(s)


def _matches(name: str, suffixes) -> bool:
    n = name.lower()
    return any(n.endswith(s) and len(n) > len(s) for s in suffixes)


def _project_dir_of(root, path) -> str | None:
    """The project folder (results, runs, ...) that ``path`` lies in, or None."""
    p = os.path.normcase(os.path.abspath(path))
    for d in PROJECT_DIRS:
        top = os.path.normcase(os.path.abspath(os.path.join(root, d)))
        if p == top or p.startswith(top + os.sep):
            return d
    return None


def _plain(s: str) -> str:
    """A path string without the Windows extended-length prefix that fsutil.long_path adds."""
    unc, local = "\\\\?\\UNC\\", "\\\\?\\"          # the two prefixes of fsutil.long_path
    if s.startswith(unc):
        return "\\\\" + s[len(unc):]
    if s.startswith(local):
        return s[len(local):]
    return s


def _exists(path) -> bool:
    return os.path.exists(long_path(path))


def path_too_long(path, opener: str = "the image readers") -> str | None:
    """Why ``path`` cannot be used, on Windows without long paths, when it exists only in its extended-length form
    (fsutil.long_path): past the 260 characters Windows allows, which the image readers and the model loader do
    not lift. None elsewhere, and whenever the plain path works."""
    if os.name != "nt" or os.path.exists(path) or not _exists(path):
        return None
    return (f"the path is {len(os.path.abspath(path))} characters long, past the 260 characters Windows allows, so "
            f"{opener} cannot open it: move the folder to a shorter path (or map a drive letter to it), or turn on "
            f"Windows long paths (LongPathsEnabled)")


def discover_folder(folder, suffixes, root) -> list[Path]:
    """The image files under ``folder`` (module docstring), sorted by path. The walk goes through
    fsutil.long_path, so a deep NAS tree past Windows' 260-character limit is searched too."""
    skip = {os.path.normcase(os.path.abspath(os.path.join(root, d))) for d in PROJECT_DIRS}
    out: list[Path] = []
    for long_dir, dirnames, filenames in os.walk(long_path(os.path.abspath(folder))):
        dirpath = _plain(long_dir)
        keep = []
        for d in sorted(dirnames):
            full = os.path.join(dirpath, d)
            if os.path.normcase(os.path.abspath(full)) in skip:
                continue
            if _matches(d, [s for s in suffixes if s in _DIR_SUFFIXES]):
                out.append(Path(full))                       # an OME-Zarr image: one image, not searched
                continue
            keep.append(d)
        dirnames[:] = keep
        for f in filenames:
            if _matches(f, suffixes) and not is_output_name(f):
                out.append(Path(dirpath) / f)
    return sorted(out, key=lambda p: str(p))


def input_suffixes(cfg) -> tuple[str, ...]:
    """The recipe's io.input_suffixes, or every suffix germquant reads when it sets none. A value that
    recipe_keys.input_suffixes_problem refuses (every run of the recipe and `germquant batch` refuse it too: the
    rule is a non-empty list of the suffixes germquant reads, in lower case with the dot) is a ProjectError with
    that message, so a folder is never searched with suffixes the recipe cannot run."""
    raw = cfg.get("io.input_suffixes")
    if raw is None:
        return tuple(SUPPORTED_SUFFIXES)
    problem = input_suffixes_problem(raw)
    if problem:
        from .manifest import ProjectError

        raise ProjectError(problem)
    return tuple(raw)


def candidates(root, paths, suffixes) -> list[tuple[str, Path | None, str | None]]:
    """``[(as given, path, refusal)]`` in the order of ``paths``: the files each entry names, expanded (a folder is
    searched, a glob pattern matched); a pattern that matches nothing is one refusal (path None)."""
    found: list[tuple[str, Path | None, str | None]] = []
    for item in paths:
        s = os.fspath(item)
        p = Path(s)
        if _GLOB_CHARS & set(s) and not p.exists():
            hits = sorted(glob.glob(s, recursive=True))
            if not hits:
                found.append((s, None, "no file matches the pattern"))
            for h in hits:
                hp = Path(h)
                if hp.is_dir() and not _matches(hp.name, _DIR_SUFFIXES):
                    found += [(str(f), f, None) for f in discover_folder(hp, suffixes, root)]
                else:
                    found.append((h, hp, None))
            continue
        if os.path.isdir(long_path(p)) and not _matches(p.name, _DIR_SUFFIXES):
            found += [(str(f), f, None) for f in discover_folder(p, suffixes, root)]
        else:
            found.append((s, p, None))
    return found


def check_candidate(root, path: Path):
    """``(header, None)`` for a file that may be added, or ``(None, reason)``."""
    from .headers import read_header

    if not _exists(path):
        return None, f"not found: {path}"
    too_long = path_too_long(path)
    if too_long is not None:
        return None, too_long
    where = _project_dir_of(root, path)
    if where is not None:
        return None, f"inside the project's {where}/ folder, which germquant writes or keeps for itself"
    out = is_output_name(path.name)
    if out is not None:
        return None, f"a germquant output file ({out}), not an input image"
    try:
        image_format(path)
    except ValueError as e:
        return None, f"unsupported: {e}"
    try:
        header = read_header(path)
    except Exception as e:  # noqa: BLE001 - every reader failure is a refusal with its message
        return None, f"header unreadable: {type(e).__name__}: {e}"
    n = int(header.get("n_images", 1) or 1)
    if header.get("reader") == "nd2" and n > 1:
        return None, (f"an .nd2 with {n} positions (P = {n}); germquant reads one position per file: export each "
                      f"position as its own .nd2")
    if n > 1:
        return None, (f"a multi-scene file: {n} images in one file; germquant reads one image per file: export or "
                      f"split it into one file per scene")
    return header, None


def image_entry(root, path: Path, header: dict) -> dict:
    """The manifest record of an added image: ``{id, path, size, mtime[, timepoints], excluded, reason}``."""
    st = os.stat(long_path(path))
    entry = {"id": image_stem(path), "path": stored_path(root, path),
             "size": None if os.path.isdir(long_path(path)) else int(st.st_size), "mtime": float(st.st_mtime)}
    n_t = int((header.get("sizes") or {}).get("T", 1) or 1)
    if n_t > 1:
        entry["timepoints"] = n_t
    entry.update({"excluded": False, "reason": None})
    return entry


_DEVICE_NAMES = frozenset({"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)),
                           *(f"LPT{i}" for i in range(1, 10))})


def windows_device_name(name: str) -> bool:
    """Whether ``name`` (before any dot) is a name Windows reserves for a device, such as CON or COM1."""
    return str(name).split(".")[0].strip().upper() in _DEVICE_NAMES


def id_reason(image_id: str) -> str | None:
    """Why ``image_id`` cannot name a results folder on every system, or None."""
    if not image_id.strip():
        return "the file name leaves an empty image id"
    if image_id != image_id.rstrip(". "):
        return (f"the image id {image_id!r} ends with a dot or a space, which Windows drops from folder names, so its "
                f"results folder would not carry its id: rename the file")
    if windows_device_name(image_id):
        return f"the image id {image_id!r} is a Windows device name, which no results folder may have: rename the file"
    return None


def duplicate_reason(image_id: str, entries: list[dict]) -> str | None:
    """Why ``image_id`` cannot join ``entries`` (an id that cannot name a folder, an equal id, or one that differs
    only in letter case), or None."""
    bad = id_reason(image_id)
    if bad is not None:
        return bad
    for e in entries:
        if e["id"] == image_id:
            return f"duplicate image id {image_id}: the project already holds {e['path']}"
    for e in entries:
        if e["id"].lower() == image_id.lower():
            return (f"duplicate image id {image_id}: the project holds {e['id']} ({e['path']}), and ids that differ "
                    f"only in letter case would share one results folder on Windows and macOS")
    return None
