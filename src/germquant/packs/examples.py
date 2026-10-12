"""Pack examples: the synthetic example inside a pack, and real public images fetched on demand
(`germquant pack fetch-example`; docs/PACKS.md; docs/DESIGN_system_packs.md, sections 3.3 and 6.9, package P3).

A pack ships a tiny synthetic example (``example/<name>.ome.tif``, under 2 MB, ``{id, image}``) and may name
real public images, which are never redistributed. A real example lists its source files at a public
address, each pinned by size and sha256, and how to assemble them into one OME-TIFF (germquant.assemble:
declared axes, channel names, voxel size, crop), pinned by pixels_sha256:

    examples:
      - id: sbiad3750_sut509_p1
        recipe: cells_spots                         # the recipe it is meant for (default: the pack's default)
        base_url: "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/750/S-BIAD3750/Files/"
        sources:
          - {path: "SpotMAX_gt_dataset/.../<stem>_Bright.tif", size: 14701884, sha256: "<64 hex>"}
          - {path: "SpotMAX_gt_dataset/.../<stem>_Cy3.tif", size: 14701884, sha256: "<64 hex>"}
        assemble: {axes: ZYX, channels: [Bright, Cy3], spacing_zyx_um: [0.24, 0.0720635, 0.0720635]}
        pixels_sha256: "<64 hex>"
        licence: CC0
        licence_spdx: CC0-1.0                       # pack lint holds it to PERMISSIVE_DATA_LICENCES
        licence_url: "https://creativecommons.org/publicdomain/zero/1.0/"
        attribution: "Padovani et al. 2026, Science Advances; BioImage Archive S-BIAD3750"
        expect: {n_objects: [lo, hi], spots_per_object: [lo, hi]}

A source inside a zip archive names the archive as its ``path`` and adds ``member`` (and optionally
``zip_size``, the archive's size, and ``crc32``); only that member is read, with HTTP Range requests
(datafetch.fetch_zip_member). Digests and the CRC-32 are quoted strings in pack.yaml (an unquoted value of
digits only would be read as a number). Each source is saved under its file name (the member's, for a zip
member), which must be a plain file name, unique within the example. ``expect`` gives the ranges a run of the
example's recipe should land in: the image_summary's n_nuclei (objects, whatever the pack's object_kind) and
mean_spots (spots per object).

fetch_example writes ``<out>/<pack id>/<example id>/source/<file>`` (each checked by size and sha256; a file
already there with the right digest is kept), assembles ``<out>/<pack id>/<example id>/<example id>.ome.tif``,
checks its pixels_sha256 and writes ``<example id>.json`` beside it: where every byte came from, the licence,
what was changed and the attribution. A second call finds the image with its pinned pixels, channel names and
voxel size in place and downloads nothing (the source files may then be deleted to save space); an image whose
pixels match but whose channel names or voxel size no longer match the pack (a later pack fix) is written again
from its pixels, and a record that no longer matches is rewritten. ``out`` defaults to ``$GERMQUANT_EXAMPLES``,
else the per-user folder (``%LOCALAPPDATA%/germquant/examples``, ``~/Library/Application Support/germquant/
examples`` or ``~/.local/share/germquant/examples``); a folder inside a pack root or a pack is refused, because
public images never go into a pack. A pack that is unusable (its manifest has problems) fetches nothing.

The address of a real example is https. Plain http to this machine (127.0.0.1, localhost, [::1]) is accepted
only while ``$GERMQUANT_ALLOW_LOOPBACK_EXAMPLES`` is 1, for the stand-in servers of tests; `pack lint` refuses it
in any case.

Torch-free: datafetch (standard library) and germquant.assemble (numpy, tifffile) only, so `python -m
germquant.packs fetch-example` never loads torch, Cellpose, SpotMAX or CuPy (gate G6).
"""
from __future__ import annotations

import csv
import json
import math
import os
import re
import sys
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

EXAMPLES_ENV = "GERMQUANT_EXAMPLES"
ALLOW_LOOPBACK_ENV = "GERMQUANT_ALLOW_LOOPBACK_EXAMPLES"
EXAMPLE_ID_RE = re.compile(r"^[a-z0-9][a-z0-9_]{0,63}$")
# a saved source's file name: plain characters, no leading dot or dash, no Windows device name, never a name
# the fetcher uses for its own partial files
SOURCE_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._()+-]{0,200}$")
_DEVICE_NAME_RE = re.compile(r"^(con|prn|aux|nul|com[0-9]|lpt[0-9])(\..*)?$", re.IGNORECASE)
_SPDX_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9.+-]{0,63}$")
EXPECT_KEYS = ("n_objects", "spots_per_object")


class ExampleError(RuntimeError):
    """An example that cannot be fetched or does not assemble to its pinned pixels."""


def examples_dir() -> Path:
    """Where fetched examples go by default: ``$GERMQUANT_EXAMPLES``, else the per-user folder."""
    env = os.environ.get(EXAMPLES_ENV, "").strip()
    if env:
        return Path(env)
    if os.name == "nt":
        base = os.environ.get("LOCALAPPDATA") or str(Path.home() / "AppData" / "Local")
        return Path(base) / "germquant" / "examples"
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support" / "germquant" / "examples"
    base = os.environ.get("XDG_DATA_HOME") or str(Path.home() / ".local" / "share")
    return Path(base) / "germquant" / "examples"


# ---- the manifest's examples -----------------------------------------------------------------------------

def _plain_relative(path: str) -> bool:
    pp = PurePosixPath(str(path))
    return (bool(str(path)) and not pp.is_absolute() and ".." not in pp.parts and "\\" not in str(path)
            and not re.match(r"^[A-Za-z]:", str(path)))


def _range_ok(rng) -> bool:
    return (isinstance(rng, list) and len(rng) == 2
            and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in rng) and rng[0] <= rng[1])


def _hex_problem(value, n_digits: int, what: str) -> str | None:
    """A digest or CRC-32 written in pack.yaml: a quoted string of ``n_digits`` lowercase hex digits. An
    unquoted value of digits only is read by YAML as a number, and a CRC-32 such as 78066861 would then be taken
    as decimal: refused, with the fix."""
    if isinstance(value, str) and re.fullmatch(f"[0-9a-f]{{{n_digits}}}", value):
        return None
    if isinstance(value, int) and not isinstance(value, bool):
        return (f"{what} must be {n_digits} lowercase hex digits in quotes; pack.yaml gives the number {value} "
                f"(YAML reads an unquoted value of digits only as a number): quote it")
    return f"{what} must be {n_digits} lowercase hex digits"


def is_remote(spec: dict) -> bool:
    return isinstance(spec, dict) and bool(spec.get("sources"))


def example_problems(examples, recipe_ids=None) -> list[str]:
    """Problems with a pack.yaml ``examples`` list (module docstring): an id that breaks its pattern or repeats;
    an example that is neither a local ``image`` (a relative path inside the pack) nor remote ``sources``; a
    remote example without an https base_url (plain http only to this machine, and only while
    $GERMQUANT_ALLOW_LOOPBACK_EXAMPLES is 1), sized and digested sources saved under plain, unique file names, an
    assemble block (axes, channel names, three positive voxel sizes, an optional crop), pixels_sha256, a
    licence or an attribution; a licence_spdx or licence_url that is malformed; an ``expect`` range that is not
    [low, high]; a ``recipe`` the pack does not declare."""
    probs: list[str] = []
    if examples is None:
        return probs
    if not isinstance(examples, list):
        return ["examples is not a list"]
    seen = set()
    for i, ex in enumerate(examples):
        if not isinstance(ex, dict):
            probs.append(f"examples[{i}] is not a mapping")
            continue
        eid = str(ex.get("id", ""))
        where = f"example {eid!r}" if eid else f"examples[{i}]"
        if not EXAMPLE_ID_RE.fullmatch(eid):
            probs.append(f"{where}: id does not match {EXAMPLE_ID_RE.pattern}")
        if eid in seen:
            probs.append(f"{where} is declared twice")
        seen.add(eid)
        if recipe_ids is not None and ex.get("recipe") is not None and str(ex.get("recipe")) not in recipe_ids:
            probs.append(f"{where}: recipe {ex.get('recipe')!r} is not a recipe of the pack")
        has_image, has_sources = ex.get("image") is not None, ex.get("sources") is not None
        if has_image == has_sources:
            probs.append(f"{where}: give either image (a file inside the pack) or sources (files to fetch)")
            continue
        if has_image:
            if not _plain_relative(str(ex["image"])):
                probs.append(f"{where}: image {ex['image']!r} must be a relative path inside the pack")
            continue
        probs += [f"{where}: {p}" for p in _remote_problems(ex)]
    return probs


def is_loopback_url(url: str) -> bool:
    """Whether ``url`` is plain http to this machine (127.0.0.1, localhost or [::1]): the stand-in servers of
    tests, which a loader accepts only while $GERMQUANT_ALLOW_LOOPBACK_EXAMPLES is 1 and `pack lint` refuses."""
    host = urllib.parse.urlsplit(str(url)).hostname or ""
    return str(url).startswith("http://") and host in ("127.0.0.1", "localhost", "::1")


def loopback_allowed() -> bool:
    """Whether plain http to this machine is accepted as an example address (tests set the variable)."""
    return os.environ.get(ALLOW_LOOPBACK_ENV, "").strip() == "1"


def source_name_problem(name: str) -> str | None:
    """Why a source's saved file name is refused, or None (module docstring)."""
    if not SOURCE_NAME_RE.fullmatch(name):
        return (f"saves as {name!r}, which is not a plain file name (letters, digits and . _ ( ) + -, starting "
                f"with a letter or digit)")
    if _DEVICE_NAME_RE.fullmatch(name):
        return f"saves as {name!r}, a Windows device name"
    if name.lower().endswith((".part", ".zpart")):
        return f"saves as {name!r}, which ends like the fetcher's partial files"
    return None


def _remote_problems(ex: dict) -> list[str]:
    out: list[str] = []
    base = str(ex.get("base_url", "") or "")
    if not (base.startswith("https://") or (is_loopback_url(base) and loopback_allowed())):
        hint = (f" (plain http to this machine is for test servers, with {ALLOW_LOOPBACK_ENV}=1)"
                if is_loopback_url(base) else "")
        out.append(f"base_url {base!r} must be an https:// address{hint}")
    srcs = ex.get("sources")
    if not isinstance(srcs, list) or not srcs:
        return out + ["sources must be a non-empty list"]
    names: list[str] = []
    for j, s in enumerate(srcs):
        if not isinstance(s, dict):
            out.append(f"sources[{j}] is not a mapping")
            continue
        path_ok = _plain_relative(str(s.get("path", "")))
        if not path_ok:
            out.append(f"sources[{j}]: path {s.get('path')!r} must be a relative path under base_url")
        size = s.get("size")
        if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
            out.append(f"sources[{j}]: size must be the file's size in bytes")
        p = _hex_problem(s.get("sha256"), 64, "sha256")
        if p:
            out.append(f"sources[{j}]: {p}")
        member_ok = True
        if s.get("member") is not None and not _plain_relative(str(s["member"])):
            member_ok = False
            out.append(f"sources[{j}]: member {s['member']!r} must be a path inside the archive")
        if s.get("crc32") is not None:
            p = _hex_problem(s["crc32"], 8, "crc32")
            if p:
                out.append(f"sources[{j}]: {p}")
        zs = s.get("zip_size")
        if zs is not None and (isinstance(zs, bool) or not isinstance(zs, int) or zs <= 0):
            out.append(f"sources[{j}]: zip_size must be the archive's size in bytes")
        if (s.get("crc32") is not None or zs is not None) and s.get("member") is None:
            out.append(f"sources[{j}]: crc32 and zip_size belong to a zip member (member)")
        if path_ok and member_ok:
            name = source_name(s)
            p = source_name_problem(name)
            if p:
                out.append(f"sources[{j}] {p}")
            names.append(name)
    folded = [n.lower() for n in names]
    dups = sorted({n for n in names if folded.count(n.lower()) > 1})
    if dups:
        out.append(f"sources save under the same file name {', '.join(dups)} (names are compared ignoring case, as "
                   f"Windows and macOS do); each source of an example needs its own file name")
    asm = ex.get("assemble")
    if not isinstance(asm, dict):
        out.append("assemble must be a mapping with axes, channels and spacing_zyx_um")
    else:
        axes = str(asm.get("axes", ""))
        if not re.fullmatch(r"[CZYX]+", axes) or len(set(axes)) != len(axes) or "Y" not in axes or "X" not in axes:
            out.append(f"assemble.axes {axes!r} must use the letters C, Z, Y and X once each at most, with Y and X")
        chans = asm.get("channels")
        if (not isinstance(chans, list) or not chans or not all(isinstance(c, str) and c.strip() for c in chans)
                or len(set(chans)) != len(chans)):
            out.append("assemble.channels must list unique channel names")
        sp = asm.get("spacing_zyx_um")
        if not (isinstance(sp, list) and len(sp) == 3
                and all(isinstance(v, (int, float)) and not isinstance(v, bool) and v > 0 for v in sp)):
            out.append("assemble.spacing_zyx_um must be three positive numbers (dz, dy, dx in micrometres)")
        crop = asm.get("crop")
        if crop is not None:
            if not isinstance(crop, dict) or not crop:
                out.append("assemble.crop must map z, y or x to [start, stop)")
            else:
                for k, rng in crop.items():
                    if str(k) not in ("z", "y", "x") or not (
                            isinstance(rng, list) and len(rng) == 2
                            and all(isinstance(v, int) and not isinstance(v, bool) for v in rng)
                            and 0 <= rng[0] < rng[1]):
                        out.append(f"assemble.crop {k}: give [start, stop) as whole numbers, start < stop")
    p = _hex_problem(ex.get("pixels_sha256"), 64, "pixels_sha256")
    if p:
        out.append(f"{p} (germquant.assemble.pixels_sha256)")
    if not str(ex.get("licence", "") or "").strip():
        out.append("licence is missing (the data licence, checked at build time)")
    spdx = ex.get("licence_spdx")
    if spdx is not None and not (isinstance(spdx, str) and _SPDX_RE.fullmatch(spdx)):
        out.append(f"licence_spdx {spdx!r} must be an SPDX licence id such as CC0-1.0 or CC-BY-4.0")
    lurl = ex.get("licence_url")
    if lurl is not None and not (isinstance(lurl, str) and lurl.startswith("https://") and len(lurl) > len("https://")):
        out.append(f"licence_url {lurl!r} must be an https:// address (the licence text)")
    if not str(ex.get("attribution", "") or "").strip():
        out.append("attribution is missing (who made the data, and where it is deposited)")
    exp = ex.get("expect")
    if exp is not None:
        if not isinstance(exp, dict):
            out.append("expect must be a mapping")
        else:
            for k, rng in exp.items():
                if k not in EXPECT_KEYS:
                    out.append(f"expect.{k}: unknown (known: {', '.join(EXPECT_KEYS)})")
                elif not _range_ok(rng):
                    out.append(f"expect.{k} must be [low, high], got {rng!r}")
    return out


def list_examples(pack) -> list[dict]:
    """The examples a pack declares, each with ``kind`` (local or remote) and the recipe it is meant for."""
    out = []
    for ex in pack.data.get("examples") or []:
        if not isinstance(ex, dict):
            continue
        rec = dict(ex)
        rec["kind"] = "remote" if is_remote(ex) else "local"
        rec["recipe"] = str(ex.get("recipe") or pack.default_recipe().id)
        out.append(rec)
    return out


def find_example(pack, example: str | None = None) -> dict:
    """The example ``example`` names, else the pack's first remote example, else its first example.
    FileNotFoundError naming the examples."""
    exs = list_examples(pack)
    if example:
        for ex in exs:
            if ex.get("id") == example:
                return ex
        raise FileNotFoundError(f"pack {pack.ref} has no example {example!r}; its examples: "
                                f"{', '.join(str(e.get('id')) for e in exs) or 'none'}")
    remote = [e for e in exs if e["kind"] == "remote"]
    if remote:
        return remote[0]
    if exs:
        return exs[0]
    raise FileNotFoundError(f"pack {pack.ref} declares no example")


def source_url(spec: dict, src: dict) -> str:
    base = str(spec["base_url"])
    if not base.endswith("/"):
        base += "/"
    return base + urllib.parse.quote(str(src["path"]), safe="/")


def source_name(src: dict) -> str:
    """The file name a source is saved under: the member's (for a zip member) or the path's last part."""
    return PurePosixPath(str(src.get("member") or src["path"])).name


# ---- fetching --------------------------------------------------------------------------------------------

def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _verified(path: Path, size: int, sha256: str) -> bool:
    from ..datafetch import sha256_file

    return path.is_file() and path.stat().st_size == int(size) and sha256_file(path) == str(sha256).lower()


def _fetch_source(spec: dict, src: dict, dest: Path, *, download: bool, say) -> dict:
    from .. import datafetch as D

    url = source_url(spec, src)
    if _verified(dest, src["size"], src["sha256"]):
        say(f"  {dest.name}: already here, size and sha256 checked")
        return _source_record(spec, src, dest, "present")
    if dest.exists():
        say(f"  {dest.name}: present but not the pinned file; fetching it again")
        dest.unlink()
    if not download:
        raise ExampleError(f"{dest} is missing and downloads are switched off (--no-download)")
    progress = _progress(say, dest.name)
    if src.get("member"):
        say(f"  {dest.name}: reading member {src['member']} ({int(src['size']):,} bytes) of {url} with range requests")
        D.fetch_zip_member(url, str(src["member"]), dest, size=int(src["size"]), sha256=str(src["sha256"]),
                           crc32=src.get("crc32"), zip_size=src.get("zip_size"), progress=progress)
    else:
        say(f"  {dest.name}: downloading {int(src['size']):,} bytes from {url}")
        D.download(url, dest, size=int(src["size"]), sha256=str(src["sha256"]), progress=progress)
    return _source_record(spec, src, dest, "downloaded")


def _progress(say, name: str):
    import time

    state = {"t": time.time()}

    def report(have: int, total: int | None) -> None:
        now = time.time()
        if now - state["t"] >= 30 and total:
            state["t"] = now
            say(f"  {name}: {have / 1e6:.1f} of {total / 1e6:.1f} MB")
    return report


def _inside(path: Path, folder: Path) -> bool:
    return path == folder or folder in path.parents


def _examples_root(out, pack) -> Path:
    """The examples folder, resolved; ExampleError when it lies inside a pack root, the engine's built-in or
    test pack folders, or the pack itself (fetched public images never go into a pack: git tracks every file
    under packs/, and pack verify and build-sums would list them)."""
    from ..config import builtin_roots
    from .discovery import pack_roots

    root = (Path(out) if out is not None else examples_dir()).expanduser().resolve()
    guarded = [(kind, Path(folder)) for kind, folder in pack_roots()]
    guarded += [("built-in test packs", Path(r) / "tests" / "fixtures" / "packs") for r in builtin_roots()]
    guarded.append((f"pack {pack.ref}", Path(pack.folder)))
    for kind, folder in guarded:
        f = folder.expanduser().resolve()
        if _inside(root, f):
            raise ExampleError(f"the examples folder {root} lies inside {f} ({kind}); fetched public images never go "
                               f"into a pack folder (git tracks the files there, and pack verify would list them), "
                               f"so choose another folder with --out, such as the default {examples_dir()}")
    return root


def _metadata_differences(got: dict, asm: dict) -> list[str]:
    """What an assembled image's channel names and voxel size get wrong against the example's assemble block."""
    out = []
    if [str(n) for n in got.get("channel_names") or []] != [str(c) for c in asm["channels"]]:
        out.append(f"channel names {got.get('channel_names')} (the pack says {list(asm['channels'])})")
    sp, want = got.get("spacing_zyx_um"), [float(v) for v in asm["spacing_zyx_um"]]
    if (sp is None or len(sp) != 3
            or not all(math.isclose(float(a), b, rel_tol=1e-9, abs_tol=0.0) for a, b in zip(sp, want))):
        out.append(f"voxel size {sp} um (the pack says {want})")
    return out


def fetch_example(pack, example: str | None = None, out: str | Path | None = None, *, download: bool = True,
                  say=None) -> dict:
    """Fetch and assemble one example of a pack (module docstring). Returns its record: ``id``, ``kind``,
    ``image`` (the OME-TIFF to run), ``recipe``, ``status`` (local; fetched: assembled now, or its channel
    names and voxel size written again from its pinned pixels; present: in place and current) and, for a
    remote example, ``pixels_sha256``, ``sources``, ``licence``, ``licence_url``, ``attribution``, ``expect`` and
    ``record`` (the JSON beside the image). ExampleError when the pack is unusable, when the examples folder lies
    inside a pack, or when a file or the assembled pixels do not match their pins (nothing wrong is kept)."""
    from ..assemble import assemble, read_ome_tiff, write_ome_tiff
    from . import load_pack

    say = say or (lambda *a: None)
    p = load_pack(pack)
    if p.problems:
        raise ExampleError(f"pack {p.ref} is unusable, so no example is fetched: " + "; ".join(p.problems))
    spec = find_example(p, example)
    base = {"pack": p.id, "version": p.version, "id": spec["id"], "kind": spec["kind"], "recipe": spec["recipe"],
            "expect": spec.get("expect")}
    if spec["kind"] == "local":
        image = p.folder / str(spec["image"])
        if not image.is_file():
            raise ExampleError(f"pack {p.ref}: the example image {spec['image']} is missing")
        return {**base, "image": str(image), "status": "local"}
    root = _examples_root(out, p)
    folder = (root / p.id / str(spec["id"])).resolve()
    if folder.parent.parent != root:
        raise ExampleError(f"pack {p.ref}: example {spec['id']!r} would be written to {folder}, outside {root}")
    image = folder / f"{spec['id']}.ome.tif"
    record_path = folder / f"{spec['id']}.json"
    want = str(spec["pixels_sha256"])
    asm = spec["assemble"]
    rec = {**base, "image": str(image), "record": str(record_path), "pixels_sha256": want,
           "licence": spec.get("licence"), "licence_url": spec.get("licence_url"),
           "attribution": spec.get("attribution")}
    srcs = spec["sources"]
    paths = [folder / "source" / source_name(s) for s in srcs]
    if (any(d.parent != folder / "source" or source_name_problem(d.name) for d in paths)
            or len({d.name.lower() for d in paths}) != len(paths)):                 # the loader refuses these too
        raise ExampleError(f"pack {p.ref}: example {spec['id']!r} names sources that do not save as distinct plain "
                           f"files in {folder / 'source'}")
    if image.is_file():
        try:
            got = read_ome_tiff(image)
        except Exception:  # noqa: BLE001 - an unreadable leftover is assembled again
            got = None
        if got is not None and got["pixels_sha256"] == want:     # the pixels are pinned: the sources may be gone
            stale = _metadata_differences(got, asm)
            status = "present"
            if stale:                                            # a pack fix of names or voxel size, same pixels
                say(f"{p.id}/{spec['id']}: {image.name} holds the pinned pixels, but its " + "; ".join(stale)
                    + "; writing it again from its pixels (nothing is downloaded)")
                write_ome_tiff(image, got["data"], asm["channels"], asm["spacing_zyx_um"])
                got = read_ome_tiff(image)
                if got["pixels_sha256"] != want or _metadata_differences(got, asm):
                    raise ExampleError(f"{p.id}/{spec['id']}: {image} did not keep its pixels or take the pack's "
                                       f"channel names and voxel size when written again")
                status = "fetched"
            else:
                say(f"{p.id}/{spec['id']}: already assembled at {image} (pixels_sha256, channel names and voxel "
                    f"size checked)")
            info = {"shape": list(got["data"].shape), "dtype": str(got["data"].dtype),
                    "channel_names": got["channel_names"], "spacing_zyx_um": got["spacing_zyx_um"]}
            sources = _current_sources(record_path, spec, srcs, paths)
            if stale or _record_stale(record_path, p, spec, image, sources, info):
                _write_record(record_path, p, spec, image, sources, info)
            return {**rec, "status": status, "sources": sources}
    say(f"{p.id}/{spec['id']}: {len(srcs)} source file(s); licence {spec.get('licence')}; {spec.get('attribution')}")
    fetched = [_fetch_source(spec, s, d, download=download, say=say) for s, d in zip(srcs, paths)]
    say(f"  assembling {image.name}: axes {asm['axes']}, channels {asm['channels']}, voxel "
        f"{asm['spacing_zyx_um']} um" + (f", crop {asm['crop']}" if asm.get("crop") else ""))
    res = assemble(paths, asm["channels"], asm["spacing_zyx_um"], image, axes=str(asm["axes"]), crop=asm.get("crop"))
    if res["pixels_sha256"] != want:
        image.unlink(missing_ok=True)
        record_path.unlink(missing_ok=True)          # a record of an earlier assembly would describe other pixels
        raise ExampleError(f"{p.id}/{spec['id']}: the assembled pixels have sha256 {res['pixels_sha256'][:12]}..., the "
                           f"pack pins {want[:12]}...; the image was discarded (the source files matched their "
                           f"digests, so the assemble settings or the reader differ from the pack's)")
    _write_record(record_path, p, spec, image, fetched, res)
    say(f"  -> {image}")
    return {**rec, "status": "fetched", "sources": fetched}


def _source_record(spec: dict, src: dict, dest: Path, status: str) -> dict:
    out = {"url": source_url(spec, src), "file": dest.name, "size": int(src["size"]), "sha256": str(src["sha256"])}
    if src.get("member"):
        out["member"] = str(src["member"])
    return {**out, "status": status}


_SOURCE_KEYS = ("url", "file", "size", "sha256", "member")


def _current_sources(record_path: Path, spec: dict, srcs: list, paths: list) -> list:
    """The record's sources when they describe the example's sources (their statuses are kept), else records
    built from the spec (kept when the file is still there)."""
    want = [_source_record(spec, s, d, "") for s, d in zip(srcs, paths)]
    old = _read_record(record_path).get("sources") or []
    if (isinstance(old, list) and len(old) == len(want) and all(isinstance(o, dict) for o in old)
            and all({k: o.get(k) for k in _SOURCE_KEYS} == {k: w.get(k) for k in _SOURCE_KEYS}
                    for o, w in zip(old, want))):
        return old
    return [{**w, "status": "kept" if d.is_file() else "not kept"} for w, d in zip(want, paths)]


def changes_text(spec: dict) -> str:
    """What fetch-example did to the public files, in plain words (CC BY 4.0 asks a shared copy to say so)."""
    asm = spec["assemble"]
    srcs = spec["sources"]
    parts = []

    def archive(path: str) -> str:              # Zenodo's API path ends in /content after the zip's name
        pp = PurePosixPath(path).parts
        return next((p for p in reversed(pp) if p.lower().endswith(".zip")), str(path))

    members = [f"member {s['member']} read out of the archive {archive(str(s['path']))}" for s in srcs
               if s.get("member")]
    files = [PurePosixPath(str(s["path"])).name for s in srcs if not s.get("member")]
    if members:
        parts.append("; ".join(members))
    if files:
        parts.append(("file " if len(files) == 1 else "files ") + ", ".join(files)
                     + (" stacked as channels" if len(srcs) > 1 else ""))
    parts.append(f"pages read with the declared axes {asm['axes']}")
    if asm.get("crop"):
        parts.append("cropped to " + ", ".join(f"{k} [{v[0]}, {v[1]})" for k, v in asm["crop"].items()))
    parts.append("channels named " + ", ".join(str(c) for c in asm["channels"]))
    sp = asm["spacing_zyx_um"]
    parts.append(f"voxel size {sp[0]} x {sp[1]} x {sp[2]} um recorded")
    parts.append("written as a zlib-compressed OME-TIFF by germquant assemble")
    return "; ".join(parts)


def _record_data(pack, spec: dict, image: Path, sources: list, image_info: dict) -> dict:
    asm = spec["assemble"]
    return {"format": 1, "pack": pack.id, "version": pack.version, "example": spec["id"], "recipe": spec["recipe"],
            "image": image.name, "pixels_sha256": str(spec["pixels_sha256"]), "shape_czyx": image_info.get("shape"),
            "dtype": image_info.get("dtype"), "channel_names": image_info.get("channel_names"),
            "spacing_zyx_um": image_info.get("spacing_zyx_um"), "crop": asm.get("crop"),
            "declared_axes": str(asm["axes"]), "licence": spec.get("licence"),
            "licence_spdx": spec.get("licence_spdx"), "licence_url": spec.get("licence_url"),
            "attribution": spec.get("attribution"), "changes": changes_text(spec), "sources": sources}


def _write_record(path: Path, pack, spec: dict, image: Path, sources: list, image_info: dict) -> None:
    """``<example>.json`` beside the image: where every byte came from, the licence (with its address), what
    was changed and the attribution."""
    from .. import __version__

    record = {**_record_data(pack, spec, image, sources, image_info), "written_at": _now(), "germquant": __version__}
    path.write_text(json.dumps(record, indent=1) + "\n", encoding="ascii")


def _record_stale(path: Path, pack, spec: dict, image: Path, sources: list, image_info: dict) -> bool:
    """Whether ``<example>.json`` is missing or no longer says what a record written now would say (another
    pack version, licence, attribution, channel names, voxel size, sources ...)."""
    old = _read_record(path)
    if not old:
        return True
    for key, value in _record_data(pack, spec, image, sources, image_info).items():
        have = old.get(key)
        if key == "spacing_zyx_um" and isinstance(value, list) and isinstance(have, list) and len(have) == len(value):
            if not all(math.isclose(float(a), float(b), rel_tol=1e-9, abs_tol=0.0) for a, b in zip(have, value)):
                return True
        elif json.loads(json.dumps(value)) != have:
            return True
    return False


def _read_record(record_path: Path) -> dict:
    try:
        rec = json.loads(record_path.read_text(encoding="ascii"))
        return rec if isinstance(rec, dict) else {}
    except (OSError, ValueError):
        return {}


# ---- the expect ranges -------------------------------------------------------------------------------------

def _summary_row(results: Path) -> dict:
    files = sorted(Path(results).glob("*__image_summary.csv"))
    if len(files) != 1:
        raise FileNotFoundError(f"{results}: expected one <image_id>__image_summary.csv, found {len(files)}")
    with open(files[0], newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if len(rows) != 1:
        raise ValueError(f"{files[0]}: expected one row, found {len(rows)}")
    return rows[0]


def check_expect(spec: dict, results: str | Path) -> dict:
    """Compare a results folder of the example's run with its ``expect`` ranges: n_objects is the image
    summary's n_nuclei, spots_per_object its mean_spots. Returns ``{"ok", "checks": [{key, value, range, ok}]}``
    (ok is True when the example sets no ranges)."""
    row = _summary_row(Path(results))
    values = {"n_objects": row.get("n_nuclei"), "spots_per_object": row.get("mean_spots")}
    checks = []
    for key, rng in (spec.get("expect") or {}).items():
        raw = values.get(key)
        try:
            val = float(raw)
        except (TypeError, ValueError):
            val = float("nan")
        ok = bool(rng[0] <= val <= rng[1])
        checks.append({"key": key, "value": val, "range": list(rng), "ok": ok})
    return {"ok": all(c["ok"] for c in checks), "checks": checks}
