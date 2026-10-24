"""The 2D path: objects segmented on one 2D image (docs/DESIGN_system_packs.md section 4, docs/PLANAR.md).

Projected runs (package P2a). ``segmentation.nuclei.project`` (none | max | mean | plane) makes the segment
stage segment a 3D stack's objects on ONE 2D image of the segmentation channel: its maximum or mean
projection, or one plane picked by ``segmentation.nuclei.project_plane`` (mid | sharpest | an index). The
(1, Y, X) labels are then copied through every plane and everything downstream stays 3D: spots are counted
in 3D and count for the object whose column they fall in, intensities are means over that column. This is
how the yeast smFISH deposit of S-BIAD3750 was made (Cell-ACDC segmented one plane of the transmitted
channel and SpotMAX stacked the 2D masks through z), and how a thin monolayer stack is segmented.

What a projected run records, and only a projected run: the flag ``segment:projected_<mode>``,
``segment_projection`` {mode, plane, z} in the run geometry (run_manifest.json and <image_id>__stages.json),
``volume_um3`` as NA (a prism is not a cell volume) and the nuclei columns ``area_um2`` and
``equivalent_diameter_um`` from the 2D labels (schema.PLANAR_COLUMNS, written after the declared columns).
The montage's outline panel sits on the 2D image the segmenter saw.

Refuse, never skip silently: ``check_planar`` runs right after the read and raises PlanarConfigError (a
ValueError, fatal like a failed read) when a 3D-only stage (Stage.planar_ok False) is switched on with its
channel roles bound, or when the segmentation backend does not declare ``planar_ok``. An index out of range
raises before segmentation.

Nothing here is reachable from a frozen config: ``segmentation.nuclei.project`` is read with the code
default ``none``, the other keys are read only in a projected run, and every existing run takes the
existing path. Single-plane images (io.planar, package P2b) reuse segment_planar.
"""
from __future__ import annotations

import copy
import logging
import math
import re

import numpy as np

log = logging.getLogger(__name__)

PROJECT_MODES = ("none", "max", "mean", "plane")
PLANE_SELECTORS = ("mid", "sharpest")
SHARPNESS_SIGMA_PX = 1.0     # the Gaussian smoothing before the Laplacian of the sharpest-plane score


class PlanarConfigError(ValueError):
    """A 2D-path setting that cannot run: a 3D-only stage switched on with its roles bound, a backend
    without planar_ok, an unknown projection mode or plane selector, a plane index out of range."""


def planar_ok(fn) -> bool:
    """Whether a backend callable (or a Stage) declares that it handles 2D objects (its planar_ok attribute,
    False when absent: the input_roles pattern of the spot backends)."""
    return bool(getattr(fn, "planar_ok", False))


# ---- settings ---------------------------------------------------------------------------------------

def parse_plane_selector(value) -> str | int:
    """``mid``, ``sharpest`` or a plane index (an int, or a string of digits) from
    segmentation.nuclei.project_plane. Whether an index lies inside the stack is checked by resolve_plane,
    which knows the number of planes."""
    if isinstance(value, bool):
        raise PlanarConfigError(f"segmentation.nuclei.project_plane must be mid, sharpest or a plane index; "
                                f"got {value!r}")
    if isinstance(value, (int, np.integer)):
        return int(value)
    text = str(value).strip().lower()
    if text in PLANE_SELECTORS:
        return text
    if re.fullmatch(r"-?[0-9]+", text):
        return int(text)
    raise PlanarConfigError(f"segmentation.nuclei.project_plane must be mid, sharpest or a plane index; "
                            f"got {value!r}")


def projection_settings(cfg) -> dict | None:
    """``{"mode", "plane"}`` for a projected run, None for every other run (the code default ``none``).
    ``plane`` is the parsed segmentation.nuclei.project_plane (code default ``mid``) in mode ``plane``. The
    projections use every plane, so for them the key is not read at all and ``plane`` is None (docs/PLANAR.md:
    project_plane is ignored by max and mean). An unknown mode, or in mode ``plane`` an unknown selector,
    raises PlanarConfigError."""
    raw = cfg.get("segmentation.nuclei.project", "none")
    if raw is None or raw is False:           # YAML null, or an unquoted no / off
        return None
    mode = str(raw).strip().lower()
    if mode == "none":
        return None
    if mode not in PROJECT_MODES:
        raise PlanarConfigError(f"segmentation.nuclei.project must be one of {', '.join(PROJECT_MODES)}; "
                                f"got {raw!r}")
    if mode != "plane":
        return {"mode": mode, "plane": None}
    return {"mode": mode, "plane": parse_plane_selector(cfg.get("segmentation.nuclei.project_plane", "mid"))}


def min_area_setting(cfg) -> float:
    """segmentation.nuclei.min_area_um2, the size gate of 2D objects (code default 0.0: no gate)."""
    value = float(cfg.get("segmentation.nuclei.min_area_um2", 0.0) or 0.0)
    if not value >= 0.0:
        raise PlanarConfigError(f"segmentation.nuclei.min_area_um2 must be 0 or more; got {value!r}")
    return value


# ---- the 2D image ------------------------------------------------------------------------------------

def sharpness(plane) -> float:
    """The focus score of one plane: the variance of the Laplacian of the plane smoothed by a Gaussian of
    sigma 1 px, both filters in float32 (the variance is summed in float64)."""
    from scipy import ndimage as ndi

    p = np.asarray(plane, dtype=np.float32)
    lap = ndi.laplace(ndi.gaussian_filter(p, sigma=SHARPNESS_SIGMA_PX))
    return float(np.var(lap, dtype=np.float64))


def plane_scores(img) -> np.ndarray:
    """sharpness() of every plane of a (Z, Y, X) image; a plane whose score is not finite scores -inf."""
    arr = _zyx(img)
    scores = np.array([sharpness(arr[z]) for z in range(arr.shape[0])], dtype=np.float64)
    return np.where(np.isfinite(scores), scores, -np.inf)


def resolve_plane(img, selector) -> int:
    """The plane index a selector picks in a (Z, Y, X) image: ``mid`` is Z // 2, ``sharpest`` the plane with
    the highest sharpness() (ties go to the lowest z), an index itself. An index outside 0 to Z - 1 raises
    PlanarConfigError."""
    arr = _zyx(img)
    nz = int(arr.shape[0])
    sel = parse_plane_selector(selector)
    if sel == "mid":
        return nz // 2
    if sel == "sharpest":
        return int(np.argmax(plane_scores(arr)))          # argmax returns the first maximum: the lowest z
    if not 0 <= sel < nz:
        raise PlanarConfigError(f"segmentation.nuclei.project_plane {sel} is out of range: the stack has {nz} "
                                f"plane(s), indices 0 to {nz - 1}")
    return int(sel)


def project_image(img, mode: str, plane=None) -> np.ndarray:
    """The (Y, X) image a projected run segments: ``max`` the maximum projection (the input dtype), ``mean``
    the mean projection (summed in float64, returned as float32), ``plane`` the plane that ``plane`` (a
    selector or an index, resolve_plane) picks, as a copy."""
    arr = _zyx(img)
    if mode == "max":
        return arr.max(axis=0)
    if mode == "mean":
        return arr.mean(axis=0, dtype=np.float64).astype(np.float32)
    if mode == "plane":
        return np.array(arr[resolve_plane(arr, "mid" if plane is None else plane)])
    raise PlanarConfigError(f"projection mode must be max, mean or plane; got {mode!r}")


def _zyx(img) -> np.ndarray:
    arr = np.asarray(img)
    if arr.ndim != 3:
        raise ValueError(f"expected a (Z, Y, X) image, got shape {arr.shape}")
    return arr


# ---- segmentation ------------------------------------------------------------------------------------

def segment_planar(
    img_1yx,
    spacing,
    *,
    method: str = "auto",
    cellpose_model: str = "cpsam",
    diameter_um: float = 3.0,
    min_area_um2: float = 0.0,
    device: str | None = None,
    backend: str | None = None,
    cfg_section: dict | None = None,
) -> tuple[np.ndarray, str]:
    """Segment ONE plane: ``(labels (1, Y, X) int32, method_used)``. ``img_1yx`` is (1, Y, X) or (Y, X);
    ``spacing`` is (dz, dy, dx) (dz plays no part in a 2D result).

    The control flow of segment.segment_nuclei, which it mirrors without touching: method ``cellpose`` and
    ``auto`` run the backend that ``backend`` names (an unknown name raises UnknownBackendError, never a
    fallback), ``auto`` falls back to the classical watershed when the backend fails, ``cellpose`` re-raises,
    ``classical`` never calls the backend, and a backend's output is held to the contract (an integer label
    image the shape of the input). A backend that declares ``planar_ok`` gets the keyword ``planar=True``
    (cellpose then runs its 2D call); one that does not is called as for a 3D volume of one plane, as the
    input_roles pattern of the spot backends does, and check_planar refuses such a backend in a pipeline
    run. The size gate is drop_small_area (``min_area_um2``; 0 = no gate), followed by a consecutive
    relabel; segmentation.nuclei.min_volume_um3 does not apply to 2D objects."""
    from .segment import nuclei
    from .segment.backends import DEFAULT_BACKEND, get_backend

    img = np.asarray(img_1yx)
    if img.ndim == 2:
        img = img[None]
    if img.ndim != 3 or img.shape[0] != 1:
        raise ValueError(f"segment_planar takes one plane, (1, Y, X) or (Y, X); got shape {np.shape(img_1yx)}")
    if method in ("cellpose", "auto"):
        name = (str(backend).strip() if backend is not None else "") or DEFAULT_BACKEND
        run = get_backend(name)          # a misspelt backend is a config error: raise, never fall back
        what = "Cellpose (2D)" if name == DEFAULT_BACKEND else f"segmentation backend {name!r} (2D)"
        section = copy.deepcopy(cfg_section) if cfg_section else {}
        kw = {"planar": True} if planar_ok(run) else {}
        try:
            labels = run(img, spacing, model=cellpose_model, diameter_um=diameter_um, device=device,
                         cfg_section=section, **kw)
            labels = nuclei._check_labels(labels, img.shape, name)
            return drop_small_area(labels, spacing, min_area_um2), name
        except Exception as e:
            if method == "cellpose":
                raise
            log.warning("%s unavailable (%s); using classical watershed fallback.", what, e)
    labels = nuclei._classical(img, spacing, diameter_um)
    return drop_small_area(labels, spacing, min_area_um2), "classical"


def drop_small_area(labels, spacing, min_area_um2: float) -> np.ndarray:
    """Drop 2D objects below a physical area, then relabel consecutively: an object with fewer than
    int(min_area_um2 / (dy * dx)) pixels goes, the rule segment.nuclei._drop_small applies to volumes.
    ``spacing`` is (dz, dy, dx) or (dy, dx); 0 keeps every object (the relabel still runs)."""
    from skimage.segmentation import relabel_sequential

    lab = np.array(labels, dtype=np.int32)
    dy, dx = float(spacing[-2]), float(spacing[-1])
    min_px = max(1, int(float(min_area_um2) / (dy * dx)))
    counts = np.bincount(lab.ravel())
    small = np.where(counts < min_px)[0]
    small = small[small != 0]
    if small.size:
        lab[np.isin(lab, small)] = 0
    lab, _, _ = relabel_sequential(lab)
    return lab.astype(np.int32)


def repeat_through_z(labels_1yx, nz: int) -> np.ndarray:
    """The (nz, Y, X) label image of a projected run: the (1, Y, X) labels copied into every plane."""
    lab = np.asarray(labels_1yx)
    if lab.ndim == 2:
        lab = lab[None]
    if lab.ndim != 3 or lab.shape[0] != 1:
        raise ValueError(f"expected (1, Y, X) labels, got shape {lab.shape}")
    return np.repeat(lab.astype(np.int32), int(nz), axis=0)


# ---- measurement -------------------------------------------------------------------------------------

def area_columns(nucleus_ids, labels_2d, spacing) -> tuple[np.ndarray, np.ndarray]:
    """(area_um2, equivalent_diameter_um) of each id: its pixel count in the 2D labels times dy * dx, and the
    diameter of the circle of that area. An id absent from the labels has area 0."""
    lab = np.asarray(labels_2d)
    if lab.ndim == 3:
        lab = lab[0]
    ids = np.asarray(list(nucleus_ids), dtype=np.int64)
    counts = np.bincount(lab.ravel().astype(np.int64), minlength=int(ids.max()) + 1 if ids.size else 1)
    dy, dx = float(spacing[-2]), float(spacing[-1])
    area = counts[ids].astype(np.float64) * dy * dx if ids.size else np.zeros(0)
    return area, 2.0 * np.sqrt(area / math.pi)


def projected_measurements(nuclei, labels_2d, spacing):
    """The nuclei table of a projected run, from measure_objects' table on the copied labels (whose
    intensities are means over each 3D column and whose centroid_z_um is the stack middle, never NaN):
    volume_um3 becomes NA, because a prism is not a cell volume, and area_um2 and equivalent_diameter_um
    are added from the 2D labels. An empty table is returned as it is."""
    if nuclei.empty or "nucleus_id" not in nuclei.columns:
        return nuclei
    df = nuclei.copy()
    area, diameter = area_columns(df["nucleus_id"].to_numpy(), labels_2d, spacing)
    df["volume_um3"] = np.nan
    df["area_um2"] = area
    df["equivalent_diameter_um"] = diameter
    return df


# ---- gating ------------------------------------------------------------------------------------------

def check_planar(cfg, role_to_idx: dict, *, planar: bool = False) -> None:
    """Refuse, before any pixel is processed, what cannot run on 2D objects (design 4.4). Raises
    PlanarConfigError naming every problem, with the key that switches it off:

    * a stage without planar_ok (germline, axis, envelope, audit, staging, sc_trace, granule, coloc,
      partition, granule_tail) that is switched on and whose channel roles are all bound; a stage whose
      roles are missing skips itself as it always has;
    * a segmentation backend without planar_ok (cellpose and classical declare it, nuclei3d does not), unless
      segmentation.nuclei.method is classical, which never calls the backend, or segmentation.nuclei.labels_from
      is set, which segments nothing (docs/LABELS.md);
    * with ``planar=True`` (single-plane runs, P2b) also a spot backend without planar_ok, for the main spots
      instance and every spots.instances entry whose role is bound. Projected runs count spots in 3D, so
      any spot backend may run there.
    """
    from .segment import backends as seg_backends
    from .stages import STAGES, missing_roles, stage_enabled

    what = "single-plane run" if planar else "projected run"
    problems: list[str] = []
    for st in STAGES:
        if st.planar_ok or not stage_enabled(cfg, st.name) or missing_roles(role_to_idx, st.name):
            continue
        how = f"{st.config_key}: false" + (f" or --no-{st.name}" if st.cli_switchable else "")
        problems.append(f"the {st.name} stage is 3D-only and is switched on with its channel roles bound "
                        f"(switch it off with {how})")
    if (str(cfg.get("segmentation.nuclei.method", "auto")) != "classical"
            and cfg.get("segmentation.nuclei.labels_from") is None):
        name = str(cfg.get("segmentation.nuclei.backend", seg_backends.DEFAULT_BACKEND) or "").strip() \
            or seg_backends.DEFAULT_BACKEND
        fn = seg_backends.get_backend(name)      # an unknown name raises here, as the segment stage would
        if not planar_ok(fn):
            problems.append(f"the segmentation backend {name!r} does not declare planar_ok, so it cannot segment "
                            f"a 2D image (cellpose and classical can)")
    if planar and stage_enabled(cfg, "spots"):
        from .spots import backends as spot_backends

        main = cfg.get("spots.backend", spot_backends.DEFAULT_BACKEND)
        todo = [("spots", main, "foci")]
        for inst in cfg.get("spots.instances", []) or []:
            if isinstance(inst, dict):
                todo.append((f"spots_{inst.get('name')}", inst.get("backend", main), inst.get("role")))
        for where, name, role in todo:
            if role is None or role_to_idx.get(str(role)) is None:
                continue                         # the instance skips itself for lack of its channel
            if not planar_ok(spot_backends.get_backend(str(name))):
                problems.append(f"{where}: the spot backend {name!r} does not declare planar_ok")
    if problems:
        raise PlanarConfigError(f"{what} refused before any pixel was processed: " + "; ".join(problems))


# ---- the projected run in process_image -------------------------------------------------------------

def prepare_projection(cfg, stack, role_to_idx: dict) -> dict | None:
    """What process_image needs for a projected run, right after the read and before any stage runs: the
    settings, check_planar's refusals and, in mode ``plane``, the plane index ``z`` resolved on the channel
    the segment stage will segment (the segmentation role, or channel 0 when it does not resolve, as the
    segment stage falls back; nothing is resolved when the role is missing and
    segmentation.nuclei.require_channel is true, because the segment stage then fails the image). None
    for every other run, which reads, checks and records nothing new."""
    settings = projection_settings(cfg)
    if settings is None:
        return None
    check_planar(cfg, role_to_idx)
    z = None
    if settings["mode"] == "plane":
        seg_role = str(cfg.get("segmentation.nuclei.channel", "dna") or "dna")
        img = stack.channel(role_to_idx.get(seg_role))
        if img is None and not bool(cfg.get("segmentation.nuclei.require_channel", False)):
            img = stack.data[0]
        if img is not None:
            z = resolve_plane(img, settings["plane"])
    return {"mode": settings["mode"], "plane": settings["plane"], "z": z,
            "record": {"mode": settings["mode"], "plane": settings["plane"], "z": z},
            "flag": f"segment:projected_{settings['mode']}", "min_area_um2": min_area_setting(cfg)}
