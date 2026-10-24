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

Single-plane runs (package P2b). ``io.planar`` (false | auto | true) makes an image with one plane a 2D
image, where every quantity is 2D. decide() takes the decision right after the read: auto is planar when the
header has one plane (no Z axis, or Z of 1) and the stack read has one plane, so a 3D file cut to one plane
with --z-range stays 3D; true requires one plane and refuses more, naming --z-range and
segmentation.nuclei.project. dy and dx come from ``io.voxel_size_um``, else Stack.spacing_yx (a one-plane
file's own xy calibration, bioio reader), else the stack's spacing when it was readable, else 1 um per pixel
with the flag ``voxel:UNREADABLE_xy_assumed_1um``. Array code receives (PLANAR_DZ, dy, dx): the placeholder
reaches only code whose result does not depend on dz for one plane. The object is segmented on the plane
(segment_planar, the area gate), measured by measure_planar (volume_um3 NA, centroid_z_um 0.0, area_um2 and
equivalent_diameter_um), its spots are detected with ``planar=True`` (only spotmax declares planar_ok) and its
label, spot and mask images are written as YX (write_labels_yx, write_spots_yx, write_mask_yx).

What a single-plane run records, and only such a run: the flag ``read:planar``, voxel_dz_um NaN in every
table, ``planar: true`` and ``spacing_source`` in the run geometry, ``"planar": true`` in the read, segment,
measure and spots sub-hashes and the done-marker extra ``{"planar": true}``. ``io.voxel_size_um`` applies in
3D and projected runs too (voxel_override): spacing_source ``override`` and the flag ``voxel:override``, plus
``voxel:override_differs_from_file`` when the file's own calibration is more than 1 percent away. Time series
of single planes are refused in this release (PLANAR_SERIES_REFUSED).

Nothing here is reachable from a frozen config: ``segmentation.nuclei.project`` is read with the code
default ``none``, ``io.planar`` with the code default false and ``io.voxel_size_um`` with null, the other
keys are read only in a projected or single-plane run, and every existing run takes the existing path.
"""
from __future__ import annotations

import copy
import logging
import math
import re

import numpy as np

from .fsutil import long_path

log = logging.getLogger(__name__)

PROJECT_MODES = ("none", "max", "mean", "plane")
PLANE_SELECTORS = ("mid", "sharpest")
SHARPNESS_SIGMA_PX = 1.0     # the Gaussian smoothing before the Laplacian of the sharpest-plane score
PLANAR_MODES = ("false", "auto", "true")
PLANAR_DZ = 1.0              # the z spacing array code receives in a single-plane run; the tables write NaN
OVERRIDE_TOLERANCE = 0.01    # voxel:override_differs_from_file: the override is more than 1 percent off the file
UM = chr(0x00B5) + "m"       # the OME unit of length
PLANAR_SERIES_REFUSED = (
    "single-plane time series are refused in this release: io.planar is true, or auto on an image with one plane, "
    "and timepoints were asked for (--t, or t=). Run the series with io.planar false, which processes each "
    "timepoint as a stack of one plane, as before; or save the timepoints as separate single-plane images and run "
    "those (docs/PLANAR.md, 'Time series')")


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


def planar_setting(cfg) -> str:
    """io.planar as ``"false"`` (the code default, also for null), ``"auto"`` or ``"true"``. Any other value raises
    PlanarConfigError (recipe_keys refuses it first, with the same rule: a quoted "true" is text)."""
    raw = cfg.get("io.planar", None)
    if raw is None or raw is False:
        return "false"
    if raw is True:
        return "true"
    if isinstance(raw, str) and raw == "auto":
        return "auto"
    raise PlanarConfigError(f"io.planar must be false, auto or true (or null for false); got {raw!r}")


def voxel_size_setting(cfg) -> tuple[float | None, float, float] | None:
    """io.voxel_size_um as (dz or None, dy, dx) in micrometres, the file's own voxel size before xy_stride: a
    number is dy = dx, [dy, dx] leaves dz to the file, [dz, dy, dx] sets all three. None when the key is unset
    (null, the code default). A malformed value raises PlanarConfigError (recipe_keys refuses it first)."""
    raw = cfg.get("io.voxel_size_um", None)
    if raw is None:
        return None
    from .recipe_keys import voxel_size_problem

    problem = voxel_size_problem(raw)
    if problem:
        raise PlanarConfigError(problem)
    if isinstance(raw, (list, tuple)):
        vals = [float(v) for v in raw]
        return (None, vals[0], vals[1]) if len(vals) == 2 else (vals[0], vals[1], vals[2])
    return None, float(raw), float(raw)


# ---- single-plane runs: the decision and the spacing ------------------------------------------------

def one_plane_header(meta) -> bool:
    """Whether a header (read_image_metadata) describes one plane: its sizes have no Z axis (the nd2 reader leaves
    Z out of a 2D file) or a Z of 1 (bioio). is_2d is not used: it means "no Z axis" in one reader and "Z is 1" in
    the other. A header without sizes is not taken for one plane."""
    sizes = meta.get("sizes") if isinstance(meta, dict) else None
    if not sizes:
        return False
    z = sizes.get("Z")
    return z is None or int(z) == 1


def _file_yx(stack) -> tuple[float, float] | None:
    """The file's own (dy, dx) before xy_stride, or None when the file carries no xy calibration."""
    s = float(getattr(stack, "downsample_xy", 1) or 1)
    yx = getattr(stack, "spacing_yx", None)
    if yx is not None:
        return float(yx[0]) / s, float(yx[1]) / s
    if stack.spacing_ok:
        return float(stack.spacing[1]) / s, float(stack.spacing[2]) / s
    return None


def _differs(value: float, file_value: float) -> bool:
    return abs(float(value) - float(file_value)) > OVERRIDE_TOLERANCE * abs(float(file_value))


def voxel_override(cfg, stack) -> dict | None:
    """The spacing io.voxel_size_um gives a 3D or projected run (decision D6 of the Train B plan), or None when the
    key is unset, which is every existing run. ``{"spacing": (dz, dy, dx), "flags": [...], "spacing_source":
    "override"}``: dy and dx are the override times xy_stride, dz the override's own or else the file's. The flags
    are ``voxel:override``, plus ``voxel:override_differs_from_file`` when an axis the override sets is more than
    1 percent away from the file's own calibration of it, plus today's ``voxel:UNREADABLE_assumed_isotropic_1um``
    when dz is still the reader's 1 um guess (an override without dz on a file without a voxel size). A
    single-plane run takes its override through decide()."""
    setting = voxel_size_setting(cfg)
    if setting is None:
        return None
    dz_o, dy_o, dx_o = setting
    s = float(getattr(stack, "downsample_xy", 1) or 1)
    file_yx = _file_yx(stack)
    pairs = [] if file_yx is None else [(dy_o, file_yx[0]), (dx_o, file_yx[1])]
    if dz_o is not None and stack.spacing_ok:
        pairs.append((dz_o, float(stack.spacing[0])))
    flags = ["voxel:override"]
    if any(_differs(o, f) for o, f in pairs):
        flags.append("voxel:override_differs_from_file")
    if dz_o is None and not stack.spacing_ok:
        flags.append("voxel:UNREADABLE_assumed_isotropic_1um")      # dz stays the reader's 1 um guess
    dz = float(dz_o) if dz_o is not None else float(stack.spacing[0])
    return {"spacing": (dz, dy_o * s, dx_o * s), "flags": flags, "spacing_source": "override"}


def _xy_spacing(cfg, stack) -> tuple[float, float, str, list[str]]:
    """(dy, dx, spacing_source, flags) of a single-plane run, after xy_stride (design 4.3)."""
    s = float(getattr(stack, "downsample_xy", 1) or 1)
    setting = voxel_size_setting(cfg)
    if setting is not None:
        _, dy, dx = setting                  # a dz in the override has no place in a 2D run
        flags = ["voxel:override"]
        file_yx = _file_yx(stack)
        if file_yx is not None and (_differs(dy, file_yx[0]) or _differs(dx, file_yx[1])):
            flags.append("voxel:override_differs_from_file")
        return dy * s, dx * s, "override", flags
    yx = getattr(stack, "spacing_yx", None)
    if yx is not None:
        return float(yx[0]), float(yx[1]), "file_yx", []
    if stack.spacing_ok:
        return float(stack.spacing[1]), float(stack.spacing[2]), "file", []
    return s, s, "assumed", ["voxel:UNREADABLE_xy_assumed_1um"]


def decide(cfg, meta, stack) -> dict | None:
    """The single-plane decision, right after the read (design 4.2 and 4.3). None unless io.planar is auto or true,
    so every existing run takes the existing path. Otherwise ``{"planar", "dy", "dx", "spacing_source", "flags"}``:

    * auto: planar when the header has one plane (one_plane_header) and the stack read has one plane, so a 3D file
      cut to one plane with --z-range stays 3D; true: planar, and more than one plane read raises
      PlanarConfigError naming --z-range and segmentation.nuclei.project;
    * dy, dx and spacing_source (``override``, ``file_yx``, ``file`` or ``assumed``): io.voxel_size_um (times
      xy_stride), else Stack.spacing_yx, else stack.spacing[1:] when spacing_ok, else 1 um per file pixel;
    * flags, for a planar decision only: ``read:planar``, then ``voxel:override`` (and
      ``voxel:override_differs_from_file``) or ``voxel:UNREADABLE_xy_assumed_1um``. A 3D decision carries none;
      its run records what any 3D run records."""
    mode = planar_setting(cfg)
    if mode == "false":
        return None
    nz = int(np.shape(stack.data)[1])
    if mode == "true" and nz != 1:
        raise PlanarConfigError(
            f"io.planar is true but the image has {nz} planes as read, and a single-plane run needs one plane: read "
            f"one plane with --z-range Z Z+1 (z_range), or set io.planar to false or auto and segment a projection "
            f"of the stack with segmentation.nuclei.project (docs/PLANAR.md)")
    planar = nz == 1 and (mode == "true" or one_plane_header(meta))
    dy, dx, source, flags = _xy_spacing(cfg, stack)
    return {"planar": planar, "dy": dy, "dx": dx, "spacing_source": source,
            "flags": ["read:planar", *flags] if planar else []}


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


def measure_planar(labels, intensity: dict, spacing):
    """The object table of a single-plane run (design 4.6): measure_objects, unchanged, on the (1, Y, X) labels with
    spacing (PLANAR_DZ, dy, dx), then volume_um3 NA (one plane has no volume), centroid_z_um 0.0, and area_um2
    (the pixel count times dy * dx) and equivalent_diameter_um (the diameter of the circle of that area) added.
    The id column keeps measure_objects' name (``label``); an empty table is returned as it is."""
    from .measure import measure_objects

    df = measure_objects(labels, intensity, spacing, compute_surface=False)
    if df.empty or "label" not in df.columns:
        return df
    area, diameter = area_columns(df["label"].to_numpy(), labels, spacing)
    df["volume_um3"] = np.nan
    df["centroid_z_um"] = 0.0
    df["area_um2"] = area
    df["equivalent_diameter_um"] = diameter
    return df


# ---- single-plane writers (design 4.8; the 3D writers of pipeline.py are not edited) -------------------

def _yx(arr) -> np.ndarray:
    a = np.asarray(arr)
    if a.ndim == 3 and a.shape[0] == 1:
        a = a[0]
    if a.ndim != 2:
        raise ValueError(f"expected a (1, Y, X) or (Y, X) image, got shape {a.shape}")
    return a


def write_labels_yx(labels, path, spacing_yx) -> None:
    """int32 label image of a single-plane run as an OME-TIFF with axes YX and PhysicalSizeX and Y only (no Z), so
    Fiji, Imaris and bioio load it on the right grid. A failed tagged write falls back to a plain TIFF of the same
    array, as the 3D writer does; a failure of both is logged, never fatal."""
    try:
        import tifffile

        arr = _yx(labels).astype(np.int32)
        dy, dx = (float(s) for s in spacing_yx)
        try:
            tifffile.imwrite(long_path(path), arr, compression="zlib", ome=True, resolution=(1 / dx, 1 / dy),
                             metadata={"axes": "YX", "PhysicalSizeX": dx, "PhysicalSizeXUnit": UM,
                                       "PhysicalSizeY": dy, "PhysicalSizeYUnit": UM})
            return
        except Exception as e:  # noqa: BLE001 - fall back to the untagged file, never a stub
            log.warning("tagged label write failed (%s); writing a plain TIFF", e)
        tifffile.imwrite(long_path(path), arr, compression="zlib")
    except Exception as e:  # noqa: BLE001
        log.warning("could not save label image: %s", e)


def write_spots_yx(spots, shape_yx, spacing_yx, path, radius_um: float = 0.3) -> None:
    """uint16 (Y, X) image of a single-plane run with a disk of ``radius_um`` (at least one pixel) at each spot,
    value 65535, as an ImageJ TIFF calibrated in y and x (unit um, no z step). ``spots`` needs y_um and x_um."""
    try:
        import tifffile

        h, w = (int(s) for s in shape_yx[-2:])
        dy, dx = (float(s) for s in spacing_yx)
        img = np.zeros((h, w), np.uint16)
        if spots is not None and len(spots):
            ry, rx = max(1, int(round(radius_um / dy))), max(1, int(round(radius_um / dx)))
            oy, ox = np.ogrid[-ry:ry + 1, -rx:rx + 1]
            offs = np.argwhere((oy / ry) ** 2 + (ox / rx) ** 2 <= 1.0) - np.array([ry, rx])
            yc = np.clip(np.round(spots["y_um"].to_numpy(dtype=float) / dy).astype(int), 0, h - 1)
            xc = np.clip(np.round(spots["x_um"].to_numpy(dtype=float) / dx).astype(int), 0, w - 1)
            for y, x in zip(yc, xc):
                img[np.clip(y + offs[:, 0], 0, h - 1), np.clip(x + offs[:, 1], 0, w - 1)] = 65535
        tifffile.imwrite(long_path(path), img, compression="zlib", imagej=True, resolution=(1 / dx, 1 / dy),
                         metadata={"unit": "um", "axes": "YX"})
    except Exception as e:  # noqa: BLE001
        log.warning("could not save spots image: %s", e)


def write_mask_yx(mask, spacing_yx, path) -> None:
    """Binary mask of a single-plane run as a uint8 (Y, X) ImageJ TIFF (0 or 255), calibrated in y and x."""
    try:
        import tifffile

        img = (_yx(mask) > 0).astype(np.uint8) * 255
        dy, dx = (float(s) for s in spacing_yx)
        tifffile.imwrite(long_path(path), img, compression="zlib", imagej=True, resolution=(1 / dx, 1 / dy),
                         metadata={"unit": "um", "axes": "YX"})
    except Exception as e:  # noqa: BLE001
        log.warning("could not save mask image: %s", e)


# ---- gating ------------------------------------------------------------------------------------------

def check_planar(cfg, role_to_idx: dict, *, planar: bool = False) -> None:
    """Refuse, before any pixel is processed, what cannot run on 2D objects (design 4.4). Raises
    PlanarConfigError naming every problem, with the key that switches it off:

    * a stage without planar_ok (germline, axis, envelope, audit, staging, sc_trace, granule, coloc,
      partition, granule_tail) that is switched on and whose channel roles are all bound; a stage whose
      roles are missing skips itself as it always has;
    * a segmentation backend without planar_ok (cellpose and classical declare it, nuclei3d does not), unless
      segmentation.nuclei.method is classical, which never calls the backend;
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
    if str(cfg.get("segmentation.nuclei.method", "auto")) != "classical":
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


# ---- the single-plane run in process_image -----------------------------------------------------------

def prepare_planar(cfg, role_to_idx: dict) -> dict:
    """What process_image needs for a single-plane run, after decide() and before any stage runs: a
    segmentation.nuclei.project value is still checked and is then ignored, because one plane needs no projection
    (decision D7: the segment stage raises segment:project_ignored_planar); check_planar's refusals with
    ``planar=True`` (also a spot backend without planar_ok, the main instance and every spots.instances entry);
    and the area gate. ``{"min_area_um2", "project_ignored"}``."""
    ignored = projection_settings(cfg) is not None       # raises for an unknown mode or plane selector
    check_planar(cfg, role_to_idx, planar=True)
    return {"min_area_um2": min_area_setting(cfg), "project_ignored": ignored}
