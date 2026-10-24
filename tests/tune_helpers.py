"""Shared helpers for the spot tuning tests (germquant.tune; docs/TUNING.md; Train B unit P5-TUNE).

The planted examples: small two-channel fields (DAPI, foci) on the test pack's grid (0.5 x 0.2 x 0.2 um), with
gaussian nuclei and three kinds of foci on a flat background of 60 counts:

* true foci (amplitude 140 to 160, a LoG peak over background of about 3.5; a few at 260 to 270, about 5.4),
  which are the truth points;
* faint decoys (amplitude about 105, a LoG ratio of about 2.75) and fainter ones (about 72, about 2.2), which
  are not.

The test pack's basic recipe (classical nuclei, LoG spots at min_snr 2.0) segments every nucleus, and on the
default LoG grid {1.5, 2, 2.5, 3, 4, 5} only min_snr 3 finds exactly the true foci: the planted optimum.
The images are fed to process_image and to the tuners through the .nd2 reader seam (``planted_reader``), so
no image reader is needed; ``write_ome_tiff`` writes the same pixels as an OME-TIFF for the fresh-interpreter
tests, which read it through bioio.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from germquant import pipeline
from germquant.config import load_config
from germquant.io.nd2_reader import Stack

ROOT = Path(__file__).resolve().parents[1]
TEST_PACK = ROOT / "tests" / "fixtures" / "packs" / "test_pack"
BASIC = TEST_PACK / "config" / "profiles" / "basic.yaml"
SPACING = (0.5, 0.2, 0.2)
SHAPE = (12, 64, 128)
CHANNELS = ["DAPI", "foci"]
MICRON = chr(0xB5) + "m"            # the OME unit name (the schema refuses "um")
LOG_OPTIMUM = 3.0                   # spots.log.min_snr that finds exactly the true foci
T, D = True, False                  # a planted focus is truth (T) or a decoy (D)

# nucleus centre (z, y, x voxels) and its foci (dz, dy, dx voxels from the centre, amplitude, truth or decoy)
FIELD_A = [
    ((6, 18, 18), [(0, -4, -4, 150, T), (0, 4, 4, 140, T), (2, 0, 0, 105, D)]),
    ((6, 18, 52), [(0, 0, 0, 160, T)]),
    ((6, 20, 86), [(0, -4, -4, 150, T), (0, 4, 4, 155, T), (-2, -4, 4, 145, T), (2, 4, -4, 270, T)]),
    ((6, 20, 116), [(0, -4, 0, 150, T), (0, 4, 0, 72, D), (2, 0, 0, 72, D)]),
    ((6, 46, 30), [(0, -4, -4, 145, T), (0, 4, 4, 150, T), (-2, 0, 0, 108, D)]),
    ((6, 46, 66), [(0, 0, -4, 150, T), (0, 0, 4, 260, T), (-2, -4, 0, 160, T)]),
    ((6, 46, 104), [(0, 0, 0, 155, T), (2, 0, 0, 75, D)]),
]
FIELD_B = [
    ((6, 18, 20), [(0, -4, -4, 150, T), (0, 4, 4, 150, T), (0, -4, 4, 150, T)]),
    ((6, 20, 60), [(0, 0, 0, 155, T), (2, 4, 4, 106, D)]),
    ((6, 18, 100), [(0, -4, -4, 145, T), (0, 4, 4, 150, T), (-2, 0, 0, 72, D)]),
    ((6, 46, 40), [(0, 0, 0, 265, T), (0, -4, -4, 150, T)]),
    ((6, 46, 90), [(0, -4, -4, 150, T), (0, 4, 4, 150, T), (0, 4, -4, 148, T), (0, -4, 4, 152, T)]),
]
FIELDS = {"planted_a": (FIELD_A, 1), "planted_b": (FIELD_B, 2)}
# three nuclei of FIELD_A on a smaller field: the cheap example of the geometry and fresh-interpreter tests
SHAPE_TINY = (12, 40, 104)
TINY = (FIELD_A[:3], 7, SHAPE_TINY)
# single-plane fields for the 2D case (package P2b's io.planar runs; the LoG backend is refused there). On one
# plane SpotMAX's effect sizes of the 105-count decoys reach those of the faintest true foci, so these fields
# carry only the fainter decoys (72 to 75 counts): their effect sizes stay below 18.1, the true foci's above
# 24.1 (measured with P2b merged in a preview)
SHAPE_2D = (1, 64, 128)
FIELD_2D_A = [
    ((0, 18, 18), [(0, -4, -4, 150, T), (0, 4, 4, 140, T), (0, 4, -4, 72, D)]),
    ((0, 18, 52), [(0, 0, 0, 160, T)]),
    ((0, 20, 86), [(0, -4, -4, 150, T), (0, 4, 4, 155, T), (0, -4, 4, 145, T), (0, 4, -4, 270, T)]),
    ((0, 20, 116), [(0, -4, 0, 150, T), (0, 4, 0, 72, D)]),
    ((0, 46, 30), [(0, -4, -4, 145, T), (0, 4, 4, 150, T), (0, 4, -4, 75, D)]),
    ((0, 46, 66), [(0, 0, -4, 150, T), (0, 0, 4, 260, T), (0, -4, 0, 160, T)]),
    ((0, 46, 104), [(0, 0, 0, 155, T), (0, 4, 4, 75, D)]),
]
FIELD_2D_B = [
    ((0, 18, 20), [(0, -4, -4, 150, T), (0, 4, 4, 150, T), (0, -4, 4, 150, T)]),
    ((0, 20, 60), [(0, 0, 0, 155, T), (0, 4, 4, 74, D)]),
    ((0, 18, 100), [(0, -4, -4, 145, T), (0, 4, 4, 150, T), (0, -4, 4, 72, D)]),
    ((0, 46, 40), [(0, 0, 0, 265, T), (0, -4, -4, 150, T)]),
    ((0, 46, 90), [(0, -4, -4, 150, T), (0, 4, 4, 150, T), (0, 4, -4, 148, T), (0, -4, 4, 152, T)]),
]


def planted_pixels(field, seed: int, shape=SHAPE, spacing=SPACING) -> tuple[np.ndarray, list, list]:
    """(C, Z, Y, X) uint16 pixels of a planted field, the true foci and the decoys as (z, y, x) um."""
    rng = np.random.default_rng(seed)
    zz, yy, xx = np.indices(shape).astype(np.float32)
    dna = np.full(shape, 100.0, np.float32)
    foci = np.full(shape, 60.0, np.float32)
    truth, decoys = [], []
    for (cz, cy, cx), planted in field:
        r2 = ((zz - cz) * spacing[0]) ** 2 + ((yy - cy) * spacing[1]) ** 2 + ((xx - cx) * spacing[2]) ** 2
        dna += 1500.0 * np.exp(-r2 / (2 * 1.4 ** 2))
        for dz, dy, dx, amp, true in planted:
            fz, fy, fx = cz + dz, cy + dy, cx + dx
            f2 = (((zz - fz) * spacing[0]) / 0.35) ** 2 + (((yy - fy) * spacing[1]) / 0.15) ** 2 \
                + (((xx - fx) * spacing[2]) / 0.15) ** 2
            foci += amp * np.exp(-f2 / 2)
            (truth if true else decoys).append((fz * spacing[0], fy * spacing[1], fx * spacing[2]))
    dna += rng.normal(0, 8, shape)
    foci += rng.normal(0, 1.5, shape)
    return np.clip(np.stack([dna, foci]), 0, 65535).astype(np.uint16), truth, decoys


def planted_stack(name: str) -> tuple[Stack, list, list]:
    """The planted field ``name`` (planted_a, planted_b or tiny) as a Stack, its true foci and its decoys."""
    field, seed, shape = TINY if name == "tiny" else (*FIELDS[name], SHAPE)
    pix, truth, decoys = planted_pixels(field, seed, shape)
    stack = Stack(data=pix.astype(np.float32), spacing=SPACING, channel_names=list(CHANNELS),
                  path=Path(f"{name}.nd2"), spacing_ok=True)
    return stack, truth, decoys


def _strided(stack: Stack, xy_stride=1, z_range=None) -> Stack:
    data = stack.data
    dz, dy, dx = stack.spacing
    if z_range is not None:
        data = data[:, z_range[0]:z_range[1]]
    if xy_stride > 1:
        data = data[:, :, ::xy_stride, ::xy_stride]
        dy, dx = dy * xy_stride, dx * xy_stride
    return Stack(data=np.ascontiguousarray(data), spacing=(dz, dy, dx), channel_names=list(stack.channel_names),
                 path=stack.path, downsample_xy=xy_stride, spacing_ok=stack.spacing_ok)


def patch_planted_reader(mp, stacks: dict[str, Stack]) -> None:
    """Feed ``stacks`` (by file name) to process_image and to the tuners through the .nd2 reader seam of
    germquant.pipeline (``mp`` is a pytest MonkeyPatch); the stride and the z range are applied as
    read_stack applies them."""
    def meta(path, **kw):
        s = stacks[Path(path).name]
        c, z, y, x = s.data.shape
        return {"sizes": {"C": c, "Z": z, "Y": y, "X": x}, "spacing": s.spacing, "spacing_ok": True,
                "channel_names": list(s.channel_names), "dtype": "float32", "is_2d": False}

    def read(path, *, xy_stride=1, z_range=None, **kw):
        return _strided(stacks[Path(path).name], xy_stride, z_range)

    mp.setattr(pipeline, "read_nd2_metadata", meta)
    mp.setattr(pipeline, "read_stack", read)


def basic_recipe(**overrides):
    """The test pack's basic recipe (hash 02bb83cf74f9 unchanged when no override is given)."""
    cfg = load_config(BASIC)
    for key, value in overrides.items():
        cfg.set(key, value)
    return cfg


def finished_run(image: Path, cfg, out_dir: Path, **geometry) -> dict:
    """process_image of ``image`` into ``out_dir`` (one image per folder, as a project's results/<id>/)."""
    res = pipeline.process_image(image, cfg, out_dir, prov=None, **geometry)
    assert not any(":FAILED_" in f for f in res["qc_flags"]), res["qc_flags"]
    return res


def write_truth(points_path: Path, truth_um, result_dir: Path, image_id: str, *, xy_stride=1, z_range=None):
    """Write the planted truth (file-frame um) as the truth-point file a user would save on the run."""
    import tifffile

    from germquant.tune.points import write_points

    labels = tifffile.imread(result_dir / f"{image_id}__nuclei_labels.tif")
    z0 = z_range[0] if z_range else 0
    zyx = np.asarray(truth_um, float).reshape(-1, 3).copy()
    zyx[:, 0] -= z0 * SPACING[0]                      # napari's world coordinates of the run's layers
    spacing = (SPACING[0], SPACING[1] * xy_stride, SPACING[2] * xy_stride)
    return write_points(points_path, zyx, labels=labels, spacing=spacing, xy_stride=xy_stride, z_range=z_range)


def write_ome_tiff(path: Path, pixels: np.ndarray, spacing=SPACING, names=CHANNELS) -> Path:
    import tifffile

    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(path, pixels, ome=True, photometric="minisblack", compression="zlib",
                     metadata={"axes": "CZYX", "Channel": {"Name": list(names)},
                               "PhysicalSizeZ": spacing[0], "PhysicalSizeY": spacing[1], "PhysicalSizeX": spacing[2],
                               "PhysicalSizeZUnit": MICRON, "PhysicalSizeYUnit": MICRON, "PhysicalSizeXUnit": MICRON})
    return path


def counting_backend(mp, name: str) -> list:
    """Wrap the registered spot backend ``name`` so every call is counted; returns the list of calls."""
    from germquant.spots import backends as B

    B.discover_backends()
    orig = B._REGISTRY[name]
    calls: list = []

    def wrapper(*a, **k):
        calls.append(dict(k.get("cfg_section") or {}))
        return orig(*a, **k)

    for attr in ("input_roles", "provenance"):
        if hasattr(orig, attr):
            setattr(wrapper, attr, getattr(orig, attr))
    mp.setitem(B._REGISTRY, name, wrapper)
    return calls


class FakeProject:
    """The fields of germquant.project.Project the tuning entry reads (root, images), for tests that do not
    depend on package P5-PRJ."""

    def __init__(self, root: Path, images: list[dict], cfg):
        self.root = Path(root)
        self.images = images
        self.manifest = {"images": images}
        self.pack = {}
        self.warnings = []
        self.cfg = cfg


def fake_project_api(project: FakeProject, edits: list | None = None):
    """Replacements for (open_project, project_config) and (.., write_recipe) of germquant.project."""
    def open_project(root):
        if Path(root) != project.root:
            raise FileNotFoundError(f"no project at {root}")
        return project

    def project_config(p):
        return p.cfg

    def write_recipe(p, *, edits=None, unset=(), roles=None, origin="form", plugin_version=None):
        record = {"edits": dict(edits or {}), "origin": origin, "plugin_version": plugin_version}
        if edits_log is not None:
            edits_log.append(record)
        return "c0ffee000000"

    edits_log = edits
    return open_project, project_config, write_recipe
