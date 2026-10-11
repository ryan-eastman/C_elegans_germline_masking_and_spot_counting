#!/usr/bin/env python
"""Write the generated files of the outside packs (packs/cultured_cells, packs/budding_yeast), print their recipe
hashes, and fetch, pin and run their real examples (docs/PACKS.md; docs/DESIGN_system_packs.md, section 12, P3).

    python scripts/make_pack_examples.py [--pack cultured_cells|budding_yeast|all] [--base] [--example]
    python scripts/make_pack_examples.py --pack budding_yeast --draft-real [--examples DIR]
    python scripts/make_pack_examples.py --pack budding_yeast --run-real [--examples DIR] [--results DIR]

--base        config/profiles/pack_base.yaml: the frozen base every recipe of the pack extends, written from
              config/profiles/generic.yaml as resolved (germquant.packs.authoring), with the pack's own channel
              map, every worm switch off, the neutral file-name regex and the QC settings of a field that is not
              a gonad (qc.germline_outlier_factor 1e9; qc.min_nuclei 1 for cultured cells, 10 for yeast). Frozen
              once written: an edit changes every recipe hash of the pack, and the recipe goldens catch it.
--example     example/<name>.ome.tif: the small deterministic synthetic example (under 2 MB) the recipe goldens
              run on. cultured_cells: three channels (DAPI, gH2AX, 53BP1), four nuclei with foci, 16 x 144 x 176
              voxels at 0.3 x 0.2 x 0.2 um. budding_yeast: two channels (Bright, Cy3), six cells in focus on
              plane 4 of 13 (the middle plane would be 6) with 19 smFISH-like spots and two small bright specks
              between the cells that segmentation must not count, 13 x 160 x 224 voxels at 0.24 x 0.072 x 0.072
              um (the S-BIAD3750 SUT509 voxel size).
--draft-real  a maintainer's step: fetch the pack's real example as pack.yaml describes it, accepting sources
              and pixels whose digests are not pinned yet, and print the size, sha256 and pixels_sha256 to pin.
--run-real    fetch the real example (germquant pack fetch-example), run the example's recipe on it through the
              pack (api.run, after the pack check) and compare the counts with the example's expect ranges.
              Runs on the device germquant selects (CUDA_VISIBLE_DEVICES=-1 forces the CPU).

After --base or --example: copy the printed hashes into pack.yaml, run `germquant pack golden ID --update`, then
`germquant pack build-sums packs/ID` and `germquant pack lint packs/ID`.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from germquant.config import load_config  # noqa: E402
from germquant.packs.authoring import frozen_base_data, write_frozen_base  # noqa: E402

PACKS = REPO / "packs"
MICRON = chr(0xB5) + "m"          # the OME unit name; the schema refuses "um"
SETTINGS = {
    "cultured_cells": {"map": "config/channel_maps/cultured_dna_foci.yaml", "min_nuclei": 1,
                       "recipes": ("nuclei_foci", "nuclei_foci_3d"), "example": "example/cultured_synthetic.ome.tif"},
    "budding_yeast": {"map": "config/channel_maps/yeast_transmitted_foci.yaml", "min_nuclei": 10,
                      "recipes": ("cells_spots",), "example": "example/yeast_synthetic.ome.tif"},
}

# ---- cultured cells: four flat nuclei, gamma-H2AX foci and 53BP1 foci (half of them on the same sites) --------
CULT_SPACING = (0.3, 0.2, 0.2)
CULT_SHAPE = (16, 144, 176)
CULT_CHANNELS = ["DAPI", "gH2AX", "53BP1"]
# (z, y, x) centre in voxels and (rz, ry, rx) semi-axes in um
CULT_NUCLEI = [((7.5, 36, 44), (2.0, 5.0, 6.0)), ((7.5, 40, 130), (2.2, 5.4, 5.8)),
               ((7.5, 106, 46), (1.9, 4.8, 6.2)), ((7.5, 104, 128), (2.1, 5.2, 6.0))]


def cultured_pixels() -> np.ndarray:
    """(C, Z, Y, X) uint16 of the cultured cells example (module docstring)."""
    rng = np.random.default_rng(20261012)
    zz, yy, xx = np.indices(CULT_SHAPE).astype(np.float32)
    dna = np.full(CULT_SHAPE, 150.0, np.float32)
    g = np.full(CULT_SHAPE, 100.0, np.float32)
    b = np.full(CULT_SHAPE, 100.0, np.float32)
    for k, ((cz, cy, cx), (rz, ry, rx)) in enumerate(CULT_NUCLEI):
        r = np.sqrt(((zz - cz) * CULT_SPACING[0] / rz) ** 2 + ((yy - cy) * CULT_SPACING[1] / ry) ** 2
                    + ((xx - cx) * CULT_SPACING[2] / rx) ** 2)
        inside = 1.0 / (1.0 + np.exp((r - 1.0) / 0.08))
        texture = 1.0 + 0.08 * np.sin(xx * 0.9 + k) * np.cos(yy * 0.7 - k)
        dna += 1200.0 * inside * texture
        g += 100.0 * inside
        b += 80.0 * inside
        n = 3 + k
        for j in range(n):
            ang = 2 * np.pi * j / n + 0.5 * k
            rad = 0.35 + 0.25 * ((j + k) % 3) / 2
            fy = cy + rad * ry / CULT_SPACING[1] * np.sin(ang)
            fx = cx + rad * rx / CULT_SPACING[2] * np.cos(ang)
            fz = cz + (j % 3 - 1) * 1.5
            d2 = (((zz - fz) * CULT_SPACING[0]) / 0.35) ** 2 + (((yy - fy) * CULT_SPACING[1]) / 0.15) ** 2 \
                + (((xx - fx) * CULT_SPACING[2]) / 0.15) ** 2
            spot = np.exp(-d2 / 2)
            g += (1600.0 + 120.0 * j) * spot
            if j % 2 == 0:                                   # every second gamma-H2AX focus carries 53BP1 too
                b += (1300.0 + 100.0 * j) * spot
        fy, fx = cy - 0.5 * ry / CULT_SPACING[1], cx + 0.1 * rx / CULT_SPACING[2]    # one 53BP1-only focus
        d2 = (((zz - cz) * CULT_SPACING[0]) / 0.35) ** 2 + (((yy - fy) * CULT_SPACING[1]) / 0.15) ** 2 \
            + (((xx - fx) * CULT_SPACING[2]) / 0.15) ** 2
        b += 1400.0 * np.exp(-d2 / 2)
    dna += rng.normal(0, 10, CULT_SHAPE)
    g += rng.normal(0, 8, CULT_SHAPE)
    b += rng.normal(0, 8, CULT_SHAPE)
    return np.clip(np.stack([dna, g, b]), 0, 65535).astype(np.uint16)


# ---- budding yeast: six cells in focus on one plane, smFISH-like spots, two specks of debris -----------------
YEAST_SPACING = (0.24, 0.072, 0.072)
YEAST_SHAPE = (13, 160, 224)
YEAST_CHANNELS = ["Bright", "Cy3"]
YEAST_FOCUS_Z = 4                 # the planted sharpest plane (mid would be 6)
YEAST_CELL_R_UM = 1.9
YEAST_CELLS = [(40, 40), (38, 112), (42, 184), (120, 44), (122, 114), (118, 182)]    # (y, x) voxels
YEAST_DEBRIS = [(80, 76), (80, 148)]
YEAST_DEBRIS_R_UM = 0.5           # small bright specks: the classical segmenter on the sharpest plane finds none


def yeast_pixels() -> np.ndarray:
    """(C, Z, Y, X) uint16 of the budding yeast example (module docstring)."""
    from scipy import ndimage as ndi

    rng = np.random.default_rng(20261013)
    nz, ny, nx = YEAST_SHAPE
    yy, xx = np.indices((ny, nx)).astype(np.float32)
    r_px = YEAST_CELL_R_UM / YEAST_SPACING[1]
    focus = np.full((ny, nx), 1000.0, np.float32)
    for cy, cx in YEAST_CELLS:
        inside = (yy - cy) ** 2 + (xx - cx) ** 2 <= r_px ** 2
        focus[inside] = 1600.0 + 120.0 * np.sin(xx[inside] * 0.9) * np.cos(yy[inside] * 0.8)
    for cy, cx in YEAST_DEBRIS:
        focus[(yy - cy) ** 2 + (xx - cx) ** 2 <= (YEAST_DEBRIS_R_UM / YEAST_SPACING[1]) ** 2] = 2600.0
    bright = np.zeros(YEAST_SHAPE, np.float32)
    for z in range(nz):
        s = 1.5 * abs(z - YEAST_FOCUS_Z)
        bright[z] = ndi.gaussian_filter(focus, s) if s else focus
    bright += rng.normal(0.0, 5.0, bright.shape).astype(np.float32)

    zz, yy3, xx3 = np.indices(YEAST_SHAPE).astype(np.float32)
    cy3 = np.full(YEAST_SHAPE, 100.0, np.float32)
    for k, (cy, cx) in enumerate(YEAST_CELLS):
        cy3[(yy3 - cy) ** 2 + (xx3 - cx) ** 2 <= r_px ** 2] = 300.0
        n = 2 + k % 4
        for j in range(n):
            ang = 2 * np.pi * j / n + 0.6 * k
            py, px = cy + 0.55 * r_px * np.sin(ang), cx + 0.55 * r_px * np.cos(ang)
            pz = 3 + (3 * j + k) % 7
            d2 = (((zz - pz) * YEAST_SPACING[0]) / 0.3) ** 2 + (((yy3 - py) * YEAST_SPACING[1]) / 0.1) ** 2 \
                + (((xx3 - px) * YEAST_SPACING[2]) / 0.1) ** 2
            cy3 += (1800.0 + 150.0 * j) * np.exp(-d2 / 2)
    cy3 += rng.normal(0.0, 10.0, cy3.shape).astype(np.float32)
    return np.clip(np.stack([bright, cy3]), 0, 65535).astype(np.uint16)


def write_example(path: Path, data: np.ndarray, names: list[str], spacing) -> None:
    import tifffile

    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(path, data, ome=True, photometric="minisblack", compression="zlib",
                     metadata={"axes": "CZYX", "Channel": {"Name": list(names)},
                               "PhysicalSizeZ": spacing[0], "PhysicalSizeY": spacing[1], "PhysicalSizeX": spacing[2],
                               "PhysicalSizeZUnit": MICRON, "PhysicalSizeYUnit": MICRON, "PhysicalSizeXUnit": MICRON})


def example_data(pack_id: str) -> tuple[np.ndarray, list[str], tuple]:
    if pack_id == "cultured_cells":
        return cultured_pixels(), CULT_CHANNELS, CULT_SPACING
    return yeast_pixels(), YEAST_CHANNELS, YEAST_SPACING


# ---- real examples ----------------------------------------------------------------------------------------

def _raw_example(pack_id: str, example: str | None) -> dict:
    """The example spec as pack.yaml writes it (no validation: a draft may hold digests that are not pinned)."""
    from germquant.packs.format import load_manifest

    data = load_manifest(PACKS / pack_id)
    exs = [e for e in data.get("examples") or [] if isinstance(e, dict) and e.get("sources")]
    if example:
        exs = [e for e in exs if e.get("id") == example]
    if not exs:
        raise SystemExit(f"{pack_id}: no real example{' ' + example if example else ''} in pack.yaml")
    return exs[0]


def _pinned(value) -> bool:
    import re

    return bool(re.fullmatch(r"[0-9a-f]{64}", str(value or "")))


def draft_real(pack_id: str, example: str | None, out: Path) -> dict:
    """Fetch and assemble the real example, taking each digest that is pinned and printing every digest to pin."""
    from germquant import datafetch as D
    from germquant.assemble import assemble
    from germquant.packs.examples import source_name, source_url

    spec = _raw_example(pack_id, example)
    folder = out / pack_id / spec["id"]
    paths, pins = [], []
    for src in spec["sources"]:
        url = source_url(spec, src)
        dest = folder / "source" / source_name(src)
        sha = src.get("sha256") if _pinned(src.get("sha256")) else None
        if not (dest.is_file() and dest.stat().st_size == int(src["size"])
                and (sha is None or D.sha256_file(dest) == sha)):
            t0 = time.time()
            if src.get("member"):
                D.fetch_zip_member(url, src["member"], dest, size=int(src["size"]), sha256=sha,
                                   crc32=src.get("crc32"), zip_size=src.get("zip_size"))
            else:
                D.download(url, dest, size=int(src["size"]), sha256=sha)
            print(f"fetched {dest.name} in {time.time() - t0:.1f} s")
        got = D.sha256_file(dest)
        pins.append({"path": src["path"], "member": src.get("member"), "size": dest.stat().st_size, "sha256": got})
        paths.append(dest)
    asm = spec["assemble"]
    res = assemble(paths, asm["channels"], asm["spacing_zyx_um"], folder / f"{spec['id']}.ome.tif",
                   axes=str(asm["axes"]), crop=asm.get("crop"))
    print(json.dumps({"sources": pins, "pixels_sha256": res["pixels_sha256"], "shape_czyx": res["shape"],
                      "dtype": res["dtype"], "image": res["path"]}, indent=1))
    return res


def run_real(pack_id: str, example: str | None, out: Path, results: Path | None) -> int:
    """Fetch the real example through the pack, run its recipe and check the expect ranges."""
    from germquant import api, packs
    from germquant.packs.examples import check_expect

    pack = packs.find_pack(pack_id)
    rec = packs.fetch_example(pack, example, out, say=print)
    st = packs.pack_status(pack, recipe=rec["recipe"], image=rec["image"])
    for w in st["warnings"]:
        print(f"warning: {w['what']}")
    if not st["usable"]:
        for b in st["blockers"]:
            print(f"BLOCKER: {b['what']} (fix: {b['fix']})")
        return 2
    res_dir = (results or Path(rec["image"]).parent / f"results_{rec['recipe']}")
    t0 = time.time()
    res = api.run(rec["image"], pack=pack.ref, recipe=rec["recipe"], out=res_dir)
    took = time.time() - t0
    manifest = json.loads((res_dir / "run_manifest.json").read_text(encoding="utf-8"))
    chk = check_expect(rec, res_dir)
    print(json.dumps({"pack": pack.ref, "example": rec["id"], "recipe": rec["recipe"], "seconds": round(took, 1),
                      "n_objects": res["n_nuclei"], "qc_pass": res["qc_pass"], "qc_flags": res["qc_flags"],
                      "expect": chk, "pack_record": manifest.get("pack"), "runtime": manifest.get("runtime"),
                      "segmentation_method": res["tables"]["image_summary"]["segmentation_method"].iloc[0],
                      "results": str(res_dir)}, indent=1, default=str))
    return 0 if chk["ok"] else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pack", choices=["cultured_cells", "budding_yeast", "all"], default="all")
    ap.add_argument("--base", action="store_true", help="(re)write config/profiles/pack_base.yaml")
    ap.add_argument("--example", action="store_true", help="(re)write the synthetic example")
    ap.add_argument("--draft-real", action="store_true", help="fetch the real example and print the digests to pin")
    ap.add_argument("--run-real", action="store_true", help="fetch the real example, run it and check expect")
    ap.add_argument("--real-example", default=None, help="the real example's id (default: the pack's first)")
    ap.add_argument("--examples", default=None, help="the examples folder (default: $GERMQUANT_EXAMPLES or the "
                                                     "per-user folder)")
    ap.add_argument("--results", default=None, help="--run-real: the results folder")
    a = ap.parse_args(argv)
    ids = list(SETTINGS) if a.pack == "all" else [a.pack]
    if a.draft_real or a.run_real:
        from germquant.packs.examples import examples_dir

        out = Path(a.examples) if a.examples else examples_dir()
        rc = 0
        for pid in ids:
            if a.draft_real:
                draft_real(pid, a.real_example, out)
            if a.run_real:
                rc = max(rc, run_real(pid, a.real_example, out, Path(a.results) if a.results else None))
        return rc
    for pid in ids:
        s, folder = SETTINGS[pid], PACKS / pid
        if a.base:
            data = frozen_base_data(REPO / "config" / "profiles" / "generic.yaml", channel_map=s["map"],
                                    min_nuclei=s["min_nuclei"])
            print("wrote", write_frozen_base(folder / "config" / "profiles" / "pack_base.yaml", data))
        if a.example:
            out = folder / s["example"]
            pixels, names, spacing = example_data(pid)
            write_example(out, pixels, names, spacing)
            print("wrote", out, out.stat().st_size, "bytes")
        if (folder / "config" / "profiles" / "pack_base.yaml").is_file():
            for recipe in s["recipes"]:
                print(f"{pid} {recipe}: expected_config_hash {load_config(folder / 'config' / 'profiles' / f'{recipe}.yaml').hash}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
