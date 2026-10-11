#!/usr/bin/env python
"""Write the generated files of the test pack (tests/fixtures/packs/test_pack) and print its recipe hashes.

    python scripts/make_test_pack.py [--example] [--base]

--base     config/profiles/pack_base.yaml: the frozen base, written from config/profiles/generic.yaml as resolved
           (germquant.packs.authoring), with the pack's channel map, every worm switch off, the neutral regex and
           the outside QC settings. Frozen once written: an edit changes every recipe hash of the pack.
--example  example/test_pack_synthetic.ome.tif: a small deterministic two-channel field (DAPI and foci, five
           nuclei with two puncta each, 12 x 64 x 96 voxels at 0.5 x 0.2 x 0.2 um), under 2 MB.

Then copy the printed hashes into pack.yaml, run `germquant pack golden test_pack --update` with
GERMQUANT_PACKS=tests/fixtures/packs, and `germquant pack build-sums tests/fixtures/packs/test_pack`.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from germquant.config import load_config  # noqa: E402
from germquant.packs.authoring import frozen_base_data, write_frozen_base  # noqa: E402

PACK = REPO / "tests" / "fixtures" / "packs" / "test_pack"
SPACING = (0.5, 0.2, 0.2)
SHAPE = (12, 64, 96)
CENTRES = [(6, 18, 18), (6, 18, 50), (6, 20, 80), (6, 46, 30), (6, 46, 66)]


def example_pixels() -> np.ndarray:
    """(C, Z, Y, X) uint16: DAPI blobs (sigma 1.4 um) and two foci per nucleus (sigma 0.15 um in xy, 0.35 um in z)."""
    rng = np.random.default_rng(20261011)
    z, y, x = SHAPE
    zz, yy, xx = np.indices(SHAPE).astype(np.float32)
    dna = np.full(SHAPE, 100.0, np.float32)
    foci = np.full(SHAPE, 60.0, np.float32)
    for k, (cz, cy, cx) in enumerate(CENTRES):
        r2 = ((zz - cz) * SPACING[0]) ** 2 + ((yy - cy) * SPACING[1]) ** 2 + ((xx - cx) * SPACING[2]) ** 2
        dna += 1500.0 * np.exp(-r2 / (2 * 1.4 ** 2))
        for dy, dx in ((-3, -3), (3, 4)):
            fz, fy, fx = cz, cy + dy, cx + dx
            f2 = (((zz - fz) * SPACING[0]) / 0.35) ** 2 + (((yy - fy) * SPACING[1]) / 0.15) ** 2 \
                + (((xx - fx) * SPACING[2]) / 0.15) ** 2
            foci += (2500.0 + 150.0 * k) * np.exp(-f2 / 2)
    dna += rng.normal(0, 8, SHAPE)
    foci += rng.normal(0, 5, SHAPE)
    return np.clip(np.stack([dna, foci]), 0, 65535).astype(np.uint16)


MICRON = chr(0xB5) + "m"          # the OME unit name; the schema refuses "um" (bioio-ome-tiff then drops the calibration)


def write_example(path: Path) -> None:
    import tifffile

    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(path, example_pixels(), ome=True, photometric="minisblack", compression="zlib",
                     metadata={"axes": "CZYX", "Channel": {"Name": ["DAPI", "foci"]},
                               "PhysicalSizeZ": SPACING[0], "PhysicalSizeY": SPACING[1], "PhysicalSizeX": SPACING[2],
                               "PhysicalSizeZUnit": MICRON, "PhysicalSizeYUnit": MICRON, "PhysicalSizeXUnit": MICRON})


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", action="store_true", help="(re)write config/profiles/pack_base.yaml")
    ap.add_argument("--example", action="store_true", help="(re)write the synthetic example")
    a = ap.parse_args(argv)
    if a.base:
        data = frozen_base_data(REPO / "config" / "profiles" / "generic.yaml",
                                channel_map="config/channel_maps/test_dna_foci.yaml", min_nuclei=1)
        print("wrote", write_frozen_base(PACK / "config" / "profiles" / "pack_base.yaml", data))
    if a.example:
        out = PACK / "example" / "test_pack_synthetic.ome.tif"
        write_example(out)
        print("wrote", out, out.stat().st_size, "bytes")
    for recipe in ("basic", "with_model"):
        p = PACK / "config" / "profiles" / f"{recipe}.yaml"
        if (PACK / "config" / "profiles" / "pack_base.yaml").is_file():
            print(f"{recipe}: expected_config_hash {load_config(p).hash}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
