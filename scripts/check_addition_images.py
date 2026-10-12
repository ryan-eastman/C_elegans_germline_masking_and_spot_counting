"""Protocol changelog entry 9 (d)(1) and (d)(2) for the additions: which channel holds the foci, and whether the
Imaris project holds the .nd2's own pixels. Reads only the additions' files (never a 2025-11-05 file).

    python scripts/check_addition_images.py --out DIR [--additions tests/golden/rad51_additions.json] [--tags TAG ...]
        [--planes 0 53 105]

For every addition of the record:

* (d)(1) max projections over z of the three .nd2 channels, written as one PNG per gonad (the whole field,
  downsampled, and a full-resolution crop around the median of the Imaris surface centres), each channel
  stretched to its own 1st to 99.9th percentile, for a look by eye: 545 must show discrete nuclear foci and 477
  synaptonemal-complex tracks. The default channel map (config/channel_maps/n2_dapi_syp_rad51_pgl1.yaml, match by
  name) is resolved against the .nd2 channel names and must put foci on 545.
* (d)(2) planes of channels 0 to 2 of the .ims (ResolutionLevel 0, the image part of the padded array) against the
  same .nd2 planes: equal pixel for pixel, or the number of differing pixels.

Writes DIR/image_checks.json (one record per gonad) and DIR/<tag>__mips.png. Exit code 0 when every plane is equal
and the map resolves foci to 545, else 1.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import numpy as np  # noqa: E402

from germquant.validate import benchmark_split as bs  # noqa: E402
from germquant.validate import imaris_ims as I  # noqa: E402

DEFAULT_MAP = REPO / "config" / "channel_maps" / "n2_dapi_syp_rad51_pgl1.yaml"


def compare_planes(ims_path, nd2_channels: dict, planes) -> list[dict]:
    """Planes ``planes`` of the .ims channels against the .nd2 channels (``nd2_channels``: index -> ZYX)."""
    import h5py

    out = []
    with h5py.File(str(ims_path), "r") as f:
        geom = I.read_geometry(f)
        X, Y, Z = (int(v) for v in geom["shape_xyz"])
        for c, arr in nd2_channels.items():
            if tuple(arr.shape) != (Z, Y, X):
                out.append({"channel": int(c), "equal": False, "reason": f"shape {arr.shape} against {(Z, Y, X)}"})
                continue
            ds = I._channel_dataset(f, int(c))
            for z in planes:
                if not 0 <= int(z) < Z:
                    continue
                a = ds[int(z), :Y, :X]
                diff = int((a != arr[int(z)]).sum())
                out.append({"channel": int(c), "plane": int(z), "equal": diff == 0, "n_different": diff})
    return out


def _stretch(img: np.ndarray) -> np.ndarray:
    lo, hi = np.percentile(img, [1.0, 99.9])
    return np.clip((img.astype(float) - lo) / max(hi - lo, 1e-9), 0, 1)


def write_mips(path, mips: dict, names: list[str], centre_yx, crop: int = 400) -> None:
    """One PNG: the downsampled field and a full-resolution crop per channel."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(mips)
    fig, axes = plt.subplots(2, n, figsize=(5 * n, 10))
    cy, cx = (int(v) for v in centre_yx)
    for k, (c, m) in enumerate(sorted(mips.items())):
        axes[0, k].imshow(_stretch(m[::4, ::4]), cmap="gray")
        axes[0, k].set_title(f"channel {c} ({names[c] if c < len(names) else '?'}), max over z")
        y0, x0 = max(cy - crop // 2, 0), max(cx - crop // 2, 0)
        axes[1, k].imshow(_stretch(m[y0:y0 + crop, x0:x0 + crop]), cmap="gray")
        axes[1, k].set_title(f"crop at y {y0}, x {x0} (full resolution)")
        for a in axes[:, k]:
            a.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--additions", default=str(REPO / "tests" / "golden" / "rad51_additions.json"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--tags", nargs="*", default=None)
    ap.add_argument("--planes", nargs="*", type=int, default=[0, 53, 105])
    ap.add_argument("--channel-map", default=str(DEFAULT_MAP))
    args = ap.parse_args(argv)
    import h5py

    from germquant.io.channel_map import ChannelMap

    rec = bs.load_additions(args.additions)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cmap = ChannelMap.from_yaml(args.channel_map)
    results, ok = {}, True
    for g in rec["gonads"]:
        if args.tags and g["tag"] not in args.tags:
            continue
        nd2, ims = g["nd2"]["path"], g["ims"]["path"]
        bs.refuse_sealed(nd2)
        bs.refuse_sealed(ims)
        names = list(g["nd2"].get("channel_names") or [])
        roles, flags = cmap.resolve(names)
        chans = {c: I.read_nd2_channel_zyx(nd2, c) for c in range(len(names))}
        planes = compare_planes(ims, chans, args.planes)
        with h5py.File(ims, "r") as f:
            smi = f["Scene8/Content/MegaSurfaces0/SurfaceModelInfo"][:]
            geom = I.read_geometry(f)
        com = np.c_[smi["CenterOfMassX"], smi["CenterOfMassY"], smi["CenterOfMassZ"]].astype(float)
        idx = I.imaris_to_voxel_index(np.median(com, axis=0)[None], geom)[0]
        mips = {c: a.max(axis=0) for c, a in chans.items()}
        write_mips(out / f"{g['tag']}__mips.png", mips, names, (idx[1], idx[0]))
        del chans
        foci_ok = roles.get("foci") is not None and names[roles["foci"]] == "545"
        r = {"nd2_channel_names": names, "channel_map": str(args.channel_map), "roles": roles, "map_flags": flags,
             "foci_on_545": bool(foci_ok), "planes": planes,
             "all_planes_equal": bool(planes and all(p.get("equal") for p in planes)),
             "mips_png": f"{g['tag']}__mips.png", "crop_centre_xy_voxel": [int(idx[0]), int(idx[1])],
             "by_eye": "to be filled in by the person who looks at the PNG (entry 9 (d)(1))"}
        results[g["tag"]] = r
        ok = ok and r["all_planes_equal"] and r["foci_on_545"]
        print(f"{g['tag']}: map roles {roles} (foci on 545: {foci_ok}); planes equal: {r['all_planes_equal']} "
              f"({sum(p.get('equal', False) for p in planes)}/{len(planes)}); -> {r['mips_png']}", flush=True)
    (out / "image_checks.json").write_text(json.dumps(results, indent=1, default=str) + "\n", encoding="utf-8")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
