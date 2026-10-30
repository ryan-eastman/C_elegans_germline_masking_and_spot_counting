"""Digests behind the sc_straighten v3 compatibility test (tests/test_sc_straighten_v3_compat.py) and its golden
(tests/golden/sc_straighten/v2_compat.json, made once with the package code of commit acbee94, the head of the v2
round, before any change of the v3 round).

Every declared column of every table is digested (the v2 columns included, not only those of v1), formatted as
the writer formats them, with the strips, the strip index, the number of strands and the record's class counts.
The cases run sc_straighten at its defaults and with the v2 punctum switches (protocol changelog entries 16 and
18), on registered families and on the development families, so a v3 change that leaks into a configuration
without its switch shows up here.
"""
from __future__ import annotations

import json
from pathlib import Path

import sc_v1_compat as C1

REGISTERED = tuple(("synth", f, 2001) for f in ("iso", "dense", "gaps", "foci", "pairs", "noise", "snr"))
DEVELOPMENT = (("dev", "puncta_pt", 2000), ("dev", "puncta_agg", 2001), ("dev", "dense_puncta", 2002),
               ("dev", "psf_check", 2000), ("dev", "real_dense_r2", 2001), ("dev", "real_dense_r2_brightness", 2002))
V2_PUNCTUM_LABEL = ("punctum.enabled=true",)                              # entry 18 (a label only)
V2_PUNCTUM_EXCLUDE = ("punctum.enabled=true", "punctum.exclude=true")     # entry 16
CASES = tuple((src, fam, seed, "defaults", ()) for src, fam, seed in REGISTERED + DEVELOPMENT) + \
    tuple((src, fam, seed, "v2_punctum_label", V2_PUNCTUM_LABEL) for src, fam, seed in
          (("dev", "puncta_pt", 2000), ("dev", "dense_puncta", 2002), ("synth", "gaps", 2001),
           ("dev", "psf_check", 2000))) + \
    tuple((src, fam, seed, "v2_punctum_exclude", V2_PUNCTUM_EXCLUDE) for src, fam, seed in
          (("dev", "dense_puncta", 2002), ("dev", "psf_check", 2000)))


def golden_path() -> Path:
    return Path(__file__).resolve().parent / "golden" / "sc_straighten" / "v2_compat.json"


def case_name(case) -> str:
    src, fam, seed, cfg, _sets = case
    return f"{fam}/{seed}/{cfg}"


def make_scene(src: str, fam: str, seed: int):
    from germquant.sc.straighten import synth, synth_dev

    return synth.make_scene(fam, seed) if src == "synth" else synth_dev.make_scene(fam, seed)


def run_scene(sc, sets=()):
    from germquant.sc.straighten import params as PRM
    from germquant.sc.straighten.pipeline import straighten_image

    P = PRM.resolve(None, list(sets))
    roles = {k: v for k, v in sc.channels.items() if k != "central_element"}
    spots = {"spots": sc.spots} if len(sc.spots) else None
    return straighten_image(sc.channels["central_element"], sc.labels, sc.spacing, sc.nuclei, params=P, roles=roles,
                            spot_tables=spots, image_id=f"{sc.family}_{sc.seed}",
                            wavelengths=sc.info.get("wavelengths"))


def case_digests(case, cols: dict) -> dict:
    src, fam, seed, _cfg, sets = case
    res = run_scene(make_scene(src, fam, seed), sets)
    out = C1.result_digests(res, cols)
    out["class_counts"] = json.dumps({k[6:]: int(v) for k, v in sorted(res.record["counts"].items())
                                      if k.startswith("class_")}, sort_keys=True)
    return out
