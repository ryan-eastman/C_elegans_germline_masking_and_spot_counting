"""Parse experimental design out of the filename into tidy metadata columns.

Convention (from the N2 dataset): ``<date>_<genotype>_<treatment>_<SEX>_<rep>``
e.g. ``20251105_n2_nohs_HERM _001`` / ``20251021_n2_hs_MALE_2``.
"""
from __future__ import annotations

import re
from pathlib import Path

# herm = oocytes (6 SCs), male = spermatocytes (5 SCs — X is a univalent)
_SEX_TO_GERMCELL = {"herm": "oocyte", "male": "spermatocyte"}
_TREATMENT_CANON = {"hs": "heat", "nohs": "control"}


def parse_sample(path: str | Path, regex: str, defaults: dict | None = None) -> dict:
    defaults = dict(defaults or {})
    stem = Path(path).stem.strip()
    fields = {
        "image_id": stem,
        "genotype": defaults.get("genotype", "unknown"),
        "sex": defaults.get("sex", "unknown"),
        "treatment": defaults.get("treatment", "unknown"),
        "replicate": defaults.get("replicate", "NA"),
        "acquisition_date": "NA",
        "germ_cell": "unknown",
    }
    m = re.search(regex, stem, flags=re.IGNORECASE)
    if m:
        g = {k: v for k, v in m.groupdict().items() if v}
        if "date" in g:
            fields["acquisition_date"] = g["date"]
        if "genotype" in g:
            fields["genotype"] = g["genotype"].upper()
        if "treatment" in g:
            t = g["treatment"].lower()
            fields["treatment"] = _TREATMENT_CANON.get(t, t)
        if "sex" in g:
            s = g["sex"].lower()
            fields["sex"] = s
            fields["germ_cell"] = _SEX_TO_GERMCELL.get(s, "unknown")
        if "replicate" in g:
            fields["replicate"] = g["replicate"]
    return fields


def expected_sc_count(germ_cell: str) -> int | None:
    """WT pachytene: oocyte = 6 SCs (5 autosomal + XX); spermatocyte = 5 (X univalent)."""
    return {"oocyte": 6, "spermatocyte": 5}.get(germ_cell)
