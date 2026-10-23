"""The QC flag vocabulary as data, and the rules that raise the flags.

Each entry of ``CODES`` is (code, scope, level, threshold key, description). Levels: ``info`` (a fact),
``review`` (look before trusting the number) and ``fail`` (the number is not usable for its purpose).
``qc_status`` is ``fail`` when any fail flag is present, else ``review`` when any review flag is present, else
``pass``. docs/SC_STRAIGHTEN.md carries the same table; tests/test_sc_straighten_qc.py keeps the two equal and
gives every code one positive and one negative case.
"""
from __future__ import annotations

import math

CODES: tuple[tuple[str, str, str, str, str], ...] = (
    # strand scope
    ("fragment", "strand", "info", "classes.min_complete_um", "length_um below the minimum of a whole strand"),
    ("low_snr", "strand", "review", "qc.low_snr_review", "median core SNR below the review level"),
    ("low_snr_fail", "strand", "fail", "qc.low_snr_fail", "median core SNR below the fail level (replaces low_snr)"),
    ("end_unfit", "strand", "review", "ends.fit_min_amp_se", "a free end whose half-maximum fit failed (the traced end is kept)"),
    ("contact_end", "strand", "info", "", "an end where the strand meets another strand"),
    ("truncated_z", "strand", "fail", "qc.stack_edge_um", "an end at the stack edge, or a sample within stack_edge_um of the first or last plane"),
    ("costly_link", "strand", "review", "classes.max_link_cost", "a link of the end matching costs more than max_link_cost"),
    ("contaminated", "strand", "review", "classes.max_contaminated", "more than max_contaminated of the samples lie near another strand"),
    ("merged_suspect", "strand", "review", "classes.merged_multiplicity", "integrated signal at or above merged_multiplicity x the image's complete-strand median over at least merged_run_um"),
    ("double_width", "strand", "review", "qc.double_width_ratio", "lateral or axial FWHM above double_width_ratio x the image's complete-strand median"),
    ("exits_label", "strand", "review", "classes.min_frac_in_label", "less than min_frac_in_label of the strand inside its nucleus label"),
    ("neighbour_nucleus", "strand", "fail", "qc.near_other_label_fail", "more than near_other_label_fail of the strand near another nucleus label"),
    ("steep_z", "strand", "review", "qc.steep_frac_review", "more than steep_frac_review of the strand steeper than 60 degrees to the image plane"),
    ("long", "strand", "review", "qc.long_um", "longer than long_um: two chromosomes joined end to end"),
    ("kink", "strand", "review", "qc.kink_deg", "a turn above kink_deg within kink_um"),
    ("fit_residual", "strand", "review", "qc.fit_rms_review_um", "centreline fit residual RMS above fit_rms_review_um"),
    ("gap_long", "strand", "review", "qc.gap_long_um", "a gap of gap_long_um or more"),
    ("saturated", "strand", "review", "qc.saturated_frac", "more than saturated_frac of the core samples at the detector maximum"),
    ("imported", "strand", "info", "", "geometry imported from a file (Imaris Filaments, SWC or CSV)"),
    ("verified", "strand", "info", "", "imported geometry marked verified"),
    ("error", "strand", "fail", "", "an exception while measuring the strand (error:<Type>)"),
    # nucleus scope
    ("label_merged_suspect", "nucleus", "review", "qc.merged_label_ratio", "label volume above merged_label_ratio x the median of traced nuclei"),
    ("truncated_z", "nucleus", "fail", "", "the label touches the first or last plane (totals unusable)"),
    ("edge_xy", "nucleus", "review", "", "the label touches the image border"),
    ("count_low", "nucleus", "info", "nucleus.expected", "fewer main strands than expected"),
    ("count_high", "nucleus", "info", "nucleus.expected", "more main strands than expected"),
    ("noise_mismatch", "nucleus", "review", "noise.nucleus_sd_range", "the transformed noise SD of the territory outside nucleus_sd_range"),
    ("budget_exceeded", "nucleus", "review", "budget_s", "tracing stopped at the work budget (strands so far kept)"),
    ("error", "nucleus", "fail", "", "an exception in the nucleus (error:<Type>; status failed, no zero row)"),
    # image scope
    ("sc_straighten:uncalibrated", "image", "info", "", "always in v1, until the validation protocol has scored the method"),
    ("noise_fit_poor", "image", "review", "noise.poor_rel_rms", "the noise model fits its bins with a relative RMS above poor_rel_rms (or fell back)"),
    ("psf_fallback", "image", "review", "psf.min_stretches", "fewer than min_stretches flat isolated bright stretches: the fallback PSF was used"),
    ("sc_channel_suspect", "image", "review", "qc.sc_channel_suspect_frac", "among nuclei with SC present, fewer than sc_channel_suspect_frac have a strand of 2 um or more"),
    ("foci_uncalibrated", "image", "info", "", "a spot table was used (spot detection is not calibrated for SC positions)"),
    ("geometry_from_manifest", "image", "info", "", "run geometry (xy_stride, z_range) from run_manifest.json: no done marker"),
    ("label_merged_suspect", "image", "info", "", "label_merged_suspect=N: the number of such nuclei"),
    ("nucleus_errors", "image", "info", "", "nucleus_errors=N: the number of failed nuclei"),
)

LEVEL = {(c, s): lvl for c, s, lvl, _k, _d in CODES}


def level(code: str, scope: str) -> str:
    base = code.split(":", 1)[0] if not code.startswith("sc_straighten:") else code
    base = base.split("=", 1)[0]
    return LEVEL.get((base, scope), "info")


def status(flags: list[str], scope: str) -> str:
    lv = {level(f, scope) for f in flags}
    if "fail" in lv:
        return "fail"
    if "review" in lv:
        return "review"
    return "pass"


def _f(x) -> float:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return v


def _gt(x, t) -> bool:
    v = _f(x)
    return math.isfinite(v) and v > float(t)


def _lt(x, t) -> bool:
    v = _f(x)
    return math.isfinite(v) and v < float(t)


def strand_flags(m: dict, p, image_ctx: dict | None = None) -> list[str]:
    """Flags of one strand from its measurement row ``m`` (the strands table fields). ``image_ctx``: the
    image's complete-strand medians (``fwhm_lat_med``, ``fwhm_ax_med``) for double_width, and the merged run
    (``merged_run_um`` in m) for merged_suspect."""
    ctx = image_ctx or {}
    out: list[str] = []
    if _lt(m.get("length_um"), p["classes.min_complete_um"]):
        out.append("fragment")
    snr = _f(m.get("snr_core_median"))
    if math.isfinite(snr):
        if snr < float(p["qc.low_snr_fail"]):
            out.append("low_snr_fail")
        elif snr < float(p["qc.low_snr_review"]):
            out.append("low_snr")
    ends = (m.get("end0_class", ""), m.get("end1_class", ""))
    if "free_unfit" in ends:
        out.append("end_unfit")
    if "contact" in ends:
        out.append("contact_end")
    if "stack_edge" in ends or bool(m.get("near_stack_edge", False)):
        out.append("truncated_z")
    if _gt(m.get("link_cost_max"), p["classes.max_link_cost"]):
        out.append("costly_link")
    if _gt(m.get("frac_contaminated"), p["classes.max_contaminated"]):
        out.append("contaminated")
    if _f(m.get("merged_run_um", 0.0)) >= float(p["classes.merged_run_um"]) - 1e-9 and _f(m.get("merged_run_um", 0.0)) > 0:
        out.append("merged_suspect")
    r = float(p["qc.double_width_ratio"])
    if (_gt(m.get("fwhm_lat_um"), r * _f(ctx.get("fwhm_lat_med"))) or
            _gt(m.get("fwhm_ax_um"), r * _f(ctx.get("fwhm_ax_med")))):
        out.append("double_width")
    if _lt(m.get("frac_in_label"), p["classes.min_frac_in_label"]):
        out.append("exits_label")
    if _gt(m.get("frac_near_other_label"), p["qc.near_other_label_fail"]):
        out.append("neighbour_nucleus")
    if _gt(m.get("frac_steep60"), p["qc.steep_frac_review"]):
        out.append("steep_z")
    if _gt(m.get("length_um"), p["qc.long_um"]):
        out.append("long")
    if _gt(m.get("kink_max_deg"), p["qc.kink_deg"]):
        out.append("kink")
    if _gt(m.get("fit_rms_um"), p["qc.fit_rms_review_um"]):
        out.append("fit_residual")
    if _f(m.get("gap_max_um", 0.0)) >= float(p["qc.gap_long_um"]) - 1e-9:
        out.append("gap_long")
    if _gt(m.get("frac_saturated"), p["qc.saturated_frac"]):
        out.append("saturated")
    if m.get("geometry_source") == "imported":
        out.append("imported")
    if bool(m.get("verified", False)):
        out.append("verified")
    if m.get("error"):
        out.append(f"error:{m['error']}")
    return out


def nucleus_flags(n: dict, p, median_volume: float | None = None) -> list[str]:
    out: list[str] = []
    if median_volume and _gt(n.get("volume_um3"), float(p["qc.merged_label_ratio"]) * float(median_volume)):
        out.append("label_merged_suspect")
    if bool(n.get("touches_z", False)):
        out.append("truncated_z")
    if bool(n.get("touches_xy", False)):
        out.append("edge_xy")
    cs = n.get("count_status")
    cs = cs if isinstance(cs, str) else ""
    if cs == "low":
        out.append("count_low")
    elif cs == "high":
        out.append("count_high")
    sd = _f(n.get("noise_sd"))
    lo, hi = (float(v) for v in p["noise.nucleus_sd_range"])
    if math.isfinite(sd) and (sd < lo or sd > hi):
        out.append("noise_mismatch")
    if bool(n.get("budget_exceeded", False)):
        out.append("budget_exceeded")
    if n.get("error"):
        out.append(f"error:{n['error']}")
    return out


def image_flags(info: dict, p) -> list[str]:
    out = ["sc_straighten:uncalibrated"]
    rel = _f(info.get("noise_rel_rms"))
    if bool(info.get("noise_fallback")) or (math.isfinite(rel) and rel > float(p["noise.poor_rel_rms"])):
        out.append("noise_fit_poor")
    if info.get("psf_source") == "fallback":
        out.append("psf_fallback")
    n_present = int(info.get("n_present", 0) or 0)
    n_long = int(info.get("n_present_with_2um", 0) or 0)
    if n_present >= 10 and n_long < float(p["qc.sc_channel_suspect_frac"]) * n_present:
        out.append("sc_channel_suspect")
    if info.get("spot_tables"):
        out.append("foci_uncalibrated")
    if info.get("geometry_source") == "run_manifest":
        out.append("geometry_from_manifest")
    if int(info.get("n_label_merged", 0) or 0):
        out.append(f"label_merged_suspect={int(info['n_label_merged'])}")
    if int(info.get("n_failed", 0) or 0):
        out.append(f"nucleus_errors={int(info['n_failed'])}")
    return out


def docs_table() -> str:
    """The markdown table docs/SC_STRAIGHTEN.md carries (kept equal to CODES by a test)."""
    rows = ["| code | scope | level | threshold key | fires when |", "|---|---|---|---|---|"]
    for c, s, lvl, k, d in CODES:
        rows.append(f"| {c} | {s} | {lvl} | {k or '-'} | {d} |")
    return "\n".join(rows)
