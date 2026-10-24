"""The resolved ``sc_straighten`` parameter section (docs/SC_STRAIGHTEN.md, "Config").

Precedence: the code defaults below < the ``sc_straighten`` section of ``--config`` / ``--profile`` <
``--set KEY=VALUE``. No shipped config or profile names the section, so no config_hash moves. Axis order is
(z, lateral) wherever two values appear and (z, y, x) for three. ``params_hash`` (12 hex) is the sha256 of the
canonical JSON of the resolved section; it enters the image sub-hash and every output row.
"""
from __future__ import annotations

import copy
import hashlib
import json
from typing import Any

DEFAULTS: dict[str, Any] = {
    "enabled": False,               # reserved for a later registry stage; the command ignores it
    "nuclei": "germline",           # germline | pachytene | all
    "territory_reach_um": 0.4,
    "crop_pad_um": [0.6, 0.5, 0.5],
    "profile": {"roles": ["foci"], "ds_um": 0.05, "core_um": [0.25, 0.12], "integrate_um": [0.7, 0.35],
                "background_pct": 15},
    "noise": {"planes": 12, "bins": 30, "qmax": 0.85, "poor_rel_rms": 0.15, "nucleus_sd_range": [0.7, 1.4],
              "var_floor": 1.0},
    "fields": {"sigma_um": [0.2, 0.108, 0.108], "null_quantile": 0.999, "null_shape": [40, 96, 96],
               "null_seed": 0, "crest_offset_max": 0.6},
    "follow": {"seed_min": 2.5, "lo": 1.5, "step_um": 0.1, "momentum": 0.5, "max_correction_um": 0.15,
               "lineness_max": 0.5, "weak_run_um": 0.6, "claim_vox": [1, 2, 2], "contact_cos": 0.85,
               "loop_back_um": 0.6, "min_length_um": 0.5, "max_length_um": 30.0, "stack_edge_um": 0.4,
               "contact_um": [0.3, 0.15], "max_turn_deg": 30.0, "contact_free_um": [0.6, 0.4]},
    "link": {"max_gap_um": 1.0, "max_gap_contact_um": 1.5, "max_angle_deg": 50.0, "max_skew_um": None,
             "gap_cost_per_um": 40.0, "skew_cost_per_um": 100.0, "unlinked_cost": 60.0, "tip_um": 0.3,
             "dir_window_um": [0.3, 1.2]},
    "dedupe": {"lateral_um": 0.15, "axial_um": 0.3, "min_run_um": 0.4, "parallel_cos": 0.9,
               "min_remainder_um": 0.5, "ghost_max_um": 1.5, "ghost_tz_min": 0.7, "ghost_reach_um": [1.0, 0.35],
               "shadow_frac": 0.8},
    "refine": {"ds_um": 0.1, "window_um": [0.30, 0.15], "grid_um": [0.6, 0.4], "max_shift_um": 0.3,
               "iterations": 4, "knot_um": 0.5, "sigma_loc_um": [0.06, 0.03, 0.03], "lambda": 0.1,
               "lambda_imported": 0.01, "se_curve_um": 0.10},
    "ends": {"anchor_um": 0.6, "inside_um": 0.6, "outside_um": 1.4, "step_um": 0.05, "max_recentre_um": 0.2,
             "sigma_bounds": [0.6, 1.6], "fit_min_amp_se": 5.0, "fit_max_rel_rms": 0.25, "max_extension_um": 0.5},
    "psf": {"mode": "auto", "fallback_sigma_um": [0.365, 0.15], "min_stretches": 20, "flat_tz_max": 0.3,
            "isolation_um": 1.0, "min_snr": 10.0, "min_stretch_um": 1.0,
            # v2 (protocol changelog entries 12 to 15): the foci estimate of the image PSF ("off" gives v1); withdrawn
            # by entry 16 (its confirmation failed at sigma_z 0.33 um), so the default is off
            "foci": "off", "foci_min": 15, "foci_window_um": [1.2, 0.54], "foci_isolation_um": 1.0,
            "foci_max_sigma_xy_um": 0.25, "foci_max_ratio": 1.5, "foci_max_rel_rms": 0.2, "foci_min_snr": 8.0},
    "gaps": {"q": 0.6, "window_frac": 0.9, "min_um": 0.25, "z": 3.0, "smooth_um": 0.1, "ring_um": 0.45,
             "tail_sigmas": 1.0, "plateau_um": 0.4},
    "contamination_um": [0.7, 0.35],
    "strips": {"write": True, "ds_um": 0.05, "du_um": 0.05, "dv_um": 0.1, "half_um": [1.0, 0.6],
               "extend_um": 0.5, "flat_v_um": 0.3, "axial_u_um": 0.15, "roles": ["central_element", "foci"]},
    "samples": {"step_um": 0.1},
    "foci": {"tables": "auto", "max_um": [0.6, 0.3], "ambiguity_ratio": 1.25},
    "classes": {"min_complete_um": 1.0, "max_link_cost": 60.0, "max_contaminated": 0.3,
                "merged_multiplicity": 1.7, "merged_run_um": 0.5, "min_frac_in_label": 0.8},
    "qc": {"low_snr_review": 4.0, "low_snr_fail": 2.5, "steep_frac_review": 0.5, "long_um": 10.0,
           "kink_deg": 60.0, "kink_um": 0.3, "fit_rms_review_um": 0.12, "gap_long_um": 0.6,
           "stack_edge_um": 0.4, "merged_label_ratio": 1.6, "double_width_ratio": 1.6,
           "near_other_label_fail": 0.5, "saturated_frac": 0.01, "sc_channel_suspect_frac": 0.2,
           "merged_label_window_um": 15.0},
    "nucleus": {"present_min_um": 1.0, "present_min_snr": 3.0, "expected": {"oocyte": 6, "spermatocyte": 5},
                "male_x_rule": False},
    "budget_s": 60,
    "workers": 1,
    "overview": True,
}

# The value of every v2 switch that gives the behaviour of v1 (commit 5cb3cb6), one per change of the v2 round
# (protocol changelog entries 12 on); tests/test_sc_straighten_v2_compat.py runs them all against v1's outputs.
V1_VALUES: dict[str, Any] = {
    "psf.foci": "off",
    "psf.min_stretches": 20,
}


def v1_sets() -> list[str]:
    """``--set`` assignments that switch every v2 change back to v1."""
    return [f"{k}={json.dumps(v)}" for k, v in V1_VALUES.items()]


def switch_on(value) -> bool:
    """A v2 switch value: "on" (or YAML true) is on; "off" (YAML reads an unquoted off as false) is off."""
    return value is True or str(value).strip().lower() in ("on", "true", "yes")


# keys a development sweep may tune (DESIGN section 8); everything else is fixed unless the protocol changelog
# says otherwise
TUNABLE_PREFIXES = ("follow.seed_min", "follow.lo", "link.", "refine.knot_um", "refine.lambda",
                    "refine.sigma_loc_um", "gaps.", "classes.")


class Params:
    """A resolved, read-only view of the section: ``p["follow.seed_min"]`` or ``p.get("follow.seed_min")``."""

    def __init__(self, data: dict):
        self._d = copy.deepcopy(data)

    def get(self, dotted: str, default: Any = None) -> Any:
        node: Any = self._d
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node

    def __getitem__(self, dotted: str) -> Any:
        sentinel = object()
        v = self.get(dotted, sentinel)
        if v is sentinel:
            raise KeyError(dotted)
        return v

    def as_dict(self) -> dict:
        return copy.deepcopy(self._d)

    @property
    def hash(self) -> str:
        return params_hash(self._d)

    def __getstate__(self):
        return self._d

    def __setstate__(self, state):
        self._d = state


def canonical_json(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True, default=str)


def params_hash(d: dict) -> str:
    return hashlib.sha256(canonical_json(d).encode()).hexdigest()[:12]


def _merge(base: dict, overlay: dict, path: str = "") -> dict:
    out = copy.deepcopy(base)
    for k, v in (overlay or {}).items():
        key = f"{path}{k}"
        if k not in out:
            raise KeyError(f"unknown sc_straighten key {key!r} (docs/SC_STRAIGHTEN.md lists the keys)")
        if isinstance(out[k], dict):
            if not isinstance(v, dict):
                raise ValueError(f"sc_straighten.{key} is a section; got {v!r}")
            out[k] = _merge(out[k], v, key + ".")
        else:
            out[k] = copy.deepcopy(v)
    return out


def parse_value(text: str) -> Any:
    """A ``--set`` value: YAML scalar or flow sequence (``3.5``, ``true``, ``[1, 2, 2]``, ``null``)."""
    import yaml

    return yaml.safe_load(text)


def apply_set(d: dict, assignment: str) -> dict:
    """One ``KEY=VALUE`` (the key relative to the section, for example ``follow.seed_min=3.5``)."""
    if "=" not in assignment:
        raise ValueError(f"--set wants KEY=VALUE, got {assignment!r}")
    key, text = assignment.split("=", 1)
    key = key.strip()
    if key.startswith("sc_straighten."):
        key = key[len("sc_straighten."):]
    parts = key.split(".")
    node = d
    for part in parts[:-1]:
        if part not in node or not isinstance(node[part], dict):
            raise KeyError(f"unknown sc_straighten key {key!r}")
        node = node[part]
    if parts[-1] not in node:
        raise KeyError(f"unknown sc_straighten key {key!r}")
    if isinstance(node[parts[-1]], dict):
        raise ValueError(f"sc_straighten.{key} is a section, not a value")
    node[parts[-1]] = parse_value(text)
    return d


def resolve(section: dict | None = None, sets: list[str] | tuple[str, ...] = ()) -> Params:
    """Code defaults, then the config's ``sc_straighten`` section, then each ``--set``."""
    d = _merge(DEFAULTS, section or {})
    for s in sets or ():
        d = apply_set(d, s)
    return Params(d)


def from_config(cfg, sets=()) -> Params:
    """The section of a loaded germquant Config (``cfg.get("sc_straighten")``), or the defaults for None."""
    section = cfg.get("sc_straighten", None) if cfg is not None else None
    return resolve(section if isinstance(section, dict) else None, sets)
