"""Strict checks of the opt-in recipe keys of tranche 5 package P2a (docs/DESIGN_system_packs.md sections 3.10,
3.12, 4.2 and 4.7; docs/PLANAR.md, "Recipe keys for packs").

process_image calls check(cfg) once per image, right after the read and before any stage runs, so a malformed
value stops the image before any pixel is processed and before any file is written. This closes two ways a
recipe typo could change the numbers without a word:

* a loose reading that switches a feature on when the recipe meant it off: bool("false") is True, so a quoted
  "false" for require_channel, sharpen or remove_hot_pixels would turn the behaviour on;
* an error raised inside the Cellpose backend, which method auto swallows and replaces with the classical
  watershed (segment_nuclei, segment_planar), so a bad object_diameter_um would silently swap Cellpose masks
  for watershed masks. Design 4.4: refuse, never skip silently.

Each check runs only when its key is set to something other than null. No shipped config or profile sets any
of these keys, so a default run checks nothing and changes nothing.

    segmentation.nuclei.cellpose                     null or a mapping
    segmentation.nuclei.cellpose.object_diameter_um  null or a positive number (micrometres)
    segmentation.nuclei.require_channel              null, true or false
    segmentation.nuclei.labels_from                  null or a path template (P2c, docs/LABELS.md)
    segmentation.nuclei.labels_frame                 null or a whole number, 0 or more (P2c)
    qc.min_nuclei                                    null or a whole number, 0 or more
    spots.spot_radius_z_um                           null or a positive number (micrometres)
    spots.remove_hot_pixels, spots.sharpen           null, true or false

The three spots keys are checked in the spots: section and in every spots.instances entry, which can set them
for its own instance. The accepted values are the types YAML gives: true or false, whole numbers and decimals.
A quoted value ("false", "4") is text, not a boolean or a number, and is refused. So is a numpy scalar other
than a float64, which a config could only hold through cfg.set and which the run manifest could not record
(Config.hash and as_dict serialise the config as JSON). The projection keys (segmentation.nuclei.project,
project_plane, min_area_um2) are checked by germquant.planar.prepare_projection in projected runs, where they
are read.

The module imports only the standard library (no numpy, backend, torch, Cellpose or SpotMAX), so a recipe
linter can use it.
"""
from __future__ import annotations

import math
import string

# the optional SpotMAX keys, by kind (spots.detect.SPOTMAX_OPTIONAL_KEYS holds the same three names; a test
# keeps the two in step without importing the spots package here)
SPOT_LENGTH_KEYS = ("spot_radius_z_um",)
SPOT_FLAG_KEYS = ("remove_hot_pixels", "sharpen")

CHECKED_KEYS = (
    "segmentation.nuclei.cellpose",
    "segmentation.nuclei.cellpose.object_diameter_um",
    "segmentation.nuclei.require_channel",
    "qc.min_nuclei",
    *(f"spots.{k}" for k in SPOT_LENGTH_KEYS + SPOT_FLAG_KEYS),
    "segmentation.nuclei.labels_from", "segmentation.nuclei.labels_frame",
)
# the fields of a labels_from template (germquant.io.labels formats it)
LABELS_TEMPLATE_FIELDS = ("image_dir", "image_stem", "image_id")


def _is_int(value) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)       # YAML true is an int in Python


def _quoted(value) -> str:
    return " (a quoted value is text: write it without quotes)" if isinstance(value, str) else ""


def flag_problem(key: str, value) -> str | None:
    """None when ``value`` is null, true or false; otherwise the message naming ``key``."""
    if value is None or isinstance(value, bool):
        return None
    return f"{key} must be true or false (or null for the code default); got {value!r}{_quoted(value)}"


def positive_length_problem(key: str, value) -> str | None:
    """None when ``value`` is null or a finite number above 0 (a whole number or a decimal, not a boolean);
    otherwise the message naming ``key``."""
    if value is None:
        return None
    if _is_int(value) or isinstance(value, float):
        try:
            v = float(value)
        except OverflowError:                 # a whole number too large for a float
            v = math.inf
        if math.isfinite(v) and v > 0:
            return None
    return f"{key} must be a positive number of micrometres (or null); got {value!r}{_quoted(value)}"


def count_problem(key: str, value) -> str | None:
    """None when ``value`` is null or a whole number of 0 or more (not a boolean, not 2.0); otherwise the
    message naming ``key``."""
    if value is None or (_is_int(value) and value >= 0):
        return None
    return f"{key} must be a whole number, 0 or more (or null for the code default); got {value!r}{_quoted(value)}"


def labels_template_problem(key: str, value) -> str | None:
    """None when ``value`` is null or a path template: a non-empty string whose fields are only {image_dir},
    {image_stem} and {image_id}, with no format spec or conversion ({{ and }} are literal braces). Otherwise the
    message naming ``key``."""
    if value is None:
        return None
    want = "a path template with the fields {image_dir}, {image_stem} and {image_id} (or null)"
    if not isinstance(value, str) or not value.strip():
        return f"{key} must be {want}; got {value!r}"
    try:
        parts = list(string.Formatter().parse(value))
    except ValueError as e:                    # an unbalanced brace
        return f"{key} must be {want}; {value!r} is not a template ({e})"
    for _text, field, spec, conversion in parts:
        if field is None:
            continue
        if field not in LABELS_TEMPLATE_FIELDS or spec or conversion:
            shown = "{" + field + (f"!{conversion}" if conversion else "") + (f":{spec}" if spec else "") + "}"
            return f"{key} must be {want}; got the field {shown} in {value!r}"
    return None


def object_diameter_problem(cfg_section) -> str | None:
    """The problem with segmentation.nuclei.cellpose.object_diameter_um in a segmentation.nuclei section, or
    None (unset or valid)."""
    sub = (cfg_section or {}).get("cellpose")
    if sub is None:
        return None
    if not isinstance(sub, dict):
        return f"segmentation.nuclei.cellpose must be a mapping of Cellpose keys (or null); got {sub!r}"
    return positive_length_problem("segmentation.nuclei.cellpose.object_diameter_um", sub.get("object_diameter_um"))


def object_diameter_um(cfg_section) -> float | None:
    """segmentation.nuclei.cellpose.object_diameter_um from a segmentation.nuclei section: None when unset, the
    value in micrometres when valid, ValueError naming the key otherwise (segment.backends.object_diameter_px
    turns it into Cellpose's diameter in pixels)."""
    problem = object_diameter_problem(cfg_section)
    if problem:
        raise ValueError(problem)
    value = ((cfg_section or {}).get("cellpose") or {}).get("object_diameter_um")
    return None if value is None else float(value)


def problems(cfg) -> list[str]:
    """Every malformed opt-in key of a config (``cfg.get(dotted, default)``), one message per key; empty when
    none is set or every set key is valid."""
    found: list[str] = []
    seg = cfg.get("segmentation.nuclei", {}) or {}
    if isinstance(seg, dict):
        found.append(object_diameter_problem(seg))
        found.append(flag_problem("segmentation.nuclei.require_channel", seg.get("require_channel")))
        found.append(labels_template_problem("segmentation.nuclei.labels_from", seg.get("labels_from")))
        found.append(count_problem("segmentation.nuclei.labels_frame", seg.get("labels_frame")))
    found.append(count_problem("qc.min_nuclei", cfg.get("qc.min_nuclei")))
    spots = cfg.get("spots", {}) or {}
    sections = [("spots", spots)] if isinstance(spots, dict) else []
    instances = (spots.get("instances") or []) if isinstance(spots, dict) else []
    # an entry that is not a mapping is refused by pipeline.spot_instances, with its own message
    sections += [(f"spots.instances[{i}]", inst) for i, inst in enumerate(instances) if isinstance(inst, dict)]
    for where, section in sections:
        found += [positive_length_problem(f"{where}.{k}", section.get(k)) for k in SPOT_LENGTH_KEYS]
        found += [flag_problem(f"{where}.{k}", section.get(k)) for k in SPOT_FLAG_KEYS]
    return [p for p in found if p]


def check(cfg) -> None:
    """Raise ValueError listing every malformed opt-in key (see the module docstring); return quietly when
    there is none, which is every run that sets none of these keys."""
    found = problems(cfg)
    if found:
        raise ValueError("recipe refused before any stage ran: " + "; ".join(found))
