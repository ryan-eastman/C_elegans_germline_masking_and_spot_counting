"""Acquisition metadata from the .nd2 text description: the per-channel camera, filter and laser record.

The Nikon description (``nd2.ND2File(path).text_info["description"]``) lists one "Plane #n:" block per
channel. Each block carries the channel name, the camera settings (exposure, bit depth, binning, camera
mode, conversion gain, readout speed), the modality, the X-Light filter wheels (emission, excitation,
dichroic) and the laser line that was on (``ExW:<line>; Power: <value>``). `parse_description` turns
every block into one record per channel name with the keys in `PLANE_FIELDS`; a line the description
does not carry gives None.

Consumers:
  * the ``acquisition`` stage (config key acquisition.read_exposures, off by default) maps the records
    onto the resolved channel roles by channel name and writes ``exp_ms_<role>`` and
    ``laser_pow_<role>`` exactly as the stage always has (ported from
    analysis/coloc/scripts/build_acquisition_metadata.py), then ``em_filter_<role>``,
    ``camera_mode_<role>`` and ``bit_depth_<role>``, and one ``acquisition_json`` column holding the
    complete per-channel record (channel names as keys) as compact JSON;
  * ``germquant info`` prints the per-channel table for an .nd2;
  * the expected-settings check (config section ``acquisition.expected``, code default absent = no
    check) compares each role's record with the settings a study expects and raises the advisory QC
    flags ``acq:<role>_em_filter_mismatch``, ``acq:<role>_exposure_above_max``,
    ``acq:<role>_power_above_max`` and ``acq:<role>_camera_mode_mismatch``. docs/QC.md has the
    worked example. The section is absent from every shipped config, so no default run can see them.
"""
from __future__ import annotations

import dataclasses
import json
import logging
import re
from pathlib import Path

log = logging.getLogger(__name__)

# the keys of one channel's record, in the order they are reported
PLANE_FIELDS: tuple[str, ...] = (
    "exposure_ms", "power", "em_filter", "ex_filter", "dichroic", "camera_mode", "bit_depth",
    "binning", "conversion_gain", "readout_speed", "modality",
)

# the keys an acquisition.expected entry may carry (one entry per channel role)
EXPECTED_KEYS: tuple[str, ...] = ("em_filter", "exposure_ms_max", "power_max", "camera_mode")

_PLANE_SPLIT = re.compile(r"Plane #\d+:")
# the plane name sits on its own line; "Camera Name: ..." in the header never starts a line with "Name:"
_NAME_LINE = re.compile(r"^[ \t]*Name:[ \t]*([^\r\n]+)", re.MULTILINE)
_NAME_LOOSE = re.compile(r"Name:\s*(\S+)")
_EXPOSURE = re.compile(r"Exposure:\s*([\d.]+)\s*(ms|s)\b")
_LASER_LINE = re.compile(r"ExW:([^;\r\n]+);\s*Power:\s*([\d.]+)")
_BIT_DEPTH = re.compile(r"Bit Depth:\s*(\d+)")
_STRING_FIELDS = {
    "camera_mode": re.compile(r"Camera Mode:[ \t]*([^\r\n]+)"),
    "binning": re.compile(r"Binning:[ \t]*(\S+)"),
    "conversion_gain": re.compile(r"Conversion Gain:[ \t]*([^\r\n]+)"),
    "readout_speed": re.compile(r"Readout Speed:[ \t]*([^\r\n]+)"),
    "modality": re.compile(r"Modality:[ \t]*([^\r\n]+)"),
}
# the X-Light wheels: "X-Light, FilterChanger(EM Wheel): 3 (595/31)"
_WHEELS = {
    "em_filter": re.compile(r"FilterChanger\(EM Wheel\):[ \t]*([^\r\n]+)"),
    "ex_filter": re.compile(r"FilterChanger\(EX Wheel\):[ \t]*([^\r\n]+)"),
    "dichroic": re.compile(r"FilterChanger\(DM Wheel\):[ \t]*([^\r\n]+)"),
}
_POSITION_LABEL = re.compile(r"^\d+\s*\((.*)\)$", re.DOTALL)


def _filter_label(raw: str) -> str | None:
    """The filter label of a wheel entry: ``"3 (595/31)"`` -> ``"595/31"`` (the wheel position is
    dropped; the label is what a person recognises and what acquisition.expected names). An entry with
    no position in front is returned as written."""
    raw = raw.strip()
    m = _POSITION_LABEL.match(raw)
    label = m.group(1).strip() if m else raw
    return label or None


def _string(rx: re.Pattern, block: str) -> str | None:
    m = rx.search(block)
    if not m:
        return None
    return m.group(1).strip() or None


def _plane_name(block: str) -> str | None:
    m = _NAME_LINE.search(block) or _NAME_LOOSE.search(block)
    return m.group(1).strip() if m and m.group(1).strip() else None


def _exposure_ms(block: str) -> float | None:
    m = _EXPOSURE.search(block)
    if not m:
        return None
    value = float(m.group(1))
    return value * 1000.0 if m.group(2) == "s" else value


def _power(block: str, name: str) -> float | None:
    """The power of the laser line named like the channel (``ExW:545; Power: 96.0``); None when no line
    carries the channel's name (NaN in laser_pow_<role> is the signal that the line and channel names
    disagree, the semantics every existing table was written with)."""
    m = re.search(rf"ExW:{re.escape(name)}; Power:\s*([\d.]+)", block)
    if m:
        return float(m.group(1))
    return None


def parse_description(desc: str) -> dict[str, dict]:
    """{channel name: record} from the nd2 description text, one record per "Plane #n:" block.

    Every record has every key of `PLANE_FIELDS`: ``exposure_ms`` and ``power`` are floats,
    ``bit_depth`` an int, the rest strings; a value the block does not carry is None. The filter keys
    hold the filter label (``"595/31"``), not the wheel position.
    """
    out: dict[str, dict] = {}
    for block in _PLANE_SPLIT.split(desc or "")[1:]:
        name = _plane_name(block)
        if not name:
            continue
        rec: dict = {"exposure_ms": _exposure_ms(block), "power": _power(block, name)}
        for key, rx in _WHEELS.items():
            raw = _string(rx, block)
            rec[key] = _filter_label(raw) if raw else None
        for key, rx in _STRING_FIELDS.items():
            rec[key] = _string(rx, block)
        bd = _BIT_DEPTH.search(block)
        rec["bit_depth"] = int(bd.group(1)) if bd else None
        out[name] = {k: rec.get(k) for k in PLANE_FIELDS}
    return out


def read_description(nd2_path: str | Path) -> str:
    import nd2

    with nd2.ND2File(str(nd2_path)) as h:
        return str(h.text_info.get("description", "") or "")


def role_records(per_name: dict[str, dict], channel_names: list[str], role_to_idx: dict) -> dict[str, dict]:
    """{role: record} for every resolved role: the record of the channel the role resolved to, matched
    by channel name (exact first, then a substring match either way, as the stage always did); an empty
    dict when the description has no plane for that channel."""
    out: dict[str, dict] = {}
    for role, idx in role_to_idx.items():
        if idx is None:
            continue
        name = str(channel_names[idx]) if idx < len(channel_names) else ""
        rec: dict = {}
        if name:
            rec = per_name.get(name) or next((v for k, v in per_name.items() if k and (k in name or name in k)), {})
        out[role] = dict(rec)
    return out


def _number(value) -> float:
    return float("nan") if value is None else value


def summary_fields(per_name: dict[str, dict], per_role: dict[str, dict]) -> dict:
    """The image_summary fields. ``exp_ms_<role>`` / ``laser_pow_<role>`` come first, for every role, as
    they always have (NaN when absent); then ``em_filter_<role>``, ``camera_mode_<role>`` (strings, None
    when absent, which writes an empty cell like every other missing stage value) and ``bit_depth_<role>``
    (NaN when absent); then ``acquisition_json``, the complete per-channel record as compact JSON."""
    fields: dict = {}
    for role, rec in per_role.items():
        fields[f"exp_ms_{role}"] = _number(rec.get("exposure_ms"))
        fields[f"laser_pow_{role}"] = _number(rec.get("power"))
    for role, rec in per_role.items():
        fields[f"em_filter_{role}"] = rec.get("em_filter")
        fields[f"camera_mode_{role}"] = rec.get("camera_mode")
        fields[f"bit_depth_{role}"] = _number(rec.get("bit_depth"))
    fields["acquisition_json"] = json.dumps(per_name, separators=(",", ":"))
    return fields


@dataclasses.dataclass
class AcquisitionRecord:
    per_channel: dict[str, dict]     # channel name -> record (PLANE_FIELDS), as the description lists them
    per_role: dict[str, dict]        # resolved role -> record ({} when the description has none for it)
    fields: dict                     # the image_summary fields (summary_fields)


def read_acquisition(nd2_path: str | Path, channel_names: list[str], role_to_idx: dict) -> AcquisitionRecord:
    """Read and parse the .nd2 description and map it onto the resolved roles. An unreadable description
    is logged and gives an empty record (every field NA): a metadata read must never fail the run."""
    try:
        per_name = parse_description(read_description(nd2_path))
    except Exception as e:  # noqa: BLE001 - a metadata read must never fail the run
        log.warning("acquisition metadata unreadable (%s); fields left NA", e)
        per_name = {}
    per_role = role_records(per_name, channel_names, role_to_idx)
    return AcquisitionRecord(per_channel=per_name, per_role=per_role, fields=summary_fields(per_name, per_role))


def acquisition_fields(nd2_path: str | Path, channel_names: list[str], role_to_idx: dict) -> dict:
    """image_summary fields for every resolved role (see `summary_fields`)."""
    return read_acquisition(nd2_path, channel_names, role_to_idx).fields


# ---- expected settings -------------------------------------------------------------------------------

def parse_expected(section, roles: set[str] | None = None) -> dict[str, dict]:
    """Validate an ``acquisition.expected`` section: {role: {em_filter: str | [str], exposure_ms_max:
    number, power_max: number, camera_mode: str | [str]}}. Returns {} for an absent section and
    {role: {em_filter: [str], camera_mode: [str], exposure_ms_max: float, power_max: float}} (only the
    keys given) otherwise. A role the channel map does not define, an unknown key or a value of the
    wrong type raises ValueError, so a typo fails the run before any pixel is read instead of silently
    checking nothing."""
    if section is None:
        return {}
    if not isinstance(section, dict):
        raise ValueError("acquisition.expected must be a mapping of channel role -> expected settings")
    out: dict[str, dict] = {}
    for role, spec in section.items():
        role = str(role)
        if roles is not None and role not in roles:
            raise ValueError(f"acquisition.expected: {role!r} is not a role of the channel map "
                             f"(roles: {sorted(roles)})")
        if not isinstance(spec, dict) or not spec:
            raise ValueError(f"acquisition.expected.{role} must be a mapping with at least one of {EXPECTED_KEYS}")
        unknown = sorted(set(spec) - set(EXPECTED_KEYS))
        if unknown:
            raise ValueError(f"acquisition.expected.{role}: unknown key(s) {unknown}; allowed: {list(EXPECTED_KEYS)}")
        clean: dict = {}
        for key in ("em_filter", "camera_mode"):
            if key in spec:
                values = spec[key] if isinstance(spec[key], (list, tuple)) else [spec[key]]
                if not values or any(v is None or isinstance(v, (bool, dict, list)) for v in values):
                    raise ValueError(f"acquisition.expected.{role}.{key} must be a string or a list of strings")
                clean[key] = [str(v) for v in values]
        for key in ("exposure_ms_max", "power_max"):
            if key in spec:
                v = spec[key]
                if isinstance(v, bool) or not isinstance(v, (int, float)):
                    raise ValueError(f"acquisition.expected.{role}.{key} must be a number")
                clean[key] = float(v)
        out[role] = clean
    return out


def expected_settings(cfg) -> dict[str, dict]:
    """The validated ``acquisition.expected`` section of a config ({} when absent, the code default)."""
    return parse_expected(cfg.get("acquisition.expected"), roles=set(cfg.channel_map.roles))


def _normalise(text) -> str:
    return re.sub(r"\s+", "", str(text)).lower()


def _matches(observed, accepted: list[str]) -> bool:
    """String match after removing whitespace, ignoring case. An absent value matches nothing: the
    expectation could not be confirmed, which is reported like a mismatch."""
    if observed is None:
        return False
    o = _normalise(observed)
    return any(o == _normalise(a) for a in accepted)


def expected_flags(per_role: dict[str, dict], expected: dict[str, dict]) -> list[str]:
    """The QC flags of the expected-settings check: ``acq:<role>_em_filter_mismatch``,
    ``acq:<role>_camera_mode_mismatch`` (string match after stripping whitespace, case-insensitive; an
    absent value is a mismatch), ``acq:<role>_exposure_above_max`` and ``acq:<role>_power_above_max``
    (silent when the value is absent: nothing is above a maximum). A role that did not resolve in this
    image is not checked. The flags are advisory: qc_pass does not change."""
    flags: list[str] = []
    for role, spec in expected.items():
        rec = per_role.get(role)
        if rec is None:
            continue
        if "em_filter" in spec and not _matches(rec.get("em_filter"), spec["em_filter"]):
            flags.append(f"acq:{role}_em_filter_mismatch")
        exposure = rec.get("exposure_ms")
        if "exposure_ms_max" in spec and exposure is not None and exposure > spec["exposure_ms_max"]:
            flags.append(f"acq:{role}_exposure_above_max")
        power = rec.get("power")
        if "power_max" in spec and power is not None and power > spec["power_max"]:
            flags.append(f"acq:{role}_power_above_max")
        if "camera_mode" in spec and not _matches(rec.get("camera_mode"), spec["camera_mode"]):
            flags.append(f"acq:{role}_camera_mode_mismatch")
    return flags


def format_table(per_name: dict[str, dict]) -> str:
    """The per-channel table ``germquant info`` prints: exposure, laser power, emission filter and camera
    mode per channel, NA where the description has no value."""
    def num(v):
        return "NA" if v is None else f"{v:g}"

    lines = [f"  {'channel':<10} {'exposure_ms':>11} {'power':>6}  {'em_filter':<12} camera_mode"]
    for name, rec in per_name.items():
        lines.append(f"  {name:<10} {num(rec.get('exposure_ms')):>11} {num(rec.get('power')):>6}  "
                     f"{(rec.get('em_filter') or 'NA'):<12} {rec.get('camera_mode') or 'NA'}")
    return "\n".join(lines)
