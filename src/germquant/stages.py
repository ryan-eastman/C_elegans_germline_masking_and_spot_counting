"""The stage registry ("registry-lite", docs/ROADMAP_modular_pipeline.md step 1).

`STAGES` is the single, explicitly ORDERED list of what `pipeline.process_image` runs. It is data, not
a resolver: the pipeline body still calls the stages in this order by hand, and nothing re-sorts them.
Each entry says how the stage is switched on or off, which channel roles it needs, and whether the CLI
grows a generated ``--no-<stage>`` switch for it. New quantifications (envelope, staging, sc_trace,
partition, ...) are added here and called from process_image through `run_stage`, never as inline
bodies.

Design rules (they protect the validated RAD-51 and ccw77 numbers):
  * every switch is read with ``cfg.get(key, code_default)``: the config files are not edited, so their
    whole-text ``config_hash`` keeps reproducing;
  * a stage failure is recorded and flagged as ``<stage>:FAILED_<Exception>`` and the run continues
    (``fatal=True`` stages re-raise, exactly as they did before the registry existed);
  * ``--no-coloc`` keeps setting only ``coloc.enabled`` (an override enters the hash); the granule stage
    defaults to ``coloc.enabled`` so the legacy flag still yields no granule outputs;
  * failure contract for the split granule / coloc pair: a granule failure flags ``granule:FAILED_X`` and
    coloc is skipped ("granule stage failed"); a coloc failure flags ``coloc:FAILED_X`` and KEEPS the
    granule products (granules table with nucleus_id, ``__granules_labels.tif``, n_granules on nuclei)
    while image_summary gets no coloc fields, so its shape equals a run without coloc. Consumers read
    ``<image_id>__stages.json`` or the FAILED flag, never the presence of n_granules.
"""
from __future__ import annotations

import logging
import time
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class Stage:
    name: str
    config_key: str | None            # dotted key that switches the stage; None = always on
    code_default: bool                # value assumed when the key is absent from the config file
    required_roles: tuple[str, ...]   # channel-map roles the stage needs (empty = none)
    cli_switchable: bool              # generate --no-<name> on `germquant run` / `batch`
    help: str
    default_from: str | None = None   # another stage's key whose value is the default when ours is absent
    legacy_cli_help: str | None = None
    config_sections: tuple[str, ...] = ()   # top-level config sections whose text feeds the stage sub-hash


STAGES: tuple[Stage, ...] = (
    Stage("read", None, True, ("dna",), False, "read the .nd2, resolve channel roles, parse the sample",
          config_sections=("io", "metadata")),
    Stage("acquisition", "acquisition.read_exposures", False, (), True,
          "per-channel exposure and laser power from the .nd2 description -> image_summary",
          config_sections=("acquisition",)),
    Stage("segment", None, True, (), False, "Cellpose 3D nucleus segmentation (falls back to classical)",
          config_sections=("segmentation",)),
    Stage("measure", None, True, (), False, "3D regionprops + per-role mean intensities per nucleus"),
    Stage("signal_qc", "signal_qc.enabled", True, (), True,
          "per-channel background, noise, SNR, saturation and exposure-normalized levels (image_summary sig_*) "
          "and pairwise channel cross-talk in nuclei / perinuclear ring / cytoplasm (channel_crosstalk table); "
          "its flags are off unless signal_qc.flags.enabled", config_sections=("signal_qc",)),
    Stage("germline", "germline.enabled", True, ("central_element",), True,
          "drop nuclei outside the gonad (SYP-seeded connected component)", config_sections=("germline",)),
    Stage("envelope", "envelope.enabled", False, ("lamin",), True,
          "lamin-envelope nucleus masks (seeded watershed), ring test, territories", config_sections=("envelope",)),
    Stage("audit", "audit.enabled", False, ("lamin",), True,
          "mask audit: lamin-only nuclei missed by the labels, labels with no envelope (+ overlay)",
          config_sections=("audit", "envelope")),
    Stage("axis", "axis.enabled", True, (), True,
          "principal-curve germline axis -> distal->proximal position per nucleus", config_sections=("axis",)),
    Stage("staging", "staging.enabled", False, (), True,
          "pachytene zones (early / mid / late) from the hand-traced axis in staging.traces_file",
          config_sections=("staging",)),
    Stage("spots", "spots.enabled", True, ("foci",), True,
          "SpotMAX spot counting in the `foci` channel (RAD-51); further instances (COSA-1 crossover foci) "
          "via spots.instances, each on its own channel role with its own table and column",
          legacy_cli_help="segmentation only: skip RAD-51/SpotMAX spot detection (fast, never wedges)",
          config_sections=("spots",)),
    Stage("sc_trace", "sc.trace.enabled", False, ("central_element",), True,
          "per-nucleus SC tracing: SC length, fragment lower bound, fragmentation index (needs skan)",
          config_sections=("sc",)),
    Stage("granule", "granule.enabled", True, ("granule", "central_element"), True,
          "PGL-1 granule segmentation in the perinuclear region (+ lamin shell)", default_from="coloc.enabled",
          config_sections=("granule", "coloc")),
    Stage("coloc", "coloc.enabled", True, ("granule", "central_element"), True,
          "SYP <-> PGL-1 colocalization (shell voxel, partition coefficient, operands)",
          legacy_cli_help="skip PGL-1 granule surfacing + SYP<->PGL-1 colocalization stage",
          config_sections=("coloc", "sc", "granule")),
    Stage("partition", "partition.enabled", False, ("granule", "central_element"), True,
          "SYP-3 partition coefficient into P granules (distance-matched, rotation null, z-shift floor; "
          "per zone and pooled pachytene when staged)", config_sections=("partition", "envelope", "granule", "staging")),
    Stage("granule_tail", "granule_tail.enabled", False, ("granule", "central_element"), True,
          "per-granule SYP-3 excess in nuclear units and the lit fraction (no-envelope nuclei dropped)",
          config_sections=("granule_tail", "envelope", "granule")),
    Stage("qc", None, True, (), False, "QC flags and pass/fail", config_sections=("qc",)),
    Stage("render", "render.montage", True, (), False, "montage PNG", config_sections=("render",)),
    Stage("write", None, True, (), False, "tables, label images, masks, manifest", config_sections=("output",)),
)

BY_NAME: dict[str, Stage] = {s.name: s for s in STAGES}


def stage_enabled(cfg, name: str) -> bool:
    """The stage's switch as the config resolves it (code default when the key is absent)."""
    st = BY_NAME[name]
    if st.config_key is None:
        return True
    default = st.code_default
    if st.default_from is not None:
        default = bool(cfg.get(st.default_from, default))
    return bool(cfg.get(st.config_key, default))


def missing_roles(role_to_idx: dict, name: str) -> list[str]:
    return [r for r in BY_NAME[name].required_roles if role_to_idx.get(r) is None]


@dataclass
class Outcome:
    status: str                      # ran | skipped | failed
    reason: str = ""
    elapsed_s: float = 0.0
    error: str = ""
    flags: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return {"status": self.status, "reason": self.reason, "elapsed_s": round(self.elapsed_s, 3),
                "error": self.error, "flags": list(self.flags)}


def run_stage(name: str, fn: Callable[[], Any], *, flags: list[str], outcomes: dict[str, Outcome],
              enabled: bool = True, skip_reason: str = "", fatal: bool = False, fail_value: Any = None):
    """Run one stage body with the registry's uniform gate / try / flag / record behaviour.

    Returns ``fn()``'s value, or ``fail_value`` when the stage was skipped or failed (non-fatal).
    The failure flag text ``<name>:FAILED_<ExceptionType>`` is the one the pipeline has always used.
    """
    if not enabled:
        outcomes[name] = Outcome("skipped", reason=skip_reason or "disabled")
        return fail_value
    t0 = time.perf_counter()
    n_before = len(flags)
    try:
        result = fn()
    except Exception as e:  # one stage failing must not kill the batch
        if fatal:
            outcomes[name] = Outcome("failed", elapsed_s=time.perf_counter() - t0,
                                     error=f"{type(e).__name__}: {e}", flags=flags[n_before:])
            raise
        log.warning("%s stage failed (%s: %s); continuing without it.", name, type(e).__name__, e)
        flags.append(f"{name}:FAILED_{type(e).__name__}")
        outcomes[name] = Outcome("failed", elapsed_s=time.perf_counter() - t0,
                                 error=f"{type(e).__name__}: {e}", flags=flags[n_before:])
        return fail_value
    outcomes[name] = Outcome("ran", elapsed_s=time.perf_counter() - t0, flags=flags[n_before:])
    return result


def stage_hashes(cfg, role_to_idx: dict, extras: dict[str, str], model_sha: str = "",
                 spot_backend: dict | None = None, segment_backend: dict | None = None) -> dict[str, str]:
    """Per-stage provenance sub-hash: sha256 (12 hex) of the stage's resolved config sections, the
    channel roles it consumes, the versions of the tools it depends on, and (segment) the Cellpose model
    file digest. Lets two runs with different whole-config hashes be recognised as identical for one stage
    (e.g. the RAD-51 counts of a run that only switched coloc off). Written to run_manifest.json and
    <image_id>__stages.json; never a table column.

    ``spot_backend`` is the provenance a non-default spot backend declares (net3d: its model file's
    sha256 and the extra channel roles it reads); it enters the spots payload only when not empty, so
    the hashes of a default run are unchanged. ``segment_backend`` is the same for the segmentation
    backend (nuclei3d: its model file's sha256), entering the segment payload only when not empty."""
    import hashlib
    import json

    out = {}
    for st in STAGES:
        payload = {
            "config": {sec: cfg.get(sec, None) for sec in st.config_sections},
            "enabled": stage_enabled(cfg, st.name),
            "roles": {r: role_to_idx.get(r) for r in st.required_roles},
            "tools": {k: extras.get(k, "not-installed") for k in _STAGE_TOOLS.get(st.name, ())},
        }
        if st.name == "segment":
            payload["model_sha256"] = model_sha
        if st.name == "segment" and segment_backend:
            payload["segment_backend"] = segment_backend
        if st.name == "spots" and spot_backend:
            payload["spot_backend"] = spot_backend
        out[st.name] = hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()[:12]
    return out


_STAGE_TOOLS = {
    "read": ("nd2",), "segment": ("cellpose", "torch"), "measure": ("scikit-image",), "signal_qc": ("scipy",),
    "spots": ("spotmax", "cellacdc", "cupy"), "sc_trace": ("skan", "scikit-image", "scipy"),
    "granule": ("scikit-image", "scipy"),
    "coloc": ("scikit-image", "scipy"), "partition": ("scipy",), "granule_tail": ("scipy",),
}


def cli_switches() -> list[Stage]:
    return [s for s in STAGES if s.cli_switchable]


def apply_cli_switches(cfg, args) -> list[str]:
    """Apply ``--no-<stage>`` flags to the config (each is an override, so it enters ``cfg.hash``
    exactly as ``--no-spots`` always has). Returns the messages to print."""
    msgs = []
    for st in cli_switches():
        if getattr(args, f"no_{st.name}", False):
            cfg.set(st.config_key, False)
            if st.name == "spots":
                msgs.append("segmentation only: skipping spot detection (--no-spots)")
            elif st.name == "coloc":
                msgs.append("skipping PGL-1 granule surfacing + colocalization (--no-coloc)")
            else:
                msgs.append(f"skipping the {st.name} stage (--no-{st.name})")
    return msgs
