"""SYP-3 partition coefficient into P granules and the per-granule lit fraction (stages `partition`
and `granule_tail`; docs/ROADMAP_modular_pipeline.md step 10), ported from the August 2026 ccw77
analysis (analysis/coloc/scripts/pc_lamin_worker.py, pc_zone_worker.py, granule_tail.py).

Partition coefficient PC = (SYP in granule voxels - bkg) / (SYP in non-granule cytoplasm at the SAME
distance from the nuclear envelope - bkg), averaged over 0.25 um distance bins weighted by granule
voxels, inside a 2.5 um cytoplasm shell; bkg = 3rd percentile of the SYP crop. Controls: a 180 degree
rotation null (the granule mask rotated in-plane, measured against its own outside set) and a z-shift
floor (the mask shifted 6 planes, so it samples the same blur but no granule); PC_specific = PC / z-shift.
Regions: the whole shell, each hand-traced pachytene zone (a voxel inherits the zone of the nucleus
whose envelope is nearest) and the three zones pooled ("pach", the contamination-free whole-gonad
measure). The whole-gonad gates (50 / 500 voxels per region, 20 / 50 per bin) and the zone gates
(30 / 200, 15 / 40) are the workers' own.

Lit fraction (granule_tail): per granule, excess = mean SYP in the granule minus the non-granule
cytoplasm at the same distance bin, in units of the nuclear SYP level (N = mean SYP inside the
envelope minus bkg, bkg here = 0.5th percentile far from the envelope); per gonad the fraction of
granules with excess above 0.25 / 0.5 / 1.0 and percentiles. Version 4 of the analysis: nuclei with no
lamin envelope are removed from the envelope set before the shell is built and the granules
re-segmented inside it, exactly as granule_tail.py did.

Measurability (measurable.py, docs/PARTITION_QC.md): after the stage has run, the pipeline reads the SYP
signal-to-noise the signal_qc stage measured (image_summary field sig_central_element_snr) into
partition_syp_snr and partition_measurable; with partition.min_syp_snr set in a profile (code default
0.0, off) an image below it is flagged partition:low_syp_snr while every number is still written. A SYP
channel with no live fluorophore gives a partition coefficient floor above 1 that no control removes.
"""
from .measurable import (
    LOW_SNR_FLAG,
    MEASURABILITY_FIELDS,
    MIN_SYP_SNR_DEFAULT,
    SYP_SNR_FIELD,
    partition_measurability,
    syp_snr_from_summary,
)
from .pc import DBINS_UM, partition_metrics, pc_dm, run_partition, zone_voxel_map
from .tail import run_granule_tail

__all__ = ["DBINS_UM", "partition_metrics", "pc_dm", "run_partition", "zone_voxel_map", "run_granule_tail",
           "LOW_SNR_FLAG", "MEASURABILITY_FIELDS", "MIN_SYP_SNR_DEFAULT", "SYP_SNR_FIELD",
           "partition_measurability", "syp_snr_from_summary"]
