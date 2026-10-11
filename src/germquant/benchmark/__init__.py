"""Benchmark harness: our pipeline output against the Imaris expert truth, gonad by gonad.

Two benchmarks share this package:

* ``spots``: the RAD-51 spot detector. Expert spots are mapped into our nucleus labels for a
  per-nucleus count comparison (Pearson, Lin's CCC, bias), and our spots are matched one to one with
  the expert's inside an anisotropic cylinder for precision, recall and F1.
* ``nuclei3d``: the 3D nucleus segmentation. Our label image is matched to the expert's partial label
  image by IoU (F1 at IoU 0.5, AP averaged over IoU 0.5 to 0.95, splits, merges, z-extent error), with
  a partial-annotation rule so untraced nuclei far from the expert's region are not false positives
  (``coverage="protocol"`` is the pre-registered scope of docs/BENCHMARK_PROTOCOL.md section 5.3).

``stats`` holds the agreement statistics and the cluster bootstrap by gonad; ``matching`` holds the
one-to-one assignment that both benchmarks use (the Hungarian algorithm on the connected components of
the sparse candidate graph, so no dense 30000 by 30000 matrix is ever allocated); ``manifest`` reads
the benchmark manifest, the pipeline output folder and the truth files.

Manifest shape (built by scripts/benchmark_manifest.py, read by ``manifest.load_manifest``): a JSON
list of entries, or an object with an ``entries`` list. Every entry has

    tag          the golden tag, e.g. "n2rad51_noHS_herm_001" (used as the row key everywhere)
    split        a grouping label, e.g. "noHS" or "HS", or a cross-validation fold
    results_dir  the pipeline output folder: <image>__nuclei.csv, <image>__spots.csv,
                 <image>__nuclei_labels.tif
    truth_dir    the folder holding the truth files for this gonad (see below)
    nd2          the source .nd2 (used only for the voxel size when the tables lack it)
    ims          the Imaris project the truth was exported from (its stem names the truth files)

Truth file contract (produced by the Imaris reader package, consumed here): for an .ims named
<stem>.ims the truth folder holds ``<stem>__imaris_spots.csv``, ``<stem>__imaris_nuclei.csv``,
``<stem>__imaris_labels.tif`` (optional) and ``<stem>__imaris_params.json``. Coordinates in the
``x_um``/``y_um``/``z_um`` columns are already in OUR frame: Imaris micrometres were turned into voxel
indices with Imaris's own (rounded) voxel size and then multiplied by the .nd2 voxel size, which is
how the pipeline's own coordinates are formed. The integer ``xi``/``yi``/``zi`` columns are those
voxel indices. docs/BENCHMARK.md documents every column and output file.

Nothing in here touches a pipeline stage; the package reads outputs, it never produces them.
"""
from .manifest import (
    Entry,
    PipelineOutputs,
    Truth,
    find_pipeline_outputs,
    load_labels,
    load_manifest,
    load_truth,
    our_voxel_zyx,
    roi_box_from_params,
    truth_file_names,
)
from .matching import Assignment, assign_dense, assign_greedy, assign_sparse
from .stats import (
    agreement_metrics,
    bias,
    cluster_bootstrap,
    lin_ccc,
    paired_comparison,
    pearson,
    precision_recall_f1,
)

__all__ = [
    "Entry", "PipelineOutputs", "Truth", "find_pipeline_outputs", "load_labels", "load_manifest",
    "load_truth", "our_voxel_zyx", "roi_box_from_params", "truth_file_names",
    "Assignment", "assign_dense", "assign_greedy", "assign_sparse",
    "agreement_metrics", "bias", "cluster_bootstrap", "lin_ccc", "paired_comparison", "pearson",
    "precision_recall_f1",
]
