"""nuclei3d: an EXPERIMENTAL in-house 3D nucleus segmentation model, a compact residual 3D U-Net with three
heads (foreground, distance to the nucleus boundary, closeness to the nearest neighbouring nucleus; the
cell and neighbour distances of Scherr et al. 2020), trained on the lab's Imaris nucleus surfaces with a
partial-annotation loss (docs/DESIGN_inhouse_components.md, Component 2; docs/BENCHMARK.md, "nuclei3d").

Modules: ``model`` (the network), ``data`` (targets and per-voxel weights from the recovered Imaris masks
and the deployed Cellpose labels, patches from the bench cache), ``train`` (``germquant train-nuclei3d``),
``infer`` (normalisation, the model file, deterministic tiled inference, the h-maxima seeds and the seeded
watershed) and ``backend`` (the ``nuclei3d`` segmentation backend, selected only by
``segmentation.nuclei.backend: nuclei3d``). Importing this package does not import torch; ``model``,
``train`` and the inference functions do.
"""
