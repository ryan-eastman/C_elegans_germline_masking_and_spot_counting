"""net3d: an EXPERIMENTAL learned RAD-51 foci detector, a compact residual 3D U-Net that predicts a
per-voxel heatmap of focus centres from the RAD-51 and DAPI channels, trained on the lab's Imaris spots.

Modules: ``model`` (the network), ``data`` (patches, targets and loss weights from the bench cache),
``train`` (``germquant train-spots``), ``infer`` (normalisation, the model file, deterministic tiled
inference), ``peaks`` (peak finding and nucleus assignment, numpy only) and ``backend`` (the ``net3d``
spot backend, selected only by ``spots.backend: net3d``). Importing this package does not import torch;
``model``, ``train`` and the inference functions do. docs/BENCHMARK.md ("net3d") has the training
record and the validation numbers; docs/BACKENDS.md the backend's keys.
"""
