"""germquant — reproducible 3D quantification of the C. elegans meiotic germline.

Pipeline: read .nd2 -> resolve channel roles -> segment nuclei (3D Cellpose) -> measure
(spacing-aware) -> isolate the germline -> linearize the germline axis -> count spots
(RAD-51 etc. via SpotMAX) -> tidy CSV/Parquet + QC montage. Cross-validated vs Imaris
(see docs/SPOTMAX_VALIDATION.md).

The cardinal rule: voxel spacing (dz, dy, dx) in microns is read once from each .nd2 and
passed to EVERY 3D operation. Never hardcode it.
"""

__version__ = "0.1.0"
