"""germquant — reproducible 3D quantification of the C. elegans meiotic germline.

Pipeline: read .nd2 -> resolve channel roles -> segment nuclei (3D) -> linearize
germline axis -> call zones (TZ / pachytene) -> trace SC length & fragmentation ->
detect RAD-51 foci -> measure (spacing-aware) -> tidy CSV/Parquet -> renders/QC.

The cardinal rule (see docs/ARCHITECTURE.md): voxel spacing (dz, dy, dx) in microns
is read once from each .nd2 and passed to EVERY 3D operation. Never hardcode it.
"""

__version__ = "0.1.0"
