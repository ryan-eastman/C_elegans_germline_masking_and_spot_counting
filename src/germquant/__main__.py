"""``python -m germquant ...``: the full command line, exactly as the germquant executable.

It enters cli.main, and cli.py's first import is torch (the cuBLAS load-order guard: torch must claim its
cuBLAS before any SpotMAX or CuPy import, or Cellpose silently drops to the watershed on the RTX 5090), so
this entry keeps that order. Light work (listing and checking packs) goes through ``python -m
germquant.packs`` instead, which never loads torch.
"""
from .cli import main

if __name__ == "__main__":
    raise SystemExit(main())
