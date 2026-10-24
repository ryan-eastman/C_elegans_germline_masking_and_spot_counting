"""``python -m germquant.project new|add|edit|use-model|rounds|unlock ...``.

The torch-free entry (docs/PROJECTS.md, "The light entry"; PLAN.md decisions D28 and D39): it imports neither
germquant.cli nor germquant.pipeline, so a project write never loads torch, Cellpose, SpotMAX, CuPy, bioio or
dask on an .nd2 or TIFF image. The plugin's runner makes every project change through it, with the engine
interpreter, so the code that writes a project is the code that runs it.
"""
import sys

from ..cli_projects import main

if __name__ == "__main__":
    sys.exit(main())
