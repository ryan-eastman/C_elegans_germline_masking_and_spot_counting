"""``python -m germquant.packs list|show|check|fetch-models|install|verify|build-sums|lint|--version``.

The torch-free entry (docs/DESIGN_system_packs.md, sections 3.9 and 8.1): it imports neither germquant.cli
nor germquant.pipeline, so listing packs and checking them never loads torch, Cellpose, SpotMAX or CuPy. The
plugin's subprocess runner lists and checks packs through it with the engine interpreter, and compares
``--version`` (engine version, git sha, python) with its own germquant.
"""
import sys

from ..cli_packs import main

if __name__ == "__main__":
    sys.exit(main())
