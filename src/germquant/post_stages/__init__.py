"""Post-stages: measurement steps that an outside lab adds as its own installed package (docs/POST_STAGES.md;
design 5.6).

A post-stage is a distribution with an entry point in the group ``germquant.post_stages`` and a declaration file,
germquant_post_stage.yaml, in its top-level package: what it reads (channel roles, tables, the label image, the
image), the tables and columns it writes, and its typed parameters, declared like recipe tunables. A run switches
post-stages on with a separate file, post_stages.yaml (a list of {name, params}), which is never part of
config_hash. After process_image has finished, ``run_post_stages`` calls each one on what the run wrote and
writes only new tables, <image_id>__post_<name>_<table>.csv, plus the record <image_id>__post.json with each
stage's sub-hash; a changed sub-hash reruns that post-stage only, and a failure is recorded, never failing the
image. Post-stages install only into .venv_packs or a lab's own environment, never into the validated .venv
(rule R6).

Importing this package, listing post-stages and checking or writing post_stages.yaml read metadata, the
declaration files and YAML only: no post-stage module, torch, CuPy, SpotMAX, Cellpose, bioio, numpy or pandas
is loaded (gate G6). ``run_post_stages`` imports germquant.pipeline first, so torch claims its cuBLAS before a
reader or the post-stage does (D28). The torch-free command entry is ``python -m germquant.post_stages list
[--json] [--all]`` and ``write DIR --entries JSON [--json]``. germquant.post_stages.testing holds the fixture
post-stage that the engine and plugin tests share.

The run commands take post-stages (docs/POST_STAGES.md, "Running post-stages"): `germquant run` and `germquant
batch` with ``--post-stages FILE``, api.run and api.run_loaded with ``post_stages=``, and `germquant project run`
with the project's post_stages.yaml; a pack recipe's post_stages/<recipe>.yaml is in force unless one of those
replaces it (post_stages_in_force). A run of a recipe file or a profile without them imports none of this
package; a pack run imports it to look for the pack's file.
"""
from __future__ import annotations

from .config import (
    PACK_POST_STAGES_DIR,
    POST_STAGES_FILE,
    dump_post_stages,
    entries_of,
    load_post_stages_file,
    pack_post_stages_path,
    post_stages_in_force,
    post_stages_sha256,
    validate_post_stages,
    write_post_stages_file,
)
from .declaration import DECLARATION_FILE, PARAM_TYPES, declaration_problems
from .discovery import (
    ENTRY_POINT_GROUP,
    PostStageError,
    declarations,
    list_post_stages,
    load_post_stage,
    register_post_stage,
    unregister_post_stage,
)
from .runner import (
    SERIES_REFUSAL,
    PostInputs,
    PostStageContractError,
    post_record_path,
    post_table_name,
    post_table_problems,
    read_post_record,
    run_post_stages,
)

__all__ = [
    "DECLARATION_FILE", "ENTRY_POINT_GROUP", "PARAM_TYPES", "POST_STAGES_FILE", "SERIES_REFUSAL", "PostInputs",
    "PostStageContractError", "PostStageError", "declaration_problems", "declarations", "dump_post_stages",
    "list_post_stages", "load_post_stage", "load_post_stages_file", "post_record_path", "post_stages_sha256",
    "post_table_name", "read_post_record", "register_post_stage", "run_post_stages", "unregister_post_stage",
    "validate_post_stages", "write_post_stages_file",
    # the run routes (phase WIRE): which list a run uses, and which post tables collect gathers
    "PACK_POST_STAGES_DIR", "entries_of", "pack_post_stages_path", "post_stages_in_force", "post_table_problems",
]
