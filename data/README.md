# data/  (gitignored — nothing here is committed)

Drop inputs here for local checks. The whole folder is in `.gitignore` (except this README),
so raw microscopy and ground truth never get committed.

## raw_examples/
Real `.nd2` germline stacks for end-to-end checks (ideally one HERM/oocyte and one MALE/spermatocyte;
control + heat if easy).

    germquant info  data/raw_examples/<file>.nd2          # channels + voxel size
    germquant run   data/raw_examples/<file>.nd2 --config config/config.yaml --out results_test/

## ground_truth/
Imaris ground truth for validating the spot counts. The readers in `germquant.validate` consume the
lab's native Imaris exports directly — no reformatting needed:
- **`.xlsx`** Imaris Statistics export (the coloc method) — `germquant.validate.imaris_xlsx`.
- **`.ims`** Imaris project (Spots + Surfaces) — `germquant.validate.imaris_ims`.

Validate against them with:

    python scripts/imaris_gt_counts.py      # per-gonad RAD-51/nucleus from the xlsx coloc method
    python scripts/validate_same_image.py   # pipeline output vs the .ims for one image
    # cross-gonad CV: scripts/cv_run.py -> build_cv_manifest.py -> cv_detect.py  (see docs/RUNBOOK.md)

Imaris masks only a *subset* of nuclei, so only **recall** and **per-nucleus RAD-51** are valid
comparisons (not nucleus count / precision / F1). See [docs/SPOTMAX_VALIDATION.md](../docs/SPOTMAX_VALIDATION.md).
