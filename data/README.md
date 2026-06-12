# data/  (gitignored — nothing here is committed)

Drop inputs here for local verification. The whole folder is in `.gitignore`
(except this README), so raw microscopy and ground truth never get committed.

## raw_examples/
Real `.nd2` germline stacks for end-to-end checks. A couple of representative files is
plenty (ideally one HERM/oocyte and one MALE/spermatocyte; control + heat if easy).

    germquant info  data/raw_examples/<file>.nd2          # channels + voxel size
    germquant run   data/raw_examples/<file>.nd2 --config config/config.yaml \
                    --out results_test/ --xy-stride 4      # CPU smoke (classical seg)

## ground_truth/
Imaris (or Fiji/SNT) hand quantifications, as CSV, for `germquant validate` /
`scripts/run_and_validate.sh`. Two files (names matter — the script looks for them):

**`per_nucleus.csv`** — one row per nucleus with a join key:

| column                  | meaning                                              |
|-------------------------|------------------------------------------------------|
| `image_id`              | matches the `.nd2` stem (the pipeline's image_id)    |
| `nucleus_id`            | per-nucleus id (any consistent numbering)            |
| `sc_total_length_um`    | Imaris filament length per nucleus (µm)              |
| `n_fragments`           | number of SC filaments/fragments per nucleus         |
| `n_foci`                | RAD-51 foci per nucleus                               |

**`zones.csv`** — one row per image:

| column                  | meaning                                              |
|-------------------------|------------------------------------------------------|
| `image_id`              | matches the `.nd2` stem                              |
| `transition_zone_um`    | TZ length (µm)                                        |
| `pachytene_um`          | pachytene length (µm)                                |

If Imaris nucleus numbering can't be matched to the pipeline's, that's fine — provide
per-image values and we validate at the distribution level (means/medians, CCC) instead.
Don't have all columns? Provide what you have; the script skips the rest.
