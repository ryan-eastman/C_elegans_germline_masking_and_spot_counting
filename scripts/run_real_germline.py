"""Full-res real confirmation of the new path: native-scale combined model + germline isolation +
colored montage, on the Cahoon control gonad. Confirms the integrated pipeline on real data.

    python scripts/run_real_germline.py
"""
import pandas as pd

from germquant import pipeline
from germquant.config import load_config

REAL = r"data/raw_examples/madeleline images/20251105_N2_noHS/20251105_N2_nohs_HERM _001.nd2"
OUT = "results_germline"


def main():
    cfg = load_config("config/config.yaml")
    cfg._data["segmentation"]["nuclei"]["cellpose_model"] = "models/models/germline_nuclei_combined"
    res = pipeline.process_image(REAL, cfg, OUT, prov=None)
    s = pd.read_csv(f"{OUT}/{res['image_id']}__image_summary.csv")
    print(f"DONE image={res['image_id']}")
    print(f"  n_nuclei(total)={int(s['n_nuclei'].iloc[0])}  n_germline={int(s['n_germline_nuclei'].iloc[0])}"
          f"  -> dropped {int(s['n_nuclei'].iloc[0]) - int(s['n_germline_nuclei'].iloc[0])} off-gonad")
    print(f"  mean_foci(germline)={float(s['mean_foci'].iloc[0]):.2f}  montage={OUT}/{res['image_id']}__montage.png")
    print(f"  qc_flags={s['qc_flags'].iloc[0]}")


if __name__ == "__main__":
    main()
