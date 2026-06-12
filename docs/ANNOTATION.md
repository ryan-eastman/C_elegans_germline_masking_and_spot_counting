# Annotation & Cellpose fine-tuning

Zero-shot Cellpose under-segments crowded pachytene nuclei. Fine-tuning on ~50–100 hand-
annotated images is the single biggest accuracy lever (precedent: germline Cellpose model,
Jaccard 0.47 → 0.78, microPub Biology 2023). One-time effort; reused for every run after.

## 1. Export slices to annotate
```bash
germquant prep-training /nas/madeleine_images --config config/config.yaml \
    --out training/raw --n-slices 3
```
Exports 3 evenly-spaced **DAPI z-slices per gonad** (full-res, contrast-stretched 16-bit TIF)
into `training/raw/`. Use a representative spread: herm + male, HS + noHS, distal + proximal.
Aim for **~50–100 slices** total.

> 2D slices (not full 3D volumes) keep annotation tractable. Train a 2D model, then apply it
> in 3D with `do_3D=True, anisotropy=dz/dxy` — the pipeline already does this.

## 2. Annotate the nuclei
Easiest: the **Cellpose GUI** (`pip install "germquant[gpu]"` gives you `cellpose`).
```bash
python -m cellpose            # opens the GUI
```
- Open each `training/raw/*.tif`, paint/curate nucleus masks (right-drag to draw; it
  saves `<name>_seg.npy` next to the image).
- Or annotate in **napari** + `napari-cellpose` and export label masks as `<name>_masks.tif`.

Curate, don't draw from scratch: run the generalist model first, then fix the merges/splits.

## 3. Fine-tune
```bash
germquant finetune training/raw --out-model models/germline_nuclei --epochs 100
# (add --print-only to just see the cellpose command without running it)
```
This trains from `cpsam` and writes the model to `models/germline_nuclei`. GPU only.
Cellpose CLI flags vary by version — if it errors, run the printed command manually.

## 4. Use the fine-tuned model
Point the config at it — no code change:
```yaml
# config/config.yaml
segmentation:
  nuclei:
    method: cellpose
    cellpose_model: models/germline_nuclei
    diameter_um: 3.0      # set to your measured mean pachytene-nucleus diameter
```

## 5. Verify it actually helped
Hold out a few annotated gonads, then:
```bash
germquant validate --seg-pred results/.../X__nuclei_labels.tif \
                   --seg-truth training/truth/X_masks.tif --iou 0.5 --out validation/
```
Reports detection **F1 + mean IoU** (Metrics Reloaded style). Track this before/after fine-tuning.
