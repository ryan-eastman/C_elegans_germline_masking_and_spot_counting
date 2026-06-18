# Counting RAD-51 spots in *C. elegans* germlines — a step-by-step guide

This program takes a confocal image of a worm germline (a `.nd2` file) and, **all by itself**, finds
every nucleus and counts the **RAD-51 spots** inside each one. It writes the answers into spreadsheet
files (`.csv`) you open in Excel. No clicking through Imaris, no counting by hand.

**You do not need to know how to code.** This guide assumes you have *never* used a "terminal" before.
If a word is confusing, check **"Words you might not know"** near the bottom.

---

## 🟢 The easy way (one image, drag-and-drop)

1. In this folder, find the file **`quantify.bat`**.
2. **Drag your `.nd2` image on top of `quantify.bat`** and let go.
3. A black window opens and starts working. **Leave it open** — it takes about **15–25 minutes**.
4. When it prints **DONE**, your answers are in a new folder right next to your image, named
   `yourimagename_results`.

That's the whole thing. Everything below explains the results, how to do many images, and what to do
if it misbehaves.

---

## Before you start (what you need)

- **This computer** — the one with the big NVIDIA graphics card. It's already set up. (On a different
  computer? See *"Setting up a new computer"* at the bottom, or ask whoever manages the lab computers.)
- **Your image**: a per-gonad `.nd2` file with DAPI / SYP / RAD-51 channels. The whole-slide overview
  files (names with **`10x`** or **`largeimage`**) will not work — use the individual gonad files.
- **About 20 minutes per image**, and patience. The computer does the work; you just wait.

---

## 🟡 The typing way (if drag-and-drop won't work, or you want to choose where results go)

1. **Open PowerShell:** click the **Start** menu (Windows icon, bottom-left of the screen), type
   **`PowerShell`**, and click **"Windows PowerShell"**. A black window opens.
2. **Go to the program's folder** — copy this line, paste it into the window (right-click pastes),
   and press **Enter**:
   ```
   cd "C:\Users\ryane\C_elegans_germline_masking_and_spot_counting"
   ```
   *(If the code lives in a different folder on your computer, use that folder instead.)*
3. **Run one image** — copy the line below, but **replace the two paths with yours**, then Enter:
   ```
   .venv\Scripts\germquant.exe run "C:\path\to\YOUR_IMAGE.nd2" --config config\config.yaml --out "C:\path\to\WHERE_RESULTS_GO"
   ```
   - The first quotes hold **your image file**. *Tip:* instead of typing it, **drag the `.nd2` file
     into the window** — it pastes the correct path for you.
   - After `--out` is **any empty folder** where you want the results saved.
4. Wait ~15–25 minutes. When the blinking prompt comes back (with a ✓), it's finished.

---

## 🟠 Doing a whole folder at once

To process **every** `.nd2` in a folder (it automatically skips the `10x`/`largeimage` overviews):
```
.venv\Scripts\germquant.exe batch "C:\path\to\FOLDER_OF_IMAGES" --config config\config.yaml --out "C:\path\to\RESULTS_FOLDER"
```
It makes one results sub-folder per image, plus a `batch_summary.csv` listing them all. Many images
can take hours — it does them one at a time.

---

## 📊 Your results — where they are and what they mean

Inside the results folder, each image produces files that start with the image's name:

| File (open in **Excel**) | What it is |
|---|---|
| **`..._nuclei.csv`** | **The main result.** One row per nucleus. The **`n_spots`** column = RAD-51 spots in that nucleus. The **`in_germline`** column is `True` for real germline nuclei and `False` for gut/debris (ignore the `False` ones). |
| `..._spots.csv` | One row per individual RAD-51 spot (its 3D position + brightness). For deeper analysis. |
| `..._image_summary.csv` | One row of totals for the whole image (number of nuclei, average spots, pass/fail QC). |
| `..._montage.png` | A **picture to eyeball that it worked** — nucleus outlines with the spots marked. Open it and glance at it every time. |
| `..._nuclei_labels.tif` | The 3D nucleus map. Load it into Imaris as **Surfaces** to double-check the segmentation. |
| `..._spots.tif` | The detected spots as a 3D image, on the same grid as the nuclei map. Load it into Imaris as a **Channel** (Edit → Add Channels), or run Imaris **Spots** detection on it (diameter ~0.4 µm) to get Spots objects next to your Surfaces. |

**To get "average RAD-51 spots per nucleus":** open `..._nuclei.csv` in Excel → filter so
`in_germline = True` → average the **`n_spots`** column.

> ⚠️ **Read this before trusting the numbers.** The spot counts are **validated against Imaris for N2
> (wild-type) worms — both no-heat-shock and heat-shock** (they match closely, ~1:1, across 4–21
> foci/nucleus). For **other genotypes (mutants)** they have **not** been checked yet — re-check
> against Imaris first, and ask whoever maintains this tool before reporting absolute counts on a new
> genotype. (The full why is in [docs/SPOTMAX_VALIDATION.md](docs/SPOTMAX_VALIDATION.md).)

---

## 🔧 If something breaks

- **Red error about the `.venv` folder missing** → the program isn't installed on this computer. See
  *Setting up a new computer*.
- **"running scripts is disabled on this system" / Windows blocked the file** → use the **typing way**
  (PowerShell) instead of double-clicking, or ask a lab tech to allow scripts.
- **"file not found"** → the image path is wrong. Drag the `.nd2` into the window to paste it
  correctly. Paths with spaces are fine **as long as they're inside quotes**.
- **It found 0 nuclei or used the wrong channels** → your image's channels may differ from
  DAPI/SYP/RAD-51. Tell whoever maintains this; the channel setup lives in `config\channel_maps\`.
- **The black window vanished instantly** → run it the **typing way** so you can read the error.
- **"out of memory" (GPU)** → close other heavy programs (**especially Imaris**) and try again. Only
  one big image can run at a time.
- **Still stuck?** Copy the red text and send it to whoever maintains this tool.

---

## 📖 Words you might not know

- **`.nd2`** — the raw image file the Nikon confocal saves.
- **Terminal / PowerShell / "the black window"** — where you type commands. Scrolling text is normal;
  it's just the program thinking out loud.
- **Path** — a file's full address, e.g. `C:\Users\you\images\worm1.nd2`. Drag a file into the window
  to paste its path.
- **`.csv`** — a spreadsheet file; double-click to open it in Excel.
- **Segmentation** — the computer outlining each nucleus in 3D.
- **RAD-51 / focus (plural foci) / "spot"** — the protein whose dots you're counting; they mark DNA
  double-strand breaks during meiosis.

---

## 🛠 Setting up a new computer (one time — for a lab tech)

Requirements: 64-bit Windows, an NVIDIA GPU (RTX 30/40/50-series), ~10 GB free.
1. Install **Python 3.11** from [python.org](https://www.python.org/downloads/) — during install,
   **tick "Add Python to PATH"**.
2. Install the **uv** helper: in PowerShell, run `pip install uv`.
3. Get this code: download the ZIP from GitHub and unzip it (or `git clone` it).
4. In PowerShell, `cd` into the folder, then run these one at a time:
   ```
   uv venv
   uv pip install -e .
   uv pip install spotmax cellacdc
   uv pip install torch --index-url https://download.pytorch.org/whl/cu128
   ```
5. Put the trained model at `models\models\germline_nuclei_combined` (ask Ryan for it, ~1.2 GB).
6. Confirm the GPU is seen: `.venv\Scripts\germquant.exe check-gpu` (it should print your card).
7. Test the whole thing: drag a small `.nd2` onto `quantify.bat`.

---

## For developers

Pipeline: `read .nd2 → segment nuclei (Cellpose) → isolate germline → linearize axis → count spots
(SpotMAX) → tidy CSV/Parquet + QC montage`. The `spots:` block in `config/config.yaml` holds the
cross-validated detection parameters. Validation + cross-validation tooling is in `scripts/cv_*.py`
with Imaris readers in `germquant.validate`; the calibration story and honest caveats are in
[docs/SPOTMAX_VALIDATION.md](docs/SPOTMAX_VALIDATION.md). Tests: `pytest` (CPU, no GPU/data needed).
