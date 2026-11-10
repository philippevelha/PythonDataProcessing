# PGC LVM data analysis

Python/Jupyter port of `LVM_reader.m` for the LabVIEW `.lvm` files recorded with the PGC program.

- `lvm_reader.py`: the LVM parser (handles multi-segment files and files with no extension) plus the analysis functions (periodogram PSD, low-pass detrend, minimum block std, phase unwrap)
- `notebooks/LVM_analysis.ipynb`: analysis of a single file and comparison of several files
- `site/browser.html`: web trace browser that lets you pick any recording and zoom through its time series
- `build_traces.py`: generates the browser data (`site/traces/`) from the raw files

## Run locally

```bash
pip install -r requirements.txt
cd notebooks
jupyter lab
```

Set `DATA_DIR`, `SINGLE_FILE` and `MULTI_FILES` in the configuration cell, then run all cells.
The plots are interactive Plotly figures, and they stay zoomable on the published page. Long traces are min/max decimated to `MAX_POINTS` per trace so the page stays small (about 5 MB). To see full resolution, set `T_WINDOW` or raise `MAX_POINTS`.

## Updating the trace browser

After adding or changing measurement files, run

```bash
python build_traces.py            # or: --root <data folder>  --bins 20000
```

The script writes a min/max envelope of every LVM time-series recording (10,000 bins per channel by default, `ModulatingSignal` left out) to `site/traces/`. It skips duplicates and files in other formats, such as wavelength scans. Commit `site/traces/` and push.
Recordings with both `Phase1` and `Phase2` also get a `Phase1-Phase2` channel (push-pull differential phase), computed sample by sample before the envelope is taken. It has its own panel and can be switched off like any other channel.
To preview locally, run `python -m http.server -d site` and open http://localhost:8000/browser.html.

### Event markers

In the browser, press **M** (or "Add markers") and click on the plot to place a vertical marker. Drag a marker to move it, and type its label in the list below the controls. **Save markers** writes the markers of every recording to `markers.json`. In Chrome or Edge, pick `site/markers.json` the first time; later saves in the same session overwrite it directly. Other browsers download the file, which you then move to `site/`. Commit `site/markers.json` so the published page shows the markers. Edits that are not saved yet are kept in this browser and restored on the next visit.

`markers.json` format: `{"markers": {"<recording id>": {"file": "<folder>/<name>", "markers": [{"t": <seconds>, "label": "..."}]}}}`.

The notebook reads this file (`MARKERS_FILE`) and draws the markers as dashed lines on the single-file time traces and on the chronogram. In Python: `lvm.load_markers("site/markers.json")` returns a table (`id`, `file`, `t`, `label`) and `lvm.markers_for(markers, path)` selects one recording's markers.

## Publishing

The raw data is **not** in the repository (see `.gitignore`). Run the notebook locally, save it with its outputs and push to `main`.
The `Publish notebooks` workflow then renders it to HTML on GitHub Pages without re-executing it.
In the repository, set **Settings → Pages → Source: GitHub Actions** once.
