# PGC LVM data analysis

Python/Jupyter port of `LVM_reader.m` for the LabVIEW `.lvm` files recorded with the PGC program.

- `lvm_reader.py`: the LVM parser (handles multi-segment files and files with no extension) plus the analysis functions (periodogram PSD, low-pass detrend, minimum block std, phase unwrap)
- `notebooks/LVM_analysis.ipynb`: analysis of a single file and comparison of several files

## Run locally

```bash
pip install -r requirements.txt
cd notebooks
jupyter lab
```

Set `DATA_DIR`, `SINGLE_FILE` and `MULTI_FILES` in the configuration cell, then run all cells.
The plots are interactive Plotly figures, and they stay zoomable on the published page. Long traces are min/max decimated to `MAX_POINTS` per trace so the page stays small (about 5 MB). To see full resolution, set `T_WINDOW` or raise `MAX_POINTS`.

## Publishing

The raw data is **not** in the repository (see `.gitignore`). Run the notebook locally, save it with its outputs and push to `main`.
The `Publish notebooks` workflow then renders it to HTML on GitHub Pages without re-executing it.
In the repository, set **Settings → Pages → Source: GitHub Actions** once.
