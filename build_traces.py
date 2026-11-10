"""Build the data for the web trace browser (site/browser.html).

Scans DATA_ROOT for LVM time-series recordings and writes, for each one, a compact
min/max envelope (BINS bins per channel, float32) to site/traces/<id>.bin, plus
site/traces/manifest.json describing every recording.

    python build_traces.py                 # scan this folder, default settings
    python build_traces.py --bins 20000    # finer zoom resolution (bigger files)
    python build_traces.py --root E:/MeasurementsDATA

Binary layout of <id>.bin (little-endian float32), for each channel in manifest order:
    raw files (samples <= bins): values[n]
    otherwise:                   min[n], max[n]
Bin k starts at t0 + k * bin_dt (seconds).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
import warnings
from pathlib import Path

import numpy as np

import lvm_reader as lvm

SKIP_DIRS = {".git", ".github", "notebooks", "site", "__pycache__", "_site"}
SKIP_SUFFIXES = {".mat", ".fig", ".jpg", ".png", ".pdf", ".m", ".py", ".ipynb", ".md", ".txt"}
DEFAULT_EXCLUDE_CHANNELS = ["ModulatingSignal"]   # PGC carrier: its envelope is just a constant band


def file_fingerprint(p: Path) -> tuple[int, str]:
    with open(p, "rb") as f:
        return p.stat().st_size, hashlib.md5(f.read(1 << 16)).hexdigest()


def find_recordings(root: Path) -> list[Path]:
    out, seen = [], {}
    for p in sorted(root.rglob("*")):
        rel = p.relative_to(root)
        if not p.is_file() or rel.parts[0] in SKIP_DIRS or p.suffix.lower() in SKIP_SUFFIXES:
            continue
        if not lvm.is_lvm(p):
            continue
        fp = file_fingerprint(p)
        if fp in seen:              # e.g. "spike01" and "spike01.lvm" are the same recording
            print(f"  skip duplicate {rel} (= {seen[fp]})")
            continue
        seen[fp] = rel
        out.append(p)
    return out


def envelope(y: np.ndarray, step: int, nb: int) -> tuple[np.ndarray, np.ndarray]:
    pad = np.full(nb * step, np.nan)
    pad[: len(y)] = y
    blocks = pad.reshape(nb, step)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # all-NaN bins -> NaN (gap in the plot)
        return np.nanmin(blocks, axis=1), np.nanmax(blocks, axis=1)


def recording_id(rel: Path) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", "__".join(rel.parts))


def build(root: Path, out_dir: Path, bins: int, exclude: list[str]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    files = find_recordings(root)
    print(f"{len(files)} recordings found")
    manifest, written = [], set()
    t_start = time.time()

    for i, p in enumerate(files, 1):
        rel = p.relative_to(root)
        try:
            d = lvm.read_lvm(p)
        except Exception as e:  # noqa: BLE001 - report and carry on with the other files
            print(f"  [{i}/{len(files)}] FAILED {rel}: {e}")
            continue
        n = len(d.data)
        raw = n <= bins
        step = 1 if raw else int(np.ceil(n / bins))
        nb = int(np.ceil(n / step))
        channels = [c for c in d.channels if c not in exclude]

        parts, ch_meta = [], []
        for c in channels:
            y = d.data[c].to_numpy(dtype=float)
            valid = y[~np.isnan(y)]
            if raw:
                parts.append(y.astype("<f4"))
            else:
                lo, hi = envelope(y, step, nb)
                parts += [lo.astype("<f4"), hi.astype("<f4")]
            ch_meta.append(dict(name=c,
                                min=float(valid.min()) if len(valid) else None,
                                max=float(valid.max()) if len(valid) else None))

        rid = recording_id(rel)
        blob = np.concatenate(parts).tobytes()
        (out_dir / f"{rid}.bin").write_bytes(blob)
        written.add(f"{rid}.bin")
        manifest.append(dict(
            id=rid, folder=str(rel.parent).replace("\\", "/"), name=p.name,
            date=d.header.get("Date", ""), time=d.header.get("Time", "")[:8],
            operator=d.header.get("Operator", ""),
            fs=d.fs, samples=n, duration=float(d.data.time.iloc[-1] - d.data.time.iloc[0] + d.dt),
            segments=d.n_segments, t0=float(d.data.time.iloc[0]), bin_dt=step * d.dt,
            bins=nb, raw=raw, channels=ch_meta, file=f"traces/{rid}.bin", bytes=len(blob),
        ))
        print(f"  [{i}/{len(files)}] {rel}  {n:,} samples -> {nb:,} bins  ({len(blob)/1e3:.0f} kB)")
        del d

    # remove .bin files of recordings that no longer exist
    for old in out_dir.glob("*.bin"):
        if old.name not in written:
            old.unlink()

    manifest.sort(key=lambda m: (m["date"], m["time"], m["name"]))
    (out_dir / "manifest.json").write_text(json.dumps(
        dict(generated=time.strftime("%Y-%m-%d %H:%M"), bins=bins, recordings=manifest), indent=1))
    total = sum(m["bytes"] for m in manifest)
    print(f"done: {len(manifest)} recordings, {total/1e6:.1f} MB, {time.time()-t_start:.0f} s")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=Path(__file__).parent, help="folder containing the measurements")
    ap.add_argument("--out", type=Path, default=Path(__file__).parent / "site" / "traces")
    ap.add_argument("--bins", type=int, default=10_000, help="envelope bins per channel (zoom resolution)")
    ap.add_argument("--exclude", nargs="*", default=DEFAULT_EXCLUDE_CHANNELS, help="channels to leave out")
    a = ap.parse_args()
    build(a.root.resolve(), a.out.resolve(), a.bins, a.exclude)
