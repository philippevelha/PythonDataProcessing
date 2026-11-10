"""Reader and analysis helpers for LabVIEW Measurement (.lvm) files written by the PGC program.

Python port of LVM_reader.m.

File layout handled here (Multi_Headings = Yes, X_Columns = Multi):
    - a file header terminated by ``***End_of_Header***``
    - one or more *segments*, each with its own header (Channels, Samples, Delta_X, ...)
      terminated by ``***End_of_Header***``, a column-name line
      (``X_Value  Phase1  X_Value  Phase2 ...  Comment``) and the numeric rows.
      The time column restarts at X0 (=0) in every segment.

Segments are concatenated and the time axis is made continuous, as in the MATLAB
``lastvalue`` logic (here the offset is ``previous_end + dt`` so no timestamp is duplicated).
"""

from __future__ import annotations

import io
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import signal

END_OF_HEADER = "***End_of_Header***"
_NUMERIC_START = set("0123456789-+.")


@dataclass
class LVMData:
    path: Path
    header: dict[str, str]            # file-level header (Date, Time, Operator, ...)
    segment_header: dict[str, list]   # first segment header (Delta_X, Samples, ...)
    data: pd.DataFrame                # columns: time, <channel names>
    n_segments: int
    comments: list[str] = field(default_factory=list)

    @property
    def dt(self) -> float:
        return float(self.segment_header["Delta_X"][0])

    @property
    def fs(self) -> float:
        return 1.0 / self.dt

    @property
    def channels(self) -> list[str]:
        return [c for c in self.data.columns if c != "time"]

    def __repr__(self) -> str:
        return (f"LVMData({self.path.name!r}, channels={self.channels}, "
                f"samples={len(self.data)}, segments={self.n_segments}, fs={self.fs:g} Hz)")


def is_lvm(path: str | Path) -> bool:
    """True if the file starts with the LabVIEW Measurement signature (extension is irrelevant)."""
    try:
        with open(path, "rb") as f:
            return f.read(19) == b"LabVIEW Measurement"
    except OSError:
        return False


def _parse_kv(line: str) -> tuple[str, list[str]]:
    parts = line.rstrip("\r\n").split("\t")
    return parts[0], [p for p in parts[1:] if p != ""]


def read_lvm(path: str | Path) -> LVMData:
    """Read a (possibly multi-segment) LVM file into an :class:`LVMData`."""
    path = Path(path)
    header: dict[str, str] = {}
    seg_header: dict[str, list] = {}
    names: list[str] | None = None
    comments: list[str] = []
    n_headers = 0
    numeric = io.StringIO()

    with open(path, "r", encoding="latin-1") as f:
        for line in f:
            if not line or line[0] in _NUMERIC_START:
                numeric.write(line)
                continue
            if line.startswith(END_OF_HEADER):
                n_headers += 1
                continue
            key, vals = _parse_kv(line)
            if key == "X_Value":
                if names is None:
                    names = vals
                continue
            if n_headers == 0:
                header[key] = vals[0] if vals else ""
            elif n_headers == 1 and key:
                seg_header[key] = vals
            # comments on data rows start with a digit so they are handled below

    if names is None:
        raise ValueError(f"{path} does not look like a multi-heading LVM file")

    # names = ['Phase1', 'X_Value', 'Phase2', ..., 'PD', 'Comment'] -> Y channel names
    ch_names = [n for n in names if n not in ("X_Value", "Comment")]
    n_cols = 2 * len(ch_names)
    numeric.seek(0)
    raw = pd.read_csv(numeric, sep="\t", header=None,
                      names=list(range(n_cols + 1)), na_values=["NaN"],
                      dtype={n_cols: str}, engine="c")
    comments = raw[n_cols].dropna().tolist()
    raw = raw.iloc[:, :n_cols].to_numpy(dtype=float)

    t = raw[:, 0].copy()
    dt = float(seg_header.get("Delta_X", [np.median(np.diff(t[:100]))])[0])
    # make time continuous across segments (time restarts at 0 in each segment)
    starts = np.flatnonzero(np.diff(t) < 0) + 1
    offset = np.zeros_like(t)
    for s in starts:
        offset[s:] = t[s - 1] + offset[s - 1] + dt - t[s]
    t += offset

    df = pd.DataFrame({"time": t})
    for i, name in enumerate(ch_names):
        df[name] = raw[:, 2 * i + 1]

    return LVMData(path=path, header=header, segment_header=seg_header, data=df,
                   n_segments=len(starts) + 1, comments=comments)


# --------------------------------------------------------------------------- analysis

def unwrap_phase(phase: np.ndarray, period: float = 360.0, discont: float | None = None) -> np.ndarray:
    """Unwrap a phase signal. Default assumes degrees (period 360). NaNs are preserved."""
    out = np.array(phase, dtype=float)
    ok = ~np.isnan(out)
    out[ok] = np.unwrap(out[ok], discont=discont, period=period)
    return out


def lowpass_trend(y: np.ndarray, cutoff: float = 1e-4, order: int = 12) -> np.ndarray:
    """Zero-phase Butterworth low-pass (MATLAB: designfilt lowpassiir + filtfilt).

    ``cutoff`` is normalised to Nyquist (MATLAB HalfPowerFrequency). Second-order
    sections are used because a 12th-order filter at 1e-4 is unstable in b/a form.
    """
    sos = signal.butter(order, cutoff, output="sos")
    return signal.sosfiltfilt(sos, y)


def periodogram(y: np.ndarray, fs: float) -> tuple[np.ndarray, np.ndarray]:
    """One-sided PSD via FFT, identical to the MATLAB code: psd = |X|^2 / (fs N), doubled except DC/Nyquist."""
    y = np.asarray(y, dtype=float)
    n = len(y)
    xdft = np.fft.rfft(y)
    psd = np.abs(xdft) ** 2 / (fs * n)
    psd[1:-1 if n % 2 == 0 else None] *= 2
    freq = np.fft.rfftfreq(n, d=1 / fs)
    return freq, psd


def min_block_std(y: np.ndarray, block: int = 100) -> float:
    """Smallest standard deviation over consecutive blocks of ``block`` samples (MATLAB ``minstd``)."""
    y = np.asarray(y, dtype=float)
    nb = len(y) // block
    if nb == 0:
        return float(np.std(y, ddof=1))
    return float(np.std(y[: nb * block].reshape(nb, block), axis=1, ddof=1).min())


def fs_from_filename(name: str, default: float = 1000.0) -> float:
    """MATLAB logic: '<d>kHz' in the name -> d*1000 Hz ('0kHz' means 10 kHz), else ``default``."""
    m = re.search(r"(\d)kHz", name)
    if not m:
        return default
    d = int(m.group(1))
    return 10_000.0 if d == 0 else d * 1000.0


def decimate_for_plot(x: np.ndarray, y: np.ndarray, max_points: int = 50_000):
    """Min/max decimation so long traces plot fast without hiding spikes.

    Returns float32 arrays: half the size when embedded in an HTML plot, and still
    ~1 ms time resolution up to t = 8192 s (2 ms up to 4.5 h), far finer
    than the spacing of decimated points.
    """
    x, y = np.asarray(x), np.asarray(y)
    n = len(x)
    if n <= max_points:
        return x.astype(np.float32), y.astype(np.float32)
    step = int(np.ceil(2 * n / max_points))
    m = n // step
    xs = x[: m * step].reshape(m, step)
    ys = y[: m * step].reshape(m, step)
    i_min = np.nanargmin(np.where(np.isnan(ys), np.inf, ys), axis=1)
    i_max = np.nanargmax(np.where(np.isnan(ys), -np.inf, ys), axis=1)
    rows = np.arange(m)
    idx = np.sort(np.stack([i_min, i_max], axis=1), axis=1)
    return (xs[rows[:, None], idx].ravel().astype(np.float32),
            ys[rows[:, None], idx].ravel().astype(np.float32))


def log_bin(freq: np.ndarray, psd: np.ndarray, bins_per_decade: int = 50):
    """Average a PSD in log-spaced frequency bins (for readable overlays of long periodograms)."""
    f, p = freq[1:], psd[1:]                      # drop DC
    edges = np.logspace(np.log10(f[0]), np.log10(f[-1]),
                        int(np.log10(f[-1] / f[0]) * bins_per_decade) + 1)
    idx = np.digitize(f, edges)
    sums = np.bincount(idx, weights=p, minlength=len(edges) + 1)
    counts = np.bincount(idx, minlength=len(edges) + 1)
    fsum = np.bincount(idx, weights=f, minlength=len(edges) + 1)
    ok = counts > 0
    return fsum[ok] / counts[ok], sums[ok] / counts[ok]
