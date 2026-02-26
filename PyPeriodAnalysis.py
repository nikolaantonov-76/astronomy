"""
Lomb-Scargle period analysis for AAVSO Extended format files.

Usage:
1. Run: `python astronomy/PyPeriodAnalysis.py`
2. Select an AAVSO Extended file (.txt/.csv).
3. Set period search bounds and optional filter.
4. Click "Analyze" to compute Lomb-Scargle and inspect plots.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tkinter as tk
import tkinter.font as tkfont
from astropy.timeseries import LombScargle
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import filedialog, messagebox, ttk


FOOTER_TEXT = "Nikola Antonov (nikola.antonov@iaps.institute), https://astro.iaps.instiute"

# Shared visual theme (same as PyAAVSOGenerator / PyAAVSOPlanner)
_BG = "#f0f4f8"
_BORDER = "#b8c8d8"
_ACCENT = "#1f4d7a"
_TEXT = "#2c3e50"


def _apply_theme(root) -> None:
    style = ttk.Style(root)
    try:
        style.theme_use("clam")
    except Exception:
        pass

    base_font = ("Segoe UI", 10)
    bold_font = ("Segoe UI", 10, "bold")
    fixed_font = ("Consolas", 10)
    for fname in (
        "TkDefaultFont",
        "TkTextFont",
        "TkMenuFont",
        "TkHeadingFont",
        "TkCaptionFont",
        "TkTooltipFont",
    ):
        try:
            tkfont.nametofont(fname).configure(family="Segoe UI", size=10)
        except Exception:
            pass
    try:
        tkfont.nametofont("TkFixedFont").configure(family="Consolas", size=10)
    except Exception:
        pass

    root.option_add("*Font", base_font, "interactive")
    root.option_add("*Text.Font", fixed_font, "interactive")
    root.option_add("*Label.Font", base_font, "interactive")
    root.option_add("*Button.Font", base_font, "interactive")
    root.option_add("*Entry.Font", base_font, "interactive")
    root.option_add("*Menu.Font", base_font, "interactive")

    style.configure(".", font=base_font, background=_BG, foreground=_TEXT)
    style.configure("TFrame", background=_BG)
    style.configure("TLabel", font=base_font, background=_BG, foreground=_TEXT)
    style.configure("TCheckbutton", font=base_font, background=_BG, foreground=_TEXT)
    style.configure(
        "TLabelframe",
        background=_BG,
        bordercolor=_BORDER,
        relief="groove",
        borderwidth=1,
    )
    style.configure("TLabelframe.Label", font=bold_font, background=_BG, foreground=_ACCENT)
    style.configure(
        "TButton",
        font=base_font,
        background="#dce8f5",
        foreground="#1a3a5c",
        relief="flat",
        borderwidth=1,
        padding=(8, 4),
    )
    style.map(
        "TButton",
        background=[("active", "#b8d0ea"), ("pressed", "#90b8e0")],
        relief=[("pressed", "flat")],
    )
    style.configure(
        "TEntry",
        font=base_font,
        fieldbackground="white",
        bordercolor=_BORDER,
        lightcolor=_BORDER,
        darkcolor=_BORDER,
    )
    style.configure("TCombobox", font=base_font, fieldbackground="white", bordercolor=_BORDER)
    root.configure(bg=_BG)


def _parse_delim(delim_value: str) -> str:
    token = (delim_value or "").strip().lower()
    if token in {",", "comma"}:
        return ","
    if token in {"tab", r"\t"}:
        return "\t"
    if token in {"space", " "}:
        return " "
    return ","


def _to_float_or_nan(value) -> float:
    try:
        return float(str(value).strip())
    except Exception:
        return math.nan


def read_aavso_extended(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"File not found: {path}")

    header_lines: list[str] = []
    data_lines: list[str] = []
    with p.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                header_lines.append(line)
            else:
                data_lines.append(line)

    if not data_lines:
        raise ValueError("No data rows found.")

    delim = ","
    for line in header_lines:
        if line.upper().startswith("#DELIM="):
            delim = _parse_delim(line.split("=", 1)[1])
            break

    column_line = None
    for line in reversed(header_lines):
        if line.upper().startswith("#NAME"):
            column_line = line[1:]
            break
    if not column_line:
        raise ValueError("Column definition header '#NAME,...' is missing.")

    cols = next(csv.reader([column_line], delimiter=delim))
    rows = list(csv.reader(data_lines, delimiter=delim))
    if not rows:
        raise ValueError("Could not parse any rows.")

    ncols = len(cols)
    clipped = [r[:ncols] + [""] * max(0, ncols - len(r)) for r in rows]
    df = pd.DataFrame(clipped, columns=cols)

    if "DATE" not in df.columns or "MAG" not in df.columns:
        raise ValueError("Required columns DATE and MAG are not present.")

    df["DATE"] = df["DATE"].map(_to_float_or_nan)
    df["MAG_RAW"] = df["MAG"].astype(str).str.strip()
    df["UPPER_LIMIT"] = df["MAG_RAW"].str.startswith("<")
    df["MAG_VALUE"] = df["MAG_RAW"].str.replace("<", "", regex=False).str.replace(">", "", regex=False).map(_to_float_or_nan)
    if "MERR" in df.columns:
        df["MERR_VALUE"] = df["MERR"].map(_to_float_or_nan)
    else:
        df["MERR_VALUE"] = math.nan
    if "FILT" not in df.columns:
        df["FILT"] = ""
    return df


def run_lomb_scargle(
    df: pd.DataFrame,
    filt: str,
    min_period: float,
    max_period: float,
    samples_per_peak: int,
    nyquist_factor: float,
    use_errors: bool,
    exclude_limits: bool,
    ls_method: str,
    ls_normalization: str,
    fit_mean: bool,
    center_data: bool,
    nterms: int,
    fap_method: str,
):
    work = df.copy()
    if filt and filt.upper() != "ALL":
        work = work[work["FILT"].astype(str).str.upper() == filt.upper()]
    if exclude_limits:
        work = work[~work["UPPER_LIMIT"]]

    work = work[np.isfinite(work["DATE"]) & np.isfinite(work["MAG_VALUE"])]
    if len(work) < 5:
        raise ValueError("Not enough usable points for Lomb-Scargle (need at least 5).")

    t = work["DATE"].to_numpy(dtype=float)
    y = work["MAG_VALUE"].to_numpy(dtype=float)
    t_ref = float(np.mean(t))
    t_rel = t - t_ref

    dy = None
    if use_errors:
        merr = work["MERR_VALUE"].to_numpy(dtype=float)
        if np.isfinite(merr).sum() >= 5:
            dy = np.where(np.isfinite(merr) & (merr > 0), merr, np.nanmedian(merr[np.isfinite(merr) & (merr > 0)]))
        else:
            dy = None

    ls = LombScargle(
        t_rel,
        y,
        dy=dy,
        fit_mean=fit_mean,
        center_data=center_data,
        nterms=int(nterms),
        normalization=ls_normalization,
    )
    min_freq = 1.0 / max_period
    max_freq = 1.0 / min_period
    freq, power = ls.autopower(
        minimum_frequency=min_freq,
        maximum_frequency=max_freq,
        samples_per_peak=int(samples_per_peak),
        nyquist_factor=float(nyquist_factor),
        method=ls_method,
    )
    if len(freq) == 0:
        raise ValueError("Frequency grid is empty. Adjust the period bounds.")

    best_idx = int(np.argmax(power))
    best_freq = float(freq[best_idx])
    best_period = 1.0 / best_freq

    t_rel_start = float(np.min(t_rel))
    model_t_rel_one_cycle = t_rel_start + np.linspace(0.0, best_period, 2000)
    model_mag_one_cycle = ls.model(model_t_rel_one_cycle, best_freq)
    t0_rel = float(model_t_rel_one_cycle[int(np.argmin(model_mag_one_cycle))])
    t0_phase = t0_rel + t_ref

    phase = ((t_rel - t0_rel) * best_freq) % 1.0
    order = np.argsort(phase)
    model_phase = np.linspace(0.0, 1.0, 400)
    model_t_rel = t0_rel + (model_phase / best_freq)
    model_mag = ls.model(model_t_rel, best_freq)

    period = 1.0 / freq
    per_order = np.argsort(period)
    period_sorted = period[per_order]
    power_sorted = power[per_order]

    try:
        fap = float(ls.false_alarm_probability(float(power[best_idx]), method=fap_method))
    except Exception:
        fap = math.nan
    fap_levels = {}
    for alpha in (0.1, 0.05, 0.01):
        try:
            fap_levels[alpha] = float(ls.false_alarm_level(alpha, method=fap_method))
        except Exception:
            pass

    return {
        "work": work,
        "t": t,
        "t_ref": t_ref,
        "y": y,
        "freq": freq,
        "period": period,
        "power": power,
        "period_sorted": period_sorted,
        "power_sorted": power_sorted,
        "best_freq": best_freq,
        "best_period": best_period,
        "best_period_min": best_period * 24.0 * 60.0,
        "t0_phase": t0_phase,
        "fap": fap,
        "fap_levels": fap_levels,
        "phase": phase[order],
        "y_phase": y[order],
        "model_phase": model_phase,
        "model_mag": model_mag,
    }


class PeriodAnalysisGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("AAVSO Lomb-Scargle Period Analysis")
        _apply_theme(root)

        self.vars = {
            "input_file": tk.StringVar(value=""),
            "filter": tk.StringVar(value="ALL"),
            "min_period": tk.StringVar(value="0.03"),
            "max_period": tk.StringVar(value="10.0"),
            "samples_per_peak": tk.StringVar(value="8"),
            "nyquist_factor": tk.StringVar(value="5"),
            "ls_method": tk.StringVar(value="auto"),
            "ls_normalization": tk.StringVar(value="standard"),
            "fap_method": tk.StringVar(value="baluev"),
            "fit_mean": tk.BooleanVar(value=True),
            "center_data": tk.BooleanVar(value=True),
            "nterms": tk.StringVar(value="1"),
            "use_errors": tk.BooleanVar(value=True),
            "exclude_limits": tk.BooleanVar(value=True),
        }
        self.df = None
        self._build_ui()
        self._fit_window()
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

    def _fit_window(self):
        self.root.update_idletasks()
        width = min(max(self.root.winfo_reqwidth(), 980), 1240)
        height = min(max(self.root.winfo_reqheight(), 700), 900)
        self.root.geometry(f"{width}x{height}")
        self.root.minsize(960, 680)

    def _build_ui(self):
        ttk.Label(
            self.root,
            text=FOOTER_TEXT,
            anchor="center",
            justify="center",
            foreground="#888888",
        ).pack(side="bottom", fill="x", padx=14, pady=(0, 6))

        action_bar = ttk.Frame(self.root)
        action_bar.pack(side="bottom", fill="x", padx=14, pady=(0, 8))
        ttk.Button(action_bar, text="Analyze", command=self._analyze).pack(side="right")

        outer = ttk.Frame(self.root)
        outer.pack(fill="both", expand=True, padx=14, pady=(14, 8))

        controls = ttk.LabelFrame(outer, text="Input and Search Settings", padding=10)
        controls.pack(fill="x")

        pad = {"padx": 8, "pady": 5}
        ttk.Label(controls, text="AAVSO Extended File").grid(row=0, column=0, sticky="e", **pad)
        ttk.Entry(controls, textvariable=self.vars["input_file"], width=70).grid(row=0, column=1, sticky="we", **pad)
        ttk.Button(controls, text="Browse", command=self._browse).grid(row=0, column=2, sticky="w", **pad)

        ttk.Label(controls, text="Filter").grid(row=1, column=0, sticky="e", **pad)
        self.filter_combo = ttk.Combobox(
            controls,
            textvariable=self.vars["filter"],
            values=["ALL"],
            width=16,
            state="readonly",
        )
        self.filter_combo.grid(row=1, column=1, sticky="w", **pad)

        ttk.Label(controls, text="Min Period (days)").grid(row=2, column=0, sticky="e", **pad)
        ttk.Entry(controls, textvariable=self.vars["min_period"], width=18).grid(row=2, column=1, sticky="w", **pad)
        ttk.Label(controls, text="Max Period (days)").grid(row=2, column=2, sticky="e", **pad)
        ttk.Entry(controls, textvariable=self.vars["max_period"], width=18).grid(row=2, column=3, sticky="w", **pad)

        ttk.Label(controls, text="Samples/Peak").grid(row=3, column=0, sticky="e", **pad)
        ttk.Entry(controls, textvariable=self.vars["samples_per_peak"], width=18).grid(row=3, column=1, sticky="w", **pad)
        ttk.Label(controls, text="Nyquist Factor").grid(row=3, column=2, sticky="e", **pad)
        ttk.Entry(controls, textvariable=self.vars["nyquist_factor"], width=18).grid(row=3, column=3, sticky="w", **pad)

        ttk.Label(controls, text="LS Method").grid(row=4, column=0, sticky="e", **pad)
        ttk.Combobox(
            controls,
            textvariable=self.vars["ls_method"],
            values=["auto", "fast", "slow", "cython", "chi2", "fastchi2", "scipy"],
            width=16,
            state="readonly",
        ).grid(row=4, column=1, sticky="w", **pad)
        ttk.Label(controls, text="Normalization").grid(row=4, column=2, sticky="e", **pad)
        ttk.Combobox(
            controls,
            textvariable=self.vars["ls_normalization"],
            values=["standard", "model", "log", "psd"],
            width=16,
            state="readonly",
        ).grid(row=4, column=3, sticky="w", **pad)

        ttk.Label(controls, text="FAP Method").grid(row=5, column=0, sticky="e", **pad)
        ttk.Combobox(
            controls,
            textvariable=self.vars["fap_method"],
            values=["baluev", "davies", "naive", "bootstrap"],
            width=16,
            state="readonly",
        ).grid(row=5, column=1, sticky="w", **pad)
        ttk.Label(controls, text="nterms").grid(row=5, column=2, sticky="e", **pad)
        ttk.Entry(controls, textvariable=self.vars["nterms"], width=18).grid(row=5, column=3, sticky="w", **pad)

        ttk.Checkbutton(controls, text="Use MERR as uncertainty", variable=self.vars["use_errors"]).grid(
            row=6, column=0, columnspan=2, sticky="w", padx=8, pady=5
        )
        ttk.Checkbutton(controls, text="fit_mean", variable=self.vars["fit_mean"]).grid(
            row=6, column=2, sticky="w", padx=8, pady=5
        )
        ttk.Checkbutton(controls, text="center_data", variable=self.vars["center_data"]).grid(
            row=6, column=3, sticky="w", padx=8, pady=5
        )
        ttk.Checkbutton(controls, text="Exclude upper-limit points (<MAG)", variable=self.vars["exclude_limits"]).grid(
            row=7, column=1, columnspan=3, sticky="w", padx=8, pady=(2, 8)
        )

        controls.columnconfigure(1, weight=1)

        self.summary_label = ttk.Label(
            outer,
            text="Load a file, then run analysis.",
            foreground="#304860",
        )
        self.summary_label.pack(fill="x", pady=(8, 4))

        plot_frame = ttk.LabelFrame(outer, text="Lomb-Scargle Results", padding=8)
        plot_frame.pack(fill="both", expand=True)

        self.fig, self.axes = plt.subplots(2, 1, figsize=(10, 7), dpi=100)
        self.fig.patch.set_facecolor("#f8fbff")
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self._draw_placeholder()
        self.root.bind("<Return>", lambda _event: self._analyze())

    def _draw_placeholder(self):
        ax1, ax2 = self.axes
        for ax in self.axes:
            ax.clear()
            ax.grid(True, alpha=0.3)
        ax1.set_title("Periodogram")
        ax1.set_xlabel("Period (days)")
        ax1.set_ylabel("Power")
        ax2.set_title("Phase-folded Light Curve")
        ax2.set_xlabel("Phase")
        ax2.set_ylabel("Magnitude")
        ax2.invert_yaxis()
        self.fig.tight_layout()
        self.canvas.draw_idle()

    def _browse(self):
        selected = filedialog.askopenfilename(
            title="Select AAVSO Extended file",
            filetypes=[("Text/CSV files", "*.txt *.csv"), ("All files", "*.*")],
        )
        if not selected:
            return
        self.vars["input_file"].set(selected)
        try:
            self.df = read_aavso_extended(selected)
            filters = sorted({str(v).strip().upper() for v in self.df["FILT"] if str(v).strip()})
            combo_vals = ["ALL"] + filters
            self.filter_combo["values"] = combo_vals
            self.vars["filter"].set("ALL")
            self.summary_label.config(text=f"Loaded {len(self.df)} rows. Available filters: {', '.join(filters) if filters else 'None'}")
        except Exception as exc:
            self.df = None
            self.filter_combo["values"] = ["ALL"]
            self.vars["filter"].set("ALL")
            messagebox.showerror("Load Error", str(exc))

    def _analyze(self):
        try:
            if self.df is None:
                input_path = self.vars["input_file"].get().strip()
                if not input_path:
                    raise ValueError("Select an input file first.")
                self.df = read_aavso_extended(input_path)

            min_period = float(self.vars["min_period"].get().strip())
            max_period = float(self.vars["max_period"].get().strip())
            spp = int(float(self.vars["samples_per_peak"].get().strip()))
            nyquist_factor = float(self.vars["nyquist_factor"].get().strip())
            nterms = int(float(self.vars["nterms"].get().strip()))
            if min_period <= 0 or max_period <= 0:
                raise ValueError("Period bounds must be positive.")
            if min_period >= max_period:
                raise ValueError("Min period must be smaller than max period.")
            if spp < 2:
                raise ValueError("Samples/Peak should be at least 2.")
            if nyquist_factor <= 0:
                raise ValueError("Nyquist factor must be positive.")
            if nterms < 1:
                raise ValueError("nterms must be >= 1.")

            res = run_lomb_scargle(
                df=self.df,
                filt=self.vars["filter"].get().strip(),
                min_period=min_period,
                max_period=max_period,
                samples_per_peak=spp,
                nyquist_factor=nyquist_factor,
                use_errors=bool(self.vars["use_errors"].get()),
                exclude_limits=bool(self.vars["exclude_limits"].get()),
                ls_method=self.vars["ls_method"].get().strip(),
                ls_normalization=self.vars["ls_normalization"].get().strip(),
                fit_mean=bool(self.vars["fit_mean"].get()),
                center_data=bool(self.vars["center_data"].get()),
                nterms=nterms,
                fap_method=self.vars["fap_method"].get().strip(),
            )
            self._plot_results(res)
            self.summary_label.config(
                text=(
                    f"Best period: {res['best_period']:.8f} d  "
                    f"({res['best_period_min']:.2f} min), "
                    f"T0(phase=0): JD {res['t0_phase']:.6f}, "
                    f"FAP={res['fap']:.3e}, points={len(res['work'])}"
                )
            )
        except Exception as exc:
            messagebox.showerror("Analysis Error", str(exc))

    def _plot_results(self, res):
        ax1, ax2 = self.axes
        ax1.clear()
        ax2.clear()

        ax1.plot(res["period_sorted"], res["power_sorted"], color="#1f4d7a", lw=1.2)
        ax1.axvline(res["best_period"], color="#d64545", ls="--", lw=1.0)
        x_right = float(np.nanmax(res["period_sorted"]))
        for alpha, level in sorted(res["fap_levels"].items(), key=lambda x: x[0], reverse=True):
            ax1.axhline(level, color="#8a8a8a", ls=":", lw=0.9)
            ax1.text(
                x_right,
                level,
                f"FAP {alpha:.0%}",
                ha="right",
                va="bottom",
                color="#666666",
                fontsize=8,
            )
        ax1.set_title("Lomb-Scargle Periodogram")
        ax1.set_xlabel("Period (days)")
        ax1.set_ylabel("Power")
        ax1.grid(True, alpha=0.3)

        phase2 = np.concatenate([res["phase"], res["phase"] + 1.0])
        y2 = np.concatenate([res["y_phase"], res["y_phase"]])
        model_phase2 = np.concatenate([res["model_phase"], res["model_phase"] + 1.0])
        model_mag2 = np.concatenate([res["model_mag"], res["model_mag"]])
        ax2.scatter(phase2, y2, s=18, color="#2e7d32", alpha=0.8, edgecolors="none")
        ax2.plot(model_phase2, model_mag2, color="#9c2f2f", lw=1.2)
        ax2.set_xlim(0.0, 2.0)
        ax2.set_title(
            f"Phased Light Curve (P = {res['best_period']:.8f} d = {res['best_period_min']:.2f} min, "
            f"T0 = JD {res['t0_phase']:.6f})"
        )
        ax2.set_xlabel("Phase")
        ax2.set_ylabel("Magnitude")
        ax2.grid(True, alpha=0.3)
        ax2.invert_yaxis()

        self.fig.tight_layout()
        self.canvas.draw_idle()

    def _on_close(self):
        try:
            plt.close(self.fig)
            plt.close("all")
        except Exception:
            pass
        try:
            self.root.quit()
        except Exception:
            pass
        self.root.destroy()


def main():
    root = tk.Tk()
    PeriodAnalysisGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
