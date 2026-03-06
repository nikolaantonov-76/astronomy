"""
Lomb-Scargle period analysis for AAVSO Extended format files.

Usage:
1. Run: `python astronomy/PyPeriodAnalysis.py`
2. Select an AAVSO Extended file (.txt/.csv).
3. Set period search bounds and optional filter.
4. Click "Analyze" to compute Lomb-Scargle and inspect plots.
"""

from __future__ import annotations

import calendar
import csv
import datetime as dt
import json
import math
import re
from pathlib import Path
from urllib.parse import urljoin

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import tkinter as tk
import tkinter.font as tkfont
from astropy.timeseries import LombScargle
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Rectangle
from tkinter import filedialog, messagebox, ttk


FOOTER_TEXT = (
    "Nikola Antonov (nikola.antonov@iaps.institute), https://astro.iaps.institute\n"
    "0.25-m Newtonian f/4.8, ASI553MM Pro CMOS\n"
    "Meshtitsa, Bulgaria"
)

# Shared visual theme (same as PyAAVSOGenerator / PyAAVSOPlanner)
_BG = "#f0f4f8"
_BORDER = "#b8c8d8"
_ACCENT = "#1f4d7a"
_TEXT = "#2c3e50"
AAVSO_WEBOBS_URL = "https://apps.aavso.org/webobs/results/"

# Approximate spectral-color mapping for common AAVSO filters.
AAVSO_FILTER_COLORS = {
    "U": "#5b2ca0",     # ultraviolet
    "B": "#1f5fd0",     # blue
    "V": "#2e8b57",     # visual (green)
    "R": "#c43a31",     # red
    "I": "#6a1b1a",     # near-IR
    "CV": "#333333",    # clear visual
    "C": "#444444",     # clear
    "TG": "#1c8e5f",
    "TB": "#2f65c8",
    "TR": "#c14135",
    "SG": "#2ca25f",
    "SR": "#de2d26",
    "SI": "#7f0000",
    "SZ": "#4a1486",
    "HA": "#e41a1c",
    "OIII": "#17becf",
    "VIS": "#2e8b57",
}


def _axis_limits_with_margin(values, frac=0.08, min_pad=1e-6):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return None
    vmin = float(np.min(arr))
    vmax = float(np.max(arr))
    span = vmax - vmin
    if span <= 0:
        pad = max(abs(vmin) * frac, min_pad)
    else:
        pad = max(span * frac, min_pad)
    return vmin - pad, vmax + pad


def _jd_to_calendar_date_text(jd_value) -> str:
    try:
        jd_float = float(jd_value)
        dt_value = dt.datetime(1970, 1, 1) + dt.timedelta(days=jd_float - 2440588.0)
        return dt_value.strftime("%Y-%m-%d")
    except Exception:
        return "n/a"


def _set_optimal_xy_limits(ax, x_values, y_values, x_frac=0.04, y_frac=0.08):
    x_lim = _axis_limits_with_margin(x_values, frac=x_frac)
    y_lim = _axis_limits_with_margin(y_values, frac=y_frac)
    if x_lim is not None:
        ax.set_xlim(*x_lim)
    if y_lim is not None:
        ax.set_ylim(*y_lim)


def _parse_filter_selection(value: str):
    raw = str(value or "").strip().upper()
    if not raw or raw in {"ALL", "*"}:
        return []
    parts = re.split(r"[,\s;+/|]+", raw)
    tokens = [p.strip().upper() for p in parts if p.strip()]
    # Preserve order and uniqueness.
    out = []
    seen = set()
    for t in tokens:
        if t in seen:
            continue
        seen.add(t)
        out.append(t)
    return out


def _strip_html_tags(raw: str) -> str:
    text = re.sub(r"<[^>]+>", " ", raw or "", flags=re.S)
    return re.sub(r"\s+", " ", text).strip()


def _to_float_or_none(value):
    try:
        return float(str(value).strip())
    except Exception:
        return None


def _is_jd_like(value: str) -> bool:
    jd = _to_float_or_none(value)
    return jd is not None and 2400000.0 <= jd <= 2600000.0


def _is_mag_like(value: str) -> bool:
    token = str(value).strip().replace("<", "").replace(">", "")
    mag = _to_float_or_none(token)
    return mag is not None and -5.0 <= mag <= 30.0


def _is_err_like(value: str) -> bool:
    token = str(value).strip().replace("<", "").replace(">", "")
    err = _to_float_or_none(token)
    return err is not None and 0.0 <= err <= 2.0


def _find_jd_in_cells(values):
    for cell in values:
        if _is_jd_like(cell):
            return str(cell).strip()
    return ""


def _find_mag_in_cells(values):
    for cell in values:
        if _is_mag_like(cell):
            return str(cell).strip()
    return ""


def _find_band_in_cells(values):
    allowed = {
        "U", "B", "V", "R", "I", "CV", "C", "TG", "TB", "TR", "SG", "SR", "SI", "SZ", "HA", "OIII"
    }
    for cell in values:
        token = str(cell).strip().upper()
        if token in allowed:
            return token
    return ""


def _find_observer_in_cells(values):
    for cell in values:
        token = str(cell).strip().upper()
        if re.fullmatch(r"[A-Z0-9]{3,8}", token):
            return token
    return ""


def _find_err_in_cells(values):
    for cell in values:
        if _is_err_like(cell):
            return str(cell).strip()
    return ""


def _normalize_observer_code(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9]+", "", str(value or "").upper())
    return token


def _extract_rows_from_next_data(page_html: str):
    script_match = re.search(
        r'<script[^>]*id=["\']__NEXT_DATA__["\'][^>]*>(.*?)</script>',
        page_html or "",
        flags=re.I | re.S,
    )
    if not script_match:
        return []

    try:
        payload = json.loads(script_match.group(1).strip())
    except Exception:
        return []

    found = []

    def walk(node):
        if isinstance(node, dict):
            lowered = {str(k).lower(): v for k, v in node.items()}
            jd_val = lowered.get("jd", "")
            mag_val = lowered.get("magnitude") or lowered.get("mag") or lowered.get("value") or ""
            err_val = (
                lowered.get("merr")
                or lowered.get("error")
                or lowered.get("uncertainty")
                or lowered.get("magerror")
                or ""
            )
            if _is_jd_like(str(jd_val)) and _is_mag_like(str(mag_val)):
                found.append(
                    {
                        "date": str(lowered.get("date", "") or lowered.get("obs_date", "") or ""),
                        "jd": str(jd_val),
                        "mag": str(mag_val),
                        "merr": str(err_val),
                        "band": str(lowered.get("band", "") or lowered.get("filter", "") or "").upper(),
                        "observer": str(lowered.get("observer", "") or lowered.get("obscode", "") or "").upper(),
                    }
                )
            for value in node.values():
                walk(value)
            return
        if isinstance(node, list):
            for item in node:
                walk(item)

    walk(payload)
    return found


def _extract_rows_from_html(page_html: str):
    rows_out = _extract_rows_from_next_data(page_html)
    tables = re.findall(r"<table\b.*?>.*?</table>", page_html or "", flags=re.I | re.S)
    for table in tables:
        rows = re.findall(r"<tr\b.*?>.*?</tr>", table, flags=re.I | re.S)
        if len(rows) < 2:
            continue

        header_cells = re.findall(r"<th\b.*?>(.*?)</th>", rows[0], flags=re.I | re.S)
        headers = [_strip_html_tags(cell).lower() for cell in header_cells]

        idx_jd = next((i for i, h in enumerate(headers) if h == "jd"), None)
        idx_mag = next((i for i, h in enumerate(headers) if "magnitude" in h or h == "mag"), None)
        idx_err = next((i for i, h in enumerate(headers) if "error" in h or "err" in h or "uncert" in h), None)
        idx_band = next((i for i, h in enumerate(headers) if "band" in h or "filter" in h), None)
        idx_obs = next((i for i, h in enumerate(headers) if "observer" in h or "obs" in h), None)
        idx_date = next((i for i, h in enumerate(headers) if "date" in h), None)
        if idx_jd is None or idx_mag is None:
            continue

        for row in rows[1:]:
            cells = re.findall(r"<td\b.*?>(.*?)</td>", row, flags=re.I | re.S)
            if not cells:
                continue
            values = [_strip_html_tags(cell) for cell in cells]
            jd_val = values[idx_jd] if idx_jd < len(values) else _find_jd_in_cells(values)
            mag_val = values[idx_mag] if idx_mag < len(values) else _find_mag_in_cells(values)
            if not _is_jd_like(jd_val) or not _is_mag_like(mag_val):
                continue
            rows_out.append(
                {
                    "date": values[idx_date] if idx_date is not None and idx_date < len(values) else "",
                    "jd": jd_val,
                    "mag": mag_val,
                    "merr": values[idx_err] if idx_err is not None and idx_err < len(values) else _find_err_in_cells(values),
                    "band": values[idx_band] if idx_band is not None and idx_band < len(values) else _find_band_in_cells(values),
                    "observer": values[idx_obs] if idx_obs is not None and idx_obs < len(values) else _find_observer_in_cells(values),
                }
            )

    loose_rows = re.findall(r"<tr\b.*?>.*?</tr>", page_html or "", flags=re.I | re.S)
    for row in loose_rows:
        cells = re.findall(r"<t[dh]\b.*?>(.*?)</t[dh]>", row, flags=re.I | re.S)
        if len(cells) < 3:
            continue
        values = [_strip_html_tags(cell) for cell in cells]
        jd_val = _find_jd_in_cells(values)
        mag_val = _find_mag_in_cells(values)
        if not _is_jd_like(jd_val) or not _is_mag_like(mag_val):
            continue
        rows_out.append(
            {
                "date": "",
                "jd": jd_val,
                "mag": mag_val,
                "merr": _find_err_in_cells(values),
                "band": _find_band_in_cells(values),
                "observer": _find_observer_in_cells(values),
            }
        )

    unique = {}
    for row in rows_out:
        key = (
            str(row.get("jd", "")).strip(),
            str(row.get("mag", "")).strip(),
            str(row.get("merr", "")).strip(),
            str(row.get("band", "")).strip().upper(),
            str(row.get("observer", "")).strip().upper(),
        )
        unique[key] = row

    return sorted(unique.values(), key=lambda item: float(item["jd"]))


def _extract_next_page_url(page_html: str, current_url: str) -> str:
    if not page_html:
        return ""
    # Prefer explicit rel=next links when present.
    rel_next = re.search(
        r'<a[^>]*rel=["\']next["\'][^>]*href=["\']([^"\']+)["\']',
        page_html,
        flags=re.I | re.S,
    )
    if rel_next:
        return urljoin(current_url, rel_next.group(1).strip())

    # Fallback: links with visible "Next" text.
    next_text = re.search(
        r'<a[^>]*href=["\']([^"\']+)["\'][^>]*>\s*(?:Next|›|»)\s*</a>',
        page_html,
        flags=re.I | re.S,
    )
    if next_text:
        return urljoin(current_url, next_text.group(1).strip())
    return ""


def _utc_day_to_jd_bounds(start_date: str, end_date: str):
    start_day = dt.date.fromisoformat(start_date)
    end_day = dt.date.fromisoformat(end_date)
    if end_day < start_day:
        raise ValueError("End date must be the same or after start date.")
    start_ts = pd.Timestamp(start_day.isoformat(), tz="UTC")
    end_ts = pd.Timestamp((end_day + dt.timedelta(days=1)).isoformat(), tz="UTC")
    return float(start_ts.to_julian_date()), float(end_ts.to_julian_date())


def fetch_aavso_rows(
    object_name: str,
    limit: int = 5000,
    observer: str = "",
    filt: str = "",
    start_date: str = "",
    end_date: str = "",
):
    object_name = object_name.strip()
    if not object_name:
        raise ValueError("Object name is required for AAVSO fetch.")

    headers = {
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
        "User-Agent": (
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/120.0.0.0 Safari/537.36"
        ),
    }
    page_size = 200
    target_limit = max(200, int(limit))
    obs_norm = _normalize_observer_code(observer)
    filt_norm = (filt or "").strip().upper()
    if filt_norm in {"ALL", "*"}:
        filt_norm = ""
    jd_start = -float("inf")
    jd_end = float("inf")
    if start_date and end_date:
        try:
            jd_start, jd_end = _utc_day_to_jd_bounds(start_date, end_date)
        except Exception:
            jd_start, jd_end = -float("inf"), float("inf")

    last_error = "No response from AAVSO."
    for key_name in ("star", "target"):
        try:
            base_params = {key_name: object_name, "num_results": str(page_size)}
            # IMPORTANT: WebObs uses obscode for observer filtering.
            if obs_norm:
                base_params["obscode"] = obs_norm
            # Optional hint; we still filter locally to guarantee correctness.
            if filt_norm:
                base_params["band"] = filt_norm

            rows_agg = []
            seen = set()
            max_pages = max(25, int(math.ceil(target_limit / page_size)) + 10)
            repeated_signature_count = 0
            last_signature = None

            for page_no in range(1, max_pages + 1):
                params = dict(base_params)
                params["page"] = str(page_no)
                resp = requests.get(AAVSO_WEBOBS_URL, params=params, timeout=(10, 40), headers=headers)
                resp.raise_for_status()
                rows_page = _extract_rows_from_html(resp.text)
                if not rows_page:
                    break

                jd_values = []
                for row in rows_page:
                    jd = _to_float_or_none(row.get("jd", ""))
                    if jd is not None:
                        jd_values.append(jd)

                if jd_values:
                    sig = (round(min(jd_values), 6), round(max(jd_values), 6), len(jd_values))
                    if sig == last_signature:
                        repeated_signature_count += 1
                    else:
                        repeated_signature_count = 0
                    last_signature = sig
                    if repeated_signature_count >= 2:
                        break

                for row in rows_page:
                    key = (
                        str(row.get("jd", "")).strip(),
                        str(row.get("mag", "")).strip(),
                        str(row.get("merr", "")).strip(),
                        str(row.get("band", "")).strip().upper(),
                        _normalize_observer_code(row.get("observer", "")),
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                    rows_agg.append(row)

                # Stop once we are already fully older than requested start date.
                if jd_values and max(jd_values) < jd_start:
                    break
                if len(rows_agg) >= target_limit:
                    break

            if rows_agg:
                rows_agg.sort(key=lambda item: float(item["jd"]))
                # Keep only requested date interval here to avoid "latest-only" bias.
                trimmed = []
                for row in rows_agg:
                    jd = _to_float_or_none(row.get("jd", ""))
                    if jd is None:
                        continue
                    if jd_start <= jd < jd_end:
                        trimmed.append(row)
                return trimmed if trimmed else rows_agg[:target_limit]
            last_error = "No rows parsed from AAVSO response."
        except requests.exceptions.Timeout:
            last_error = "AAVSO request timeout."
        except Exception as exc:
            last_error = str(exc)
    raise RuntimeError(f"Could not fetch AAVSO data: {last_error}")


def build_df_from_aavso_rows(rows, filt: str, observer: str, start_date: str, end_date: str, object_name: str):
    filt_tokens = _parse_filter_selection(filt)
    obs_norm = _normalize_observer_code(observer)
    if obs_norm in {"", "ALL", "*"}:
        obs_norm = ""
    jd_start = -float("inf")
    jd_end = float("inf")
    if start_date and end_date:
        jd_start, jd_end = _utc_day_to_jd_bounds(start_date, end_date)

    out_rows = []
    for row in rows:
        jd = _to_float_or_none(row.get("jd", ""))
        if jd is None:
            continue
        if not (jd_start <= jd < jd_end):
            continue

        band = str(row.get("band", "")).strip().upper()
        obs_code = _normalize_observer_code(row.get("observer", ""))
        # If band is missing in a parsed row, keep it to avoid false negatives from partial HTML parsing.
        if filt_tokens and band and band not in filt_tokens:
            continue
        # Observer filter is strict by observer code.
        if obs_norm and obs_code != obs_norm:
            continue

        mag_raw = str(row.get("mag", "")).strip()
        if not mag_raw:
            continue

        out_rows.append(
            {
                "NAME": object_name,
                "DATE": jd,
                "MAG": mag_raw,
                "MERR": str(row.get("merr", "")).strip(),
                "FILT": band,
                "OBSERVER": obs_code,
            }
        )

    if not out_rows:
        return pd.DataFrame()

    df = pd.DataFrame(out_rows)
    df["DATE"] = df["DATE"].map(_to_float_or_nan)
    df["MAG_RAW"] = df["MAG"].astype(str).str.strip()
    df["UPPER_LIMIT"] = df["MAG_RAW"].str.startswith("<")
    df["MAG_VALUE"] = (
        df["MAG_RAW"]
        .str.replace("<", "", regex=False)
        .str.replace(">", "", regex=False)
        .map(_to_float_or_nan)
    )
    df["MERR_VALUE"] = df["MERR"].map(_to_float_or_nan)
    if "FILT" not in df.columns:
        df["FILT"] = ""
    return df


class CalendarDialog(tk.Toplevel):
    def __init__(self, parent, initial_date: dt.date | None = None):
        super().__init__(parent)
        self.title("Select Date")
        self.resizable(False, False)
        self.transient(parent)
        self.grab_set()
        self.result = None

        today = dt.date.today()
        base = initial_date or today
        self.year = base.year
        self.month = base.month
        self._selected_day = base.day

        body = ttk.Frame(self, padding=8)
        body.pack(fill="both", expand=True)

        nav = ttk.Frame(body)
        nav.pack(fill="x", pady=(0, 6))
        ttk.Button(nav, text="<", width=3, command=self._prev_month).pack(side="left")
        self.month_label = ttk.Label(nav, text="", anchor="center")
        self.month_label.pack(side="left", fill="x", expand=True)
        ttk.Button(nav, text=">", width=3, command=self._next_month).pack(side="right")

        self.grid_frame = ttk.Frame(body)
        self.grid_frame.pack()

        actions = ttk.Frame(body)
        actions.pack(fill="x", pady=(8, 0))
        ttk.Button(actions, text="Today", command=self._set_today).pack(side="left")
        ttk.Button(actions, text="Cancel", command=self._cancel).pack(side="right")

        self.bind("<Escape>", lambda _e: self._cancel())
        self._build_grid()

    def _build_grid(self):
        for child in self.grid_frame.winfo_children():
            child.destroy()

        self.month_label.config(text=f"{calendar.month_name[self.month]} {self.year}")
        week_headers = ["Mo", "Tu", "We", "Th", "Fr", "Sa", "Su"]
        for col, txt in enumerate(week_headers):
            ttk.Label(self.grid_frame, text=txt, width=4, anchor="center").grid(row=0, column=col, padx=1, pady=1)

        for row_idx, week in enumerate(calendar.monthcalendar(self.year, self.month), start=1):
            for col_idx, day in enumerate(week):
                if day == 0:
                    ttk.Label(self.grid_frame, text=" ", width=4).grid(row=row_idx, column=col_idx, padx=1, pady=1)
                    continue
                ttk.Button(
                    self.grid_frame,
                    text=str(day),
                    width=4,
                    command=lambda d=day: self._select_day(d),
                ).grid(row=row_idx, column=col_idx, padx=1, pady=1)

    def _prev_month(self):
        if self.month == 1:
            self.month = 12
            self.year -= 1
        else:
            self.month -= 1
        self._build_grid()

    def _next_month(self):
        if self.month == 12:
            self.month = 1
            self.year += 1
        else:
            self.month += 1
        self._build_grid()

    def _select_day(self, day: int):
        self._selected_day = day
        self.result = dt.date(self.year, self.month, self._selected_day)
        self.destroy()

    def _set_today(self):
        today = dt.date.today()
        self.result = today
        self.destroy()

    def _cancel(self):
        self.result = None
        self.destroy()


class FilterSelectDialog(tk.Toplevel):
    def __init__(self, parent, title: str, options: list[str], selected: list[str]):
        super().__init__(parent)
        self.title(title)
        self.resizable(False, False)
        self.transient(parent)
        self.grab_set()
        self.result = None

        selected_set = {str(x).strip().upper() for x in selected if str(x).strip()}
        body = ttk.Frame(self, padding=10)
        body.pack(fill="both", expand=True)
        ttk.Label(body, text="Select one or more filters").pack(anchor="w", pady=(0, 6))

        self._vars = {}
        grid = ttk.Frame(body)
        grid.pack(fill="both", expand=True)
        for idx, opt in enumerate(options):
            var = tk.BooleanVar(value=(opt in selected_set))
            self._vars[opt] = var
            ttk.Checkbutton(grid, text=opt, variable=var).grid(
                row=idx // 6, column=idx % 6, sticky="w", padx=6, pady=3
            )

        actions = ttk.Frame(body)
        actions.pack(fill="x", pady=(8, 0))
        ttk.Button(actions, text="All", command=self._set_all).pack(side="left")
        ttk.Button(actions, text="None", command=self._set_none).pack(side="left", padx=(6, 0))
        ttk.Button(actions, text="Cancel", command=self._cancel).pack(side="right")
        ttk.Button(actions, text="Apply", command=self._apply).pack(side="right", padx=(0, 6))

        self.bind("<Escape>", lambda _e: self._cancel())

    def _set_all(self):
        for v in self._vars.values():
            v.set(True)

    def _set_none(self):
        for v in self._vars.values():
            v.set(False)

    def _apply(self):
        self.result = [k for k, v in self._vars.items() if bool(v.get())]
        self.destroy()

    def _cancel(self):
        self.result = None
        self.destroy()


def _apply_theme(root) -> None:
    style = ttk.Style(root)
    try:
        style.theme_use("clam")
    except Exception:
        pass

    base_font = ("Segoe UI", 12)
    bold_font = ("Segoe UI", 12, "bold")
    fixed_font = ("Consolas", 12)
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
    steps: int,
    use_errors: bool,
    exclude_limits: bool,
    ls_method: str,
    ls_normalization: str,
    fit_mean: bool,
    center_data: bool,
    nterms: int,
    fap_method: str,
    fold_period_override: float | None,
    fold_t0_override: float | None,
):
    work = df.copy()
    filt_tokens = _parse_filter_selection(filt)
    if filt_tokens:
        work = work[work["FILT"].astype(str).str.upper().isin(filt_tokens)]
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
    if steps < 2:
        raise ValueError("Steps must be at least 2.")
    freq_grid = np.linspace(min_freq, max_freq, int(steps), dtype=float)
    used_ls_method = ls_method
    try:
        power = ls.power(freq_grid, method=ls_method)
        freq = freq_grid
    except ModuleNotFoundError as exc:
        if "scipy" not in str(exc).lower():
            raise
        used_ls_method = "fast"
        power = ls.power(freq_grid, method=used_ls_method)
        freq = freq_grid
    if len(freq) == 0:
        raise ValueError("Frequency grid is empty. Adjust the period bounds.")

    best_idx = int(np.argmax(power))
    best_freq = float(freq[best_idx])
    best_period = 1.0 / best_freq

    t_rel_start = float(np.min(t_rel))
    model_t_rel_one_cycle = t_rel_start + np.linspace(0.0, best_period, 2000)
    model_mag_one_cycle = ls.model(model_t_rel_one_cycle, best_freq)
    # Classical convention requested: T0 at minimum brightness
    # (in magnitudes this is the maximum model magnitude).
    t0_idx = int(np.argmax(model_mag_one_cycle))
    t0_rel = float(model_t_rel_one_cycle[t0_idx])
    t0_phase = t0_rel + t_ref

    fold_period = float(fold_period_override) if fold_period_override and fold_period_override > 0 else best_period
    fold_freq = 1.0 / fold_period
    fold_t0 = float(fold_t0_override) if fold_t0_override and np.isfinite(fold_t0_override) else t0_phase
    fold_t0_rel = fold_t0 - t_ref

    # Explicit astronomy fold formula: phase = frac((JD - T0) / P)
    phase = np.mod((t - fold_t0) / fold_period, 1.0)
    order = np.argsort(phase)
    model_phase = np.linspace(0.0, 1.0, 400)
    model_t_rel = fold_t0_rel + (model_phase / fold_freq)
    try:
        model_mag = ls.model(model_t_rel, fold_freq)
    except Exception:
        model_mag = np.full_like(model_phase, np.nan, dtype=float)

    cycles_covered = (float(np.max(t)) - float(np.min(t))) / fold_period if fold_period > 0 else math.nan

    period = 1.0 / freq
    per_order = np.argsort(period)
    period_sorted = period[per_order]
    power_sorted = power[per_order]

    fap_error = ""
    fap_method_used = ""
    fap_candidates = [fap_method] + [m for m in ("baluev", "davies", "naive", "bootstrap") if m != fap_method]
    fap = math.nan
    for method_name in fap_candidates:
        try:
            fap = float(ls.false_alarm_probability(float(power[best_idx]), method=method_name))
            fap_method_used = method_name
            break
        except Exception as exc:
            fap_error = f"{type(exc).__name__}: {exc}"
    fap_levels = {}
    if fap_method_used:
        for alpha in (0.1, 0.05, 0.01):
            try:
                fap_levels[alpha] = float(ls.false_alarm_level(alpha, method=fap_method_used))
            except Exception:
                pass

    fap_text = f"{fap:.3e}" if np.isfinite(fap) else "n/a"

    return {
        "work": work,
        "t": t,
        "t_ref": t_ref,
        "y": y,
        "dy": dy,
        "freq": freq,
        "period": period,
        "power": power,
        "period_sorted": period_sorted,
        "power_sorted": power_sorted,
        "best_freq": best_freq,
        "best_period": best_period,
        "best_period_min": best_period * 24.0 * 60.0,
        "used_ls_method": used_ls_method,
        "t0_phase": t0_phase,
        "fold_period": fold_period,
        "fold_period_min": fold_period * 24.0 * 60.0,
        "fold_t0": fold_t0,
        "cycles_covered": cycles_covered,
        "fap": fap,
        "fap_text": fap_text,
        "fap_method_used": fap_method_used,
        "fap_error": fap_error,
        "fap_levels": fap_levels,
        "phase": phase[order],
        "y_phase": y[order],
        "model_phase": model_phase,
        "model_mag": model_mag,
        "search_min_period": float(min_period),
        "search_max_period": float(max_period),
    }


class PeriodAnalysisGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("AAVSO Lomb-Scargle Period Analysis")
        _apply_theme(root)
        self._scroll_canvas = None
        self._scroll_frame = None

        self.vars = {
            "input_file": tk.StringVar(value=""),
            "periodogram_filter": tk.StringVar(value="ALL"),
            "min_period": tk.StringVar(value="0.03"),
            "max_period": tk.StringVar(value="10.0"),
            "steps": tk.StringVar(value="5000"),
            "ls_method": tk.StringVar(value="auto"),
            "ls_normalization": tk.StringVar(value="standard"),
            "fap_method": tk.StringVar(value="baluev"),
            "fit_mean": tk.BooleanVar(value=True),
            "center_data": tk.BooleanVar(value=True),
            "nterms": tk.StringVar(value="1"),
            "fold_period": tk.StringVar(value=""),
            "fold_t0": tk.StringVar(value=""),
            "use_errors": tk.BooleanVar(value=True),
            "exclude_limits": tk.BooleanVar(value=True),
            "aavso_object": tk.StringVar(value=""),
            "aavso_filter": tk.StringVar(value="ALL"),
            "aavso_observer": tk.StringVar(value=""),
            "aavso_start_date": tk.StringVar(value=(dt.date.today() - dt.timedelta(days=30)).isoformat()),
            "aavso_end_date": tk.StringVar(value=dt.date.today().isoformat()),
            "aavso_all_dates": tk.BooleanVar(value=False),
            "aavso_limit": tk.StringVar(value="5000"),
        }
        self.df = None
        self.last_result = None
        self.aavso_start_entry = None
        self.aavso_start_button = None
        self.aavso_end_entry = None
        self.aavso_end_button = None
        self._build_ui()
        self._fit_window()
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

    def _fit_window(self):
        self.root.update_idletasks()
        width = min(max(self.root.winfo_reqwidth(), 1060), 1360)
        height = min(max(self.root.winfo_reqheight(), 720), 1100)
        self.root.geometry(f"{width}x{height}")
        self.root.minsize(980, 640)

    def _on_scroll_frame_configure(self, _event=None):
        if self._scroll_canvas is None or self._scroll_frame is None:
            return
        self._scroll_canvas.configure(scrollregion=self._scroll_canvas.bbox("all"))

    def _on_scroll_canvas_configure(self, event):
        if self._scroll_canvas is None or self._scroll_frame is None:
            return
        self._scroll_canvas.itemconfigure(self._scroll_window_id, width=event.width)

    def _bind_scroll_events(self):
        self.root.bind_all("<MouseWheel>", self._on_mousewheel, add="+")
        self.root.bind_all("<Button-4>", self._on_mousewheel_linux, add="+")
        self.root.bind_all("<Button-5>", self._on_mousewheel_linux, add="+")

    def _on_mousewheel(self, event):
        if self._scroll_canvas is None:
            return
        delta = int(-event.delta / 120) if event.delta else 0
        if delta != 0:
            self._scroll_canvas.yview_scroll(delta, "units")

    def _on_mousewheel_linux(self, event):
        if self._scroll_canvas is None:
            return
        if event.num == 4:
            self._scroll_canvas.yview_scroll(-1, "units")
        elif event.num == 5:
            self._scroll_canvas.yview_scroll(1, "units")

    def _build_ui(self):
        shell = ttk.Frame(self.root)
        shell.pack(fill="both", expand=True)
        self._scroll_canvas = tk.Canvas(shell, bg=_BG, highlightthickness=0, borderwidth=0)
        y_scroll = ttk.Scrollbar(shell, orient="vertical", command=self._scroll_canvas.yview)
        self._scroll_canvas.configure(yscrollcommand=y_scroll.set)
        self._scroll_canvas.pack(side="left", fill="both", expand=True)
        y_scroll.pack(side="right", fill="y")
        self._scroll_frame = ttk.Frame(self._scroll_canvas)
        self._scroll_window_id = self._scroll_canvas.create_window((0, 0), window=self._scroll_frame, anchor="nw")
        self._scroll_frame.bind("<Configure>", self._on_scroll_frame_configure)
        self._scroll_canvas.bind("<Configure>", self._on_scroll_canvas_configure)
        self._bind_scroll_events()

        ttk.Label(
            self._scroll_frame,
            text=FOOTER_TEXT,
            anchor="center",
            justify="center",
            foreground="#888888",
        ).pack(side="bottom", fill="x", padx=14, pady=(0, 6))

        action_bar = ttk.Frame(self._scroll_frame)
        action_bar.pack(side="bottom", fill="x", padx=14, pady=(0, 8))
        ttk.Button(action_bar, text="Save Outputs", command=self._save_outputs).pack(side="right", padx=(0, 8))
        ttk.Button(action_bar, text="Reset", command=self._reset_loaded_data).pack(side="right", padx=(0, 8))
        ttk.Button(action_bar, text="Analyze", command=self._analyze).pack(side="right")

        outer = ttk.Frame(self._scroll_frame)
        outer.pack(fill="both", expand=True, padx=14, pady=(14, 8))

        controls = ttk.LabelFrame(outer, text="Input and Search Settings", padding=10)
        controls.pack(fill="x")

        pad = {"padx": 8, "pady": 5}
        ttk.Label(controls, text="AAVSO Extended File").grid(row=0, column=0, sticky="e", **pad)
        ttk.Entry(controls, textvariable=self.vars["input_file"], width=70).grid(row=0, column=1, sticky="we", **pad)
        ttk.Button(controls, text="Browse", command=self._browse).grid(row=0, column=2, sticky="w", **pad)

        ttk.Label(controls, text="Periodogram Filter").grid(row=1, column=0, sticky="e", **pad)
        self.period_filter_combo = ttk.Combobox(
            controls,
            textvariable=self.vars["periodogram_filter"],
            values=["ALL"],
            width=16,
            state="readonly",
        )
        self.period_filter_combo.grid(row=1, column=1, sticky="w", **pad)

        ttk.Label(controls, text="Min Period (days)").grid(row=2, column=0, sticky="e", **pad)
        ttk.Entry(controls, textvariable=self.vars["min_period"], width=18).grid(row=2, column=1, sticky="w", **pad)
        ttk.Label(controls, text="Max Period (days)").grid(row=2, column=2, sticky="e", **pad)
        ttk.Entry(controls, textvariable=self.vars["max_period"], width=18).grid(row=2, column=3, sticky="w", **pad)

        ttk.Label(controls, text="Steps").grid(row=3, column=0, sticky="e", **pad)
        ttk.Entry(controls, textvariable=self.vars["steps"], width=18).grid(row=3, column=1, sticky="w", **pad)

        ttk.Label(controls, text="LS Method").grid(row=3, column=2, sticky="e", **pad)
        ttk.Combobox(
            controls,
            textvariable=self.vars["ls_method"],
            values=["auto", "fast", "slow", "cython", "chi2", "fastchi2"],
            width=16,
            state="readonly",
        ).grid(row=3, column=3, sticky="w", **pad)
        ttk.Label(controls, text="Normalization").grid(row=4, column=2, sticky="e", **pad)
        ttk.Combobox(
            controls,
            textvariable=self.vars["ls_normalization"],
            values=["standard", "model", "log", "psd"],
            width=16,
            state="readonly",
        ).grid(row=4, column=3, sticky="w", **pad)

        ttk.Label(controls, text="FAP Method").grid(row=4, column=0, sticky="e", **pad)
        ttk.Combobox(
            controls,
            textvariable=self.vars["fap_method"],
            values=["baluev", "davies", "naive", "bootstrap"],
            width=16,
            state="readonly",
        ).grid(row=4, column=1, sticky="w", **pad)
        ttk.Label(controls, text="nterms").grid(row=5, column=2, sticky="e", **pad)
        ttk.Entry(controls, textvariable=self.vars["nterms"], width=18).grid(row=5, column=3, sticky="w", **pad)

        ttk.Label(controls, text="Fold Period Override (d)").grid(row=6, column=0, sticky="e", **pad)
        ttk.Entry(controls, textvariable=self.vars["fold_period"], width=18).grid(row=6, column=1, sticky="w", **pad)
        ttk.Label(controls, text="Fold T0 Override (JD)").grid(row=6, column=2, sticky="e", **pad)
        ttk.Entry(controls, textvariable=self.vars["fold_t0"], width=18).grid(row=6, column=3, sticky="w", **pad)

        ttk.Checkbutton(controls, text="Use MERR as uncertainty", variable=self.vars["use_errors"]).grid(
            row=7, column=0, columnspan=2, sticky="w", padx=8, pady=5
        )
        ttk.Checkbutton(controls, text="fit_mean", variable=self.vars["fit_mean"]).grid(
            row=7, column=2, sticky="w", padx=8, pady=5
        )
        ttk.Checkbutton(controls, text="center_data", variable=self.vars["center_data"]).grid(
            row=7, column=3, sticky="w", padx=8, pady=5
        )
        ttk.Checkbutton(controls, text="Exclude upper-limit points (<MAG)", variable=self.vars["exclude_limits"]).grid(
            row=8, column=1, columnspan=3, sticky="w", padx=8, pady=(2, 8)
        )
        ttk.Label(
            controls,
            text="Tip: Periodogram is computed on one selected filter",
            foreground="#666666",
        ).grid(row=9, column=0, columnspan=4, sticky="w", padx=8, pady=(0, 6))

        controls.columnconfigure(1, weight=1)

        online = ttk.LabelFrame(outer, text="AAVSO Online Load (WebObs)", padding=10)
        online.pack(fill="x", pady=(8, 0))
        opad = {"padx": 8, "pady": 5}

        ttk.Label(online, text="Object").grid(row=0, column=0, sticky="e", **opad)
        ttk.Entry(online, textvariable=self.vars["aavso_object"], width=24).grid(row=0, column=1, sticky="w", **opad)

        ttk.Label(online, text="Filter").grid(row=0, column=2, sticky="e", **opad)
        ttk.Entry(online, textvariable=self.vars["aavso_filter"], width=22, state="readonly").grid(
            row=0, column=3, sticky="w", **opad
        )
        ttk.Button(online, text="Select...", command=self._pick_aavso_filters).grid(
            row=0, column=4, sticky="w", **opad
        )

        ttk.Label(online, text="Observer").grid(row=0, column=5, sticky="e", **opad)
        ttk.Entry(online, textvariable=self.vars["aavso_observer"], width=12).grid(row=0, column=6, sticky="w", **opad)

        ttk.Label(online, text="Start date").grid(row=1, column=0, sticky="e", **opad)
        self.aavso_start_entry = ttk.Entry(online, textvariable=self.vars["aavso_start_date"], width=14, state="readonly")
        self.aavso_start_entry.grid(row=1, column=1, sticky="w", **opad)
        self.aavso_start_button = ttk.Button(online, text="Pick", command=lambda: self._pick_date("aavso_start_date"))
        self.aavso_start_button.grid(row=1, column=2, sticky="w", **opad)

        ttk.Label(online, text="End date").grid(row=1, column=3, sticky="e", **opad)
        self.aavso_end_entry = ttk.Entry(online, textvariable=self.vars["aavso_end_date"], width=14, state="readonly")
        self.aavso_end_entry.grid(row=1, column=4, sticky="w", **opad)
        self.aavso_end_button = ttk.Button(online, text="Pick", command=lambda: self._pick_date("aavso_end_date"))
        self.aavso_end_button.grid(row=1, column=6, sticky="w", **opad)
        ttk.Checkbutton(
            online,
            text="All dates",
            variable=self.vars["aavso_all_dates"],
            command=self._toggle_aavso_date_controls,
        ).grid(row=1, column=5, sticky="w", padx=8, pady=5)

        ttk.Label(online, text="Rows to fetch").grid(row=2, column=0, sticky="e", **opad)
        ttk.Entry(online, textvariable=self.vars["aavso_limit"], width=14).grid(row=2, column=1, sticky="w", **opad)
        ttk.Label(online, text="Use Select... to choose multiple filters", foreground="#666666").grid(
            row=2, column=3, columnspan=3, sticky="w", padx=8, pady=5
        )
        ttk.Button(online, text="Load from AAVSO", command=self._load_from_aavso).grid(
            row=2, column=6, sticky="e", padx=8, pady=(8, 5)
        )
        online.columnconfigure(6, weight=1)
        self._toggle_aavso_date_controls()

        self.summary_label = ttk.Label(
            outer,
            text="Load a file, then run analysis.",
            foreground="#304860",
        )
        self.summary_label.pack(fill="x", pady=(8, 4))

        plot_frame = ttk.LabelFrame(outer, text="Lomb-Scargle Results", padding=8)
        plot_frame.pack(fill="both", expand=True)

        self.plot_tabs = ttk.Notebook(plot_frame)
        self.plot_tabs.pack(fill="both", expand=True)
        self.tab_light = ttk.Frame(self.plot_tabs)
        self.tab_period = ttk.Frame(self.plot_tabs)
        self.plot_tabs.add(self.tab_light, text="Light Curve")
        self.plot_tabs.add(self.tab_period, text="Periodogram")

        self.fig_light, self.ax_light = plt.subplots(1, 1, figsize=(10.8, 6.4), dpi=100)
        self.fig_light.patch.set_facecolor("#f8fbff")
        self.canvas_light = FigureCanvasTkAgg(self.fig_light, master=self.tab_light)
        self.canvas_light.get_tk_widget().pack(fill="both", expand=True)

        self.fig_period, self.axes_period = plt.subplots(2, 1, figsize=(10.8, 8.6), dpi=100)
        self.fig_period.patch.set_facecolor("#f8fbff")
        self.canvas_period = FigureCanvasTkAgg(self.fig_period, master=self.tab_period)
        self.canvas_period.get_tk_widget().pack(fill="both", expand=True)
        self._draw_placeholder()
        self.root.bind("<Return>", lambda _event: self._analyze())

    def _apply_figure_layout(self, fig, hspace, top=0.965, bottom=0.075):
        fig.subplots_adjust(
            left=0.08,
            right=0.985,
            top=top,
            bottom=bottom,
            hspace=hspace,
        )

    def _save_figure_with_footer(self, fig, out_path: Path, title: str = "", include_title: bool = True):
        # Temporarily overlay export title and bottom footer ribbon during PNG save.
        subplot_state = (
            fig.subplotpars.left,
            fig.subplotpars.right,
            fig.subplotpars.top,
            fig.subplotpars.bottom,
            fig.subplotpars.wspace,
            fig.subplotpars.hspace,
        )
        orig_size = fig.get_size_inches().copy()
        artists = []
        try:
            # Keep export readable regardless of current UI tab size.
            export_w = max(float(orig_size[0]), 12.8)
            export_h = max(float(orig_size[1]), 7.8)
            fig.set_size_inches(export_w, export_h, forward=False)

            footer_height = 0.125
            footer_gap = 0.08
            export_top = min(subplot_state[2], 0.92)
            fig.subplots_adjust(
                top=export_top,
                bottom=max(subplot_state[3], footer_height + footer_gap),
                left=max(subplot_state[0], 0.09),
                right=min(subplot_state[1], 0.97),
            )
            footer = Rectangle(
                (0, 0),
                1,
                footer_height,
                transform=fig.transFigure,
                facecolor="#2f3b4d",
                edgecolor="none",
                alpha=0.98,
                zorder=6,
            )
            artists.append(footer)
            fig.add_artist(footer)
            if include_title and title.strip():
                title_y = 0.5 * (1.0 + export_top)
                title_artist = fig.text(
                    0.5,
                    title_y,
                    title,
                    ha="center",
                    va="center",
                    fontsize=21,
                    fontweight="bold",
                    color="#111111",
                    zorder=8,
                )
                artists.append(title_artist)
            footer_artist = fig.text(
                0.5,
                footer_height / 2.0,
                FOOTER_TEXT,
                ha="center",
                va="center",
                fontsize=10.0,
                color="white",
                fontweight="bold",
                linespacing=1.25,
                zorder=7,
            )
            artists.append(footer_artist)
            fig.savefig(out_path, dpi=180)
        finally:
            for artist in artists:
                try:
                    artist.remove()
                except Exception:
                    pass
            fig.subplots_adjust(
                left=subplot_state[0],
                right=subplot_state[1],
                top=subplot_state[2],
                bottom=subplot_state[3],
                wspace=subplot_state[4],
                hspace=subplot_state[5],
            )
            fig.set_size_inches(orig_size[0], orig_size[1], forward=False)

    def _export_filters_label(self) -> str:
        if self.df is None or self.df.empty or "FILT" not in self.df.columns:
            return "ALL"
        filters = []
        for value in self.df["FILT"].tolist():
            token = str(value).strip().upper()
            if token and token not in filters:
                filters.append(token)
        return ",".join(filters) if filters else "ALL"

    def _export_calendar_date_label(self) -> str:
        if self.df is None or self.df.empty or "DATE" not in self.df.columns:
            return "n/a"
        jd_values = pd.to_numeric(self.df["DATE"], errors="coerce").to_numpy(dtype=float)
        jd_values = jd_values[np.isfinite(jd_values)]
        if jd_values.size == 0:
            return "n/a"
        return _jd_to_calendar_date_text(float(np.nanmin(jd_values)))

    def _draw_placeholder(self):
        self.ax_light.clear()
        self.ax_light.grid(True, alpha=0.3)
        self.ax_light.set_title("Light Curve - load data")
        self.ax_light.set_xlabel("Julian Date", fontsize=14)
        self.ax_light.set_ylabel("Magnitude", fontsize=14)
        self.ax_light.tick_params(axis="both", labelsize=12)
        self.ax_light.invert_yaxis()

        ax1, ax3 = self.axes_period
        self._reset_period_views(ax1, ax3)

        self._apply_figure_layout(self.fig_light, hspace=0.25, top=0.90, bottom=0.11)
        self._apply_figure_layout(self.fig_period, hspace=0.34, top=0.95, bottom=0.09)
        self.canvas_light.draw_idle()
        self.canvas_period.draw_idle()

    def _reset_loaded_data(self):
        self.df = None
        self.last_result = None
        self.vars["input_file"].set("")
        self.period_filter_combo["values"] = ["ALL"]
        self.vars["periodogram_filter"].set("ALL")
        self.summary_label.config(text="Load a file, then run analysis.")
        self._draw_placeholder()
        self.plot_tabs.select(self.tab_light)

    def _reset_period_views(self, ax1=None, ax3=None):
        if ax1 is None or ax3 is None:
            ax1, ax3 = self.axes_period
        for ax in self.axes_period:
            ax.clear()
            ax.grid(True, alpha=0.3)
        ax1.set_title("Lomb-Scargle Periodogram")
        ax1.set_xlabel("Period (days)")
        ax1.set_ylabel("Power")
        ax3.set_title("Phase-folded Light Curve")
        ax3.set_xlabel("Phase")
        ax3.set_ylabel("Magnitude")
        ax3.invert_yaxis()
        ax1.text(
            0.5, 0.5,
            "Press Analyze to compute periodogram.",
            transform=ax1.transAxes,
            ha="center",
            va="center",
            color="#666666",
        )

    def _refresh_loaded_light_curve(self):
        if self.df is None or self.df.empty:
            return
        self._plot_light_curve_only()

    def _plot_light_curve_only(self):
        if self.df is None or self.df.empty:
            return

        work = self.df.copy()
        loaded_filters = sorted({str(v).strip().upper() for v in work["FILT"] if str(v).strip()})

        work = work[np.isfinite(work["DATE"]) & np.isfinite(work["MAG_VALUE"])]
        object_label = "Object"
        if "NAME" in work.columns:
            names = [str(v).strip() for v in work["NAME"].dropna().unique() if str(v).strip()]
            if names:
                object_label = names[0]
        self.ax_light.clear()
        self.ax_light.grid(True, alpha=0.3)

        if len(work) == 0:
            self.ax_light.set_title("Light Curve")
            self.ax_light.set_xlabel("Julian Date", fontsize=14)
            self.ax_light.set_ylabel("Magnitude", fontsize=14)
            self.ax_light.tick_params(axis="both", labelsize=12)
            self.ax_light.text(
                0.5,
                0.5,
                "No loaded points.",
                transform=self.ax_light.transAxes,
                ha="center",
                va="center",
                color="#666666",
            )
            self.ax_light.invert_yaxis()
            self._apply_figure_layout(self.fig_light, hspace=0.25, top=0.90, bottom=0.11)
            self.canvas_light.draw_idle()
            return

        has_err = "MERR_VALUE" in work.columns
        used_filters = []
        for band, grp in work.groupby(work["FILT"].astype(str).str.strip().str.upper(), dropna=False):
            grp = grp[np.isfinite(grp["DATE"]) & np.isfinite(grp["MAG_VALUE"])]
            if len(grp) == 0:
                continue
            jd_order = np.argsort(grp["DATE"].to_numpy(dtype=float))
            t_jd = grp["DATE"].to_numpy(dtype=float)[jd_order]
            y_jd = grp["MAG_VALUE"].to_numpy(dtype=float)[jd_order]
            dy_jd = None
            if has_err:
                dy_jd = grp["MERR_VALUE"].to_numpy(dtype=float)[jd_order]

            b = band or "UNK"
            color = AAVSO_FILTER_COLORS.get(b, "#666666")
            used_filters.append((b, color))

            if dy_jd is not None and np.isfinite(dy_jd).any():
                self.ax_light.errorbar(
                    t_jd,
                    y_jd,
                    yerr=dy_jd,
                    fmt="o",
                    ms=4,
                    color=color,
                    ecolor=color,
                    elinewidth=0.8,
                    capsize=2,
                    alpha=0.9,
                    label=b,
                )
            else:
                self.ax_light.scatter(
                    t_jd,
                    y_jd,
                    s=20,
                    color=color,
                    alpha=0.85,
                    edgecolors="none",
                    label=b,
                )

        t_all = work["DATE"].to_numpy(dtype=float)
        y_all = work["MAG_VALUE"].to_numpy(dtype=float)
        shown_filters = ",".join(loaded_filters) if loaded_filters else "ALL"
        date_label = _jd_to_calendar_date_text(float(np.nanmin(t_all)))
        self.ax_light.set_title(
            f"{object_label} | {date_label} | Filters: {shown_filters}",
            fontsize=14,
            pad=8,
        )
        self.ax_light.set_xlabel("Julian Date", fontsize=14)
        self.ax_light.set_ylabel("Magnitude", fontsize=14)
        self.ax_light.tick_params(axis="both", labelsize=12)
        self.ax_light.ticklabel_format(style="plain", useOffset=False, axis="x")
        _set_optimal_xy_limits(self.ax_light, t_all, y_all, x_frac=0.03, y_frac=0.12)
        if used_filters:
            # De-duplicate legend entries while preserving order.
            handles, labels = self.ax_light.get_legend_handles_labels()
            seen = set()
            kept_handles = []
            kept_labels = []
            for h, l in zip(handles, labels):
                if l in seen:
                    continue
                seen.add(l)
                kept_handles.append(h)
                kept_labels.append(l)
            self.ax_light.legend(kept_handles, kept_labels, title="Filter", fontsize=9, loc="best")
        self.ax_light.invert_yaxis()
        self._apply_figure_layout(self.fig_light, hspace=0.25, top=0.90, bottom=0.11)
        self.canvas_light.draw_idle()

    def _pick_date(self, var_key: str):
        current = self.vars[var_key].get().strip()
        initial_date = None
        try:
            if current:
                initial_date = dt.date.fromisoformat(current)
        except ValueError:
            initial_date = None

        dlg = CalendarDialog(self.root, initial_date=initial_date)
        self.root.wait_window(dlg)
        if dlg.result is not None:
            self.vars[var_key].set(dlg.result.isoformat())

    def _pick_aavso_filters(self):
        options = ["U", "B", "V", "R", "I", "CV", "C", "VIS", "TG", "TB", "TR", "SG", "SR", "SI", "SZ", "HA", "OIII"]
        selected = _parse_filter_selection(self.vars["aavso_filter"].get())
        dlg = FilterSelectDialog(self.root, "Select AAVSO Filters", options=options, selected=selected)
        self.root.wait_window(dlg)
        if dlg.result is None:
            return
        if not dlg.result:
            self.vars["aavso_filter"].set("ALL")
            return
        self.vars["aavso_filter"].set(",".join(dlg.result))

    def _toggle_aavso_date_controls(self):
        all_dates = bool(self.vars["aavso_all_dates"].get())
        entry_state = "disabled" if all_dates else "readonly"
        button_state = "disabled" if all_dates else "normal"
        for widget in (self.aavso_start_entry, self.aavso_end_entry):
            if widget is not None:
                widget.configure(state=entry_state)
        for widget in (self.aavso_start_button, self.aavso_end_button):
            if widget is not None:
                widget.configure(state=button_state)

    def _load_from_aavso(self):
        try:
            object_name = self.vars["aavso_object"].get().strip()
            if not object_name:
                raise ValueError("AAVSO object is required.")

            all_dates = bool(self.vars["aavso_all_dates"].get())
            start_date = self.vars["aavso_start_date"].get().strip()
            end_date = self.vars["aavso_end_date"].get().strip()
            if all_dates:
                start_date = ""
                end_date = ""
            else:
                _utc_day_to_jd_bounds(start_date, end_date)

            row_limit = int(float(self.vars["aavso_limit"].get().strip()))
            if row_limit < 200:
                row_limit = 200

            selected_filter = self.vars["aavso_filter"].get().strip()
            selected_observer = self.vars["aavso_observer"].get().strip()
            rows = fetch_aavso_rows(
                object_name=object_name,
                limit=row_limit,
                observer=selected_observer,
                filt=selected_filter,
                start_date=start_date,
                end_date=end_date,
            )
            rows_in_scope = rows
            if start_date and end_date:
                jd_start, jd_end = _utc_day_to_jd_bounds(start_date, end_date)
                rows_in_scope = []
                for row in rows:
                    jd = _to_float_or_none(row.get("jd", ""))
                    if jd is not None and jd_start <= jd < jd_end:
                        rows_in_scope.append(row)
            df_online = build_df_from_aavso_rows(
                rows=rows,
                filt=selected_filter,
                observer=selected_observer,
                start_date=start_date,
                end_date=end_date,
                object_name=object_name,
            )
            if df_online.empty:
                bands = sorted({str(r.get("band", "")).strip().upper() for r in rows_in_scope if str(r.get("band", "")).strip()})
                observers = sorted({
                    _normalize_observer_code(r.get("observer", ""))
                    for r in rows_in_scope
                    if _normalize_observer_code(r.get("observer", ""))
                })
                scope_label = "all dates" if all_dates else "selected date range"
                raise ValueError(
                    "No AAVSO rows matched the selected filters/scope.\n"
                    f"Rows in scope ({scope_label}): {len(rows_in_scope)}\n"
                    f"Available bands in scope: {', '.join(bands) if bands else 'n/a'}\n"
                    f"Available observers in scope: {', '.join(observers[:20]) if observers else 'n/a'}"
                )

            # AAVSO load overrides any existing file-loaded dataset.
            self._reset_loaded_data()
            self.df = df_online
            filters = sorted({str(v).strip().upper() for v in self.df["FILT"] if str(v).strip()})
            self.period_filter_combo["values"] = filters if filters else ["ALL"]

            self.vars["periodogram_filter"].set(filters[0] if filters else "ALL")

            date_scope = "all dates" if all_dates else f"{start_date} to {end_date}"
            self.summary_label.config(
                text=(
                    f"Loaded {len(self.df)} rows from AAVSO for {object_name} ({date_scope}). "
                    f"JD span: {float(self.df['DATE'].min()):.5f} .. {float(self.df['DATE'].max()):.5f}"
                )
            )
            self._plot_light_curve_only()
            self._reset_period_views()
            self._apply_figure_layout(self.fig_period, hspace=0.34, top=0.95, bottom=0.09)
            self.canvas_period.draw_idle()
            self.plot_tabs.select(self.tab_light)
        except Exception as exc:
            messagebox.showerror("AAVSO Load Error", str(exc))

    def _browse(self):
        selected = filedialog.askopenfilename(
            title="Select AAVSO Extended file",
            filetypes=[("Text/CSV files", "*.txt *.csv"), ("All files", "*.*")],
        )
        if not selected:
            return
        try:
            loaded_df = read_aavso_extended(selected)
            # File load overrides any existing AAVSO-loaded dataset.
            self._reset_loaded_data()
            self.vars["input_file"].set(selected)
            self.df = loaded_df
            filters = sorted({str(v).strip().upper() for v in self.df["FILT"] if str(v).strip()})
            self.period_filter_combo["values"] = filters if filters else ["ALL"]
            self.vars["periodogram_filter"].set(filters[0] if filters else "ALL")
            self.summary_label.config(text=f"Loaded {len(self.df)} rows. Available filters: {', '.join(filters) if filters else 'None'}")
            self._plot_light_curve_only()
            self._reset_period_views()
            self._apply_figure_layout(self.fig_period, hspace=0.34, top=0.95, bottom=0.09)
            self.canvas_period.draw_idle()
            self.plot_tabs.select(self.tab_light)
        except Exception as exc:
            messagebox.showerror("Load Error", str(exc))

    def _analyze(self):
        try:
            if self.df is None:
                input_path = self.vars["input_file"].get().strip()
                if not input_path:
                    raise ValueError("Select an input file first.")
                self.df = read_aavso_extended(input_path)

            selected_filter = self.vars["periodogram_filter"].get().strip().upper()
            if any(sep in selected_filter for sep in [",", "+", ";", "|", " "]):
                raise ValueError("Periodogram Filter must be a single filter (e.g. B, V, R).")
            if selected_filter in {"", "ALL"}:
                available_filters = sorted(
                    {
                        str(v).strip().upper()
                        for v in self.df["FILT"]
                        if str(v).strip()
                    }
                )
                if len(available_filters) > 1:
                    raise ValueError(
                        "For periodogram analysis select a specific filter (not ALL).\n"
                        f"Available filters: {', '.join(available_filters)}"
                    )

            min_period = float(self.vars["min_period"].get().strip())
            max_period = float(self.vars["max_period"].get().strip())
            steps = int(float(self.vars["steps"].get().strip()))
            nterms = int(float(self.vars["nterms"].get().strip()))
            fold_period_raw = self.vars["fold_period"].get().strip()
            fold_t0_raw = self.vars["fold_t0"].get().strip()
            fold_period = float(fold_period_raw) if fold_period_raw else None
            fold_t0 = float(fold_t0_raw) if fold_t0_raw else None
            if min_period <= 0 or max_period <= 0:
                raise ValueError("Period bounds must be positive.")
            if min_period >= max_period:
                raise ValueError("Min period must be smaller than max period.")
            if steps < 2:
                raise ValueError("Steps must be at least 2.")
            if nterms < 1:
                raise ValueError("nterms must be >= 1.")
            if fold_period is not None and fold_period <= 0:
                raise ValueError("Fold period override must be positive.")

            res = run_lomb_scargle(
                df=self.df,
                filt=selected_filter,
                min_period=min_period,
                max_period=max_period,
                steps=steps,
                use_errors=bool(self.vars["use_errors"].get()),
                exclude_limits=bool(self.vars["exclude_limits"].get()),
                ls_method=self.vars["ls_method"].get().strip(),
                ls_normalization=self.vars["ls_normalization"].get().strip(),
                fit_mean=bool(self.vars["fit_mean"].get()),
                center_data=bool(self.vars["center_data"].get()),
                nterms=nterms,
                fap_method=self.vars["fap_method"].get().strip(),
                fold_period_override=fold_period,
                fold_t0_override=fold_t0,
            )
            self.last_result = res
            self._plot_results(res)
            self.plot_tabs.select(self.tab_period)
            self.summary_label.config(
                text=(
                    f"Periodogram filter: {selected_filter} | "
                    f"Best period: {res['best_period']:.8f} d  "
                    f"({res['best_period_min']:.2f} min), "
                    f"Fold: P={res['fold_period']:.8f} d ({res['fold_period_min']:.2f} min), "
                    f"T0(min brightness)=JD {res['fold_t0']:.6f}, Cycles~{res['cycles_covered']:.2f}, "
                    f"FAP={res['fap_text']}, points={len(res['work'])}"
                )
            )
            edge_eps = max((max_period - min_period) * 0.01, min_period * 0.02)
            if abs(res["best_period"] - min_period) <= edge_eps or abs(res["best_period"] - max_period) <= edge_eps:
                messagebox.showwarning(
                    "Period At Search Edge",
                    (
                        "Best period is very close to a search boundary.\n"
                        "Try widening Min/Max period range for a more reliable solution."
                    ),
                )
            if res.get("used_ls_method") != self.vars["ls_method"].get().strip():
                messagebox.showwarning(
                    "Method Fallback",
                    f"Selected LS method is unavailable in this environment. Switched to '{res['used_ls_method']}'.",
                )
        except Exception as exc:
            msg = str(exc)
            if "No module named 'scipy'" in msg or 'No module named "scipy"' in msg:
                msg = (
                    "SciPy is not installed, but the current option requires it.\n\n"
                    "Switch LS Method to: auto / fast / slow / cython / chi2 / fastchi2.\n"
                    "Or install SciPy in your environment."
                )
            messagebox.showerror("Analysis Error", msg)

    def _save_outputs(self):
        try:
            if self.df is None or self.df.empty:
                raise ValueError("No data loaded. Load from file or AAVSO first.")

            out_dir = filedialog.askdirectory(title="Select output folder")
            if not out_dir:
                return
            out_path = Path(out_dir)
            out_path.mkdir(parents=True, exist_ok=True)

            object_label = "object"
            if "NAME" in self.df.columns:
                names = [str(v).strip() for v in self.df["NAME"].dropna().unique() if str(v).strip()]
                if names:
                    object_label = names[0]
            safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", object_label).strip("._") or "object"
            stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
            prefix = f"{safe_name}_{stamp}"

            # Screens
            light_png = out_path / f"{prefix}_light_curve.png"
            period_png = out_path / f"{prefix}_periodogram_phase.png"
            filters_label = self._export_filters_label()
            date_label = self._export_calendar_date_label()
            base_title = f"{object_label} | {date_label} | Filters: {filters_label}"
            # For Light Curve export, render one bold figure-level title in the top header band.
            original_light_title = self.ax_light.get_title()
            try:
                self.ax_light.set_title("")
                self._save_figure_with_footer(
                    self.fig_light,
                    light_png,
                    base_title,
                    include_title=True,
                )
            finally:
                self.ax_light.set_title(original_light_title)

            if self.last_result is not None:
                period_title = f"{base_title} | Lomb-Scargle Periodogram"
            else:
                period_title = f"{base_title} | Periodogram/Phase"
            self._save_figure_with_footer(
                self.fig_period,
                period_png,
                period_title,
            )

            # Loaded data
            loaded_csv = out_path / f"{prefix}_loaded_data.csv"
            self.df.to_csv(loaded_csv, index=False)

            saved = [str(light_png), str(period_png), str(loaded_csv)]

            # Analysis outputs if available
            if self.last_result is not None:
                res = self.last_result
                period_df = pd.DataFrame(
                    {
                        "period_days": np.asarray(res["period_sorted"], dtype=float),
                        "power": np.asarray(res["power_sorted"], dtype=float),
                    }
                )
                phased_df = pd.DataFrame(
                    {
                        "phase": np.asarray(res["phase"], dtype=float),
                        "mag": np.asarray(res["y_phase"], dtype=float),
                    }
                )
                model_df = pd.DataFrame(
                    {
                        "model_phase": np.asarray(res["model_phase"], dtype=float),
                        "model_mag": np.asarray(res["model_mag"], dtype=float),
                    }
                )
                period_csv = out_path / f"{prefix}_periodogram.csv"
                phased_csv = out_path / f"{prefix}_phased_curve.csv"
                model_csv = out_path / f"{prefix}_phased_model.csv"
                summary_txt = out_path / f"{prefix}_analysis_summary.txt"
                period_df.to_csv(period_csv, index=False)
                phased_df.to_csv(phased_csv, index=False)
                model_df.to_csv(model_csv, index=False)

                summary_lines = [
                    "Lomb-Scargle Analysis Summary",
                    f"best_period_days={res['best_period']:.10f}",
                    f"best_period_minutes={res['best_period_min']:.6f}",
                    f"fold_period_days={res['fold_period']:.10f}",
                    f"fold_period_minutes={res['fold_period_min']:.6f}",
                    f"fold_t0_jd={res['fold_t0']:.8f}",
                    f"fap={res['fap_text']}",
                    f"points={len(res['work'])}",
                ]
                summary_txt.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
                saved.extend([str(period_csv), str(phased_csv), str(model_csv), str(summary_txt)])

            messagebox.showinfo("Saved", "Saved files:\n" + "\n".join(saved))
        except Exception as exc:
            messagebox.showerror("Save Error", str(exc))

    def _plot_results(self, res):
        ax1, ax3 = self.axes_period
        ax1.clear()
        ax3.clear()

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
        p_min = float(res.get("search_min_period", np.nanmin(res["period_sorted"])))
        p_max = float(res.get("search_max_period", np.nanmax(res["period_sorted"])))
        if np.isfinite(p_min) and np.isfinite(p_max) and p_min > 0 and p_max > p_min:
            ax1.set_xlim(p_min, p_max)
            # In astronomy period scans are often visualized in log-period when range is wide.
            if (p_max / p_min) >= 20.0:
                ax1.set_xscale("log")
            else:
                ax1.set_xscale("linear")
        pwr = np.asarray(res["power_sorted"], dtype=float)
        pwr = pwr[np.isfinite(pwr)]
        if pwr.size:
            p_top = float(np.nanmax(pwr))
            ax1.set_ylim(0.0, max(0.05, p_top * 1.15))
        ax1.text(
            0.01,
            0.97,
            f"Best P = {res['best_period']:.8f} d ({res['best_period_min']:.2f} min)",
            transform=ax1.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            color="#1f4d7a",
            bbox={"boxstyle": "round,pad=0.22", "fc": "#eef4fb", "ec": "#9db5cf", "alpha": 0.9},
        )

        phase2 = np.concatenate([res["phase"], res["phase"] + 1.0])
        y2 = np.concatenate([res["y_phase"], res["y_phase"]])
        model_phase2 = np.concatenate([res["model_phase"], res["model_phase"] + 1.0])
        model_mag2 = np.concatenate([res["model_mag"], res["model_mag"]])
        ax3.scatter(phase2, y2, s=16, color="#2e7d32", alpha=0.65, edgecolors="none", label="Points")
        if np.isfinite(model_mag2).any():
            ax3.plot(model_phase2, model_mag2, color="#9c2f2f", lw=1.2, label="Model")
        ax3.set_xlim(0.0, 2.0)
        ax3.set_xticks(np.arange(0.0, 2.01, 0.2))

        # Overlay phase-binned means for a standard astronomy-style phased plot.
        bins = np.linspace(0.0, 2.0, 41)
        idx = np.digitize(phase2, bins) - 1
        bin_x = []
        bin_y = []
        bin_err = []
        for i in range(len(bins) - 1):
            mask = idx == i
            if np.count_nonzero(mask) < 3:
                continue
            yy = y2[mask]
            yy = yy[np.isfinite(yy)]
            if yy.size < 3:
                continue
            x_center = 0.5 * (bins[i] + bins[i + 1])
            bin_x.append(x_center)
            bin_y.append(float(np.mean(yy)))
            bin_err.append(float(np.std(yy)))
        if bin_x:
            ax3.errorbar(
                np.asarray(bin_x, dtype=float),
                np.asarray(bin_y, dtype=float),
                yerr=np.asarray(bin_err, dtype=float),
                fmt="o",
                ms=3.5,
                color="#1f4d7a",
                ecolor="#7fa2c7",
                capsize=2,
                elinewidth=0.7,
                alpha=0.9,
                label="Binned mean",
            )

        ax3.set_title(
            f"Phased Light Curve (Fold P = {res['fold_period']:.8f} d = {res['fold_period_min']:.2f} min, "
            f"T0 min brightness = JD {res['fold_t0']:.6f})"
        )
        ax3.set_xlabel("Phase")
        ax3.set_ylabel("Magnitude")
        ax3.grid(True, alpha=0.3)
        combined_phase_mag = y2
        if np.isfinite(model_mag2).any():
            combined_phase_mag = np.concatenate([y2, model_mag2[np.isfinite(model_mag2)]])
        y_lim = _axis_limits_with_margin(combined_phase_mag, frac=0.08)
        if y_lim is not None:
            ax3.set_ylim(*y_lim)
        ax3.invert_yaxis()
        if ax3.get_legend_handles_labels()[0]:
            ax3.legend(loc="best", fontsize=8)

        self._apply_figure_layout(self.fig_period, hspace=0.34, top=0.95, bottom=0.09)
        self.canvas_period.draw_idle()

    def _on_close(self):
        try:
            self.root.unbind_all("<MouseWheel>")
            self.root.unbind_all("<Button-4>")
            self.root.unbind_all("<Button-5>")
        except Exception:
            pass
        try:
            plt.close(self.fig_light)
            plt.close(self.fig_period)
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
