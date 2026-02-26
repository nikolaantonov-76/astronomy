"""
Apply BV/BVR transformation coefficients to AAVSO Extended File Format exports and save transformed outputs.

Usage:
1. Run: `python astronomy/PyTransformApplier.py`
2. Select BV or BVR mode and provide B/V(/R) files plus the `VPhot.ini` coefficients file.
3. Click `Execute` to create transformed `.xlsx`, `.txt`, and light-curve `.png` files.

Author: Nikola Antonov
Email: nikola.antonov@iaps.institute
"""

import configparser
import math
import re
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

aij_columns = [
    "Name", "DATE", "MAG", "MERR", "FILT", "TRANS", "MTYPE",
    "CNAME", "CMAG", "KNAME", "KMAG", "AMASS", "GROUP", "CHART", "NOTES",
]
UI_SETTINGS_FILE = Path("PyTransformApplier.settings.ini")
FOOTER_TEXT = "Nikola Antonov (nikola.antonov@iaps.institute), https://astro.iaps.instiute"

# ── Shared visual theme ────────────────────────────────────────────────────
_BG     = "#f0f4f8"
_BORDER = "#b8c8d8"
_ACCENT = "#1f4d7a"
_TEXT   = "#2c3e50"


def _apply_theme(root) -> None:
    """Apply consistent soft-blue theme to the given root window."""
    import tkinter.font as tkfont

    # 1. Activate clam theme FIRST (resets everything, so must come before customising)
    style = ttk.Style(root)
    try:
        style.theme_use("clam")
    except Exception:
        pass

    # 2. Reconfigure every named system font explicitly
    _F  = ("Segoe UI", 10)
    _FM = ("Segoe UI", 10, "bold")
    _FC = ("Consolas", 10)
    for fname in ("TkDefaultFont", "TkTextFont", "TkMenuFont",
                  "TkHeadingFont", "TkCaptionFont", "TkTooltipFont"):
        try:
            tkfont.nametofont(fname).configure(family="Segoe UI", size=10)
        except Exception:
            pass
    try:
        tkfont.nametofont("TkFixedFont").configure(family="Consolas", size=10)
    except Exception:
        pass

    # 3. option_add with "interactive" priority overrides all lower-priority theme defaults
    root.option_add("*Font",        _F,  "interactive")
    root.option_add("*Text.Font",   _FC, "interactive")
    root.option_add("*Label.Font",  _F,  "interactive")
    root.option_add("*Button.Font", _F,  "interactive")
    root.option_add("*Entry.Font",  _F,  "interactive")
    root.option_add("*Menu.Font",   _F,  "interactive")

    # 4. Configure every ttk widget type explicitly with font
    #    (style root "." alone does not cascade reliably in clam)
    style.configure(".",                 font=_F,  background=_BG, foreground=_TEXT)
    style.configure("TFrame",            background=_BG)
    style.configure("TLabel",            font=_F,  background=_BG, foreground=_TEXT)
    style.configure("TCheckbutton",      font=_F,  background=_BG, foreground=_TEXT)
    style.configure("TRadiobutton",      font=_F,  background=_BG, foreground=_TEXT)
    style.configure("TLabelframe",       background=_BG, bordercolor=_BORDER,
                                         relief="groove", borderwidth=1)
    style.configure("TLabelframe.Label", font=_FM, background=_BG, foreground=_ACCENT)
    style.configure("TButton",           font=_F,  background="#dce8f5", foreground="#1a3a5c",
                                         relief="flat", borderwidth=1, padding=(8, 4))
    style.map("TButton",
        background=[("active", "#b8d0ea"), ("pressed", "#90b8e0")],
        relief=[("pressed", "flat")])
    style.configure("TEntry",            font=_F,  fieldbackground="white",
                                         bordercolor=_BORDER, lightcolor=_BORDER, darkcolor=_BORDER)
    style.configure("TCombobox",         font=_F,  fieldbackground="white", bordercolor=_BORDER)
    style.configure("Treeview",          font=_F,  background="white",
                                         fieldbackground="white", rowheight=22)
    style.configure("Treeview.Heading",  font=_FM, background="#dce8f5", foreground="#1a3a5c")
    root.configure(bg=_BG)
# ──────────────────────────────────────────────────────────────────────────


def parse_notes_column(df):
    def extract_mags(notes):
        kv = dict(
            re.findall(r"(VMAGINS|CMAGINS|CREFMAG)=([-+]?[0-9]*\.?[0-9]+)", str(notes))
        )
        return {
            "VarInstMag": float(kv.get("VMAGINS", "nan")),
            "CompInstMag": float(kv.get("CMAGINS", "nan")),
            "CompCatMag":  float(kv.get("CREFMAG", "nan")),
        }
    return df.join(df["NOTES"].apply(extract_mags).apply(pd.Series))


def load_aij_file(filepath):
    df = pd.read_csv(filepath, names=aij_columns, comment="#")
    # DATE must be numeric (JD) to compute nearest-frame matches
    df["DATE"] = pd.to_numeric(df["DATE"], errors="coerce")
    # Remove potential trailing/leading spaces in FILT after CSV parsing
    df["FILT"] = df["FILT"].astype(str).str.strip()
    df = parse_notes_column(df)
    df["Err"] = pd.to_numeric(df["MERR"], errors="coerce")
    return df


def load_transform_coeffs(ini_file):
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(ini_file)
    return {
        "Tbv":   float(config["Coefficients"]["Tbv"]),
        "Tv_bv": float(config["Coefficients"]["Tv_bv"]),
        "Tb_bv": float(config["Coefficients"].get("Tb_bv", 0.0)),
        "Tvr":   float(config["Coefficients"].get("Tvr", 0.0)),
        "Tr_vr": float(config["Coefficients"].get("Tr_vr", 0.0)),
        "Tv_vr": float(config["Coefficients"].get("Tv_vr", 0.0)),
    }


def find_nearest(df_ref, target_date, max_minutes):
    """
    Return the nearest row in df_ref to target_date.
    If the time gap is above max_minutes, return (None, diff_minutes).
    """
    diffs = (df_ref["DATE"] - target_date).abs()
    idx = diffs.idxmin()
    diff_minutes = float(diffs[idx]) * 24 * 60
    if diff_minutes > max_minutes:
        return None, diff_minutes
    return df_ref.loc[idx], diff_minutes


def propagate_r_error(sigma_r, sigma_v, sigma_b, tvr):
    """
    Propagate photometric errors through the BVR transformation for the R band.

    The transformed R magnitude depends on the R, V and B instrumental
    measurements:
        R_std = V_std + (rc_cat - vc_cat) - Tvr * ((vs - rs) - (vc - rc))

    Taking partial derivatives and adding in quadrature gives:
        σ_R_trans = √( (Tvr·σ_r)² + σ_v² + σ_b² )

    where σ_b enters because V_std itself was derived from the BV
    transformation.  B and V errors are unchanged by their own
    transformations (Tv_bv is small and does not feed back into their
    error budgets at the level of precision reported here).

    Parameters
    ----------
    sigma_r : float  – original R MERR
    sigma_v : float  – MERR of the matched V frame
    sigma_b : float  – MERR of the matched B frame
    tvr     : float  – Tvr transformation coefficient

    Returns
    -------
    float – propagated MERR rounded to 3 decimal places
    """
    return round(math.sqrt((tvr * sigma_r) ** 2 + sigma_v ** 2 + sigma_b ** 2), 3)


def iterative_transform_bv(df_b, df_v, coeffs, max_minutes, warnings):
    results = []
    df_b_only = df_b[df_b["FILT"] == "B"]
    df_v_only = df_v[df_v["FILT"] == "V"]

    for _, v in df_v_only.iterrows():
        b_match, diff = find_nearest(df_b_only, v["DATE"], max_minutes)
        if b_match is None:
            warnings.append(
                f"V frame JD={v['DATE']:.6f}: no B frame within {max_minutes} min "
                f"(nearest is {diff:.1f} min) - saved as UNTRANSFORMED"
            )
            row = v.copy(); row["MAG_ORIG"] = row["MAG"]; row["TRANS"] = "NO"
            results.append(row); continue

        bs, bc, bc_cat = b_match["VarInstMag"], b_match["CompInstMag"], b_match["CompCatMag"]
        vs, vc, vc_cat = v["VarInstMag"], v["CompInstMag"], v["CompCatMag"]
        bs_std, vs_std = bs, vs
        for _ in range(6):
            bs_std = vs_std + (bc_cat - vc_cat) + coeffs["Tbv"] * ((bs - vs) - (bc - vc))
            vs_std = vs + (vc_cat - vc) + coeffs["Tv_bv"] * (
                (bs_std - vs_std) - (bc_cat - vc_cat)
            )
        row = v.copy(); row["MAG_ORIG"] = row["MAG"]
        row["MAG"] = round(vs_std, 3); row["TRANS"] = "YES"
        results.append(row)

    for _, b in df_b_only.iterrows():
        v_match, diff = find_nearest(df_v_only, b["DATE"], max_minutes)
        if v_match is None:
            warnings.append(
                f"B frame JD={b['DATE']:.6f}: no V frame within {max_minutes} min "
                f"(nearest is {diff:.1f} min) - saved as UNTRANSFORMED"
            )
            row = b.copy(); row["MAG_ORIG"] = row["MAG"]; row["TRANS"] = "NO"
            results.append(row); continue

        bs, bc, bc_cat = b["VarInstMag"], b["CompInstMag"], b["CompCatMag"]
        vs, vc, vc_cat = v_match["VarInstMag"], v_match["CompInstMag"], v_match["CompCatMag"]
        bs_std, vs_std = bs, vs
        for _ in range(6):
            bs_std = vs_std + (bc_cat - vc_cat) + coeffs["Tbv"] * ((bs - vs) - (bc - vc))
            vs_std = vs + (vc_cat - vc) + coeffs["Tv_bv"] * (
                (bs_std - vs_std) - (bc_cat - vc_cat)
            )
        row = b.copy(); row["MAG_ORIG"] = row["MAG"]
        row["MAG"] = round(bs_std, 3); row["TRANS"] = "YES"
        results.append(row)

    return pd.DataFrame(results).sort_values(by="DATE").reset_index(drop=True)


def iterative_transform_bvr(df_b, df_v, df_r, coeffs, max_minutes, warnings):
    results = []
    df_b_only = df_b[df_b["FILT"] == "B"]
    df_v_only = df_v[df_v["FILT"] == "V"]
    df_r_only = df_r[df_r["FILT"] == "R"]

    for _, v in df_v_only.iterrows():
        b_match, diff = find_nearest(df_b_only, v["DATE"], max_minutes)
        if b_match is None:
            warnings.append(
                f"V frame JD={v['DATE']:.6f}: no B frame within {max_minutes} min "
                f"(nearest is {diff:.1f} min) - saved as UNTRANSFORMED"
            )
            row = v.copy(); row["MAG_ORIG"] = row["MAG"]; row["TRANS"] = "NO"
            results.append(row); continue

        bs, bc, bc_cat = b_match["VarInstMag"], b_match["CompInstMag"], b_match["CompCatMag"]
        vs, vc, vc_cat = v["VarInstMag"], v["CompInstMag"], v["CompCatMag"]
        bs_std, vs_std = bs, vs
        for _ in range(6):
            bs_std = vs_std + (bc_cat - vc_cat) + coeffs["Tbv"] * ((bs - vs) - (bc - vc))
            vs_std = vs + (vc_cat - vc) + coeffs["Tv_bv"] * (
                (bs_std - vs_std) - (bc_cat - vc_cat)
            )
        row = v.copy(); row["MAG_ORIG"] = row["MAG"]
        row["MAG"] = round(vs_std, 3); row["TRANS"] = "YES"
        results.append(row)

    for _, b in df_b_only.iterrows():
        v_match, diff = find_nearest(df_v_only, b["DATE"], max_minutes)
        if v_match is None:
            warnings.append(
                f"B frame JD={b['DATE']:.6f}: no V frame within {max_minutes} min "
                f"(nearest is {diff:.1f} min) - saved as UNTRANSFORMED"
            )
            row = b.copy(); row["MAG_ORIG"] = row["MAG"]; row["TRANS"] = "NO"
            results.append(row); continue

        bs, bc, bc_cat = b["VarInstMag"], b["CompInstMag"], b["CompCatMag"]
        vs, vc, vc_cat = v_match["VarInstMag"], v_match["CompInstMag"], v_match["CompCatMag"]
        bs_std, vs_std = bs, vs
        for _ in range(6):
            bs_std = vs_std + (bc_cat - vc_cat) + coeffs["Tbv"] * ((bs - vs) - (bc - vc))
            vs_std = vs + (vc_cat - vc) + coeffs["Tv_bv"] * (
                (bs_std - vs_std) - (bc_cat - vc_cat)
            )
        row = b.copy(); row["MAG_ORIG"] = row["MAG"]
        row["MAG"] = round(bs_std, 3); row["TRANS"] = "YES"
        results.append(row)

    for _, r in df_r_only.iterrows():
        v_match, v_diff = find_nearest(df_v_only, r["DATE"], max_minutes)
        b_match, b_diff = find_nearest(df_b_only, r["DATE"], max_minutes)

        if v_match is None or b_match is None:
            missing = []
            if v_match is None: missing.append(f"V ({v_diff:.1f} min)")
            if b_match is None: missing.append(f"B ({b_diff:.1f} min)")
            warnings.append(
                f"R frame JD={r['DATE']:.6f}: no {' and '.join(missing)} frame within "
                f"{max_minutes} min - saved as UNTRANSFORMED"
            )
            row = r.copy(); row["MAG_ORIG"] = row["MAG"]; row["TRANS"] = "NO"
            results.append(row); continue

        rs, rc, rc_cat = r["VarInstMag"], r["CompInstMag"], r["CompCatMag"]
        vs, vc, vc_cat = v_match["VarInstMag"], v_match["CompInstMag"], v_match["CompCatMag"]
        bs, bc, bc_cat = b_match["VarInstMag"], b_match["CompInstMag"], b_match["CompCatMag"]
        bs_std, vs_std, rs_std = bs, vs, rs
        for _ in range(6):
            bs_std = vs_std + (bc_cat - vc_cat) + coeffs["Tbv"] * ((bs - vs) - (bc - vc))
            vs_std = vs + (vc_cat - vc) + coeffs["Tv_bv"] * (
                (bs_std - vs_std) - (bc_cat - vc_cat)
            )
            rs_std = vs_std - (vc_cat - rc_cat) - coeffs["Tvr"] * ((vs - rs) - (vc - rc))

        # Propagate photometric errors: σ_R = √((Tvr·σ_r)² + σ_v² + σ_b²)
        new_merr = propagate_r_error(
            sigma_r=float(r["Err"]),
            sigma_v=float(v_match["Err"]),
            sigma_b=float(b_match["Err"]),
            tvr=coeffs["Tvr"],
        )

        row = r.copy(); row["MAG_ORIG"] = row["MAG"]
        row["MAG"]  = round(rs_std, 3)
        row["MERR"] = str(new_merr)
        row["Err"]  = new_merr
        row["TRANS"] = "YES"
        results.append(row)

    return pd.DataFrame(results).sort_values(by="DATE").reset_index(drop=True)


def save_to_excel(df, output_path):
    df.to_excel(output_path, index=False)
    print(f"Saved to: {output_path}")


def obs_date_from_df(df):
    """Return observation date (Europe/Sofia) from the minimum JD in the series."""
    from datetime import timedelta
    min_jd = float(df["DATE"].min())
    unix_ts = (min_jd - 2440587.5) * 86400.0
    dt_utc = datetime.utcfromtimestamp(unix_ts)
    try:
        from zoneinfo import ZoneInfo
        return dt_utc.replace(tzinfo=ZoneInfo("UTC")).astimezone(ZoneInfo("Europe/Sofia")).strftime("%Y-%m-%d")
    except Exception:
        return (dt_utc + timedelta(hours=2)).strftime("%Y-%m-%d")


def save_transformed_txt(df, v_header_source, filters):
    object_name = str(df["Name"].dropna().unique()[0]).replace(" ", "_")
    obs_date = obs_date_from_df(df)
    filename = f"{object_name}_{obs_date}_{filters}.txt"
    with open(v_header_source, "r", encoding="utf-8") as f:
        header_lines = [line.strip() for line in f if line.startswith("#")]
    with open(filename, "w", encoding="utf-8") as f:
        for line in header_lines:
            f.write(line + "\n")
        df[aij_columns].to_csv(f, index=False, header=False)
    print(f"Transformed file saved as: {filename}")
    return filename


def plot_light_curve(df, output_filename, show_plot=False):
    object_name = df["Name"].dropna().unique()[0]
    filters = sorted(df["FILT"].unique())
    filter_colors = {"B": "blue", "V": "green", "R": "red"}

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(12, 8))
    for filt in filters:
        df_filt = df[df["FILT"] == filt]
        color = filter_colors.get(filt, "black")
        ax.errorbar(
            df_filt["DATE"], df_filt["MAG"], yerr=df_filt["Err"],
            fmt="o", ecolor=color, color=color,
            label=filt, capsize=2, markersize=4,
        )

    ax.invert_yaxis()
    ax.set_title(f"Light Curve for {object_name}", fontsize=16)
    ax.set_xlabel("Julian Date", fontsize=12)
    ax.set_ylabel("Magnitude", fontsize=12)
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    if show_plot:
        plt.show()
    else:
        plt.close(fig)
    print(f"Light curve plot saved as: {output_filename}")


def run_transform(scenario, b_file, v_file, r_file, ini_file, max_minutes):
    df_b = load_aij_file(b_file)
    df_v = load_aij_file(v_file)
    coeffs = load_transform_coeffs(ini_file)
    warnings = []

    if scenario == "BV":
        df_transformed = iterative_transform_bv(df_b, df_v, coeffs, max_minutes, warnings)
        filters = "BV"
    else:
        df_r = load_aij_file(r_file)
        df_transformed = iterative_transform_bvr(df_b, df_v, df_r, coeffs, max_minutes, warnings)
        filters = "BVR"

    if df_transformed.empty:
        raise ValueError("No transformed data was produced.")

    # All output files use the observation date derived from the first frame
    object_name = str(df_transformed["Name"].dropna().unique()[0]).replace(" ", "_")
    obs_date = obs_date_from_df(df_transformed)
    output_excel = f"{object_name}_{obs_date}_{filters}.xlsx"
    output_plot  = f"{object_name}_{obs_date}_{filters}_light_curve.png"

    save_to_excel(df_transformed, output_excel)
    output_txt = save_transformed_txt(df_transformed, v_header_source=v_file, filters=filters)
    plot_light_curve(df_transformed, output_plot, show_plot=True)
    return output_excel, output_txt, output_plot, warnings


class SettingsDialog(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Transform Applier Settings")
        self.resizable(True, True)
        _apply_theme(self)

        self.scenario_var    = tk.StringVar(value="BVR")
        self.b_file_var      = tk.StringVar(value="B.txt")
        self.v_file_var      = tk.StringVar(value="V.txt")
        self.r_file_var      = tk.StringVar(value="R.txt")
        self.ini_file_var    = tk.StringVar(value="VPhot.ini")
        self.max_minutes_var = tk.StringVar(value="15")
        self.coeffs_status_var = tk.StringVar(value="Coefficients: not loaded")
        self.coeffs_table = None

        self._build_ui()
        self._load_saved_settings()
        self._refresh_coefficients_table()
        self._toggle_r_controls()
        self.update_idletasks()
        w = min(max(self.winfo_reqwidth(),  700), 860)
        h = min(max(self.winfo_reqheight(), 520), 660)
        self.geometry(f"{w}x{h}")
        self.minsize(w, h)

    def _build_ui(self):
        frame = ttk.LabelFrame(self, text="Transformation Settings", padding=12)
        frame.grid(row=0, column=0, sticky="nsew", padx=12, pady=12)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        ttk.Label(frame, text="Scenario").grid(row=0, column=0, sticky="w", pady=4)
        scenario_combo = ttk.Combobox(
            frame, textvariable=self.scenario_var,
            values=["BV", "BVR"], state="readonly", width=12,
        )
        scenario_combo.grid(row=0, column=1, sticky="w", pady=4)
        scenario_combo.bind("<<ComboboxSelected>>", lambda _: self._toggle_r_controls())

        self._add_file_row(frame, 1, "B file", self.b_file_var)
        self._add_file_row(frame, 2, "V file", self.v_file_var)
        self.r_entry, self.r_button = self._add_file_row(frame, 3, "R file", self.r_file_var)
        self._add_file_row(frame, 4, "INI file", self.ini_file_var)
        ttk.Button(frame, text="Save INI", command=self._save_ini_setting).grid(
            row=4, column=3, sticky="w", pady=4
        )

        ttk.Label(frame, text="Max match (min)").grid(row=5, column=0, sticky="w", pady=4)
        ttk.Entry(frame, textvariable=self.max_minutes_var, width=12).grid(
            row=5, column=1, sticky="w", pady=4)
        ttk.Label(frame, text="Frames outside threshold -> TRANS=NO", foreground="#777777").grid(
            row=5, column=2, sticky="w", padx=(8, 0))

        ttk.Label(frame, text="Loaded coefficients").grid(row=6, column=0, sticky="nw", pady=(10, 4))
        table_frame = ttk.Frame(frame)
        table_frame.grid(row=6, column=1, columnspan=3, sticky="we", pady=(10, 4))
        columns = ("coefficient", "value", "error", "r2")
        self.coeffs_table = ttk.Treeview(table_frame, columns=columns, show="headings", height=5)
        self.coeffs_table.heading("coefficient", text="Coefficient")
        self.coeffs_table.heading("value", text="Value")
        self.coeffs_table.heading("error", text="Error")
        self.coeffs_table.heading("r2", text="R²")
        self.coeffs_table.column("coefficient", width=120, anchor="w")
        self.coeffs_table.column("value", width=85, anchor="e")
        self.coeffs_table.column("error", width=85, anchor="e")
        self.coeffs_table.column("r2", width=85, anchor="e")
        self.coeffs_table.grid(row=0, column=0, sticky="nsew")
        table_scroll = ttk.Scrollbar(table_frame, orient="vertical", command=self.coeffs_table.yview)
        self.coeffs_table.configure(yscrollcommand=table_scroll.set)
        table_scroll.grid(row=0, column=1, sticky="ns")
        ttk.Label(frame, textvariable=self.coeffs_status_var, foreground="#888888").grid(
            row=7, column=1, columnspan=3, sticky="w"
        )

        ttk.Button(frame, text="Execute", command=self._execute).grid(
            row=8, column=0, columnspan=4, sticky="ew", pady=(10, 0))
        ttk.Label(
            frame,
            text=FOOTER_TEXT,
            anchor="center",
            justify="center",
            foreground="#888888",
        ).grid(row=9, column=0, columnspan=4, sticky="ew", pady=(8, 0))

    def _add_file_row(self, parent, row, label, variable):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=4)
        entry = ttk.Entry(parent, textvariable=variable, width=36)
        entry.grid(row=row, column=1, sticky="we", pady=4, padx=(8, 8))
        button = ttk.Button(parent, text="Browse", command=lambda: self._browse_file(variable))
        button.grid(row=row, column=2, sticky="w", pady=4)
        return entry, button

    def _browse_file(self, target_var):
        selected = filedialog.askopenfilename(
            title="Select file", initialdir=str(Path.cwd()),
            filetypes=[("Text and INI files", "*.txt *.ini"), ("All files", "*.*")],
        )
        if selected:
            target_var.set(selected)
            if target_var is self.ini_file_var:
                self._refresh_coefficients_table()

    def _refresh_coefficients_table(self):
        for row_id in self.coeffs_table.get_children():
            self.coeffs_table.delete(row_id)

        ini_path = self.ini_file_var.get().strip()
        if not ini_path:
            self.coeffs_status_var.set("Coefficients: no INI file selected")
            return
        if not Path(ini_path).is_file():
            self.coeffs_status_var.set("Coefficients: INI file not found")
            return

        try:
            config = configparser.ConfigParser()
            config.optionxform = str
            config.read(ini_path, encoding="utf-8")
            coeffs = dict(config.items("Coefficients"))
            errors = dict(config.items("Error")) if config.has_section("Error") else {}
            r2_vals = (
                dict(config.items("R Squared Values"))
                if config.has_section("R Squared Values")
                else {}
            )

            for key, coeff_raw in coeffs.items():
                coeff_val = float(coeff_raw)
                err_raw = errors.get(key, "")
                r2_raw = r2_vals.get(key, "")
                err_val = f"{float(err_raw):.3f}" if err_raw.strip() else "n/a"
                r2_val = f"{float(r2_raw):.3f}" if r2_raw.strip() else "n/a"
                self.coeffs_table.insert(
                    "",
                    "end",
                    values=(key, f"{coeff_val:.3f}", err_val, r2_val),
                )

            self.coeffs_status_var.set(f"Coefficients loaded from: {ini_path}")
        except Exception as exc:
            self.coeffs_status_var.set(f"Coefficients: failed to parse INI ({exc})")

    def _load_saved_settings(self):
        if not UI_SETTINGS_FILE.is_file():
            return
        config = configparser.ConfigParser()
        try:
            config.read(UI_SETTINGS_FILE, encoding="utf-8")
            saved_ini = config.get("Settings", "ini_file", fallback="").strip()
            if saved_ini:
                self.ini_file_var.set(saved_ini)
                self._refresh_coefficients_table()
        except Exception:
            # Ignore settings read errors and keep defaults.
            return

    def _save_ini_setting(self):
        ini_path = self.ini_file_var.get().strip()
        if not ini_path:
            messagebox.showerror("Invalid input", "INI file path is empty.")
            return

        config = configparser.ConfigParser()
        config["Settings"] = {"ini_file": ini_path}
        try:
            with open(UI_SETTINGS_FILE, "w", encoding="utf-8") as f:
                config.write(f)
            self._refresh_coefficients_table()
            messagebox.showinfo("Saved", f"INI setting saved to:\n{UI_SETTINGS_FILE}")
        except Exception as exc:
            messagebox.showerror("Save failed", str(exc))

    def _toggle_r_controls(self):
        state = "normal" if self.scenario_var.get() == "BVR" else "disabled"
        self.r_entry.configure(state=state)
        self.r_button.configure(state=state)

    def _execute(self):
        scenario = self.scenario_var.get()
        b_file   = self.b_file_var.get().strip()
        v_file   = self.v_file_var.get().strip()
        r_file   = self.r_file_var.get().strip()
        ini_file = self.ini_file_var.get().strip()

        try:
            max_minutes = float(self.max_minutes_var.get().strip())
            if max_minutes <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror("Invalid input", "Max match must be a positive number (in minutes).")
            return

        required_paths = {"B file": b_file, "V file": v_file, "INI file": ini_file}
        if scenario == "BVR":
            required_paths["R file"] = r_file

        missing = [name for name, path in required_paths.items()
                   if not path or not Path(path).is_file()]
        if missing:
            messagebox.showerror("Invalid input",
                "Missing or invalid files:\n" + "\n".join(f"- {n}" for n in missing))
            return

        try:
            output_excel, output_txt, output_plot, warns = run_transform(
                scenario=scenario, b_file=b_file, v_file=v_file,
                r_file=r_file, ini_file=ini_file, max_minutes=max_minutes,
            )
        except Exception as exc:
            messagebox.showerror("Execution failed", str(exc))
            return

        msg = f"Scenario: {scenario}\nSaved:\n- {output_excel}\n- {output_txt}\n- {output_plot}"
        if warns:
            msg += f"\n\nWarning: {len(warns)} frame(s) saved as UNTRANSFORMED:\n"
            msg += "\n".join(f"• {w}" for w in warns)
        messagebox.showinfo("Done", msg)


if __name__ == "__main__":
    SettingsDialog().mainloop()