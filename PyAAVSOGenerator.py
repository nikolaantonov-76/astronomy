"""
Generate AAVSO Extended Format reports from AstroImageJ (AIJ) photometry data.
"""

import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
import math
import re
import json
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.image as mpimg

# --- CONFIGURATION ---
STAR_NAME = "KR Aur"
FILTER = "V"
OBSERVER_CODE = "ANNA"
INPUT_TBL_FILE = "KR Aur 2026-02-18 V.tbl"
CHART_ID = "X40996ANY"

# --- Star configuration ---
SINGLE_COMP_TBL_ID = "na"
SINGLE_COMP_AAVSO_NAME = "ENSEMBLE"
COMP_STAR_STD_MAG_ERR = 0.05

# --- Check-star configuration ---
CHECK_STAR_TBL_ID = "T4"
CHECK_STAR_AAVSO_NAME = "133"
CHECK_STAR_STD_MAG = 13.308
CHECK_STAR_STD_MAG_ERR = 0.029

TRANSFORMED = "NO"

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
    _F  = ("Segoe UI", 8)
    _FM = ("Segoe UI", 8, "bold")
    _FC = ("Consolas", 8)
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

# --- END CONFIGURATION ---
FOOTER_TEXT = "Nikola Antonov (nikola.antonov@iaps.institute), https://astro.iaps.instiute"


def _to_float_or_none(value):
    try:
        val = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(val) or math.isinf(val):
        return None
    return val


def _fmt_or_na(value, decimals=3):
    num = _to_float_or_none(value)
    if num is None:
        return "na"
    return f"{num:.{decimals}f}"


def _series_local_date_str_from_min_jd(df, tz_name="Europe/Sofia"):
    try:
        min_jd = float(df["JD_UTC"].min())
        unix_ts = (min_jd - 2440587.5) * 86400.0
        dt_utc = datetime.utcfromtimestamp(unix_ts)
        try:
            from zoneinfo import ZoneInfo
            dt_local = dt_utc.replace(tzinfo=ZoneInfo("UTC")).astimezone(ZoneInfo(tz_name))
        except Exception:
            dt_local = dt_utc + timedelta(hours=2)
        return dt_local.strftime("%Y-%m-%d")
    except Exception:
        return datetime.now().strftime("%Y-%m-%d")


def get_instrumental_mag(row, obj_id):
    """
    Compute instrumental magnitude normalized by exposure time,
    following the VPhot convention:
        m_inst = -2.5 * log10(flux / EXPTIME)
    If EXPTIME is missing or 0, use an unnormalized ZP=0 fallback.
    """
    flux_col = f"Source-Sky_{obj_id}"
    if flux_col not in row or pd.isna(row[flux_col]):
        return None

    flux_val = _to_float_or_none(row[flux_col])
    if flux_val is None or flux_val <= 0:
        return None

    exptime = _to_float_or_none(row.get("EXPTIME", None))
    if exptime and exptime > 0:
        return -2.5 * math.log10(flux_val / exptime)
    else:
        # Fallback: unnormalized (ZP=0)
        return -2.5 * math.log10(flux_val)


def build_notes_string(row, single_comp_id, check_id, check_std_mag, comp_err, check_err):
    parts = []
    use_ensemble = (not single_comp_id or single_comp_id.lower() == "na")

    # 1. CHECK STAR
    kmag_measured = _to_float_or_none(row.get(f"Source_AMag_{check_id}", None))
    kmag_ins = get_instrumental_mag(row, check_id)
    if kmag_measured is not None:
        parts.append(f"|KMAGSTD={kmag_measured:.3f}")
    if kmag_ins is not None:
        parts.append(f"|KMAGINS={kmag_ins:.3f}")
    if check_std_mag is not None:
        parts.append(f"|KREFMAG={check_std_mag:.3f}")
    if check_err is not None:
        parts.append(f"|KREFERR={check_err:.3f}")

    # 2. VARIABLE STAR (T1)
    vmag_ins = get_instrumental_mag(row, "T1")
    if vmag_ins is not None:
        parts.append(f"|VMAGINS={vmag_ins:.3f}")

    # 3. COMPARISON STAR(S)
    if use_ensemble:
        comp_ids = sorted(
            set([re.search(r'C(\d+)', col).group(0)
                 for col in row.index if re.search(r'Source-Sky_C(\d+)', col)]),
            key=lambda x: int(x[1:])
        )

        # CMAGINS = -2.5 * log10(total ensemble flux / EXPTIME)
        # Physically correct: the ensemble is treated as a single star
        # with combined flux, matching AIJ's internal calculation.
        exptime = _to_float_or_none(row.get("EXPTIME", None))
        total_flux = 0.0
        for cid in comp_ids:
            flux_col = f"Source-Sky_{cid}"
            fv = _to_float_or_none(row.get(flux_col, None))
            if fv is not None and fv > 0:
                total_flux += fv

        if total_flux > 0:
            if exptime and exptime > 0:
                ensemble_cmagins = -2.5 * math.log10(total_flux / exptime)
            else:
                ensemble_cmagins = -2.5 * math.log10(total_flux)
            parts.append(f"|CMAGINS={ensemble_cmagins:.3f}")

        cref_vals = []
        for cid in comp_ids:
            cref = _to_float_or_none(row.get(f"Source_AMag_{cid}", None))
            if cref is not None:
                cref_vals.append(cref)
        if cref_vals:
            avg_crefmag = sum(cref_vals) / len(cref_vals)
            parts.append(f"|CREFMAG={avg_crefmag:.3f}")

        parts.append("|ENSTYPE=1")
    else:
        cmag_ins = get_instrumental_mag(row, single_comp_id)
        cref_mag = _to_float_or_none(row.get(f"Source_AMag_{single_comp_id}", None))
        if cmag_ins is not None:
            parts.append(f"|CMAGINS={cmag_ins:.3f}")
        if cref_mag is not None:
            parts.append(f"|CREFMAG={cref_mag:.3f}")

    if comp_err is not None:
        parts.append(f"|CREFERR={comp_err:.3f}")

    return "".join(parts) if parts else "na"


def plot_light_curve(df, star_name, filt, check_star_id, check_star_name):
    if df.empty or len(df) < 2:
        print("Not enough data to generate the plot.")
        return None

    filter_colors = {'U': 'darkviolet', 'B': 'blue', 'V': 'green', 'R': 'red'}
    data_color = filter_colors.get(filt.upper(), 'black')

    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(14, 8))

    ax.errorbar(
        df['JD_UTC'], df['Source_AMag_T1'], yerr=df['Source_AMag_Err_T1'],
        fmt='o', color=data_color, ecolor='dimgray',
        capsize=3, markersize=6, label=f'{star_name}'
    )

    kmag_col = f'Source_AMag_{check_star_id}'
    kerr_col = f'Source_AMag_Err_{check_star_id}'
    if kmag_col in df.columns and kerr_col in df.columns:
        ax.errorbar(
            df['JD_UTC'], df[kmag_col], yerr=df[kerr_col],
            fmt='^', color='gray', ecolor='silver',
            capsize=3, markersize=6, label=f'Check: {check_star_name}'
        )

    min_jd = float(df['JD_UTC'].min())
    max_jd = float(df['JD_UTC'].max())
    jd_span = max_jd - min_jd
    x_pad = max(jd_span * 0.08, 5.0 / (24.0 * 60.0))
    x_left = min_jd - x_pad
    x_right = max_jd + x_pad
    ax.set_xlim(x_left, x_right)
    ax.set_xticks(np.linspace(start=x_left, stop=x_right, num=7))

    ax.invert_yaxis()
    ax.set_title(f'Light Curve for {star_name} - Filter {filt}', fontsize=16)
    ax.set_ylabel('Magnitude', fontsize=12)
    ax.legend(fontsize=12)

    jd_int_part = int(df['JD_UTC'].iloc[0])
    ax.set_xlabel(f'JD (UTC) - {jd_int_part}', fontsize=12)
    ax.ticklabel_format(style='plain', useOffset=False, axis='x')
    tick_locations = ax.get_xticks()
    ax.set_xticklabels([f"{loc - jd_int_part:.4f}" for loc in tick_locations], rotation=45, ha='right')

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    instrument_credit = "N. Antonov, AO Meshtitsa | 0.25-m f/4.8 Newtonian | ASI533MM Pro"
    plt.figtext(0.5, 0.01, instrument_credit, ha='center', va='bottom', fontsize=10, style='italic', color='gray')

    today_str = _series_local_date_str_from_min_jd(df)
    plot_filename = f"{star_name.replace(' ', '_')}_{today_str}_{filt}_light_curve.png"

    plt.savefig(plot_filename, dpi=300)
    plt.close(fig)
    print(f"Plot created successfully: '{plot_filename}'")
    return plot_filename


def create_aavso_report(input_file, star_name, filt, obs_code,
                        single_comp_id, single_comp_name,
                        check_star_id, check_star_name, check_star_std_mag,
                        comp_err, check_err, transformed, chart_id):
    check_star_id = (check_star_id or "").strip()
    try:
        df = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        print(f"Error reading file: {e}")
        return None, None, []

    if not check_star_id:
        raise ValueError("Check Star TBL ID is required.")

    required_cols = [
        "JD_UTC", "Source_AMag_T1", "Source_AMag_Err_T1", "AIRMASS",
        f"Source_AMag_{check_star_id}",
    ]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in .tbl file: {', '.join(missing_cols)}")

    plot_filename = plot_light_curve(df, star_name, filt, check_star_id, check_star_name)

    aavso_data = []
    for _, row in df.iterrows():
        use_ensemble = (not single_comp_id or single_comp_id.lower() == "na")

        cname = "ENSEMBLE" if use_ensemble else single_comp_name

        # CMAG: instrumental magnitude of comp/ensemble (consistent with VPhot)
        if use_ensemble:
            exptime = _to_float_or_none(row.get("EXPTIME", None))
            comp_ids = sorted(
                set([re.search(r'C(\d+)', col).group(0)
                     for col in row.index if re.search(r'Source-Sky_C(\d+)', col)]),
                key=lambda x: int(x[1:])
            )
            total_flux = sum(
                _to_float_or_none(row.get(f"Source-Sky_{cid}", None)) or 0
                for cid in comp_ids
            )
            if total_flux > 0:
                if exptime and exptime > 0:
                    cmag_ins = -2.5 * math.log10(total_flux / exptime)
                else:
                    cmag_ins = -2.5 * math.log10(total_flux)
                cmag_str = _fmt_or_na(cmag_ins)
            else:
                cmag_str = "na"
        else:
            cmag_str = _fmt_or_na(get_instrumental_mag(row, single_comp_id))

        # KMAG: instrumental magnitude of the check star (consistent with VPhot)
        kmag_str = _fmt_or_na(get_instrumental_mag(row, check_star_id))

        notes = build_notes_string(row, single_comp_id, check_star_id, check_star_std_mag, comp_err, check_err)

        aavso_row = [
            star_name, _fmt_or_na(row.get("JD_UTC"), 8), _fmt_or_na(row.get("Source_AMag_T1"), 3),
            _fmt_or_na(row.get("Source_AMag_Err_T1"), 3), filt, transformed, "STD",
            cname, cmag_str, check_star_name, kmag_str,
            _fmt_or_na(row.get("AIRMASS"), 4), "0", chart_id, notes
        ]
        aavso_data.append(",".join(map(str, aavso_row)))

    today_str = _series_local_date_str_from_min_jd(df)
    output_filename = f"{star_name.replace(' ', '_')}_{today_str}_{filt}.txt"

    header = [
        "#TYPE=EXTENDED",
        f"#OBSCODE={obs_code}",
        "#SOFTWARE=AstroImageJ to AAVSO Report Generator v1.1",
        "#DELIM=,",
        "#DATE=JD",
        "#OBSTYPE=CCD",
        "#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES"
    ]

    with open(output_filename, 'w') as f:
        f.write("\n".join(header))
        f.write("\n")
        f.write("\n".join(aavso_data))

    print(f"AAVSO report created successfully: '{output_filename}'")
    preview_lines = header + aavso_data
    return output_filename, plot_filename, preview_lines


class AavsoGeneratorGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("AAVSO Report Generator")
        _apply_theme(root)
        self.config_dir = Path.cwd() / "aavso_configs"
        self.config_dir.mkdir(parents=True, exist_ok=True)

        self.vars = {
            "star_name": tk.StringVar(value=STAR_NAME),
            "filt": tk.StringVar(value=FILTER),
            "obs_code": tk.StringVar(value=OBSERVER_CODE),
            "input_file": tk.StringVar(value=INPUT_TBL_FILE),
            "chart_id": tk.StringVar(value=CHART_ID),
            "single_comp_id": tk.StringVar(value=SINGLE_COMP_TBL_ID),
            "single_comp_name": tk.StringVar(value=SINGLE_COMP_AAVSO_NAME),
            "comp_err": tk.StringVar(value=str(COMP_STAR_STD_MAG_ERR)),
            "check_id": tk.StringVar(value=CHECK_STAR_TBL_ID),
            "check_name": tk.StringVar(value=CHECK_STAR_AAVSO_NAME),
            "check_std_mag": tk.StringVar(value=str(CHECK_STAR_STD_MAG)),
            "check_err": tk.StringVar(value=str(CHECK_STAR_STD_MAG_ERR)),
            "transformed": tk.StringVar(value=TRANSFORMED),
            "profile_name": tk.StringVar(value=STAR_NAME.replace(" ", "_")),
        }
        self._build_ui()
        self._refresh_profiles()
        self._fit_window_to_content()

    def _fit_window_to_content(self):
        self.root.update_idletasks()
        width  = min(max(self.root.winfo_reqwidth(),  820), 980)
        height = min(max(self.root.winfo_reqheight(), 620), 740)
        self.root.geometry(f"{width}x{height}")
        self.root.minsize(width, height)

    def _build_ui(self):
        pad = {"padx": 8, "pady": 5}
        ttk.Label(
            self.root,
            text=FOOTER_TEXT,
            anchor="center",
            justify="center",
            foreground="#888888",
        ).pack(side="bottom", fill="x", padx=14, pady=(0, 6))
        action_bar = ttk.Frame(self.root)
        action_bar.pack(side="bottom", fill="x", padx=14, pady=(0, 6))
        ttk.Button(action_bar, text="Execute", command=self._generate).pack(side="right")

        outer = ttk.Frame(self.root)
        outer.pack(fill="both", expand=True, padx=14, pady=(14, 8))

        frame = ttk.LabelFrame(outer, text="Observation Settings", padding=10)
        frame.pack(fill="both", expand=True)

        fields = [
            ("Star Name", "star_name"),
            ("Filter", "filt"),
            ("Observer Code", "obs_code"),
            ("Input .tbl File", "input_file"),
            ("Chart ID", "chart_id"),
            ("Single Comp TBL ID", "single_comp_id"),
            ("Single Comp AAVSO Name", "single_comp_name"),
            ("Comp Std Mag Err", "comp_err"),
            ("Check Star TBL ID", "check_id"),
            ("Check Star AAVSO Name", "check_name"),
            ("Check Star Std Mag", "check_std_mag"),
            ("Check Star Std Mag Err", "check_err"),
            ("Transformed", "transformed"),
        ]

        profile_row = ttk.Frame(frame)
        profile_row.grid(row=0, column=0, columnspan=3, sticky="we", padx=8, pady=(0, 8))
        profile_row.columnconfigure(1, weight=1)
        ttk.Label(profile_row, text="Configuration").grid(row=0, column=0, sticky="w", padx=(0, 8))
        self.profile_combo = ttk.Combobox(profile_row, textvariable=self.vars["profile_name"])
        self.profile_combo.grid(row=0, column=1, sticky="we", padx=(0, 8))
        ttk.Button(profile_row, text="Load Settings", command=self._load_profile).grid(row=0, column=2, padx=(0, 6))
        ttk.Button(profile_row, text="Save Settings", command=self._save_profile).grid(row=0, column=3)

        for row_idx, (label, key) in enumerate(fields, start=1):
            ttk.Label(frame, text=label).grid(row=row_idx, column=0, sticky="e", **pad)
            if key == "filt":
                widget = ttk.Combobox(frame, textvariable=self.vars[key], values=["U", "B", "V", "R", "I"], width=28, state="readonly")
            elif key == "transformed":
                widget = ttk.Combobox(frame, textvariable=self.vars[key], values=["NO", "YES"], width=28, state="readonly")
            else:
                widget = ttk.Entry(frame, textvariable=self.vars[key], width=32)
            widget.grid(row=row_idx, column=1, sticky="we", **pad)

            if key == "input_file":
                ttk.Button(frame, text="Browse", command=self._browse_file).grid(row=row_idx, column=2, sticky="w", **pad)

        frame.columnconfigure(1, weight=1)

        hint = (
            "Use 'na' in Single Comp TBL ID to enable ensemble mode.\n"
            "The report and light-curve image are saved in the current folder."
        )
        ttk.Label(frame, text=hint, foreground="#555555").grid(
            row=len(fields) + 1, column=0, columnspan=3, sticky="w", padx=8, pady=(10, 12))

        self.root.bind("<Return>", lambda _event: self._generate())

    def _browse_file(self):
        selected = filedialog.askopenfilename(
            title="Select .tbl file",
            filetypes=[("TBL files", "*.tbl"), ("All files", "*.*")]
        )
        if selected:
            self.vars["input_file"].set(selected)

    def _read_float(self, key, label):
        val = _to_float_or_none(self.vars[key].get().strip())
        if val is None:
            raise ValueError(f"Invalid numeric value for '{label}'.")
        return val

    def _safe_profile_name(self, raw_name):
        name = raw_name.strip()
        if not name:
            name = self.vars["star_name"].get().strip() or "default"
        name = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._")
        return name or "default"

    def _profile_path(self, profile_name):
        return self.config_dir / f"{profile_name}.json"

    def _refresh_profiles(self):
        profiles = sorted(path.stem for path in self.config_dir.glob("*.json"))
        self.profile_combo["values"] = profiles

    def _current_settings(self):
        keys = [
            "star_name", "filt", "obs_code", "input_file", "chart_id",
            "single_comp_id", "single_comp_name", "comp_err", "check_id",
            "check_name", "check_std_mag", "check_err", "transformed"
        ]
        return {key: self.vars[key].get().strip() for key in keys}

    def _apply_settings(self, settings):
        for key, value in settings.items():
            if key in self.vars and key != "profile_name":
                self.vars[key].set(str(value))

    def _save_profile(self):
        try:
            profile_name = self._safe_profile_name(self.vars["profile_name"].get())
            self.vars["profile_name"].set(profile_name)
            profile_path = self._profile_path(profile_name)
            with open(profile_path, "w", encoding="utf-8") as f:
                json.dump(self._current_settings(), f, indent=2)
            self._refresh_profiles()
            messagebox.showinfo("Saved", f"Configuration saved:\n{profile_path}")
        except Exception as exc:
            messagebox.showerror("Error", f"Could not save configuration:\n{exc}")

    def _load_profile(self):
        try:
            profile_name = self._safe_profile_name(self.vars["profile_name"].get())
            profile_path = self._profile_path(profile_name)
            if not profile_path.exists():
                raise ValueError(f"Configuration not found:\n{profile_path}")
            with open(profile_path, "r", encoding="utf-8") as f:
                settings = json.load(f)
            if not isinstance(settings, dict):
                raise ValueError("Invalid configuration format.")
            self._apply_settings(settings)
            self.vars["profile_name"].set(profile_name)
            messagebox.showinfo("Loaded", f"Configuration loaded:\n{profile_path}")
        except Exception as exc:
            messagebox.showerror("Error", f"Could not load configuration:\n{exc}")

    def _generate(self):
        try:
            input_file = self.vars["input_file"].get().strip()
            if not input_file:
                raise ValueError("Input .tbl file is required.")
            if not Path(input_file).exists():
                raise ValueError(f"Input file not found: {input_file}")

            output_file, plot_file, generated_values = create_aavso_report(
                input_file=input_file,
                star_name=self.vars["star_name"].get().strip(),
                filt=self.vars["filt"].get().strip().upper(),
                obs_code=self.vars["obs_code"].get().strip(),
                single_comp_id=self.vars["single_comp_id"].get().strip(),
                single_comp_name=self.vars["single_comp_name"].get().strip(),
                check_star_id=self.vars["check_id"].get().strip(),
                check_star_name=self.vars["check_name"].get().strip(),
                check_star_std_mag=self._read_float("check_std_mag", "Check Star Std Mag"),
                comp_err=self._read_float("comp_err", "Comp Std Mag Err"),
                check_err=self._read_float("check_err", "Check Star Std Mag Err"),
                transformed=self.vars["transformed"].get().strip().upper(),
                chart_id=self.vars["chart_id"].get().strip(),
            )

            if not output_file:
                raise RuntimeError("Report generation failed.")

            msg = f"Report created:\n{output_file}"
            if plot_file:
                msg += f"\n\nPlot created:\n{plot_file}"
            messagebox.showinfo("Success", msg)
            if plot_file and Path(plot_file).exists():
                self._show_plot_window(plot_file)
            self._show_values_window(generated_values)
        except Exception as exc:
            messagebox.showerror("Error", str(exc))

    def _show_plot_window(self, plot_file):
        preview = tk.Toplevel(self.root)
        preview.title(f"Light Curve Preview - {Path(plot_file).name}")
        preview.geometry("980x640")

        fig, ax = plt.subplots(figsize=(10, 6), dpi=100)
        img = mpimg.imread(plot_file)
        ax.imshow(img)
        ax.axis("off")
        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=preview)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

        def _on_close():
            plt.close(fig)
            preview.destroy()

        preview.protocol("WM_DELETE_WINDOW", _on_close)

    def _show_values_window(self, values):
        viewer = tk.Toplevel(self.root)
        viewer.title("Generated AAVSO Values")
        viewer.geometry("1180x520")

        ttk.Label(
            viewer,
            text="Generated lines (same values written in the output file):"
        ).pack(anchor="w", padx=10, pady=(10, 4))

        text = ScrolledText(viewer, wrap="none")
        text.pack(fill="both", expand=True, padx=10, pady=(0, 10))
        text.insert("1.0", "\n".join(values) if values else "No values generated.")
        text.configure(state="disabled")


if __name__ == "__main__":
    root = tk.Tk()
    app = AavsoGeneratorGUI(root)
    root.mainloop()