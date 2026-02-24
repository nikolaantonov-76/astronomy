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
    # DATE трябва да е число (JD) за аритметика при намиране на най-близък кадър
    df["DATE"] = pd.to_numeric(df["DATE"], errors="coerce")
    # Премахване на евентуални интервали от FILT след CSV парсване
    df["FILT"] = df["FILT"].astype(str).str.strip()
    df = parse_notes_column(df)
    df["Err"] = pd.to_numeric(df["MERR"], errors="coerce")
    return df


def load_transform_coeffs(ini_file):
    config = configparser.ConfigParser()
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
    Връща най-близкия ред от df_ref до target_date.
    Ако разликата надвишава max_minutes, връща (None, diff_minutes).
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
                f"V кадър JD={v['DATE']:.6f}: няма B кадър в рамките на {max_minutes} мин "
                f"(най-близкият е {diff:.1f} мин) — записан НЕТРАНСФОРМИРАН"
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
                f"B кадър JD={b['DATE']:.6f}: няма V кадър в рамките на {max_minutes} мин "
                f"(най-близкият е {diff:.1f} мин) — записан НЕТРАНСФОРМИРАН"
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
                f"V кадър JD={v['DATE']:.6f}: няма B кадър в рамките на {max_minutes} мин "
                f"(най-близкият е {diff:.1f} мин) — записан НЕТРАНСФОРМИРАН"
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
                f"B кадър JD={b['DATE']:.6f}: няма V кадър в рамките на {max_minutes} мин "
                f"(най-близкият е {diff:.1f} мин) — записан НЕТРАНСФОРМИРАН"
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
            if v_match is None: missing.append(f"V ({v_diff:.1f} мин)")
            if b_match is None: missing.append(f"B ({b_diff:.1f} мин)")
            warnings.append(
                f"R кадър JD={r['DATE']:.6f}: няма {' и '.join(missing)} кадър в рамките на "
                f"{max_minutes} мин — записан НЕТРАНСФОРМИРАН"
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
    """Връща датата на наблюдението (Europe/Sofia) от минималния JD в серията."""
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

    # Всички изходни файлове използват датата на наблюдението от първия кадър
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
        self.resizable(False, False)

        self.scenario_var    = tk.StringVar(value="BVR")
        self.b_file_var      = tk.StringVar(value="B.txt")
        self.v_file_var      = tk.StringVar(value="V.txt")
        self.r_file_var      = tk.StringVar(value="R.txt")
        self.ini_file_var    = tk.StringVar(value="VPhot.ini")
        self.max_minutes_var = tk.StringVar(value="15")

        self._build_ui()
        self._toggle_r_controls()

    def _build_ui(self):
        frame = ttk.Frame(self, padding=12)
        frame.grid(row=0, column=0, sticky="nsew")

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

        ttk.Label(frame, text="Max match (min)").grid(row=5, column=0, sticky="w", pady=4)
        ttk.Entry(frame, textvariable=self.max_minutes_var, width=12).grid(
            row=5, column=1, sticky="w", pady=4)
        ttk.Label(frame, text="Кадри извън прага → TRANS=NO", foreground="#777777").grid(
            row=5, column=2, sticky="w", padx=(8, 0))

        ttk.Button(frame, text="Execute", command=self._execute).grid(
            row=6, column=0, columnspan=3, sticky="ew", pady=(10, 0))

    def _add_file_row(self, parent, row, label, variable):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=4)
        entry = ttk.Entry(parent, textvariable=variable, width=45)
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
            messagebox.showerror("Invalid input", "Max match трябва да е положително число (в минути).")
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
            msg += f"\n\n⚠️  {len(warns)} кадър(а) записани НЕТРАНСФОРМИРАНИ:\n"
            msg += "\n".join(f"• {w}" for w in warns)
        messagebox.showinfo("Done", msg)


if __name__ == "__main__":
    SettingsDialog().mainloop()