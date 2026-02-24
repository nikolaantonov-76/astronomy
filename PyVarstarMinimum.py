"""Variable-star minima calculator with VSX lookup."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
import json
import math
import tkinter as tk
from tkinter import messagebox, ttk
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET


VSX_API_URL = "https://www.aavso.org/vsx/index.php"
DATETIME_FMT = "%Y-%m-%d %H:%M"


# --- JD <-> datetime conversion ---
def datetime_to_jd(dt: datetime) -> float:
    if dt.tzinfo is None or dt.utcoffset() != timedelta(0):
        raise ValueError("Datetime must be timezone-aware UTC.")

    y, m = dt.year, dt.month
    if m <= 2:
        y -= 1
        m += 12

    a = math.floor(y / 100)
    b = 2 - a + math.floor(a / 4)
    day = dt.day + (
        dt.hour
        + dt.minute / 60
        + dt.second / 3600
        + dt.microsecond / 3.6e9
    ) / 24

    return (
        math.floor(365.25 * (y + 4716))
        + math.floor(30.6001 * (m + 1))
        + day
        + b
        - 1524.5
    )


def jd_to_datetime(jd: float) -> datetime:
    jd += 0.5
    z = int(jd)
    f = jd - z

    if z < 2299161:
        a = z
    else:
        alpha = int((z - 1867216.25) / 36524.25)
        a = z + 1 + alpha - int(alpha / 4)

    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)

    day = b - d - int(30.6001 * e) + f
    month = e - 1 if e < 14 else e - 13
    year = c - 4716 if month > 2 else c - 4715

    day_int = int(day)
    frac_day = day - day_int
    sec_total = frac_day * 86400

    hour = int(sec_total // 3600)
    minute = int((sec_total % 3600) // 60)
    sec = int(sec_total % 60)
    usec = int((sec_total - int(sec_total)) * 1e6)

    return datetime(year, month, day_int, hour, minute, sec, usec, tzinfo=timezone.utc)


def generate_minima(jd_epoch: float, period: float, start_utc: datetime, end_utc: datetime):
    jd_start = datetime_to_jd(start_utc)
    jd_end = datetime_to_jd(end_utc)

    n = math.ceil((jd_start - jd_epoch) / period)
    while True:
        jd_min = jd_epoch + n * period
        if jd_min > jd_end:
            break
        yield jd_to_datetime(jd_min)
        n += 1


# --- VSX integration ---
def _pick_first(mapping: dict[str, str], *candidates: str) -> str:
    for key in candidates:
        value = mapping.get(key.lower(), "").strip()
        if value:
            return value
    return ""


def _flatten_json_strings(value, out: dict[str, str], parent_key: str = ""):
    if isinstance(value, dict):
        for key, val in value.items():
            _flatten_json_strings(val, out, key)
        return
    if isinstance(value, list):
        for item in value:
            _flatten_json_strings(item, out, parent_key)
        return
    if value is None:
        return

    text = str(value).strip()
    if text and parent_key:
        out[parent_key.lower()] = text


def _parse_vsx_payload(raw_text: str) -> dict[str, str]:
    raw_text = raw_text.strip()
    if not raw_text:
        raise ValueError("Empty response from VSX.")

    # JSON branch
    if raw_text.startswith("{") or raw_text.startswith("["):
        parsed = json.loads(raw_text)
        flat = {}
        _flatten_json_strings(parsed, flat)
    else:
        # XML branch
        root = ET.fromstring(raw_text)
        flat = {}
        for node in root.iter():
            tag = node.tag
            if "}" in tag:
                tag = tag.split("}", 1)[1]
            text = (node.text or "").strip()
            if text:
                flat[tag.lower()] = text

    # Common error indicators in API responses
    err = _pick_first(flat, "error", "errormsg", "message", "status")
    if err and "ok" not in err.lower() and "success" not in err.lower():
        if "not found" in err.lower() or "no" in err.lower():
            raise ValueError(f"VSX: {err}")

    result = {
        "name": _pick_first(flat, "name", "obj", "object", "ident"),
        "type": _pick_first(flat, "vartype", "type", "variabilitytype"),
        "period": _pick_first(flat, "period", "periodvalue"),
        "epoch": _pick_first(flat, "epoch", "epoch_hjd", "hjd0", "minima", "m0"),
        "ra": _pick_first(flat, "ra", "ra2000", "raj2000", "raj"),
        "dec": _pick_first(flat, "dec", "declination2000", "de2000", "dej2000", "declination"),
        "max_mag": _pick_first(flat, "max", "maxmag", "maxmagnitude"),
        "min_mag": _pick_first(flat, "min", "minmag", "minmagnitude"),
    }

    if not any(result.values()):
        raise ValueError("VSX did not return recognized object fields.")

    return result


def fetch_vsx_info(star_name: str, timeout: int = 15) -> dict[str, str]:
    star_name = star_name.strip()
    if not star_name:
        raise ValueError("Enter a star name first.")

    params = urlencode({"view": "api.object", "ident": star_name})
    url = f"{VSX_API_URL}?{params}"
    req = Request(url, headers={"User-Agent": "PyVarstarMinimum/1.0"})

    try:
        with urlopen(req, timeout=timeout) as response:
            body = response.read().decode("utf-8", errors="replace")
    except HTTPError as exc:
        raise RuntimeError(f"VSX HTTP error: {exc.code}") from exc
    except URLError as exc:
        raise RuntimeError(f"VSX connection error: {exc.reason}") from exc

    return _parse_vsx_payload(body)


class VarstarMinimumApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("Variable Star Minima")
        self.root.geometry("920x620")
        self.root.minsize(820, 560)

        now_utc = datetime.now(timezone.utc)
        default_start = (now_utc - timedelta(days=2)).strftime(DATETIME_FMT)
        default_end = (now_utc + timedelta(days=2)).strftime(DATETIME_FMT)

        self.vars = {
            "star_name": tk.StringVar(value=""),
            "jd_epoch": tk.StringVar(value="2444293.0235"),
            "period": tk.StringVar(value="0.158432492"),
            "start_utc": tk.StringVar(value=default_start),
            "end_utc": tk.StringVar(value=default_end),
            "tz_offset": tk.StringVar(value="3"),
            "vsx_summary": tk.StringVar(value="VSX: not queried"),
        }

        self._build_menu()
        self._build_ui()

    def _build_menu(self):
        menu = tk.Menu(self.root)
        file_menu = tk.Menu(menu, tearoff=0)
        file_menu.add_command(label="Exit", command=self.root.quit)
        menu.add_cascade(label="File", menu=file_menu)

        tools_menu = tk.Menu(menu, tearoff=0)
        tools_menu.add_command(label="Fetch From AAVSO VSX", command=self._fetch_vsx)
        menu.add_cascade(label="Tools", menu=tools_menu)

        self.root.config(menu=menu)

    def _build_ui(self):
        frame = ttk.Frame(self.root, padding=12)
        frame.pack(fill="both", expand=True)

        fields = [
            ("Variable Star Name", "star_name"),
            ("JD Epoch", "jd_epoch"),
            ("Period (days)", "period"),
            (f"Start UTC ({DATETIME_FMT})", "start_utc"),
            (f"End UTC ({DATETIME_FMT})", "end_utc"),
            ("Local UTC Offset (hours)", "tz_offset"),
        ]

        for idx, (label, key) in enumerate(fields):
            ttk.Label(frame, text=label).grid(row=idx, column=0, sticky="e", padx=8, pady=6)
            ttk.Entry(frame, textvariable=self.vars[key], width=36).grid(
                row=idx, column=1, sticky="we", padx=8, pady=6
            )

            if key == "star_name":
                ttk.Button(frame, text="VSX Lookup", command=self._fetch_vsx).grid(
                    row=idx, column=2, sticky="w", padx=8, pady=6
                )

        ttk.Label(frame, textvariable=self.vars["vsx_summary"], foreground="#1f4d7a").grid(
            row=len(fields), column=0, columnspan=3, sticky="w", padx=8, pady=(4, 10)
        )

        ttk.Button(frame, text="Calculate Minima", command=self._calculate).grid(
            row=len(fields) + 1, column=0, columnspan=3, sticky="we", padx=8, pady=(0, 10)
        )

        self.output = tk.Text(frame, wrap="none", height=18)
        self.output.grid(row=len(fields) + 2, column=0, columnspan=3, sticky="nsew", padx=8, pady=(0, 8))

        frame.columnconfigure(1, weight=1)
        frame.rowconfigure(len(fields) + 2, weight=1)

        self.root.bind("<Return>", lambda _event: self._calculate())

    def _parse_float(self, key: str, label: str) -> float:
        raw = self.vars[key].get().strip()
        try:
            return float(raw)
        except ValueError as exc:
            raise ValueError(f"Invalid {label}: {raw}") from exc

    def _parse_utc_datetime(self, key: str, label: str) -> datetime:
        raw = self.vars[key].get().strip()
        try:
            dt = datetime.strptime(raw, DATETIME_FMT)
        except ValueError as exc:
            raise ValueError(f"Invalid {label}. Use format: {DATETIME_FMT}") from exc
        return dt.replace(tzinfo=timezone.utc)

    def _fetch_vsx(self):
        try:
            info = fetch_vsx_info(self.vars["star_name"].get())

            if info["name"]:
                self.vars["star_name"].set(info["name"])

            if info["period"]:
                self.vars["period"].set(info["period"])

            if info["epoch"]:
                self.vars["jd_epoch"].set(info["epoch"])

            summary = (
                f"VSX: {info['name'] or 'n/a'} | Type: {info['type'] or 'n/a'} | "
                f"Period: {info['period'] or 'n/a'} | Epoch: {info['epoch'] or 'n/a'} | "
                f"RA: {info['ra'] or 'n/a'} | Dec: {info['dec'] or 'n/a'}"
            )
            self.vars["vsx_summary"].set(summary)
        except Exception as exc:
            messagebox.showerror("VSX Lookup", str(exc))

    def _calculate(self):
        try:
            jd_epoch = self._parse_float("jd_epoch", "JD epoch")
            period = self._parse_float("period", "period")
            if period <= 0:
                raise ValueError("Period must be positive.")

            start_utc = self._parse_utc_datetime("start_utc", "start UTC")
            end_utc = self._parse_utc_datetime("end_utc", "end UTC")
            if end_utc < start_utc:
                raise ValueError("End UTC must be after start UTC.")

            tz_hours = self._parse_float("tz_offset", "local UTC offset")
            local_tz = timezone(timedelta(hours=tz_hours))

            minima = list(generate_minima(jd_epoch, period, start_utc, end_utc))

            lines = []
            star = self.vars["star_name"].get().strip() or "(unnamed)"
            lines.append(f"Star: {star}")
            lines.append(f"JD Epoch: {jd_epoch}")
            lines.append(f"Period: {period} d")
            lines.append(f"Window: {start_utc.strftime('%Y-%m-%d %H:%M:%S')} UTC -> {end_utc.strftime('%Y-%m-%d %H:%M:%S')} UTC")
            lines.append("")
            lines.append(f"{'N':>3} | {'JD':>14} | {'UTC':<19} | {'Local':<19}")
            lines.append("-" * 72)

            for idx, utc_dt in enumerate(minima):
                jd_val = datetime_to_jd(utc_dt)
                local_dt = utc_dt.astimezone(local_tz)
                lines.append(
                    f"{idx:>3} | {jd_val:>14.6f} | {utc_dt.strftime('%Y-%m-%d %H:%M:%S')} | {local_dt.strftime('%Y-%m-%d %H:%M:%S')}"
                )

            if not minima:
                lines.append("No minima in the selected interval.")

            self.output.delete("1.0", tk.END)
            self.output.insert("1.0", "\n".join(lines))
        except Exception as exc:
            messagebox.showerror("Calculation Error", str(exc))


def main():
    root = tk.Tk()
    VarstarMinimumApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
