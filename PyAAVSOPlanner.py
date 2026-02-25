"""
Time/Altitude planner with settings menu and graphical visualization.

Usage:
1. Run: `python astronomy/PyAAVSOPlanner.py`
2. Enter target name, site coordinates, date, altitude threshold, and sampling step.
3. Click `Calculate` to plot altitude and view rise/set and recent AAVSO observations.

Author: Nikola Antonov
Email: nikola.antonov@iaps.institute
"""

import datetime as dt
import html
import json
from pathlib import Path
import queue
import re
import threading
import zoneinfo

import astropy.units as u
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import requests
import tkinter as tk
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time, TimeDelta
from astropy.utils import iers
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import messagebox, ttk
from tkinter.scrolledtext import ScrolledText


DEFAULTS = {
    "object_name": "AM Her",
    "latitude": "42.67",
    "longitude": "23.0",
    "obs_date": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d"),
    "target_altitude": "30.0",
    "sampling_min": "1",
    "timezone": "Europe/Sofia",
    "show_rise_set": True,
}

_SESAME = [
    "https://cds.unistra.fr/cgi-bin/nph-sesame/-oJ?{name}",
    "https://cds.u-strasbg.fr/cgi-bin/nph-sesame/-oJ?{name}",
    "https://cds.unistra.fr/cgi-bin/nph-sesame/-oI?{name}",
    "https://cds.u-strasbg.fr/cgi-bin/nph-sesame/-oI?{name}",
]

_AAVSO_WEBOBS_URL = "https://apps.aavso.org/webobs/results/"

SETTINGS_FILE = Path("PyTimeAltitude.settings.json")
# Prevent long blocking attempts to download IERS data in offline/limited networks.
iers.conf.auto_download = False
iers.conf.auto_max_age = None
iers.conf.remote_timeout = 3.0


def _json_ra_dec(txt):
    data = json.loads(txt)
    return float(data["target"]["resolver"]["jradeg"]), float(data["target"]["resolver"]["jdedeg"])


def _ascii_ra_dec(txt):
    match = re.search(r"^%J\s+([+-]?\d+(?:\.\d*)?)\s+([+-]?\d+(?:\.\d*)?)", txt, re.M)
    if not match:
        raise ValueError("No %J line found.")
    return float(match[1]), float(match[2])


def get_coord(name: str) -> SkyCoord:
    for url in _SESAME:
        try:
            response = requests.get(url.format(name=requests.utils.quote(name)), timeout=4)
            response.raise_for_status()
            ra, dec = (_json_ra_dec if "-oJ" in url else _ascii_ra_dec)(response.text)
            return SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
        except Exception:
            continue
    raise RuntimeError(f"Sesame could not resolve coordinates for '{name}'.")


def altitude_deg(moment: Time, coord: SkyCoord, loc: EarthLocation) -> float:
    return coord.transform_to(AltAz(obstime=moment, location=loc)).alt.deg


def find_crossings(coord, loc, date_iso, h_deg, step_min):
    # Search in the same time window as the plot: -12h..+12h around midnight
    midnight = Time(f"{date_iso} 00:00:00", scale="utc")
    times = midnight + np.arange(-720, 720 + step_min, step_min) * u.min
    altitudes = np.array([altitude_deg(t, coord, loc) for t in times])
    diff = altitudes - h_deg
    change = np.where(np.diff(np.sign(diff)) != 0)[0]

    events = []
    for idx in change:
        lo, hi = times[idx], times[idx + 1]
        for _ in range(22):
            mid = lo + (hi - lo) / 2
            if (hi - lo) < TimeDelta(1 * u.s):
                break
            if (altitude_deg(mid, coord, loc) - h_deg) * (altitude_deg(lo, coord, loc) - h_deg) <= 0:
                hi = mid
            else:
                lo = mid
        event_time = lo + (hi - lo) / 2
        direction = "up" if altitudes[idx] < h_deg else "down"
        events.append((event_time, direction))
    return sorted(events, key=lambda x: x[0].jd)


def build_day_curve_centered_midnight(coord, loc, date_iso, step_min):
    midnight = Time(f"{date_iso} 00:00:00", scale="utc")
    times = midnight + np.arange(-720, 720 + step_min, step_min) * u.min
    altitudes = np.array([altitude_deg(t, coord, loc) for t in times])
    return times, altitudes


def build_sun_alt_curve(times: Time, loc: EarthLocation):
    sun_alt = get_sun(times).transform_to(AltAz(obstime=times, location=loc)).alt.deg
    return np.array(sun_alt, dtype=float)


def pretty_time(when: Time, tz_name=""):
    utc_dt = when.to_datetime(timezone=dt.timezone.utc)
    if not tz_name.strip():
        return f"{utc_dt.strftime('%Y-%m-%d %H:%M:%S')} UTC"
    try:
        tz = zoneinfo.ZoneInfo(tz_name)
        local_dt = utc_dt.astimezone(tz)
        return (
            f"{utc_dt.strftime('%Y-%m-%d %H:%M:%S')} UTC"
            f"  |  {local_dt.strftime('%Y-%m-%d %H:%M:%S')} {local_dt.tzname()}"
        )
    except Exception:
        return f"{utc_dt.strftime('%Y-%m-%d %H:%M:%S')} UTC"


def _strip_html_tags(raw: str) -> str:
    text = re.sub(r"<[^>]+>", " ", raw, flags=re.S)
    text = html.unescape(text)
    return re.sub(r"\s+", " ", text).strip()


def _is_jd_like(value: str) -> bool:
    try:
        jd = float(value.strip())
    except Exception:
        return False
    return 2400000.0 <= jd <= 2600000.0


def _to_float_or_none(value: str):
    try:
        return float(str(value).strip())
    except Exception:
        return None


def _has_magnitude_like_value(value: str) -> bool:
    cleaned = value.strip()
    if not cleaned:
        return False
    if cleaned in {"—", "-", "na", "n/a"}:
        return False
    return bool(re.search(r"\d", cleaned))


def _is_mag_like(value: str) -> bool:
    cleaned = str(value).strip().replace("<", "").replace(">", "")
    mag = _to_float_or_none(cleaned)
    if mag is None:
        return False
    return -5.0 <= mag <= 30.0


def _find_jd_in_cells(values):
    for cell in values:
        if _is_jd_like(cell):
            return cell
    return ""


def _find_mag_in_cells(values):
    for cell in values:
        if _is_mag_like(cell):
            return cell
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


def _find_date_in_cells(values):
    for cell in values:
        token = str(cell).strip()
        if re.search(r"\d{4}[-/ ]\d{2}[-/ ]\d{2}", token):
            return token
        if re.search(r"\d{4}\s+[A-Za-z]{3,}\.?\s+\d{2}\.\d+", token):
            return token
    return ""


def _extract_rows_from_html_loose(page_html: str):
    found = []
    rows = re.findall(r"<tr\b.*?>.*?</tr>", page_html, flags=re.I | re.S)
    for row in rows:
        cells = re.findall(r"<t[dh]\b.*?>(.*?)</t[dh]>", row, flags=re.I | re.S)
        if len(cells) < 3:
            continue
        values = [_strip_html_tags(cell) for cell in cells]
        jd_value = _find_jd_in_cells(values)
        mag_value = _find_mag_in_cells(values)
        if not _is_jd_like(jd_value):
            continue
        if not _is_mag_like(mag_value):
            continue
        found.append(
            {
                "date": _find_date_in_cells(values),
                "jd": jd_value,
                "mag": mag_value,
                "band": _find_band_in_cells(values),
                "observer": _find_observer_in_cells(values),
            }
        )
    return found


def _extract_rows_from_next_data(page_html: str):
    script_match = re.search(
        r'<script[^>]*id=["\']__NEXT_DATA__["\'][^>]*>(.*?)</script>',
        page_html,
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
            mag_val = (
                lowered.get("magnitude")
                or lowered.get("mag")
                or lowered.get("value")
                or ""
            )
            if _is_jd_like(str(jd_val)) and _has_magnitude_like_value(str(mag_val)):
                found.append(
                    {
                        "date": str(lowered.get("date", "") or lowered.get("obs_date", "") or ""),
                        "jd": str(jd_val),
                        "mag": str(mag_val),
                        "band": str(lowered.get("band", "") or lowered.get("filter", "") or ""),
                        "observer": str(lowered.get("observer", "") or lowered.get("obscode", "") or ""),
                    }
                )
            for val in node.values():
                walk(val)
            return
        if isinstance(node, list):
            for item in node:
                walk(item)

    walk(payload)
    return found


def _parse_webobs_rows(page_html: str):
    all_rows = _extract_rows_from_next_data(page_html)
    tables = re.findall(r"<table\b.*?>.*?</table>", page_html, flags=re.I | re.S)
    for table in tables:
        rows = re.findall(r"<tr\b.*?>.*?</tr>", table, flags=re.I | re.S)
        if len(rows) < 2:
            continue

        header_cells = re.findall(r"<th\b.*?>(.*?)</th>", rows[0], flags=re.I | re.S)
        headers = [_strip_html_tags(cell).lower() for cell in header_cells]
        if not headers:
            continue

        idx_jd = next((i for i, h in enumerate(headers) if h == "jd"), None)
        idx_date = next((i for i, h in enumerate(headers) if "date" in h), None)
        idx_mag = next((i for i, h in enumerate(headers) if "magnitude" in h or h == "mag"), None)
        idx_band = next((i for i, h in enumerate(headers) if "band" in h or "filter" in h), None)
        idx_obs = next((i for i, h in enumerate(headers) if "observer" in h or h == "obs"), None)
        if idx_jd is None or idx_mag is None:
            continue

        for row in rows[1:]:
            cells = re.findall(r"<td\b.*?>(.*?)</td>", row, flags=re.I | re.S)
            if not cells:
                continue
            values = [_strip_html_tags(cell) for cell in cells]
            if idx_mag is not None and idx_mag >= len(values):
                continue
            jd_value = values[idx_jd] if idx_jd is not None and idx_jd < len(values) else _find_jd_in_cells(values)
            mag_value = values[idx_mag] if idx_mag is not None else _find_mag_in_cells(values)
            if not _is_jd_like(jd_value):
                continue
            if not _has_magnitude_like_value(mag_value):
                continue
            all_rows.append(
                {
                    "date": values[idx_date] if idx_date is not None and idx_date < len(values) else "",
                    "jd": jd_value,
                    "mag": mag_value,
                    "band": values[idx_band] if idx_band is not None and idx_band < len(values) else "",
                    "observer": values[idx_obs] if idx_obs is not None and idx_obs < len(values) else "",
                }
            )

    all_rows.extend(_extract_rows_from_html_loose(page_html))

    if not all_rows:
        return []

    unique_by_jd = {}
    for row in all_rows:
        key = (
            row.get("jd", ""),
            str(row.get("mag", "")).strip(),
            str(row.get("observer", "")).strip(),
            str(row.get("band", "")).strip(),
        )
        unique_by_jd[key] = row

    return sorted(unique_by_jd.values(), key=lambda item: float(item["jd"]), reverse=True)


def fetch_latest_aavso_observations(object_name: str, limit: int = 5):
    n = str(max(200, limit))
    headers = {
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
        "User-Agent": (
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/120.0.0.0 Safari/537.36"
        ),
    }
    # Try WebObs using both query parameter variants (star= and target=)
    candidates = [
        (_AAVSO_WEBOBS_URL, {"star": object_name, "num_results": n}),
        (_AAVSO_WEBOBS_URL, {"target": object_name, "num_results": n}),
    ]
    last_error = "No response from AAVSO."
    for url, params in candidates:
        try:
            response = requests.get(url, params=params, timeout=(8, 30), headers=headers)
            response.raise_for_status()
            rows = _parse_webobs_rows(response.text)
            if rows:
                unique_by_jd = {}
                for row in rows:
                    key = (
                        row.get("jd", ""),
                        str(row.get("mag", "")).strip(),
                        str(row.get("observer", "")).strip(),
                        str(row.get("band", "")).strip(),
                    )
                    unique_by_jd[key] = row
                ordered = sorted(unique_by_jd.values(), key=lambda item: float(item["jd"]), reverse=True)
                return ordered[:limit], None
            last_error = f"No observations parsed from {url}"
        except requests.exceptions.Timeout:
            last_error = f"Timeout: {url}"
            continue
        except Exception as exc:
            last_error = f"{url}: {exc}"
            continue
    return [], last_error


class TimeAltitudeApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("AAVSO Observation Planner")
        self.geometry("1180x760")
        self.minsize(980, 640)

        self.vars = {
            "object_name": tk.StringVar(value=DEFAULTS["object_name"]),
            "latitude": tk.StringVar(value=DEFAULTS["latitude"]),
            "longitude": tk.StringVar(value=DEFAULTS["longitude"]),
            "obs_date": tk.StringVar(value=DEFAULTS["obs_date"]),
            "target_altitude": tk.StringVar(value=DEFAULTS["target_altitude"]),
            "sampling_min": tk.StringVar(value=DEFAULTS["sampling_min"]),
            "timezone": tk.StringVar(value=DEFAULTS["timezone"]),
            "show_rise_set": tk.BooleanVar(value=DEFAULTS["show_rise_set"]),
            "coord_info": tk.StringVar(value="Coordinates: not resolved"),
        }

        self.figure, self.ax = plt.subplots(figsize=(8.2, 5.2), dpi=100)
        self.canvas = None
        self.run_button = None
        self._query_running = False
        self._query_token = 0
        self._result_queue = queue.Queue()
        self._poll_after_id = None

        self._build_menu()
        self._build_ui()
        self._load_settings_from_disk(show_feedback=False)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _build_menu(self):
        main_menu = tk.Menu(self)

        file_menu = tk.Menu(main_menu, tearoff=0)
        file_menu.add_command(label="Run", command=self._run)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self._on_close)
        main_menu.add_cascade(label="File", menu=file_menu)

        settings_menu = tk.Menu(main_menu, tearoff=0)
        settings_menu.add_command(label="Save Settings", command=self._save_settings_to_disk)
        settings_menu.add_command(label="Load Settings", command=self._load_settings_from_disk)
        settings_menu.add_command(label="Reset Defaults", command=self._reset_defaults)
        settings_menu.add_command(label="Set Date = Today (UTC)", command=self._set_today_utc)
        settings_menu.add_command(label="Set Date = Tomorrow (UTC)", command=self._set_tomorrow_utc)
        main_menu.add_cascade(label="Settings", menu=settings_menu)

        help_menu = tk.Menu(main_menu, tearoff=0)
        help_menu.add_command(label="About", command=self._show_about)
        main_menu.add_cascade(label="Help", menu=help_menu)

        self.config(menu=main_menu)

    def _build_ui(self):
        root_frame = ttk.Frame(self, padding=12)
        root_frame.pack(fill="both", expand=True)

        left = ttk.LabelFrame(root_frame, text="Query Settings", padding=10)
        left.pack(side="left", fill="y", padx=(0, 10))

        right = ttk.Frame(root_frame)
        right.pack(side="right", fill="both", expand=True)

        fields = [
            ("Object Name", "object_name"),
            ("Latitude (deg)", "latitude"),
            ("Longitude (deg)", "longitude"),
            ("Date UTC (YYYY-MM-DD)", "obs_date"),
            ("Target Altitude (deg)", "target_altitude"),
            ("Sampling (min)", "sampling_min"),
            ("Timezone", "timezone"),
        ]

        for row, (label, key) in enumerate(fields):
            ttk.Label(left, text=label).grid(row=row, column=0, sticky="w", pady=4)
            ttk.Entry(left, textvariable=self.vars[key], width=28).grid(row=row, column=1, sticky="we", pady=4, padx=(8, 0))

        ttk.Checkbutton(
            left,
            text="Include rise/set (0 deg)",
            variable=self.vars["show_rise_set"],
        ).grid(row=len(fields), column=0, columnspan=2, sticky="w", pady=(8, 6))

        self.run_button = ttk.Button(left, text="Run Query", command=self._run)
        self.run_button.grid(row=len(fields) + 1, column=0, columnspan=2, sticky="we", pady=(6, 2))

        ttk.Label(left, textvariable=self.vars["coord_info"], foreground="#1f4d7a", wraplength=290, justify="left").grid(
            row=len(fields) + 2, column=0, columnspan=2, sticky="w", pady=(10, 0)
        )

        left.columnconfigure(1, weight=1)
        left.rowconfigure(len(fields) + 2, weight=1)

        plot_frame = ttk.LabelFrame(right, text="Altitude Chart", padding=8)
        plot_frame.pack(fill="both", expand=True)

        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        text_frame = ttk.LabelFrame(right, text="Events", padding=8)
        text_frame.pack(fill="both", expand=False, pady=(10, 0))
        self.output = ScrolledText(text_frame, wrap="word", height=11)
        self.output.pack(fill="both", expand=True)

    def _reset_defaults(self):
        for key, value in DEFAULTS.items():
            if key == "show_rise_set":
                self.vars[key].set(bool(value))
            else:
                self.vars[key].set(str(value))
        self.vars["coord_info"].set("Coordinates: not resolved")

    def _snapshot_settings(self):
        data = {}
        for key in DEFAULTS:
            if key == "show_rise_set":
                data[key] = bool(self.vars[key].get())
            else:
                data[key] = str(self.vars[key].get())
        return data

    def _apply_settings(self, data):
        for key in DEFAULTS:
            if key not in data:
                continue
            if key == "show_rise_set":
                self.vars[key].set(bool(data[key]))
            else:
                self.vars[key].set(str(data[key]))

    def _save_settings_to_disk(self):
        try:
            SETTINGS_FILE.write_text(
                json.dumps(self._snapshot_settings(), indent=2, ensure_ascii=True) + "\n",
                encoding="utf-8",
            )
            messagebox.showinfo("Settings", f"Settings saved to {SETTINGS_FILE}")
        except Exception as exc:
            messagebox.showerror("Settings Error", f"Could not save settings: {exc}")

    def _load_settings_from_disk(self, show_feedback=True):
        if not SETTINGS_FILE.exists():
            if show_feedback:
                messagebox.showinfo("Settings", f"No saved settings found at {SETTINGS_FILE}")
            return
        try:
            data = json.loads(SETTINGS_FILE.read_text(encoding="utf-8"))
            if not isinstance(data, dict):
                raise ValueError("Settings file is not a JSON object.")
            self._apply_settings(data)
            if show_feedback:
                messagebox.showinfo("Settings", f"Settings loaded from {SETTINGS_FILE}")
        except Exception as exc:
            if show_feedback:
                messagebox.showerror("Settings Error", f"Could not load settings: {exc}")

    def _set_today_utc(self):
        self.vars["obs_date"].set(dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d"))

    def _set_tomorrow_utc(self):
        date_val = dt.datetime.now(dt.timezone.utc) + dt.timedelta(days=1)
        self.vars["obs_date"].set(date_val.strftime("%Y-%m-%d"))

    def _show_about(self):
        messagebox.showinfo(
            "About",
            "Computes altitude crossings and visualizes altitude vs time for one UTC date.\n"
            "Coordinates are resolved through CDS Sesame.",
        )

    def _read_inputs(self):
        object_name = self.vars["object_name"].get().strip()
        if not object_name:
            raise ValueError("Object name is required.")

        try:
            latitude = float(self.vars["latitude"].get().strip())
            longitude = float(self.vars["longitude"].get().strip())
            target_altitude = float(self.vars["target_altitude"].get().strip())
            sampling_min = int(self.vars["sampling_min"].get().strip())
            if sampling_min <= 0:
                raise ValueError
        except ValueError as exc:
            raise ValueError("Latitude/Longitude/Target/Sampling must be valid numbers.") from exc

        if abs(latitude) > 90:
            raise ValueError("Latitude must be in [-90, 90].")
        if abs(longitude) > 180:
            raise ValueError("Longitude must be in [-180, 180].")

        obs_date = self.vars["obs_date"].get().strip()
        dt.datetime.strptime(obs_date, "%Y-%m-%d")

        tz_name = self.vars["timezone"].get().strip()
        if tz_name:
            try:
                zoneinfo.ZoneInfo(tz_name)
            except Exception as exc:
                raise ValueError(f"Invalid timezone: {tz_name}") from exc

        return object_name, latitude, longitude, obs_date, target_altitude, sampling_min, tz_name

    def _run(self):
        if self._query_running:
            return
        try:
            (
                object_name,
                latitude,
                longitude,
                obs_date,
                target_altitude,
                sampling_min,
                tz_name,
            ) = self._read_inputs()
            self._set_query_running(True, "Resolving coordinates and computing altitude...")
            self._query_token += 1
            token = self._query_token
            threading.Thread(
                target=self._run_worker,
                args=(
                    token,
                    object_name,
                    latitude,
                    longitude,
                    obs_date,
                    target_altitude,
                    sampling_min,
                    tz_name,
                    self.vars["show_rise_set"].get(),
                ),
                daemon=True,
            ).start()
            self._start_polling_results()
        except Exception as exc:
            messagebox.showerror("Query Error", str(exc))

    def _set_query_running(self, running: bool, status_text: str = ""):
        self._query_running = running
        if self.run_button is not None:
            self.run_button.configure(state=("disabled" if running else "normal"))
        self.config(cursor=("watch" if running else ""))
        if status_text:
            self.output.delete("1.0", tk.END)
            self.output.insert("1.0", status_text + "\n")
        self.update_idletasks()

    def _start_polling_results(self):
        if self._poll_after_id is not None:
            self.after_cancel(self._poll_after_id)
            self._poll_after_id = None
        self._poll_after_id = self.after(120, self._poll_results)

    def _poll_results(self):
        self._poll_after_id = None
        handled = False
        while True:
            try:
                kind, token, payload = self._result_queue.get_nowait()
            except queue.Empty:
                break
            if token != self._query_token:
                continue
            handled = True
            if kind == "error":
                self._set_query_running(False)
                messagebox.showerror("Query Error", payload)
            else:
                coord = payload["coord"]
                self.vars["coord_info"].set(
                    f"Coordinates: RA={coord.ra.deg:.6f} deg, Dec={coord.dec.deg:.6f} deg"
                )
                self._render_text(
                    payload["object_name"],
                    payload["obs_date"],
                    payload["latitude"],
                    payload["longitude"],
                    payload["target_altitude"],
                    payload["tz_name"],
                    payload["rise_set_events"],
                    payload["target_events"],
                    payload["aavso_rows"],
                    payload["aavso_error"],
                )
                self._render_plot(
                    payload["object_name"],
                    payload["obs_date"],
                    payload["target_altitude"],
                    payload["tz_name"],
                    payload["curve_times"],
                    payload["curve_alt"],
                    payload["sun_altitudes"],
                    payload["rise_set_events"],
                    payload["target_events"],
                )
                self._set_query_running(False)
        if self._query_running and not handled:
            self._poll_after_id = self.after(120, self._poll_results)

    def _run_worker(
        self,
        token,
        object_name,
        latitude,
        longitude,
        obs_date,
        target_altitude,
        sampling_min,
        tz_name,
        include_rise_set,
    ):
        try:
            coord = get_coord(object_name)
            site = EarthLocation(lat=latitude * u.deg, lon=longitude * u.deg)
            target_events = find_crossings(coord, site, obs_date, target_altitude, sampling_min)
            rise_set_events = (
                find_crossings(coord, site, obs_date, 0.0, sampling_min) if include_rise_set else []
            )
            curve_step = max(1, min(sampling_min, 5))
            curve_times, curve_alt = build_day_curve_centered_midnight(coord, site, obs_date, curve_step)
            sun_altitudes = build_sun_alt_curve(curve_times, site)
            aavso_rows, aavso_error = fetch_latest_aavso_observations(object_name, limit=5)
            payload = {
                "object_name": object_name,
                "latitude": latitude,
                "longitude": longitude,
                "obs_date": obs_date,
                "target_altitude": target_altitude,
                "tz_name": tz_name,
                "coord": coord,
                "rise_set_events": rise_set_events,
                "target_events": target_events,
                "curve_times": curve_times,
                "curve_alt": curve_alt,
                "sun_altitudes": sun_altitudes,
                "aavso_rows": aavso_rows,
                "aavso_error": aavso_error,
            }
            self._result_queue.put(("ok", token, payload))
        except Exception as exc:
            self._result_queue.put(("error", token, str(exc)))

    def _render_text(
        self,
        object_name,
        obs_date,
        latitude,
        longitude,
        target_altitude,
        tz_name,
        rise_set_events,
        target_events,
        aavso_rows,
        aavso_error,
    ):
        lines = [
            f"Object:   {object_name}",
            f"Date:     {obs_date} (UTC)",
            f"Site:     {latitude:+.3f} deg, {longitude:+.3f} deg",
            "",
        ]

        if self.vars["show_rise_set"].get():
            lines.append("0 deg crossings (rise/set):")
            if not rise_set_events:
                lines.append("  - No 0 deg crossings on this date.")
            else:
                for when, direction in rise_set_events:
                    label = "rise" if direction == "up" else "set"
                    lines.append(f"  - {label:<4}: {pretty_time(when, tz_name)}")
            lines.append("")

        lines.append(f"{target_altitude:.1f} deg crossings:")
        if not target_events:
            lines.append("  - No target-altitude crossings on this date.")
        else:
            for when, direction in target_events:
                label = "up" if direction == "up" else "down"
                lines.append(f"  - {label:<4}: {pretty_time(when, tz_name)}")

        lines.append("")
        lines.append("Latest 5 AAVSO observations:")
        if aavso_rows:
            for row in aavso_rows:
                date_label = row["date"] or "n/a"
                jd_label = row["jd"] or "n/a"
                mag_label = row["mag"] or "n/a"
                band_label = row["band"] or "n/a"
                obs_label = row["observer"] or "n/a"
                lines.append(
                    f"  - Date: {date_label} | JD: {jd_label} | Mag: {mag_label} | Band: {band_label} | Obs: {obs_label}"
                )
        else:
            lines.append(f"  - Could not fetch observations: {aavso_error}")

        self.output.delete("1.0", tk.END)
        self.output.insert("1.0", "\n".join(lines))

    def _to_plot_datetimes(self, times, tz_name):
        utc_datetimes = [t.to_datetime(timezone=dt.timezone.utc) for t in times]
        if not tz_name:
            return utc_datetimes, "UTC"
        tz = zoneinfo.ZoneInfo(tz_name)
        return [d.astimezone(tz) for d in utc_datetimes], tz_name

    def _render_plot(
        self,
        object_name,
        obs_date,
        target_altitude,
        tz_name,
        curve_times,
        curve_alt,
        sun_altitudes,
        rise_set_events,
        target_events,
    ):
        x_datetimes, axis_tz = self._to_plot_datetimes(curve_times, tz_name)
        axis_tzinfo = zoneinfo.ZoneInfo(axis_tz) if axis_tz != "UTC" else dt.timezone.utc
        self.ax.clear()
        self.figure.patch.set_facecolor("#eff3f8")
        self.ax.set_facecolor("#fbfdff")
        for spine in self.ax.spines.values():
            spine.set_color("#9aa4b2")

        y_min, y_max = -20, 90
        twilight_mask = (sun_altitudes < 0.0) & (sun_altitudes > -18.0)
        astro_dark_mask = sun_altitudes <= -18.0
        self.ax.fill_between(
            x_datetimes,
            y_min,
            y_max,
            where=twilight_mask,
            color="#f6c453",
            alpha=0.16,
            interpolate=True,
            zorder=0,
        )
        self.ax.fill_between(
            x_datetimes,
            y_min,
            y_max,
            where=astro_dark_mask,
            color="#0d1b3d",
            alpha=0.22,
            interpolate=True,
            zorder=0,
        )

        self.ax.plot(x_datetimes, curve_alt, color="#1f5aa6", linewidth=2.4, label="Object altitude")
        self.ax.axhline(0.0, color="#616161", linestyle="--", linewidth=1.1, label="Horizon (0 deg)")
        self.ax.axhline(
            target_altitude,
            color="#bc2f2f",
            linestyle="--",
            linewidth=1.3,
            label=f"Target ({target_altitude:.1f} deg)",
        )

        for when, direction in target_events:
            dt_event = self._to_plot_datetimes([when], tz_name)[0][0]
            marker = "^" if direction == "up" else "v"
            self.ax.scatter([dt_event], [target_altitude], marker=marker, color="#d62728", s=70, zorder=5)

        if self.vars["show_rise_set"].get():
            for when, direction in rise_set_events:
                dt_event = self._to_plot_datetimes([when], tz_name)[0][0]
                marker = "^" if direction == "up" else "v"
                self.ax.scatter([dt_event], [0.0], marker=marker, color="#444444", s=55, zorder=5)

        max_idx = int(np.argmax(curve_alt))
        max_alt = float(curve_alt[max_idx])
        max_time = x_datetimes[max_idx]
        self.ax.scatter([max_time], [max_alt], color="#1f8f4d", s=65, zorder=6, label="Daily max alt")
        self.ax.annotate(
            f"Max {max_alt:.1f} deg",
            xy=(max_time, max_alt),
            xytext=(10, 10),
            textcoords="offset points",
            fontsize=9,
            color="#214c36",
            bbox={"boxstyle": "round,pad=0.25", "fc": "#e9f5ee", "ec": "#8fbba2", "alpha": 0.9},
        )

        self.ax.set_title(f"Altitude Plan: {object_name} ({obs_date} UTC)", fontsize=13, color="#19365e", pad=10)
        self.ax.set_ylabel("Altitude (deg)")
        self.ax.set_xlabel(f"Time ({axis_tz})")
        self.ax.set_ylim(y_min, y_max)
        self.ax.grid(True, linestyle=":", alpha=0.45, color="#7e8ea3")

        self.ax.legend_.remove() if self.ax.legend_ else None

        self.ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=axis_tzinfo))
        self.figure.autofmt_xdate()
        self.figure.tight_layout()
        self.canvas.draw()

    def _on_close(self):
        if self._poll_after_id is not None:
            try:
                self.after_cancel(self._poll_after_id)
            except Exception:
                pass
            self._poll_after_id = None
        try:
            plt.close(self.figure)
        except Exception:
            pass
        self.quit()
        self.destroy()


def main():
    app = TimeAltitudeApp()
    app.mainloop()


if __name__ == "__main__":
    main()
