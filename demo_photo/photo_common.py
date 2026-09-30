"""Shared photo metadata, cropping, identification, and attitude fitting."""

import json
import math
import os
import shutil
import subprocess
from datetime import datetime, timedelta

import piexif
import piexif.helper
from PIL import Image, ImageOps

CROP_SIZE = 560
MAX_RMS_PIXELS = 3.0
MIN_INLIERS = 3
CAMERA_FRAME = "x_right_y_down_z_forward"


def calculate_gha_aries(utc_time):
    year = utc_time.year
    month = utc_time.month
    day = utc_time.day
    hour = utc_time.hour
    minute = utc_time.minute
    second = utc_time.second + utc_time.microsecond / 1_000_000.0

    a = (14 - month) // 12
    y = year + 4800 - a
    m = month + 12 * a - 3
    jdn = day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045
    jd = jdn + (hour - 12) / 24.0 + minute / 1440.0 + second / 86400.0

    t = (jd - 2451545.0) / 36525.0
    gmst = (
        280.46061837
        + 360.98564736629 * (jd - 2451545.0)
        + 0.000387933 * t * t
        - (t * t * t) / 38710000.0
    )
    return jd, gmst % 360.0


def decimal_year(utc_time):
    start = datetime(utc_time.year, 1, 1)
    end = datetime(utc_time.year + 1, 1, 1)
    return utc_time.year + (utc_time - start).total_seconds() / (end - start).total_seconds()


def parse_rational(value):
    numerator, denominator = value
    if denominator == 0:
        raise ValueError("zero EXIF rational denominator")
    return float(numerator) / float(denominator)


def get_exif_gps(gps_ifd):
    if not gps_ifd:
        return None, None
    try:
        lat_ref = gps_ifd[piexif.GPSIFD.GPSLatitudeRef].decode("ascii").upper()
        lon_ref = gps_ifd[piexif.GPSIFD.GPSLongitudeRef].decode("ascii").upper()
        lat_val = gps_ifd[piexif.GPSIFD.GPSLatitude]
        lon_val = gps_ifd[piexif.GPSIFD.GPSLongitude]
        lat = parse_rational(lat_val[0]) + parse_rational(lat_val[1]) / 60.0
        lat += parse_rational(lat_val[2]) / 3600.0
        lon = parse_rational(lon_val[0]) + parse_rational(lon_val[1]) / 60.0
        lon += parse_rational(lon_val[2]) / 3600.0
    except (KeyError, TypeError, ValueError, UnicodeDecodeError):
        return None, None

    if lat_ref not in ("N", "S") or lon_ref not in ("E", "W"):
        return None, None
    lat = -lat if lat_ref == "S" else lat
    lon = -lon if lon_ref == "W" else lon
    if not -90.0 <= lat <= 90.0 or not -180.0 <= lon <= 180.0:
        return None, None
    return lat, lon


def get_user_comment(exif_ifd):
    if not exif_ifd:
        return None
    try:
        raw = exif_ifd[piexif.ExifIFD.UserComment]
        value = json.loads(piexif.helper.UserComment.load(raw))
        return value if isinstance(value, dict) else None
    except (KeyError, TypeError, ValueError, UnicodeDecodeError, json.JSONDecodeError):
        return None


def parse_exif_utc(exif_dict):
    exif_ifd = exif_dict["Exif"]
    date_value = exif_ifd.get(piexif.ExifIFD.DateTimeOriginal)
    offset_tags = (36881, 36880, 36882)
    if date_value is None:
        date_value = exif_dict["0th"].get(piexif.ImageIFD.DateTime)
        offset_tags = (36880, 36881, 36882)
    if date_value is None:
        raise ValueError("no EXIF timestamp")

    date_text = date_value.decode("ascii") if isinstance(date_value, bytes) else str(date_value)
    dt = datetime.strptime(date_text, "%Y:%m:%d %H:%M:%S")
    subsec = exif_ifd.get(piexif.ExifIFD.SubSecTimeOriginal)
    if subsec:
        subsec_text = subsec.decode("ascii") if isinstance(subsec, bytes) else str(subsec)
        digits = "".join(ch for ch in subsec_text if ch.isdigit())[:6]
        if digits:
            dt = dt.replace(microsecond=int(digits.ljust(6, "0")))

    offset_hours = None
    for tag in offset_tags:
        value = exif_ifd.get(tag)
        if not value:
            continue
        text = value.decode("ascii") if isinstance(value, bytes) else str(value)
        if len(text) != 6 or text[0] not in "+-" or text[3] != ":":
            continue
        try:
            hours = int(text[1:3])
            minutes = int(text[4:6])
        except ValueError:
            continue
        if hours > 14 or minutes > 59 or (hours == 14 and minutes != 0):
            continue
        sign = 1 if text[0] == "+" else -1
        offset_hours = sign * (hours + minutes / 60.0)
        break

    if offset_hours is None:
        configured = os.environ.get("STARFIX_UTC_OFFSET_HOURS")
        try:
            offset_hours = float(configured) if configured is not None else None
        except ValueError:
            offset_hours = None
    if offset_hours is None or not math.isfinite(offset_hours) or not -14.0 <= offset_hours <= 14.0:
        raise ValueError(
            "EXIF timestamp has no valid UTC offset; set STARFIX_UTC_OFFSET_HOURS explicitly"
        )
    return dt - timedelta(hours=offset_hours)


def crop_center_to_pgm(image_path, output_path, crop_size=CROP_SIZE):
    with Image.open(image_path) as source:
        image = ImageOps.exif_transpose(source)
        width, height = image.size
        if crop_size < 1 or width < crop_size or height < crop_size:
            raise ValueError(f"image is smaller than the required {crop_size}x{crop_size} crop")
        left = (width - crop_size) // 2
        top = (height - crop_size) // 2
        image.crop((left, top, left + crop_size, top + crop_size)).convert("L").save(output_path)
        return width, height


def get_cli_path():
    for candidate in ("./starfix_cli", "../starfix_cli"):
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return os.path.abspath(candidate)
    resolved = shutil.which("starfix_cli")
    if resolved:
        return os.path.abspath(resolved)
    raise FileNotFoundError("starfix_cli was not found; run `make starfix_cli`")


def get_db_path():
    for candidate in ("data/starfix_db.bin", "../data/starfix_db.bin"):
        if os.path.isfile(candidate):
            return os.path.abspath(candidate)
    raise FileNotFoundError("data/starfix_db.bin was not found")


def estimate_crop_fov(exif_dict, width, height, crop_size=CROP_SIZE):
    focal_35 = exif_dict["Exif"].get(piexif.ExifIFD.FocalLengthIn35mmFilm, 24)
    if isinstance(focal_35, bytes):
        focal_35 = focal_35.decode("ascii")
    try:
        focal_35 = float(focal_35)
    except (TypeError, ValueError):
        focal_35 = 24.0
    if not math.isfinite(focal_35) or focal_35 <= 0.0:
        focal_35 = 24.0

    # 35 mm equivalent focal lengths preserve diagonal FOV, including on 4:3 sensors.
    focal_pixels = math.hypot(width, height) * focal_35 / math.hypot(36.0, 24.0)
    return 2.0 * math.degrees(math.atan((crop_size / 2.0) / focal_pixels))


def run_command(command, work_dir):
    return subprocess.run(command, cwd=work_dir, capture_output=True, text=True, check=False)


def parse_solver_output(result):
    parsed = {"rms": math.inf, "rotation": None, "inliers": 0, "matches": 0}
    if result.returncode != 0:
        return parsed
    for line in result.stdout.splitlines():
        if line.startswith("RMS_ERROR:"):
            try:
                parsed["rms"] = float(line.split()[1])
            except (IndexError, ValueError):
                pass
        elif line.startswith("INLIERS:"):
            try:
                inliers, matches = line.split()[1].split("/")
                parsed["inliers"] = int(inliers)
                parsed["matches"] = int(matches)
            except (IndexError, ValueError):
                pass
        elif line.startswith("ROTATION_MATRIX:"):
            try:
                values = [float(value) for value in line.split()[1:]]
            except ValueError:
                continue
            if len(values) == 9 and all(math.isfinite(value) for value in values):
                parsed["rotation"] = [values[0:3], values[3:6], values[6:9]]
    return parsed


def solution_is_acceptable(solution):
    return (
        solution["rotation"] is not None
        and math.isfinite(solution["rms"])
        and solution["rms"] <= MAX_RMS_PIXELS
        and solution["inliers"] >= MIN_INLIERS
        and solution["matches"] >= MIN_INLIERS
        and solution["inliers"] / solution["matches"] >= 0.6
    )


def solve_plate(image_path, exif_dict, cli_bin, work_dir, utc_dt,
                zenith_camera=(0.0, 0.0, 1.0), fit_attitude=None):
    """Identify stars and fit attitude; leave the selected crop/matches in work_dir/data."""
    with Image.open(image_path) as source:
        width, height = ImageOps.exif_transpose(source).size
    if min(width, height) < CROP_SIZE:
        raise ValueError(f"image is smaller than the required {CROP_SIZE}x{CROP_SIZE} crop")
    year = decimal_year(utc_dt)
    data_dir = os.path.join(work_dir, "data")
    matches_path = os.path.join(data_dir, "identified_stars.json")
    centroids_path = os.path.join(data_dir, "detected_centroids.json")
    sizes = sorted({min(min(width, height), int(CROP_SIZE * factor))
                    for factor in (1, 1.5, 2, 2.5, 3, 4)})

    for crop_size in sizes:
        crop_center_to_pgm(image_path, os.path.join(data_dir, "mock_image.pgm"), crop_size)
        fov_estimate = estimate_crop_fov(exif_dict, width, height, crop_size)

        def identify(fov):
            return run_command([cli_bin, "--identify", str(crop_size), str(crop_size),
                                f"{fov:.9f}", f"{year:.9f}"], work_dir)

        def solve(fov):
            if fit_attitude is not None:
                return fit_attitude(fov, crop_size)
            return parse_solver_output(run_command([
                cli_bin, "--fit-photo", *(f"{value:.12f}" for value in zenith_camera),
                f"{fov:.9f}", f"{year:.9f}", str(crop_size)], work_dir))

        for threshold in (45, 35):
            detected_result = run_command([cli_bin, "--centroid", str(threshold)], work_dir)
            if detected_result.returncode != 0:
                continue
            with open(centroids_path, encoding="utf-8") as source:
                detected = json.load(source)
            # Reject edge fragments and low-flux noise before selecting the brightest quads.
            centroids = [c for c in detected if c["flux"] >= 100 and
                         15 <= c["u"] <= crop_size - 15 and 15 <= c["v"] <= crop_size - 15]
            print(f"  Crop {crop_size}px, threshold {threshold}: {len(centroids)} usable stars")
            if len(centroids) < 5:
                continue  # four seed stars plus an independent verification star
            with open(centroids_path, "w", encoding="utf-8") as output:
                json.dump(centroids, output)

            low_search = max(2.0, fov_estimate * 0.75)
            high_search = min(90.0, fov_estimate * 1.25)
            step = (high_search - low_search) / 40.0
            best = None
            for i in range(41):
                fov = low_search + i * step
                if identify(fov).returncode != 0:
                    continue
                candidate = solve(fov)
                if not solution_is_acceptable(candidate):
                    continue
                score = (candidate["inliers"], -candidate["rms"])
                if best is None or score > best[0]:
                    with open(matches_path, encoding="utf-8") as source:
                        best = (score, fov, json.load(source), candidate)
            if best is None:
                continue
            with open(matches_path, "w", encoding="utf-8") as output:
                json.dump(best[2], output)

            # Keep the correspondence fixed while refining focal length.
            half_width = max(step, best[1] * 0.025)
            low, high = max(1.0, best[1] - half_width), min(100.0, best[1] + half_width)
            ratio = (math.sqrt(5.0) - 1.0) / 2.0
            x1, x2 = high - ratio * (high - low), low + ratio * (high - low)
            solution1, solution2 = solve(x1), solve(x2)
            for _ in range(18):
                if solution1["rms"] <= solution2["rms"]:
                    high, x2, solution2 = x2, x1, solution1
                    x1 = high - ratio * (high - low)
                    solution1 = solve(x1)
                else:
                    low, x1, solution1 = x1, x2, solution2
                    x2 = low + ratio * (high - low)
                    solution2 = solve(x2)
            fov, solution = min(((best[1], best[3]), (x1, solution1), (x2, solution2)),
                                key=lambda item: item[1]["rms"])
            if solution_is_acceptable(solution):
                return {"crop_size": crop_size, "fov_deg": fov, "solution": solution,
                        "centroids": centroids, "identified": best[2],
                        "width": width, "height": height}
    raise RuntimeError("no verified plate solution after searching central and wider crops")
