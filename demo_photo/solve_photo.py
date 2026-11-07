#!/usr/bin/env python3
"""Solve a photo using a measured camera zenith and UTC time."""

import argparse
import glob
import math
import os
import sys
import tempfile

import piexif
from PIL import Image, ImageDraw, ImageOps

if __package__:
    from . import photo_common
else:
    import photo_common


def get_camera_zenith(tilts):
    frame = tilts.get("camera_frame")
    if frame is not None and frame != photo_common.CAMERA_FRAME:
        raise ValueError("camera zenith must use +X right, +Y down, +Z forward")
    zenith_camera = tilts.get("zenith_camera")
    if isinstance(zenith_camera, list) and len(zenith_camera) == 3:
        g_x, g_y, g_z = (float(value) for value in zenith_camera)
        g_norm = math.sqrt(g_x * g_x + g_y * g_y + g_z * g_z)
        if not math.isfinite(g_norm) or abs(g_norm - 1.0) > 1e-3:
            raise ValueError("EXIF zenith_camera is not a unit vector")
        return g_x / g_norm, g_y / g_norm, g_z / g_norm
    if zenith_camera is not None:
        raise ValueError("zenith_camera must contain three components")
    roll, pitch = float(tilts.get("roll_deg", 0.0)), float(tilts.get("pitch_deg", 0.0))
    if not math.isfinite(roll) or not math.isfinite(pitch):
        raise ValueError("non-finite legacy tilt component")
    g_x, g_y = math.sin(math.radians(roll)), math.sin(math.radians(pitch))
    g_z_sq = 1.0 - g_x * g_x - g_y * g_y
    if g_z_sq < -1e-9:
        raise ValueError("legacy EXIF tilt components do not form a unit vector")
    return g_x, g_y, math.sqrt(max(0.0, g_z_sq))


def refract_camera_ray(camera, zenith):
    cos_theta = sum(c * z for c, z in zip(camera, zenith))
    cos_theta = max(-1.0, min(1.0, cos_theta))
    theta_t = math.acos(cos_theta)
    alt_t = math.pi / 2.0 - theta_t
    if theta_t <= 1e-6 or alt_t <= 0.0:
        return camera
    alt_deg = math.degrees(alt_t)
    if alt_deg >= 15.0:
        refraction = 7.888888e-5 * 1010.0 / (283.0 * math.tan(alt_t))
    else:
        num = ((2e-5 * alt_deg + 1.96e-2) * alt_deg + 1.594e-1) * 1010.0
        den = 283.0 * ((8.45e-2 * alt_deg + 5.05e-1) * alt_deg + 1.0)
        refraction = math.radians(num / den)
    theta_a = max(0.0, theta_t - refraction)
    tangent = [c - cos_theta * z for c, z in zip(camera, zenith)]
    norm = math.sqrt(sum(t * t for t in tangent))
    if norm <= 1e-12:
        return camera
    return [math.cos(theta_a) * z + math.sin(theta_a) * (t / norm)
            for z, t in zip(zenith, tangent)]


def solve_photo_overlay(image_path, zenith_camera=None):
    image_path = os.path.abspath(image_path)
    exif_dict = piexif.load(image_path)
    utc_dt = photo_common.parse_exif_utc(exif_dict)
    _, gha_aries = photo_common.calculate_gha_aries(utc_dt)
    metadata = photo_common.get_user_comment(exif_dict["Exif"]) or {}
    if zenith_camera is not None:
        metadata = {"zenith_camera": list(zenith_camera)}
    elif not any(key in metadata for key in ("zenith_camera", "roll_deg", "pitch_deg")):
        print("  No tilt measurement: assuming the image center points exactly at zenith")
    zenith = get_camera_zenith(metadata)
    cli_bin, db_file = photo_common.get_cli_path(), photo_common.get_db_path()
    print(f"\nProcessing {os.path.basename(image_path)} at UTC {utc_dt.isoformat()}")
    with tempfile.TemporaryDirectory(prefix="starfix-photo-") as work_dir:
        os.mkdir(os.path.join(work_dir, "data"))
        os.symlink(db_file, os.path.join(work_dir, "data", "starfix_db.bin"))
        plate = photo_common.solve_plate(image_path, exif_dict, cli_bin, work_dir,
                                     utc_dt, zenith_camera=zenith)
        command = [cli_bin, "--solve-photo", f"{gha_aries:.9f}",
                   *(f"{component:.12f}" for component in zenith),
                   f"{plate['fov_deg']:.9f}",
                   f"{photo_common.decimal_year(utc_dt):.9f}", str(plate["crop_size"])]
        result = photo_common.run_command(command, work_dir)
        solution = photo_common.parse_solver_output(result)
        if not photo_common.solution_is_acceptable(solution):
            raise RuntimeError(f"photo position failed quality checks:\n{result.stdout}{result.stderr}")
        position = {}
        for line in result.stdout.splitlines():
            if "Solved Latitude" in line:
                position["latitude"] = float(line.split(":")[1].split()[0])
            elif "Solved Longitude" in line:
                position["longitude"] = float(line.split(":")[1].split()[0])
        if len(position) != 2:
            raise RuntimeError("solver did not return a complete position")
        print(result.stdout)

    with Image.open(image_path) as source:
        image = ImageOps.exif_transpose(source).convert("RGBA")
    width, height = image.size
    overlay = Image.new("RGBA", image.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    size = plate["crop_size"]
    left, top = (width - size) // 2, (height - size) // 2
    draw.rectangle((left, top, left + size, top + size), outline="lime", width=3)
    focal = size / (2.0 * math.tan(math.radians(plate["fov_deg"]) / 2.0))
    rotation = solution["rotation"]
    for match in plate["identified"]:
        centroid = plate["centroids"][match["centroid_idx"]]
        x, y = left + centroid["u"], top + centroid["v"]
        draw.ellipse((x - 12, y - 12, x + 12, y + 12), outline="red", width=2)
        ra, dec = match["ra"], match["dec"]
        celestial = (math.cos(dec)*math.cos(ra), math.cos(dec)*math.sin(ra), math.sin(dec))
        camera = refract_camera_ray([sum(rotation[row][column]*celestial[column]
                                         for column in range(3)) for row in range(3)], zenith)
        if camera[2] <= 0:
            continue
        x = left + size/2 + focal * camera[0]/camera[2]
        y = top + size/2 + focal * camera[1]/camera[2]
        draw.ellipse((x - 4, y - 4, x + 4, y + 4), fill="cyan")
        draw.text((x + 16, y), f"HIP {match['hip']}", fill="white")
    text = ["STARFIX NAVIGATION", os.path.basename(image_path), f"UTC: {utc_dt.isoformat()}",
            f"Camera zenith (+X right, +Y down, +Z forward): {zenith}",
            f"Latitude: {position['latitude']:.6f}  Longitude: {position['longitude']:.6f}",
            f"Crop: {size}px  RMS: {solution['rms']:.3f}px  Inliers: {solution['inliers']}/{solution['matches']}"]
    draw.rectangle((20, 20, min(width - 20, 1050), 180), fill=(0, 0, 0, 210))
    for i, line in enumerate(text):
        draw.text((30, 30 + 23*i), line, fill="white")
    out_jpeg = os.path.splitext(image_path)[0] + "_overlay.jpg"
    Image.alpha_composite(image, overlay).convert("RGB").save(out_jpeg, quality=85)
    print(f"Saved {out_jpeg}")
    return out_jpeg


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("images", nargs="*")
    parser.add_argument("--zenith", type=float, nargs=3, metavar=("X", "Y", "Z"),
                        help="measured unit zenith vector: +X right, +Y down, +Z forward")
    arguments = parser.parse_args(argv)
    search_dir = "demo_photo" if os.path.isdir("demo_photo") else "."
    images = arguments.images or [path for path in sorted(glob.glob(os.path.join(search_dir, "*.jpg")))
                                  if not path.endswith("_overlay.jpg")]
    if not images:
        parser.error("no JPEG images found")
    failed = 0
    for path in images:
        try:
            solve_photo_overlay(path, zenith_camera=arguments.zenith)
        except (OSError, ValueError, RuntimeError) as error:
            print(f"Failed {path}: {error}", file=sys.stderr)
            failed += 1
    return 2 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
