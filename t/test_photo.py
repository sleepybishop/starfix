#!/usr/bin/env python3
import importlib.util
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from datetime import datetime
from pathlib import Path

import piexif
import piexif.helper
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "demo_photo"))


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module("photo_common", "demo_photo/photo_common.py")
photo = load_module("solve_photo", "demo_photo/solve_photo.py")


class CompletedResult:
    def __init__(self, stdout, returncode=0):
        self.stdout = stdout
        self.returncode = returncode


class PhotoTests(unittest.TestCase):
    def test_fov_is_invariant_under_portrait_landscape_rotation(self):
        exif = {"Exif": {piexif.ExifIFD.FocalLengthIn35mmFilm: 24}}
        portrait = common.estimate_crop_fov(exif, 3000, 4000)
        landscape = common.estimate_crop_fov(exif, 4000, 3000)
        self.assertAlmostEqual(portrait, landscape, places=12)

    def test_real_photo_geometry_and_mirror_rejection(self):
        with open(ROOT / "t/fixtures/photo_centroids.json", encoding="utf-8") as source:
            fixture = json.load(source)
        latitude, longitude = fixture["known_position"]
        for case in fixture["photos"]:
            with self.subTest(image=case["name"]), tempfile.TemporaryDirectory() as directory:
                data = Path(directory) / "data"
                data.mkdir()
                os.symlink(ROOT / "data/starfix_db.bin", data / "starfix_db.bin")
                centroids = case["centroids"]
                (data / "detected_centroids.json").write_text(json.dumps(centroids))
                utc = datetime.fromisoformat(case["utc"].replace("Z", ""))
                _, gha = common.calculate_gha_aries(utc)
                year = common.decimal_year(utc)

                def run(*arguments):
                    return subprocess.run([str(ROOT / "starfix_cli"), *map(str, arguments)],
                                          cwd=directory, capture_output=True, text=True)

                result = run("--identify", case["size"], case["size"], case["fov_deg"], year)
                self.assertEqual(result.returncode, 0, result.stdout)
                arguments = ("--fit-photo", 0, 0, 1,
                             case["fov_deg"], year, case["size"])
                result = run(*arguments)
                self.assertEqual(result.returncode, 0, result.stdout)
                solution = common.parse_solver_output(result)
                self.assertTrue(common.solution_is_acceptable(solution), result.stdout)
                self.assertGreaterEqual(solution["inliers"], 5)
                rotation = solution["rotation"]
                # Compare to independent full-image astrometry, not the simulator's projection.
                lat, lon = case["reference_center_lat_lon"]
                dec, ra = math.radians(lat), math.radians(lon + gha)
                expected = (math.cos(dec)*math.cos(ra), math.cos(dec)*math.sin(ra), math.sin(dec))
                dot = sum(a*b for a, b in zip(rotation[2], expected))
                error = math.degrees(math.acos(max(-1.0, min(1.0, dot))))
                self.assertLess(error, 0.08)
                dec, ra = math.radians(latitude), math.radians(longitude + gha)
                zenith = (math.cos(dec)*math.cos(ra), math.cos(dec)*math.sin(ra), math.sin(dec))
                z_camera = sum(a*b for a, b in zip(rotation[2], zenith))
                tilt = math.degrees(math.acos(max(-1.0, min(1.0, z_camera))))
                self.assertAlmostEqual(tilt, case["reference_tilt_deg"], delta=0.08)

                # The same matched stars reflected vertically cannot be explained by a rotation.
                mirrored = [dict(c, v=case["size"] - c["v"]) for c in centroids]
                (data / "detected_centroids.json").write_text(json.dumps(mirrored))
                result = run(*arguments)
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertNotIn("ROTATION_MATRIX:", result.stdout)

    def test_solver_output_requires_complete_quality_data(self):
        result = CompletedResult(
            "RMS_ERROR: 0.25 pixels\n"
            "INLIERS: 4/5\n"
            "ROTATION_MATRIX: 1 0 0 0 1 0 0 0 1\n"
        )
        solution = common.parse_solver_output(result)
        self.assertTrue(common.solution_is_acceptable(solution))
        solution["inliers"] = 2
        self.assertFalse(common.solution_is_acceptable(solution))
        solution["inliers"] = 4
        solution["rms"] = math.inf
        self.assertFalse(common.solution_is_acceptable(solution))

    def test_exif_time_uses_offset_and_subseconds(self):
        exif = {
            "Exif": {
                piexif.ExifIFD.DateTimeOriginal: b"2024:07:09 22:30:15",
                piexif.ExifIFD.SubSecTimeOriginal: b"25",
                36881: b"-07:00",
            },
            "0th": {},
        }
        self.assertEqual(
            common.parse_exif_utc(exif), datetime(2024, 7, 10, 5, 30, 15, 250000)
        )

    def test_negative_zenith_hemisphere_round_trips(self):
        expected = [0.2, -0.3, -math.sqrt(0.87)]
        actual = photo.get_camera_zenith({"zenith_camera": expected})
        for expected_value, actual_value in zip(expected, actual):
            self.assertAlmostEqual(expected_value, actual_value, places=12)

    def test_exif_orientation_is_applied_before_crop(self):
        with tempfile.TemporaryDirectory() as directory:
            source_path = os.path.join(directory, "portrait.jpg")
            output_path = os.path.join(directory, "crop.pgm")
            exif = {"0th": {piexif.ImageIFD.Orientation: 6}, "Exif": {}, "GPS": {}, "1st": {}}
            Image.new("RGB", (600, 800), "black").save(source_path, exif=piexif.dump(exif))
            width, height = common.crop_center_to_pgm(source_path, output_path)
            self.assertEqual((width, height), (800, 600))
            with Image.open(output_path) as crop:
                self.assertEqual(crop.size, (common.CROP_SIZE, common.CROP_SIZE))

    def test_sparse_center_uses_wider_crop(self):
        with open(ROOT / "t/fixtures/photo_centroids.json", encoding="utf-8") as source:
            fixture = json.load(source)
        case = fixture["photos"][0]
        with tempfile.TemporaryDirectory() as directory:
            # Render measured star positions; the central 560px region is too sparse.
            image = Image.new("L", (case["size"], case["size"]), 0)
            draw = ImageDraw.Draw(image)
            for centroid in case["centroids"]:
                x, y = round(centroid["u"]), round(centroid["v"])
                draw.ellipse((x-2, y-2, x+2, y+2), fill=220)
            path = Path(directory) / "sparse.jpg"
            utc = datetime.fromisoformat(case["utc"].replace("Z", ""))
            focal_mm = round(case["size"] / (2*math.tan(math.radians(case["fov_deg"])/2))
                             * math.hypot(36, 24) / math.hypot(case["size"], case["size"]))
            exif = {"0th": {}, "Exif": {
                piexif.ExifIFD.DateTimeOriginal: utc.strftime("%Y:%m:%d %H:%M:%S").encode(),
                36881: b"+00:00", piexif.ExifIFD.FocalLengthIn35mmFilm: focal_mm},
                "GPS": {}, "1st": {}}
            image.save(path, quality=100, exif=piexif.dump(exif))
            original = path.read_bytes()
            data = Path(directory) / "data"
            data.mkdir()
            os.symlink(ROOT / "data/starfix_db.bin", data / "starfix_db.bin")
            result = common.solve_plate(path, exif, str(ROOT / "starfix_cli"), directory, utc)
            self.assertEqual(result["crop_size"], 840)
            self.assertEqual(path.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
