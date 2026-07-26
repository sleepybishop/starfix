#!/usr/bin/env python3
import os
import sys
import argparse
import json
import math
import struct
import subprocess
from datetime import datetime
from PIL import Image, ImageDraw, ImageFont
import piexif
import piexif.helper

def calculate_gha_aries(utc_time):
    year = utc_time.year
    month = utc_time.month
    day = utc_time.day
    hour = utc_time.hour
    minute = utc_time.minute
    second = utc_time.second
    
    a = (14 - month) // 12
    y = year + 4800 - a
    m = month + 12 * a - 3
    jdn = day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045
    jd = jdn + (hour - 12) / 24.0 + minute / 1440.0 + second / 86400.0
    
    T = (jd - 2451545.0) / 36525.0
    gmst = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * T*T - (T*T*T) / 38710000.0
    gmst = gmst % 360.0
    if gmst < 0:
        gmst += 360.0
    return jd, gmst

def parse_rational(r):
    return float(r[0]) / float(r[1])

def get_exif_gps(gps_ifd):
    if not gps_ifd:
        return None, None
    try:
        lat_ref = gps_ifd[piexif.GPSIFD.GPSLatitudeRef].decode('utf-8')
        lat_val = gps_ifd[piexif.GPSIFD.GPSLatitude]
        lon_ref = gps_ifd[piexif.GPSIFD.GPSLongitudeRef].decode('utf-8')
        lon_val = gps_ifd[piexif.GPSIFD.GPSLongitude]
        
        lat = parse_rational(lat_val[0]) + parse_rational(lat_val[1])/60.0 + parse_rational(lat_val[2])/3600.0
        if lat_ref == 'S':
            lat = -lat
            
        lon = parse_rational(lon_val[0]) + parse_rational(lon_val[1])/60.0 + parse_rational(lon_val[2])/3600.0
        if lon_ref == 'W':
            lon = -lon
            
        return lat, lon
    except KeyError:
        return None, None

def get_exif_tilts(exif_ifd):
    if not exif_ifd:
        return None
    try:
        user_comment_raw = exif_ifd[piexif.ExifIFD.UserComment]
        comment_str = piexif.helper.UserComment.load(user_comment_raw)
        tilts = json.loads(comment_str)
        return tilts
    except (KeyError, ValueError, json.JSONDecodeError):
        return None

def crop_center_to_pgm(image_path, out_pgm_path, crop_size=560):
    img = Image.open(image_path)
    w, h = img.size
    left = (w - crop_size) // 2
    top = (h - crop_size) // 2
    right = left + crop_size
    bottom = top + crop_size
    cropped = img.crop((left, top, right, bottom)).convert('L')
    cropped.save(out_pgm_path)

def load_and_precess_catalog(db_path, year=2026.56):
    with open(db_path, "rb") as f:
        magic = f.read(4)
        if magic != b"SFIX":
            raise ValueError("Invalid database format")
        n_stars, n_entries, b_factor = struct.unpack("<III", f.read(12))
        
        stars = []
        # Each catalog entry is 12 bytes
        for _ in range(n_stars):
            data = f.read(12)
            hip, ra_q, dec_q, mag_q = struct.unpack("<iHHB", data[:9])
            
            # Dequantize J2000
            ra0 = ra_q * (2.0 * math.pi) / 65535.0
            dec0 = dec_q * math.pi / 65535.0 - math.pi / 2.0
            mag = mag_q * 10.0 / 255.0 - 2.0
            
            # Propagate precession
            T = (year - 2000.0) / 100.0
            m = 46.1244 * math.pi / (3600.0 * 180.0)
            n = 20.0431 * math.pi / (3600.0 * 180.0)
            
            d_dec = T * n * math.cos(ra0)
            tan_dec = math.tan(dec0) if abs(dec0) < 1.55 else 0.0
            d_ra = T * (m + n * math.sin(ra0) * tan_dec)
            
            ra = (ra0 + d_ra) % (2.0 * math.pi)
            dec = max(-math.pi/2.0, min(math.pi/2.0, dec0 + d_dec))
            
            stars.append({"hip": hip, "ra": ra, "dec": dec, "mag": mag})
    return stars
def main():
    parser = argparse.ArgumentParser(description="StarFix Photographic Pose & Position Solver")
    parser.add_argument("image", nargs="?", default="pa.jpg", help="Path to raw JPEG photograph")
    args = parser.parse_args()

    image_path = args.image
    if not os.path.exists(image_path):
        print(f"Error: {image_path} not found!")
        sys.exit(1)
        
    print(f"====================================================")
    print(f"STARFIX PIPELINE & FULL-FRAME PLATE SOLVER OVERLAY")
    print(f"====================================================")
    
    # 1. Open image to get dimensions
    with Image.open(image_path) as img:
        w_orig, h_orig = img.size
        
    # 2. Load EXIF metadata
    exif_dict = piexif.load(image_path)
    
    # 2. Extract timestamp and calculate GHA Aries
    date_str = exif_dict["0th"].get(piexif.ImageIFD.DateTime)
    if not date_str:
        date_str = exif_dict["Exif"].get(piexif.ExifIFD.DateTimeOriginal)
        
    if not date_str:
        print("Error: No timestamp found in EXIF!")
        sys.exit(1)
        
    date_str = date_str.decode('utf-8')
    dt = datetime.strptime(date_str, "%Y:%m:%d %H:%M:%S")
    
    # Parse offset (PDT: UTC-7)
    offset_str = exif_dict["Exif"].get(42033)
    offset_hours = -7.0
    if offset_str:
        try:
            offset_str = offset_str.decode('utf-8')
            sign = 1 if offset_str[0] == '+' else -1
            h = int(offset_str[1:3])
            m = int(offset_str[4:6])
            offset_hours = sign * (h + m/60.0)
        except Exception:
            pass
            
    import datetime as dt_mod
    utc_dt = dt - dt_mod.timedelta(hours=offset_hours)
    jd, gha_aries = calculate_gha_aries(utc_dt)
    
    print(f"Time : {date_str} (UTC{'-' if offset_hours < 0 else '+'}{abs(offset_hours):g})")
    
    # 3. Extract GPS coordinates (use as Assumed Position)
    ap_lat, ap_lon = get_exif_gps(exif_dict["GPS"])
    if ap_lat is None or ap_lon is None:
        print("Error: No GPS coordinates found in EXIF!")
        sys.exit(1)
    
    # 4. Extract phone tilts and build gravity vector
    tilts = get_exif_tilts(exif_dict["Exif"])
    if not tilts:
        print("Error: No tilt metadata found in EXIF UserComment!")
        sys.exit(1)
        
    roll = tilts.get("roll_deg", 0.0)
    pitch = tilts.get("pitch_deg", 0.0)
    
    roll_rad = math.radians(roll)
    pitch_rad = math.radians(pitch)
    
    g_x = math.sin(roll_rad)
    g_y = math.sin(pitch_rad)
    g_z = math.sqrt(max(0.0, 1.0 - g_x*g_x - g_y*g_y))
    
    # 4. Extract Focal Length from EXIF (tag 41989 in Exif IFD)
    focal_35 = exif_dict["Exif"].get(piexif.ExifIFD.FocalLengthIn35mmFilm, 24)
    if isinstance(focal_35, bytes):
        focal_35 = int(focal_35.decode('utf-8'))
    elif not isinstance(focal_35, int):
        focal_35 = 24
        
    f_pixels = w_orig * (focal_35 / 36.0)
    fov_deg = 2.0 * math.degrees(math.atan(280.0 / f_pixels))
    print(f"Focal Length (35mm Equiv): {focal_35} mm")
    print(f"Calculated Crop FOV      : {fov_deg:.4f} degrees")
    
    # 5. Crop central region temporarily to data/mock_image.pgm
    os.makedirs("data", exist_ok=True)
    temp_pgm = "data/mock_image.pgm"
    temp_centroids_json = "data/detected_centroids.json"
    
    print(f"\n[1/4] Cropping central 560x560 region to temp file...")
    crop_center_to_pgm(image_path, temp_pgm)
    
    # 6. Run C Centroiding step
    print(f"[2/4] Extracting star centroids using starfix_cli...")
    subprocess.run(["../starfix_cli", "--centroid"], check=True)
    
    # Auto-bootstrap matching stars dynamically using blind identification grid sweep
    print("Auto-bootstrapping lens FOV and star identification blindly...")
    best_fov = None
    low = 0.8 * fov_deg
    high = 1.2 * fov_deg
    # Sweep FOV in 0.1 degree steps
    for test_fov in [low + i * 0.1 for i in range(int((high - low) / 0.1) + 1)]:
        if os.path.exists("data/identified_stars.json"):
            os.remove("data/identified_stars.json")
        cmd = ["../starfix_cli", "--identify", "560", "560", f"{test_fov:.2f}"]
        subprocess.run(cmd, capture_output=True)
        if os.path.exists("data/identified_stars.json"):
            with open("data/identified_stars.json", "r") as f:
                matches = json.load(f)
            if len(matches) >= 3:
                print(f"  Auto-bootstrap succeeded at {test_fov:.2f} degrees (matched {len(matches)} stars)")
                best_fov = test_fov
                break
    if best_fov is None:
        print("Error: Could not blindly identify matching stars in crop. Auto-bootstrap failed.")
        sys.exit(1)
    fov_deg = best_fov
    
    # Optimize lens focal length and FOV dynamically to minimize star projection RMS error
    print("Optimizing lens focal length and FOV dynamically...")
    def run_solver_for_fov(test_fov):
        cmd = [
            "../starfix_cli",
            "--solve-photo",
            f"{gha_aries:.6f}",
            f"{g_x:.6f}",
            f"{g_y:.6f}",
            f"{g_z:.6f}",
            f"{ap_lat:.6f}",
            f"{ap_lon:.6f}",
            f"{test_fov:.6f}"
        ]
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in res.stdout.splitlines():
            if line.startswith("RMS_ERROR:"):
                return float(line.split()[1])
        return float('inf')

    low = 0.8 * fov_deg
    high = 1.2 * fov_deg
    for _ in range(25):
        m1 = low + (high - low) / 3.0
        m2 = high - (high - low) / 3.0
        err1 = run_solver_for_fov(m1)
        err2 = run_solver_for_fov(m2)
        if err1 < err2:
            high = m2
        else:
            low = m1
    optimal_fov = (low + high) / 2.0
    optimal_rms = run_solver_for_fov(optimal_fov)
    print(f"Optimal FOV solved       : {optimal_fov:.4f} degrees (RMS: {optimal_rms:.4f} pixels)")
    fov_deg = optimal_fov

    # 7. Run C solver with EXIF inputs
    print(f"[3/4] Solving camera attitude & terrestrial position...")
    cmd = [
        "../starfix_cli",
        "--solve-photo",
        f"{gha_aries:.6f}",
        f"{g_x:.6f}",
        f"{g_y:.6f}",
        f"{g_z:.6f}",
        f"{ap_lat:.6f}",
        f"{ap_lon:.6f}",
        f"{fov_deg:.6f}"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    print(result.stdout)
    
    # 8. Parse solved attitude matrix R and solved coordinates
    R_flat = None
    solved_lat_str, solved_lon_str = "", ""
    for line in result.stdout.splitlines():
        if line.startswith("ROTATION_MATRIX:"):
            R_flat = list(map(float, line.split()[1:]))
        elif "Solved Latitude" in line:
            solved_lat_str = line.split(":")[-1].strip()
        elif "Solved Longitude" in line:
            solved_lon_str = line.split(":")[-1].strip()
            
    if not R_flat:
        print("Error: Could not parse rotation matrix from C solver output!")
        # Clean up temp files
        if os.path.exists(temp_pgm): os.remove(temp_pgm)
        if os.path.exists(temp_centroids_json): os.remove(temp_centroids_json)
        sys.exit(1)
        
    # Reconstruct 3x3 rotation matrix R
    R = [R_flat[0:3], R_flat[3:6], R_flat[6:9]]
    
    # 9. Generate full-frame visualization
    print(f"[4/4] Generating full-frame ($4000 \\times 3000$) plate overlay...")
    img = Image.open(image_path).convert('RGBA')
    w_orig, h_orig = img.size # should be 4000x3000
    
    # Create transparent overlay for alpha drawing
    overlay = Image.new('RGBA', img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    
    # Crop boundaries
    crop_size = 560
    left = (w_orig - crop_size) // 2
    top = (h_orig - crop_size) // 2
    right = left + crop_size
    bottom = top + crop_size
    
    # Draw green crop box rectangle
    draw.rectangle([left, top, right, bottom], outline=(0, 255, 0, 255), width=5)
    draw.text((left + 10, top + 10), "Crop region (560x560 px)", fill=(0, 255, 0, 255))
    
    # Load and draw detected centroids
    with open(temp_centroids_json, "r") as f:
        centroids = json.load(f)
        
    for idx, c in enumerate(centroids):
        u_crop, v_crop = c['u'], c['v']
        # Map to original full-frame coordinates
        u_orig = left + u_crop
        v_orig = top + v_crop
        # Draw red circle around centroids
        draw.ellipse([u_orig - 15, v_orig - 15, u_orig + 15, v_orig + 15], outline=(255, 0, 0, 255), width=3)
        draw.text((u_orig + 20, v_orig - 10), f"#{idx}", fill=(255, 255, 0, 255))

    # Load and project catalog stars
    print("  Loading star catalog and projecting coordinates...")
    stars = load_and_precess_catalog("../data/starfix_db.bin")
    
    # Focal length of phone lens in pixels
    # (using dynamic fov_deg calculated above)
    fov_rad = math.radians(fov_deg)
    f_len = crop_size / (2.0 * math.tan(fov_rad / 2.0))
    
    for star in stars:
        # 1. Celestial unit vector
        v_cel = [
            math.cos(star["dec"]) * math.cos(star["ra"]),
            math.cos(star["dec"]) * math.sin(star["ra"]),
            math.sin(star["dec"])
        ]
        
        # 2. Camera frame unit vector: v_cam = R * v_cel
        v_cam = [
            R[0][0]*v_cel[0] + R[0][1]*v_cel[1] + R[0][2]*v_cel[2],
            R[1][0]*v_cel[0] + R[1][1]*v_cel[1] + R[1][2]*v_cel[2],
            R[2][0]*v_cel[0] + R[2][1]*v_cel[1] + R[2][2]*v_cel[2]
        ]
        
        # 3. Project if in front of camera
        if v_cam[2] > 0.0:
            x_proj = v_cam[0] / v_cam[2]
            y_proj = v_cam[1] / v_cam[2]
            
            # Map to full image coordinates (origin at center: 2000, 1500)
            u_proj = w_orig / 2.0 + x_proj * f_len
            v_proj = h_orig / 2.0 + y_proj * f_len
            
            if 0.0 <= u_proj <= w_orig and 0.0 <= v_proj <= h_orig:
                # Draw cyan dots (4px wide -> radius 2px) at 50% transparency (alpha 128)
                if star["mag"] <= 6.5: # only draw visible catalog stars
                    draw.ellipse([u_proj - 2, v_proj - 2, u_proj + 2, v_proj + 2], fill=(0, 255, 255, 128))

    # Draw telemetry overlay box
    overlay_w, overlay_h = 950, 480
    box_x, box_y = 50, 50
    draw.rectangle([box_x, box_y, box_x + overlay_w, box_y + overlay_h], fill=(0, 0, 0, 180), outline=(255, 255, 255, 255), width=4)
    
    text_y = box_y + 30
    draw.text((box_x + 30, text_y), "STARFIX NAVIGATION TELEMETRY", fill=(0, 255, 0, 255))
    text_y += 50
    draw.text((box_x + 30, text_y), f"Image Name       : pa.jpg", fill=(255, 255, 255, 255))
    text_y += 40
    draw.text((box_x + 30, text_y), f"UTC Timestamp    : {utc_dt}", fill=(255, 255, 255, 255))
    text_y += 40
    draw.text((box_x + 30, text_y), f"GHA Aries        : {gha_aries:.6f} deg", fill=(255, 255, 255, 255))
    text_y += 40
    draw.text((box_x + 30, text_y), f"Camera Tilts     : Roll = {roll:.4f} deg, Pitch = {pitch:.4f} deg", fill=(255, 255, 255, 255))
    text_y += 50
    draw.text((box_x + 30, text_y), f"Solved Latitude  : {solved_lat_str}", fill=(0, 255, 255, 255))
    text_y += 40
    draw.text((box_x + 30, text_y), f"Solved Longitude : {solved_lon_str}", fill=(0, 255, 255, 255))

    # Alpha composite the overlay onto the original image
    final_img = Image.alpha_composite(img, overlay).convert('RGB')

    # Save output files
    out_jpeg = "solve_overlay.jpg"
    vis_artifact = "/home/joe/.gemini/antigravity-cli/brain/c5b3e84e-183a-494f-854b-e1459382e39b/solve_overlay.jpg"
    
    final_img.save(out_jpeg, quality=85)
    final_img.save(vis_artifact, quality=85)
    
    print(f"\nSuccessfully generated full-frame overlay image!")
    print(f"  Saved locally to: {out_jpeg}")
    print(f"  Saved as conversation artifact to: {vis_artifact}")
    
    # 10. Clean up temporary files in data/
    print("\nCleaning up temporary files from data/ directory...")
    if os.path.exists(temp_pgm): os.remove(temp_pgm)
    if os.path.exists(temp_centroids_json): os.remove(temp_centroids_json)
    print("Done!")

if __name__ == "__main__":
    main()
