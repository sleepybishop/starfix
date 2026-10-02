# StarFix

[![CI](https://github.com/sleepybishop/starfix/actions/workflows/ci.yml/badge.svg)](https://github.com/sleepybishop/starfix/actions/workflows/ci.yml)

StarFix is a high-performance, bare-metal, star tracker and celestial navigation pipeline designed for real-time attitude determination and position estimation on embedded devices for hobbyists, sailors, and off-grid navigators.

## Quick Start

### Build the Software
Build the photo CLI and reusable `libstarfix.a`:
```bash
make
```

The default build contains centroiding, star identification, attitude fitting, and
position solving. Navigation support is optional:
```bash
make navigation   # EKF, factor graph, NMEA library, and navigation example
make simulations  # simulated camera/position and navigation examples
```
`make check` builds and tests both the core and optional examples. Simulated
attitude/position commands live in `examples/starfix_sim`; trajectory demos live
in `examples/starfix_navigation`. The photo CLI does not emit NMEA files.

### Rebuild the Database
Download the HYG catalog and generate the optimized guide database:
```bash
make db
```

### Run Tests and Demos
Execute the test runner:
```bash
make check
```
Run the included photo through the solver and generate a diagnostic overlay:
```bash
make demo
```

### Photo solving

Install the photo tools' dependencies with `python3 -m pip install pillow piexif`.
The camera frame is right-handed: **+X image right, +Y image down, +Z forward**,
after applying EXIF orientation. IMU zenith measurements must be transformed into
this same camera frame.

The included `demo_photo/demo.jpg` was taken at 63.78402778° N, 18.01597222° W
on September 21, 2026 at 23:47:51.074 UTC. Run the photo demo with `make demo`.
Its supplied zenith vector was inferred from that known location, so recovering
the location demonstrates consistency, not independent geolocation accuracy.

To calculate a position using an independently measured camera zenith:
```bash
python3 demo_photo/solve_photo.py --zenith 0 0 1 photo.jpg
```
The position is computed directly from the star attitude, camera zenith, and UTC
time. No initial position or EXIF GPS is required. Longitude is wrapped to
[-180°, 180°); at an exact pole its undefined value is reported as zero.
With no tilt metadata or `--zenith`, the tool explicitly assumes the image center
points straight up. The timestamp must include a UTC offset, or
`STARFIX_UTC_OFFSET_HOURS` must be set.

The tool searches wider central crops when the default 560-pixel crop cannot be
solved, verifies each catalog quad with at least a fifth star, and reports residuals
in original image pixels. It writes an overlay next to the source JPEG.

`make check` includes compact centroid fixtures from two real photos, independent
sky-direction references, reflected-image rejection, and a sparse-center crop
fallback test. The original JPEGs are not needed to run these tests.

## Prior Art

StarFix's pattern matching algorithms build upon several excellent open-source projects:
* **TETRA**[^1]: The original Python-based 4-star quad matching algorithm developed by the European Space Agency.
* **Cedar-Solve**[^2]: An optimized fork of TETRA introducing uniform database pattern density and separate Rust centroid filtering.
* **Tetra3rs**[^3]: A high-performance compiled Rust port of TETRA.

[^1]: [ESA TETRA3 (Python)](https://github.com/esa/tetra3)
[^2]: [Cedar-Solve (Python/Rust)](https://github.com/smroid/cedar-solve)
[^3]: [Tetra3rs (Rust)](https://github.com/ssmichael1/tetra3rs)
