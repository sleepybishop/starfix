# StarFix

StarFix is a high-performance, bare-metal, zero-allocation C star tracker and celestial navigation pipeline designed for real-time attitude determination and position estimation on embedded devices (e.g. Raspberry Pi Zero W), tailormade for hobbyists, sailors, and off-grid navigators.

## Key Features

* **Zero-Allocation Runtime**: Implemented in pure C99 with a strict memory arena allocator (`arena`) limited to a deterministic 10MB RAM footprint.
* **Dual Star-ID Engine**:
  * **Lost-in-Space Mode**: Accelerated 4-star quad ratio pattern matching.
  * **Tracking Mode**: Fast nearest-neighbor matching using attitude hints to bypass hashing overhead.
* **Physical Correction Suite**: Relativistic stellar aberration, atmospheric refraction, and epoch precession corrections.
* **Integrated State Estimation**: Real-time trajectory smoothing combining a high-rate Extended Kalman Filter (EKF) and a Factor Graph Optimizer (using `nanoqsp`).
* **Terrestrial Position Solver**: Resolves exact Earth coordinates (latitude and longitude) by combining camera pointing vectors with gravity vectors.

## Quick Start

### Build the Software
Compile the CLI and test suite:
```bash
make
```

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
Run the end-to-end simulated camera pipeline and generate diagnostic visualizations:
```bash
make demo
```

## Prior Art Footnotes

StarFix's pattern matching algorithms build upon several excellent open-source projects:
* **TETRA**[^1]: The original Python-based 4-star quad matching algorithm developed by the European Space Agency.
* **Cedar-Solve**[^2]: An optimized fork of TETRA introducing uniform database pattern density and separate Rust centroid filtering.
* **Tetra3rs**[^3]: A high-performance compiled Rust port of TETRA.

[^1]: [ESA TETRA3 (Python)](https://github.com/esa/tetra3)
[^2]: [Cedar-Solve (Python/Rust)](https://github.com/smroid/cedar-solve)
[^3]: [Tetra3rs (Rust)](https://github.com/ssmichael1/tetra3rs)

## License
Licensed under the MIT License. See [LICENSE](LICENSE) for details.
