# StarFix

[![CI](https://github.com/sleepybishop/starfix/actions/workflows/ci.yml/badge.svg)](https://github.com/sleepybishop/starfix/actions/workflows/ci.yml)

StarFix is a high-performance, bare-metal, star tracker and celestial navigation pipeline designed for real-time attitude determination and position estimation on embedded devices for hobbyists, sailors, and off-grid navigators.

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

## Prior Art

StarFix's pattern matching algorithms build upon several excellent open-source projects:
* **TETRA**[^1]: The original Python-based 4-star quad matching algorithm developed by the European Space Agency.
* **Cedar-Solve**[^2]: An optimized fork of TETRA introducing uniform database pattern density and separate Rust centroid filtering.
* **Tetra3rs**[^3]: A high-performance compiled Rust port of TETRA.

[^1]: [ESA TETRA3 (Python)](https://github.com/esa/tetra3)
[^2]: [Cedar-Solve (Python/Rust)](https://github.com/smroid/cedar-solve)
[^3]: [Tetra3rs (Rust)](https://github.com/ssmichael1/tetra3rs)
