# 2D-DWE

2D-DWE is a two-dimensional local-inertial overland-flow solver on a uniform raster.
The solver core is C++17 with OpenMP and has no third-party runtime dependency. GDAL
support adds georeferenced terrain, roughness and GeoTIFF result I/O for real-world
and configuration-driven cases.

**Status: Scheme A implemented and verified.**

Start with **[DESIGN.md](DESIGN.md)**; the implementation log is Appendix A.

## Build and verify

Build the dependency-free solver and run the two self-contained analytical checks:

```sh
make -j
make test
```

To enable raster and configuration-file cases, install GDAL development files and build with:

```sh
make clean
make -j WITH_GDAL=1
```

The Environment Agency Test 4 benchmark additionally needs Python, NumPy and Rasterio:

```sh
make test-ea WITH_GDAL=1
```

The qmake project remains available for Qt Creator users:

```sh
qmake6 2D-DWE.pro CONFIG+=with_gdal && make -j
```

## Verification status

| test | criterion | result |
|---|---|---|
| tilted plane | Manning normal depth | 0.2276 m expected, **0.2276 m** (+0.00 %) |
| closed depression | nothing over the rim | **0.00000 m** |
| closed depression | drains to the notch | 8.43 → 8.14 → **8.04 m** at 2/10/40 ks |
| both | mass closure | **< 1e-14** |
| EA 2D benchmark Test 4 (Neelz and Pender 2013) | depths and arrival times at points 1, 3, 5, 6 within the 17 packages' range | **pass** (`cases/ea_test4/RESULT.md`) |

## The problem in one paragraph

Figueredo Wash runs about 3,000 ft north-west of the Coyote Canyon community and sits
*higher* than it, separated by a natural rise. A 1-D HEC-RAS study established that the
100-year water surface stands **3.3 ft above the lowest point of that rise**, over 786 ft of
the 2,620 ft community frontage. What a 1-D model cannot produce is the volume that crosses,
the depth on the community bench, or how far the spill travels — which is what a mitigation
decision needs. This solver produces those, cheaply enough to sweep berm heights and
sensitivities.

It is **not** the model of record. HEC-RAS 7.0.1 is, because a reviewer can audit it. This is
a second opinion and a sweep engine, validated against RAS on the same domain.

## Key design decisions

| | |
|---|---|
| Scheme | **Local inertial** (Bates et al. 2010) — explicit, Courant-limited, GPU-friendly, and better than pure diffusive wave at the Fr ≈ 1 stations. Semi-implicit DW retained as an alternative. |
| Grid | Uniform structured, cell-centred finite volume. Develop at 3 m, produce at **1 m** — at native DEM resolution the whole crest-representation problem disappears. |
| Terrain | 1 m LiDAR, **unsmoothed**. A moving average lowers ridges and raises inverts; on a 3.3 ft margin that bias is not affordable. |
| Wet/dry | Face-based, depth over the higher bed, with a positivity limiter and hysteresis. |
| Acceptance | Mass balance closes to < 0.1 % — or the run is discarded. |

## Layout

```
DESIGN.md     the design note — read this first
src/          (not yet written)
tests/        (not yet written)
cases/c6/     Coyote Canyon configuration
```

Inputs are prepared by `NHA/PythonScripts/ras2d/C6-2d-inputs.py` and live in
`WS3_GIS/NM15-097/outputs_2D/`.
