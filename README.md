# 2D-DWE

2D-DWE is a two-dimensional local-inertial overland-flow solver on a uniform raster.
The solver core is C++17 with OpenMP and has no third-party runtime dependency. GDAL
support adds georeferenced terrain, roughness and GeoTIFF result I/O for real-world
cases.

**Status: the solver core and two analytical verification cases are implemented.**

Start with **[DESIGN.md](DESIGN.md)** for the numerical formulation and design record.

## Build and verify

Build the dependency-free solver and run both analytical checks:

```sh
make -j
make test
```

The individual checks can also be run directly:

```sh
./bin/2D-DWE plane
./bin/2D-DWE pond
```

To enable real-world raster cases, install GDAL development files and build with:

```sh
make clean
make -j WITH_GDAL=1
```

The qmake project remains available for Qt Creator users:

```sh
qmake6 2D-DWE.pro && make -j
```

Results are written as VTK ImageData (`.vti`) without linking VTK and can be opened
in ParaView. GDAL-enabled cases also write GeoTIFF summaries.

## Verification status

| test | criterion | result |
|---|---|---|
| tilted plane | Manning normal depth | 0.2276 m expected, **0.2276 m** (+0.00 %) |
| closed depression | nothing over the rim | **0.00000 m** |
| closed depression | drains to the notch | 8.43 → 8.14 → **8.04 m** at 2/10/40 ks |
| both | mass closure | **< 1e-14** |

## Original application

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

## Repository layout

```
DESIGN.md     numerical design and implementation record
src/          solver, command-line driver and optional GDAL I/O
tests/        executable analytical verification checks
cases/c6/     Coyote Canyon input manifest and roughness lookup
tools/        case preparation, plotting and animation utilities
```

Large terrain, forcing and result rasters are intentionally not stored in Git. The
Coyote Canyon manifest documents the external inputs used for that application.

## Reproducibility scope

The `plane` and `pond` cases are self-contained and exercised in continuous integration.
The site cases currently require external GIS and hydrologic inputs and retain some
project-specific paths in the command-line driver and tools. Making those paths fully
configuration-driven is planned work; published studies should archive their exact input
bundle separately and identify the solver commit used.
