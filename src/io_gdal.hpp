// =============================================================================
// io_gdal.hpp -- raster I/O. Compiled only with CONFIG+=with_gdal.
//
// Kept behind a flag and behind a narrow interface so that GDAL never appears
// inside a timestep. Rasters are read once into the flat arrays the solver
// already uses, and results are written on an output cadence. Nothing in the
// hot loops knows GDAL exists (DESIGN 9).
// =============================================================================
#pragma once

#ifdef DWE_WITH_GDAL

#include <string>
#include <vector>

#include "grid.hpp"

namespace dwe {

// How a fine source raster is collapsed onto a coarser model cell.
//
// This choice is not cosmetic for a barrier problem. Measured on the C6 divide,
// aggregating 1 m terrain by MEAN lowers the crest by 0.26 ft (5 m cells) to
// 0.95 ft (15 m cells) and raises channel inverts by a similar amount -- a
// one-directional bias, always toward more spill. MIN preserves channels and
// destroys ridges; MAX does the reverse. Running at the native 1 m removes the
// question, which is the reason the design says to produce results there.
enum class Agg { Mean, Min, Max };

struct RasterMeta
{
    int    nx = 0, ny = 0;
    double dx = 0.0;
    double x0 = 0.0, y0 = 0.0;   // north-west corner
    double nodata = 0.0;
    bool   has_nodata = false;
    std::string crs;
};

bool raster_meta(const std::string& path, RasterMeta& m);

// Read `path` and aggregate it onto a model grid of cell size `target_dx`.
// The grid is DERIVED from the raster (extent and origin), so the model domain
// is exactly the raster's footprint. `valid[k] = 0` marks nodata cells, which
// the caller should treat as permanently dry walls.
bool load_terrain(const std::string& path,
                  double target_dx,
                  Agg agg,
                  Grid& g,
                  std::vector<double>& z,
                  std::vector<char>& valid,
                  std::string* err = nullptr);

// Sample a categorical raster (NLCD codes) BY COORDINATE onto an existing grid
// and map each code to a Manning n via `csv` (columns: code,name,n).
// Coordinate sampling, not index arithmetic: the land-cover raster has its own
// resolution and origin and will not line up with the model grid.
bool load_manning(const std::string& landcover_path,
                  const std::string& csv_path,
                  const Grid& g,
                  std::vector<double>& n,
                  double n_default,
                  std::string* err = nullptr);

// Sample a CONTINUOUS Manning-n raster (float, n values directly -- not land-cover
// codes) BY COORDINATE onto an existing grid. Used when the roughness field is
// built alongside a carved terrain and already carries the channel lining, which
// no land-cover product can supply.
bool load_manning_direct(const std::string& n_path,
                         const Grid& g,
                         std::vector<double>& n,
                         double n_default,
                         std::string* err = nullptr);

// Write a cell field as a GeoTIFF on the model grid.
// Any single-band float raster sampled onto the grid at cell centres (nearest),
// with no range check; cells outside it or nodata get `fill`. Used for the
// initial depth of a hot start.
bool load_field(const std::string& path,
                const Grid& g,
                std::vector<double>& a,
                double fill,
                std::string* err = nullptr);

bool write_gtiff(const std::string& path,
                 const Grid& g,
                 const std::vector<double>& a,
                 const std::string& crs_wkt_or_epsg = "EPSG:26912",
                 double nodata = -9999.0);

} // namespace dwe

#endif // DWE_WITH_GDAL
