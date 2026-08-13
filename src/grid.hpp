// =============================================================================
// grid.hpp -- uniform structured grid, cell-centred, and the face indexing.
//
// ORIENTATION
//   Cells are stored in raster order to match the GeoTIFF the terrain comes
//   from: i increases EAST, j increases SOUTH, and (x0, y0) is the NORTH-WEST
//   corner of cell (0,0). The solver itself is orientation-agnostic; keeping
//   raster order simply removes a class of flip bugs at the I/O boundary.
//
// FACE INDEXING (staggered / MAC)
//   qx lives on WEST-EAST faces:  (nx+1) x ny
//       fx(i,j) is the face on the WEST edge of cell (i,j).
//       fx(0,j)  = domain west boundary,  fx(nx,j) = domain east boundary.
//       q positive  =>  flow from cell (i-1,j) to cell (i,j)   [eastward]
//
//   qy lives on NORTH-SOUTH faces: nx x (ny+1)
//       fy(i,j) is the face on the NORTH edge of cell (i,j).
//       fy(i,0)  = domain north boundary, fy(i,ny) = domain south boundary.
//       q positive  =>  flow from cell (i,j-1) to cell (i,j)   [southward]
//
//   So a cell's continuity is
//       dh/dt * dx = -( qx(i+1,j) - qx(i,j) ) - ( qy(i,j+1) - qy(i,j) )
// =============================================================================
#pragma once

#include <cstddef>

namespace dwe {

struct Grid
{
    int    nx = 0;
    int    ny = 0;
    double dx = 1.0;     // m, square cells
    double x0 = 0.0;     // m, west edge of column 0
    double y0 = 0.0;     // m, north edge of row 0

    Grid() = default;
    Grid(int nx_, int ny_, double dx_, double x0_ = 0.0, double y0_ = 0.0)
        : nx(nx_), ny(ny_), dx(dx_), x0(x0_), y0(y0_) {}

    // --- cells ----------------------------------------------------------
    inline std::size_t ncell() const { return std::size_t(nx) * std::size_t(ny); }
    inline std::size_t k(int i, int j) const { return std::size_t(j) * nx + i; }
    inline double area() const { return dx * dx; }

    inline double xc(int i) const { return x0 + (i + 0.5) * dx; }
    inline double yc(int j) const { return y0 - (j + 0.5) * dx; }

    // --- faces ----------------------------------------------------------
    inline std::size_t nfx() const { return std::size_t(nx + 1) * std::size_t(ny); }
    inline std::size_t nfy() const { return std::size_t(nx) * std::size_t(ny + 1); }

    inline std::size_t fx(int i, int j) const { return std::size_t(j) * (nx + 1) + i; }
    inline std::size_t fy(int i, int j) const { return std::size_t(j) * nx + i; }

    inline bool in_range(int i, int j) const
    { return i >= 0 && i < nx && j >= 0 && j < ny; }
};

} // namespace dwe
