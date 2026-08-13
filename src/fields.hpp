// =============================================================================
// fields.hpp -- state, stored as a structure of arrays.
//
// Flat std::vector<double>, no cell objects, no virtual dispatch anywhere in
// the inner loops. The same loop bodies become GPU kernel bodies later without
// restructuring (DESIGN 9).
// =============================================================================
#pragma once

#include <vector>
#include <algorithm>
#include <cmath>

#include "grid.hpp"
#include "params.hpp"

namespace dwe {

struct Fields
{
    Grid g;

    std::vector<double> z;    // m, bed elevation      [ncell]
    std::vector<double> h;    // m, water depth        [ncell]
    std::vector<double> n;    // Manning n             [ncell]
    std::vector<double> qx;   // m^2/s, unit discharge [nfx]
    std::vector<double> qy;   // m^2/s, unit discharge [nfy]

    // scratch, allocated once
    std::vector<double> alpha;   // per-cell flux limiter               [ncell]
    std::vector<char>   active;  // 1 if the cell is wet                [ncell]

    explicit Fields(const Grid& grid, double n_default)
        : g(grid),
          z(grid.ncell(), 0.0),
          h(grid.ncell(), 0.0),
          n(grid.ncell(), n_default),
          qx(grid.nfx(), 0.0),
          qy(grid.nfy(), 0.0),
          alpha(grid.ncell(), 1.0),
          active(grid.ncell(), 0)
    {}

    // Water-surface elevation. Not stored: it is z+h and keeping a third array
    // in sync is a bug waiting to happen.
    inline double H(std::size_t k) const { return z[k] + h[k]; }

    double max_depth() const
    {
        double m = 0.0;
        for (double v : h) m = std::max(m, v);
        return m;
    }

    double volume() const
    {
        // Kahan summation: at 10^7 cells a naive sum loses enough precision to
        // muddy a 0.1 % mass-balance test.
        double s = 0.0, c = 0.0;
        for (double v : h) { const double y = v - c; const double t = s + y; c = (t - s) - y; s = t; }
        return s * g.area();
    }
};

// -----------------------------------------------------------------------------
// THE FACE-DEPTH RULE  (DESIGN 4.1)
//
// The single most important function in the code. Never average depths across a
// face: use the depth over the HIGHER of the two beds. This is what stops water
// climbing a rise. An averaged face depth leaks water onto dry high ground and
// the leak is invisible in the output -- it looks like a plausible flood.
// -----------------------------------------------------------------------------
inline double face_depth(double zL, double hL, double zR, double hR)
{
    const double zf = std::max(zL, zR);
    const double Hf = std::max(zL + hL, zR + hR);
    return Hf - zf;
}

} // namespace dwe
