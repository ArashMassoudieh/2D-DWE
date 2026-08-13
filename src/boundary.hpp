// =============================================================================
// boundary.hpp -- domain-edge conditions (DESIGN 6).
//
// Three kinds:
//   WALL          closed face, q = 0.                    (default everywhere)
//   INFLOW        total discharge Q(t) spread over the face cells.
//   NORMAL_DEPTH  free outfall, q = (1/n) h^(5/3) sqrt(S_bed).
//
// A boundary is a run of edge faces on one side of the rectangle. Faces are
// addressed through the staggered indexing in grid.hpp; nothing here touches
// interior faces.
// =============================================================================
#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include "grid.hpp"

namespace dwe {

enum class Side  { West, East, North, South };
enum class BcType { Wall, Inflow, NormalDepth };

// A time series, linearly interpolated. Constant outside its range -- a
// hydrograph that runs out should hold its last value, not fall to zero and
// silently drain the model.
struct Series
{
    std::vector<double> t;   // s
    std::vector<double> v;

    static Series constant(double value) { Series s; s.t = {0.0}; s.v = {value}; return s; }

    double at(double tq) const
    {
        if (t.empty()) return 0.0;
        if (tq <= t.front()) return v.front();
        if (tq >= t.back())  return v.back();
        const auto it = std::upper_bound(t.begin(), t.end(), tq);
        const std::size_t b = std::size_t(it - t.begin());
        const std::size_t a = b - 1;
        const double w = (tq - t[a]) / (t[b] - t[a]);
        return v[a] + w * (v[b] - v[a]);
    }
};

struct Boundary
{
    Side    side  = Side::West;
    BcType  type  = BcType::Wall;

    // index range along the side, [begin, end). For West/East these are j
    // (rows); for North/South they are i (columns). Empty => the whole side.
    int     begin = -1;
    int     end   = -1;

    Series  Q;              // m^3/s, total across the boundary (Inflow)
    double  slope = 0.001;  // bed slope for NormalDepth

    std::string name = "bc";

    void resolve(const Grid& g)
    {
        const int n = (side == Side::West || side == Side::East) ? g.ny : g.nx;
        if (begin < 0) begin = 0;
        if (end   < 0) end   = n;
        begin = std::max(0, std::min(begin, n));
        end   = std::max(begin, std::min(end, n));
    }

    int count() const { return end - begin; }
};

} // namespace dwe
