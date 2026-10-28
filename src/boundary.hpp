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
#include <fstream>
#include <string>
#include <algorithm>
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

    // Read a two-column CSV "t_s,value"; '#' comments and a header line are
    // skipped. Returns an empty Series if the file cannot be read, so a caller
    // can fall back rather than silently run with no inflow at all.
    static Series from_csv(const std::string& path)
    {
        Series s;
        std::ifstream f(path);
        if (!f) return s;
        std::string line;
        while (std::getline(f, line)) {
            if (line.empty() || line[0] == '#') continue;
            const std::size_t c = line.find(',');
            if (c == std::string::npos) continue;
            try {
                const double a = std::stod(line.substr(0, c));
                const double b = std::stod(line.substr(c + 1));
                s.t.push_back(a); s.v.push_back(b);
            } catch (...) { continue; }        // the header line lands here
        }
        return s;
    }

    double peak()  const { return v.empty() ? 0.0 : *std::max_element(v.begin(), v.end()); }
    double span()  const { return t.empty() ? 0.0 : t.back(); }
    double t_peak() const
    {
        if (v.empty()) return 0.0;
        return t[std::size_t(std::max_element(v.begin(), v.end()) - v.begin())];
    }

    // Trapezoidal volume [m^3] -- lets the run report how much water the
    // hydrograph should deliver, which the mass balance can then be checked against.
    double volume() const
    {
        double s2 = 0.0;
        for (std::size_t i = 1; i < t.size(); ++i)
            s2 += 0.5 * (v[i] + v[i - 1]) * (t[i] - t[i - 1]);
        return s2;
    }

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

    // Inflow only: distribute Q across the run in proportion to the local
    // conveyance, max(H_edge - z, 0)^(5/3), with H_edge the highest water surface
    // on the run, instead of uniformly. A uniform split forces the whole flood
    // through every cell of the run alike -- with a run cut to the low-flow
    // channel (Rock Creek: 20 m) the overbank flow is pushed through the channel
    // and mounds up a metre at the edge; with a run spanning the valley it is
    // poured onto dry terraces. Weighted, the inflow follows the cross-section.
    bool    conveyance = false;

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
