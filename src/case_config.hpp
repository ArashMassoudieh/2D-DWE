// =============================================================================
// case_config.hpp -- a case read from a text file instead of compiled in.
//
// Added for the Rock Creek digital-twin demonstration (2026-10): the twin writes
// a fresh config and hydrographs for every forecast that crosses the flood
// trigger, so the case cannot be a hard-coded function. Included from main.cpp
// after the Case struct.
//
// Format: one directive per line, '#' starts a comment. Paths run to the end of
// the line, so they may contain spaces.
//
//   name      rockcreek
//   dem       /path/dem_3m.tif            terrain GeoTIFF, used at its native cell size
//   manning   /path/manning_3m.tif        Manning n GeoTIFF on the same grid
//   crs       EPSG:26918                  CRS written to the output GeoTIFFs
//   t_end     172800                      s
//   t_out     900                         s, frame cadence
//   t_gauge   300                         s, gauge / hydrograph cadence
//   cfl 0.7 | theta 0.9 | n_default 0.06 | dt_max 5
//   inflow    north 551 558 /path/q.csv   side, [begin,end) cells, Q(t) CSV (t_s,Q_cms), split uniformly
//   inflow_conveyance north 500 600 /path/q.csv   the same, split by local conveyance (boundary.hpp)
//   h0        /path/depth.tif             initial depth (hot start), sampled onto the grid; default dry
//   outfall   south 0 -1 0.0042           side, [begin,end) (-1 = to the end), bed slope
//   point     Q_luzon 453 481 /path/q.csv name, column i, row j, Q(t) CSV
//   gauge     Joyce 345 498               name, column i, row j: depth and WSE every t_gauge
//   frames    depth speed                 GeoTIFF frames written every t_out (frames/)
//   vti       0                           1 = also write .vti frames
//   flood_h   0.15                        m, depth counted as "flooded" in the area series
//   time0     2014-04-30T00:00:00Z        label only: the wall-clock time of t = 0
// =============================================================================
#pragma once

#ifdef DWE_WITH_GDAL

#include <fstream>
#include <sstream>
#include <string>

static std::string cfg_trim(const std::string& s)
{
    const auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    const auto b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

static Side cfg_side(const std::string& w)
{
    if (w == "west")  return Side::West;
    if (w == "east")  return Side::East;
    if (w == "north") return Side::North;
    if (w == "south") return Side::South;
    std::printf("*** config: unknown side '%s'\n", w.c_str());
    std::exit(2);
}

static Series cfg_series(const std::string& path, const std::string& what)
{
    Series s = Series::from_csv(path);
    if (s.t.size() < 2) {
        std::printf("*** config: could not read %s hydrograph %s\n", what.c_str(), path.c_str());
        std::exit(2);     // a missing inflow must stop the run, not silently drop the water
    }
    return s;
}

static Case make_config(const std::string& path)
{
    Case c;
    c.is_config = true;
    c.name = "config";
    std::ifstream f(path);
    if (!f) { std::printf("*** config: cannot open %s\n", path.c_str()); std::exit(2); }

    std::string dem;
    struct InflowSpec { Side side; int b, e; std::string csv; bool conv; };
    std::string h0file;
    struct OutSpec    { Side side; int b, e; double slope; };
    std::vector<InflowSpec> ins;
    std::vector<OutSpec>    outs;
    struct PtSpec { std::string name; int i, j; std::string csv; };
    std::vector<PtSpec> pts;

    std::string line;
    while (std::getline(f, line)) {
        const auto hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        line = cfg_trim(line);
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string key; ss >> key;
        std::string rest;
        auto tail = [&]() { std::getline(ss, rest); return cfg_trim(rest); };

        if      (key == "name")      c.name = tail();
        else if (key == "dem")       dem = tail();
        else if (key == "manning")   c.manning = tail();
        else if (key == "crs")       c.crs = tail();
        else if (key == "t_end")     ss >> c.p.t_end;
        else if (key == "t_out")     ss >> c.p.t_out;
        else if (key == "t_gauge")   ss >> c.t_gauge;
        else if (key == "cfl")       ss >> c.p.cfl;
        else if (key == "theta")     ss >> c.p.theta;
        else if (key == "n_default") ss >> c.p.n_default;
        else if (key == "dt_max")    ss >> c.p.dt_max;
        else if (key == "flood_h")   ss >> c.flood_h;
        else if (key == "time0")     c.time0 = tail();
        else if (key == "vti")       { int v = 0; ss >> v; c.p.write_vti = (v != 0); }
        else if (key == "frames")    { std::string w; while (ss >> w) c.frames.push_back(w); }
        else if (key == "inflow" || key == "inflow_conveyance")
                                     { InflowSpec s; std::string w; ss >> w >> s.b >> s.e;
                                       s.side = cfg_side(w); s.csv = tail(); s.conv = (key != "inflow");
                                       ins.push_back(s); }
        else if (key == "h0")        h0file = tail();
        else if (key == "outfall")   { OutSpec s; std::string w; ss >> w >> s.b >> s.e >> s.slope;
                                       s.side = cfg_side(w); outs.push_back(s); }
        else if (key == "point")     { PtSpec s; ss >> s.name >> s.i >> s.j; s.csv = tail(); pts.push_back(s); }
        else if (key == "gauge")     { Case::Gauge gg; ss >> gg.name >> gg.i >> gg.j; c.gauges.push_back(gg); }
        else { std::printf("*** config: unknown directive '%s'\n", key.c_str()); std::exit(2); }
    }

    std::vector<char> valid;
    std::string err;
    if (dem.empty() || !load_terrain(dem, 0.0, Agg::Mean, c.g, c.z, valid, &err)) {
        std::printf("*** terrain: %s\n", dem.empty() ? "no 'dem' in the config" : err.c_str());
        std::exit(2);
    }
    std::printf("config        %s\n", path.c_str());
    std::printf("terrain       %s\n              %d x %d @ %.2f m\n", dem.c_str(), c.g.nx, c.g.ny, c.g.dx);
    c.h0.assign(c.g.ncell(), 0.0);
    if (!h0file.empty()) {
        if (!load_field(h0file, c.g, c.h0, 0.0, &err)) {
            std::printf("*** h0: %s\n", err.c_str());
            std::exit(2);
        }
        double v = 0.0; std::size_t nw = 0;
        for (auto& d : c.h0) { if (!(d > 0.0)) d = 0.0; v += d; nw += (d > 0.002); }
        std::printf("initial depth %s\n              %zu wet cells, %.0f m3\n", h0file.c_str(), nw, v * c.g.area());
    }

    for (const auto& s : ins) {
        Boundary b; b.side = s.side; b.type = BcType::Inflow; b.begin = s.b; b.end = s.e;
        b.Q = cfg_series(s.csv, "inflow"); b.name = "inflow"; b.conveyance = s.conv;
        std::printf("inflow        cells [%d,%d) %s  peak %.2f m3/s  volume %.0f m3  %s\n",
                    s.b, s.e, s.conv ? "conveyance-weighted" : "uniform",
                    b.Q.peak(), b.Q.volume(), s.csv.c_str());
        c.bcs.push_back(b);
    }
    for (const auto& s : outs) {
        Boundary b; b.side = s.side; b.type = BcType::NormalDepth; b.begin = s.b; b.end = s.e;
        b.slope = s.slope; b.name = "outfall";
        std::printf("outfall       normal depth, S = %.5f\n", s.slope);
        c.bcs.push_back(b);
    }
    for (const auto& s : pts) {
        if (!c.g.in_range(s.i, s.j)) { std::printf("*** point %s outside the grid\n", s.name.c_str()); std::exit(2); }
        SchemeInertial::PointSource q; q.i = s.i; q.j = s.j; q.name = s.name;
        q.Q = cfg_series(s.csv, s.name);
        std::printf("point         %-14s (%d,%d)  peak %.2f m3/s  volume %.0f m3\n",
                    s.name.c_str(), s.i, s.j, q.Q.peak(), q.Q.volume());
        c.points.push_back(q);
    }
    for (const auto& gg : c.gauges)
        if (!c.g.in_range(gg.i, gg.j)) { std::printf("*** gauge %s outside the grid\n", gg.name.c_str()); std::exit(2); }
    return c;
}

#endif // DWE_WITH_GDAL
