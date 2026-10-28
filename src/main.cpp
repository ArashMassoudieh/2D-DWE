// =============================================================================
// main.cpp -- driver.
//
// Two cases are wired in so the solver can be exercised the moment it builds,
// before any GDAL terrain is involved. Both are from the validation ladder in
// DESIGN 10, in the order that ladder specifies.
//
//   plane      steady uniform flow down a tilted plane. Must converge to the
//              Manning normal depth, which is known analytically. Tests the
//              friction term, the face rule and the boundaries.
//
//   pond       a closed depression. Water must fill to the spill elevation and
//              NO FURTHER. This is the one that exposes a leaky face rule
//              immediately, and it is directly analogous to the C6 divide.
//
//   c6         the real thing: Figueredo Wash and the Coyote Canyon community,
//              on 1 m LiDAR. Requires CONFIG+=with_gdal.
//
// Usage:  2D-DWE [plane|pond|c6] [--out DIR] [--tend S] [--dx M] [--agg mean|min|max]
// =============================================================================
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <filesystem>

#include "params.hpp"
#include "grid.hpp"
#include "fields.hpp"
#include "boundary.hpp"
#include "massbalance.hpp"
#include "scheme_inertial.hpp"
#include "vti_writer.hpp"
#ifdef DWE_WITH_GDAL
#include "io_gdal.hpp"
#endif

using namespace dwe;

// -----------------------------------------------------------------------------
// Manning normal depth for wide sheet flow: q = (1/n) h^(5/3) sqrt(S)
// -----------------------------------------------------------------------------
static double normal_depth(double q, double n, double S)
{
    return std::pow(q * n / std::sqrt(S), 3.0 / 5.0);
}

// -----------------------------------------------------------------------------
struct Case
{
    std::string          name;
    Grid                 g;
    Params               p;
    std::vector<Boundary> bcs;
    std::vector<double>  z;
    std::vector<double>  h0;
    double               expect = -1.0;   // analytical answer, if there is one
    std::string          expect_what;
    // Which design storm this case was built for. Kept explicitly rather than
    // inferred from t_end, because --tend overrides t_end on the command line
    // and would then silently select the wrong rainfall file.
    std::string          storm = "6hr";

    // --- config-file cases (case_config.hpp) ---------------------------------
    struct Gauge { std::string name; int i = 0, j = 0; };
    bool                 is_config = false;
    std::string          manning;              // Manning n GeoTIFF on the terrain grid
    std::string          crs = "EPSG:26912";   // CRS stamped on the output GeoTIFFs
    std::vector<SchemeInertial::PointSource> points;
    std::vector<Gauge>   gauges;
    double               t_gauge = 300.0;      // s
    std::vector<std::string> frames;           // "depth", "speed": GeoTIFF frames every t_out
    double               flood_h = 0.15;       // m, threshold for the flooded-area series
    std::string          time0;                // wall-clock label of t = 0
};

static Case make_plane()
{
    Case c;
    c.name = "plane";

    const int    nx = 200, ny = 20;
    const double dx = 5.0;
    const double S  = 0.005;          // 0.5 % down to the east
    const double n  = 0.030;
    const double Q  = 20.0;           // m^3/s total across the west face

    c.g = Grid(nx, ny, dx, 0.0, ny * dx);
    c.p.t_end       = 4000.0;
    c.p.t_out       = 500.0;
    c.p.n_default   = n;
    c.p.report_every = 2000;

    c.z.assign(c.g.ncell(), 0.0);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            c.z[c.g.k(i, j)] = (nx - 1 - i) * dx * S;   // high in the west

    c.h0.assign(c.g.ncell(), 0.0);

    Boundary in;  in.side = Side::West;  in.type = BcType::Inflow;
    in.Q = Series::constant(Q);          in.name = "inflow";
    Boundary out; out.side = Side::East; out.type = BcType::NormalDepth;
    out.slope = S;                       out.name = "outfall";
    c.bcs = {in, out};

    // q per unit width across the 20-cell (100 m) face
    const double q_unit = Q / (ny * dx);
    c.expect = normal_depth(q_unit, n, S);
    c.expect_what = "Manning normal depth";
    return c;
}

static Case make_pond()
{
    Case c;
    c.name = "pond";

    const int    nx = 120, ny = 120;
    const double dx = 2.0;

    c.g = Grid(nx, ny, dx, 0.0, ny * dx);
    c.p.t_end        = 2000.0;
    c.p.t_out        = 200.0;
    c.p.n_default    = 0.050;
    c.p.report_every = 2000;

    // A bowl inside a rim, with ONE low notch on the east side.
    //
    //   floor 6.0 m      notch 8.0 m      rim 10.0 m
    //
    // The bowl is PRE-FILLED to a water surface of 9.0 m -- above the notch,
    // below the rim -- and then left alone. No inflow. The only way out is the
    // notch, so the correct answer is: drain through the notch, settle at
    // 8.0 m, and never put a drop over the 10.0 m rim.
    //
    // This isolates the face-depth rule (DESIGN 4.1) with no boundary forcing
    // to hide behind, and it is the direct analogue of the C6 divide: a body of
    // water held by a rise, with one low point that governs.
    const double rim = 10.0, floor_z = 6.0, notch = 8.0, fill = 9.0;
    c.z.assign(c.g.ncell(), rim);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const double r = std::hypot(i - nx / 2.0, j - ny / 2.0);
            if (r < 40.0) c.z[c.g.k(i, j)] = floor_z;
            else if (r < 45.0)
                c.z[c.g.k(i, j)] = floor_z + (r - 40.0) / 5.0 * (rim - floor_z);
        }
    }
    // the notch: a 6-cell gap cut through the rim, due east of centre
    for (int j = ny / 2 - 3; j < ny / 2 + 3; ++j)
        for (int i = nx / 2 + 40; i < nx; ++i)
            c.z[c.g.k(i, j)] = notch;

    // pre-fill to WSE = 9.0 m wherever the bed is below it
    c.h0.assign(c.g.ncell(), 0.0);
    for (std::size_t k = 0; k < c.g.ncell(); ++k)
        if (c.z[k] < fill) c.h0[k] = fill - c.z[k];

    Boundary out; out.side = Side::East; out.type = BcType::NormalDepth;
    out.slope = 0.01; out.name = "notch outfall";
    c.bcs = {out};                     // no inflow: the pond drains and stops

    c.expect = notch;
    c.expect_what = "spill (notch) elevation -- must settle here, never over the rim";
    return c;
}

#ifdef DWE_WITH_GDAL
// -----------------------------------------------------------------------------
// The C6 case. Terrain, roughness and the domain all come from the Phase-1
// products in WS3_GIS/NM15-097/outputs_2D/ (see cases/c6/inputs_manifest.txt).
//
// Boundaries, as established by the Phase-1 geometry work:
//   inflow   the wash crosses the domain edge exactly ONCE, on the WEST side.
//            The reach has one inflow and no laterals, so the discharge is the
//            same wherever the face is placed: 32,628 cfs = 924 m3/s for the
//            governing 6-hour storm.
//   outflow  the FULL EAST edge. Sampled ground along it is 6,130-6,138 ft, the
//            low side of the domain, and both the wash and any spill onto the
//            community bench drain that way.
// -----------------------------------------------------------------------------
static const char* C6_DIR =
    "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097/outputs_2D";

// NOTE the path: outputs_2D/inputs, not output_2D. This site briefly had BOTH
// "output_2D" (model inputs) and "outputs_2D" (maps, animations) -- two directories in
// one site differing by a single letter, which is a trap for anyone looking for a file
// and found one immediately. Everything 2D now lives under outputs_2D: inputs/, maps/,
// animations/.
static const char* C4B_DIR =
    "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/AZ12-100/outputs_2D/inputs";

// Cluster C4 ANALYSIS A -- the MAIN WASH at the Indian Route 4 crossing, with the
// culvert COMPLETELY PLUGGED. Separate from C4B_DIR: two 2D models of one cluster
// on two different washes, one directory each, named for the analysis.
//
// The plugged case needs NO terrain edit. The culvert is a buried pipe and is
// invisible to LiDAR, so the bare-earth DEM already carries an unbroken road
// embankment -- 6060.4 ft over the road against 6041-6045 ft in the channel on
// either side. This model is therefore the one condition at this crossing that
// 2D-DWE can represent EXACTLY: no 1D structure is needed because there is no
// opening. Do not reuse it for the culvert-flowing case, which it cannot do.
static const char* C4A_DIR =
    "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/AZ12-100/outputs_2D_A/inputs";

static const char* NM141_DIR =
    "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-141Local/output_2D";
static const char* NM141_LC =
    "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-141Local/landcover/nlcd_2023_NM15-141Local.tif";

static const char* AZ140_SITE =
    "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/AZ12-140";

// Cluster C7 Analysis B, EXISTING condition. The local hillslope catchment above
// the Tohajiilee subdivisions (NM15-135).
static const char* C7B_SITE =
    "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-135";

static Case make_c6(double dx_target, Agg agg, const std::string& storm = "6hr",
                    const std::string& dem_override = std::string())
{
    Case c;
    c.name = "c6";

    // The engineer-supplied clip: 4372 x 2106 @ 1 m, EPSG:26912, no nodata.
    // The DEM defines the domain -- load_terrain derives the grid from it -- so
    // there is no separate boundary polygon to keep in step.
    // --dem lets a modified terrain (e.g. a carved berm) be substituted without
    // touching the case; everything else -- inflow face, roughness, storm -- is
    // unchanged, so the two runs are directly comparable.
    const std::string terr = dem_override.empty()
                           ? std::string(C6_DIR) + "/cliped_dem.tif"
                           : dem_override;
    const std::string land = std::string(C6_DIR) + "/landcover_2d.tif";
    const std::string csv  = std::string(C6_DIR) + "/manning_n_2d.csv";

    std::vector<char> valid;
    std::string err;
    if (!load_terrain(terr, dx_target, agg, c.g, c.z, valid, &err)) {
        std::printf("*** terrain: %s\n", err.c_str());
        std::exit(2);
    }
    std::printf("terrain       %s\n", terr.c_str());
    std::printf("              %d x %d @ %.2f m  (requested %.2f m)\n",
                c.g.nx, c.g.ny, c.g.dx, dx_target);
    std::printf("storm         %s 100-yr\n", storm == "24hr" ? "24-hour" : "6-hour");

    const bool s24  = (storm == "24hr");
    c.storm         = storm;
    c.p.n_default   = 0.050;
    // 24-hr HMS record runs 36 h: LOWER peak (29,759 vs 32,431 cfs) but far MORE
    // volume (6,600 vs 4,720 ac-ft). Which governs stage is the open question.
    c.p.t_end       = (s24 ? 36.0 : 9.0) * 3600.0;
    c.p.t_out       = 600.0;
    c.p.report_every = 500;
    c.p.dt_max      = 5.0;

    c.h0.assign(c.g.ncell(), 0.0);

    // --- inflow: 32,628 cfs on the west edge, on the wash ----------------
    //
    // The wash crosses the west edge at y = 3,959,857 (Phase 1). Convert that to
    // a row index on whatever grid the terrain gave us, and spread the discharge
    // over a face centred on it.
    // Find the wash on the west edge FROM THE TERRAIN rather than trusting a
    // hard-coded northing: the incised channel is the low ground, so take the
    // west column's minimum and keep every row within BAND_M of it. This
    // survives a change of DEM extent, which a fixed coordinate does not.
    const double BAND_M = 3.0;
    double zmin = 1e30;
    int    jmin = 0;
    for (int j = 0; j < c.g.ny; ++j) {
        const double zz = c.z[c.g.k(0, j)];
        if (zz < zmin) { zmin = zz; jmin = j; }
    }
    int jlo = jmin, jhi = jmin;
    while (jlo > 0            && c.z[c.g.k(0, jlo - 1)] <= zmin + BAND_M) --jlo;
    while (jhi < c.g.ny - 1   && c.z[c.g.k(0, jhi + 1)] <= zmin + BAND_M) ++jhi;

    Boundary in;
    in.side  = Side::West;
    in.type  = BcType::Inflow;
    in.begin = jlo;
    in.end   = jhi + 1;
    // TRANSIENT INFLOW -- the routed HMS hydrograph, not a steady peak.
    // Holding the peak steady forever is not conservative, it is simply wrong:
    // it delivers unbounded volume, so inundation keeps growing with run length
    // and no reported extent is a property of the storm. The hydrograph carries
    // 4,720 ac-ft and is done in 9 hours.
    //
    // The element is Reach_58, NOT Junction_28. Junction_28 holds the 32,628 cfs
    // headline peak but lies 1,190 m west of this domain; Reach_58 discharges
    // 115 m from the inflow face, attenuated to 32,431 cfs. See tools/hms_hydrograph.py.
    const std::string qcsv = std::string(C6_DIR)
        + (s24 ? "/inflow_24hr_100yr.csv" : "/inflow_6hr_100yr.csv");
    in.Q = Series::from_csv(qcsv);
    if (in.Q.t.size() < 2) {
        std::printf("*** could not read %s -- falling back to a STEADY peak.\n"
                    "*** Run tools/hms_hydrograph.py to regenerate it.\n", qcsv.c_str());
        in.Q = Series::constant(32431.0 * 0.028316847);
    }
    in.name  = "Figueredo Wash, 6-hr 100-yr";

    Boundary out;
    out.side  = Side::East;
    out.type  = BcType::NormalDepth;
    out.slope = 0.005;
    out.name  = "east edge outfall";

    c.bcs = {in, out};

    std::printf("inflow        HMS Reach_58, %s 100-yr TRANSIENT hydrograph\n"
                "              peak %.1f m3/s (%.0f cfs) at %.2f h, volume %.0f ac-ft,"
                " %.1f h long\n"
                "              over rows %d-%d on the west edge\n",
                s24 ? "24-hr" : "6-hr",
                in.Q.peak(), in.Q.peak() / 0.028316847, in.Q.t_peak() / 3600.0,
                in.Q.volume() / 1233.4818, in.Q.span() / 3600.0,
                in.begin, in.end);
    std::printf("              located from the terrain: west-edge low %.2f m at row %d,\n"
                "              face width %.0f m (%.0f ft) within %.1f m of it\n",
                zmin, jmin, (in.end - in.begin) * c.g.dx,
                (in.end - in.begin) * c.g.dx * 3.280839895, BAND_M);
    std::printf("outflow       normal depth, full east edge, S = %.4f\n", out.slope);

    // --- roughness -------------------------------------------------------
    // load_manning fills a vector on this grid; stash it in h0's sibling by
    // returning through the Case. Fields are built by the caller, so hand the
    // roughness back via a static and let main() apply it.
    return c;
}
#endif // DWE_WITH_GDAL

// -----------------------------------------------------------------------------

// =============================================================================
// C10 -- second Coyote Canyon cluster, with a culverted road crossing.
//
// GEOMETRY: the wash runs SOUTH to NORTH, entering the south edge near easting
// 711,290 at 1876.27 m and leaving the north edge near 711,285 at 1859.19 m --
// a 17 m drop over 1,633 m. Inflow is therefore on the SOUTH edge and the
// outfall on the NORTH, the opposite arrangement to C6.
//
// THE CULVERT is represented as a carved notch, not as a structure, and the DEM
// is supplied on the command line (--dem) precisely so the notch width can be
// swept. The crossing carries 3,870 cfs in the 100-year event while any
// realistic barrel passes a few percent of it, so the road overtops and the
// EMBANKMENT controls. Whether the culvert is represented at all should
// therefore change little -- the sweep is how that gets demonstrated rather
// than asserted.
//
// Do NOT feed this case outputs/dem_carved.tif from the GIS project. That DEM is
// carved full-width for D8 delineation, which needs hydrologic connectivity;
// used as 2D terrain it deletes the road as a control and can never overtop.
// =============================================================================
static Case make_nm141(double dx_target, Agg agg, const std::string& dem_override,
                       const std::string& storm)
{
    Case c;
    c.name = "nm141";

    const std::string terr = dem_override.empty()
                           ? std::string(NM141_DIR) + "/dem_2d_carved_2m.tif"
                           : dem_override;
    std::vector<char> valid;
    std::string err;
    if (!load_terrain(terr, dx_target, agg, c.g, c.z, valid, &err)) {
        std::printf("*** terrain: %s\n", err.c_str());
        std::exit(2);
    }
    std::printf("terrain       %s\n", terr.c_str());
    std::printf("              %d x %d @ %.2f m  (requested %.2f m)\n",
                c.g.nx, c.g.ny, c.g.dx, dx_target);

    c.p.n_default    = 0.050;
    // The 24-hr storm's HMS record runs 36 h, four times the 6-hr window, and it
    // carries MORE volume (1,238 ac-ft against 916) despite a LOWER peak.
    const bool s24 = (storm == "24hr");
    c.storm = storm;
    c.p.t_end        = (s24 ? 36.0 : 9.0) * 3600.0;
    c.p.t_out        = s24 ? 7200.0 : 1800.0;
    c.p.report_every = 500;
    c.p.dt_max       = 5.0;
    c.h0.assign(c.g.ncell(), 0.0);

    // --- THREE inflows: the wash and both tributaries -----------------------
    //
    // The HMS topology has three reaches crossing into this domain, and the
    // first version of this case used only one of them. Reach_122 from the west
    // is nearly as large as the main wash and was missing entirely, so the model
    // was carrying 51 % of the peak and 43 % of the volume:
    //
    //   Reach_147  SOUTH edge  3,870 cfs  394 ac-ft   the main wash
    //   Reach_122  WEST  edge  3,374 cfs  399 ac-ft   <-- was omitted
    //   Reach_99   EAST  edge    424 cfs   78 ac-ft   <-- was omitted
    //   local Subbasins 1,2,3,57  265 cfs  43 ac-ft   not represented (5 %)
    //   HMS Outlet total        7,619 cfs  915 ac-ft
    //
    // The three reaches carry 95 % of the outlet volume, so boundary inflows are
    // the right representation here and direct rainfall is a minor refinement.
    // That is the OPPOSITE of C6, where 94 % of the domain lay outside any
    // subwatershed and the local rainfall term was the only local contribution.
    //
    // Each face is located from the terrain NEAR THE KNOWN NHD CROSSING rather
    // than at the edge's global minimum. On the west edge the tributary is only
    // 0.4 m below its surroundings -- a global-minimum search finds it by luck,
    // not by construction, and would silently move if the DEM changed.
    struct InflowSpec { Side side; double target; const char* csv; const char* name; };
    std::printf("storm         %s 100-yr\n", s24 ? "24-hour" : "6-hour");
    const InflowSpec specs[3] = {
        { Side::South, 711299.0, s24 ? "/inflow_south_24hr_100yr.csv" : "/inflow_south_6hr_100yr.csv", "Reach_147, main wash (south)" },
        { Side::West,  3960293.0, s24 ? "/inflow_west_24hr_100yr.csv"  : "/inflow_west_6hr_100yr.csv",  "Reach_122, west tributary"    },
        { Side::East,  3959922.0, s24 ? "/inflow_east_24hr_100yr.csv"  : "/inflow_east_6hr_100yr.csv",  "Reach_99, east tributary"     },
    };
    const double BAND_M   = 3.0;
    const double SEARCH_M = 60.0;
    // MAX_FACE_M caps the inflow face. The main wash is incised 5.4 m at the
    // south edge and the band walk finds it cleanly, but NEITHER tributary is
    // incised where it crosses: the west edge is flat to within 0.8 m over
    // 400 m and the east edge simply rises monotonically. With nothing to stop
    // it the walk ran to 708 m and 906 m, which would spread 3,374 cfs into an
    // unphysical sheet a few inches deep.
    //
    // The cap is an ASSUMPTION, not a measurement, and it matters: at 60 m the
    // west tributary enters about 3 ft deep, at 200 m about 1.5 ft. It is set
    // here so it can be found and tested rather than buried.
    const double MAX_FACE_M = 60.0;

    std::vector<Boundary> bcs;
    for (const auto& sp : specs) {
        const bool horiz = (sp.side == Side::South || sp.side == Side::North);
        const int  n     = horiz ? c.g.nx : c.g.ny;
        // index of the target coordinate along this edge
        const int ctr = horiz ? int((sp.target - c.g.x0) / c.g.dx)
                              : int((c.g.y0 - sp.target) / c.g.dx);
        const int rad = std::max(1, int(SEARCH_M / c.g.dx));

        auto zedge = [&](int m) {
            switch (sp.side) {
            case Side::South: return c.z[c.g.k(m, c.g.ny - 1)];
            case Side::North: return c.z[c.g.k(m, 0)];
            case Side::West:  return c.z[c.g.k(0, m)];
            default:          return c.z[c.g.k(c.g.nx - 1, m)];
            }
        };
        double zmin = 1e30; int mmin = std::min(std::max(ctr, 0), n - 1);
        for (int m = std::max(0, ctr - rad); m <= std::min(n - 1, ctr + rad); ++m)
            if (zedge(m) < zmin) { zmin = zedge(m); mmin = m; }
        int lo = mmin, hi = mmin;
        const int maxhalf = std::max(1, int(0.5 * MAX_FACE_M / c.g.dx));
        while (lo > 0     && zedge(lo - 1) <= zmin + BAND_M && mmin - lo < maxhalf) --lo;
        while (hi < n - 1 && zedge(hi + 1) <= zmin + BAND_M && hi - mmin < maxhalf) ++hi;
        const bool capped = (mmin - lo >= maxhalf) || (hi - mmin >= maxhalf);

        Boundary in;
        in.side  = sp.side;
        in.type  = BcType::Inflow;
        in.begin = lo;
        in.end   = hi + 1;
        in.Q     = Series::from_csv(std::string(NM141_DIR) + sp.csv);
        in.name  = sp.name;
        if (in.Q.t.size() < 2) {
            std::printf("*** could not read %s%s -- that inflow is MISSING\n",
                        NM141_DIR, sp.csv);
            continue;
        }
        const char* sd = (sp.side == Side::South) ? "SOUTH"
                       : (sp.side == Side::West)  ? "WEST" : "EAST";
        std::printf("inflow        %-28s %-5s edge, peak %7.1f m3/s (%5.0f cfs)"
                    " at %.2f h\n"
                    "              %5.0f ac-ft, face %.0f m wide, edge low %.2f m\n",
                    sp.name, sd, in.Q.peak(), in.Q.peak() / 0.028316847,
                    in.Q.t_peak() / 3600.0, in.Q.volume() / 1233.4818,
                    (in.end - in.begin) * c.g.dx, zmin);
        if (capped)
            std::printf("              NOTE: face CAPPED at %.0f m -- no incised channel"
                        " on this edge, width is assumed\n", MAX_FACE_M);
        bcs.push_back(in);
    }

    Boundary out;
    out.side  = Side::North;
    out.type  = BcType::NormalDepth;
    out.slope = 0.010;
    out.name  = "north edge outfall";
    bcs.push_back(out);
    c.bcs = bcs;

    double qtot = 0.0, vtot = 0.0;
    for (const auto& b : c.bcs)
        if (b.type == BcType::Inflow) { qtot += b.Q.peak(); vtot += b.Q.volume(); }
    std::printf("              --- %zu inflows, peaks sum to %.0f cfs, volume %.0f ac-ft\n",
                c.bcs.size() - 1, qtot / 0.028316847, vtot / 1233.4818);
    std::printf("              (HMS Outlet: %s -- peaks are not concurrent)\n",
                s24 ? "6,782 cfs, 1,238 ac-ft" : "7,619 cfs, 915 ac-ft");
    std::printf("outflow       normal depth, full north edge, S = %.3f\n", out.slope);
    return c;
}


// =============================================================================
// NM15-032 (Ojo Amarillo, cluster C9) -- the NORTH-SIDE CHANNEL of the subdivision.
//
// GIS project:  WS3_GIS/NM15-032LocalE_North_Channel   (inputs: outputs_2D/inputs,
// built by outputs_2D/make_2d_inputs.py; HMS: WS3_HMS/NM15_032LocalE_North_Channel).
//
// GEOMETRY: a road-side ditch runs NORTH along the west side of N 3005, passes
// culvert 1 (under N 4080) and culvert 2 (under N 3005), then runs WEST along the
// excavated channel on the subdivision's north edge, through culvert 3 (Ojo
// Amarillo Rd) and culvert 4 (Elementary School Rd), and on WEST in a shallow
// swale to the watershed outlet. Domain E 734450-736100, N 4063700-4064750.
//
// TERRAIN: dem_2d_culverts.tif (default) carries a >= 1.5 m conveyance notch at
// each culvert; dem_2d_raw.tif is every culvert blocked. The culvert sizes are
// ASSUMED (2x48 in at 1-2, 2x24 in at 3, 2x18 in at 4: the largest standard CMP
// that fits under each road) -- NHA has not supplied them.
//
// INFLOWS: Junction_9 (the ditch, SOUTH edge, 4.1 km2) and Subbasin_4 (EAST
// edge). Subbasins 1+2+3 are a uniform volume-preserving rain source. A
// separate ~1.5 km2 drainage enters the south edge WEST of the subdivision; it is
// outside this watershed and NOT represented.
//
// ROUGHNESS: run with  --manning <inputs>/manning_existing_1m.tif  (channel full
// of shrubs, n 0.080) or manning_cleared_1m.tif (n 0.035). Without --manning it
// falls back to NLCD (landcover_2d.tif + manning_n_2d.csv), which does not see
// the channel vegetation at all.
// =============================================================================
static const char* NM032N_DIR =
    "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-032LocalE_North_Channel/outputs_2D/inputs";

static Case make_nm032n(double dx_target, Agg agg, const std::string& dem_override,
                        const std::string& storm)
{
    Case c;
    c.name = "nm032n";
    const std::string terr = dem_override.empty()
                           ? std::string(NM032N_DIR) + "/dem_2d_culverts.tif"
                           : dem_override;
    std::vector<char> valid;
    std::string err;
    if (!load_terrain(terr, dx_target, agg, c.g, c.z, valid, &err)) {
        std::printf("*** terrain: %s\n", err.c_str());
        std::exit(2);
    }
    std::printf("terrain       %s\n", terr.c_str());
    std::printf("              %d x %d @ %.2f m  (requested %.2f m)\n",
                c.g.nx, c.g.ny, c.g.dx, dx_target);

    c.p.n_default    = 0.050;
    const bool s24 = (storm == "24hr");
    c.storm = storm;
    c.p.t_end        = (s24 ? 36.0 : 9.0) * 3600.0;     // the HMS record length
    c.p.t_out        = s24 ? 7200.0 : 1800.0;
    c.p.report_every = 500;
    c.p.dt_max       = 5.0;
    c.h0.assign(c.g.ncell(), 0.0);

    struct InflowSpec { Side side; double target; const char* csv; const char* name; };
    std::printf("storm         %s 100-yr\n", s24 ? "24-hour" : "6-hour");
    const InflowSpec specs[2] = {
        { Side::South, 735908.0,  s24 ? "/inflow_south_24hr_100yr.csv" : "/inflow_south_6hr_100yr.csv",
          "Junction_9, road-side ditch (south)" },
        { Side::East,  4063808.0, s24 ? "/inflow_east_24hr_100yr.csv"  : "/inflow_east_6hr_100yr.csv",
          "Subbasin_4, NE subbasin (east)" },
    };
    // The ditch is a narrow V beside a 1-2 m road embankment: a 3 m band (NM15-141)
    // would run the face up onto the road. 1 m and a 30 m cap keep it in the ditch.
    const double BAND_M = 1.0, SEARCH_M = 30.0, MAX_FACE_M = 30.0;

    std::vector<Boundary> bcs;
    for (const auto& sp : specs) {
        const bool horiz = (sp.side == Side::South || sp.side == Side::North);
        const int  n     = horiz ? c.g.nx : c.g.ny;
        const int ctr = horiz ? int((sp.target - c.g.x0) / c.g.dx)
                              : int((c.g.y0 - sp.target) / c.g.dx);
        const int rad = std::max(1, int(SEARCH_M / c.g.dx));
        auto zedge = [&](int m) {
            switch (sp.side) {
            case Side::South: return c.z[c.g.k(m, c.g.ny - 1)];
            case Side::North: return c.z[c.g.k(m, 0)];
            case Side::West:  return c.z[c.g.k(0, m)];
            default:          return c.z[c.g.k(c.g.nx - 1, m)];
            }
        };
        double zmin = 1e30; int mmin = std::min(std::max(ctr, 0), n - 1);
        for (int m = std::max(0, ctr - rad); m <= std::min(n - 1, ctr + rad); ++m)
            if (zedge(m) < zmin) { zmin = zedge(m); mmin = m; }
        int lo = mmin, hi = mmin;
        const int maxhalf = std::max(1, int(0.5 * MAX_FACE_M / c.g.dx));
        while (lo > 0     && zedge(lo - 1) <= zmin + BAND_M && mmin - lo < maxhalf) --lo;
        while (hi < n - 1 && zedge(hi + 1) <= zmin + BAND_M && hi - mmin < maxhalf) ++hi;
        const bool capped = (mmin - lo >= maxhalf) || (hi - mmin >= maxhalf);

        Boundary in;
        in.side  = sp.side;
        in.type  = BcType::Inflow;
        in.begin = lo;
        in.end   = hi + 1;
        in.Q     = Series::from_csv(std::string(NM032N_DIR) + sp.csv);
        in.name  = sp.name;
        if (in.Q.t.size() < 2) {
            // Both inflows are required: refuse rather than run a model that is
            // silently missing one (the NM15-141 case carries on, and prints a
            // "MISSING" line nobody reads).
            std::printf("*** could not read %s%s -- refusing to run without it\n",
                        NM032N_DIR, sp.csv);
            std::exit(2);
        }
        std::printf("inflow        %-36s %-5s edge, peak %6.2f m3/s (%4.0f cfs) at %.2f h\n"
                    "              %5.1f ac-ft, face %.0f m wide, edge low %.2f m\n",
                    sp.name, sp.side == Side::South ? "SOUTH" : "EAST",
                    in.Q.peak(), in.Q.peak() / 0.028316847, in.Q.t_peak() / 3600.0,
                    in.Q.volume() / 1233.4818, (in.end - in.begin) * c.g.dx, zmin);
        if (capped)
            std::printf("              NOTE: face CAPPED at %.0f m\n", MAX_FACE_M);
        bcs.push_back(in);
    }
    // Outfalls: the swale leaves through the WEST edge (bed slope ~0.006); a small
    // area north of the channel drains off the NORTH edge.
    Boundary w;  w.side = Side::West;  w.type = BcType::NormalDepth; w.slope = 0.006;
    w.name = "west edge outfall (watershed outlet)";  bcs.push_back(w);
    Boundary nn; nn.side = Side::North; nn.type = BcType::NormalDepth; nn.slope = 0.010;
    nn.name = "north edge outfall";                  bcs.push_back(nn);
    c.bcs = bcs;
    std::printf("outflow       normal depth: west edge S = 0.006, north edge S = 0.010\n");
    std::printf("              HMS Outlet peak %s\n", s24 ? "204 cfs" : "207 cfs");
    return c;
}


// -----------------------------------------------------------------------------
// C4 ANALYSIS B -- the South Tributary Wash at Cottonwood / Tselani, AZ
// (NHA projects AZ12-100 and AZ12-116, cluster C4).
//
// WHY THIS CASE EXISTS.  The 1-D HEC-RAS model of this wash cannot answer the
// question the study is actually asking.  The flow is shallow and unconfined --
// 0.7 to 1.7 ft deep over top widths of 114 to 230 ft -- spreading across a
// nearly level bench, and the exposure of the AZ12-116 homes turns entirely on
// whether an unnamed road embankment between the wash and the subdivision
// overtops, and where.  A 1-D model crosses that embankment as ordinary ground
// at every section: it can neither hold water against it nor route water ALONG
// it to the low point and spill there.  The 1-D floodplain is a bounded lower
// estimate and the report says so.
//
// The embankment is resolved cell-for-cell at 1 m -- 5.6 to 6.2 ft above the
// wash bed at every station checked -- so no breakline and no mesh is involved.
// That is the whole reason a structured-grid solver suits this site.
//
// ONE INFLOW, ON THE SOUTH EDGE.  The wash crosses the south edge in a flat
// swale, X 599970-600040 at 6087 ft, shoulders rising to 6094 west and 6127
// east.  There is NO incised channel anywhere on this reach -- 2.4 ft of relief
// over 105 m -- so the face is found by the same edge-low walk nm141 uses and
// then CAPPED, exactly as its two unincised tributaries are.  70 m is the swale
// and it matches the 114-230 ft top widths the 1-D model computes; the cap is
// an assumption, and it is here where it can be found and tested.
//
// The hydrograph is Junction_2 plus the RS 706 lateral, 230 cfs.  The lateral
// joins 100 m inside the domain and its catchment lies west, so a west-edge
// inflow looks natural -- but the MAIN WASH crosses that edge 30 m away and
// 21 ft lower, and an edge-low search would put 32 cfs into the wrong drainage.
// Folding it into the south edge moves it about 110 m upstream on a 320 m
// reach, which is the smaller and the more conservative error.
//
// TWO OUTFALLS.  The domain drains NORTH into the Main Wash at 1.14 %; the east
// edge is flat to within 0.001 and only takes water that gets that far across
// the bench.  Both are normal-depth so nothing ponds against a domain edge.
//
// The east edge sits at X 600400, 13 m past the last AZ12-100 structure. An
// earlier version ran to X 601000 on the reasoning that the bench falls east
// toward the Main Wash lowland -- true on a transect at Y 3992420, but the
// domain-wide gradient settles it: north 1.14 %, east -0.0008. The water leaves
// NORTH into the Main Wash, and the 600 m of eastern ground was carrying no
// flow and 55 % of the cells. Trimmed. The east edge stays normal-depth, so
// anything that does reach it still leaves rather than ponding on a wall.
//
// WHAT THIS DOMAIN DOES NOT CONTAIN, DELIBERATELY.  The Indian Route 4 culvert
// is Analysis A's structure on the MAIN WASH; the north edge stops 55 m short of
// it.  A stretch of the Main Wash channel does fall inside the rectangle and
// carries no flow here, so the tributary drains into an empty channel that in
// reality holds the Analysis A flood.  The Main Wash 24-hr surface stands 0.63 ft
// below this tributary's outlet invert with the culvert clean and 0.03 ft below
// it as found silted -- free-discharging in both, but only just.  The exposure
// this model exists to resolve is ponding behind the road some 200 m UPSTREAM of
// that outlet, driven by the tributary's own 230 cfs, so the outlet condition
// does not decide it.  That is the argument for the simplification.  It is NOT a
// claim that the two washes are independent, and it is the first thing to revisit
// if the answer comes out close.
// -----------------------------------------------------------------------------
static Case make_c4b(double dx_target, Agg agg, const std::string& dem_override,
                     const std::string& storm)
{
    Case c;
    c.name = "c4b";

    const std::string terr = dem_override.empty()
                           ? std::string(C4B_DIR) + "/dem_2d_raw.tif"
                           : dem_override;
    std::vector<char> valid;
    std::string err;
    if (!load_terrain(terr, dx_target, agg, c.g, c.z, valid, &err)) {
        std::printf("*** terrain: %s\n", err.c_str());
        std::exit(2);
    }
    std::printf("terrain       %s\n", terr.c_str());
    std::printf("              %d x %d @ %.2f m  (requested %.2f m)\n",
                c.g.nx, c.g.ny, c.g.dx, dx_target);
    std::printf("              RAW LiDAR -- NOT outputs/dem_carved.tif, which is cut\n"
                "              full width through the road for D8 delineation and would\n"
                "              delete the very control this model exists to represent.\n");

    c.p.n_default    = 0.050;                 // Shrub/Scrub, 85% of the domain
    const bool s24   = (storm == "24hr");
    c.storm          = storm;
    c.p.t_end        = (s24 ? 36.0 : 9.0) * 3600.0;   // the HMS record length
    c.p.t_out        = s24 ? 3600.0 : 900.0;
    c.p.report_every = 500;
    c.p.dt_max       = 5.0;
    c.h0.assign(c.g.ncell(), 0.0);

    std::printf("storm         %s 100-yr\n", s24 ? "24-hour" : "6-hour");

    // --- the single south-edge inflow ---------------------------------------
    const double TARGET_X   = 600005.0;   // centre of the swale on the south edge
    const double SEARCH_M   = 60.0;
    const double MAX_FACE_M = 80.0;       // the swale is ~70 m; see the note above
    const std::string csv = std::string(C4B_DIR)
                          + (s24 ? "/inflow_south_24hr_100yr.csv"
                                 : "/inflow_south_6hr_100yr.csv");

    const int n   = c.g.nx;
    const int ctr = int((TARGET_X - c.g.x0) / c.g.dx);
    const int rad = std::max(1, int(SEARCH_M / c.g.dx));
    auto zedge = [&](int m) { return c.z[c.g.k(m, c.g.ny - 1)]; };

    double zmin = 1e30; int mmin = std::min(std::max(ctr, 0), n - 1);
    for (int m = std::max(0, ctr - rad); m <= std::min(n - 1, ctr + rad); ++m)
        if (zedge(m) < zmin) { zmin = zedge(m); mmin = m; }

    const double BAND_M  = 3.0;
    const int    maxhalf = std::max(1, int(0.5 * MAX_FACE_M / c.g.dx));
    int lo = mmin, hi = mmin;
    while (lo > 0            && zedge(lo - 1) < zmin + BAND_M && (mmin - lo) < maxhalf) --lo;
    while (hi < n - 1        && zedge(hi + 1) < zmin + BAND_M && (hi - mmin) < maxhalf) ++hi;
    const bool capped = ((mmin - lo) >= maxhalf) || ((hi - mmin) >= maxhalf);

    Boundary in;
    in.side  = Side::South;
    in.type  = BcType::Inflow;
    in.begin = lo;
    in.end   = hi + 1;
    in.name  = "Junction_2 + RS 706 lateral, south edge";
    in.Q     = Series::from_csv(csv);
    if (in.Q.t.size() < 2) {
        std::printf("*** inflow: could not read %s\n", csv.c_str());
        std::exit(2);
    }
    std::printf("inflow        %s\n", in.name.c_str());
    std::printf("              peak %6.1f m3/s (%6.0f cfs) at %.2f h, %.0f ac-ft\n",
                in.Q.peak(), in.Q.peak() / 0.028316847, in.Q.t_peak() / 3600.0,
                in.Q.volume() / 1233.4818);
    std::printf("              face %.0f m wide, edge low %.2f m%s\n",
                (in.end - in.begin) * c.g.dx, zmin,
                capped ? "  -- CAPPED, width is assumed" : "");

    std::vector<Boundary> bcs{in};

    Boundary outN; outN.side = Side::North; outN.type = BcType::NormalDepth;
    outN.slope = 0.011; outN.name = "north edge outfall (to the Main Wash)";
    bcs.push_back(outN);
    Boundary outE; outE.side = Side::East;  outE.type = BcType::NormalDepth;
    outE.slope = 0.010; outE.name = "east edge outfall (across the bench)";
    bcs.push_back(outE);
    c.bcs = bcs;

    std::printf("outflow       normal depth on the NORTH edge (S = %.3f) and the\n"
                "              EAST edge (S = %.3f)\n", outN.slope, outE.slope);
    std::printf("              (HMS Outlet: %s)\n",
                s24 ? "233 cfs, 25.3 ac-ft" : "232 cfs, 12.5 ac-ft");
    return c;
}


// -----------------------------------------------------------------------------
// C4 ANALYSIS A -- Main Wash at Indian Route 4, CULVERT COMPLETELY PLUGGED.
//
// WHAT THIS ANSWERS THAT THE 1D MODEL COULD NOT
//   1. The blockage condition itself. The field form records both barrels 95%
//      blocked. HEC-RAS still passes flow through a culvert opening, so the 1D
//      study bounds that condition with a deck held level at 6059.95 ft. Here
//      the embankment simply is solid, straight from the LiDAR.
//   2. The lateral spill. The 1D model reports nine sections carrying water at
//      or above a cut-line end, all upstream of the crossing, and has no
//      mechanism to shed that water sideways. The domain runs west to X 599450,
//      past the NW end of the road at (599850, 3992945), so the solver decides
//      for itself whether the pool escapes around the embankment.
//
// INFLOW.  One face on the SOUTH edge, in the incised notch at X 599628 (invert
// 6069 ft, banks 6097 ft, 6 m wide at invert+6 ft). All three HMS inflows are
// summed into it -- see C4A-2d-inputs.py for why that is legitimate here.
//
// NO RAINFALL TERM.  1.35 x 1.22 km inside a 12.69 km2 watershed; the boundary
// hydrograph already carries the routed runoff and a uniform excess on top of it
// would double count. Deliberately NOT added to the rainfall block below.
// -----------------------------------------------------------------------------
static Case make_c4a(double dx_target, Agg agg, const std::string& dem_override,
                     const std::string& storm)
{
    Case c;
    c.name = "c4a";

    const std::string terr = dem_override.empty()
                           ? std::string(C4A_DIR) + "/dem_2d_raw.tif"
                           : dem_override;
    std::vector<char> valid;
    std::string err;
    if (!load_terrain(terr, dx_target, agg, c.g, c.z, valid, &err)) {
        std::printf("*** terrain: %s\n", err.c_str());
        std::exit(2);
    }
    std::printf("terrain       %s\n", terr.c_str());
    std::printf("              %d x %d @ %.2f m  (requested %.2f m)\n",
                c.g.nx, c.g.ny, c.g.dx, dx_target);
    std::printf("              RAW LiDAR. The culvert is a buried pipe, invisible to\n"
                "              LiDAR, so this terrain IS the completely plugged case --\n"
                "              no notch carved, no deck assumed, no 1D structure needed.\n");

    c.p.n_default    = 0.050;                 // Shrub/Scrub, 92% of the domain
    const bool s24   = (storm == "24hr");
    c.storm          = storm;
    c.p.t_end        = (s24 ? 36.0 : 9.0) * 3600.0;   // the HMS record length
    c.p.t_out        = s24 ? 3600.0 : 900.0;
    c.p.report_every = 500;
    c.p.dt_max       = 5.0;
    c.h0.assign(c.g.ncell(), 0.0);

    std::printf("storm         %s 100-yr\n", s24 ? "24-hour" : "6-hour");

    // --- the single south-edge inflow, in the incised notch -----------------
    const double TARGET_X   = 599628.0;
    const double SEARCH_M   = 80.0;
    const double BAND_M     = 4.0;    // the notch is steep; 4 m picks up the bed
    const double MAX_FACE_M = 40.0;
    const std::string csv = std::string(C4A_DIR)
                          + (s24 ? "/inflow_mainwash_24hr_100yr.csv"
                                 : "/inflow_mainwash_6hr_100yr.csv");

    const int n   = c.g.nx;
    const int ctr = int((TARGET_X - c.g.x0) / c.g.dx);
    const int rad = std::max(1, int(SEARCH_M / c.g.dx));
    // Side::South is the LOW-Y edge, which in a north-up raster is the LAST row,
    // not row 0. Getting this wrong puts the inflow on the north edge -- where the
    // ground is 15 ft lower -- and the model quietly runs the flood in backwards.
    auto zedge = [&](int m) { return c.z[c.g.k(m, c.g.ny - 1)]; };

    double zmin = 1e30; int mmin = std::min(std::max(ctr, 0), n - 1);
    for (int m = std::max(0, ctr - rad); m <= std::min(n - 1, ctr + rad); ++m)
        if (zedge(m) < zmin) { zmin = zedge(m); mmin = m; }

    const int maxhalf = std::max(1, int(0.5 * MAX_FACE_M / c.g.dx));
    int lo = mmin, hi = mmin;
    while (lo > 0     && zedge(lo - 1) < zmin + BAND_M && (mmin - lo) < maxhalf) --lo;
    while (hi < n - 1 && zedge(hi + 1) < zmin + BAND_M && (hi - mmin) < maxhalf) ++hi;
    const bool capped = ((mmin - lo) >= maxhalf) || ((hi - mmin) >= maxhalf);

    Boundary in;
    in.side  = Side::South;
    in.type  = BcType::Inflow;
    in.begin = lo;
    in.end   = hi + 1;
    in.name  = "Main Wash: Reach_77 + Reach_82 laterals + local subbasins";
    in.Q     = Series::from_csv(csv);
    if (in.Q.t.size() < 2) {
        std::printf("*** inflow: could not read %s\n", csv.c_str());
        std::exit(2);
    }
    std::printf("inflow        %s\n", in.name.c_str());
    std::printf("              peak %6.1f m3/s (%6.0f cfs) at %.2f h, %.0f ac-ft\n",
                in.Q.peak(), in.Q.peak() / 0.028316847, in.Q.t_peak() / 3600.0,
                in.Q.volume() / 1233.4818);
    std::printf("              face %.0f m wide in the notch, edge low %.2f m%s\n",
                (in.end - in.begin) * c.g.dx, zmin,
                capped ? "  -- CAPPED, width is assumed" : "");

    std::vector<Boundary> bcs{in};

    // Outfalls. The ground falls consistently to the NE; the north and east
    // edges are 99% and 86% below the 1D pool elevation. The WEST edge is
    // mostly high ground, but its northern 15% lies BEYOND the NW end of the
    // road -- water reaching it has already escaped the impoundment, so holding
    // it in would be the model inventing a control the terrain does not have.
    Boundary outN; outN.side = Side::North; outN.type = BcType::NormalDepth;
    outN.slope = 0.010; outN.name = "north edge outfall (downstream, ground falls NE)";
    bcs.push_back(outN);
    Boundary outE; outE.side = Side::East;  outE.type = BcType::NormalDepth;
    outE.slope = 0.010; outE.name = "east edge outfall";
    bcs.push_back(outE);
    Boundary outW; outW.side = Side::West;  outW.type = BcType::NormalDepth;
    outW.slope = 0.008; outW.name = "west edge outfall (past the NW end of the road)";
    bcs.push_back(outW);
    c.bcs = bcs;

    std::printf("outflow       normal depth on NORTH (S=%.3f), EAST (S=%.3f) and\n"
                "              WEST (S=%.3f) edges; SOUTH closed by terrain\n",
                outN.slope, outE.slope, outW.slope);
    std::printf("              (1D comparison: HW 6060.56 ft, 17.0 ac, 9 structures)\n");
    return c;
}


// -----------------------------------------------------------------------------
// AZ12-140 -- Tuba City tributary, reach 89.  An INDEPENDENT CHECK on the
// HEC-RAS 2-D model built by hecras-qt, not a study in its own right.
//
// It is deliberately set up to match that model rather than to improve on it:
// the same DEM, the same NLCD roughness, the same 24-hr 100-yr hydrograph on
// the same edge, and the same free-flow treatment of the other three edges.
// Where the two solvers disagree, the disagreement is the finding.
//
// WHAT THE COMPARISON IS WORTH.  The RAS run this checks closed its volume
// accounting to only 69 %, with the east boundary discharging 125 % of the
// inflow.  A second solver on the same terrain either reproduces that or does
// not, and 2D-DWE reports its own mass balance to 1e-14, so it can tell the two
// apart: a real hydraulic result that both models find, or an artifact of one.
//
// NO RAINFALL TERM.  Unlike c6 and nm141 this case applies none.  The domain is
// a 1.9 x 2.0 km reach inside a much larger watershed and the boundary
// hydrograph already carries the routed upstream runoff; a uniform excess on
// top of it would double count.  The local contribution is a real but small
// omission, and it biases this check LOW against the RAS model, which makes the
// same omission.
// -----------------------------------------------------------------------------
static Case make_az140(double dx_target, Agg agg, const std::string& dem_override,
                       const std::string& storm)
{
    Case c;
    c.name = "az140";

    const std::string terr = dem_override.empty()
                           ? std::string(AZ140_SITE) + "/demlr/cliped_utm_for_ras2d.tif"
                           : dem_override;
    std::vector<char> valid;
    std::string err;
    if (!load_terrain(terr, dx_target, agg, c.g, c.z, valid, &err)) {
        std::printf("*** terrain: %s\n", err.c_str());
        std::exit(2);
    }
    std::printf("terrain       %s\n", terr.c_str());
    std::printf("              %d x %d @ %.2f m  (requested %.2f m)\n",
                c.g.nx, c.g.ny, c.g.dx, dx_target);

    const bool s24 = (storm != "6hr");
    c.storm = s24 ? "24hr" : "6hr";
    c.p.n_default    = 0.045;
    c.p.t_end        = 72.0 * 3600.0;      // the RAS plan's duration, to compare
    c.p.t_out        = 3600.0;             // hourly, matching the RAS output steps
    c.p.report_every = 500;
    c.p.dt_max       = 5.0;
    c.h0.assign(c.g.ncell(), 0.0);

    // --- inflow on the NORTH edge ------------------------------------------
    // The face is anchored on the CHANNEL CROSSING, not on the edge's global
    // minimum.  The first version of this case used the minimum, and on this
    // domain the north edge runs downhill to its east end, so the minimum IS the
    // corner: the inflow landed beside the east outlet and the water left 40 m
    // from where it entered, wetting 0.4 acres in 72 hours.  HEC-RAS made the
    // identical mistake -- 99 % of its inflow went into the last two faces of
    // 99 -- which is where its 69-153 % volume error came from.
    //
    // The crossing is taken from the watershed flow-accumulation raster
    // (outputs/flow_acc.tif): the north edge's channel carries 549,214 cells of
    // upstream area at easting 479966.5, 63 % along the edge.  This is the same
    // rule nm141 uses, and for the same reason its comments give.
    const double INFLOW_X   = 479966.5;   // UTM 12N easting of the north crossing
    const double SEARCH_M   = 60.0;       // refine to the local low within this
    const double BAND_M     = 3.0;
    const double MAX_FACE_M = 400.0;
    std::vector<Boundary> bcs;
    {
        auto zedge = [&](int i) { return c.z[c.g.k(i, 0)]; };   // north row
        const int ctr = int((INFLOW_X - c.g.x0) / c.g.dx);
        const int rad = std::max(1, int(SEARCH_M / c.g.dx));
        double zmin = 1e30; int imin = std::min(std::max(ctr, 0), c.g.nx - 1);
        for (int i = std::max(0, ctr - rad); i <= std::min(c.g.nx - 1, ctr + rad); ++i)
            if (zedge(i) < zmin) { zmin = zedge(i); imin = i; }
        int lo = imin, hi = imin;
        const int maxhalf = std::max(1, int(0.5 * MAX_FACE_M / c.g.dx));
        while (lo > 0          && zedge(lo - 1) <= zmin + BAND_M && imin - lo < maxhalf) --lo;
        while (hi < c.g.nx - 1 && zedge(hi + 1) <= zmin + BAND_M && hi - imin < maxhalf) ++hi;

        Boundary in;
        in.side  = Side::North;
        in.type  = BcType::Inflow;
        in.begin = lo;
        in.end   = hi + 1;
        in.name  = "24-hr 100-yr hydrograph (north)";
        in.Q     = Series::from_csv(std::string(AZ140_SITE) +
                       (s24 ? "/outputs_2D/dwe_inflow_24hr_100yr.csv"
                            : "/outputs_2D/dwe_inflow_6hr_100yr.csv"));
        if (in.Q.t.size() < 2) {
            std::printf("*** could not read the inflow CSV -- run the converter first\n");
            std::exit(2);
        }
        std::printf("storm         %s 100-yr\n", s24 ? "24-hour" : "6-hour");
        std::printf("inflow        %-30s NORTH edge, peak %6.2f m3/s (%5.0f cfs) at %.2f h\n"
                    "              %5.0f ac-ft, face %.0f m wide, channel low %.2f m\n",
                    in.name.c_str(), in.Q.peak(), in.Q.peak() / 0.028316847,
                    in.Q.t_peak() / 3600.0, in.Q.volume() / 1233.4818,
                    (in.end - in.begin) * c.g.dx, zmin);
        bcs.push_back(in);
    }

    // --- free flow out of the other three edges -----------------------------
    // The RAS model puts a Normal Depth BC on south, east and west, each with a
    // slope measured from the terrain draining through it.  The same three
    // slopes are used here so the two models are given the same outlet.
    struct Out { Side side; double slope; const char* name; };
    const Out outs[3] = {
        { Side::South, 0.02307, "south, normal depth" },
        { Side::East,  0.04512, "east, normal depth"  },
        { Side::West,  0.02307, "west, normal depth"  },
    };
    for (const auto& o : outs) {
        Boundary b;
        b.side  = o.side;
        b.type  = BcType::NormalDepth;
        b.slope = o.slope;
        b.name  = o.name;
        bcs.push_back(b);
        std::printf("outflow       %-30s S = %.5f\n", o.name, o.slope);
    }
    c.bcs = bcs;
    return c;
}


// -----------------------------------------------------------------------------
// C7B -- the local hillslope catchment above the Tohajiilee subdivisions,
// EXISTING condition (no interceptor channel).
//
// WHY RAIN ON GRID AND NO INFLOW BOUNDARY AT ALL.  The mechanism here is not a
// channel arriving at a boundary: it is 54 acres of hillside delivering sheet
// and shallow concentrated flow across 97 separate cells along 2,658 ft of the
// subdivisions' eastern frontage.  Imposing that as boundary hydrographs would
// mean discretising a line into points, and measurement shows that cannot be
// done: pour points on the two largest draws capture 27 % of the contributing
// area, and even 50 of them reach only 96 %.
//
// So the WHOLE 274-acre catchment is put inside the domain and forced with
// excess precipitation.  Every one of those 97 entry flows is then generated by
// the model rather than imposed on it, and there is nothing to discretise.
//
// EXCESS, NOT GROSS.  The solver has no infiltration model.  The series is the
// HMS PRECIP-EXCESS, area-weighted over all 28 subbasins -- the same SCS
// curve-number loss, the same CNs and the same hyetograph as the rest of the
// study.  See PythonScripts/c7b/C7B-rain-excess.py.
//
// NO DOUBLE COUNTING is possible here, unlike c6 and nm141: there is no
// boundary hydrograph to double count against.
//
// OUTFLOW ON THE EAST EDGE ONLY.  Measured, not assumed: tracing D8 downstream
// from the catchment outlet leaves the domain at its south-east corner, and of
// the four edges only the east has ground falling outward over its low quartile
// (S = 0.00572; south is -0.00109, i.e. flat).  A normal-depth BC needs a
// positive slope, so the south edge is left closed and the run is CHECKED
// afterwards for water piling against it.
// -----------------------------------------------------------------------------
static Case make_c7b(double dx_target, Agg agg, const std::string& dem_override,
                     const std::string& storm)
{
    Case c;
    c.name = "c7b";

    const std::string terr = dem_override.empty()
                           ? std::string(C7B_SITE) + "/outputs_2D/terrain_2d.tif"
                           : dem_override;
    std::vector<char> valid;
    std::string err;
    if (!load_terrain(terr, dx_target, agg, c.g, c.z, valid, &err)) {
        std::printf("*** terrain: %s\n", err.c_str());
        std::exit(2);
    }
    std::printf("terrain       %s\n", terr.c_str());
    std::printf("              %d x %d @ %.2f m  (requested %.2f m)\n",
                c.g.nx, c.g.ny, c.g.dx, dx_target);

    const bool s24 = (storm != "6hr");
    c.storm = s24 ? "24hr" : "6hr";
    c.p.n_default    = 0.045;
    // The excess series spans 36 h (24-hr storm) and 9 h (6-hr). Run past the
    // end of the rain so the recession is captured, as at the other sites.
    c.p.t_end        = (s24 ? 48.0 : 18.0) * 3600.0;
    c.p.t_out        = 900.0;
    c.p.report_every = 500;
    c.p.dt_max       = 5.0;
    c.h0.assign(c.g.ncell(), 0.0);

    std::printf("storm         %s 100-yr, EXISTING condition (no interceptor)\n",
                s24 ? "24-hour" : "6-hour");
    std::printf("inflow        NONE -- direct rainfall on the domain\n");

    // EVERY EDGE CARRIES A BC. The first version left the other three closed on
    // the reasoning that only the east falls outward. A closed edge is a
    // REFLECTIVE WALL, and a real terrain hollow straddles the west edge at
    // 318,000 / 3,882,600 -- 277 m OUTSIDE the catchment, in the margin. Dammed
    // there, it filled to 3.9 ft (the domain maximum) and then sloshed, throwing
    // a 164 ft/s velocity. Neither figure is physical and both were in the
    // domain-wide statistics. AZ12-140 met the same thing in HEC-RAS.
    //
    // The east slope is measured from the terrain over the low quartile of that
    // edge. The other three have ground RISING outward, so no slope can be
    // measured from them; they get a nominal free drain, which is closer to the
    // truth than a wall and cannot affect the catchment interior 277 m away.
    struct Out { Side side; double slope; const char* name; };
    const Out outs[4] = {
        { Side::East,  0.00572, "east, normal depth (measured)"  },
        { Side::South, 0.00500, "south, normal depth (nominal)"  },
        { Side::West,  0.00500, "west, normal depth (nominal)"   },
        { Side::North, 0.00500, "north, normal depth (nominal)"  },
    };
    std::vector<Boundary> bcs;
    for (const auto& o : outs) {
        Boundary b;
        b.side = o.side; b.type = BcType::NormalDepth; b.slope = o.slope; b.name = o.name;
        bcs.push_back(b);
        std::printf("outflow       %-32s S = %.5f\n", o.name, o.slope);
    }
    c.bcs = bcs;
    return c;
}

#include "case_config.hpp"

int main(int argc, char** argv)
{
    // Line-buffer stdout. Redirected to a file it is block-buffered by default,
    // so a multi-hour run shows an EMPTY log for its first several minutes and
    // looks hung when it is not. Progress output is worthless if it arrives 4 KB
    // at a time.
    std::setvbuf(stdout, nullptr, _IOLBF, 0);

    std::string which = "plane";
    std::string outdir;
    double tend = -1.0;
    double tout = -1.0, cfl = -1.0, theta = -1.0;
    double dxt  = 3.0;                  // model cell size for the c6 case
    std::string aggs = "mean";
    std::string demov;
    std::string nrast;   // --manning: a direct n raster, overriding land cover
    std::string storm = "6hr";
    std::string cfgfile;
    for (int a = 1; a < argc; ++a) {
        if (std::strcmp(argv[a], "--out") == 0 && a + 1 < argc) outdir = argv[++a];
        else if (std::strcmp(argv[a], "--tend") == 0 && a + 1 < argc) tend = std::atof(argv[++a]);
        else if (std::strcmp(argv[a], "--dx")   == 0 && a + 1 < argc) dxt  = std::atof(argv[++a]);
        else if (std::strcmp(argv[a], "--agg")  == 0 && a + 1 < argc) aggs = argv[++a];
        // Output cadence matters more than it looks: .vti is ASCII, so at 3 m a
        // frame is ~60 MB and the 60 s default would write 10 GB over a 3 h run.
        else if (std::strcmp(argv[a], "--tout") == 0 && a + 1 < argc) tout = std::atof(argv[++a]);
        else if (std::strcmp(argv[a], "--cfl")  == 0 && a + 1 < argc) cfl  = std::atof(argv[++a]);
        else if (std::strcmp(argv[a], "--theta") == 0 && a + 1 < argc) theta = std::atof(argv[++a]);
        // --dem runs one case against several terrains, which is what turns the
        // culvert-notch sensitivity study into a loop rather than a rebuild.
        else if (std::strcmp(argv[a], "--dem")  == 0 && a + 1 < argc) demov = argv[++a];
        else if (std::strcmp(argv[a], "--manning") == 0 && a + 1 < argc) nrast = argv[++a];
        else if (std::strcmp(argv[a], "--storm")== 0 && a + 1 < argc) storm = argv[++a];
        else if (std::strcmp(argv[a], "--config")== 0 && a + 1 < argc) { cfgfile = argv[++a]; which = "config"; }
        else which = argv[a];
    }

    Case c;
    if (which == "config") {
#ifdef DWE_WITH_GDAL
        c = make_config(cfgfile);
#else
        std::printf("*** --config needs a GDAL build (CONFIG+=with_gdal)\n");
        return 2;
#endif
    }
    else if (which == "pond")      c = make_pond();
    else if (which == "nm141") {
#ifdef DWE_WITH_GDAL
        const Agg agg = (aggs == "min") ? Agg::Min
                      : (aggs == "max") ? Agg::Max : Agg::Mean;
        c = make_nm141(dxt, agg, demov, storm);
#else
        std::printf("*** the nm141 case needs GDAL: rebuild with  qmake6 CONFIG+=with_gdal\n");
        return 2;
#endif
    }
    else if (which == "nm032n") {
#ifdef DWE_WITH_GDAL
        const Agg agg = (aggs == "min") ? Agg::Min
                      : (aggs == "max") ? Agg::Max : Agg::Mean;
        c = make_nm032n(dxt, agg, demov, storm);
#else
        std::printf("*** the nm032n case needs GDAL: rebuild with  qmake6 CONFIG+=with_gdal\n");
        return 2;
#endif
    }
    else if (which == "c4a") {
#ifdef DWE_WITH_GDAL
        const Agg agg = (aggs == "min") ? Agg::Min
                      : (aggs == "max") ? Agg::Max : Agg::Mean;
        c = make_c4a(dxt, agg, demov, storm);
#else
        std::printf("*** the c4a case needs GDAL: rebuild with  qmake6 CONFIG+=with_gdal\n");
        return 2;
#endif
    }
    else if (which == "c4b") {
#ifdef DWE_WITH_GDAL
        const Agg agg = (aggs == "min") ? Agg::Min
                      : (aggs == "max") ? Agg::Max : Agg::Mean;
        c = make_c4b(dxt, agg, demov, storm);
#else
        std::printf("*** the c4b case needs GDAL: rebuild with  qmake6 CONFIG+=with_gdal\n");
        return 2;
#endif
    }
    else if (which == "az140") {
#ifdef DWE_WITH_GDAL
        const Agg agg = (aggs == "min") ? Agg::Min
                      : (aggs == "max") ? Agg::Max : Agg::Mean;
        c = make_az140(dxt, agg, demov, storm);
#else
        std::printf("*** the az140 case needs GDAL: rebuild with  qmake6 CONFIG+=with_gdal\n");
        return 2;
#endif
    }
    else if (which == "c7b") {
#ifdef DWE_WITH_GDAL
        const Agg agg = (aggs == "min") ? Agg::Min
                      : (aggs == "max") ? Agg::Max : Agg::Mean;
        c = make_c7b(dxt, agg, demov, storm);
#else
        std::printf("*** the c7b case needs GDAL: rebuild with  qmake6 CONFIG+=with_gdal\n");
        return 2;
#endif
    }
    else if (which == "c6") {
#ifdef DWE_WITH_GDAL
        const Agg agg = (aggs == "min") ? Agg::Min
                      : (aggs == "max") ? Agg::Max : Agg::Mean;
        if (agg == Agg::Mean && dxt > 1.5)
            std::printf("NOTE: aggregating 1 m terrain by MEAN at %.1f m lowers ridge\n"
                        "      crests and raises channel inverts -- a one-directional bias,\n"
                        "      always toward more spill. See DESIGN 4.3.\n\n", dxt);
        c = make_c6(dxt, agg, storm, demov);
#else
        std::printf("*** the c6 case needs GDAL: rebuild with  qmake6 CONFIG+=with_gdal\n");
        return 2;
#endif
    }
    else c = make_plane();
    if (tend > 0.0) c.p.t_end = tend;
    if (tout > 0.0) c.p.t_out = tout;
    if (cfl  > 0.0) c.p.cfl   = cfl;
    if (theta > 0.0) {
        if (theta > 1.0) { std::printf("*** --theta must be in (0,1]\n"); return 2; }
        c.p.theta = theta;
    }
    if (!outdir.empty()) c.p.out_dir = outdir;
    else                 c.p.out_dir = "out_" + c.name;
    std::filesystem::create_directories(c.p.out_dir);

    Fields F(c.g, c.p.n_default);
    F.z = c.z;
    F.h = c.h0;

#ifdef DWE_WITH_GDAL
    if (nrast.empty() && !c.manning.empty()) nrast = c.manning;
    if (!nrast.empty()) {
        // A roughness raster supplied directly. This is the carved-channel case:
        // the lining is 1.2 m wide and no land-cover product resolves it, so the
        // n field is built with the terrain and read straight in.
        std::string err;
        if (!load_manning_direct(nrast, c.g, F.n, c.p.n_default, &err)) {
            std::printf("*** roughness: %s\n", err.c_str());
            std::exit(2);          // asked for a specific n field: do not guess
        }
        double nmin = 1e9, nmax = -1e9, nsum = 0.0;
        std::size_t nlin = 0;
        for (double v : F.n) { nmin = std::min(nmin, v); nmax = std::max(nmax, v);
                               nsum += v; if (v < 0.020) ++nlin; }
        std::printf("roughness     n = %.3f - %.3f (mean %.4f) from %s\n",
                    nmin, nmax, nsum / F.n.size(), nrast.c_str());
        std::printf("              %zu cells at n < 0.020 -- the concrete lining\n", nlin);
    }
    else if (c.name == "c6" || c.name == "nm141" || c.name == "az140"
             || c.name == "c7b" || c.name == "c4b" || c.name == "c4a" || c.name == "nm032n") {
        std::string err;
        const bool is141 = (c.name == "nm141");
        const bool is140 = (c.name == "az140");
        const bool is7b  = (c.name == "c7b");
        const bool is4b  = (c.name == "c4b");
        const bool is4a  = (c.name == "c4a");
        const bool is32n = (c.name == "nm032n");
        const std::string lcras =
            is32n ? std::string(NM032N_DIR) + "/landcover_2d.tif"
          : is4a  ? std::string(C4A_DIR) + "/landcover_2d.tif"
          : is4b  ? std::string(C4B_DIR) + "/landcover_2d.tif"
          : is7b  ? std::string(C7B_SITE) + "/outputs_2D/landcover_2d.tif"
          : is140 ? std::string(AZ140_SITE) + "/landcover/nlcd_2023_AZ12-140_utm12n.tif"
          : is141 ? std::string(NM141_LC)
                  : std::string(C6_DIR) + "/landcover_2d.tif";
        const std::string lccsv =
            is32n ? std::string(NM032N_DIR) + "/manning_n_2d.csv"
          : is4a  ? std::string(C4A_DIR) + "/manning_n_2d.csv"
          : is4b  ? std::string(C4B_DIR) + "/manning_n_2d.csv"
          : is7b  ? std::string(C7B_SITE) + "/outputs_2D/manning_n_2d.csv"
          : is140 ? std::string(AZ140_SITE) + "/outputs_2D/manning_n_2d.csv"
          : is141 ? std::string(NM141_DIR) + "/manning_n_2d.csv"
                  : std::string(C6_DIR) + "/manning_n_2d.csv";
        if (!load_manning(lcras, lccsv, c.g, F.n, c.p.n_default, &err)) {
            std::printf("*** roughness: %s -- falling back to a uniform n = %.3f\n",
                        err.c_str(), c.p.n_default);
        } else {
            double nmin = 1e9, nmax = -1e9, nsum = 0.0;
            for (double v : F.n) { nmin = std::min(nmin, v); nmax = std::max(nmax, v); nsum += v; }
            std::printf("roughness     n = %.3f - %.3f (mean %.4f) from NLCD\n",
                        nmin, nmax, nsum / F.n.size());
        }
    }
#endif

    SchemeInertial solver(F, c.p, c.bcs);
    solver.points = c.points;

    // ---- DIRECT RAINFALL ON THE DOMAIN ------------------------------------
    // Excess precipitation, applied uniformly. This is NOT double counting the
    // boundary hydrograph: only 147 of the domain's 2,275 acres lie inside any
    // delineated subwatershed, and that subbasin (70) drains straight to the
    // Outlet without passing through Reach_58. The community bench is outside
    // the watershed altogether -- the controlling divide IS the watershed
    // boundary. See tools/hms_rainfall.py.
    if (c.name == "c6" || c.name == "nm141" || c.name == "c7b" || c.name == "c4b"
        || c.name == "nm032n") {
        const bool is141 = (c.name == "nm141");
        const bool is32n = (c.name == "nm032n");
        const bool is7b  = (c.name == "c7b");
        const bool is4b  = (c.name == "c4b");
        const std::string rcsv = is32n
            ? std::string(NM032N_DIR) + (c.storm == "24hr"
                  ? "/rain_local_24hr_100yr.csv"
                  : "/rain_local_6hr_100yr.csv")
          : is4b
            ? std::string(C4B_DIR) + (c.storm == "24hr"
                  ? "/rain_domain_24hr_100yr.csv"
                  : "/rain_domain_6hr_100yr.csv")
          : is7b
            ? std::string(C7B_SITE) + (c.storm == "24hr"
                  ? "/outputs_2D/rain_excess_24hr_100yr.csv"
                  : "/outputs_2D/rain_excess_6hr_100yr.csv")
          : is141
            ? std::string(NM141_DIR) + (c.storm == "24hr"
                  ? "/rain_excess_local_24hr_100yr.csv"
                  : "/rain_excess_local_6hr_100yr.csv")
            : std::string(C6_DIR) + (c.storm == "24hr"
                  ? "/rain_excess_24hr_100yr.csv"
                  : "/rain_excess_6hr_100yr.csv");
        solver.source = Series::from_csv(rcsv);
        if (solver.source.t.size() < 2) {
            std::printf("rainfall      NONE -- could not read %s\n", rcsv.c_str());
        } else {
            double dep = 0.0;
            for (std::size_t i = 1; i < solver.source.t.size(); ++i)
                dep += solver.source.v[i - 1] * (solver.source.t[i] - solver.source.t[i - 1]);
            std::printf("rainfall      %s\n"
                        "              %.4f in (%.1f ac-ft), peak %.3f in/h, %.1f h long\n",
                        is32n ? "LOCAL runoff of HMS Subbasins 1+2+3 as a uniform source"
                      : is141 ? "LOCAL subbasin runoff (Sub 1,2,3,57) as a uniform source"
                              : is7b ? "excess precipitation, area-weighted over 28 subbasins"
                                     : "excess precipitation, uniform over the domain",
                        dep / 0.0254, dep * c.g.ncell() * c.g.area() / 1233.4818,
                        solver.source.peak() * 3600.0 / 0.0254,
                        solver.source.span() / 3600.0);
            if (is141)
                std::printf("              volume-preserving: Subbasins 5 and 7 are EXCLUDED --\n"
                            "              they route through Reach_99 and are already in the east inflow\n");
        }
    }
    MassBalance mb;
    mb.set_initial(F.volume());

    std::printf("=====================================================================\n");
    std::printf("2D-DWE  --  Scheme A, local inertial (Bates et al. 2010)\n");
    std::printf("=====================================================================\n");
    std::printf("case          %s\n", c.name.c_str());
    std::printf("grid          %d x %d  @ %.2f m   (%zu cells)\n",
                c.g.nx, c.g.ny, c.g.dx, c.g.ncell());
    std::printf("duration      %.0f s,  CFL %.2f\n", c.p.t_end, c.p.cfl);
    if (c.p.theta < 1.0)
        std::printf("momentum      de Almeida (2012) smoothing, theta = %.3f\n"
                    "              (theta = 1 is Bates 2010, the default)\n", c.p.theta);
    std::printf("thresholds    h_dry %.4f  h_wet %.4f  h_face %.4f m\n",
                c.p.h_dry, c.p.h_wet, c.p.h_face_min);
    if (c.expect > 0)
        std::printf("expected      %.4f m  (%s)\n", c.expect, c.expect_what.c_str());
    std::printf("---------------------------------------------------------------------\n");

    double t = 0.0, t_next_out = 0.0;
    long   step = 0, nout = 0;

    // Running maximum depth, updated EVERY step.
    //
    // This used to write F.h at the end of the run and call it depth_max. On a
    // steady inflow the final state is the maximum, so the file looked right and
    // the error stayed hidden. With a real hydrograph it is not even close: the
    // 10 m run peaks at 6.62 m near t = 4 h and has fallen to 3.40 m by the end,
    // so the "maximum" was reporting barely half the flood.
    //
    // Sampling at the output cadence would not fix it either -- the peak does not
    // land on an output step. Envelope quantities must be accumulated per step.
    std::vector<double> h_max(c.g.ncell(), 0.0);
    std::vector<double> wse_max(c.g.ncell(), -9999.0);

    // ---- config cases: speed envelope, flood arrival, gauges, GeoTIFF frames ----
    // Cell speed from the four face discharges: |(qx_w+qx_e)/2, (qy_n+qy_s)/2| / h.
    auto cell_speed = [&](int i, int j) {
        const std::size_t k = c.g.k(i, j);
        if (F.h[k] <= c.p.h_dry) return 0.0;
        const double u = 0.5 * (F.qx[c.g.fx(i, j)] + F.qx[c.g.fx(i + 1, j)]);
        const double v = 0.5 * (F.qy[c.g.fy(i, j)] + F.qy[c.g.fy(i, j + 1)]);
        return std::hypot(u, v) / std::max(F.h[k], 0.05);   // floor: no spurious speed on a film
    };
    // Discharge through a boundary run, + into the domain.
    auto edge_flow = [&](const Boundary& b) {
        double q = 0.0;
        const int e = (b.end < 0) ? ((b.side == Side::West || b.side == Side::East) ? c.g.ny : c.g.nx) : b.end;
        for (int m = std::max(b.begin, 0); m < e; ++m) {
            switch (b.side) {
            case Side::West:  q +=  F.qx[c.g.fx(0, m)];         break;
            case Side::East:  q += -F.qx[c.g.fx(c.g.nx, m)];    break;
            case Side::North: q +=  F.qy[c.g.fy(m, 0)];         break;
            case Side::South: q += -F.qy[c.g.fy(m, c.g.ny)];    break;
            }
        }
        return q * c.g.dx;
    };
    std::vector<double> spd_max, t_flood;
    FILE*  gfile    = nullptr;
    double t_next_g = 0.0;
#ifdef DWE_WITH_GDAL
    if (c.is_config) {
        spd_max.assign(c.g.ncell(), 0.0);
        t_flood.assign(c.g.ncell(), -9999.0);
        if (!c.frames.empty()) std::filesystem::create_directories(c.p.out_dir + "/frames");
        gfile = std::fopen((c.p.out_dir + "/gauges.csv").c_str(), "w");
        std::fprintf(gfile, "# 2D-DWE gauges, case %s, t = 0 at %s\n", c.name.c_str(),
                     c.time0.empty() ? "(not given)" : c.time0.c_str());
        std::fprintf(gfile, "t_s");
        for (const auto& gg : c.gauges) std::fprintf(gfile, ",%s_depth_m,%s_wse_m", gg.name.c_str(), gg.name.c_str());
        std::fprintf(gfile, ",Q_in_edge_cms,Q_in_points_cms,Q_out_cms,flooded_area_m2,volume_m3\n");
    }
#endif

    while (t < c.p.t_end && step < c.p.max_steps) {
        if (t >= t_next_out && (c.p.write_vti || !c.frames.empty())) {
            char path[512];
            if (c.p.write_vti) {
                std::snprintf(path, sizeof(path), "%s/%s_%04ld.vti",
                              c.p.out_dir.c_str(), c.name.c_str(), nout);
                write_vti(path, F);
            }
#ifdef DWE_WITH_GDAL
            for (const auto& fr : c.frames) {
                std::vector<double> a(c.g.ncell(), 0.0);
                if (fr == "depth") a = F.h;
                else if (fr == "speed") {
#pragma omp parallel for schedule(static)
                    for (int j = 0; j < c.g.ny; ++j)
                        for (int i = 0; i < c.g.nx; ++i) a[c.g.k(i, j)] = cell_speed(i, j);
                }
                else continue;
                std::snprintf(path, sizeof(path), "%s/frames/%s_%04ld.tif",
                              c.p.out_dir.c_str(), fr.c_str(), nout);
                write_gtiff(path, c.g, a, c.crs);
            }
#endif
            ++nout;
            t_next_out += c.p.t_out;
        }
        if (gfile && t >= t_next_g) {
            std::fprintf(gfile, "%.1f", t);
            for (const auto& gg : c.gauges) {
                const std::size_t k = c.g.k(gg.i, gg.j);
                std::fprintf(gfile, ",%.4f,%.4f", F.h[k], F.H(k));
            }
            double qin = 0.0, qout = 0.0, qpt = 0.0;
            for (const auto& b : c.bcs) {
                if (b.type == BcType::Inflow) qin += edge_flow(b);
                else if (b.type == BcType::NormalDepth) qout -= edge_flow(b);
            }
            for (const auto& pnt : c.points) qpt += std::max(0.0, pnt.Q.at(t));
            std::size_t nfl = 0;
#pragma omp parallel for schedule(static) reduction(+:nfl)
            for (std::size_t k = 0; k < c.g.ncell(); ++k) nfl += (F.h[k] > c.flood_h);
            std::fprintf(gfile, ",%.4f,%.4f,%.4f,%.1f,%.1f\n", qin, qpt, qout,
                         nfl * c.g.area(), F.volume());
            std::fflush(gfile);
            t_next_g += c.t_gauge;
        }

        const double dt = solver.step(t, mb);
#pragma omp parallel for schedule(static)
        for (std::size_t k = 0; k < c.g.ncell(); ++k) {
            if (F.h[k] > h_max[k]) {
                h_max[k]   = F.h[k];
                wse_max[k] = F.H(k);
            }
        }
        if (c.is_config) {
            const double tn = t + dt;
#pragma omp parallel for schedule(static)
            for (int j = 0; j < c.g.ny; ++j)
                for (int i = 0; i < c.g.nx; ++i) {
                    const std::size_t k = c.g.k(i, j);
                    if (F.h[k] <= c.p.h_dry) continue;
                    const double sp = cell_speed(i, j);
                    if (sp > spd_max[k]) spd_max[k] = sp;
                    if (t_flood[k] < 0.0 && F.h[k] > c.flood_h) t_flood[k] = tn;
                }
        }
        t += dt;
        ++step;

        if (step % c.p.report_every == 0) {
            std::printf("  step %8ld  t %9.2f s  dt %8.5f  hmax %7.4f m  vol %12.3f m3\n",
                        step, t, dt, F.max_depth(), F.volume());
        }
        if (dt <= c.p.dt_min * 1.0001) {
            std::printf("*** timestep collapsed to the floor (%g s) -- stopping\n", dt);
            break;
        }
    }

    std::printf("---------------------------------------------------------------------\n");
    std::printf("finished      %ld steps, t = %.2f s\n", step, t);
    if (gfile) std::fclose(gfile);
    mb.report(t, F.volume(), c.p.mass_tol);

    // --- case-specific acceptance ---------------------------------------
    if (c.name == "plane") {
        // sample the middle of the plane, away from both boundaries
        double sum = 0.0; int cnt = 0;
        for (int j = 0; j < c.g.ny; ++j)
            for (int i = c.g.nx / 3; i < 2 * c.g.nx / 3; ++i) { sum += F.h[c.g.k(i, j)]; ++cnt; }
        const double got = sum / cnt;
        const double err = (got - c.expect) / c.expect;
        std::printf("  normal depth  expected %.4f m, got %.4f m   (%+.2f %%)  %s\n",
                    c.expect, got, 100.0 * err,
                    std::fabs(err) < 0.05 ? "OK" : "*** OUT OF TOLERANCE ***");
    }
    else if (c.name == "pond") {
        // settled level inside the bowl, and any water that got onto the rim
        double wse_bowl = 0.0, h_on_rim = 0.0;
        for (int j = 0; j < c.g.ny; ++j) {
            for (int i = 0; i < c.g.nx; ++i) {
                const std::size_t k = c.g.k(i, j);
                const double r = std::hypot(i - c.g.nx / 2.0, j - c.g.ny / 2.0);
                if (r < 38.0) wse_bowl = std::max(wse_bowl, F.H(k));
                else if (F.z[k] > 9.5) h_on_rim = std::max(h_on_rim, F.h[k]);
            }
        }
        // The pond drains asymptotically -- discharge over the notch tends to
        // zero as the head does -- so it never reaches the spill elevation in
        // finite time. The test is that it approaches from ABOVE and never sits
        // below: settling low would mean the notch was leaking.
        const double over = wse_bowl - c.expect;
        std::printf("  settled water surface in the bowl %.4f m  (notch %.4f m, +%.4f)  %s\n",
                    wse_bowl, c.expect, over,
                    (over > -0.02 && over < 0.50) ? "OK -- draining to the notch"
                                                  : "*** WRONG LEVEL ***");
        std::printf("  greatest depth anywhere on the 10.0 m rim  %.5f m   %s\n",
                    h_on_rim,
                    h_on_rim < 0.01 ? "OK -- rim held" : "*** RIM LEAKED ***");
    }

#ifdef DWE_WITH_GDAL
    // c7b added 2026-09-02. Without it the C7 runs wrote no depth_max, and the
    // maximum had to be rebuilt from the sampled .vti frames -- which UNDERSTATES
    // it, because a flashy 30-minute-Tc peak can fall between two frames. h_max
    // here is the running maximum over every timestep.
    if (c.name == "c6" || c.name == "nm141" || c.name == "az140"
        || c.name == "c7b" || c.name == "c4b" || c.name == "c4a"
        || c.name == "nm032n" || c.is_config) {
        if (c.is_config) {
            write_gtiff(c.p.out_dir + "/speed_max.tif", c.g, spd_max, c.crs);
            std::vector<double> th(c.g.ncell(), -9999.0);
            for (std::size_t k = 0; k < c.g.ncell(); ++k)
                if (t_flood[k] >= 0.0) th[k] = t_flood[k] / 3600.0;
            write_gtiff(c.p.out_dir + "/t_flood_h.tif", c.g, th, c.crs);   // hours after t = 0
        }
        for (std::size_t k = 0; k < c.g.ncell(); ++k)
            if (h_max[k] <= c.p.h_dry) wse_max[k] = -9999.0;
        write_gtiff(c.p.out_dir + "/depth_max.tif", c.g, h_max, c.crs);
        write_gtiff(c.p.out_dir + "/wse_max.tif",   c.g, wse_max, c.crs);

        std::vector<double> wse_f(c.g.ncell());
        for (std::size_t k = 0; k < c.g.ncell(); ++k)
            wse_f[k] = (F.h[k] > c.p.h_dry) ? F.H(k) : -9999.0;
        write_gtiff(c.p.out_dir + "/depth_final.tif", c.g, F.h, c.crs);
        write_gtiff(c.p.out_dir + "/wse_final.tif",   c.g, wse_f, c.crs);

        double hm = 0.0;
        for (std::size_t k = 0; k < c.g.ncell(); ++k) hm = std::max(hm, h_max[k]);
        std::printf("  peak depth anywhere over the whole run  %.3f m (%.2f ft)\n",
                    hm, hm * 3.280839895);

        // ---- DEPTH HAZARD BANDS -------------------------------------------
        // Thresholds are stated in FEET because that is the unit the report and
        // the client use, but the solver is SI throughout, so they are converted
        // here rather than anywhere near the physics.
        //
        // The classified raster stores the THRESHOLD ITSELF (0/1/3/6), not an
        // ordinal 0/1/2/3. A reviewer opening it in QGIS then reads the legend
        // directly instead of having to find a key in a document.
        static const double BAND_FT[3] = { 1.0, 3.0, 6.0 };
        const double FT   = 3.280839895;
        const double acre = c.g.dx * c.g.dx / 4046.8564224;

        std::vector<double> band(c.g.ncell(), 0.0);
        std::size_t nb[3] = { 0, 0, 0 };
        for (std::size_t k = 0; k < c.g.ncell(); ++k) {
            const double dft = h_max[k] * FT;
            for (int q = 0; q < 3; ++q)
                if (dft > BAND_FT[q]) { ++nb[q]; band[k] = BAND_FT[q]; }
        }
        std::printf("  depth bands (max over the run, area above each threshold):\n");
        for (int q = 0; q < 3; ++q)
            std::printf("      > %.0f ft : %9zu cells  %10.1f acres  %8.3f sq mi\n",
                        BAND_FT[q], nb[q], nb[q] * acre, nb[q] * acre / 640.0);
        write_gtiff(c.p.out_dir + "/depth_bands.tif", c.g, band, c.crs);
        std::printf("  wrote depth_max / wse_max / depth_final / wse_final / depth_bands"
                    " .tif to %s/\n", c.p.out_dir.c_str());
    }
#endif
    std::printf("  wrote %ld .vti frames to %s/\n", nout, c.p.out_dir.c_str());
    const bool ok = mb.closes(F.volume(), c.p.mass_tol);
    std::printf("=====================================================================\n");
    return ok ? 0 : 1;
}
