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

static Case make_c6(double dx_target, Agg agg)
{
    Case c;
    c.name = "c6";

    const std::string terr = std::string(C6_DIR) + "/terrain_2d.tif";
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

    c.p.n_default   = 0.050;
    c.p.t_end       = 9.0 * 3600.0;    // the 6-hour storm's control window
    c.p.t_out       = 600.0;
    c.p.report_every = 500;
    c.p.dt_max      = 5.0;

    c.h0.assign(c.g.ncell(), 0.0);

    // --- inflow: 32,628 cfs on the west edge, on the wash ----------------
    //
    // The wash crosses the west edge at y = 3,959,857 (Phase 1). Convert that to
    // a row index on whatever grid the terrain gave us, and spread the discharge
    // over a face centred on it.
    const double y_wash   = 3959857.0;
    const double half_m   = 150.0;
    const int j_wash = int((c.g.y0 - y_wash) / c.g.dx);
    const int jhalf  = std::max(1, int(half_m / c.g.dx));

    Boundary in;
    in.side  = Side::West;
    in.type  = BcType::Inflow;
    in.begin = std::max(0, j_wash - jhalf);
    in.end   = std::min(c.g.ny, j_wash + jhalf);
    in.Q     = Series::constant(32628.0 * 0.028316847);   // cfs -> m3/s
    in.name  = "Figueredo Wash, 6-hr 100-yr";

    Boundary out;
    out.side  = Side::East;
    out.type  = BcType::NormalDepth;
    out.slope = 0.005;
    out.name  = "east edge outfall";

    c.bcs = {in, out};

    std::printf("inflow        %.1f m3/s (32,628 cfs) over rows %d-%d on the west edge\n",
                in.Q.at(0.0), in.begin, in.end);
    std::printf("outflow       normal depth, full east edge, S = %.4f\n", out.slope);

    // --- roughness -------------------------------------------------------
    // load_manning fills a vector on this grid; stash it in h0's sibling by
    // returning through the Case. Fields are built by the caller, so hand the
    // roughness back via a static and let main() apply it.
    return c;
}
#endif // DWE_WITH_GDAL

// -----------------------------------------------------------------------------
int main(int argc, char** argv)
{
    std::string which = "plane";
    std::string outdir;
    double tend = -1.0;
    double dxt  = 3.0;                  // model cell size for the c6 case
    std::string aggs = "mean";
    for (int a = 1; a < argc; ++a) {
        if (std::strcmp(argv[a], "--out") == 0 && a + 1 < argc) outdir = argv[++a];
        else if (std::strcmp(argv[a], "--tend") == 0 && a + 1 < argc) tend = std::atof(argv[++a]);
        else if (std::strcmp(argv[a], "--dx")   == 0 && a + 1 < argc) dxt  = std::atof(argv[++a]);
        else if (std::strcmp(argv[a], "--agg")  == 0 && a + 1 < argc) aggs = argv[++a];
        else which = argv[a];
    }

    Case c;
    if (which == "pond")      c = make_pond();
    else if (which == "c6") {
#ifdef DWE_WITH_GDAL
        const Agg agg = (aggs == "min") ? Agg::Min
                      : (aggs == "max") ? Agg::Max : Agg::Mean;
        if (agg == Agg::Mean && dxt > 1.5)
            std::printf("NOTE: aggregating 1 m terrain by MEAN at %.1f m lowers ridge\n"
                        "      crests and raises channel inverts -- a one-directional bias,\n"
                        "      always toward more spill. See DESIGN 4.3.\n\n", dxt);
        c = make_c6(dxt, agg);
#else
        std::printf("*** the c6 case needs GDAL: rebuild with  qmake6 CONFIG+=with_gdal\n");
        return 2;
#endif
    }
    else c = make_plane();
    if (tend > 0.0) c.p.t_end = tend;
    if (!outdir.empty()) c.p.out_dir = outdir;
    else                 c.p.out_dir = "out_" + c.name;
    std::filesystem::create_directories(c.p.out_dir);

    Fields F(c.g, c.p.n_default);
    F.z = c.z;
    F.h = c.h0;

#ifdef DWE_WITH_GDAL
    if (c.name == "c6") {
        std::string err;
        if (!load_manning(std::string(C6_DIR) + "/landcover_2d.tif",
                          std::string(C6_DIR) + "/manning_n_2d.csv",
                          c.g, F.n, c.p.n_default, &err)) {
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
    MassBalance mb;
    mb.set_initial(F.volume());

    std::printf("=====================================================================\n");
    std::printf("2D-DWE  --  Scheme A, local inertial (Bates et al. 2010)\n");
    std::printf("=====================================================================\n");
    std::printf("case          %s\n", c.name.c_str());
    std::printf("grid          %d x %d  @ %.2f m   (%zu cells)\n",
                c.g.nx, c.g.ny, c.g.dx, c.g.ncell());
    std::printf("duration      %.0f s,  CFL %.2f\n", c.p.t_end, c.p.cfl);
    std::printf("thresholds    h_dry %.4f  h_wet %.4f  h_face %.4f m\n",
                c.p.h_dry, c.p.h_wet, c.p.h_face_min);
    if (c.expect > 0)
        std::printf("expected      %.4f m  (%s)\n", c.expect, c.expect_what.c_str());
    std::printf("---------------------------------------------------------------------\n");

    double t = 0.0, t_next_out = 0.0;
    long   step = 0, nout = 0;
    std::vector<double> wse(c.g.ncell());

    while (t < c.p.t_end && step < c.p.max_steps) {
        if (c.p.write_vti && t >= t_next_out) {
            for (std::size_t k = 0; k < c.g.ncell(); ++k) wse[k] = F.H(k);
            char path[512];
            std::snprintf(path, sizeof(path), "%s/%s_%04ld.vti",
                          c.p.out_dir.c_str(), c.name.c_str(), nout);
            write_vti(path, c.g, F.h, F.z, wse);
            ++nout;
            t_next_out += c.p.t_out;
        }

        const double dt = solver.step(t, mb);
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
    if (c.name == "c6") {
        write_gtiff(c.p.out_dir + "/depth_max.tif", c.g, F.h);
        std::vector<double> wse_f(c.g.ncell());
        for (std::size_t k = 0; k < c.g.ncell(); ++k)
            wse_f[k] = (F.h[k] > c.p.h_dry) ? F.H(k) : -9999.0;
        write_gtiff(c.p.out_dir + "/wse_final.tif", c.g, wse_f);
        std::printf("  wrote depth_max.tif and wse_final.tif to %s/\n", c.p.out_dir.c_str());
    }
#endif
    std::printf("  wrote %ld .vti frames to %s/\n", nout, c.p.out_dir.c_str());
    const bool ok = mb.closes(F.volume(), c.p.mass_tol);
    std::printf("=====================================================================\n");
    return ok ? 0 : 1;
}
