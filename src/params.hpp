// =============================================================================
// params.hpp -- run configuration and the physical/numerical constants.
//
// Every threshold that can change an answer lives here, not buried in a loop.
// =============================================================================
#pragma once

#include <string>

namespace dwe {

constexpr double G = 9.80665;   // m/s^2

struct Params
{
    // --- time -----------------------------------------------------------
    double t_end      = 3600.0;  // s, simulated duration
    double cfl        = 0.7;     // alpha in  dt <= alpha*dx/sqrt(g*h)
    double dt_min     = 1.0e-4;  // s, floor -- below this the run is stopped
    double dt_max     = 10.0;    // s, ceiling regardless of what CFL allows
    double t_out      = 60.0;    // s, output cadence
    long   max_steps  = 50000000;

    // --- wetting and drying ---------------------------------------------
    //
    // Two thresholds, deliberately unequal. Equal thresholds make cells
    // flicker wet/dry every step and the solution never settles (DESIGN 5.3).
    double h_dry      = 0.002;   // m, below this a cell is dry
    double h_wet      = 0.005;   // m, must exceed this to be counted newly wet
    double h_face_min = 0.001;   // m, a face carrying less than this is closed

    // --- numerics -------------------------------------------------------
    //
    // The local-inertial friction term divides by hf^(7/3). On a very thin
    // film that denominator explodes and q is driven to zero anyway, but the
    // intermediate can overflow, so hf is floored before it is used.
    double h_fric_min = 0.001;   // m

    // --- momentum smoothing (de Almeida et al. 2012) ---------------------
    //
    // The Bates (2010) update above uses the face's OWN q^n in both the
    // numerator and the friction denominator. That is the theta = 1 case, and
    // it is the variant most prone to the chequerboard oscillation that appears
    // when the flow approaches Fr ~ 1 or spills as a thin fast film beside a
    // deep pool -- exactly the C4 Analysis A crossing, where 110 m3/s in a steep
    // incised channel is arrested against a vertical 19 ft embankment and then
    // spills over the road.
    //
    // de Almeida, Bates, Freer & Souvignet (2012), WRR 48, W05528, replace q^n
    // with a weighted average along the same axis:
    //
    //     q_hat = theta*q_i + (1-theta)/2 * (q_{i-1} + q_{i+1})
    //
    // theta = 1 reproduces Bates (2010) EXACTLY, and is the default so that no
    // existing result moves. theta ~ 0.7-0.9 damps the oscillation at small
    // cost in accuracy; lower values damp harder and start to smear real
    // gradients. Reducing the timestep does NOT substitute for this: the
    // instability is in the spatial discretisation of the momentum equation,
    // not in dt, and at C4 a 2.8x cut in CFL moved the peak depth by 1 %.
    double theta      = 1.0;     // 1.0 = Bates (2010); <1 = de Almeida smoothing

    // --- roughness ------------------------------------------------------
    double n_default  = 0.050;   // Manning n if no roughness raster is given

    // --- reporting ------------------------------------------------------
    int    report_every = 200;   // steps between console lines
    double mass_tol     = 1.0e-3;// 0.1 % -- the acceptance test (DESIGN 7)

    // --- output ---------------------------------------------------------
    std::string out_dir  = "out";
    bool        write_vti = true;
};

} // namespace dwe
