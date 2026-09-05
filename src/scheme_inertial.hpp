// =============================================================================
// scheme_inertial.hpp -- Scheme A, the local-inertial formulation.
//
//   Bates, Horritt & Fewtrell (2010), J. Hydrol. 387, 33-45.
//
// Per face, explicit in the unit discharge q [m^2/s]:
//
//                q^n  -  g*hf*dt*(dH/dx)
//     q^{n+1} = --------------------------------
//               1 + g*dt*n^2*|q^n| / hf^(7/3)
//
// then cell continuity
//
//     h^{n+1} = h^n - (dt/dx)*( qx_east - qx_west + qy_south - qy_north ) + dt*s
//
// WHY THIS AND NOT STRICT DIFFUSIVE WAVE
//   Retaining the local acceleration term makes the scheme EXPLICIT with a
//   Courant-type step, dt <= alpha*dx/sqrt(g*h). Strict DW has no inertia, so
//   its stability limit is diffusive, dt <= dx^2/(4K) with K -> infinity as the
//   water surface flattens: at 1 m cells and h = 1 m that is ~4e-5 s, which is
//   unusable. Local inertia also behaves better through Fr ~ 1, which is where
//   the steepest ~10 % of the C6 divide stations sit. See DESIGN 2 and 3.
// =============================================================================
#pragma once

#include "fields.hpp"
#include "params.hpp"
#include "boundary.hpp"
#include "massbalance.hpp"

#include <vector>

namespace dwe {

class SchemeInertial
{
public:
    SchemeInertial(Fields& f, const Params& p, std::vector<Boundary> bcs);

    // Advance one step. Returns the dt actually taken.
    double step(double t, MassBalance& mb);

    // Courant-limited timestep. Already adaptive -- recomputed every step -- but
    // it uses the WAVE + FLOW speed, not the wave speed alone, and its growth is
    // rate-limited. See the .cpp for why both matter.
    double stable_dt() const;

    // Largest (|u| + sqrt(g h)) seen in the last flux update, cached so the
    // timestep does not need its own full pass over the grid.
    double last_max_speed = 0.0;

private:
    mutable double dt_prev = 0.0;
    double src_now = 0.0;   // source rate used by the current step
public:

    // Spatially uniform source term s(t) in Equation (continuity), m/s -- direct
    // rainfall on the domain. Empty by default, which evaluates to zero.
    //
    // This is EXCESS precipitation, not gross rainfall: the solver has no
    // infiltration model, so feeding it the gross storm depth would route all of
    // it as runoff. See tools/hms_rainfall.py.
    Series source;

private:
    void update_faces(double dt);
    void apply_boundaries(double t, double dt);
    void limit_fluxes(double dt);
    void update_depths(double dt, MassBalance& mb);
    void dry_cells(MassBalance& mb);

    // Previous-step flux snapshots, used ONLY when Params::theta < 1 (the
    // de Almeida momentum smoothing). Left empty at theta = 1 so the default
    // path allocates nothing.
    std::vector<double>   qx_old, qy_old;

    Fields&               F;
    const Params&         P;
    std::vector<Boundary> BC;
};

} // namespace dwe
