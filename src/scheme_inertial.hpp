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

    // Courant-limited timestep from the current depth field.
    double stable_dt() const;

    // Optional uniform source (rainfall etc.), m/s. Zero by default.
    double source = 0.0;

private:
    void update_faces(double dt);
    void apply_boundaries(double t, double dt);
    void limit_fluxes(double dt);
    void update_depths(double dt, MassBalance& mb);
    void dry_cells(MassBalance& mb);

    Fields&               F;
    const Params&         P;
    std::vector<Boundary> BC;
};

} // namespace dwe
