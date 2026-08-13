// =============================================================================
// massbalance.hpp -- the acceptance test, not a debugging aid (DESIGN 7).
//
//     inflow - outflow - dStorage - discarded  <  tol * inflow
//
// A solver that produces a plausible-looking depth map with a 5 % mass error is
// worse than one that crashes, because the plausible map gets used.
// =============================================================================
#pragma once

#include <cstdio>
#include <cmath>
#include <algorithm>

namespace dwe {

struct MassBalance
{
    double v_in        = 0.0;   // m^3, cumulative through inflow boundaries
    double v_out       = 0.0;   // m^3, cumulative through outflow boundaries
    double v_rain      = 0.0;   // m^3, cumulative source term
    double v_discarded = 0.0;   // m^3, round-off negatives clamped to zero
    double v_film      = 0.0;   // m^3, water currently held below h_dry (NOT lost)
    std::size_t n_film = 0;     // cells holding that film
    double v_initial   = 0.0;   // m^3, storage at t=0

    void set_initial(double v) { v_initial = v; }

    // signed closure error, m^3
    double error(double v_now) const
    {
        return (v_initial + v_in + v_rain) - (v_out + v_discarded + v_now);
    }

    double relative_error(double v_now) const
    {
        const double scale = std::max({std::fabs(v_initial + v_in + v_rain), 1.0});
        return error(v_now) / scale;
    }

    bool closes(double v_now, double tol) const
    {
        return std::fabs(relative_error(v_now)) < tol;
    }

    void report(double t, double v_now, double tol) const
    {
        const double e = error(v_now);
        const double r = relative_error(v_now);
        std::printf("  mass balance @ t=%10.2f s :"
                    " in %12.3f  rain %10.3f  out %12.3f  stored %12.3f  discarded %10.4f"
                    "  film %10.1f (%zu cells)"
                    "  err %11.4e m3 (%8.3e)  %s\n",
                    t, v_in, v_rain, v_out, v_now, v_discarded, v_film, n_film, e, r,
                    std::fabs(r) < tol ? "OK" : "*** FAILS ***");
    }
};

} // namespace dwe
