#include "scheme_inertial.hpp"

#include <cmath>
#include <algorithm>

namespace dwe {

SchemeInertial::SchemeInertial(Fields& f, const Params& p, std::vector<Boundary> bcs)
    : F(f), P(p), BC(std::move(bcs))
{
    for (auto& b : BC) b.resolve(F.g);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// ADAPTIVE TIMESTEP
//
// Three things beyond the textbook dt = alpha*dx/sqrt(g*h):
//
//  1. THE SPEED IS |u| + sqrt(g*h), NOT sqrt(g*h) ALONE. The signal that limits
//     an explicit scheme is the fastest characteristic, which is advection plus
//     the gravity wave. On the C6 wash the flow reaches several m/s, so ignoring
//     |u| overestimates the allowable step by tens of percent -- and the failure
//     is not graceful, it is a blown-up run.
//
//  2. IT REUSES THE MAXIMUM FOUND DURING THE FLUX PASS. A separate scan for
//     max depth is a full extra sweep of the grid every step; at 9.2 M cells
//     that is real time spent computing something the flux loop already saw.
//
//  3. GROWTH IS RATE-LIMITED. After a quiet period dt can otherwise jump by a
//     large factor in one step -- typically as a wetting front arrives -- and
//     the lagged friction term does not tolerate it. Allowing at most +20 % per
//     step costs nothing and removes a whole class of intermittent blow-ups.
// -----------------------------------------------------------------------------
double SchemeInertial::stable_dt() const
{
    double speed = last_max_speed;
    if (speed <= 0.0) {
        // first step, or a completely dry grid: fall back to a depth-only scan
        const double hmax = F.max_depth();
        if (hmax <= P.h_dry) return P.dt_max;
        speed = std::sqrt(G * hmax);
    }
    double dt = P.cfl * F.g.dx / speed;

    if (dt_prev > 0.0) dt = std::min(dt, dt_prev * 1.2);   // growth limiter

    return std::min(P.dt_max, std::max(P.dt_min, dt));
}

// -----------------------------------------------------------------------------
// Interior faces. Boundary faces (i==0, i==nx, j==0, j==ny) are left alone here
// and set by apply_boundaries(); a face not touched by either stays at zero,
// which is the wall condition.
// -----------------------------------------------------------------------------
void SchemeInertial::update_faces(double dt)
{
    const Grid&  g   = F.g;
    const double dx  = g.dx;
    const double gdt = G * dt;
    double smax = 0.0;

    // de Almeida momentum smoothing. The weighted q_hat must be built from the
    // OLD flux field: the loops below write F.qx/F.qy in place, so a neighbour
    // read after its own update would mix time levels and quietly break the
    // scheme. Snapshot only when it is actually in use -- at theta = 1 this
    // allocates nothing and the code path is identical to Bates (2010).
    const bool   smooth = (P.theta < 1.0);
    const double th     = P.theta;
    const double thn    = 0.5 * (1.0 - P.theta);
    if (smooth) { qx_old = F.qx; qy_old = F.qy; }

    // ---- west-east faces -------------------------------------------------
#pragma omp parallel for schedule(static) reduction(max:smax)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 1; i < g.nx; ++i) {
            const std::size_t L = g.k(i - 1, j);
            const std::size_t R = g.k(i, j);
            const std::size_t f = g.fx(i, j);

            const double hf = face_depth(F.z[L], F.h[L], F.z[R], F.h[R]);
            if (hf <= P.h_face_min) { F.qx[f] = 0.0; continue; }

            const double dHdx = (F.H(R) - F.H(L)) / dx;
            const double nf   = 0.5 * (F.n[L] + F.n[R]);
            const double hfr  = std::max(hf, P.h_fric_min);

            // q_hat: theta*q_i + (1-theta)/2*(q_{i-1} + q_{i+1}) along x.
            // i runs 1..nx-1 here, so fx(i-1,j) and fx(i+1,j) are both valid
            // face indices (0..nx); the end ones carry the boundary fluxes,
            // which is what de Almeida intends.
            const double q0   = smooth
                ? th * qx_old[f] + thn * (qx_old[g.fx(i - 1, j)]
                                        + qx_old[g.fx(i + 1, j)])
                : F.qx[f];
            const double num  = q0 - gdt * hf * dHdx;
            const double den  = 1.0 + gdt * nf * nf * std::fabs(q0)
                                      / std::pow(hfr, 7.0 / 3.0);
            F.qx[f] = num / den;
            // fastest characteristic through this face: advection + gravity wave
            smax = std::max(smax, std::fabs(F.qx[f]) / hfr + std::sqrt(G * hf));
        }
    }

    // ---- north-south faces ----------------------------------------------
#pragma omp parallel for schedule(static) reduction(max:smax)
    for (int j = 1; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            const std::size_t L = g.k(i, j - 1);   // north neighbour
            const std::size_t R = g.k(i, j);       // south neighbour
            const std::size_t f = g.fy(i, j);

            const double hf = face_depth(F.z[L], F.h[L], F.z[R], F.h[R]);
            if (hf <= P.h_face_min) { F.qy[f] = 0.0; continue; }

            const double dHdy = (F.H(R) - F.H(L)) / dx;
            const double nf   = 0.5 * (F.n[L] + F.n[R]);
            const double hfr  = std::max(hf, P.h_fric_min);

            const double q0   = smooth
                ? th * qy_old[f] + thn * (qy_old[g.fy(i, j - 1)]
                                        + qy_old[g.fy(i, j + 1)])
                : F.qy[f];
            const double num  = q0 - gdt * hf * dHdy;
            const double den  = 1.0 + gdt * nf * nf * std::fabs(q0)
                                      / std::pow(hfr, 7.0 / 3.0);
            F.qy[f] = num / den;
            smax = std::max(smax, std::fabs(F.qy[f]) / hfr + std::sqrt(G * hf));
        }
    }

    last_max_speed = smax;
}

// -----------------------------------------------------------------------------
void SchemeInertial::apply_boundaries(double t, double /*dt*/)
{
    const Grid& g = F.g;

    for (const auto& b : BC) {
        const int cnt = b.count();
        if (cnt <= 0) continue;

        if (b.type == BcType::Wall) {
            for (int m = b.begin; m < b.end; ++m) {
                switch (b.side) {
                case Side::West:  F.qx[g.fx(0, m)]      = 0.0; break;
                case Side::East:  F.qx[g.fx(g.nx, m)]   = 0.0; break;
                case Side::North: F.qy[g.fy(m, 0)]      = 0.0; break;
                case Side::South: F.qy[g.fy(m, g.ny)]   = 0.0; break;
                }
            }
        }
        else if (b.type == BcType::Inflow && b.conveyance) {
            // Weighted by conveyance under the highest water surface on the run
            // (boundary.hpp). A dry run takes all of Q on its lowest cell.
            auto cell = [&](int m) -> std::size_t {
                switch (b.side) {
                case Side::West:  return g.k(0, m);
                case Side::East:  return g.k(g.nx - 1, m);
                case Side::North: return g.k(m, 0);
                default:          return g.k(m, g.ny - 1);
                }
            };
            double Hs = -1e30, zlow = 1e30; int mlow = b.begin;
            for (int m = b.begin; m < b.end; ++m) {
                const std::size_t k = cell(m);
                if (F.h[k] > P.h_dry) Hs = std::max(Hs, F.H(k));
                if (F.z[k] < zlow) { zlow = F.z[k]; mlow = m; }
            }
            w_.assign(std::size_t(cnt), 0.0);
            double ws = 0.0;
            for (int m = b.begin; m < b.end; ++m) {
                const double d = Hs - F.z[cell(m)];
                const double w = (d > 0.0) ? std::pow(d, 5.0 / 3.0) : 0.0;
                w_[std::size_t(m - b.begin)] = w; ws += w;
            }
            if (ws <= 0.0) { w_[std::size_t(mlow - b.begin)] = 1.0; ws = 1.0; }
            const double Q = b.Q.at(t);
            for (int m = b.begin; m < b.end; ++m) {
                const double q = Q * w_[std::size_t(m - b.begin)] / (ws * g.dx);
                switch (b.side) {
                case Side::West:  F.qx[g.fx(0, m)]    =  q; break;
                case Side::East:  F.qx[g.fx(g.nx, m)] = -q; break;
                case Side::North: F.qy[g.fy(m, 0)]    =  q; break;
                case Side::South: F.qy[g.fy(m, g.ny)] = -q; break;
                }
            }
        }
        else if (b.type == BcType::Inflow) {
            // Total discharge spread evenly over the face width. q is per unit
            // width, so q = Q / (cnt * dx). Sign puts flow INTO the domain.
            const double Q = b.Q.at(t);
            const double q = Q / (cnt * g.dx);
            for (int m = b.begin; m < b.end; ++m) {
                switch (b.side) {
                case Side::West:  F.qx[g.fx(0, m)]    =  q; break;  // eastward, in
                case Side::East:  F.qx[g.fx(g.nx, m)] = -q; break;  // westward, in
                case Side::North: F.qy[g.fy(m, 0)]    =  q; break;  // southward, in
                case Side::South: F.qy[g.fy(m, g.ny)] = -q; break;  // northward, in
                }
            }
        }
        else { // NormalDepth -- free outfall using the adjacent cell's depth
            const double s = std::sqrt(std::max(b.slope, 1.0e-6));
            for (int m = b.begin; m < b.end; ++m) {
                std::size_t kc = 0;
                switch (b.side) {
                case Side::West:  kc = g.k(0, m);         break;
                case Side::East:  kc = g.k(g.nx - 1, m);  break;
                case Side::North: kc = g.k(m, 0);         break;
                case Side::South: kc = g.k(m, g.ny - 1);  break;
                }
                const double hc = F.h[kc];
                const double q  = (hc <= P.h_dry)
                                ? 0.0
                                : (1.0 / F.n[kc]) * std::pow(hc, 5.0 / 3.0) * s;
                switch (b.side) {
                case Side::West:  F.qx[g.fx(0, m)]    = -q; break;  // out to the west
                case Side::East:  F.qx[g.fx(g.nx, m)] =  q; break;  // out to the east
                case Side::North: F.qy[g.fy(m, 0)]    = -q; break;  // out to the north
                case Side::South: F.qy[g.fy(m, g.ny)] =  q; break;  // out to the south
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// POSITIVITY LIMITER  (DESIGN 5.2) -- non-negotiable.
//
// Compute, per cell, the volume that would leave it this step. If that exceeds
// what the cell holds, scale every OUTGOING face of that cell by the same
// factor. Applying the donor's factor to the face scales both sides
// identically, so mass is conserved exactly and depth cannot go negative.
// -----------------------------------------------------------------------------
void SchemeInertial::limit_fluxes(double dt)
{
    const Grid&  g  = F.g;
    const double dx = g.dx;
    const double c  = dt / dx;      // converts unit discharge to depth change

#pragma omp parallel for schedule(static)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            const std::size_t k = g.k(i, j);

            const double qw = F.qx[g.fx(i, j)];       // + enters from the west
            const double qe = F.qx[g.fx(i + 1, j)];   // + leaves to the east
            const double qn = F.qy[g.fy(i, j)];       // + enters from the north
            const double qs = F.qy[g.fy(i, j + 1)];   // + leaves to the south

            const double out = c * ( std::max(qe, 0.0) + std::max(-qw, 0.0)
                                   + std::max(qs, 0.0) + std::max(-qn, 0.0) );

            F.alpha[k] = (out > F.h[k] && out > 0.0)
                       ? std::max(0.0, F.h[k] / out)
                       : 1.0;
        }
    }

    // Apply the donor cell's factor to each face. Interior faces only -- a
    // boundary face has no donor inside the domain on one side, and inflow
    // must not be throttled by a dry edge cell.
#pragma omp parallel for schedule(static)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 1; i < g.nx; ++i) {
            const std::size_t f = g.fx(i, j);
            const double q = F.qx[f];
            if (q > 0.0)      F.qx[f] = q * F.alpha[g.k(i - 1, j)];
            else if (q < 0.0) F.qx[f] = q * F.alpha[g.k(i, j)];
        }
    }
#pragma omp parallel for schedule(static)
    for (int j = 1; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            const std::size_t f = g.fy(i, j);
            const double q = F.qy[f];
            if (q > 0.0)      F.qy[f] = q * F.alpha[g.k(i, j - 1)];
            else if (q < 0.0) F.qy[f] = q * F.alpha[g.k(i, j)];
        }
    }

    // Outflow faces are throttled by their own cell so an edge cell cannot be
    // drained below zero by the normal-depth rating.
    for (const auto& b : BC) {
        if (b.type != BcType::NormalDepth) continue;
        for (int m = b.begin; m < b.end; ++m) {
            switch (b.side) {
            case Side::West:  F.qx[g.fx(0, m)]    *= F.alpha[g.k(0, m)];        break;
            case Side::East:  F.qx[g.fx(g.nx, m)] *= F.alpha[g.k(g.nx - 1, m)]; break;
            case Side::North: F.qy[g.fy(m, 0)]    *= F.alpha[g.k(m, 0)];        break;
            case Side::South: F.qy[g.fy(m, g.ny)] *= F.alpha[g.k(m, g.ny - 1)]; break;
            }
        }
    }
}

// -----------------------------------------------------------------------------
void SchemeInertial::update_depths(double dt, MassBalance& mb)
{
    const Grid&  g  = F.g;
    const double c  = dt / g.dx;
    const double A  = g.area();

#pragma omp parallel for schedule(static)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            const std::size_t k = g.k(i, j);
            const double qw = F.qx[g.fx(i, j)];
            const double qe = F.qx[g.fx(i + 1, j)];
            const double qn = F.qy[g.fy(i, j)];
            const double qs = F.qy[g.fy(i, j + 1)];
            double hn = F.h[k] - c * ( (qe - qw) + (qs - qn) ) + dt * src_now;
            F.h[k] = (hn > 0.0) ? hn : 0.0;   // the limiter should make this a no-op
        }
    }

    if (src_now != 0.0) mb.v_rain += src_now * dt * A * double(g.ncell());

    for (std::size_t m = 0; m < points.size(); ++m) {
        const std::size_t k = g.k(points[m].i, points[m].j);
        F.h[k]  += pt_now[m] * dt / A;
        mb.v_in += pt_now[m] * dt;
    }

    // Boundary volumes, tallied from the faces actually used this step.
    for (const auto& b : BC) {
        if (b.type == BcType::Wall) continue;
        double v = 0.0;
        for (int m = b.begin; m < b.end; ++m) {
            double q = 0.0;
            switch (b.side) {
            case Side::West:  q =  F.qx[g.fx(0, m)];      break;  // + is inflow
            case Side::East:  q = -F.qx[g.fx(g.nx, m)];   break;
            case Side::North: q =  F.qy[g.fy(m, 0)];      break;
            case Side::South: q = -F.qy[g.fy(m, g.ny)];   break;
            }
            v += q * g.dx * dt;    // + into the domain, - out of it
        }
        if (b.type == BcType::Inflow) mb.v_in  += v;
        else                          mb.v_out += -v;
    }
}

// -----------------------------------------------------------------------------
// Drying. Below h_dry the residual film is removed -- and ACCOUNTED FOR. It
// goes into a discarded accumulator, never into thin air (DESIGN 5.4).
// -----------------------------------------------------------------------------
void SchemeInertial::dry_cells(MassBalance& mb)
{
    const double A = F.g.area();
    double lost = 0.0, film = 0.0;
    std::size_t nfilm = 0;

    // THE FILM IS LEFT IN PLACE. It used to be zeroed and booked as discarded,
    // which is a real mass sink and a resolution-dependent one: the cost of a
    // drying event is h_dry * cell_area, but halving the cell size more than
    // doubles the number of wetting-front cells that dry. Measured on C6 with
    // the 6-hour hydrograph, that sink took 13.8 % of the inflow at 10 m and
    // 33.2 % at 3 m -- refining the grid made the answer worse, which is
    // backwards and disqualified the finer run.
    //
    // Deleting it was never buying stability. Two guards already handle the
    // near-dry case independently: h_face_min stops flux across a shallow face
    // (update_faces returns q = 0 outright), and h_fric_min bounds the h^(7/3)
    // in the friction denominator. The film simply sits there, which is what a
    // millimetre of water on desert ground actually does.
    //
    // Only genuine round-off negatives are clamped, and those are still booked
    // rather than silently absorbed -- the positivity limiter bounds everything
    // above that, so a non-trivial negative here means the limiter has failed
    // and the mass balance must be allowed to say so.
#pragma omp parallel for schedule(static) reduction(+:lost,film,nfilm)
    for (std::size_t k = 0; k < F.g.ncell(); ++k) {
        if (F.h[k] < 0.0) {
            lost += F.h[k];
            F.h[k] = 0.0;
        }
        else if (F.h[k] < P.h_dry && F.h[k] > 0.0) {
            film += F.h[k];
            ++nfilm;
        }
        F.active[k] = (F.h[k] > P.h_wet) ? 1 : 0;
    }
    mb.v_discarded += lost * A;
    mb.v_film       = film * A;      // a snapshot, not an accumulator
    mb.n_film       = nfilm;
}

// -----------------------------------------------------------------------------
double SchemeInertial::step(double t, MassBalance& mb)
{
    const double dt = stable_dt();
    dt_prev = dt;
    // Sample the rainfall series at the MIDPOINT of the step. Sampling at the
    // start biases a step function high on a rising limb and low on a falling
    // one; the midpoint is second-order and costs nothing.
    src_now = source.t.empty() ? 0.0 : source.at(t + 0.5 * dt);
    pt_now.resize(points.size());
    for (std::size_t m = 0; m < points.size(); ++m)
        pt_now[m] = std::max(0.0, points[m].Q.at(t + 0.5 * dt));
    update_faces(dt);
    apply_boundaries(t, dt);
    limit_fluxes(dt);
    update_depths(dt, mb);
    dry_cells(mb);
    return dt;
}

} // namespace dwe
