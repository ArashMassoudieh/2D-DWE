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
double SchemeInertial::stable_dt() const
{
    const double hmax = F.max_depth();
    if (hmax <= P.h_dry) return P.dt_max;
    const double dt = P.cfl * F.g.dx / std::sqrt(G * hmax);
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

    // ---- west-east faces -------------------------------------------------
#pragma omp parallel for schedule(static)
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

            const double q0   = F.qx[f];
            const double num  = q0 - gdt * hf * dHdx;
            const double den  = 1.0 + gdt * nf * nf * std::fabs(q0)
                                      / std::pow(hfr, 7.0 / 3.0);
            F.qx[f] = num / den;
        }
    }

    // ---- north-south faces ----------------------------------------------
#pragma omp parallel for schedule(static)
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

            const double q0   = F.qy[f];
            const double num  = q0 - gdt * hf * dHdy;
            const double den  = 1.0 + gdt * nf * nf * std::fabs(q0)
                                      / std::pow(hfr, 7.0 / 3.0);
            F.qy[f] = num / den;
        }
    }
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
            double hn = F.h[k] - c * ( (qe - qw) + (qs - qn) ) + dt * source;
            F.h[k] = (hn > 0.0) ? hn : 0.0;   // the limiter should make this a no-op
        }
    }

    if (source != 0.0) mb.v_rain += source * dt * A * double(g.ncell());

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
    double lost = 0.0;

#pragma omp parallel for schedule(static) reduction(+:lost)
    for (std::size_t k = 0; k < F.g.ncell(); ++k) {
        if (F.h[k] > 0.0 && F.h[k] < P.h_dry) {
            lost += F.h[k];
            F.h[k] = 0.0;
        }
        F.active[k] = (F.h[k] > P.h_wet) ? 1 : 0;
    }
    mb.v_discarded += lost * A;
}

// -----------------------------------------------------------------------------
double SchemeInertial::step(double t, MassBalance& mb)
{
    const double dt = stable_dt();
    update_faces(dt);
    apply_boundaries(t, dt);
    limit_fluxes(dt);
    update_depths(dt, mb);
    dry_cells(mb);
    return dt;
}

} // namespace dwe
