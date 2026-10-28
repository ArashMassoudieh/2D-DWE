#!/usr/bin/env python3
"""
Culvert / crossing sizing for Navajo Route 9 at the C10 wash crossing.

FHWA HDS-5 inlet control + outlet control, English units. The governing
headwater is the larger of the two, per HDS-5.

Design flow is the converged HEC-HMS 100-year 6-hour peak (the governing storm).
"""
import math

G = 32.174

# ---- site ------------------------------------------------------------------
Q100      = 11922.0     # cfs, converged HMS 6-hr 100-yr outlet peak
INV_IN    = 6112.8      # ft, inlet invert (report Table 7, from 1 m LiDAR)
INV_OUT   = 6110.6      # ft, outlet invert
ROAD_LOW  = 6120.6      # ft, road low point (report Table 7)
CREST_CRN = 6122.7      # ft, embankment crown (derived, tools/road_overtopping.py)
L         = 130.0       # ft, barrel length
S0        = (INV_IN - INV_OUT) / L
FREEBOARD = 1.0         # ft below the controlling crest

# ---- HDS-5 inlet-control coefficients --------------------------------------
# submerged form:  HW/D = c (Q/(A D^0.5))^2 + Y - 0.5 S
INLET = {
    "CMP, headwall":            dict(c=0.0379, Y=0.69, ke=0.50, n=0.024, shape="circ"),
    "CMP, mitered to slope":    dict(c=0.0463, Y=0.75, ke=0.70, n=0.024, shape="circ"),
    "RCP, groove end w/ HW":    dict(c=0.0292, Y=0.74, ke=0.20, n=0.012, shape="circ"),
    "RCB, 90 deg HW chamfered": dict(c=0.0347, Y=0.81, ke=0.50, n=0.012, shape="box"),
}


def area_circ(D):   return math.pi * D * D / 4.0
def area_box(W, H): return W * H


def hw_inlet(Q, A, D, c, Y, S):
    """HDS-5 submerged inlet control, ft of head above the inlet invert."""
    x = Q / (A * math.sqrt(D))
    return D * (c * x * x + Y - 0.5 * S)


def hw_outlet(Q, A, R, D, ke, n, L, S0, TW):
    """HDS-5 outlet control: HW = TW + H - S0 L, with H the total head loss."""
    V = Q / A
    H = (1.0 + ke + (29.0 * n * n * L) / R ** (4.0 / 3.0)) * V * V / (2.0 * G)
    tw = max(TW, (D + critical_depth_proxy(D)) / 2.0)   # HDS-5 ho convention
    return tw + H - S0 * L


def critical_depth_proxy(D):
    return 0.7 * D          # adequate for the near-full barrels screened here


def required_cells(Q, geom, kind, hw_allow):
    """Smallest whole number of identical barrels that keeps HW <= hw_allow."""
    k = INLET[kind]
    for N in range(1, 400):
        q = Q / N
        if k["shape"] == "circ":
            D = geom
            A = area_circ(D); R = D / 4.0
        else:
            W, H = geom
            D = H
            A = area_box(W, H); R = A / (2.0 * (W + H))
        hi = hw_inlet(q, A, D, k["c"], k["Y"], S0)
        ho = hw_outlet(q, A, R, D, k["ke"], k["n"], L, S0, INV_OUT - INV_IN + 2.0)
        hw = max(hi, ho)
        if hw <= hw_allow:
            return N, hw, hi, ho, A * N
    return None, None, None, None, None


def main():
    print("C10 / Navajo Route 9 crossing -- culvert sizing (HDS-5)\n")
    print("design flow Q100 (6-hr, converged HMS) : %8.0f cfs" % Q100)
    print("inlet invert                           : %8.1f ft" % INV_IN)
    print("barrel slope S0                        : %8.4f  (%.2f%%)" % (S0, 100 * S0))
    for lbl, crest in (("road low point", ROAD_LOW), ("embankment crown", CREST_CRN)):
        print("allowable HW to %-18s: %8.2f ft  (crest %.1f - invert %.1f - %.1f ft freeboard)"
              % (lbl, crest - INV_IN - FREEBOARD, crest, INV_IN, FREEBOARD))
    hw_allow = ROAD_LOW - INV_IN - FREEBOARD
    print("\nsizing against the road low point: HW_allow = %.2f ft\n" % hw_allow)

    print("%-28s %-16s %6s %9s %9s %9s %11s" %
          ("entrance type", "barrel", "N", "HW (ft)", "HWi", "HWo", "area (sf)"))
    print("-" * 96)
    rows = []
    for kind in INLET:
        sizes = ([(W, H) for W, H in ((8, 6), (10, 8), (12, 10))]
                 if INLET[kind]["shape"] == "box" else [4.0, 6.0, 8.0, 10.0])
        for g in sizes:
            N, hw, hi, ho, atot = required_cells(Q100, g, kind, hw_allow)
            gl = ("%g x %g ft box" % g) if isinstance(g, tuple) else ("%g ft dia" % g)
            if N is None:
                print("%-28s %-16s %6s %9s %9s %9s %11s" % (kind, gl, ">400", "-", "-", "-", "-"))
                continue
            width = (N * g[0]) if isinstance(g, tuple) else (N * g)
            rows.append((kind, gl, N, hw, atot, width))
            print("%-28s %-16s %6d %9.2f %9.2f %9.2f %11.0f" % (kind, gl, N, hw, hi, ho, atot))

    print("\ntotal waterway width needed (barrel width x N, no pier allowance):")
    for kind, gl, N, hw, atot, width in sorted(rows, key=lambda r: r[5])[:6]:
        print("   %-28s %-16s N=%-4d %7.0f ft of opening" % (kind, gl, N, width))

    print("\nfor comparison, the structure that is there now:")
    for D, nb, lbl in ((40 / 12.0, 2, "as tabulated for culvert 8"),
                       (85 / 12.0, 2, "as photographed 17 Aug 2026")):
        k = INLET["CMP, headwall"]; A = area_circ(D)
        # capacity at the allowable head, inverting the inlet-control form
        x = math.sqrt(max((hw_allow / D - k["Y"] + 0.5 * S0) / k["c"], 0.0))
        q = x * A * math.sqrt(D)
        print("   %d x %.2f ft CMP (%-28s): %6.0f cfs = %4.1f%% of Q100"
              % (nb, D, lbl, nb * q, 100 * nb * q / Q100))


if __name__ == "__main__":
    main()
