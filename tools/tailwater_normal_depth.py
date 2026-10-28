#!/usr/bin/env python3
"""
Tailwater at the C10 Route 9 crossing, at the OUTLET section.

    python3 tools/tailwater_normal_depth.py [Q_cfs]

Reproduces the tailwater of Report_for_C10/culvert_sizing_C10.tex and of the
Tailwater sheet of C10_crossing_sizing.xlsx. Default Q is the converged HMS
6-hour 100-year peak.

WHY NOT TAKE IT FROM THE 2-D MODEL
  The 2-D domain ends 96 m north of the road and its outflow boundary is itself
  a normal-depth condition with an ASSUMED friction slope. Sampling the 2-D water
  surface just upstream of that boundary mostly recovers the boundary assumption.
  It is also an envelope (wse_max), so upstream and downstream maxima need not be
  simultaneous. This computes the tailwater from terrain and Manning/critical
  flow instead; nothing here depends on the 2-D run. The 2-D value is printed
  only as an independent corroboration.

TWO CORRECTIONS OVER THE EARLIER VERSION OF THIS SCRIPT (see .BAK)
  1. SECTION. The old version took a section 150 m DOWNSTREAM of the road, whose
     invert (6,097.1 ft) is ~13 ft below the culvert outlet invert (6,110.6 ft),
     and then compared that water surface against the crossing's inlet invert.
     That mixes datums and understates the tailwater. The controlling section is
     the outlet section itself, at y = 3,960,500.
  2. BED SLOPE. The old version least-squares fitted the thalweg over a 900 m
     reach inside an 800 m-wide window. In that window the "lowest cell per row"
     jumps between the wash and unconnected depressions, so the fit returned
     S0 = 0.00063 although the same thalweg falls 19.6 ft over the reach
     (~0.0066). The slope used here, 0.0121, is measured along the TRACKED
     thalweg (27.4 ft over 2,264 ft) and is the value carried in the report.

  Together those two errors put the tailwater at +1.53 ft above the inlet invert
  instead of +7.42 ft, which reverses the drowned-outlet finding.
"""
import sys, math
import numpy as np
from osgeo import gdal
gdal.UseExceptions()

FT, G = 3.280839895, 32.174
DEM = "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-141Local/demhr/cliped_utm.tif"
XC, Y_OUT = 711237.0, 3960500.0      # crossing easting; outlet section northing
HALF_M, DX_M, NPT = 250.0, 2.0, 250  # section extent and sampling

S_DS  = 0.0121      # ft/ft, tracked thalweg, 27.4 ft over 2,264 ft
N_CH  = 0.045       # channel Manning n (NLCD-derived range 0.035-0.050)
INV_IN, INV_OUT, ROAD_LOW = 6112.80, 6110.60, 6120.60
FREEBOARD = 1.0
HW_MAX = ROAD_LOW - INV_IN - FREEBOARD          # 6.80 ft above inlet invert
TW_2D  = 6119.89                                 # 2-D wse_max, corroboration only

Q = float(sys.argv[1]) if len(sys.argv) > 1 else 11922.0


def outlet_section():
    """Bed elevations (ft) across the outlet section, straight off the 1 m LiDAR."""
    ds = gdal.Open(DEM); GT = ds.GetGeoTransform(); band = ds.GetRasterBand(1)
    j = int((GT[3] - Y_OUT) / abs(GT[5]))
    xs = [XC - HALF_M + DX_M * i for i in range(NPT)]
    i0 = int((min(xs) - GT[0]) / GT[1]); i1 = int((max(xs) - GT[0]) / GT[1]) + 2
    row = band.ReadAsArray(i0, j, i1 - i0, 1).astype(float)[0]
    return np.array([row[int((x - GT[0]) / GT[1]) - i0] * FT for x in xs])


Z = outlet_section()
BED = float(np.nanmin(Z))
DX_FT = DX_M * FT


def geom(d):
    """Area (sf), top width (ft), wetted perimeter (ft), hydraulic radius (ft)."""
    wse = BED + d
    wet = Z < wse
    if not wet.any():
        return 0.0, 0.0, 0.0, 0.0
    A = float(np.sum(wse - Z[wet])) * DX_FT
    T = float(np.sum(wet)) * DX_FT
    P = T + 2.0 * d
    return A, T, P, (A / P if P > 0 else 0.0)


def bisect(f, lo=0.05, hi=40.0, it=100):
    for _ in range(it):
        m = 0.5 * (lo + hi)
        if f(m) < 0: lo = m
        else:        hi = m
    return 0.5 * (lo + hi)


def critical_depth(Q):
    """Q^2 T / (g A^3) = 1."""
    return bisect(lambda d: (lambda A, T, P, R: 1.0 - Q * Q * T / (G * A ** 3) if A > 0 else -1.0)(*geom(d)))


def normal_depth(Q, n=N_CH, S=S_DS):
    return bisect(lambda d: (lambda A, T, P, R: (1.49 / n) * A * R ** (2.0 / 3.0) * math.sqrt(S) - Q)(*geom(d)))


print("Outlet section: %d points at %.1f m over +/-%.0f m, y = %.0f" % (NPT, DX_M, HALF_M, Y_OUT))
print("  invert %.2f ft   bed slope S_ds = %.4f (tracked thalweg)   n = %.3f\n" % (BED, S_DS, N_CH))
print("Design discharge Q = %,.0f cfs".replace("%,", "%") % Q if False else "Design discharge Q = {:,.0f} cfs\n".format(Q))

dc = critical_depth(Q); Ac, Tc, _, _ = geom(dc)
print("critical depth : d_c = %5.2f ft   A = %6.0f sf  T = %5.0f ft   WSE = %8.2f ft" % (dc, Ac, Tc, BED + dc))
for n in (0.035, N_CH, 0.050):
    dn = normal_depth(Q, n=n)
    tag = "  <- carried" if abs(n - N_CH) < 1e-9 else ""
    print("normal depth   : n = %.3f  d_n = %5.2f ft            WSE = %8.2f ft%s" % (n, dn, BED + dn, tag))

dn = normal_depth(Q)
wse_c, wse_n = BED + dc, BED + dn
gov, which = (wse_n, "normal") if dn > dc else (wse_c, "critical")
ho_in, ho_out = gov - INV_IN, gov - INV_OUT

print("\n--- comparison at the crossing ---")
print("  inlet invert %8.2f ft    outlet invert %8.2f ft    road low %8.2f ft" % (INV_IN, INV_OUT, ROAD_LOW))
print("  governing tailwater (%s depth) %8.2f ft" % (which, gov))
print("     = %+6.2f ft above the inlet invert   (h_o = %+.2f ft above the outlet invert)" % (ho_in, ho_out))
print("  allowable headwater HW_max          %+6.2f ft" % HW_MAX)
if ho_in > HW_MAX:
    print("  >> DROWNED OUTLET: the tailwater alone exceeds HW_max by %.2f ft." % (ho_in - HW_MAX))
    print("     No opening of any size satisfies the headwater limit.")
else:
    print("  >> NOT drowned: %.2f ft of head remains to drive the barrel," % (HW_MAX - ho_in))
    print("     so outlet control governs the sizing.")
dch = HW_MAX / 1.5
q_unit = math.sqrt(G * dch ** 3)
print("\n  inlet-control floor at HW_max: d_c = %.2f ft, q = %.1f cfs/ft -> opening >= %,.0f ft"
      .replace("%,", "%") % (dch, q_unit, Q / q_unit))
print("\n  corroboration: 2-D wse_max downstream of the road = %.2f ft (%+.2f ft vs critical)"
      % (TW_2D, TW_2D - wse_c))
