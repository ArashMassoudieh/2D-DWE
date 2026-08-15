#!/usr/bin/env python3
"""
carve_culvert_notch.py -- represent a culvert in 2D terrain as a bounded-width
notch whose CONVEYANCE matches the barrel, not as an open trench.

    python3 tools/carve_culvert_notch.py DEM_IN LINE DEM_OUT --diameter-m 2.0

WHY NOT REUSE carve_dem_by_line.py
  That tool lowers every cell the line passes through, producing a full-width
  ramp. That is correct for D8 delineation, where the only requirement is
  hydrologic connectivity so the upstream area is not lost. It is wrong for 2D
  hydraulics, because it deletes the road embankment as a control: the model
  then conveys the entire flood through a trench and never ponds or overtops,
  which at this site is the actual flood mechanism.

WHAT THIS DOES INSTEAD
  1. Finds the true bed elevation at each end of the culvert line by taking the
     DEM minimum within a search radius, rather than sampling the line's vertex.
     A digitised endpoint often lands on the bank or the fill; at NM15-141 one
     of them sat 3 m up the embankment.
  2. Interpolates an invert linearly between those two bed elevations.
  3. Lowers only cells within HALF-WIDTH of the centreline, where the half-width
     is solved so the notch's Manning conveyance equals the culvert's capacity.
  4. Never raises ground. A cell already below the invert profile is left alone.

SIZING
  Culvert capacity uses the orifice form Q = Cd*A*sqrt(2 g H) with Cd = 0.6 and
  H measured to the barrel centroid -- appropriate for a submerged inlet, which
  is the governing condition when the design flood far exceeds the barrel.

  Notch conveyance uses Manning with the true hydraulic radius A/P. Using R ~ y
  for a narrow deep notch overestimates conveyance several-fold and will size
  the notch far too small.

  The match is at ONE head. A carved notch is an open channel: above the barrel
  crown its capacity grows as H^(5/3) while a real culvert's grows as sqrt(H),
  so the notch over-conveys at higher stages. Where the design flood greatly
  exceeds the barrel -- the usual case for arid-wash crossings -- this matters
  little, because the road overtops and the embankment, not the barrel, controls.
"""
import argparse
import math
import numpy as np
from osgeo import gdal, ogr, osr

gdal.UseExceptions(); ogr.UseExceptions()
FT = 3.280839895
G_FPS2 = 32.174


def culvert_capacity_cfs(D_m, head_ft, Cd=0.6):
    """Submerged-inlet orifice capacity of a circular barrel, cfs."""
    D = D_m * FT
    A = math.pi * D * D / 4.0
    H = head_ft - D / 2.0
    if H <= 0:
        return 0.0, A
    return Cd * A * math.sqrt(2.0 * G_FPS2 * H), A


def notch_Q_cfs(w_m, y_ft, n, S):
    """Manning conveyance of a rectangular notch, cfs (US units)."""
    w = w_m * FT
    A = w * y_ft
    P = w + 2.0 * y_ft
    R = A / P if P > 0 else 0.0
    return (1.486 / n) * A * R ** (2.0 / 3.0) * math.sqrt(S)


def solve_width(Q_target, y_ft, n, S):
    lo, hi = 0.01, 200.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if notch_Q_cfs(mid, y_ft, n, S) < Q_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dem_in")
    ap.add_argument("line")
    ap.add_argument("dem_out")
    ap.add_argument("--diameter-m", type=float, default=None,
                    help="size the notch to match this barrel diameter")
    ap.add_argument("--n", type=float, default=0.030)
    ap.add_argument("--search-r", type=float, default=6.0,
                    help="radius (m) to search for the true bed at each end")
    ap.add_argument("--width-m", type=float, default=None,
                    help="override the solved notch width")
    A = ap.parse_args()

    ds = gdal.Open(A.dem_in)
    b = ds.GetRasterBand(1)
    gt = ds.GetGeoTransform()
    nx, ny = ds.RasterXSize, ds.RasterYSize
    px = abs(gt[1])
    z = b.ReadAsArray().astype("float64")
    nodata = b.GetNoDataValue()

    def rc(x, y):
        return int((x - gt[0]) / gt[1]), int((gt[3] - y) / abs(gt[5]))

    lds = ogr.Open(A.line); lly = lds.GetLayer()
    lsr = lly.GetSpatialRef(); dsr = osr.SpatialReference(wkt=ds.GetProjection())
    feat = lly.GetNextFeature()
    geom = feat.GetGeometryRef().Clone()
    if lsr and dsr and not lsr.IsSame(dsr):
        lsr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        dsr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        geom.Transform(osr.CoordinateTransformation(lsr, dsr))
    pts = [(p[0], p[1]) for p in geom.GetPoints()]
    P0, P1 = pts[0], pts[-1]

    # --- true bed at each end: DEM minimum within search radius --------------
    def bed_min(P):
        i0, j0 = rc(P[0], P[1])
        r = int(math.ceil(A.search_r / px))
        sub = z[max(0, j0 - r):j0 + r + 1, max(0, i0 - r):i0 + r + 1]
        if nodata is not None:
            sub = sub[sub != nodata]
        return float(np.nanmin(sub))

    zA, zB = bed_min(P0), bed_min(P1)
    L = math.hypot(P1[0] - P0[0], P1[1] - P0[1])

    # --- road crest along the line, for the available head -------------------
    ux, uy = (P1[0] - P0[0]) / L, (P1[1] - P0[1]) / L
    crest = -1e30
    for s in np.arange(0, L + 0.01, px / 2):
        i, j = rc(P0[0] + ux * s, P0[1] + uy * s)
        if 0 <= i < nx and 0 <= j < ny:
            crest = max(crest, z[j, i])
    inv_hi = max(zA, zB)
    head_ft = (crest - inv_hi) * FT
    S = abs(zA - zB) / L if L > 0 else 0.001
    S = max(S, 1e-4)

    y_ft = head_ft
    if A.width_m is None and A.diameter_m is None:
        raise SystemExit("give --width-m or --diameter-m")

    if A.width_m is not None:
        # Width given: report the barrel it is equivalent to, by inverting the
        # capacity relation. This is what makes a width sweep interpretable --
        # a notch width on its own says nothing about the structure it stands for.
        w = A.width_m
        Qn = notch_Q_cfs(w, y_ft, A.n, S)
        # Capacity is NOT monotonic in diameter: once D/2 exceeds the available
        # head the orifice head goes negative and Q collapses to zero, so a plain
        # bisection walks off the end and returns the search bound. Scan the
        # physically admissible range D < 2*head and take the SMALLEST diameter
        # that reaches the notch conveyance; report the peak if none does.
        Dmax_ft = 2.0 * head_ft
        grid = np.linspace(0.05, Dmax_ft / FT * 0.999, 4000)
        caps = np.array([culvert_capacity_cfs(dd, head_ft)[0] for dd in grid])
        Qpeak = float(caps.max()); Dpeak = float(grid[int(caps.argmax())])
        idx = np.where(caps >= Qn)[0]
        if idx.size:
            D_eq = float(grid[idx[0]])
            equiv = ("  notch conveys %.0f cfs == a %.2f m (%.1f ft) barrel"
                     % (Qn, D_eq, D_eq * FT))
        else:
            D_eq = Dpeak
            equiv = ("  notch conveys %.0f cfs -- MORE than any barrel fits under this head\n"
                     "  (peak capacity %.0f cfs at %.2f m dia; the notch over-represents the culvert)"
                     % (Qn, Qpeak, Dpeak))
        Qc, Abar = Qn, math.pi * (D_eq * FT) ** 2 / 4.0
    else:
        Qc, Abar = culvert_capacity_cfs(A.diameter_m, head_ft)
        w = solve_width(Qc, y_ft, A.n, S)
        equiv = "  culvert %.2f m dia, area %.1f sf -> capacity %.0f cfs" % (A.diameter_m, Abar, Qc)

    print("CULVERT -> NOTCH SIZING")
    print("  line length              %8.1f m" % L)
    print("  bed at end A / end B     %8.2f / %.2f m" % (zA, zB))
    print("  road crest on the line   %8.2f m" % crest)
    print("  head above higher invert %8.2f m (%.1f ft)" % (crest - inv_hi, head_ft))
    print("  barrel slope             %8.3f%%" % (100 * S))
    print(equiv)
    print("  notch width              %8.2f m  (%.1f cells at %.1f m)"
          % (w, w / px, px))

    # --- carve ---------------------------------------------------------------
    half = w / 2.0
    step = px / 3.0
    nmod = 0
    zmin_carved = 1e30
    rad = int(math.ceil(half / px)) + 1
    for s in np.arange(0, L + 0.01, step):
        x = P0[0] + ux * s; y = P0[1] + uy * s
        inv = zA + (zB - zA) * (s / L)
        i0, j0 = rc(x, y)
        for jj in range(j0 - rad, j0 + rad + 1):
            for ii in range(i0 - rad, i0 + rad + 1):
                if not (0 <= ii < nx and 0 <= jj < ny):
                    continue
                cx = gt[0] + (ii + 0.5) * gt[1]
                cy = gt[3] - (jj + 0.5) * abs(gt[5])
                # perpendicular distance from the centreline
                d = abs((cx - x) * uy - (cy - y) * ux)
                if d > half:
                    continue
                if z[jj, ii] > inv:
                    z[jj, ii] = inv
                    nmod += 1
                    zmin_carved = min(zmin_carved, inv)

    drv = gdal.GetDriverByName("GTiff")
    out = drv.Create(A.dem_out, nx, ny, 1, gdal.GDT_Float32,
                     ["COMPRESS=DEFLATE", "TILED=YES"])
    out.SetGeoTransform(gt); out.SetProjection(ds.GetProjection())
    ob = out.GetRasterBand(1)
    if nodata is not None:
        ob.SetNoDataValue(nodata)
    ob.WriteArray(z)
    out.FlushCache(); out = None

    print("  cells lowered            %8d" % nmod)
    print("  wrote %s" % A.dem_out)


if __name__ == "__main__":
    main()
