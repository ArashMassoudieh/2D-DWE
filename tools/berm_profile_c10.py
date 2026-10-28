#!/usr/bin/env python3
"""
Verify the existing berm around the C10 sub-community, measure its height, and
locate where it has deteriorated.

    python3 tools/berm_profile_c10.py [--offset-ft 5] [--fig out.png]

WHAT IS MEASURED
  Three parallel profiles are cut along each digitised line: the line itself and
  the same line offset +/- the --offset-ft distance perpendicular to it. A berm
  is a built ridge, so all three should ride high together; if only the centre
  profile is high the "berm" is a digitising artefact or a one-cell spike, and
  if all three sit at ground level there is no berm there at all. That is the
  confirmation the offset profiles provide.

  Height is measured as the crest ABOVE THE LOCAL GRADE, obtained by detrending
  each transect: a straight line is drawn between the two toe elevations at
  -/+ TOE_FT and the crest is measured above that line. The detrend is what
  makes the number mean anything here. Line 2 runs along a hillside, and taking
  the ground as the higher toe -- the obvious first choice -- reported the north
  berm as 0.17 ft BELOW ground and found no berm at all, because the uphill toe
  sits above the crest. Against the local grade the same berm reads as a real,
  if low, ridge. The crest is located as the highest detrended point within
  +/- CREST_SEARCH_FT of the digitised line rather than assumed to lie on it.

DETERIORATION is a crest that falls away from its own neighbours, not one that
is merely low: the crest is compared with a running median along the line, and a
run of stations sitting more than --drop-ft below it is reported as a breach.
The entrance road crosses the north line and shows up the same way, so road
crossings are located from TIGER and labelled separately -- an opening that is a
road is not a defect.
"""
import argparse, os
import numpy as np
from osgeo import gdal, ogr, osr

gdal.UseExceptions(); ogr.UseExceptions()
FT = 3.280839895
NHA = "/home/arash/Dropbox/Chloeta/NHA"
SITE = NHA + "/WS3_GIS/NM15-141Local"

ap = argparse.ArgumentParser()
ap.add_argument("--berm", default=SITE + "/outputs/existing_berm.shp")
ap.add_argument("--dem",  default=SITE + "/demhr/cliped_sd.tif")
ap.add_argument("--offset-ft", type=float, default=5.0)
ap.add_argument("--toe-ft", type=float, default=39.0, help="half-width of the detrend transect")
ap.add_argument("--crest-search-ft", type=float, default=16.0)
ap.add_argument("--step", type=float, default=1.0, help="station spacing (m)")
ap.add_argument("--drop-ft", type=float, default=1.0)
ap.add_argument("--min-run", type=float, default=4.0, help="m; shorter dips are noise")
ap.add_argument("--smooth", type=float, default=41.0, help="running-median window (m)")
ap.add_argument("--fig", default=None)
ap.add_argument("--csv", default=None)
A = ap.parse_args()

_ds = gdal.Open(A.dem)
GT = _ds.GetGeoTransform()
Z = _ds.GetRasterBand(1).ReadAsArray().astype(float)
_nd = _ds.GetRasterBand(1).GetNoDataValue()
if _nd is not None:
    Z = np.where(Z == _nd, np.nan, Z)
NY, NX = Z.shape


def samp(x, y):
    i = np.clip(((np.asarray(x) - GT[0]) / GT[1]).astype(int), 0, NX - 1)
    j = np.clip(((GT[3] - np.asarray(y)) / abs(GT[5])).astype(int), 0, NY - 1)
    return Z[j, i]


def runmed(v, win):
    k = max(1, int(win) // 2)
    return np.array([np.nanmedian(v[max(0, i - k):i + k + 1]) for i in range(len(v))])


def lines():
    ds = ogr.Open(A.berm); ly = ds.GetLayer()
    out = []
    for f in ly:
        g = f.GetGeometryRef()
        p = np.array(g.GetPoints())[:, :2]
        out.append((f.GetField("id"), p))
    return out


def road_crossings(px, py):
    """Stations where a TIGER road passes within 12 m of the line."""
    ds = ogr.Open(SITE + "/roads/tl_2025_35031_roads/tl_2025_35031_roads.shp")
    ly = ds.GetLayer(); sr = ly.GetSpatialRef()
    tg = osr.SpatialReference(); tg.ImportFromEPSG(26912)
    sr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tg.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = None if sr.IsSame(tg) else osr.CoordinateTransformation(sr, tg)
    hit = np.zeros(len(px), bool)
    for f in ly:
        g = f.GetGeometryRef()
        if g is None:
            continue
        g = g.Clone()
        if tr:
            g.Transform(tr)
        e = g.GetEnvelope()
        if e[1] < px.min() - 50 or e[0] > px.max() + 50 or e[3] < py.min() - 50 or e[2] > py.max() + 50:
            continue
        for k in range(len(px)):
            pt = ogr.Geometry(ogr.wkbPoint); pt.AddPoint(float(px[k]), float(py[k]))
            if g.Distance(pt) < 12.0:
                hit[k] = True
    return hit


def runs(mask, s, min_len):
    out = []; i = 0
    while i < len(mask):
        if mask[i]:
            j = i
            while j + 1 < len(mask) and mask[j + 1]:
                j += 1
            if s[j] - s[i] >= min_len:
                out.append((i, j))
            i = j + 1
        else:
            i += 1
    return out


NAMES = {0: "Line 1 - west boundary", 1: "Line 2 - north boundary"}
off = A.offset_ft / FT

results = []

for idx, (fid, pts) in enumerate(lines()):
    name = NAMES.get(idx, "Line %d" % (idx + 1))
    seg = np.r_[0, np.cumsum(np.hypot(*np.diff(pts, axis=0).T))]
    s = np.arange(0, seg[-1], A.step)
    px = np.interp(s, seg, pts[:, 0]); py = np.interp(s, seg, pts[:, 1])
    dx = np.gradient(px); dy = np.gradient(py)
    L = np.hypot(dx, dy); ux, uy = -dy / L, dx / L      # left normal

    zc = samp(px, py)
    zl = samp(px + off * ux, py + off * uy)
    zr = samp(px - off * ux, py - off * uy)

    # detrended transect: straight line between the two toes, crest above it
    T = A.toe_ft / FT
    offs = np.arange(-T, T + 1e-9, 1.0)
    grid = np.array([samp(px + o * ux, py + o * uy) for o in offs])   # (noff, nstn)
    zlo, zhi = grid[0], grid[-1]
    trend = zlo[None, :] + (zhi - zlo)[None, :] * ((offs[:, None] + T) / (2 * T))
    resid = grid - trend
    win = np.abs(offs) <= A.crest_search_ft / FT
    h = np.nanmax(resid[win], axis=0) * FT                 # berm height above local grade
    kbest = np.argmax(np.where(win[:, None], resid, -np.inf), axis=0)
    crest_off = offs[kbest]
    crest = grid[kbest, np.arange(grid.shape[1])]
    ground = crest - h / FT
    # the two flanking profiles, also detrended, for the confirmation test
    hl = np.nanmax(resid[win & (offs >= 0)], axis=0) * FT
    hr = np.nanmax(resid[win & (offs <= 0)], axis=0) * FT

    base = runmed(h, A.smooth)
    road = road_crossings(px, py)
    deter = (h < base - A.drop_ft) & ~road

    print("=" * 78)
    print("%s   %.0f m (%.0f ft), %d stations at %.1f m" % (name, s[-1], s[-1] * FT, len(s), A.step))
    print("  crest height above local ground (ft): mean %.2f  median %.2f  min %.2f  max %.2f"
          % (np.nanmean(h), np.nanmedian(h), np.nanmin(h), np.nanmax(h)))
    print("  offset profiles at +/-%.0f ft (detrended): centre %.2f | +%.0f ft %.2f | -%.0f ft %.2f"
          % (A.offset_ft, np.nanmean((zc - ground) * FT),
             A.offset_ft, np.nanmean((zl - ground) * FT),
             A.offset_ft, np.nanmean((zr - ground) * FT)))
    print("  crest offset from the digitised line: median %+.1f m, |offset| < 3 m at %.0f%% of stations"
          % (np.median(crest_off), 100 * np.mean(np.abs(crest_off) < 3.0)))
    print("  both flanks raised (a real ridge, not a slope break) at %.0f%% of stations"
          % (100 * np.mean((hl > 0.3) & (hr > 0.3))))
    print("  stations reading as a raised ridge (>1 ft): %.0f%%" % (100 * np.mean(h > 1.0)))
    if road.any():
        rr = runs(road, s, 0.0)
        print("  road crossings (excluded from defects): %s"
              % ", ".join("s=%.0f-%.0f m" % (s[i], s[j]) for i, j in rr))
    dr = runs(deter, s, A.min_run)
    if dr:
        print("  DETERIORATED / LOW SEGMENTS (crest >%.1f ft below its own running median):" % A.drop_ft)
        for i, j in dr:
            print("     s=%6.0f-%6.0f m (%5.0f ft long)  crest height %.2f ft vs %.2f ft expected"
                  " -- deficit %.2f ft  @ (%.0f, %.0f)"
                  % (s[i], s[j], (s[j] - s[i]) * FT, np.nanmin(h[i:j + 1]),
                     np.nanmean(base[i:j + 1]), np.nanmean(base[i:j + 1]) - np.nanmin(h[i:j + 1]),
                     px[(i + j) // 2], py[(i + j) // 2]))
    else:
        print("  no segment falls more than %.1f ft below its running median" % A.drop_ft)
    results.append(dict(name=name, s=s, px=px, py=py, h=h, base=base, zc=zc, zl=zl, zr=zr,
                        ground=ground, road=road, deter=deter, runs=dr))

if A.csv:
    with open(A.csv, "w") as fh:
        fh.write("line,station_m,east,north,crest_ft,ground_ft,height_ft,expected_ft,road,deteriorated\n")
        for r in results:
            for k in range(len(r["s"])):
                fh.write("%s,%.1f,%.1f,%.1f,%.2f,%.2f,%.2f,%.2f,%d,%d\n" %
                         (r["name"], r["s"][k], r["px"][k], r["py"][k],
                          np.nanmax([r["zc"][k], r["zl"][k], r["zr"][k]]) * FT,
                          r["ground"][k] * FT, r["h"][k], r["base"][k],
                          int(r["road"][k]), int(r["deter"][k])))
    print("\nwrote " + A.csv)

if A.fig:
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource
    fig = plt.figure(figsize=(15, 11), dpi=150)
    gs = fig.add_gridspec(3, 1, height_ratios=[1.45, 1, 1], hspace=0.33)
    ax = fig.add_subplot(gs[0])
    hs = LightSource(azdeg=315, altdeg=45).hillshade(np.nan_to_num(Z, nan=np.nanmin(Z)),
                                                     vert_exag=3, dx=1, dy=1)
    ext = [GT[0], GT[0] + NX * GT[1], GT[3] + NY * GT[5], GT[3]]
    ax.imshow(hs, cmap="gray", extent=ext, origin="upper", zorder=1)
    for r, col in zip(results, ("#1a4f95", "#7a1020")):
        ax.plot(r["px"], r["py"], c=col, lw=2.2, zorder=4, label=r["name"])
        m = r["deter"]
        ax.scatter(r["px"][m], r["py"][m], s=16, c="#e8112d", zorder=6)
        m2 = r["road"]
        ax.scatter(r["px"][m2], r["py"][m2], s=16, c="#00b050", zorder=6)
    ax.scatter([], [], s=16, c="#e8112d", label="Deteriorated / low")
    ax.scatter([], [], s=16, c="#00b050", label="Road crossing")
    ax.set_xlim(min(r["px"].min() for r in results) - 80, max(r["px"].max() for r in results) + 80)
    ax.set_ylim(min(r["py"].min() for r in results) - 80, max(r["py"].max() for r in results) + 80)
    ax.set_aspect("equal"); ax.legend(loc="lower right", fontsize=9, framealpha=.93)
    ax.set_title("C10 — existing berm: alignment, crest condition and defects", fontsize=12.5)
    ax.set_xlabel("UTM Zone 12N Easting (m)"); ax.set_ylabel("Northing (m)")
    for k, r in enumerate(results):
        a2 = fig.add_subplot(gs[k + 1])
        sf = r["s"] * FT
        a2.plot(sf, r["h"], c="#1a4f95", lw=1.5, label="Crest height above local ground")
        a2.plot(sf, r["base"], c="#888", lw=1.2, ls="--", label="Running median (expected crest)")
        a2.fill_between(sf, r["h"], r["base"], where=r["deter"], color="#e8112d", alpha=.5,
                        label="Deteriorated / low")
        if r["road"].any():
            a2.fill_between(sf, 0, np.nanmax(r["h"]), where=r["road"], color="#00b050", alpha=.18,
                            label="Road crossing")
        a2.axhline(0, c="k", lw=.8)
        a2.set_xlim(0, sf[-1]); a2.set_xlabel("Distance along %s (ft)" % r["name"])
        a2.set_ylabel("Height (ft)"); a2.grid(alpha=.3)
        a2.legend(fontsize=8, loc="upper right", ncol=2)
        a2.set_title(r["name"], fontsize=10.5)
    fig.savefig(A.fig, bbox_inches="tight")
    print("wrote " + A.fig)
