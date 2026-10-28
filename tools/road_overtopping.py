#!/usr/bin/env python3
"""
Road-crest overtopping for C10 (NM15-141): maximum head over the Navajo Route 9
crest and the fraction of the crest overtopped, for each 2-D run.

    python3 tools/road_overtopping.py                      # all four post-C10 runs
    python3 tools/road_overtopping.py --qc-figure out.png

METHOD -- the C6 divide calculation, transplanted
  C6 sampled a HAND-DIGITISED barrier polyline (natural_barrier.shp), densified
  at 5 m, and reported max(depth on the line) and the share of stations with
  more than 0.5 ft on them.  The same two statistics are produced here and the
  0.5 ft threshold is carried over unchanged.

  C10 HAS NO DIGITISED ROAD-CREST LINE, so one is derived.  That derivation is
  the only methodological difference from C6 and it is the main uncertainty in
  the answer -- see CREST DERIVATION.  A digitised crest would replace it
  directly: pass --crest-line and the derivation is skipped.

CREST DERIVATION
  The TIGER/Line centreline for Navajo Route 9 is a cartographic centreline, not
  a survey: over this 1,102 m it lies a median 12 m off the embankment crown and
  0.82 m BELOW it (p90 1.90 m, max 3.38 m).  Sampling it directly would
  understate the crest by about that much and overstate every head by the same
  amount -- on the order of the reported number itself, so it cannot be ignored.

  At each 1 m station the crest is taken as the highest ground within +/-20 m
  perpendicular to the centreline.  20 m is chosen because it is the widest
  window at which the search almost never runs to its own edge (3 % of stations,
  against 18 % at 15 m and 18 % at 25 m); running to the edge means the window,
  not the terrain, is choosing the crest.  The resulting offsets are then
  running-median filtered (--smooth, default 101 m) because over two stretches
  -- s = 200-300 m and s > 1000 m -- a natural high south of the road outranks
  the embankment and the raw snap jumps off the road; the filter holds the line
  on the embankment through those stretches.  Both effects are reported by
  --sensitivity rather than asserted.

WHY HEAD IS WSE MINUS CREST, NOT depth_max
  C6 read depth off the line, which works because C6's terrain has no structure
  cut into it.  Here the culvert-8 notch is a 3 m x 41 m north-south strip and
  the road-perpendicular transect runs ALONG it, so on the carved runs a depth
  sample at the crossing can be depth over the notch invert -- water standing in
  the opening, not over the road.  Notch cells (z_raw - z_carved > 0.05 m) are
  therefore excluded from crest candidacy, the crest elevation is taken from the
  RAW DEM so it is identical for all four runs, and head is WSE - crest.  Open
  and blocked are then differenced against one common crest.
"""
import argparse, os, sys
import numpy as np
from osgeo import gdal, ogr, osr

gdal.UseExceptions(); ogr.UseExceptions()
FT = 3.280839895
SITE = "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-141Local"
RES  = SITE + "/output_2D/2D_Results"

RUNS = [("6hr",  "open",    "out_nm141_1m_6hr_carved_2m_postC10"),
        ("6hr",  "blocked", "out_nm141_1m_6hr_raw_postC10"),
        ("24hr", "open",    "out_nm141_1m_24hr_carved_2m_postC10"),
        ("24hr", "blocked", "out_nm141_1m_24hr_raw_postC10")]

ap = argparse.ArgumentParser()
ap.add_argument("--raw-dem",    default=SITE + "/output_2D/dem_2d_raw.tif")
ap.add_argument("--carved-dem", default=SITE + "/output_2D/dem_2d_carved_2m.tif")
ap.add_argument("--roads",      default=SITE + "/roads/tl_2025_35031_roads/tl_2025_35031_roads.shp")
ap.add_argument("--linearid",   default="1106088001728", help="TIGER LINEARID of Navajo Route 9")
ap.add_argument("--crest-line", default=None, help="digitised crest polyline; skips the derivation")
ap.add_argument("--window",    type=float, default=20.0, help="crest search half-width (m)")
ap.add_argument("--smooth",    type=float, default=101.0, help="running-median window on the offset (m)")
ap.add_argument("--threshold", type=float, default=0.5,  help="overtopping threshold (ft), C6 used 0.5")
ap.add_argument("--step",      type=float, default=1.0,  help="station spacing (m)")
ap.add_argument("--qc-figure", default=None)
ap.add_argument("--sensitivity", action="store_true")
ap.add_argument("--csv", default=None)
A = ap.parse_args()

# ------------------------------------------------------------------ rasters --
_keep = []
def read(path):
    ds = gdal.Open(path); _keep.append(ds)
    b = ds.GetRasterBand(1)
    a = b.ReadAsArray().astype(float)
    nd = b.GetNoDataValue()
    if nd is not None:
        a = np.where(a == nd, np.nan, a)
    return a, ds.GetGeoTransform()

zraw, GT = read(A.raw_dem)
zcar, _  = read(A.carved_dem)
NY, NX = zraw.shape
EXT = [GT[0], GT[0] + NX * GT[1], GT[3] + NY * GT[5], GT[3]]
NOTCH = (zraw - zcar) > 0.05

def samp(a, x, y):
    i = np.clip(((np.asarray(x) - GT[0]) / GT[1]).astype(int), 0, NX - 1)
    j = np.clip(((GT[3] - np.asarray(y)) / abs(GT[5])).astype(int), 0, NY - 1)
    return a[j, i]

# ----------------------------------------------------------------- alignment --
def domain_poly():
    r = ogr.Geometry(ogr.wkbLinearRing)
    for x, y in [(EXT[0], EXT[2]), (EXT[1], EXT[2]), (EXT[1], EXT[3]),
                 (EXT[0], EXT[3]), (EXT[0], EXT[2])]:
        r.AddPoint(x, y)
    p = ogr.Geometry(ogr.wkbPolygon); p.AddGeometry(r); return p

def load_line(path, linearid=None):
    ds = ogr.Open(path); _keep.append(ds); ly = ds.GetLayer()
    sr = ly.GetSpatialRef(); tg = osr.SpatialReference(); tg.ImportFromEPSG(26912)
    sr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tg.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = None if sr.IsSame(tg) else osr.CoordinateTransformation(sr, tg)
    dp = domain_poly(); best = None
    for f in ly:
        if linearid is not None and f.GetField("LINEARID") != linearid:
            continue
        g = f.GetGeometryRef()
        if g is None: continue
        g = g.Clone()
        if tr: g.Transform(tr)
        c = g.Intersection(dp)
        if c is None or c.IsEmpty(): continue
        parts = [c] if c.GetGeometryCount() == 0 else \
                [c.GetGeometryRef(k) for k in range(c.GetGeometryCount())]
        parts = [q for q in parts if q.GetPointCount() > 1]
        if not parts: continue
        q = max(parts, key=lambda t: t.Length())
        if best is None or q.Length() > best.Length(): best = q.Clone()
    if best is None: sys.exit("no line found in " + path)
    return np.array(best.GetPoints())[:, :2]

def stations(pts, step):
    if pts[0, 0] > pts[-1, 0]:          # west -> east
        pts = pts[::-1]
    seg = np.r_[0, np.cumsum(np.hypot(*np.diff(pts, axis=0).T))]
    s = np.arange(0, seg[-1], step)
    return s, np.interp(s, seg, pts[:, 0]), np.interp(s, seg, pts[:, 1])

def runmed(v, win):
    if win <= 1: return v
    k = max(1, int(win) // 2)
    out = np.empty_like(v)
    for i in range(len(v)):
        out[i] = np.median(v[max(0, i - k):i + k + 1])
    return out

def derive_crest(window, smooth):
    """Returns (s, x, y, z, offset, edge-hit fraction, valid-station mask).

    The search runs on the RAW DEM, which has no notch in it, so the road crown
    is found at every station including the crossing.  The notch is handled
    afterwards, on the finished line: the three stations where the culvert-8
    footprint meets the crest are an OPENING, not road, and are dropped from the
    statistics.  Excluding notch cells from the search instead -- the obvious
    first move -- punches a hole in the transect at exactly the crossing, and
    the snap then falls 20 m off the embankment onto the toe and reports the
    head over the toe as if it were head over the crest (4.17 ft against
    1.56 ft, the single largest error found in building this).
    """
    pts = load_line(A.roads, A.linearid)
    s, cx, cy = stations(pts, A.step)
    dx = np.gradient(cx); dy = np.gradient(cy)
    L = np.hypot(dx, dy); ux, uy = -dy / L, dx / L      # left normal
    offs = np.arange(-window, window + 1e-9, 1.0)
    zg = np.array([samp(zraw, cx + o * ux, cy + o * uy) for o in offs])
    raw_ob = offs[np.nanargmax(zg, axis=0)]
    ob = runmed(raw_ob, smooth)
    x = cx + ob * ux; y = cy + ob * uy
    valid = ~samp(NOTCH, x, y).astype(bool)
    return (s, x, y, samp(zraw, x, y), ob,
            float(np.mean(np.abs(raw_ob) >= window - 0.5)), valid)

if A.crest_line:
    pts = load_line(A.crest_line)
    S, X, Y = stations(pts, A.step)
    Z = samp(zraw, X, Y); OB = np.zeros_like(S); EDGE = 0.0
    VALID = ~samp(NOTCH, X, Y).astype(bool)
    SRC = "digitised: " + os.path.basename(A.crest_line)
else:
    S, X, Y, Z, OB, EDGE, VALID = derive_crest(A.window, A.smooth)
    SRC = "derived from TIGER %s, +/-%.0f m snap, %.0f m median filter" % (A.linearid, A.window, A.smooth)

# ------------------------------------------------------------------ analysis --
def head_ft(run):
    w, _ = read(os.path.join(RES, run, "wse_max.tif"))
    h = (samp(w, X, Y) - Z) * FT
    return np.where(np.isfinite(h), h, 0.0)

def stats(h, valid=None):
    v = VALID if valid is None else valid
    hv = h[v]
    ov = hv > A.threshold
    imax = int(np.flatnonzero(v)[int(np.argmax(hv))])
    return dict(maxhead=float(hv.max()), pct=100.0 * ov.mean(),
                length_ft=float(v.sum() * A.step * FT), ov_ft=float(ov.sum() * A.step * FT),
                s_at_max=float(S[imax]), x=float(X[imax]), y=float(Y[imax]))

print("C10 -- Navajo Route 9 road-crest overtopping, post-culvert-10 runs")
print("crest: %s" % SRC)
print("       %.0f m (%,.0f ft) inside the domain, %d stations at %.1f m"
      .replace("%,.0f", "%.0f") % (S[-1], S[-1] * FT, len(S), A.step))
print("       crest elevation %.2f - %.2f m, low point %.2f m; window-edge hits %.0f%%"
      % (Z.min(), Z.max(), Z.min(), 100 * EDGE))
print("threshold for 'overtopped': head > %.1f ft (C6 value)\n" % A.threshold)
hdr = "%-6s %-8s %10s %10s %12s %10s" % ("storm", "culvert-8", "max head", "% crest", "length ovt", "at station")
print(hdr); print("-" * len(hdr))
rows = []
for storm, cond, run in RUNS:
    p = os.path.join(RES, run, "wse_max.tif")
    if not os.path.exists(p):
        print("%-6s %-8s   MISSING %s" % (storm, cond, p)); continue
    st = stats(head_ft(run))
    rows.append((storm, cond, run, st))
    print("%-6s %-8s %9.2f ft %9.1f%% %9.0f ft %9.0f m" %
          (storm, cond, st["maxhead"], st["pct"], st["ov_ft"], st["s_at_max"]))

if A.csv:
    with open(A.csv, "w") as fh:
        fh.write("storm,culvert8,run,max_head_ft,pct_crest_overtopped,crest_length_ft,overtopped_length_ft,station_of_max_m,east,north\n")
        for storm, cond, run, st in rows:
            fh.write("%s,%s,%s,%.2f,%.1f,%.0f,%.0f,%.0f,%.1f,%.1f\n" %
                     (storm, cond, run, st["maxhead"], st["pct"], st["length_ft"],
                      st["ov_ft"], st["s_at_max"], st["x"], st["y"]))
    print("\nwrote " + A.csv)

# --------------------------------------------------------------- sensitivity --
if A.sensitivity:
    print("\nsensitivity of max head (ft) / %% crest overtopped to the crest derivation:")
    print("  %-26s %s" % ("crest definition", "  ".join("%-16s" % ("%s %s" % (s, c)) for s, c, _ in RUNS)))
    def line(lbl, s_, x_, y_, z_, v_=None):
        v_ = np.ones(len(s_), bool) if v_ is None else v_
        cells = []
        for storm, cond, run in RUNS:
            w, _ = read(os.path.join(RES, run, "wse_max.tif"))
            h = (samp(w, x_, y_) - z_) * FT
            h = np.where(np.isfinite(h), h, 0.0)[v_]
            cells.append("%6.2f / %5.1f%%" % (h.max(), 100 * (h > A.threshold).mean()))
        print("  %-26s %s" % (lbl, "  ".join("%-16s" % c for c in cells)))
    pts = load_line(A.roads, A.linearid)
    s0, x0, y0 = stations(pts, A.step)
    line("TIGER centreline as-is", s0, x0, y0, samp(zraw, x0, y0))
    for w_ in (10.0, 15.0, 20.0, 25.0, 30.0):
        s_, x_, y_, z_, _, e_, v_ = derive_crest(w_, A.smooth)
        line("snap +/-%2.0f m, med %.0f m" % (w_, A.smooth), s_, x_, y_, z_, v_)
    for sm in (0.0, 51.0, 201.0):
        s_, x_, y_, z_, _, e_, v_ = derive_crest(A.window, sm)
        line("snap +/-%2.0f m, med %3.0f m" % (A.window, sm), s_, x_, y_, z_, v_)

# ------------------------------------------------------------------ QC figure --
if A.qc_figure:
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource
    from matplotlib.lines import Line2D
    pts = load_line(A.roads, A.linearid)
    s0, x0, y0 = stations(pts, A.step)
    hs = LightSource(azdeg=315, altdeg=45).hillshade(np.nan_to_num(zraw, nan=np.nanmin(zraw)),
                                                     vert_exag=3, dx=1, dy=1)
    fig = plt.figure(figsize=(15, 9.5), dpi=140)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.25, 1], hspace=0.28)
    ax = fig.add_subplot(gs[0])
    pad = 90
    ax.imshow(hs, cmap="gray", extent=EXT, origin="upper", zorder=1)
    ax.plot(x0, y0, c="#1f77b4", lw=1.6, ls="--", zorder=4, label="TIGER centreline (as supplied)")
    ax.plot(X, Y, c="#e8112d", lw=2.0, zorder=5, label="Derived road crest (used)")
    h6 = head_ft(RUNS[0][2])
    ov = (h6 > A.threshold) & VALID
    ax.scatter(X[ov], Y[ov], s=9, c="#ffd400", ec="none", zorder=6,
               label="Overtopped, 6-hr culvert-8 open")
    jj, ii = np.where(NOTCH)
    ax.scatter(GT[0] + ii + .5, GT[3] - jj - .5, s=6, c="#00d0ff", zorder=7, label="Culvert-8 notch")
    ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(min(Y.min(), y0.min()) - pad, max(Y.max(), y0.max()) + pad)
    ax.set_aspect("equal"); ax.legend(loc="lower left", fontsize=8.5, framealpha=.93)
    ax.set_title("C10 — Navajo Route 9 crest derivation (hillshade: raw 1 m DEM)", fontsize=12)
    ax.set_xlabel("UTM Zone 12N Easting (m)"); ax.set_ylabel("Northing (m)")
    ax2 = fig.add_subplot(gs[1])
    sf = S * FT
    ax2.fill_between(sf, Z * FT, Z.min() * FT - 10, color="#c8b89a", zorder=2, label="Road crest (derived)")
    ax2.plot(sf, Z * FT, c="#6b5b3e", lw=1.5, zorder=3)
    ax2.plot(S * FT, samp(zraw, x0, y0) * FT, c="#1f77b4", lw=1.2, ls="--", zorder=3,
             label="Ground on the TIGER centreline")
    for (storm, cond, run), col in zip(RUNS, ("#1a4f95", "#7a1020", "#2e8b57", "#8a3800")):
        h = head_ft(run)
        ax2.plot(sf, (Z + h / FT) * FT, lw=1.5, c=col, zorder=4, label="Peak WSE — %s, %s" % (storm, cond))
    ax2.set_xlim(0, sf[-1]); ax2.set_ylim(Z.min() * FT - 10, (Z.max() + 2.5) * FT)
    ax2.set_xlabel("Distance along Navajo Route 9, west → east (ft)")
    ax2.set_ylabel("Elevation (ft NAVD88)")
    ax2.grid(alpha=.3, zorder=1); ax2.legend(fontsize=8, loc="upper left", ncol=2)
    ax2.set_title("Crest profile and peak water surfaces", fontsize=11)
    fig.savefig(A.qc_figure, bbox_inches="tight")
    print("\nwrote " + A.qc_figure)
