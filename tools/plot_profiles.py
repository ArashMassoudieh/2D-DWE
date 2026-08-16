#!/usr/bin/env python3
"""
Plan view plus ground-surface / max-WSEL profiles along a cut-line shapefile.

    python3 tools/plot_profiles.py out_nm141_1m_v3_carved_2m 1

The map carries subdivision boundaries, roads and structures; each profile
carries the same subdivisions and road crossings, so a feature can be followed
between the two views.

WHY THE PROFILE GROUND IS THE MODEL TERRAIN, NOT THE RAW DEM
  The water surface comes from a run on the carved terrain, so plotting it over
  the un-carved ground would show water standing below the bed at the culvert.
  The model terrain is used, and the notch is annotated where the cut line
  crosses it.
"""
import sys, os, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, LightSource
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from osgeo import gdal, ogr, osr

gdal.UseExceptions(); ogr.UseExceptions()
FT = 3.280839895

RUN = sys.argv[1] if len(sys.argv) > 1 else "out_nm141_1m_v3_carved_2m"
DXM = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
SITE = "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-141Local"
C6   = "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097/outputs/"
# terrain must match the run: plotting a carved profile under a blocked-run
# water surface would show water standing in a notch the run never had.
TERR = sys.argv[3] if len(sys.argv) > 3 else (
    SITE + "/output_2D/dem_2d_raw.tif" if RUN.endswith("raw")
    else SITE + "/output_2D/dem_2d_carved_2m.tif")

# ---------------------------------------------------------------- rasters ---
def band(path):
    ds = gdal.Open(path)
    a = ds.GetRasterBand(1).ReadAsArray().astype(float)
    nd = ds.GetRasterBand(1).GetNoDataValue()
    gt = ds.GetGeoTransform()
    return a, nd, gt

dep, ndd, gt = band(os.path.join(RUN, "depth_max.tif"))
dep = np.where((dep == ndd) | ~np.isfinite(dep), 0.0, dep)
wse, ndw, _ = band(os.path.join(RUN, "wse_max.tif"))
wse = np.where((wse == ndw) | ~np.isfinite(wse), np.nan, wse)
ny, nx = dep.shape
ext = [gt[0], gt[0] + nx * gt[1], gt[3] - ny * abs(gt[5]), gt[3]]

t = gdal.Warp("", TERR, format="MEM", width=nx, height=ny, resampleAlg="average")
z = t.GetRasterBand(1).ReadAsArray().astype(float)
hs = LightSource(azdeg=315, altdeg=45).hillshade(z, vert_exag=3, dx=DXM, dy=DXM)

def samp(arr, x, y):
    i = int((x - gt[0]) / gt[1]); j = int((gt[3] - y) / abs(gt[5]))
    if 0 <= i < nx and 0 <= j < ny:
        return arr[j, i]
    return np.nan

# ---------------------------------------------------------------- vectors ---
def load(path, want_epsg=26912):
    ds = ogr.Open(path); ly = ds.GetLayer()
    sr = ly.GetSpatialRef()
    tr = None
    if sr:
        tgt = osr.SpatialReference(); tgt.ImportFromEPSG(want_epsg)
        sr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        tgt.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        if not sr.IsSame(tgt):
            tr = osr.CoordinateTransformation(sr, tgt)
    out = []
    d = ly.GetLayerDefn()
    flds = [d.GetFieldDefn(i).GetName() for i in range(d.GetFieldCount())]
    for f in ly:
        g = f.GetGeometryRef()
        if g is None:
            continue
        g = g.Clone()
        if tr:
            g.Transform(tr)
        out.append((g, {k: f.GetField(k) for k in flds}))
    return out

subs   = load(C6 + "05_Subdivisions_clip.gpkg")
strucs = load(C6 + "00_NHA_Structures_clip.gpkg")
roads  = load(SITE + "/roads/tl_2025_35031_roads/tl_2025_35031_roads.shp")
cuts   = load(SITE + "/output_2D/Profile_cut_line.shp")

def inview(g):
    e = g.GetEnvelope()
    return not (e[1] < ext[0] or e[0] > ext[1] or e[3] < ext[2] or e[2] > ext[3])

subs   = [(g, a) for g, a in subs   if inview(g)]
roads  = [(g, a) for g, a in roads  if inview(g)]
sx, sy = [], []
for g, a in strucs:
    if ext[0] <= g.GetX() <= ext[1] and ext[2] <= g.GetY() <= ext[3]:
        sx.append(g.GetX()); sy.append(g.GetY())
sx, sy = np.array(sx), np.array(sy)

# order the cut lines: E-W first, then the two N-S by easting
def key(item):
    p = item[0].GetPoints()
    dx = abs(p[-1][0] - p[0][0]); dy = abs(p[-1][1] - p[0][1])
    return (0 if dx > dy else 1, p[0][0])
cuts.sort(key=key)
LBL = ["A", "B", "C", "D", "E"]

# ------------------------------------------------------------------ figure --
nprof = len(cuts)
# The domain is PORTRAIT (1,098 x 1,633 m). Giving the map a full-width cell with
# aspect="equal" letterboxes it into a narrow strip; it gets its own column and a
# summary panel fills the space beside it.
fig = plt.figure(figsize=(13.5, 11.5), dpi=150)
gs = fig.add_gridspec(1, 2, width_ratios=[1.12, 1], wspace=0.05)

ax = fig.add_subplot(gs[0, 0])
ax.imshow(hs, cmap="gray", extent=ext, alpha=.92, zorder=1)
lv = [0.5, 1, 2, 3, 5, 8, 12, 20]
cm = ListedColormap(["#c6e5f5", "#8ec9ec", "#4ba3d8", "#2272b5", "#1a4f95", "#0d2f6e", "#5a1a6e"])
im = ax.imshow(np.ma.masked_less(dep * FT, 0.5), cmap=cm, norm=BoundaryNorm(lv, cm.N),
               extent=ext, alpha=.85, zorder=2)

for g, a in roads:
    for k in range(max(1, g.GetGeometryCount())):
        gg = g.GetGeometryRef(k) if g.GetGeometryCount() else g
        p = np.array(gg.GetPoints())
        if p.size:
            ax.plot(p[:, 0], p[:, 1], c="#5a3010", lw=3.0, zorder=5, solid_capstyle="round")
            ax.plot(p[:, 0], p[:, 1], c="#d9a441", lw=1.4, zorder=5)

_SLOT = []
for g, a in subs:
    for k in range(max(1, g.GetGeometryCount())):
        gg = g.GetGeometryRef(k) if g.GetGeometryCount() else g
        ring = gg.GetGeometryRef(0) if gg.GetGeometryCount() else gg
        p = np.array(ring.GetPoints())
        if p.size:
            ax.plot(p[:, 0], p[:, 1], c="#111", lw=2.6, zorder=6)
            ax.plot(p[:, 0], p[:, 1], c="#ffd400", lw=1.4, zorder=6)
    e = g.GetEnvelope()
    _SLOT.append(str(a.get("proj_num", "")))
    ax.annotate(str(a.get("proj_num", "")),
                xy=((e[0] + e[1]) / 2, (e[2] + e[3]) / 2),
                xytext=(ext[1] - 0.34 * (ext[1] - ext[0]),
                        ext[2] + 0.42 * (ext[3] - ext[2]) - 95 * (len(_SLOT) - 1)),
                fontsize=8.5, weight="bold", zorder=12, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", lw=.8, color="k"),
                bbox=dict(fc="#ffd400", ec="k", lw=.6, pad=1.6))

# The 1 ft threshold traced explicitly. The graduated fill shows how deep, but the
# 1 ft line is what the hazard tables are keyed to, so it is drawn as a hard edge
# rather than left to be inferred from a colour break.
ax.contour(dep * FT, levels=[1.0], colors=["#e8112d"], linewidths=1.5,
           extent=ext, origin="upper", zorder=7)

ax.scatter(sx, sy, s=26, c="w", ec="#222", lw=.8, zorder=8)

for i, (g, a) in enumerate(cuts):
    p = np.array(g.GetPoints())
    ax.plot(p[:, 0], p[:, 1], c="#e8112d", lw=2.2, ls="--", zorder=9)
    ax.text(p[0, 0], p[0, 1], LBL[i], color="w", fontsize=11, weight="bold", ha="center",
            va="center", zorder=11, bbox=dict(fc="#e8112d", ec="k", boxstyle="circle,pad=0.25"))
    ax.text(p[-1, 0], p[-1, 1], LBL[i] + "'", color="w", fontsize=11, weight="bold",
            ha="center", va="center", zorder=11,
            bbox=dict(fc="#e8112d", ec="k", boxstyle="circle,pad=0.25"))

cb = fig.colorbar(im, ax=ax, shrink=.62, pad=.02, aspect=28, ticks=lv)
cb.set_label("Maximum depth over the event (ft)", fontsize=10)
cb.ax.tick_params(labelsize=9)
ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect("equal")
fig.suptitle("C10 -- maximum simulated depth, 100-yr 6-hr event\n"
             "2D-DWE, %.0f m grid, all three inflows + local runoff -- PRELIMINARY" % DXM,
             fontsize=13.5, y=0.965)
ax.set_xlabel("UTM Zone 12N Easting (m)"); ax.set_ylabel("Northing (m)")
ax.legend(handles=[Patch(fc="none", ec="#ffd400", lw=2, label="NHA subdivision"),
                   Line2D([], [], c="#d9a441", lw=2.5, label="Road"),
                   Line2D([], [], marker="o", ls="", mfc="w", mec="#222", label="Structure (%d)" % len(sx)),
                   Line2D([], [], c="#e8112d", lw=1.6,
                          label="Depth > 1 ft (%.1f ac)" % ((dep * FT > 1).sum() * DXM * DXM / 4046.8564224)),
                   Line2D([], [], c="#e8112d", lw=2, ls="--", label="Profile line")],
          loc="upper left", fontsize=8.5, framealpha=.93)
bl = ext[0] + 420
ax.plot([bl, bl + 304.8], [ext[2] + 55, ext[2] + 55], c="k", lw=4, zorder=9, solid_capstyle="butt")
ax.text(bl + 152, ext[2] + 78, "1,000 ft", ha="center", fontsize=9, zorder=9)
ax.annotate("", xy=(ext[1] - 70, ext[3] - 60), xytext=(ext[1] - 70, ext[3] - 200),
            arrowprops=dict(arrowstyle="-|>", lw=2, color="k"), zorder=9)
ax.text(ext[1] - 70, ext[3] - 40, "N", ha="center", fontsize=11, weight="bold", zorder=9)

# -------------------------------------------------------------- summary panel
axs = fig.add_subplot(gs[0, 1]); axs.axis("off")
A_ac = DXM * DXM / 4046.8564224
lines = [("C10 -- 2D results", 13, "bold"), ("", 6, "normal")]
for thr in (1.0, 3.0, 6.0):
    lines.append(("Area with max depth > %.0f ft:   %.1f acres"
                  % (thr, (dep * FT > thr).sum() * A_ac), 11, "normal"))
lines += [("", 6, "normal"),
          ("Peak depth in domain:   %.1f ft" % (dep.max() * FT), 11, "normal"),
          ("Structures in view:   %d" % len(sx), 11, "normal"),
          ("", 6, "normal"),
          ("Inflows (HMS, 6-hr 100-yr)", 11.5, "bold"),
          ("  Reach_147  south    3,870 cfs", 10, "normal"),
          ("  Reach_122  west     3,374 cfs", 10, "normal"),
          ("  Reach_99   east       424 cfs", 10, "normal"),
          ("  local runoff (rain)   43.9 ac-ft", 10, "normal"),
          ("  total volume         915.9 ac-ft", 10, "normal"),
          ("  HMS Outlet           915.0 ac-ft", 10, "normal"),
          ("", 6, "normal"),
          ("Culvert represented as a 2 m carved notch", 9.5, "normal"),
          ("(~2.6 m barrel). Road overtops regardless.", 9.5, "normal")]
yy = 0.93
for txt, fs, w in lines:
    if txt:
        axs.text(0.0, yy, txt, transform=axs.transAxes, fontsize=fs, weight=w, va="top")
    yy -= (fs + 9) / 430.0

# ------------------------------------------------------------------ save map -
out_map = os.path.join(RUN, "nm141_map.png")
plt.savefig(out_map, bbox_inches="tight")
plt.close(fig)
print("wrote", out_map)

# ------------------------------------------------------------------ profiles -
figp = plt.figure(figsize=(15, 3.5 * nprof + 0.8), dpi=150)
gsp = figp.add_gridspec(nprof, 1, hspace=0.46)
figp.suptitle("C10 -- ground surface and maximum water-surface elevation\n"
              "2D-DWE, %.0f m grid, 100-yr 6-hr event, all three inflows + local runoff"
              " -- PRELIMINARY" % DXM, fontsize=13, y=0.995)

for i, (g, a) in enumerate(cuts):
    axp = figp.add_subplot(gsp[i])
    p = np.array(g.GetPoints())[:, :2]
    seg = np.r_[0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))]
    s = np.arange(0, seg[-1], max(DXM / 2, 0.5))
    px = np.interp(s, seg, p[:, 0]); py = np.interp(s, seg, p[:, 1])
    zg = np.array([samp(z, a_, b_) for a_, b_ in zip(px, py)]) * FT
    zw = np.array([samp(wse, a_, b_) for a_, b_ in zip(px, py)]) * FT
    sf = s * FT

    axp.fill_between(sf, zg, np.nanmin(zg) - 12, color="#cbbd9a", zorder=2)
    axp.plot(sf, zg, c="#6b5b3e", lw=1.6, zorder=3, label="Ground surface (model terrain)")
    axp.plot(sf, zw, c="#1a4f95", lw=2.0, zorder=5, label="Maximum WSEL")
    # WET_FT screens out the millimetre film that spreads over almost everything;
    # without it "inundated over 82 % of its length" counts hundredths of a foot.
    WET_FT = 0.5
    wetmask = np.isfinite(zw) & ((zw - zg) > WET_FT)
    axp.fill_between(sf, zg, zw, where=wetmask,
                     color="#4ba3d8", alpha=.55, zorder=4,
                     label="Inundation (> %.1f ft)" % WET_FT)

    # subdivision extents along the line
    line = ogr.Geometry(ogr.wkbLineString)
    for xx, yy in zip(px, py):
        line.AddPoint(float(xx), float(yy))
    for gsub, asub in subs:
        if not gsub.Intersects(line):
            continue
        inside = np.array([gsub.Contains(ogr.CreateGeometryFromWkt("POINT(%f %f)" % (xx, yy)))
                           for xx, yy in zip(px, py)])
        if inside.any():
            axp.fill_between(sf, np.nanmin(zg) - 12, np.nanmax(zg) + 8, where=inside,
                             color="#ffd400", alpha=.20, zorder=1)
            k = int(np.argmax(inside)) + int(inside.sum()) // 2
            axp.text(sf[min(k, len(sf) - 1)], np.nanmax(zg) + 4, str(asub.get("proj_num", "")),
                     ha="center", fontsize=8, weight="bold", zorder=9,
                     bbox=dict(fc="#ffd400", ec="k", lw=.5, pad=1.2))

    # road crossings
    for grd, _ in roads:
        it = grd.Intersection(line)
        if it is None or it.IsEmpty():
            continue
        pts = []
        if it.GetGeometryType() in (ogr.wkbPoint, ogr.wkbPoint25D):
            pts = [(it.GetX(), it.GetY())]
        else:
            for kk in range(it.GetGeometryCount()):
                sub = it.GetGeometryRef(kk)
                pts.append((sub.GetX(), sub.GetY()))
        for xx, yy in pts:
            d0 = float(s[np.argmin((px - xx) ** 2 + (py - yy) ** 2)]) * FT
            axp.axvline(d0, color="#d9a441", lw=2.4, zorder=6)
            axp.text(d0, np.nanmax(zg) + 4, "road", rotation=90, fontsize=7.5,
                     ha="right", va="top", color="#5a3010", zorder=9)

    axp.set_xlim(0, sf[-1])
    axp.set_ylim(np.nanmin(zg) - 6, max(np.nanmax(zg), np.nanmax(zw[np.isfinite(zw)]) if np.isfinite(zw).any() else 0) + 10)
    axp.set_xlabel("Distance along %s–%s' (ft)" % (LBL[i], LBL[i]))
    axp.set_ylabel("Elevation (ft NAVD88)")
    axp.grid(alpha=.3, zorder=0)
    wet = wetmask
    axp.set_title("Profile %s-%s'   length %s ft   inundated over %.0f%% of its length"
                  "   max depth %.1f ft"
                  % (LBL[i], LBL[i], format(int(sf[-1]), ","), 100 * wet.mean(),
                     np.nanmax((zw - zg)[wet]) if wet.any() else 0.0),
                  fontsize=10.5)
    if i == 0:
        axp.legend(fontsize=8.5, loc="upper right", ncol=2, framealpha=.92)

out = os.path.join(RUN, "nm141_profiles.png")
figp.savefig(out, bbox_inches="tight")
plt.close(figp)
print("wrote", out)
for i, (g, a) in enumerate(cuts):
    p = np.array(g.GetPoints())[:, :2]
    L = float(np.hypot(*(p[-1] - p[0])))
    print("  %s-%s'  %.0f m (%.0f ft)  (%.0f,%.0f) -> (%.0f,%.0f)"
          % (LBL[i], LBL[i], L, L * FT, p[0][0], p[0][1], p[-1][0], p[-1][1]))
