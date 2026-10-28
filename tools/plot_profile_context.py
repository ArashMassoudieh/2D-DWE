#!/usr/bin/env python3
"""
Zoomed depth map around the three profile cut-lines.

    python3 tools/plot_profile_context.py

Replaces the 24-hour depth map, which repeated the 6-hour one almost exactly
(the two water surfaces differ by at most 0.4 ft on any section). What the
report lacked instead was a view at the scale the profiles are actually read
at: the domain map is 1,633 m tall and the cut-lines occupy a 550 m band of it,
so the reader cannot see where a section sits relative to the parcels, the
structures or the crossing. This frames exactly that band.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from osgeo import gdal, ogr, osr

gdal.UseExceptions(); ogr.UseExceptions()
FT = 3.280839895
NHA = "/home/arash/Dropbox/Chloeta/NHA"
S   = NHA + "/WS3_GIS/NM15-141Local"
C6  = NHA + "/WS3_GIS/NM15-097/outputs/"
RUN = "/home/arash/Projects/2D-DWE/out_nm141_1m_6hr_carved_2m_conv"
OUT = S + "/Report_for_C10/Figures/NM15-141_2d_profile_context.png"
PAD = 90.0

def to_utm(path):
    ds = ogr.Open(path); ly = ds.GetLayer()
    sr = ly.GetSpatialRef(); tg = osr.SpatialReference(); tg.ImportFromEPSG(26912)
    sr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tg.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = None if sr.IsSame(tg) else osr.CoordinateTransformation(sr, tg)
    out = []
    for f in ly:
        g = f.GetGeometryRef()
        if g is None: continue
        g = g.Clone()
        if tr: g.Transform(tr)
        out.append(g)
    return out

def rings(g):
    r = []
    if g.GetGeometryName() == "MULTIPOLYGON":
        for k in range(g.GetGeometryCount()): r += rings(g.GetGeometryRef(k))
    elif g.GetGeometryName() == "POLYGON":
        p = np.array(g.GetGeometryRef(0).GetPoints())
        if p.size: r.append(p[:, :2])
    return r

dd = gdal.Open(RUN + "/depth_max.tif"); bd = dd.GetRasterBand(1)
D = bd.ReadAsArray().astype(float); nd = bd.GetNoDataValue()
D = np.where((D == nd) | ~np.isfinite(D), 0.0, D) * FT
gt = dd.GetGeoTransform(); ny, nx = D.shape
EXT = [gt[0], gt[0] + nx * gt[1], gt[3] + ny * gt[5], gt[3]]

t = gdal.Warp("", S + "/output_2D/dem_2d_carved_2m.tif", format="MEM",
              width=nx, height=ny, resampleAlg="average")
Z = t.GetRasterBand(1).ReadAsArray().astype(float)
hs = LightSource(azdeg=315, altdeg=45).hillshade(Z, vert_exag=3, dx=1, dy=1)

cuts = [np.array(g.GetPoints())[:, :2] for g in to_utm(S + "/output_2D/Profile_cut_line.shp")]
cuts.sort(key=lambda p: (0 if abs(p[-1][0] - p[0][0]) > abs(p[-1][1] - p[0][1]) else 1, p[0][0]))
xs = np.concatenate([c[:, 0] for c in cuts]); ys = np.concatenate([c[:, 1] for c in cuts])
X0, X1 = xs.min() - PAD, xs.max() + PAD
Y0, Y1 = ys.min() - PAD, ys.max() + PAD

fig, ax = plt.subplots(figsize=(13.5, 13.5 * (Y1 - Y0) / (X1 - X0) + 1.0), dpi=160)
ax.imshow(hs, cmap="gray", extent=EXT, origin="upper", alpha=.92, zorder=1)
lv = [0.5, 1, 2, 3, 5, 8, 12, 20]
cm = ListedColormap(["#c6e5f5", "#8ec9ec", "#4ba3d8", "#2272b5", "#1a4f95", "#0d2f6e", "#5a1a6e"])
im = ax.imshow(np.ma.masked_less(D, 0.5), cmap=cm, norm=BoundaryNorm(lv, cm.N),
               extent=EXT, origin="upper", alpha=.82, zorder=2)
for g in to_utm(S + "/roads/tl_2025_35031_roads/tl_2025_35031_roads.shp"):
    e = g.GetEnvelope()
    if e[1] < X0 or e[0] > X1 or e[3] < Y0 or e[2] > Y1: continue
    for k in range(max(1, g.GetGeometryCount())):
        gg = g.GetGeometryRef(k) if g.GetGeometryCount() else g
        q = np.array(gg.GetPoints())
        if q.size:
            ax.plot(q[:, 0], q[:, 1], c="#5a3010", lw=3.0, zorder=5, solid_capstyle="round")
            ax.plot(q[:, 0], q[:, 1], c="#d9a441", lw=1.4, zorder=5)
for g in to_utm(C6 + "05_Subdivisions_clip.gpkg"):
    for q in rings(g):
        ax.plot(q[:, 0], q[:, 1], c="#111", lw=2.4, zorder=6)
        ax.plot(q[:, 0], q[:, 1], c="#ffd400", lw=1.3, zorder=6)
sx = [g.GetX() for g in to_utm(C6 + "00_NHA_Structures_clip.gpkg") if X0 <= g.GetX() <= X1 and Y0 <= g.GetY() <= Y1]
sy = [g.GetY() for g in to_utm(C6 + "00_NHA_Structures_clip.gpkg") if X0 <= g.GetX() <= X1 and Y0 <= g.GetY() <= Y1]
ax.scatter(sx, sy, s=26, c="w", ec="#222", lw=.8, zorder=8)
ax.contour(D, levels=[1.0], colors=["#e8112d"], linewidths=0.7, alpha=.75,
           extent=EXT, origin="upper", zorder=7)
for i, p in enumerate(cuts):
    L = "ABC"[i]
    ax.plot(p[:, 0], p[:, 1], c="#e8112d", lw=2.4, ls="--", zorder=9)
    for q, lab in ((p[0], L), (p[-1], L + "$'$")):
        ax.text(q[0], q[1], lab, color="w", fontsize=11, weight="bold", ha="center", va="center",
                zorder=11, bbox=dict(fc="#e8112d", ec="k", boxstyle="circle,pad=0.26"))
ax.plot([711237], [3960470], marker="*", ms=19, mfc="#e8112d", mec="k", mew=1.1, zorder=10)
ax.annotate("Route 9 crossing", xy=(711237, 3960470), xytext=(711300, 3960530),
            fontsize=9.5, weight="bold", arrowprops=dict(arrowstyle="->", lw=1.5, color="k"), zorder=12)
cb = fig.colorbar(im, ax=ax, shrink=.70, pad=.012, aspect=28, ticks=lv)
cb.set_label("Maximum depth over the 6-hour event (ft)", fontsize=10)
ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_aspect("equal")
ax.set_xlabel("UTM Zone 12N Easting (m)"); ax.set_ylabel("Northing (m)")
ax.set_title("C10 — the profile cut-lines in context\n"
             "maximum depth, 100-yr 6-hour event, 1 m grid, 2 m culvert notch", fontsize=12.5)
ax.legend(handles=[Line2D([], [], c="#e8112d", lw=2.2, ls="--", label="Profile cut-line"),
                   Line2D([], [], c="#e8112d", lw=1.6, label="Depth $>$ 1 ft"),
                   Patch(fc="none", ec="#ffd400", lw=2, label="NHA subdivision"),
                   Line2D([], [], marker="o", ls="", mfc="w", mec="#222", label="Structure (%d)" % len(sx)),
                   Line2D([], [], c="#d9a441", lw=2.2, label="Road"),
                   Line2D([], [], marker="*", ls="", mfc="#e8112d", mec="k", ms=12, label="Route 9 crossing")],
          loc="lower left", fontsize=9, framealpha=.93)
bl = X1 - 152.4 - 0.05 * (X1 - X0)   # lower right; the legend holds the lower left
ax.plot([bl, bl + 152.4], [Y0 + 0.035 * (Y1 - Y0)] * 2, c="k", lw=4, zorder=12, solid_capstyle="butt")
ax.text(bl + 76, Y0 + 0.048 * (Y1 - Y0), "500 ft", ha="center", fontsize=9, zorder=12)
nx_ = X1 - 0.06 * (X1 - X0)
ax.annotate("", xy=(nx_, Y1 - 0.04 * (Y1 - Y0)), xytext=(nx_, Y1 - 0.11 * (Y1 - Y0)),
            arrowprops=dict(arrowstyle="-|>", lw=2, color="k"), zorder=12)
ax.text(nx_, Y1 - 0.025 * (Y1 - Y0), "N", ha="center", fontsize=11, weight="bold", zorder=12)
fig.tight_layout(); fig.savefig(OUT, bbox_inches="tight")
print("wrote", OUT)
print("extent x %.0f..%.0f  y %.0f..%.0f  (%.0f x %.0f m)" % (X0, X1, Y0, Y1, X1 - X0, Y1 - Y0))
