#!/usr/bin/env python3
"""
The two berm figures for Section 4.3, built together so their section labels
always agree.

    python3 tools/berm_figures.py

  NM15-141_berm_location_map.png   berm alignment, cross-section locations,
                                   community, and the governing flood depth
  NM15-141_berm_cross_sections.png the twelve perpendicular sections

Sections are picked once and shared, so A1..A6 / B1..B6 on the map are the same
cuts plotted in the section figure. Getting that wrong is the usual failure of
a location map drawn separately from the sections it indexes.
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
FIG = S + "/Report_for_C10/Figures/"
RUN = "/home/arash/Projects/2D-DWE/out_nm141_1m_6hr_raw_conv"   # governing case
TOE_FT, NSEC = 39.0, 6

dem = gdal.Open(S + "/demhr/cliped_sd.tif"); G = dem.GetGeoTransform()
Z = dem.GetRasterBand(1).ReadAsArray().astype(float); NY, NX = Z.shape
EXT = [G[0], G[0] + NX * G[1], G[3] + NY * G[5], G[3]]
def sdem(x, y):
    i = np.clip(((np.asarray(x) - G[0]) / G[1]).astype(int), 0, NX - 1)
    j = np.clip(((G[3] - np.asarray(y)) / abs(G[5])).astype(int), 0, NY - 1)
    return Z[j, i]

dd = gdal.Open(RUN + "/depth_max.tif"); bd = dd.GetRasterBand(1)
D = bd.ReadAsArray().astype(float); nd = bd.GetNoDataValue()
D = np.where((D == nd) | ~np.isfinite(D), 0.0, D) * FT
gt = dd.GetGeoTransform(); dny, dnx = D.shape
DEXT = [gt[0], gt[0] + dnx * gt[1], gt[3] + dny * gt[5], gt[3]]

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

berm = [np.array(g.GetPoints())[:, :2] for g in to_utm(S + "/outputs/existing_berm.shp")]
NAMES = ["Line 1 — west boundary", "Line 2 — north boundary"]
COLS  = ["#1a4f95", "#7a1020"]; PRE = ["A", "B"]

# ---- pick the sections ONCE, shared by both figures ------------------------
SEC = []
for k, p in enumerate(berm):
    seg = np.r_[0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))]
    s = np.arange(0, seg[-1], 1.0)
    px = np.interp(s, seg, p[:, 0]); py = np.interp(s, seg, p[:, 1])
    dx = np.gradient(px); dy = np.gradient(py)
    L = np.hypot(dx, dy); ux, uy = -dy / L, dx / L
    for n, i in enumerate(np.linspace(8, s[-1] - 8, NSEC).astype(int)):
        SEC.append(dict(line=k, lbl="%s%d" % (PRE[k], n + 1), s=s[i],
                        x=px[i], y=py[i], ux=ux[i], uy=uy[i]))
    berm[k] = (p, s, px, py, ux, uy)

# =========================================================== 1. LOCATION MAP
fig, ax = plt.subplots(figsize=(12.5, 10.5), dpi=160)
hs = LightSource(azdeg=315, altdeg=45).hillshade(Z, vert_exag=5, dx=1, dy=1)
ax.imshow(hs, cmap="gray", extent=EXT, origin="upper", zorder=1)
lv = [0.5, 1, 2, 3, 5, 8, 12, 20]
cm = ListedColormap(["#c6e5f5", "#8ec9ec", "#4ba3d8", "#2272b5", "#1a4f95", "#0d2f6e", "#5a1a6e"])
im = ax.imshow(np.ma.masked_less(D, 0.5), cmap=cm, norm=BoundaryNorm(lv, cm.N),
               extent=DEXT, origin="upper", alpha=.60, zorder=2)
for g in to_utm(S + "/roads/tl_2025_35031_roads/tl_2025_35031_roads.shp"):
    e = g.GetEnvelope()
    if e[1] < EXT[0] or e[0] > EXT[1] or e[3] < EXT[2] or e[2] > EXT[3]: continue
    for k in range(max(1, g.GetGeometryCount())):
        gg = g.GetGeometryRef(k) if g.GetGeometryCount() else g
        q = np.array(gg.GetPoints())
        if q.size: ax.plot(q[:, 0], q[:, 1], c="#8a6a3a", lw=1.6, zorder=4)
for g in to_utm(C6 + "05_Subdivisions_clip.gpkg"):
    for q in rings(g):
        ax.plot(q[:, 0], q[:, 1], c="#111", lw=2.2, zorder=6)
        ax.plot(q[:, 0], q[:, 1], c="#ffd400", lw=1.2, zorder=6)
sx = [g.GetX() for g in to_utm(C6 + "00_NHA_Structures_clip.gpkg")
      if EXT[0] <= g.GetX() <= EXT[1] and EXT[2] <= g.GetY() <= EXT[3]]
sy = [g.GetY() for g in to_utm(C6 + "00_NHA_Structures_clip.gpkg")
      if EXT[0] <= g.GetX() <= EXT[1] and EXT[2] <= g.GetY() <= EXT[3]]
ax.scatter(sx, sy, s=20, c="w", ec="#222", lw=.7, zorder=8)
for k, (p, s, px, py, ux, uy) in enumerate(berm):
    ax.plot(px, py, c=COLS[k], lw=3.0, zorder=7)
T = TOE_FT / FT
for d in SEC:
    x0, y0 = d["x"] - T * d["ux"], d["y"] - T * d["uy"]
    x1, y1 = d["x"] + T * d["ux"], d["y"] + T * d["uy"]
    ax.plot([x0, x1], [y0, y1], c="#e8112d", lw=1.9, zorder=9)
    ax.text(x1, y1, d["lbl"], color="w", fontsize=8, weight="bold", ha="center", va="center",
            zorder=11, bbox=dict(fc="#e8112d", ec="k", lw=.5, boxstyle="circle,pad=0.20"))
ax.plot([711427], [3960327], marker="o", ms=16, mfc="none", mec="#00b050", mew=2.6, zorder=10)
ax.annotate("entrance road\nthrough Line 2", xy=(711427, 3960327), xytext=(711468, 3960243),
            fontsize=9, weight="bold", color="#00702f",
            arrowprops=dict(arrowstyle="->", lw=1.6, color="#00b050"), zorder=12)
cb = fig.colorbar(im, ax=ax, shrink=.62, pad=.012, aspect=26, ticks=lv)
cb.set_label("Maximum depth, 6-hr 100-yr, culvert blocked (ft)", fontsize=9)
ax.set_xlim(711175, 711665); ax.set_ylim(3960010, 3960400); ax.set_aspect("equal")
ax.set_xlabel("UTM Zone 12N Easting (m)"); ax.set_ylabel("Northing (m)")
ax.set_title("C10 — existing perimeter berm: alignment and cross-section locations", fontsize=12.5)
ax.legend(handles=[Line2D([], [], c=COLS[0], lw=3, label="Line 1 — west boundary (sections A1–A6)"),
                   Line2D([], [], c=COLS[1], lw=3, label="Line 2 — north boundary (sections B1–B6)"),
                   Line2D([], [], c="#e8112d", lw=2, label="Cross-section (78 ft wide)"),
                   Patch(fc="none", ec="#ffd400", lw=2, label="NHA subdivision"),
                   Line2D([], [], marker="o", ls="", mfc="w", mec="#222", label="Structure"),
                   Line2D([], [], c="#8a6a3a", lw=1.6, label="Road (TIGER 2025)")],
          loc="lower right", fontsize=8.5, framealpha=.94)
bl = 711195
ax.plot([bl, bl + 100], [3960030, 3960030], c="k", lw=4, zorder=12, solid_capstyle="butt")
ax.text(bl + 50, 3960038, "100 m", ha="center", fontsize=8.5, zorder=12)
ax.annotate("", xy=(711640, 3960385), xytext=(711640, 3960345),
            arrowprops=dict(arrowstyle="-|>", lw=2, color="k"), zorder=12)
ax.text(711640, 3960390, "N", ha="center", fontsize=11, weight="bold", zorder=12)
fig.tight_layout(); fig.savefig(FIG + "NM15-141_berm_location_map.png", bbox_inches="tight")
print("wrote " + FIG + "NM15-141_berm_location_map.png")

# ====================================================== 2. CROSS-SECTIONS
offs = np.arange(-T, T + 1e-9, 0.5)
fig, axes = plt.subplots(2, NSEC, figsize=(3.05 * NSEC, 7.6), dpi=150)
for d in SEC:
    r, c = d["line"], int(d["lbl"][1:]) - 1
    ax = axes[r, c]
    xs = d["x"] + offs * d["ux"]; ys = d["y"] + offs * d["uy"]
    z = sdem(xs, ys) * FT
    chord = z[0] + (z[-1] - z[0]) * (offs - offs[0]) / (offs[-1] - offs[0])
    resid = z - chord
    win = np.abs(offs) <= 16.0 / FT
    h = np.nanmax(resid[win]); ko = int(np.argmax(np.where(win, resid, -np.inf)))
    ax.fill_between(offs * FT, chord, z, where=z > chord, color="#c8b89a", zorder=2)
    ax.plot(offs * FT, z, c="#5a3010", lw=1.7, zorder=4)
    ax.plot(offs * FT, chord, c="#1a4f95", lw=1.3, ls="--", zorder=3)
    ax.axvline(0, c="#999", lw=.8, ls=":", zorder=1)
    ax.annotate("", xy=(offs[ko] * FT, z[ko]), xytext=(offs[ko] * FT, chord[ko]),
                arrowprops=dict(arrowstyle="<->", lw=1.4, color="#e8112d"), zorder=6)
    ax.text(offs[ko] * FT + 2, (z[ko] + chord[ko]) / 2, "%.2f ft" % h, fontsize=8.5,
            color="#8a1015", weight="bold", va="center", zorder=7)
    ax.set_title("%s   s = %d m" % (d["lbl"], d["s"]), fontsize=9.5, weight="bold",
                 color=COLS[r])
    ax.set_xlabel("offset (ft)", fontsize=8); ax.tick_params(labelsize=7); ax.grid(alpha=.25)
    if c == 0: ax.set_ylabel("elevation (ft NAVD88)", fontsize=8.5)
axes[0, 0].plot([], [], c="#5a3010", lw=1.7, label="ground (1 m LiDAR)")
axes[0, 0].plot([], [], c="#1a4f95", lw=1.3, ls="--", label="toe-to-toe chord (the profile datum)")
axes[0, 0].legend(fontsize=7.5, loc="lower left")
fig.suptitle("C10 — existing berm, perpendicular cross-sections.  Top row A1–A6 (west boundary), "
             "bottom row B1–B6 (north boundary).\nHeight is measured above the toe-to-toe chord, "
             "which removes the hillside the berm sits on.", fontsize=11.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(FIG + "NM15-141_berm_cross_sections.png", bbox_inches="tight")
print("wrote " + FIG + "NM15-141_berm_cross_sections.png")
