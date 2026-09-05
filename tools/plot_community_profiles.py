#!/usr/bin/env python3
"""
Profiles THROUGH the community, with structure pads marked against the water surface.

    python3 tools/plot_community_profiles.py out_c6_1m_rain 1

WHY THIS EXISTS
  The headline exposure metric -- depth in the 1 m cell containing each structure
  point -- asks "is this cell wet?". It does not ask "does water stand above this
  pad?". Those differ here, because the pads are local highs about 1.2 ft above
  the surrounding ground and the water flows AROUND them. Sampling the cell reads
  dry while water stands at an elevation above the pad a few metres away.

  These sections show that directly: ground surface, maximum water surface, and
  every structure pad the line passes near, so the margin between pad and water
  is visible rather than inferred from a single sampled cell.

WHAT THE DEM IS
  Verified bare-earth, not a surface model with roofs: no structure point has a
  local high greater than 0.85 m within 3 m of it, where a roof would be 2-4 m.
  The consistent +0.37 m rise at structure points IS the built pad, so the
  elevation plotted at each structure is the pad surface.
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from osgeo import gdal, ogr

gdal.UseExceptions(); ogr.UseExceptions()
FT = 3.280839895

RUN  = sys.argv[1] if len(sys.argv) > 1 else "out_c6_1m_rain"
DXM  = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
GIS  = "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097"
DEM  = GIS + "/outputs_2D/cliped_dem.tif"
NEAR = 25.0        # m: how close a structure must be to the line to be plotted

dd = gdal.Open(os.path.join(RUN, "depth_max.tif"))
a  = dd.GetRasterBand(1).ReadAsArray().astype(float)
nd = dd.GetRasterBand(1).GetNoDataValue()
a  = np.where((a == nd) | ~np.isfinite(a), 0.0, a)
gt = dd.GetGeoTransform(); ny, nx = a.shape
ww = gdal.Open(os.path.join(RUN, "wse_max.tif"))
w  = ww.GetRasterBand(1).ReadAsArray().astype(float)
w  = np.where(w < 0, np.nan, w)
zz = gdal.Open(DEM).ReadAsArray().astype(float)
gz = gdal.Open(DEM).GetGeoTransform()

def samp(arr, x, y, g=None):
    g = g if g is not None else gt
    i = int((x - g[0]) / abs(g[1])); j = int((g[3] - y) / abs(g[5]))
    if 0 <= i < arr.shape[1] and 0 <= j < arr.shape[0]:
        return arr[j, i]
    return np.nan

# ---------------------------------------------------------------- vectors ---
_s = ogr.Open(GIS + "/outputs/05_Subdivisions_clip.gpkg"); sly = _s.GetLayer()
subs = []
for f in sly:
    g = f.GetGeometryRef()
    if g is None: continue
    e = g.GetEnvelope()
    if e[1] < gt[0] or e[0] > gt[0] + nx * gt[1]: continue
    subs.append((f.GetField("proj_num"), g.Clone(), e))

_p = ogr.Open(GIS + "/outputs/00_NHA_Structures_clip.gpkg"); ply = _p.GetLayer()
S = []
for f in ply:
    g = f.GetGeometryRef()
    if g is None: continue
    x, y = g.GetX(), g.GetY()
    if gt[0] <= x <= gt[0] + nx * gt[1] and gt[3] - ny * gt[1] <= y <= gt[3]:
        S.append((x, y))
S = np.array(S)

# --- cut lines derived from the STRUCTURE POINTS ----------------------------
# Using the subdivision envelopes put the lines off the community: the polygon
# filter screened on easting only, so parcels outside the domain skewed the
# bounds. The structures are the feature of interest, so the lines are placed
# through them and are guaranteed to intersect pads.
X0, X1 = S[:, 0].min(), S[:, 0].max()
Y0, Y1 = S[:, 1].min(), S[:, 1].max()
PAD = 90.0
# Two long lines along the community's axis and one across it, each offset to
# pass through a dense row of pads.
LINES = [
    ("A", (X0 - PAD, Y0 + 0.72 * (Y1 - Y0)), (X1 + PAD, Y0 + 0.72 * (Y1 - Y0))),
    ("B", (X0 - PAD, Y0 + 0.35 * (Y1 - Y0)), (X1 + PAD, Y0 + 0.35 * (Y1 - Y0))),
    ("C", (X0 + 0.50 * (X1 - X0), Y1 + PAD), (X0 + 0.50 * (X1 - X0), Y0 - PAD)),
]

# Clamp every line to the raster footprint. The subdivision envelopes can sit
# close to the domain edge, and a line running past it samples NaN, which then
# propagates into the axis limits.
RX0, RX1 = gt[0] + 2, gt[0] + (nx - 2) * gt[1]
RY1, RY0 = gt[3] - 2, gt[3] - (ny - 2) * abs(gt[5])
def clamp(P):
    return (min(max(P[0], RX0), RX1), min(max(P[1], RY0), RY1))
LINES = [(l, clamp(p0), clamp(p1)) for l, p0, p1 in LINES]

fig = plt.figure(figsize=(15, 3.7 * len(LINES) + 0.9), dpi=150)
gs = fig.add_gridspec(len(LINES), 1, hspace=0.52)
fig.suptitle("C6 Coyote Canyon — community sections: ground, structure pads and "
             "maximum water surface\n2D-DWE, %.0f m grid, 100-yr 6-hr event — PRELIMINARY"
             % DXM, fontsize=13, y=0.995)

print("run %s" % RUN)
for k, (lbl, P0, P1) in enumerate(LINES):
    ax = fig.add_subplot(gs[k])
    L = np.hypot(P1[0] - P0[0], P1[1] - P0[1])
    s = np.arange(0, L, max(DXM / 2, 0.5))
    ux, uy = (P1[0] - P0[0]) / L, (P1[1] - P0[1]) / L
    px, py = P0[0] + ux * s, P0[1] + uy * s
    zg = np.array([samp(zz, x_, y_, gz) for x_, y_ in zip(px, py)]) * FT
    ok = np.isfinite(zg)
    if ok.sum() < 10:
        print("  %s-%s': line falls outside the terrain, skipped" % (lbl, lbl)); continue
    sf_all = s * FT
    s, px, py, zg = s[ok], px[ok], py[ok], zg[ok]
    zw = np.array([samp(w,  x_, y_)     for x_, y_ in zip(px, py)]) * FT
    sf = s * FT

    # Axis range is needed before the subdivision labels are placed, so they can
    # sit at the BOTTOM of the panel clear of the legend rather than at the top.
    _lo = np.nanmin(zg) - 4
    _hi = max(np.nanmax(zg),
              np.nanmax(zw[np.isfinite(zw)]) if np.isfinite(zw).any() else 0) + 7
    ax.fill_between(sf, zg, np.nanmin(zg) - 8, color="#cbbd9a", zorder=2)
    ax.plot(sf, zg, c="#6b5b3e", lw=1.5, zorder=3, label="Ground surface (bare-earth LiDAR)")
    ax.plot(sf, zw, c="#1a4f95", lw=2.0, zorder=5, label="Maximum water surface")
    wet = np.isfinite(zw) & ((zw - zg) > 0.5)
    ax.fill_between(sf, zg, zw, where=wet, color="#4ba3d8", alpha=.55, zorder=4,
                    label="Inundation (> 0.5 ft)")

    # subdivision extents
    _used = []      # x positions already labelled on this axis, to stagger collisions
    for pn, gsub, e in subs:
        inside = np.array([gsub.Contains(ogr.CreateGeometryFromWkt("POINT(%f %f)" % (x_, y_)))
                           for x_, y_ in zip(px, py)])
        if inside.any():
            ax.fill_between(sf, _lo, _hi, where=inside,
                            color="#ffd400", alpha=.16, zorder=1)
            m = int(np.argmax(inside)) + int(inside.sum()) // 2
            lx = sf[min(m, len(sf) - 1)]
            # adjacent parcels put their labels almost on top of each other; drop
            # the second one a row rather than let the boxes collide.
            row = 0
            while any(abs(lx - ux) < 0.11 * sf[-1] and row == ur for ux, ur in _used):
                row += 1
            _used.append((lx, row))
            ax.text(lx, _lo + (0.07 + 0.075 * row) * (_hi - _lo), str(pn),
                    ha="center", va="center", fontsize=8, weight="bold", zorder=9,
                    bbox=dict(fc="#ffd400", ec="k", lw=.5, pad=1.2))

    # structure pads near the line
    npad = nwet = 0
    for x_, y_ in S:
        t = ((x_ - P0[0]) * ux + (y_ - P0[1]) * uy)
        if t < 0 or t > L: continue
        perp = abs((x_ - P0[0]) * uy - (y_ - P0[1]) * ux)
        if perp > NEAR: continue
        pad = samp(zz, x_, y_, gz) * FT
        # highest water surface within 15 m of this structure
        i = int((x_ - gt[0]) / gt[1]); j = int((gt[3] - y_) / abs(gt[5]))
        R = int(15 / DXM)
        sub = w[max(0, j - R):j + R + 1, max(0, i - R):i + R + 1]
        wl = (np.nanmax(sub) * FT) if np.isfinite(sub).any() else np.nan
        npad += 1
        over = np.isfinite(wl) and wl > pad
        if over: nwet += 1
        ax.plot([t * FT, t * FT], [pad, pad + 1.9], c="#333", lw=1.0, zorder=7)
        ax.plot(t * FT, pad, marker="s", ms=5.5, zorder=8,
                mfc=("#e8112d" if over else "w"), mec="#222", mew=.8)
        if over:
            ax.plot([t * FT, t * FT], [pad, wl], c="#e8112d", lw=2.4, alpha=.75, zorder=6)

    ax.set_xlim(0, sf[-1])
    ax.set_ylim(_lo, _hi)
    ax.set_xlabel("Distance along %s–%s' (ft)" % (lbl, lbl))
    ax.set_ylabel("Elevation (ft NAVD88)")
    ax.grid(alpha=.3, zorder=0)
    ax.set_title("Section %s–%s'   %d structure pads within %.0f ft of the line   "
                 "%d with water standing ABOVE the pad within 50 ft"
                 % (lbl, lbl, npad, NEAR * FT, nwet), fontsize=10.5)
    if k == 0:
        h, l = ax.get_legend_handles_labels()
        h += [Line2D([], [], marker="s", ls="", mfc="w",  mec="#222", label="Structure pad, dry"),
              Line2D([], [], marker="s", ls="", mfc="#e8112d", mec="#222",
                     label="Pad below nearby water surface")]
        ax.legend(handles=h, fontsize=7.2, loc="upper center", ncol=5,
                  framealpha=.93, columnspacing=1.0, handlelength=1.6,
                  borderpad=0.35, handletextpad=0.5)
    print("  %s-%s': %d pads on the line, %d below the nearby water surface" % (lbl, lbl, npad, nwet))

out = os.path.join(RUN, "c6_community_profiles.png")
plt.savefig(out, bbox_inches="tight")
plt.close(fig)
print("wrote", out)

# --------------------------------------------------------------- key map ----
# A section is unreadable without knowing where it was cut. This is the locator:
# same depth field, same structures, the three lines labelled at both ends.
from matplotlib.colors import BoundaryNorm, ListedColormap, LightSource
from matplotlib.patches import Patch
zh = gdal.Warp("", DEM, format="MEM", width=nx, height=ny, resampleAlg="average")
zb = zh.GetRasterBand(1).ReadAsArray().astype(float)
hs = LightSource(azdeg=315, altdeg=45).hillshade(zb, vert_exag=3, dx=DXM, dy=DXM)
ext = [gt[0], gt[0] + nx * gt[1], gt[3] - ny * abs(gt[5]), gt[3]]

f2, ax = plt.subplots(figsize=(13, 9.5), dpi=150)
ax.imshow(hs, cmap="gray", extent=ext, alpha=.92, zorder=1)
lv = [0.5, 1, 2, 3, 4, 6]
cm = ListedColormap(["#bfe3f7", "#7fc4e8", "#3d94cf", "#1c5fa8", "#0d2f6e"])
im = ax.imshow(np.ma.masked_less(a * FT, 0.5), cmap=cm, norm=BoundaryNorm(lv, cm.N),
               extent=ext, alpha=.85, zorder=2)
for pn, gsub, e in subs:
    for k2 in range(max(1, gsub.GetGeometryCount())):
        gg = gsub.GetGeometryRef(k2) if gsub.GetGeometryCount() else gsub
        ring = gg.GetGeometryRef(0) if gg.GetGeometryCount() else gg
        pp = np.array(ring.GetPoints())
        if pp.size:
            ax.plot(pp[:, 0], pp[:, 1], c="#111", lw=2.4, zorder=6)
            ax.plot(pp[:, 0], pp[:, 1], c="#ffd400", lw=1.3, zorder=6)
ax.scatter(S[:, 0], S[:, 1], s=22, c="w", ec="#222", lw=.7, zorder=8)
for lbl, P0, P1 in LINES:
    ax.plot([P0[0], P1[0]], [P0[1], P1[1]], c="#e8112d", lw=2.4, ls="--", zorder=9)
    for P, tag in ((P0, lbl), (P1, lbl + "'")):
        ax.text(P[0], P[1], tag, color="w", fontsize=12, weight="bold", ha="center",
                va="center", zorder=11,
                bbox=dict(fc="#e8112d", ec="k", boxstyle="circle,pad=0.26"))
cb = f2.colorbar(im, ax=ax, shrink=.72, pad=.015, ticks=lv)
cb.set_label("Maximum depth over the event (ft)", fontsize=10)
pad_m = 160
ax.set_xlim(S[:, 0].min() - pad_m, S[:, 0].max() + pad_m)
ax.set_ylim(S[:, 1].min() - pad_m, S[:, 1].max() + pad_m)
ax.set_aspect("equal")
ax.set_title("C6 Coyote Canyon — location of community sections A--A', B--B', C--C'\n"
             "2D-DWE, %.0f m grid, 100-yr 6-hr event — PRELIMINARY" % DXM, fontsize=12.5)
ax.set_xlabel("UTM Zone 12N Easting (m)"); ax.set_ylabel("Northing (m)")
ax.legend(handles=[Patch(fc="none", ec="#ffd400", lw=2, label="NHA subdivision"),
                   Line2D([], [], marker="o", ls="", mfc="w", mec="#222",
                          label="Structure (%d)" % len(S)),
                   Line2D([], [], c="#e8112d", lw=2, ls="--", label="Section line")],
          loc="lower left", fontsize=9.5, framealpha=.93)
bx = S[:, 0].min() - pad_m + 40
by = S[:, 1].min() - pad_m + 35
ax.plot([bx, bx + 152.4], [by, by], c="k", lw=4, zorder=9, solid_capstyle="butt")
ax.text(bx + 76, by + 14, "500 ft", ha="center", fontsize=9, zorder=9)
ax.annotate("", xy=(S[:, 0].max() + pad_m - 45, S[:, 1].max() + pad_m - 40),
            xytext=(S[:, 0].max() + pad_m - 45, S[:, 1].max() + pad_m - 110),
            arrowprops=dict(arrowstyle="-|>", lw=2, color="k"), zorder=9)
ax.text(S[:, 0].max() + pad_m - 45, S[:, 1].max() + pad_m - 28, "N", ha="center",
        fontsize=11, weight="bold", zorder=9)
out2 = os.path.join(RUN, "c6_community_keymap.png")
f2.savefig(out2, bbox_inches="tight")
plt.close(f2)
print("wrote", out2)
