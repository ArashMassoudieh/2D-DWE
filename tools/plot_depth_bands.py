#!/usr/bin/env python3
"""
Three hazard-band maps (> 1 ft, > 3 ft, > 6 ft) over the run's depth envelope,
with NHA subdivision boundaries and structures.

    python3 tools/plot_depth_bands.py out_c6_3m 3

WHY THESE THRESHOLDS
  1 ft  -- nuisance/access. Standing water on roads and yards; a passenger
           vehicle is already at risk of stalling and loss of control.
  3 ft  -- structural/life-safety. Above typical finished-floor height, and
           past the depth at which an adult loses footing in moving water.
  6 ft  -- destructive. Full submergence of a single storey's lower level.

WHY THE MAPS SHOW GROUND DEPTH, NOT BUILDING DEPTH
  Structures are not represented in the DEM and no finished-floor elevations
  are available. A structure plotted inside a band has that depth of water on
  the GROUND at its location; it does not follow that the building floods. This
  is stated on every figure so the distinction cannot be lost downstream.
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
from osgeo import gdal, ogr

gdal.UseExceptions(); ogr.UseExceptions()

RUN = sys.argv[1] if len(sys.argv) > 1 else "out_c6_3m"
DXM = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
# The subdivision and structure layers are NHA-wide clips held in the C6 folder;
# they cover both sites. Only the TERRAIN is site-specific, so it is an argument.
GIS = "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097"
DEM = sys.argv[3] if len(sys.argv) > 3 else GIS + "/outputs_2D/cliped_dem.tif"
FT  = 3.280839895
BANDS = [(1.0, "#2f8fd0", "Depth > 1 ft — nuisance and access"),
         (3.0, "#1c4f9c", "Depth > 3 ft — structural / life-safety"),
         (6.0, "#3b1055", "Depth > 6 ft — destructive")]

# ---------------------------------------------------------------- rasters ---
r = gdal.Open(os.path.join(RUN, "depth_max.tif"))
b = r.GetRasterBand(1)
a = b.ReadAsArray().astype(float)
nd = b.GetNoDataValue()
a = np.where((a == nd) | ~np.isfinite(a), 0.0, a) * FT      # -> feet
gt = r.GetGeoTransform(); ny, nx = a.shape
ext = [gt[0], gt[0] + nx * gt[1], gt[3] - ny * gt[1], gt[3]]
acre = DXM * DXM / 4046.8564224

t = gdal.Warp("", DEM, format="MEM", width=nx, height=ny, resampleAlg="average")
zb = t.GetRasterBand(1).ReadAsArray().astype(float)
hs = LightSource(azdeg=315, altdeg=45).hillshade(zb, vert_exag=3, dx=DXM, dy=DXM)

# ---------------------------------------------------------------- vectors ---
_sd = ogr.Open(GIS + "/outputs/05_Subdivisions_clip.gpkg"); sly = _sd.GetLayer()
subs = []
for f in sly:
    g = f.GetGeometryRef()
    if g is None:
        continue
    e = g.GetEnvelope()
    if e[1] < ext[0] or e[0] > ext[1] or e[3] < ext[2] or e[2] > ext[3]:
        continue
    subs.append((f.GetField("proj_num"), g.Clone()))

_st = ogr.Open(GIS + "/outputs/00_NHA_Structures_clip.gpkg"); ply = _st.GetLayer()
sx, sy = [], []
for f in ply:
    g = f.GetGeometryRef()
    if g is None:
        continue
    x, y = g.GetX(), g.GetY()
    i = int((x - ext[0]) / DXM); j = int((ext[3] - y) / DXM)
    if 0 <= i < nx and 0 <= j < ny:
        sx.append(x); sy.append(y)
sx, sy = np.array(sx), np.array(sy)
sd_ft = np.array([a[int((ext[3] - y) / DXM), int((x - ext[0]) / DXM)]
                  for x, y in zip(sx, sy)])


def rasterize(geom):
    mem = ogr.GetDriverByName("Memory").CreateDataSource("m")
    l2 = mem.CreateLayer("p", sly.GetSpatialRef(), ogr.wkbPolygon)
    ft_ = ogr.Feature(l2.GetLayerDefn()); ft_.SetGeometry(geom); l2.CreateFeature(ft_)
    msk = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Byte)
    msk.SetGeoTransform(gt); msk.SetProjection(r.GetProjection())
    gdal.RasterizeLayer(msk, [1], l2, burn_values=[1])
    return msk.GetRasterBand(1).ReadAsArray().astype(bool)


submask = [(pn, rasterize(g), g) for pn, g in subs]

# zoom window on the subdivisions, for the inset
ee = [g.GetEnvelope() for _, g in subs]
zx0 = min(e[0] for e in ee) - 200; zx1 = max(e[1] for e in ee) + 200
zy0 = min(e[2] for e in ee) - 200; zy1 = max(e[3] for e in ee) + 200


_LABEL_SLOT = {}


def draw_subs(ax, lw=1.9, label=False):
    for pn, g in subs:
        for k in range(max(1, g.GetGeometryCount())):
            gg = g.GetGeometryRef(k) if g.GetGeometryCount() else g
            ring = gg.GetGeometryRef(0) if gg.GetGeometryCount() else gg
            p = np.array(ring.GetPoints())
            if p.size:
                ax.plot(p[:, 0], p[:, 1], c="#111", lw=lw + 0.9, zorder=6)
                ax.plot(p[:, 0], p[:, 1], c="#ffd400", lw=lw, zorder=6)
        if label:
            e = g.GetEnvelope()
            # The three parcels are adjacent and their centroids project onto
            # nearly the same point, so a single label row overplots itself.
            dy = 45 + 78 * (_LABEL_SLOT.setdefault(str(pn), len(_LABEL_SLOT)) % 3)
            ax.annotate(str(pn), xy=((e[0] + e[1]) / 2, e[3]),
                        xytext=((e[0] + e[1]) / 2, e[3] + dy),
                        ha="center", fontsize=9, weight="bold", zorder=10,
                        arrowprops=dict(arrowstyle="-", lw=.8, color="k"),
                        bbox=dict(fc="#ffd400", ec="k", lw=.7, pad=1.8))


print("run %s   %.0f m cells   peak %.2f ft" % (RUN, DXM, a.max()))
for thr, col, title in BANDS:
    m = a > thr
    tot_ac = m.sum() * acre
    nstr = int((sd_ft > thr).sum())

    fig = plt.figure(figsize=(15, 11), dpi=150)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.30, 1],
                          width_ratios=[1.05, 1], hspace=0.24, wspace=0.10)

    # ---- main panel: whole domain ------------------------------------
    ax = fig.add_subplot(gs[0, :])
    ax.imshow(hs, cmap="gray", extent=ext, alpha=.92, zorder=1)
    ax.imshow(np.ma.masked_where(~m, m.astype(float)),
              cmap=matplotlib.colors.ListedColormap([col]),
              extent=ext, alpha=.85, zorder=2)
    draw_subs(ax)
    rx0, rx1 = max(zx0, ext[0]), min(zx1, ext[1])
    ry0, ry1 = max(zy0, ext[2]), min(zy1, ext[3])
    ax.add_patch(Rectangle((rx0, ry0), rx1 - rx0, ry1 - ry0, fill=False,
                           ec="#e8112d", lw=2.0, ls="--", zorder=11))
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect("equal")
    ax.set_title("C6 Coyote Canyon — %s\n2D-DWE, %.0f m grid, 100-yr 6-hr transient "
                 "hydrograph — PRELIMINARY" % (title, DXM), fontsize=12.5)
    ax.set_xlabel("UTM Zone 12N Easting (m)"); ax.set_ylabel("Northing (m)")
    bl = ext[1] - 1609.34 - 1450
    ax.plot([bl, bl + 1609.34], [ext[2] + 130, ext[2] + 130], c="k", lw=4,
            zorder=9, solid_capstyle="butt")
    ax.text(bl + 804, ext[2] + 185, "1 mile", ha="center", fontsize=9, zorder=9)
    ax.annotate("", xy=(ext[1] - 230, ext[3] - 180), xytext=(ext[1] - 230, ext[3] - 480),
                arrowprops=dict(arrowstyle="-|>", lw=2, color="k"), zorder=9)
    ax.text(ext[1] - 230, ext[3] - 140, "N", ha="center", fontsize=11,
            weight="bold", zorder=9)
    ax.legend(handles=[Patch(fc=col, ec="none", label="Inundated above threshold"),
                       Patch(fc="none", ec="#ffd400", lw=2, label="NHA subdivision"),
                       Line2D([], [], color="#e8112d", ls="--", lw=1.8,
                              label="Detail below")],
              loc="lower left", fontsize=9.5, framealpha=.93)

    # ---- inset: the subdivisions -------------------------------------
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.imshow(hs, cmap="gray", extent=ext, alpha=.92, zorder=1)
    ax2.imshow(np.ma.masked_where(~m, m.astype(float)),
               cmap=matplotlib.colors.ListedColormap([col]),
               extent=ext, alpha=.85, zorder=2)
    draw_subs(ax2, lw=2.2, label=True)
    dry = sd_ft <= thr
    ax2.scatter(sx[dry], sy[dry], s=24, c="w", ec="#333", lw=.7, zorder=7)
    if (~dry).any():
        ax2.scatter(sx[~dry], sy[~dry], s=95, c="#e8112d", ec="k", lw=1.1,
                    marker="s", zorder=8)
    ax2.set_xlim(zx0, zx1); ax2.set_ylim(zy0, zy1); ax2.set_aspect("equal")
    ax2.set_xlabel("UTM Zone 12N Easting (m)"); ax2.set_ylabel("Northing (m)")
    bl2 = zx0 + 70
    ax2.plot([bl2, bl2 + 304.8], [zy0 + 55, zy0 + 55], c="k", lw=4, zorder=9,
             solid_capstyle="butt")
    ax2.text(bl2 + 152, zy0 + 78, "1,000 ft", ha="center", fontsize=8.5, zorder=9)
    ax2.legend(handles=[Line2D([], [], marker="o", ls="", mfc="w", mec="#333",
                               label="Structure below threshold (%d)" % int(dry.sum())),
                        Line2D([], [], marker="s", ls="", mfc="#e8112d", mec="k",
                               label="Structure above threshold (%d)" % nstr)],
               loc="lower left", fontsize=9, framealpha=.93)

    ax2.set_title("Detail — NHA subdivisions", fontsize=11)

    ax3 = fig.add_subplot(gs[1, 1]); ax3.axis("off")
    lines = [("Subdivision", "Acres", "Ac above", "% area"), None]
    for pn, msk_, g in submask:
        ac_tot = msk_.sum() * acre
        ac_thr = (a[msk_] > thr).sum() * acre
        pct = 100.0 * ac_thr / ac_tot if ac_tot else 0.0
        lines.append((str(pn), "%.0f" % ac_tot, "%.1f" % ac_thr, "%.1f%%" % pct))
    y = 0.95
    for row in lines:
        if row is None:
            ax3.plot([0.02, 0.98], [y + 0.012, y + 0.012], c="k", lw=1,
                     transform=ax3.transAxes); y -= 0.035; continue
        bold = row[0] == "Subdivision"
        for xx, txt, ha in zip((0.02, 0.52, 0.72, 0.95), row, ("left",)*1 + ("right",)*3):
            ax3.text(xx, y, txt, transform=ax3.transAxes, fontsize=10.5, ha=ha,
                     weight="bold" if bold else "normal")
        y -= 0.075
    y -= 0.04
    ax3.text(0.02, y, "Structures above threshold: %d of %d" % (nstr, len(sx)),
             transform=ax3.transAxes, fontsize=11, weight="bold")
    y -= 0.085
    ax3.text(0.02, y, "Domain-wide area above threshold:", transform=ax3.transAxes,
             fontsize=10.5)
    y -= 0.070
    ax3.text(0.02, y, "%.0f acres  (%.2f sq mi)" % (tot_ac, tot_ac / 640.0),
             transform=ax3.transAxes, fontsize=13, weight="bold")
    y -= 0.105
    ax3.text(0.02, y, "Peak depth anywhere in domain: %.1f ft" % a.max(),
             transform=ax3.transAxes, fontsize=10.5)
    y -= 0.070
    ax3.text(0.02, y, "Grid %.0f m  |  9-hour transient run" % DXM,
             transform=ax3.transAxes, fontsize=10.5)

    fig.text(0.5, -0.005,
             "Total area above threshold: %.0f acres (%.2f sq mi).  Depths are on the "
             "GROUND; structures are not in the DEM and no finished-floor elevations "
             "are available, so a structure above the threshold is not a flooded building."
             % (tot_ac, tot_ac / 640.0), ha="center", fontsize=9.5, style="italic")

    out = os.path.join(RUN, "c6_band_%dft.png" % int(thr))
    plt.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("  > %.0f ft : %8.0f acres (%.2f sq mi), %d/%d structures  -> %s"
          % (thr, tot_ac, tot_ac / 640.0, nstr, len(sx), out))
