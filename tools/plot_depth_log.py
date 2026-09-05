#!/usr/bin/env python3
"""
Maximum depth on a LOG colour scale, from a hair above zero, over the
subdivisions and the boundary berm.

    python3 tools/plot_depth_log.py out_c6_1m_berm_6hr [title]

The standard maps mask depths below 0.5 ft, which is the right threshold for
hazard but hides the question "is ANY water getting in, and from where?".  A log
scale spans 0.01 ft to the domain peak in one figure, so a hundredth of a foot
inside the berm and twenty feet in the wash are both legible.

Water inside the subdivisions is split into two classes, because they mean
different things and cannot be told apart by depth alone:
  * CONNECTED  -- joined to the exterior flood by a continuous wet path, so it
                  entered from outside;
  * ISOLATED   -- a separate wet region, i.e. direct rainfall ponded in place,
                  which no perimeter berm can prevent.
"""
import sys, os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LightSource
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import ndimage
from osgeo import gdal, ogr
gdal.UseExceptions(); ogr.UseExceptions()

FT   = 3.280839895
RUN  = sys.argv[1] if len(sys.argv) > 1 else "out_c6_1m_berm_6hr"
TITLE= sys.argv[2] if len(sys.argv) > 2 else "C6 Coyote Canyon — with boundary berm"
SITE = "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097"
CONN_TH = 0.05      # ft; connectivity judged above a film thickness

d0 = gdal.Open(SITE + "/outputs_2D/cliped_dem_berm.tif")
gt = d0.GetGeoTransform(); NX, NY = d0.RasterXSize, d0.RasterYSize
Z  = d0.GetRasterBand(1).ReadAsArray().astype("f8")
ds = gdal.Open(os.path.join(RUN, "depth_max.tif")); bb = ds.GetRasterBand(1)
D  = bb.ReadAsArray().astype("f8"); nd = bb.GetNoDataValue()
D  = np.where((D == nd) | ~np.isfinite(D), 0.0, D) * FT

src = ogr.Open(SITE + "/outputs/05_Subdivisions_clip.gpkg"); L = src.GetLayer(0)
mem = ogr.GetDriverByName("Memory").CreateDataSource("m")
lyr = mem.CreateLayer("s", L.GetSpatialRef(), ogr.wkbPolygon)
polys = []
for f in L:
    g = f.GetGeometryRef(); e = g.GetEnvelope()
    if f.GetField("chapter") == "COYOTE CANYON" and 706000 < e[0] < 708000 and 3959000 < e[2] < 3961000:
        polys.append((f.GetField("proj_num"), g.Clone()))
        ft_ = ogr.Feature(lyr.GetLayerDefn()); ft_.SetGeometry(g.Clone()); lyr.CreateFeature(ft_)
mr = gdal.GetDriverByName("MEM").Create("", NX, NY, 1, gdal.GDT_Byte)
mr.SetGeoTransform(gt); mr.SetProjection(d0.GetProjection())
gdal.RasterizeLayer(mr, [1], lyr, burn_values=[1])
SUB = mr.GetRasterBand(1).ReadAsArray().astype(bool)
rr = gdal.Open(SITE + "/outputs_2D/berm_raise.tif")
BERM = rr.GetRasterBand(1).ReadAsArray() > 0.05

W = D > CONN_TH
lab, n = ndimage.label(W, structure=np.ones((3, 3)))
sz = ndimage.sum(W, lab, range(1, n + 1)); main = int(np.argmax(sz)) + 1
CONN = (lab == main) & SUB
ISO  = W & SUB & ~CONN

jj, ii = np.where(BERM | SUB)
pad = 90
i0, i1 = max(ii.min() - pad, 0), min(ii.max() + pad, NX)
j0, j1 = max(jj.min() - pad, 0), min(jj.max() + pad, NY)
ext = [gt[0] + i0, gt[0] + i1, gt[3] - j1, gt[3] - j0]
sl = (slice(j0, j1), slice(i0, i1))

fig, ax = plt.subplots(figsize=(13.5, 10.5), dpi=160)
hs = LightSource(315, 45).hillshade(Z[sl], vert_exag=3, dx=1, dy=1)
ax.imshow(hs, cmap="gray", extent=ext, alpha=0.95, zorder=1)
VMIN, VMAX = 0.01, max(D[sl].max(), 1.0)
im = ax.imshow(np.ma.masked_less_equal(D[sl], VMIN), cmap="turbo",
               norm=LogNorm(vmin=VMIN, vmax=VMAX), extent=ext, alpha=0.88, zorder=2)
ax.contour(np.linspace(ext[0], ext[1], i1 - i0), np.linspace(ext[3], ext[2], j1 - j0),
           BERM[sl].astype(float), levels=[0.5], colors="#000000", linewidths=2.0, zorder=5)
for pn, g in polys:
    for k in range(max(1, g.GetGeometryCount())):
        gg = g.GetGeometryRef(k) if g.GetGeometryCount() else g
        ring = gg.GetGeometryRef(0) if gg.GetGeometryCount() else gg
        p = np.array(ring.GetPoints())
        if p.size:
            ax.plot(p[:, 0], p[:, 1], c="#111", lw=2.6, zorder=6)
            ax.plot(p[:, 0], p[:, 1], c="#ffd400", lw=1.3, zorder=6)
ax.contour(np.linspace(ext[0], ext[1], i1 - i0), np.linspace(ext[3], ext[2], j1 - j0),
           CONN[sl].astype(float), levels=[0.5], colors="#ff00ff", linewidths=1.6, zorder=7)
# mark the berm's open north-east end
bj, bi = np.where(BERM); ne_j = int(bj.min()); ne_i = int(bi[bj.argmin()])
nex, ney = gt[0] + ne_i, gt[3] - ne_j
ax.annotate("berm ends here —\nwater enters around it",
            xy=(nex, ney), xytext=(-235, -58), textcoords="offset points",
            fontsize=11, fontweight="bold", color="k", zorder=10,
            arrowprops=dict(arrowstyle="-|>", lw=2.0, color="#ff00ff"),
            bbox=dict(fc="#ffd400", ec="k", lw=0.8, alpha=0.95, pad=3))
cb = fig.colorbar(im, ax=ax, shrink=0.86, pad=0.012)
cb.set_label("Maximum depth (ft, log scale)", fontsize=10.5)
sx, sy = ext[0] + 60, ext[2] + 55
ax.plot([sx, sx + 152.4], [sy, sy], c="k", lw=4.5, zorder=9, solid_capstyle="butt")
ax.text(sx + 76, sy + 22, "500 ft", ha="center", fontsize=9.5, fontweight="bold", zorder=9)
ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("%s\nmaximum depth from 0.01 ft — log scale" % TITLE,
             fontsize=13, fontweight="bold")
ax.legend(handles=[
    Line2D([], [], color="#000000", lw=2.0, label="Berm"),
    Line2D([], [], color="#ffd400", lw=2.0, label="NHA subdivision boundary"),
    Line2D([], [], color="#ff00ff", lw=1.8,
           label="Water connected to the exterior flood (%.2f ac)" % (CONN.sum() / 4046.86)),
    Patch(fc="none", ec="none",
          label="Isolated ponding, direct rain (%.2f ac)" % (ISO.sum() / 4046.86)),
], loc="lower right", fontsize=9, framealpha=0.94)
out = os.path.join(RUN, "c6_depth_log.png")
plt.savefig(out, bbox_inches="tight")
print("wrote", out)
print("interior water above %.2f ft: connected %.3f ac | isolated %.3f ac | max %.2f ft"
      % (CONN_TH, CONN.sum() / 4046.86, ISO.sum() / 4046.86, D[SUB].max()))
