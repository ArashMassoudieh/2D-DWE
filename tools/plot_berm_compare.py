#!/usr/bin/env python3
"""
Baseline against berm, maximum depth on a LOG colour scale, over the subdivisions.

    python3 tools/plot_berm_compare.py

A log scale is used because the two panels differ by three orders of magnitude
inside the subdivisions -- feet in the baseline, hundredths of a foot with the
berm -- and a linear scale renders the berm panel uniformly blank, which reads
as "no data" rather than "no water".
"""
import os, sys
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LightSource
from matplotlib.lines import Line2D
from scipy import ndimage
from osgeo import gdal, ogr
gdal.UseExceptions(); ogr.UseExceptions()

FT = 3.280839895
SITE = "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097"
DWE  = "/home/arash/Projects/2D-DWE"
STORM = (sys.argv[1] if len(sys.argv) > 1 else "6hr").lower()
SUF   = "" if STORM == "6hr" else "_24hr"
OUT   = SITE + "/Report_for_C6/Figures/NM15-097_berm_depth_log%s.png" % SUF
_RUNS = {"6hr":  ("out_c6_1m_rain",     "out_c6_1m_berm_6hr"),
         "24hr": ("out_c6_1m_24hr_v2",  "out_c6_1m_berm_24hr")}[STORM]
PANES = [("Existing conditions", _RUNS[0], SITE + "/outputs_2D/cliped_dem.tif"),
         ("With the boundary berm", _RUNS[1], SITE + "/outputs_2D/cliped_dem_berm.tif")]
CONN_TH = 0.05

g0 = gdal.Open(SITE + "/outputs_2D/cliped_dem.tif"); gt = g0.GetGeoTransform()
NX, NY = g0.RasterXSize, g0.RasterYSize
src = ogr.Open(SITE + "/outputs/05_Subdivisions_clip.gpkg"); L = src.GetLayer(0)
mem = ogr.GetDriverByName("Memory").CreateDataSource("m")
lyr = mem.CreateLayer("s", L.GetSpatialRef(), ogr.wkbPolygon); polys = []
for f in L:
    g = f.GetGeometryRef(); e = g.GetEnvelope()
    if f.GetField("chapter") == "COYOTE CANYON" and 706000 < e[0] < 708000 and 3959000 < e[2] < 3961000:
        polys.append(g.Clone())
        ft_ = ogr.Feature(lyr.GetLayerDefn()); ft_.SetGeometry(g.Clone()); lyr.CreateFeature(ft_)
mr = gdal.GetDriverByName("MEM").Create("", NX, NY, 1, gdal.GDT_Byte)
mr.SetGeoTransform(gt); mr.SetProjection(g0.GetProjection())
gdal.RasterizeLayer(mr, [1], lyr, burn_values=[1])
SUB = mr.GetRasterBand(1).ReadAsArray().astype(bool)
rr = gdal.Open(SITE + "/outputs_2D/berm_raise.tif")
BERM = rr.GetRasterBand(1).ReadAsArray() > 0.05

jj, ii = np.where(BERM | SUB); pad = 80
i0, i1 = max(ii.min() - pad, 0), min(ii.max() + pad, NX)
j0, j1 = max(jj.min() - pad, 0), min(jj.max() + pad, NY)
ext = [gt[0] + i0, gt[0] + i1, gt[3] - j1, gt[3] - j0]
sl = (slice(j0, j1), slice(i0, i1))
xs = np.linspace(ext[0], ext[1], i1 - i0); ys = np.linspace(ext[3], ext[2], j1 - j0)

def depth(run):
    ds = gdal.Open(os.path.join(DWE, run, "depth_max.tif")); b = ds.GetRasterBand(1)
    a = b.ReadAsArray().astype("f8"); nd = b.GetNoDataValue()
    return np.where((a == nd) | ~np.isfinite(a), 0.0, a) * FT

VMIN, VMAX = 0.01, 25.0
fig, axes = plt.subplots(1, 2, figsize=(18.0, 9.6), dpi=160)
stats = []
for ax, (ttl, run, demp) in zip(axes, PANES):
    D = depth(run)
    _zd = gdal.Open(demp)          # hold the dataset; a bare band is freed
    Z = _zd.GetRasterBand(1).ReadAsArray().astype("f8")
    ax.imshow(LightSource(315, 45).hillshade(Z[sl], vert_exag=3, dx=1, dy=1),
              cmap="gray", extent=ext, alpha=0.95, zorder=1)
    im = ax.imshow(np.ma.masked_less_equal(D[sl], VMIN), cmap="turbo",
                   norm=LogNorm(vmin=VMIN, vmax=VMAX), extent=ext, alpha=0.88, zorder=2)
    if "berm" in run:
        ax.contour(xs, ys, BERM[sl].astype(float), levels=[0.5], colors="k", linewidths=2.0, zorder=5)
    for g in polys:
        for k in range(max(1, g.GetGeometryCount())):
            gg = g.GetGeometryRef(k) if g.GetGeometryCount() else g
            ring = gg.GetGeometryRef(0) if gg.GetGeometryCount() else gg
            p = np.array(ring.GetPoints())
            if p.size:
                ax.plot(p[:, 0], p[:, 1], c="#111", lw=2.4, zorder=6)
                ax.plot(p[:, 0], p[:, 1], c="#ffd400", lw=1.2, zorder=6)
    W = D > CONN_TH
    lab, n = ndimage.label(W, structure=np.ones((3, 3)))
    szs = ndimage.sum(W, lab, range(1, n + 1)); main = int(np.argmax(szs)) + 1
    conn = (lab == main) & SUB
    ax.contour(xs, ys, conn[sl].astype(float), levels=[0.5], colors="#ff00ff", linewidths=1.5, zorder=7)
    d = D[SUB]
    stats.append((ttl, (d > 1).sum() / 4046.86, d.max(), conn.sum() / 4046.86))
    ax.set_title("%s\n%.2f ac over 1 ft   |   max %.2f ft   |   connected inflow %.2f ac"
                 % (ttl, stats[-1][1], stats[-1][2], stats[-1][3]),
                 fontsize=11.5, fontweight="bold")
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    sx, sy = ext[0] + 55, ext[2] + 50
    ax.plot([sx, sx + 152.4], [sy, sy], c="k", lw=4.5, zorder=9, solid_capstyle="butt")
    ax.text(sx + 76, sy + 20, "500 ft", ha="center", fontsize=9.5, fontweight="bold", zorder=9)
cb = fig.colorbar(im, ax=axes, shrink=0.72, pad=0.012)
cb.set_label("Maximum depth (ft, log scale)", fontsize=10.5)
fig.legend(handles=[Line2D([], [], color="k", lw=2.0, label="Boundary berm"),
                    Line2D([], [], color="#ffd400", lw=2.0, label="NHA subdivision boundary"),
                    Line2D([], [], color="#ff00ff", lw=1.6,
                           label="Interior water connected to the exterior flood")],
           loc="lower center", ncol=3, fontsize=9.5, bbox_to_anchor=(0.5, 0.015))
fig.suptitle("C6 Coyote Canyon — maximum depth, 100-year %s storm" % ("6-hour" if STORM=="6hr" else "24-hour"),
             fontsize=13.5, fontweight="bold", y=0.965)
plt.savefig(OUT, bbox_inches="tight")
print("wrote", OUT)
for t, a1, mx, cn in stats:
    print("  %-24s  %6.2f ac >1 ft | max %5.2f ft | connected %5.3f ac" % (t, a1, mx, cn))
