#!/usr/bin/env python3
"""
The existing AECOM 1%-annual-chance floodplain around the C10 community.

    python3 tools/plot_aecom_c10.py

Section 3.1.1 states in words that the AECOM programme-scale mapping does not
reach this community. This draws it, because the reader's question is not
whether the statement is true but HOW FAR AWAY the mapping stops and in which
direction -- which is what decides whether the absence is a near miss or a
different drainage entirely.

DISTANCES ARE MEASURED, NOT QUOTED. The nearest-polygon distance annotated on
the map is computed here from the master AECOM layer against the subdivision
polygons, so the figure cannot drift from the layer it is drawn from.

Matplotlib, not QGIS: the QGIS map scripts in PythonScripts/maps hard-code a
Windows QGIS_PREFIX_PATH, and when that is wrong the layout silently renders
every label as a missing-glyph box and drops any online basemap. This figure is
offline and self-contained.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from osgeo import gdal, ogr, osr

gdal.UseExceptions(); ogr.UseExceptions()

NHA  = "/home/arash/Dropbox/Chloeta/NHA"
SITE = NHA + "/WS3_GIS/NM15-141Local"
C6   = NHA + "/WS3_GIS/NM15-097/outputs/"
AEC  = NHA + "/WS3_GIS/dwelling_unit_inventory/floodplains.geojson"
DEM  = SITE + "/demlr/cliped_utm.tif"
OUT  = NHA + "/WS3_GIS/NM15-141Local/Report_for_C10/Figures/NM15-141_aecom_floodplain.png"
EPSG = 26912
# A FIXED window on the community, not one sized to the AECOM polygons. Those
# polygons are long sinuous corridors tens of kilometres end to end, so letting
# their envelopes set the extent pulls the frame out to ~30 km, past the edge of
# the DEM, and shrinks the subdivisions to a speck -- which loses the one thing
# the figure is for. 5.4 km square, centred on the community and kept inside the
# terrain tile, holds the community and the nearest mapped floodplain together.
CX, CY = 711400.0, 3959900.0
HALF   = 2700.0


def to_utm(path, where=None):
    ds = ogr.Open(path); ly = ds.GetLayer()
    sr = ly.GetSpatialRef(); tg = osr.SpatialReference(); tg.ImportFromEPSG(EPSG)
    sr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tg.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = None if sr.IsSame(tg) else osr.CoordinateTransformation(sr, tg)
    out = []
    for f in ly:
        g = f.GetGeometryRef()
        if g is None:
            continue
        g = g.Clone()
        if tr:
            g.Transform(tr)
        out.append((g, f))
    return out


def rings(g):
    """Every exterior ring of a (multi)polygon, as Nx2 arrays."""
    res = []
    if g.GetGeometryName() == "MULTIPOLYGON":
        for k in range(g.GetGeometryCount()):
            res += rings(g.GetGeometryRef(k))
    elif g.GetGeometryName() == "POLYGON":
        r = g.GetGeometryRef(0)
        p = np.array(r.GetPoints())
        if p.size:
            res.append(p[:, :2])
    return res


# ------------------------------------------------------------------ features --
subs = [(f.GetField("proj_num"), g) for g, f in to_utm(C6 + "05_Subdivisions_clip.gpkg")]
subs = [(n, g) for n, g in subs
        if 710000 <= g.GetEnvelope()[0] <= 712500 and 3959500 <= g.GetEnvelope()[2] <= 3961000]
if not subs:
    raise SystemExit("no C10 subdivisions found")

aec_all = [g for g, f in to_utm(AEC)]
# keep what is near enough to be informative; the layer spans the whole Nation
X0, X1 = CX - HALF, CX + HALF
Y0, Y1 = CY - HALF, CY + HALF
win = ogr.Geometry(ogr.wkbLinearRing)
for x, y in [(X0, Y0), (X1, Y0), (X1, Y1), (X0, Y1), (X0, Y0)]:
    win.AddPoint(x, y)
WIN = ogr.Geometry(ogr.wkbPolygon); WIN.AddGeometry(win)

near = []
for g in aec_all:
    d = min(g.Distance(s) for _, s in subs)
    c = g.Intersection(WIN)
    if c is not None and not c.IsEmpty() and c.GetArea() > 0:
        near.append((d, c.Clone()))
near.sort(key=lambda t: t[0])
DMIN, GMIN = near[0]
print("AECOM polygons intersecting the window: %d ; nearest %.0f m (%.2f km)" % (len(near), DMIN, DMIN / 1000))
for n, g in subs:
    print("   %-18s nearest AECOM %.0f m" % (n, min(gg.Distance(g) for _, gg in near)))

# ------------------------------------------------------------------- terrain --
W = int((X1 - X0) / 8.0)
t = gdal.Warp("", DEM, format="MEM", outputBounds=(X0, Y0, X1, Y1),
              width=W, height=int(W * (Y1 - Y0) / (X1 - X0)), resampleAlg="average")
z = t.GetRasterBand(1).ReadAsArray().astype(float)
nd = t.GetRasterBand(1).GetNoDataValue()
if nd is not None:
    z = np.where(z == nd, np.nan, z)
z = np.where(np.isfinite(z), z, np.nanmin(z))
hs = LightSource(azdeg=315, altdeg=45).hillshade(z, vert_exag=2.2, dx=8, dy=8)

# -------------------------------------------------------------------- figure --
fig, ax = plt.subplots(figsize=(13.5, 13.5 * (Y1 - Y0) / (X1 - X0) + 1.1), dpi=170)
ax.imshow(hs, cmap="gray", extent=[X0, X1, Y0, Y1], origin="upper", alpha=.92, zorder=1)

for g, f in to_utm(SITE + "/roads/tl_2025_35031_roads/tl_2025_35031_roads.shp"):
    e = g.GetEnvelope()
    if e[1] < X0 or e[0] > X1 or e[3] < Y0 or e[2] > Y1:
        continue
    for k in range(max(1, g.GetGeometryCount())):
        gg = g.GetGeometryRef(k) if g.GetGeometryCount() else g
        p = np.array(gg.GetPoints())
        if p.size:
            ax.plot(p[:, 0], p[:, 1], c="#b0b0b0", lw=.8, zorder=3)

for d, g in near:
    for p in rings(g):
        ax.fill(p[:, 0], p[:, 1], fc="#6aa9dc", ec="#1c5fa8", lw=1.4, alpha=.55, zorder=4)

for n, g in subs:
    for p in rings(g):
        ax.plot(p[:, 0], p[:, 1], c="#111", lw=2.6, zorder=6)
        ax.plot(p[:, 0], p[:, 1], c="#ffd400", lw=1.4, zorder=6)

sx, sy = [], []
for g, f in to_utm(C6 + "00_NHA_Structures_clip.gpkg"):
    if X0 <= g.GetX() <= X1 and Y0 <= g.GetY() <= Y1:
        sx.append(g.GetX()); sy.append(g.GetY())
ax.scatter(sx, sy, s=16, c="w", ec="#222", lw=.7, zorder=7)

# the measured gap, drawn between the two closest points
_, gsub = min(subs, key=lambda t: GMIN if False else min(gg.Distance(t[1]) for _, gg in near))
pa, pb = None, None
for p in rings(gsub):
    for x, y in p:
        pt = ogr.Geometry(ogr.wkbPoint); pt.AddPoint(float(x), float(y))
        d = GMIN.Distance(pt)
        if pa is None or d < pa[0]:
            pa = (d, (x, y))
best = None
for p in rings(GMIN):
    for x, y in p:
        pt = ogr.Geometry(ogr.wkbPoint); pt.AddPoint(float(x), float(y))
        d = pt.Distance(gsub)
        if best is None or d < best[0]:
            best = (d, (x, y))
pb = best
ax.annotate("", xy=pb[1], xytext=pa[1],
            arrowprops=dict(arrowstyle="<|-|>", lw=2.0, color="#d1373a",
                            shrinkA=0, shrinkB=0), zorder=9)
mx, my = (pa[1][0] + pb[1][0]) / 2, (pa[1][1] + pb[1][1]) / 2
ax.text(mx, my + 0.018 * (Y1 - Y0),
        "%.2f km to the nearest\nAECOM 1%% floodplain" % (DMIN / 1000),
        ha="center", va="bottom", fontsize=10.5, weight="bold", color="#8a1015", zorder=10,
        bbox=dict(fc="w", ec="#d1373a", lw=1.0, alpha=.92, pad=2.6))

ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_aspect("equal")
ax.set_xlabel("UTM Zone 12N Easting (m)"); ax.set_ylabel("Northing (m)")
ax.set_title("C10 — the existing AECOM 1%-annual-chance floodplain and the NM15-141 community\n"
             "AECOM (formerly URS), Final Report, Floodplain Analysis, Phase 3, 1 February 2016 "
             "— context only", fontsize=12.5)
ax.legend(handles=[Patch(fc="#6aa9dc", ec="#1c5fa8", alpha=.6, label="AECOM 1%-annual-chance floodplain"),
                   Patch(fc="none", ec="#ffd400", lw=2, label="NHA subdivision (NM15-141 / NM15-034)"),
                   Line2D([], [], marker="o", ls="", mfc="w", mec="#222", label="NHA structure (%d)" % len(sx)),
                   Line2D([], [], c="#b0b0b0", lw=1.2, label="Road (TIGER 2025)"),
                   Line2D([], [], c="#d1373a", lw=2, label="Measured gap to nearest polygon")],
          loc="upper left", fontsize=9, framealpha=.94)

bl = X0 + 0.06 * (X1 - X0); by = Y0 + 0.045 * (Y1 - Y0)
ax.plot([bl, bl + 1000], [by, by], c="k", lw=4, zorder=9, solid_capstyle="butt")
ax.text(bl + 500, by + 0.012 * (Y1 - Y0), "1 km", ha="center", fontsize=9, zorder=9)
nx_ = X1 - 0.07 * (X1 - X0)
ax.annotate("", xy=(nx_, Y1 - 0.045 * (Y1 - Y0)), xytext=(nx_, Y1 - 0.115 * (Y1 - Y0)),
            arrowprops=dict(arrowstyle="-|>", lw=2, color="k"), zorder=9)
ax.text(nx_, Y1 - 0.03 * (Y1 - Y0), "N", ha="center", fontsize=11, weight="bold", zorder=9)

fig.tight_layout()
fig.savefig(OUT, bbox_inches="tight")
print("wrote", OUT)
