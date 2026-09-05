#!/usr/bin/env python3
"""
Structures classified by freeboard above the 2-D water surface, against the two
NHA Drainage Design Manual criteria, for the 6-hour and 24-hour 100-year storms.

    python3 tools/plot_freeboard.py

Freeboard values are read from the table generator's output, not recomputed, so
this map, the report tables and the workbook cannot drift apart:

    WS3_GIS/NM15-097/outputs_2D/structure_freeboard.csv
      <- PythonScripts/maps/C6-table-structure-freeboard.py

Classes (Manual Sec. 3.4 pad >= WSEL + 1 ft; Sec. 3.4 / 5.2.2 floor >= WSEL + 2 ft):
    < 0 ft     water surface stands ABOVE the modelled ground
    0 - 1 ft   below the pad criterion
    1 - 2 ft   meets the pad criterion, below the finished-floor criterion
    >= 2 ft    meets both

Freeboard is measured to BARE-EARTH ground, not to a surveyed pad and not to a
finished floor, so the classes are a lower bound.  The map is not a finding of
non-compliance -- see report Section 4.4.6.
"""
import csv, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, LightSource
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from osgeo import gdal, ogr
gdal.UseExceptions(); ogr.UseExceptions()

FT   = 3.280839895
SITE = "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097"
DWE  = "/home/arash/Projects/2D-DWE"
SCEN = (sys.argv[1] if len(sys.argv) > 1 else "baseline").lower()
SUF  = "" if SCEN == "baseline" else "_berm"
CSVP = os.path.join(SITE, "outputs_2D/structure_freeboard%s.csv" % SUF)
DEMP = os.path.join(SITE, "outputs_2D/cliped_dem.tif")
SUBP = os.path.join(SITE, "outputs/05_Subdivisions_clip.gpkg")
OUTP = os.path.join(SITE, "Report_for_C6/Figures/NM15-097_2d_structure_freeboard%s.png" % SUF)
_RUN = {"baseline": ("out_c6_1m_rain", "out_c6_1m_24hr_v2"),
        "berm":     ("out_c6_1m_berm_6hr", "out_c6_1m_berm_24hr")}[SCEN]
def _dep(run):
    for f in ("depth_max.tif", "depth_max_snap.tif"):
        q = os.path.join(DWE, run, f)
        if os.path.exists(q): return q
    return None
RUNS = [(k, _dep(r), t) for k, r, t in
        (("6hr", _RUN[0], "100-year 6-hour"), ("24hr", _RUN[1], "100-year 24-hour"))]
RUNS = [r for r in RUNS if r[1] is not None]

# class edges in ft, colours, labels
EDGES  = [-99.0, 0.0, 1.0, 2.0, 99.0]
COLORS = ["#d7191c", "#fdae61", "#ffe14d", "#1a9641"]
LABELS = ["Freeboard $<$ 0 ft  (water above ground)",
          "0 – 1 ft  (below the 1 ft pad criterion)",
          "1 – 2 ft  (meets pad, below the 2 ft floor criterion)",
          "$\\geq$ 2 ft  (meets both criteria)"]

rows = list(csv.DictReader(open(CSVP)))
for r in rows:
    for k in ("x", "y") + tuple(k + "_fb" for k, _, _ in RUNS):
        r[k] = float(r[k])
print("structures read: %d" % len(rows))

# subdivision polygons in the model area
sds = ogr.Open(SUBP); sly = sds.GetLayer()
polys = []
for f in sly:
    g = f.GetGeometryRef()
    if g is None:
        continue
    e = g.GetEnvelope()
    if 706000 < e[0] < 708000 and 3959000 < e[2] < 3961000:
        polys.append((f.GetField("proj_num"), g.Clone()))
env = [p[1].GetEnvelope() for p in polys]
PAD = 120.0
bx0 = min(e[0] for e in env) - PAD; bx1 = max(e[1] for e in env) + PAD
by0 = min(e[2] for e in env) - PAD; by1 = max(e[3] for e in env) + PAD

dem = gdal.Open(DEMP); gt = dem.GetGeoTransform()
nx, ny = dem.RasterXSize, dem.RasterYSize
ext = [gt[0], gt[0] + nx * gt[1], gt[3] - ny * gt[1], gt[3]]
zb = dem.GetRasterBand(1).ReadAsArray().astype(float)
hs = LightSource(azdeg=315, altdeg=45).hillshade(zb, vert_exag=3, dx=1.0, dy=1.0)

fig, axes = plt.subplots(1, len(RUNS), figsize=(8.9*len(RUNS), 7.4), dpi=150, squeeze=False)
axes = axes[0]
lv = [0.5, 1, 2, 3, 4, 6]
dcm = ListedColormap(["#cfe8f7", "#93cbe9", "#4f9ed2", "#2166ac", "#0d2f6e"])

for ax, (key, dpath, sub) in zip(axes, RUNS):
    ds = gdal.Open(dpath); bnd = ds.GetRasterBand(1)
    a = bnd.ReadAsArray().astype(float); nd = bnd.GetNoDataValue()
    a = np.where((a == nd) | ~np.isfinite(a), 0.0, a) * FT
    ax.imshow(hs, cmap="gray", extent=ext, alpha=0.9, zorder=1)
    im = ax.imshow(np.ma.masked_less(a, 0.5), cmap=dcm,
                   norm=BoundaryNorm(lv, dcm.N), extent=ext, alpha=0.80, zorder=2)
    for pn, g in polys:
        for k in range(max(1, g.GetGeometryCount())):
            gg = g.GetGeometryRef(k) if g.GetGeometryCount() else g
            ring = gg.GetGeometryRef(0) if gg.GetGeometryCount() else gg
            p = np.array(ring.GetPoints())
            if p.size:
                ax.plot(p[:, 0], p[:, 1], c="#111", lw=2.6, zorder=6)
                ax.plot(p[:, 0], p[:, 1], c="#ffd400", lw=1.3, zorder=6)

    # Anchor each label to its own polygon and stagger the callouts down the left
    # margin: the three parcels overlap in x, so centred labels land on top of
    # one another.
    _used = []
    for pn, g in polys:
        try:
            rp = g.PointOnSurface(); axp, ayp = rp.GetX(), rp.GetY()
        except Exception:
            e = g.GetEnvelope(); axp, ayp = (e[0] + e[1]) / 2, (e[2] + e[3]) / 2
        lx = bx0 + 0.13 * (bx1 - bx0)
        ly = ayp
        span = by1 - by0
        for prev in _used:
            if abs(ly - prev) < 0.085 * span:
                ly = prev - 0.085 * span
        _used.append(ly)
        ax.annotate(str(pn), xy=(axp, ayp), xytext=(lx, ly),
                    ha="center", va="center", fontsize=9.5, weight="bold", zorder=11,
                    arrowprops=dict(arrowstyle="-|>", lw=1.2, color="k",
                                    shrinkA=0, shrinkB=2),
                    bbox=dict(fc="#ffd400", ec="k", lw=0.7, pad=2.0))
    fb = np.array([r[key + "_fb"] for r in rows])
    xs = np.array([r["x"] for r in rows]); ys = np.array([r["y"] for r in rows])
    counts = []
    for c in range(4):
        m = (fb >= EDGES[c]) & (fb < EDGES[c + 1])
        counts.append(int(m.sum()))
        if m.any():
            ax.scatter(xs[m], ys[m], s=78, c=COLORS[c], ec="k", lw=0.9,
                       marker="o", zorder=8 + (3 - c))
    ax.set_xlim(bx0, bx1); ax.set_ylim(by0, by1); ax.set_aspect("equal")
    ax.set_title(sub, fontsize=11.5, weight="bold")
    ax.set_xlabel("UTM Zone 12N Easting (m)", fontsize=9)
    ax.tick_params(labelsize=8)
    bl = bx0 + 60
    ax.plot([bl, bl + 152.4], [by0 + 45, by0 + 45], c="k", lw=4.5,
            zorder=9, solid_capstyle="butt")
    ax.text(bl + 76, by0 + 62, "500 ft", ha="center", fontsize=9, weight="bold", zorder=9)
    ax.annotate("N", xy=(bx1 - 55, by1 - 55), fontsize=12, weight="bold",
                ha="center", zorder=9)
    ax.annotate("", xy=(bx1 - 55, by1 - 48), xytext=(bx1 - 55, by1 - 118),
                arrowprops=dict(fc="k", ec="k", width=2.0, headwidth=7.5), zorder=9)
    ax.set_label(key)
    print("  %-5s classes  <0: %2d | 0-1: %2d | 1-2: %2d | >=2: %2d"
          % (key, counts[0], counts[1], counts[2], counts[3]))
    if key == "6hr":
        axes[0].set_ylabel("Northing (m)", fontsize=9)
    ax._counts = counts

for _a in axes[1:]: _a.tick_params(labelleft=False)
cb = fig.colorbar(im, ax=axes, shrink=0.72, pad=0.012, ticks=lv, location="right")
cb.set_label("Maximum ground depth (ft)", fontsize=10)

handles = [Line2D([], [], marker="o", ls="", ms=10, mfc=COLORS[c], mec="k", mew=0.9,
                  label="%s\n        %s" % (LABELS[c],
                        "   ".join("%s: %d" % (RUNS[i][0], axes[i]._counts[c])
                                   for i in range(len(RUNS)))))
           for c in range(4)]
handles.append(Patch(fc="none", ec="#ffd400", lw=2, label="NHA subdivision boundary"))
fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=9,
           frameon=True, framealpha=0.95, bbox_to_anchor=(0.5, -0.085))
fig.suptitle("C6 Coyote Canyon — structure freeboard above the 2-D water surface, "
             "against the NHA Manual pad and finished-floor criteria",
             fontsize=13.5, weight="bold", y=1.005)
fig.text(0.5, 0.945,
         "Freeboard measured to BARE-EARTH ground, not to a surveyed pad or finished floor — "
         "a lower bound. Not a finding of non-compliance (report Section 4.4.6).",
         ha="center", fontsize=9.5, style="italic")
plt.savefig(OUTP, bbox_inches="tight")
print("wrote", OUTP)
