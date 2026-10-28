#!/usr/bin/env python3
"""
Perpendicular cross-sections across the C10 existing berm.

    python3 tools/berm_cross_sections.py --fig out.png

The longitudinal profile plots height above a DETRENDED local grade, and the
zero on that plot is not an elevation datum -- it is the straight chord joining
the two toe elevations at -/+ TOE_FT on the same transect. That is a per-station
reference, so it removes the hillside the berm is built on and leaves only the
built ridge. These cross-sections plot the raw elevations with that chord drawn
on top, so the reference can be seen rather than taken on trust.
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from osgeo import gdal, ogr

gdal.UseExceptions(); ogr.UseExceptions()
FT = 3.280839895
S = "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-141Local"

ap = argparse.ArgumentParser()
ap.add_argument("--dem", default=S + "/demhr/cliped_sd.tif")
ap.add_argument("--berm", default=S + "/outputs/existing_berm.shp")
ap.add_argument("--toe-ft", type=float, default=39.0)
ap.add_argument("--n", type=int, default=6, help="sections per line")
ap.add_argument("--fig", default=None)
A = ap.parse_args()

d = gdal.Open(A.dem); GT = d.GetGeoTransform()
Z = d.GetRasterBand(1).ReadAsArray().astype(float); NY, NX = Z.shape
def samp(x, y):
    i = np.clip(((np.asarray(x) - GT[0]) / GT[1]).astype(int), 0, NX - 1)
    j = np.clip(((GT[3] - np.asarray(y)) / abs(GT[5])).astype(int), 0, NY - 1)
    return Z[j, i]

ds = ogr.Open(A.berm); ly = ds.GetLayer()
lines = [np.array(f.GetGeometryRef().GetPoints())[:, :2] for f in ly]
NAMES = ["Line 1 - west boundary", "Line 2 - north boundary"]
T = A.toe_ft / FT
offs = np.arange(-T, T + 1e-9, 0.5)

fig, axes = plt.subplots(2, A.n, figsize=(3.05 * A.n, 7.4), dpi=150)
for r, p in enumerate(lines):
    seg = np.r_[0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))]
    s = np.arange(0, seg[-1], 1.0)
    px = np.interp(s, seg, p[:, 0]); py = np.interp(s, seg, p[:, 1])
    dx = np.gradient(px); dy = np.gradient(py)
    L = np.hypot(dx, dy); ux, uy = -dy / L, dx / L
    stns = np.linspace(8, s[-1] - 8, A.n).astype(int)
    for cnum, k in enumerate(stns):
        ax = axes[r, cnum]
        xs = px[k] + offs * ux[k]; ys = py[k] + offs * uy[k]
        z = samp(xs, ys) * FT
        chord = z[0] + (z[-1] - z[0]) * (offs - offs[0]) / (offs[-1] - offs[0])
        resid = z - chord
        win = np.abs(offs) <= 16.0 / FT
        h = np.nanmax(resid[win]); ko = int(np.argmax(np.where(win, resid, -np.inf)))
        ax.fill_between(offs * FT, chord, z, where=z > chord, color="#c8b89a", zorder=2)
        ax.plot(offs * FT, z, c="#5a3010", lw=1.7, zorder=4)
        ax.plot(offs * FT, chord, c="#1a4f95", lw=1.4, ls="--", zorder=3)
        ax.axvline(0, c="#888", lw=.8, ls=":", zorder=1)
        ax.annotate("", xy=(offs[ko] * FT, z[ko]), xytext=(offs[ko] * FT, chord[ko]),
                    arrowprops=dict(arrowstyle="<->", lw=1.4, color="#e8112d"), zorder=6)
        ax.text(offs[ko] * FT + 2, (z[ko] + chord[ko]) / 2, "%.2f ft" % h,
                fontsize=8.5, color="#8a1015", weight="bold", va="center", zorder=7)
        ax.set_title("%s\ns = %d m" % (NAMES[r] if cnum == 0 else "", k), fontsize=8.5)
        ax.set_xlabel("offset (ft)", fontsize=8); ax.tick_params(labelsize=7)
        if cnum == 0:
            ax.set_ylabel("elevation (ft NAVD88)", fontsize=8.5)
        ax.grid(alpha=.25)
axes[0, 0].plot([], [], c="#5a3010", lw=1.7, label="ground (1 m LiDAR)")
axes[0, 0].plot([], [], c="#1a4f95", lw=1.4, ls="--", label="toe-to-toe chord = the profile's zero")
axes[0, 0].legend(fontsize=7.5, loc="lower left")
fig.suptitle("C10 existing berm — perpendicular cross-sections (chord over %.0f ft, "
             "crest height measured above it)" % (2 * A.toe_ft), fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.965])
if A.fig:
    fig.savefig(A.fig, bbox_inches="tight"); print("wrote " + A.fig)
