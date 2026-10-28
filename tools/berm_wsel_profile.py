#!/usr/bin/env python3
"""
Longitudinal profile of the existing berm against the modelled peak water
surface, with the crest raise each line needs.

    python3 tools/berm_wsel_profile.py --fig out.png

Crest is the highest ground within +/-5 ft of the digitised line, on the 1 m
LiDAR. Water surfaces are the 2-D peak envelopes from the converged runs. The
required crest is the governing (blocked-culvert) water surface plus freeboard;
where the existing crest already clears it, nothing is required.
"""
import argparse, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from osgeo import gdal, ogr
gdal.UseExceptions(); ogr.UseExceptions()

FT = 3.280839895
S = "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-141Local"
P = "/home/arash/Projects/2D-DWE"

ap = argparse.ArgumentParser()
ap.add_argument("--freeboard-ft", type=float, default=1.0)
ap.add_argument("--fig", default=None)
ap.add_argument("--csv", default=None)
A = ap.parse_args()

dem = gdal.Open(S + "/demhr/cliped_sd.tif"); G = dem.GetGeoTransform()
Z = dem.GetRasterBand(1).ReadAsArray().astype(float); NY, NX = Z.shape
def sdem(x, y):
    i = np.clip(((x - G[0]) / G[1]).astype(int), 0, NX - 1)
    j = np.clip(((G[3] - y) / abs(G[5])).astype(int), 0, NY - 1)
    return Z[j, i]

RUNS = [("6-hr, culvert open",    "out_nm141_1m_6hr_carved_2m_conv",  "#2e8b57", "--"),
        ("6-hr, culvert blocked", "out_nm141_1m_6hr_raw_conv",        "#7a1020", "-"),
        ("24-hr, culvert blocked","out_nm141_1m_24hr_raw_conv",       "#8a3800", "-.")]
W = {}
for lbl, run, _, _ in RUNS:
    d = gdal.Open(os.path.join(P, run, "wse_max.tif")); b = d.GetRasterBand(1)
    a = b.ReadAsArray().astype(float); nd = b.GetNoDataValue()
    W[lbl] = (np.where((a == nd) | ~np.isfinite(a), np.nan, a), d.GetGeoTransform(), a.shape)
def swse(lbl, x, y):
    a, gt, (ny, nx) = W[lbl]
    i = np.clip(((x - gt[0]) / gt[1]).astype(int), 0, nx - 1)
    j = np.clip(((gt[3] - y) / abs(gt[5])).astype(int), 0, ny - 1)
    return a[j, i]

ds = ogr.Open(S + "/outputs/existing_berm.shp"); ly = ds.GetLayer()
lines = [np.array(f.GetGeometryRef().GetPoints())[:, :2] for f in ly]
NAMES = ["Line 1 — west boundary", "Line 2 — north boundary"]

fig, axes = plt.subplots(2, 1, figsize=(14, 9.2), dpi=150)
rows = []
for k, p in enumerate(lines):
    seg = np.r_[0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))]
    s = np.arange(0, seg[-1], 1.0)
    px = np.interp(s, seg, p[:, 0]); py = np.interp(s, seg, p[:, 1])
    dx = np.gradient(px); dy = np.gradient(py)
    L = np.hypot(dx, dy); ux, uy = -dy / L, dx / L
    offs = np.arange(-5, 5.1, 1.0)
    crest = np.nanmax(np.array([sdem(px + o * ux, py + o * uy) for o in offs]), axis=0) * FT
    sf = s * FT
    ax = axes[k]
    ax.fill_between(sf, crest, crest.min() - 6, color="#c8b89a", zorder=2)
    ax.plot(sf, crest, c="#5a3010", lw=2.0, zorder=5, label="Existing berm crest (1 m LiDAR)")
    gov = None
    # Stations that were NEVER wet carry nodata in wse_max, so the WSE lines
    # break there. Those breaks are a result, not missing data -- on Line 1 the
    # crest simply stays dry over 29% of its length -- so the dry reaches are
    # shaded explicitly rather than left as ambiguous gaps in the line.
    for lbl, run, col, ls in RUNS:
        w = swse(lbl, px, py) * FT
        ax.plot(sf, w, c=col, lw=1.6, ls=ls, zorder=4, label="Peak WSE — " + lbl)
        if "6-hr, culvert blocked" in lbl: gov = w
    drygov = ~np.isfinite(gov)
    if drygov.any():
        ax.fill_between(sf, crest.min() - 6, crest, where=drygov, color="#999999", alpha=.20,
                        zorder=1, label="Never inundated (no water surface)")
    over = np.isfinite(gov) & (gov > crest)
    ax.fill_between(sf, crest, np.where(np.isfinite(gov), gov, crest), where=over,
                    color="#e8112d", alpha=.42, zorder=3, label="Overtopped (governing case)")
    # Required crest FOLLOWS the water surface; a single horizontal design
    # elevation is meaningless on a line whose WSE falls 6 ft end to end.
    # Comparing max(WSE) with min(crest) -- at different stations -- gives a
    # 7 ft "raise" on Line 1 where the real per-station deficit is under 1.5 ft.
    gv = np.where(np.isfinite(gov), gov, crest)
    req_p = gv + A.freeboard_ft
    deficit = np.maximum(req_p - crest, 0.0)
    ax.plot(sf, req_p, c="#111", lw=1.4, ls=":", zorder=6,
            label="Required crest (peak WSE + %.0f ft freeboard)" % A.freeboard_ft)
    raise_ = deficit.max()
    fill_cy = np.trapz(deficit, sf) * 8.0 / 27.0     # 8 ft mean crest width, cy
    rows.append((NAMES[k], sf[-1], crest.min(), crest.max(), np.nanmax(gv),
                 100 * over.mean(), raise_, deficit[deficit > 0.05].mean() if (deficit > 0.05).any() else 0.0,
                 fill_cy))
    ax.set_xlim(0, sf[-1]); ax.set_ylim(crest.min() - 4, max(req_p.max(), crest.max()) + 1.5)
    ax.set_ylabel("Elevation (ft NAVD88)")
    ax.set_xlabel("Distance along %s (ft)" % NAMES[k].split("—")[1].strip())
    ax.grid(alpha=.3, zorder=1)
    ax.set_title("%s — overtopped along %.0f%% of its length; maximum crest deficit %.2f ft"
                 % (NAMES[k], 100 * over.mean(), raise_), fontsize=11)
    ax.legend(fontsize=7.5, loc="lower left", ncol=3, framealpha=.93)
    ax.text(0.995, 0.03, "dry along %.0f%% of this line (governing case)" % (100 * drygov.mean()),
            transform=ax.transAxes, ha="right", fontsize=8, style="italic", color="#555")
fig.suptitle("C10 — existing perimeter berm against the modelled 100-year water surface "
             "(converged hydrology)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.96])
if A.fig:
    fig.savefig(A.fig, bbox_inches="tight"); print("wrote " + A.fig)
print("\n%-26s %8s %9s %9s %8s %9s %9s %10s"
      % ("line", "len ft", "crest min", "peak WSE", "overtop", "max raise", "mean raise", "fill cy"))
for r in rows:
    print("%-26s %8.0f %9.2f %9.2f %7.0f%% %9.2f %10.2f %10.0f"
          % (r[0], r[1], r[2], r[4], r[5], r[6], r[7], r[8]))
