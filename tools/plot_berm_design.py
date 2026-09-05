#!/usr/bin/env python3
"""
Berm design profile and typical section.

    python3 tools/plot_berm_design.py

The section geometry and crest are read from tools/carve_berm.py, which is the
single definition used to build the modelled terrain.  This figure was once
hand-parameterised and drifted out of step with the carve -- an 8 ft crest on the
drawing against a 6 ft crest in the model -- so it now derives everything.
"""
import re, os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from osgeo import gdal, ogr
gdal.UseExceptions(); ogr.UseExceptions()

FT   = 3.280839895
SITE = "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097"
OUT  = SITE + "/Report_for_C6/Figures/NM15-097_berm_profile.png"

src = open(os.path.join(os.path.dirname(__file__), "carve_berm.py")).read()
def const(name, default):
    m = re.search(name + r"\s*=\s*(?:float\(sys\.argv\[\d\]\)\s*if.*?else\s*)?(-?[\d.]+)", src)
    return float(m.group(1)) if m else default
TOP_W  = const("^TOP_W", 6.0)
SIDE_Z = const("^SIDE_Z", 3.0)
GRADE  = const("^GRADE", -0.0020)
CREST0 = const("^CREST_AT_315", 6153.4)
STAREF = const("^STA_REF", 315.0)
print("read from carve_berm.py: crest %.1f ft at sta %.0f, grade %.3f%%, %.0f ft top, %.1fH:1V"
      % (CREST0, STAREF, 100 * GRADE, TOP_W, SIDE_Z))

S = "/tmp/claude-1000/-home-arash-Dropbox-Chloeta-NHA-WS3-GIS-NM15-097/0de18460-b665-429a-9baa-4a8fa3f3d0d0/scratchpad/"
t, x, y, zg, w, need = np.load(S + "berm_need.npy")
ds = np.gradient(t)
crest = np.where(t >= STAREF, CREST0 + GRADE * (t - STAREF), CREST0)
h = np.clip(crest - zg, 0, None)
flat = np.full_like(crest, np.nanmax(need))
hflat = np.clip(flat - zg, 0, None)
vol  = ((h * (TOP_W + SIDE_Z * h)) * ds).sum() / 27.0
volf = ((hflat * (TOP_W + SIDE_Z * hflat)) * ds).sum() / 27.0
hmean = h[h > 0.05].mean()
print("crest height mean %.1f ft, max %.1f ft; fill %,.0f cy" .replace(",", "") % (hmean, h.max(), vol))

fig = plt.figure(figsize=(14, 9.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1.9, 1.0], hspace=0.30)
ax = fig.add_subplot(gs[0])
ax.fill_between(t, zg.min() - 3, zg, color="#c8b89a", zorder=1)
ax.fill_between(t, zg, crest, where=crest > zg + 0.02, color="#e07000", alpha=0.5, zorder=2,
                label="Embankment")
ax.plot(t, zg, color="#5a4a32", lw=1.9, zorder=4, label="Existing ground")
ax.fill_between(t, zg, w, where=w > zg, color="#1f77b4", alpha=0.35, zorder=3)
ax.plot(t, w, color="#0050b3", lw=1.7, zorder=5, label="Flood water surface (envelope, both storms)")
ax.plot(t, np.where(w > zg, w + 2.0, np.nan), color="#9d4edd", lw=1.4, ls=":", zorder=5,
        label="Required crest = flood WSE $+$ 2 ft")
ax.plot(t, crest, color="#8a3800", lw=2.6, zorder=6,
        label="Adopted crest: flat %.1f ft, then %.3f%%" % (CREST0, 100 * GRADE))
ax.plot(t, flat, color="#c1121f", lw=1.8, ls="--", zorder=5,
        label="Level crest for comparison (%s cy)" % format(round(volf), ","))
for xv in (STAREF, 587):
    ax.axvline(xv, color="k", lw=1.0, ls=":", zorder=7)
ax.annotate("existing berm", (150, np.nanmax(need) + 1.4), ha="center", fontsize=9, style="italic")
ax.annotate("new embankment", (1850, np.nanmax(need) + 1.4), ha="center", fontsize=9.5, style="italic")
ax.set_xlim(0, t[-1]); ax.set_ylim(zg.min() - 3, np.nanmax(need) + 2.6)
ax.set_xlabel("Station along the berm centreline (ft) — south-west end at 0", fontsize=10)
ax.set_ylabel("Elevation (ft, NAVD88)", fontsize=10)
ax.grid(alpha=0.25, lw=0.5); ax.legend(loc="upper right", fontsize=8.6, framealpha=0.95)
ax.set_title("C6 boundary berm — design profile along the centreline (3,117 ft)",
             fontsize=12.5, fontweight="bold")

ax2 = fig.add_subplot(gs[1])
hh = hmean
xs = np.array([-(TOP_W / 2 + SIDE_Z * hh), -TOP_W / 2, TOP_W / 2, TOP_W / 2 + SIDE_Z * hh])
ys = np.array([0, hh, hh, 0])
ax2.plot([-34, 34], [0, 0], color="#5a4a32", lw=2.0, zorder=2)
ax2.fill_between(xs, 0, ys, color="#c8b89a", ec="#5a4a32", lw=2.2, zorder=3)
ax2.fill_between([-34, -(TOP_W / 2 + SIDE_Z * (hh - 2))], 0, hh - 2, color="#1f77b4", alpha=0.32, zorder=1)
ax2.plot([-34, -(TOP_W / 2 + SIDE_Z * (hh - 2))], [hh - 2, hh - 2], color="#0050b3", lw=2.0, zorder=4)
ax2.annotate("flood side", (-31, hh - 1.75), fontsize=10, color="#0050b3", fontweight="bold")
ax2.annotate("protected side\n(subdivision)", (17, 0.5), fontsize=10, fontweight="bold", color="#5a4a32")
ax2.annotate("", xy=(-TOP_W / 2 - 1.4, hh), xytext=(-TOP_W / 2 - 1.4, hh - 2),
             arrowprops=dict(arrowstyle="<->", color="k", lw=1.3))
ax2.annotate("2 ft\nfreeboard", (-TOP_W / 2 - 2.6, hh - 1.0), ha="right", va="center", fontsize=9.5)
ax2.annotate("", xy=(-TOP_W / 2, hh + 0.8), xytext=(TOP_W / 2, hh + 0.8),
             arrowprops=dict(arrowstyle="<->", color="k", lw=1.3))
ax2.annotate("%.0f ft crest width" % TOP_W, (0, hh + 1.15), ha="center", fontsize=10.5, fontweight="bold")
ax2.annotate("%.0f:1" % SIDE_Z, (TOP_W / 2 + SIDE_Z * hh * 0.5 + 1.8, hh * 0.45), fontsize=11, fontweight="bold")
ax2.annotate("%.0f:1" % SIDE_Z, (-(TOP_W / 2 + SIDE_Z * hh * 0.5) - 4, hh * 0.45), fontsize=11, fontweight="bold")
ax2.set_xlim(-34, 34); ax2.set_ylim(-1.8, hh + 2.8); ax2.set_aspect("equal")
ax2.set_xticks([]); ax2.set_yticks([])
for sp in ax2.spines.values(): sp.set_visible(False)
ax2.set_title("Typical section — %.0f ft crest, %.0fH:1V side slopes, at the mean height of %.1f ft (true scale)"
              % (TOP_W, SIDE_Z, hh), fontsize=11, fontweight="bold")
fig.text(0.5, 0.045,
         "Adopted crest %s cy against %s cy for a level crest. Height 0 to %.1f ft. "
         "Section and crest read from tools/carve_berm.py, the definition used to build the modelled terrain."
         % (format(round(vol), ","), format(round(volf), ","), h.max()),
         ha="center", fontsize=9.5, style="italic")
plt.savefig(OUT, bbox_inches="tight")
print("wrote", OUT)
