#!/usr/bin/env python3
"""
Animate a run's .vti frames to mp4 (or avi) with a colour bar, map axes, a north
arrow, a scale bar and a clock.

    python3 tools/make_movie.py out_c6_3m 3                 # depth, mp4
    python3 tools/make_movie.py out_c6_3m 3 --field speed
    python3 tools/make_movie.py out_c6_3m 3 --vectors --fps 4 --format avi

NOTE ON THE FILE TYPE
  The solver writes .vti (VTK ImageData), not .vtu. .vtu is the UNSTRUCTURED
  grid container -- it stores an explicit node list and cell connectivity, which
  a uniform structured grid does not need and which would multiply the file size
  for no gain. This reads .vti; nothing else is produced.

WHY MATPLOTLIB AND NOT PARAVIEW BATCH
  pvpython would render this too, but it makes the animation depend on a working
  ParaView install and on a state file that has to be kept in step with the
  solver's array names. Reading the XML directly keeps the pipeline to numpy +
  matplotlib + ffmpeg, all of which are already required elsewhere in this repo.

FRAME TIMING
  Frame times come from the run log's output cadence (--tout), passed here as
  --dt-frame, defaulting to 1800 s. The .vti container carries no time value, so
  a wrong cadence produces a plausible-looking clock that is simply wrong; the
  value used is printed and stamped in the corner.
"""
import argparse, os, re, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, LightSource
from osgeo import gdal

gdal.UseExceptions()
FT = 3.280839895

ap = argparse.ArgumentParser()
ap.add_argument("run")
ap.add_argument("dx", type=float)
ap.add_argument("--field", default="depth", choices=["depth", "speed", "head", "total_head"])
ap.add_argument("--fps", type=float, default=3.0)
ap.add_argument("--dt-frame", type=float, default=1800.0, help="seconds between frames")
ap.add_argument("--format", default="mp4", choices=["mp4", "avi"])
ap.add_argument("--vectors", action="store_true", help="overlay velocity arrows")
ap.add_argument("--dpi", type=int, default=130)
ap.add_argument("--out", default=None)
A = ap.parse_args()

DEM = "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097/outputs_2D/cliped_dem.tif"

frames = sorted(f for f in os.listdir(A.run) if re.match(r".*_\d{4}\.vti$", f))
if not frames:
    sys.exit("no .vti frames in " + A.run)


def read_vti(path, names):
    """Pull named arrays out of a .vti. Returns (dict, (nx, ny), origin, spacing)."""
    with open(path) as fh:
        s = fh.read()
    m = re.search(r'WholeExtent="0 (\d+) 0 (\d+)', s)
    nx, ny = int(m.group(1)), int(m.group(2))
    o = re.search(r'Origin="([-\d.eE+]+) ([-\d.eE+]+)', s)
    sp = re.search(r'Spacing="([-\d.eE+]+)', s)
    out = {}
    for n in names:
        k = s.split('Name="%s"' % n)
        if len(k) < 2:
            continue
        blk = k[1].split(">", 1)[1].split("</DataArray")[0]
        v = np.fromstring(blk, sep=" ")
        ncomp = 3 if 'NumberOfComponents="3"' in k[1].split(">", 1)[0] else 1
        if ncomp == 3:
            v = v.reshape(ny, nx, 3)[::-1]          # emitted bottom-up
        else:
            v = v.reshape(ny, nx)[::-1]
        out[n] = v
    return out, (nx, ny), (float(o.group(1)), float(o.group(2))), float(sp.group(1))

want = [A.field] + (["velocity"] if A.vectors else [])
d0, (nx, ny), (ox, oy), spac = read_vti(os.path.join(A.run, frames[0]), want)
ext = [ox, ox + nx * spac, oy, oy + ny * spac]

t = gdal.Warp("", DEM, format="MEM", width=nx, height=ny, resampleAlg="average")
zb = t.GetRasterBand(1).ReadAsArray().astype(float)
hs = LightSource(azdeg=315, altdeg=45).hillshade(zb, vert_exag=3, dx=A.dx, dy=A.dx)

# A FIXED colour scale across the whole animation. Rescaling per frame makes a
# receding flood look constant -- the colours stay saturated while the water
# leaves -- which is exactly backwards from what the animation is meant to show.
peak = 0.0
for f in frames:
    dd, *_ = read_vti(os.path.join(A.run, f), [A.field])
    if A.field in dd:
        peak = max(peak, float(np.nanmax(dd[A.field])))

if A.field == "depth":
    lv = [0.5, 1, 2, 3, 5, 8, 12, max(25.0, peak * FT * 1.02)]
    cm = ListedColormap(["#c6e5f5", "#8ec9ec", "#4ba3d8", "#2272b5",
                         "#1a4f95", "#0d2f6e", "#5a1a6e"])
    label, conv, floor = "Water depth (ft)", FT, 0.5
elif A.field == "speed":
    lv = [0.5, 1, 2, 4, 6, 9, 13, max(20.0, peak * FT * 1.02)]
    cm = ListedColormap(["#d8f0d8", "#a8dda8", "#f4e04d", "#f2a541",
                         "#e8712d", "#d1373a", "#7a1020"])
    label, conv, floor = "Flow speed (ft/s)", FT, 0.5
else:
    lv = list(np.linspace(zb.min() * FT, zb.max() * FT, 8))
    cm = ListedColormap(plt.get_cmap("viridis")(np.linspace(0, 1, 7)))
    label, conv, floor = label_map = ("%s (ft)" % A.field, FT, -1e30)

norm = BoundaryNorm(lv, cm.N)
out = A.out or os.path.join(A.run, "c6_%s.%s" % (A.field, A.format))
codec = ["-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20"] if A.format == "mp4" \
        else ["-c:v", "mpeg4", "-q:v", "4"]

fig, ax = plt.subplots(figsize=(15, 7.6), dpi=A.dpi)
proc = None
print("run %s | %d frames | field %s | %.0f s per frame | peak %.3f (SI)"
      % (A.run, len(frames), A.field, A.dt_frame, peak))

for i, f in enumerate(frames):
    dd, *_ = read_vti(os.path.join(A.run, f), want)
    v = dd[A.field] * conv
    ax.clear()
    ax.imshow(hs, cmap="gray", extent=ext, alpha=.9, zorder=1, origin="upper")
    im = ax.imshow(np.ma.masked_less(v, floor), cmap=cm, norm=norm, extent=ext,
                   alpha=.85, zorder=2, origin="upper")
    if A.vectors and "velocity" in dd:
        st = max(1, nx // 60)
        yy, xx = np.mgrid[0:ny:st, 0:nx:st]
        vx = dd["velocity"][::st, ::st, 0]
        vy = dd["velocity"][::st, ::st, 1]
        X = ext[0] + (xx + .5) * spac
        Y = ext[3] - (yy + .5) * spac
        msk = np.hypot(vx, vy) > 0.15
        ax.quiver(X[msk], Y[msk], vx[msk], vy[msk], color="k", alpha=.55,
                  scale=90, width=.0016, zorder=5)

    if i == 0:
        cb = fig.colorbar(im, ax=ax, shrink=.88, pad=.012, ticks=lv[:-1])
        cb.set_label(label, fontsize=10)

    tt = i * A.dt_frame
    ax.set_title("C6 Coyote Canyon — 2D-DWE, %.0f m grid, 100-yr 6-hr transient "
                 "hydrograph\nt = %2d h %02d min   (frame %d of %d)"
                 % (A.dx, int(tt // 3600), int((tt % 3600) // 60), i + 1, len(frames)),
                 fontsize=12)
    ax.set_xlabel("UTM Zone 12N Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect("equal")

    bl = ext[0] + 300
    ax.plot([bl, bl + 1609.34], [ext[2] + 140, ext[2] + 140], c="k", lw=4,
            zorder=9, solid_capstyle="butt")
    ax.text(bl + 804, ext[2] + 195, "1 mile", ha="center", fontsize=9, zorder=9)
    ax.annotate("", xy=(ext[1] - 240, ext[3] - 180), xytext=(ext[1] - 240, ext[3] - 480),
                arrowprops=dict(arrowstyle="-|>", lw=2, color="k"), zorder=9)
    ax.text(ext[1] - 240, ext[3] - 140, "N", ha="center", fontsize=11,
            weight="bold", zorder=9)
    fig.tight_layout()

    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba())[:, :, :3]
    if proc is None:
        h, w = buf.shape[:2]
        proc = subprocess.Popen(
            ["ffmpeg", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24",
             "-s", "%dx%d" % (w, h), "-r", str(A.fps), "-i", "-",
             *codec, "-r", str(max(A.fps, 10)), out],
            stdin=subprocess.PIPE, stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL)
    proc.stdin.write(buf.tobytes())
    print("  frame %2d/%d  t=%5.2f h  max %.2f" % (i + 1, len(frames), tt / 3600., v.max()))

proc.stdin.close(); proc.wait()
print("wrote %s  (%.1f MB)" % (out, os.path.getsize(out) / 1e6))
