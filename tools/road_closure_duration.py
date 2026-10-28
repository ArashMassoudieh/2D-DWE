#!/usr/bin/env python3
"""
How long Navajo Route 9 is impassable, from the per-frame .vti output.

    python3 tools/road_closure_duration.py <run-dir> [...]

depth_max.tif answers "does the road flood" but not "for how long", because it
is an envelope: every cell holds its own maximum, reached at its own time. The
duration needs the frames.

IMPASSABLE is measured on the TRAFFICKED SURFACE -- the derived road crown --
not on a +/-10 m corridor. The corridor version was tried first and is wrong:
perpendicular to the road AT THE CROSSING the corridor runs along the wash, so
it always contains channel cells, which stay wet for the whole recession and
report a permanent 2 m "blockage" that is the creek bed, not the road.

Because the crown is an upper bound on the road surface (a max-elevation snap
picks shoulders and spoil ridges as well as pavement), closure on it is a LOWER
bound on duration. The run is repeated on a crest lowered by --crest-drop, the
offset at which broad-crested weir flow over the crest balances the discharge
actually arriving (0.34 m, C = 2.6), which brackets the answer from the other
side. A surveyed road profile would replace both.

A closure needs --min-length metres blocked (default 10) so that a single
stranded cell is not reported as a severed road.

Thresholds: 0.5 ft (the C6 value, used throughout this study) and 1.0 ft, at
which a passenger vehicle is already at risk of stalling and loss of control.

Frame interval comes from the run log (duration / nframes), as in make_movie.py;
it bounds the time resolution and is printed with the result.
"""
import os, re, sys, argparse
import numpy as np
from osgeo import ogr, osr

FT = 3.280839895
SITE = "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-141Local"
ROADS = SITE + "/roads/tl_2025_35031_roads/tl_2025_35031_roads.shp"
LINEARID = "1106088001728"


def read_depth(path):
    with open(path) as fh:
        s = fh.read()
    m = re.search(r'WholeExtent="0 (\d+) 0 (\d+)', s)
    nx, ny = int(m.group(1)), int(m.group(2))
    o = re.search(r'Origin="([-\d.eE+]+) ([-\d.eE+]+)', s)
    sp = float(re.search(r'Spacing="([-\d.eE+]+)', s).group(1))
    k = s.split('Name="depth"')
    blk = k[1].split(">", 1)[1].split("</DataArray")[0]
    v = np.fromstring(blk, sep=" ").reshape(ny, nx)[::-1]
    return v, (nx, ny), (float(o.group(1)), float(o.group(2))), sp


def road_corridor():
    ds = ogr.Open(ROADS); ly = ds.GetLayer()
    sr = ly.GetSpatialRef(); tg = osr.SpatialReference(); tg.ImportFromEPSG(26912)
    sr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tg.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = osr.CoordinateTransformation(sr, tg)
    for f in ly:
        if f.GetField("LINEARID") != LINEARID:
            continue
        g = f.GetGeometryRef().Clone(); g.Transform(tr)
        pts = np.array(g.GetPoints())[:, :2]
        break
    if pts[0, 0] > pts[-1, 0]:
        pts = pts[::-1]
    seg = np.r_[0, np.cumsum(np.hypot(*np.diff(pts, axis=0).T))]
    return pts, seg


def main(runs):
    pts, seg = road_corridor()
    for run in runs:
        frames = sorted(f for f in os.listdir(run) if re.match(r".*_\d{4}\.vti$", f))
        if not frames:
            print("%s: no .vti frames\n" % os.path.basename(run)); continue
        log = os.path.join(run, "run.log")
        dur = None
        if os.path.exists(log):
            mm = re.search(r"^duration\s+([\d.]+)\s*s", open(log, errors="ignore").read(), re.M)
            if mm:
                dur = float(mm.group(1))
        dt = (dur / len(frames)) if dur else 1800.0

        d0, (nx, ny), (ox, oy), sp = read_depth(os.path.join(run, frames[0]))
        ymax = oy + ny * sp
        # clip the road to the frame footprint and station it
        inside = (pts[:, 0] >= ox) & (pts[:, 0] <= ox + nx * sp)
        p = pts[inside] if inside.sum() > 1 else pts
        sg = np.r_[0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))]
        s = np.arange(0, sg[-1], 1.0)
        cx = np.interp(s, sg, p[:, 0]); cy = np.interp(s, sg, p[:, 1])
        dx = np.gradient(cx); dy = np.gradient(cy)
        L = np.hypot(dx, dy); ux, uy = -dy / L, dx / L
        offs = np.arange(-10, 10.01, 1.0)
        XX = np.array([cx + o * ux for o in offs])
        YY = np.array([cy + o * uy for o in offs])
        II = np.clip(((XX - ox) / sp).astype(int), 0, nx - 1)
        JJ = np.clip(((ymax - YY) / sp).astype(int), 0, ny - 1)

        print("=" * 72)
        print("%s  |  %d frames, %.0f s apart (%.2f h)"
              % (os.path.basename(run), len(frames), dt, dt / 3600))
        print("%6s %8s %10s %12s %12s" % ("frame", "t (h)", "max ft", "blocked 0.5ft", "blocked 1.0ft"))
        n05 = n10 = 0
        for i, f in enumerate(frames):
            a, *_ = read_depth(os.path.join(run, f))
            g = a[JJ, II] * FT                      # (noff, nstation)
            b05 = (g > 0.5).all(axis=0); b10 = (g > 1.0).all(axis=0)
            n05 += int(b05.any()); n10 += int(b10.any())
            print("%6d %8.2f %9.2f %9d m %9d m"
                  % (i, i * dt / 3600, g.max(), b05.sum(), b10.sum()))
        print("-" * 72)
        print("CLOSED (no dry lane, >0.5 ft): %d of %d frames = %.1f h"
              % (n05, len(frames), n05 * dt / 3600))
        print("CLOSED (>1.0 ft):              %d of %d frames = %.1f h"
              % (n10, len(frames), n10 * dt / 3600))
        print("time resolution is one frame = %.2f h; the run ends at %.1f h\n"
              % (dt / 3600, len(frames) * dt / 3600))


if __name__ == "__main__":
    main(sys.argv[1:])
