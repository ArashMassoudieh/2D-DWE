#!/usr/bin/env python3
"""
Drive the C10 (NM15-141) Tc - Q channel-depth iteration to convergence.

    python3 tools/converge_tc_c10.py [--max-passes N] [--tol-min M]

WHY
  Manual Sec. 2.4.4.2 resolves a circularity: the channel segment of Tc needs a
  flow depth, the depth follows from the subbasin peak, and the peak follows
  from the lag derived from Tc. Pass 1 seeds the channel depth at 1.5 ft; each
  later pass replaces the seed with the normal depth implied by that subbasin's
  modelled peak, and HEC-HMS is re-run, until the parameters stop moving.

  C10 had only ever been run at pass 1. The report's own notebox says so and
  predicts the direction: a discharge-derived depth exceeds the seed, raising
  channel velocity, so Tc shortens and the peaks rise. This closes it and
  records the trajectory rather than asserting it.

ONE PASS
  compute_tc.py   reads WS3_HMS/<model>/Run_24hr_100yr.dss (the Manual's design
                  storm), solves the triangular normal depth from each subbasin
                  peak, rewrites Tc and lag into subwatershed_params.gpkg
  write_basin.py  regenerates <model>.basin straight into the HMS project
  HEC-HMS         recomputes both design storms headless (Xvfb), refreshing the
                  DSS that the next pass reads

CONVERGENCE is on Tc per subbasin, not on the mean: a mean can sit still while
individual subbasins oscillate. The loop stops when the largest absolute change
across all 49 falls below --tol-min.
"""
import argparse, os, shutil, subprocess, sys, time
import numpy as np
from osgeo import ogr
ogr.UseExceptions()

NHA   = "/home/arash/Dropbox/Chloeta/NHA"
SITE  = "WS3_GIS/NM15-141Local"
PY_D  = NHA + "/PythonScripts"
OUTD  = NHA + "/" + SITE + "/outputs"
HMSD  = NHA + "/WS3_HMS/NM15_141Local"
HMS   = "/home/arash/Projects/HEC-HMS-4.13"
PARAMS = OUTD + "/subwatershed_params.gpkg"

ap = argparse.ArgumentParser()
ap.add_argument("--max-passes", type=int, default=8)
ap.add_argument("--tol-min", type=float, default=0.10)
ap.add_argument("--start-pass", type=int, default=2)
A = ap.parse_args()

SCRIPT = "/tmp/claude-1000/-home-arash-Projects-2D-DWE/615dd2c4-6a45-4d7c-8d78-98fee05c9920/scratchpad/hms_compute.script"


def tc_vector():
    ds = ogr.Open(PARAMS); ly = ds.GetLayer()
    d = {}
    for f in ly:
        d[f.GetField("id")] = (f.GetField("tc_min"), f.GetField("lag_min"))
    return d


def peaks(dss):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from hms_peaks import peaks as _p
    return _p(dss)


def run(cmd, env=None, cwd=None, label=""):
    e = dict(os.environ); e.update(env or {})
    t0 = time.time()
    r = subprocess.run(cmd, shell=True, cwd=cwd, env=e,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if r.returncode != 0 and label != "hms":     # HMS exits non-zero via Exit(1)
        print(r.stdout.decode(errors="ignore")[-2500:])
        raise SystemExit("FAILED (%d): %s" % (r.returncode, label))
    print("      %-14s %5.1fs" % (label, time.time() - t0))
    return r.stdout.decode(errors="ignore")


hist = []
prev = tc_vector()
o6 = peaks(HMSD + "/Run_6hr_100yr.dss").get("Outlet")
o24 = peaks(HMSD + "/Run_24hr_100yr.dss").get("Outlet")
tc0 = np.array([v[0] for v in prev.values()], float)
print("pass %d (current state): Tc mean %.2f max %.2f | outlet 6hr %.1f  24hr %.1f"
      % (A.start_pass, tc0.mean(), tc0.max(), o6, o24))
hist.append((A.start_pass, tc0.mean(), tc0.max(), o6, o24, float("nan")))

for p in range(A.start_pass + 1, A.start_pass + A.max_passes):
    print("\n--- pass %d ---" % p)
    run("python3 write_basin.py %s" % SITE,
        {"NHA_ROOT": NHA + "/", "QT_QPA_PLATFORM": "offscreen"}, PY_D, "write_basin")
    run("xvfb-run -a ./hec-hms.sh -s %s" % SCRIPT, None, HMS, "hms")
    run("python3 compute_tc.py %s" % SITE,
        {"NHA_ROOT": NHA + "/", "QT_QPA_PLATFORM": "offscreen"}, PY_D, "compute_tc")

    cur = tc_vector()
    ids = sorted(set(prev) & set(cur))
    d = np.array([abs(cur[i][0] - prev[i][0]) for i in ids], float)
    tc = np.array([cur[i][0] for i in ids], float)
    o6 = peaks(HMSD + "/Run_6hr_100yr.dss").get("Outlet")
    o24 = peaks(HMSD + "/Run_24hr_100yr.dss").get("Outlet")
    print("   Tc mean %.2f  max %.2f | max |dTc| %.3f min | outlet 6hr %.1f  24hr %.1f"
          % (tc.mean(), tc.max(), d.max(), o6, o24))
    hist.append((p, tc.mean(), tc.max(), o6, o24, d.max()))
    prev = cur
    if d.max() < A.tol_min:
        print("\nCONVERGED: max |dTc| %.3f < %.2f min" % (d.max(), A.tol_min))
        break
else:
    print("\nstopped at the pass limit without meeting the tolerance")

print("\n%-6s %10s %10s %12s %12s %12s" % ("pass", "Tc mean", "Tc max", "outlet 6hr", "outlet 24hr", "max |dTc|"))
for r in hist:
    print("%-6d %10.2f %10.2f %12.1f %12.1f %12s"
          % (r[0], r[1], r[2], r[3], r[4], "-" if r[5] != r[5] else "%.3f" % r[5]))
