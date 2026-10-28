#!/usr/bin/env python3
"""Tc / lag summary from subwatershed_params.gpkg, for tracking the Tc-Q iteration."""
import sys
import numpy as np
from osgeo import ogr
ogr.UseExceptions()

def summary(path):
    ds = ogr.Open(path); ly = ds.GetLayer()
    tc, lag, src = [], [], {}
    for f in ly:
        tc.append(f.GetField("tc_min")); lag.append(f.GetField("lag_min"))
        s = f.GetField("q_src") if f.GetFieldIndex("q_src") >= 0 else "?"
        src[s] = src.get(s, 0) + 1
    tc = np.array([t for t in tc if t is not None], float)
    lag = np.array([l for l in lag if l is not None], float)
    return tc, lag, src

if __name__ == "__main__":
    for p in sys.argv[1:]:
        tc, lag, src = summary(p)
        print("%-58s" % p.split("/")[-1])
        print("   Tc  mean %7.2f  min %6.2f  max %7.2f" % (tc.mean(), tc.min(), tc.max()))
        print("   lag mean %7.2f  min %6.2f  max %7.2f" % (lag.mean(), lag.min(), lag.max()))
        print("   q_src: %s" % ", ".join("%s=%d" % kv for kv in sorted(src.items())))
