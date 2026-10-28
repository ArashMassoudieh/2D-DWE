#!/usr/bin/env python3
"""Peak FLOW per location from an HEC-HMS run DSS.

    python3 tools/hms_peaks.py <run.dss> [<run.dss> ...]

Reads EVERY dated block and takes the max, because reading with a blanked
D-part truncates the record and under-reports the peak (the same trap
compute_tc.py documents).
"""
import sys, os

def peaks(dss_path):
    from pydsstools.heclib.dss import HecDss
    out = {}
    with HecDss.Open(dss_path) as fh:
        try:
            allp = fh.getPathnameList("", sort=1)
        except AttributeError:
            allp = fh.path_dict().get("ts", [])
        for p in allp:
            parts = p.split("/")
            if len(parts) <= 4 or parts[3].upper() != "FLOW":
                continue
            loc = parts[2]
            try:
                ts = fh.read_ts(p)
            except Exception:
                continue
            for v in ts.values:
                try:
                    v = float(v)
                except (TypeError, ValueError):
                    continue
                if v == v and v > -1e30 and v > out.get(loc, 0.0):
                    out[loc] = v
    return out

if __name__ == "__main__":
    for d in sys.argv[1:]:
        p = peaks(d)
        subs = {k: v for k, v in p.items() if k.lower().startswith("subbasin")}
        print("%-22s outlet %10.2f cfs | %d subbasins, sum of peaks %9.1f"
              % (os.path.basename(d), p.get("Outlet", float("nan")), len(subs), sum(subs.values())))
