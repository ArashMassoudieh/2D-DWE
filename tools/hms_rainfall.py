#!/usr/bin/env python3
"""
Extract an HMS excess-precipitation time series as a CSV the solver applies as a
uniform source term (direct rainfall on the 2D domain).

    python3 tools/hms_rainfall.py [dss] [subbasin] [out.csv]

WHY EXCESS PRECIPITATION AND NOT GROSS RAINFALL
  The solver has no infiltration model. Applying the 2.03 in gross storm depth
  would route all of it as runoff and overstate the direct contribution by an
  order of magnitude. PRECIP-EXCESS is the depth remaining after the SCS
  curve-number loss HMS already applied, i.e. exactly the runoff-generating
  depth, computed with the same loss method and the same hyetograph as the rest
  of the study.

WHY Subbasin_70, AND WHY THIS IS NOT DOUBLE COUNTING
  The inflow hydrograph at the west boundary is the routed outflow of Reach_58,
  which carries everything upstream of it. Adding rain on the 2D domain would
  double count if the domain drained into that hydrograph. It does not:

    * Only 147 acres of the 2,275-acre domain fall inside ANY delineated
      subwatershed. The community bench lies OUTSIDE the watershed entirely,
      because the controlling divide IS the watershed boundary.
    * That 147 acres belongs to Subbasin_70, whose flow path is
      Subbasin_70 -> Outlet. It never passes through Reach_58.

  So none of the domain's own runoff is present in the boundary hydrograph, and
  all of it is genuinely additional.

  Subbasin_70 is used because it is the subbasin that physically overlaps the
  domain and shares its terrain and soils. Note it is near the LOW end of the
  basin: 0.205 in of excess against a basin mean of 0.85 in, reflecting locally
  high infiltration on the valley alluvium. The curve-number grid covers only
  5% of the domain (it too is clipped to the watershed), but where it does, the
  mean CN of 73 implies 0.335 in of excess -- 63% more than used here. This
  choice is therefore conservative-low and the sensitivity is documented.
"""
import sys
import numpy as np
from pydsstools.heclib.dss import HecDss

DSS  = sys.argv[1] if len(sys.argv) > 1 else \
       "/home/arash/Dropbox/Chloeta/NHA/WS3_HMS/NM15_097/Run_6hr_100yr.dss"
ELEM = sys.argv[2] if len(sys.argv) > 2 else "Subbasin_70"
OUT  = sys.argv[3] if len(sys.argv) > 3 else \
       "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097/outputs_2D/rain_excess_6hr_100yr.csv"

IN2M = 0.0254

d = HecDss.Open(DSS)
paths = [p for p in d.search_path('/*/*/PRECIP-EXCESS/*/*/*/')
         if p.split('/')[2] == ELEM]
if not paths:
    sys.exit("no PRECIP-EXCESS path for %s" % ELEM)

rows = []
interval = None
for p in paths:
    ts = d.read_ts(p)
    v  = np.array(ts.values, dtype=float)
    m  = np.array(ts.nodata, dtype=bool)
    t  = np.array([x.datetime() for x in ts.times])
    interval = ts.interval
    rows.append((t[~m], v[~m]))
d.close()

t = np.concatenate([r[0] for r in rows])
q = np.concatenate([r[1] for r in rows])      # incremental depth per interval, inches
o = np.argsort(t); t, q = t[o], q[o]
keep = np.r_[True, np.diff(t).astype('timedelta64[s]').astype(float) > 0]
t, q = t[keep], q[keep]

secs = np.array([(x - t[0]).total_seconds() for x in t])

# Incremental DEPTH over each interval -> RATE in m/s. The value at index i is
# the depth accumulated during the interval ENDING at secs[i], so it is applied
# as a constant rate over the interval that precedes it: the series is written
# as a step function sampled at the interval start.
rate = q * IN2M / interval                     # m/s
t_start = np.r_[0.0, secs[:-1]]

with open(OUT, "w") as f:
    f.write("# HMS %s excess precipitation, %s\n" % (ELEM, DSS.split('/')[-1]))
    f.write("# total %.4f in (%.4f m) over %.2f h\n"
            % (q.sum(), q.sum() * IN2M, secs[-1] / 3600.))
    f.write("# peak rate %.4f in/h (%.3e m/s)\n"
            % (rate.max() * 3600 / IN2M, rate.max()))
    f.write("# applied as a spatially uniform source; see the module docstring\n")
    f.write("t_s,rate_mps\n")
    for a, b in zip(t_start, rate):
        f.write("%.1f,%.9e\n" % (a, b))
    f.write("%.1f,0.000000000\n" % secs[-1])   # rain stops at the end of the record

DOMAIN_AC = 2275.0
vol_acft = q.sum() / 12.0 * DOMAIN_AC
print("element    %s  (%s)" % (ELEM, DSS.split('/')[-1]))
print("records    %d at %d s" % (len(t), interval))
print("total      %.4f in excess over %.2f h" % (q.sum(), secs[-1] / 3600.))
print("peak rate  %.3f in/h" % (rate.max() * 3600 / IN2M))
print("volume     %.1f ac-ft over the %.0f-acre domain" % (vol_acft, DOMAIN_AC))
print("           = %.2f%% of the 4,720 ac-ft boundary hydrograph"
      % (100 * vol_acft / 4720.0))
print("wrote      %s" % OUT)
