#!/usr/bin/env python3
"""
Extract an HMS element hydrograph from a run DSS as a CSV the solver can read.

WHY Reach_58 AND NOT Junction_28
  Junction_28 carries the 32,628 cfs peak quoted in the report, but it sits at
  (702,024, 3,960,135) -- 1,190 m WEST of the 2D domain, outside it. Reach_58
  routes that flow downstream and discharges at (703,075, 3,959,909), 115 m from
  the west-edge inflow face. Its outflow is what physically arrives at the
  boundary, attenuated to 32,431 cfs -- within 0.6 % of the 1D figure, so the two
  models stay consistent while the 2D one is fed at its own boundary.

WHY THE THREE DATE BLOCKS ARE CONCATENATED, NOT MAX'D
  HMS splits a run across calendar days: 31Dec1999 / 01Jan2000 / 02Jan2000 here.
  These are consecutive pieces of ONE event, not competing versions of it. Taking
  a max across them -- the failure that contaminated the Tc iteration -- would
  collapse a 36-hour hydrograph to a single number. Each block also carries a
  nodata mask; the padding is masked, not absent, and reading `values` without
  applying `ts.nodata` returns the DSS undefined sentinel (2^65) as live data.
"""
import sys, numpy as np
from pydsstools.heclib.dss import HecDss

DSS  = sys.argv[1] if len(sys.argv) > 1 else \
       "/home/arash/Dropbox/Chloeta/NHA/WS3_HMS/NM15_097/Run_6hr_100yr.dss"
ELEM = sys.argv[2] if len(sys.argv) > 2 else "Reach_58"
OUT  = sys.argv[3] if len(sys.argv) > 3 else \
       "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097/outputs_2D/inflow_6hr_100yr.csv"
CFS2CMS = 0.028316846592

d = HecDss.Open(DSS)
paths = [p for p in d.search_path('/*/*/FLOW/*/*/*/') if p.split('/')[2] == ELEM]
if not paths:
    sys.exit("no FLOW path for %s in %s" % (ELEM, DSS))

rows = []
for p in paths:
    ts = d.read_ts(p)
    v  = np.array(ts.values, dtype=float)
    m  = np.array(ts.nodata, dtype=bool)
    t  = np.array([x.datetime() for x in ts.times])
    rows.append((t[~m], v[~m]))
d.close()

t = np.concatenate([r[0] for r in rows])
q = np.concatenate([r[1] for r in rows])
o = np.argsort(t); t, q = t[o], q[o]
keep = np.r_[True, np.diff(t).astype('timedelta64[s]').astype(float) > 0]
t, q = t[keep], q[keep]

secs = np.array([(x - t[0]).total_seconds() for x in t])
qm   = q * CFS2CMS

with open(OUT, "w") as f:
    f.write("# HMS %s outflow, %s\n" % (ELEM, DSS.split('/')[-1]))
    f.write("# peak %.1f cfs (%.2f m3/s) at t = %.0f s (%.2f h)\n"
            % (q.max(), qm.max(), secs[q.argmax()], secs[q.argmax()] / 3600.))
    f.write("# volume %.0f ac-ft\n" % (np.trapz(q, secs) / 43560.0))
    f.write("t_s,Q_cms\n")
    for a, b in zip(secs, qm):
        f.write("%.1f,%.6f\n" % (a, b))

pk = q.argmax()
rec = np.where(q[pk:] < 0.05 * q.max())[0]
print("element   %s   (%s)" % (ELEM, DSS.split('/')[-1]))
print("records   %d at %.0f s (%.1f min) steps, span %.2f h"
      % (len(t), secs[1] - secs[0], (secs[1] - secs[0]) / 60., secs[-1] / 3600.))
print("peak      %.1f cfs = %.2f m3/s at %.2f h" % (q.max(), qm.max(), secs[pk] / 3600.))
print("volume    %.0f ac-ft" % (np.trapz(q, secs) / 43560.0))
print("recedes   to 5%% of peak at %.2f h" %
      (secs[pk + rec[0]] / 3600. if rec.size else float('nan')))
print("wrote     %s" % OUT)
