#!/usr/bin/env bash
# Block until all four corrected-depth C6 runs exit, then report.
set -u
cd /home/arash/Projects/2D-DWE
RUNS="out_c6_1m_rain out_c6_1m_berm_6hr out_c6_1m_24hr_v2 out_c6_1m_berm_24hr"
while pgrep -f "2D-DWE c6" > /dev/null 2>&1; do sleep 120; done
echo "=========== ALL C6 SOLVERS EXITED ==========="
for R in $RUNS; do
  echo "--- $R"
  if [ -f "$R/depth_max.tif" ]; then
    echo "    envelope OK: $(ls $R/*.tif | wc -l) rasters, $(ls $R/*.vti 2>/dev/null | wc -l) snapshots"
    tail -3 "$R.log" | sed 's/^/    /'
  else
    echo "    *** NO depth_max.tif -- run did not complete cleanly"
    tail -5 "$R.log" | sed 's/^/    /'
  fi
done
