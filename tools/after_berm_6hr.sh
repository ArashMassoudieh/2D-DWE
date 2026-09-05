#!/usr/bin/env bash
# Wait for the berm 6-hour run, post-process it, then queue BOTH 24-hour runs
# (no-berm baseline and with-berm scenario) in parallel at 16 threads each.
set -u
cd /home/arash/Projects/2D-DWE
RUN=out_c6_1m_berm_6hr
DEMB=/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097/outputs_2D/cliped_dem_berm.tif
LOG=tools/after_berm_6hr.out
say(){ echo "$(date '+%F %T')  $*" >> $LOG; }

: > $LOG
say "waiting for the 6-hour solver, PID 400443"
while kill -0 400443 2>/dev/null; do sleep 60; done
say "solver exited"
tail -3 ${RUN}.log >> $LOG

if [ ! -f ${RUN}/depth_max.tif ]; then
  say "*** NO ENVELOPE RASTERS -- the 6-hour run did not complete cleanly."
  say "*** Not post-processing and NOT queueing the 24-hour runs. Investigate."
  exit 1
fi
say "envelope rasters present:"; ls -la ${RUN}/*.tif >> $LOG

say "max-depth map";   python3 tools/plot_c6.py           $RUN 1 $DEMB >> $LOG 2>&1
say "hazard bands";    python3 tools/plot_depth_bands.py  $RUN 1 $DEMB >> $LOG 2>&1
say "subdivisions";    python3 tools/plot_subdivisions.py $RUN 1 $DEMB "C6 Coyote Canyon -- WITH BERM" >> $LOG 2>&1
say "depth animation"; python3 tools/make_movie.py        $RUN 1 >> $LOG 2>&1
say "speed animation"; python3 tools/make_movie.py        $RUN 1 --field speed >> $LOG 2>&1
say "post-processing done"
ls -la ${RUN}/*.png ${RUN}/*.mp4 2>/dev/null >> $LOG

# ---- both 24-hour runs, in parallel, 16 OpenMP threads each -----------------
# Parallel rather than sequential: the pair finishes in about the time one run
# would take alone at full width, and if anything kills them both have snapshots
# rather than one being complete and the other never started.
say "launching the 24-hour BASELINE (no berm)"
rm -rf out_c6_1m_24hr_v2
OMP_NUM_THREADS=16 nohup ./bin/2D-DWE c6 --dx 1 --storm 24hr \
      --tend 129600 --tout 7200 --out out_c6_1m_24hr_v2 \
      > out_c6_1m_24hr_v2.log 2>&1 &
say "  baseline pid $!"

sleep 20
say "launching the 24-hour WITH BERM"
rm -rf out_c6_1m_berm_24hr
OMP_NUM_THREADS=16 nohup ./bin/2D-DWE c6 --dx 1 --storm 24hr --dem $DEMB \
      --tend 129600 --tout 7200 --out out_c6_1m_berm_24hr \
      > out_c6_1m_berm_24hr.log 2>&1 &
say "  berm pid $!"
say "both 24-hour runs queued"
