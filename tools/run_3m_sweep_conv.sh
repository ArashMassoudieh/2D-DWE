#!/usr/bin/env bash
# Notch sweep + grid-convergence runs at the 10 ft (3 m) grid on the CONVERGED
# hydrology, to replace the pre-culvert-10 tables 10 and 11.
set -u
cd /home/arash/Projects/2D-DWE
S=/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-141Local/output_2D
L=tools/run_3m_sweep_conv.out; : > $L
say(){ echo "$(date '+%T')  $*" | tee -a $L; }
SET="dem_2d_raw.tif|none
dem_2d_carved_1m.tif|1m
dem_2d_carved_2m.tif|2m
dem_2d_carved_3m.tif|3m
dem_2d_carved_4m.tif|4m"
pids=()
while IFS='|' read -r dem tag; do
  [ -n "$dem" ] || continue
  out="out_nm141_3m_${tag}_conv"; rm -rf "$out"; mkdir -p "$out"
  say "launch $out ($dem)"
  OMP_NUM_THREADS=6 nohup ./bin/2D-DWE nm141 --dx 3 --storm 6hr --dem "$S/$dem" \
      --out "$out" --tout 1800 --cfl 0.70 > "$out.log" 2>&1 &
  pids+=($!)
done <<< "$SET"
say "5 runs launched"
for p in "${pids[@]}"; do wait "$p" || say "pid $p non-zero"; done
say "ALL DONE"
for d in out_nm141_3m_*_conv; do
  say "$d: $(grep -c '^finished' $d.log) finished | $(grep -A4 'depth bands' $d.log | grep -oE '[0-9.]+ acres' | tr '\n' ' ')"
done
