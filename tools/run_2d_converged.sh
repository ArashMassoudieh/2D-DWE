#!/usr/bin/env bash
# Re-run the four C10 2-D cases on the CONVERGED HMS hydrographs.
#
# The four *_postC10 runs were driven by the pass-1 (seeded Tc) hydrographs.
# Carrying the Tc-Q iteration to convergence raised the boundary peaks 16-23%,
# so every 2-D result derived from them is superseded. Same terrain, same
# solver, same settings -- only the inflow CSVs changed.
#
# Frames stay OUT of Dropbox (1.9 GB per run); only the GeoTIFF envelopes are
# copied over afterwards.
set -u
cd /home/arash/Projects/2D-DWE
S=/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-141Local/output_2D
L=tools/run_2d_converged.out
: > $L
say(){ echo "$(date '+%T')  $*" | tee -a $L; }

#      storm  tout   dem                    suffix
RUNS="6hr|1800|dem_2d_carved_2m.tif|carved_2m
6hr|1800|dem_2d_raw.tif|raw
24hr|7200|dem_2d_carved_2m.tif|carved_2m
24hr|7200|dem_2d_raw.tif|raw"

pids=()
while IFS='|' read -r storm tout dem suf; do
  [ -n "$storm" ] || continue
  out="out_nm141_1m_${storm}_${suf}_conv"
  rm -rf "$out"; mkdir -p "$out"
  say "launch $out  (storm $storm, tout $tout, dem $dem)"
  OMP_NUM_THREADS=7 nohup ./bin/2D-DWE nm141 --dx 1 --storm "$storm" \
      --dem "$S/$dem" --out "$out" --tout "$tout" --cfl 0.70 \
      > "$out.log" 2>&1 &
  pids+=($!)
done <<< "$RUNS"

say "4 runs launched, pids: ${pids[*]}"
for p in "${pids[@]}"; do wait "$p" || say "pid $p exited non-zero"; done
say "ALL RUNS FINISHED"
for d in out_nm141_1m_*_conv; do
  say "$d: $(ls $d/*.vti 2>/dev/null | wc -l) frames, $(du -sh $d | cut -f1)"
done
