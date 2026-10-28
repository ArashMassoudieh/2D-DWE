#!/usr/bin/env bash
# Regenerate the C10 report figures from the POST-culvert-10 runs and copy them
# into Report_for_C10/Figures under the report names.
#
# The four post-culvert-10 runs live under the GIS project, not here:
#   output_2D/2D_Results/out_nm141_1m_{6hr,24hr}_{carved_2m,raw}_postC10
# Each already holds depth_max.tif / wse_max.tif, so nothing is re-solved.
set -u
cd /home/arash/Projects/2D-DWE
S=/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-141Local
R=$S/output_2D/2D_Results
FIG=$S/Report_for_C10/Figures
L=tools/make_c10_figures.out
: > $L
say(){ echo "$(date '+%T')  $*" | tee -a $L; }

#     run-folder                             report-suffix
SET="out_nm141_1m_6hr_carved_2m_postC10|
out_nm141_1m_6hr_raw_postC10|_blocked
out_nm141_1m_24hr_carved_2m_postC10|_24hr
out_nm141_1m_24hr_raw_postC10|_24hr_blocked"

# Keep the superseded pre-culvert-10 figures rather than overwriting them away.
BK=$FIG/_pre_culvert10_$(date +%Y%m%d)
mkdir -p "$BK"
for f in "$FIG"/NM15-141_2d_map*.png "$FIG"/NM15-141_2d_profiles*.png "$FIG"/NM15-141_2d_band_*.png; do
  [ -e "$f" ] && cp -n "$f" "$BK"/ 2>/dev/null
done
say "backed up $(ls "$BK" | wc -l) superseded figures to $(basename "$BK")"

while IFS='|' read -r run suf; do
  [ -n "$run" ] || continue
  d=$R/$run
  if [ ! -s "$d/depth_max.tif" ]; then say "MISSING $d/depth_max.tif"; continue; fi
  say "=== $run  (suffix '${suf:-<none>}')"
  say "  profiles + map"
  python3 tools/plot_profiles.py    "$d" 1 >> $L 2>&1 || say "  FAILED plot_profiles $run"
  say "  depth bands"
  python3 tools/plot_depth_bands.py "$d" 1 >> $L 2>&1 || say "  FAILED plot_depth_bands $run"
  say "  subdivisions"
  python3 tools/plot_subdivisions.py "$d" 1 >> $L 2>&1 || say "  FAILED plot_subdivisions $run"

  cp -v "$d/nm141_map.png"      "$FIG/NM15-141_2d_map${suf}.png"      >> $L 2>&1
  cp -v "$d/nm141_profiles.png" "$FIG/NM15-141_2d_profiles${suf}.png" >> $L 2>&1
  for b in 1 3 6; do
    cp -v "$d/c6_band_${b}ft.png" "$FIG/NM15-141_2d_band_${b}ft${suf}.png" >> $L 2>&1
  done
done <<< "$SET"

say "ALL DONE"
