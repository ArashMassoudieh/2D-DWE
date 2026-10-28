#!/usr/bin/env bash
# Regenerate every 2-D figure and animation from the CONVERGED runs and copy
# them to the report / delivery folders. Supersedes the *_postC10 products,
# which were driven by the pass-1 (seeded Tc) hydrographs.
set -u
cd /home/arash/Projects/2D-DWE
S=/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-141Local
FIG=$S/Report_for_C10/Figures
ANI=$S/output_2D/animation
L=tools/update_2d_products_conv.out
: > $L
say(){ echo "$(date '+%T')  $*" | tee -a $L; }

TS=$(date +%Y%m%d)
BK=$FIG/_preTcIter_$TS; mkdir -p "$BK"
for f in "$FIG"/NM15-141_2d_map*.png "$FIG"/NM15-141_2d_profiles*.png \
         "$FIG"/NM15-141_2d_band_*.png "$FIG"/NM15-141_2d_subdivisions*.png; do
  [ -e "$f" ] && cp -n "$f" "$BK"/ 2>/dev/null
done
say "backed up $(ls "$BK" | wc -l) figures to $(basename "$BK")"
BKA=$ANI/_preTcIter_$TS; mkdir -p "$BKA"
for f in "$ANI"/C10_NM15-141_1m_*.mp4; do [ -e "$f" ] && cp -n "$f" "$BKA"/ 2>/dev/null; done
say "backed up $(ls "$BKA" 2>/dev/null | wc -l) animations"

#      run-dir                            report-suffix  storm  culvert-label        dem
SET="out_nm141_1m_6hr_carved_2m_conv||6hr|culvert-open|dem_2d_carved_2m.tif
out_nm141_1m_6hr_raw_conv|_blocked|6hr|culvert-blocked|dem_2d_raw.tif
out_nm141_1m_24hr_carved_2m_conv|_24hr|24hr|culvert-open|dem_2d_carved_2m.tif
out_nm141_1m_24hr_raw_conv|_24hr_blocked|24hr|culvert-blocked|dem_2d_raw.tif"

while IFS='|' read -r run suf storm culv dem; do
  [ -n "$run" ] || continue
  say "=== $run"
  python3 tools/plot_profiles.py     "$run" 1 >> $L 2>&1 || say "  FAIL profiles $run"
  python3 tools/plot_depth_bands.py  "$run" 1 >> $L 2>&1 || say "  FAIL bands $run"
  python3 tools/plot_subdivisions.py "$run" 1 >> $L 2>&1 || say "  FAIL subdiv $run"
  cp "$run/nm141_map.png"              "$FIG/NM15-141_2d_map${suf}.png"
  cp "$run/nm141_profiles.png"         "$FIG/NM15-141_2d_profiles${suf}.png"
  cp "$run/c6_subdivision_flooding.png" "$FIG/NM15-141_2d_subdivisions${suf}.png"
  for b in 1 3 6; do cp "$run/c6_band_${b}ft.png" "$FIG/NM15-141_2d_band_${b}ft${suf}.png"; done
  say "  figures copied"

  case $storm in 6hr) scap="100-yr 6-hr event";; 24hr) scap="100-yr 24-hr event";; esac
  case $culv in culvert-open) ccap="culvert open (2 m notch)";; *) ccap="culvert blocked";; esac
  for fld in depth speed; do
    name=$([ $fld = speed ] && echo velocity || echo depth)
    out="$ANI/C10_NM15-141_1m_${storm}_${culv}_${name}.mp4"
    say "  animation $fld -> $(basename "$out")"
    python3 tools/make_movie.py "$run" 1 --field $fld \
      --dem "$S/output_2D/$dem" --title "C10 Coyote Canyon — $ccap" \
      --storm "$scap" --out "$out" >> $L 2>&1 || say "  FAIL movie $run $fld"
  done
done <<< "$SET"
say "ALL DONE"
