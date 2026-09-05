#!/usr/bin/env bash
# Render the eight C10 (NM15-141) report animations.
#
# Four runs x two fields. Every file is named for the three things that
# distinguish it -- storm duration, culvert representation, and field -- so the
# eight can share one folder and be told apart WITHOUT their parent directory.
# The old names could not: three of the four runs wrote "..._depth.mp4".
#
# The hillshade DEM must match the run. Rendering a blocked-culvert run over the
# carved terrain draws water standing below the road notch, which reads as a
# solver error rather than the intended comparison.
set -u
cd /home/arash/Projects/2D-DWE
SITE=/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-141Local/output_2D
L=tools/make_c10_movies.out
: > $L

#      run                              storm  culvert   dem
RUNS="out_nm141_1m_v3_raw|6hr|culvert-blocked|dem_2d_raw.tif
out_nm141_1m_v3_carved_2m|6hr|culvert-open|dem_2d_carved_2m.tif
out_nm141_1m_24hr_raw|24hr|culvert-blocked|dem_2d_raw.tif
out_nm141_1m_24hr_carved_2m|24hr|culvert-open|dem_2d_carved_2m.tif"

say() { echo "$(date '+%T')  $*" | tee -a $L; }

while IFS='|' read -r run storm culv dem; do
  [ -n "$run" ] || continue
  n=$(ls $run/nm141_*.vti 2>/dev/null | wc -l)
  want=$(python3 -c "print($n/3.0)")          # fps 3
  case $storm in
    6hr)  scap="100-yr 6-hr event"  ;;
    24hr) scap="100-yr 24-hr event" ;;
  esac
  case $culv in
    culvert-open)    ccap="culvert open (2 m notch)" ;;
    culvert-blocked) ccap="culvert blocked" ;;
  esac
  for fld in depth speed; do
    name=$([ $fld = speed ] && echo velocity || echo depth)
    out="$run/C10_NM15-141_1m_${storm}_${culv}_${name}.mp4"
    # A non-empty mp4 is NOT proof of a finished render: a killed run leaves a
    # short, playable file. Require the duration to match the frame count.
    if [ -s "$out" ]; then
      got=$(ffprobe -v error -show_entries format=duration -of csv=p=0 "$out" 2>/dev/null || echo 0)
      if [ "$(python3 -c "print(1 if abs($got-$want)<0.5 else 0)")" = "1" ]; then
        say "ok    $out (${got}s)"; continue
      fi
      say "redo  $out (${got}s, want ${want}s)"; rm -f "$out"
    fi
    say "render $run $fld -> $(basename $out)"
    python3 tools/make_movie.py "$run" 1 --field $fld \
      --dem "$SITE/$dem" \
      --title "C10 Coyote Canyon — $ccap" \
      --storm "$scap" \
      --out "$out" >> $L 2>&1 || say "FAILED $out"
  done
done <<< "$RUNS"

say "ALL DONE"
ls -la out_nm141_1m_*/C10_*.mp4 >> $L 2>&1
