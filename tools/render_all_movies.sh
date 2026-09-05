#!/usr/bin/env bash
# Render the eight report animations, skipping only those already complete.
cd /home/arash/Projects/2D-DWE
L=tools/render_all_movies.out
: > $L
for r in out_c6_1m_rain out_c6_1m_berm_6hr out_c6_1m_24hr_v2 out_c6_1m_berm_24hr; do
  n=$(ls $r/c6_*.vti 2>/dev/null | wc -l); want=$(python3 -c "print($n/3.0)")
  for fld in depth speed; do
    f=$r/c6_coyote_canyon_$fld.mp4
    if [ -s "$f" ]; then
      got=$(ffprobe -v error -show_entries format=duration -of csv=p=0 "$f" 2>/dev/null || echo 0)
      if [ "$(python3 -c "print(1 if abs($got-$want)<0.5 else 0)")" = "1" ]; then
        echo "$(date '+%T')  ok   $f (${got}s)" >> $L; continue
      fi
      echo "$(date '+%T')  redo $f (${got}s, want ${want}s)" >> $L; rm -f "$f"
    fi
    echo "$(date '+%T')  render $r $fld" >> $L
    python3 tools/make_movie.py $r 1 --field $fld >> $L 2>&1
  done
done
echo "$(date '+%T')  ALL DONE" >> $L
