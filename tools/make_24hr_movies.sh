#!/usr/bin/env bash
cd /home/arash/Projects/2D-DWE
L=tools/make_24hr_movies.out
: > $L
for r in out_c6_1m_24hr_v2 out_c6_1m_berm_24hr; do
  for fld in depth speed; do
    f=$r/c6_coyote_canyon_$fld.mp4
    # A non-empty file is NOT proof of a complete render: a run killed mid-write
    # leaves a short, playable mp4.  Require the duration to match the frame count.
    want=$(python3 -c "import glob;print(len(glob.glob('$r/c6_*.vti'))/3.0)")
    if [ -s "$f" ]; then
      got=$(ffprobe -v error -show_entries format=duration -of csv=p=0 "$f" 2>/dev/null || echo 0)
      ok=$(python3 -c "print(1 if abs($got-$want)<0.5 else 0)")
      if [ "$ok" = "1" ]; then echo "$(date '+%T')  skip $f (complete, ${got}s)" >> $L; continue; fi
      echo "$(date '+%T')  re-render $f (truncated: ${got}s, expected ${want}s)" >> $L
      rm -f "$f"
    fi
    echo "$(date '+%T')  rendering $r $fld" >> $L
    python3 tools/make_movie.py $r 1 --field $fld >> $L 2>&1
  done
done
echo "$(date '+%T')  done" >> $L
ls -la out_c6_1m_24hr_v2/*.mp4 out_c6_1m_berm_24hr/*.mp4 >> $L 2>&1
