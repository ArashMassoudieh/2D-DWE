#!/usr/bin/env bash
set -u
cd /home/arash/Projects/2D-DWE
for r in out_c6_1m_24hr_v2 out_c6_1m_berm_24hr; do
  for fld in depth speed; do
    echo "=== $r / $fld"
    if [ "$fld" = "depth" ]; then
      python3 tools/make_movie.py "$r" 1 2>&1 | grep -viE "memory leak" | tail -2
    else
      python3 tools/make_movie.py "$r" 1 --field speed 2>&1 | grep -viE "memory leak" | tail -2
    fi
  done
done
echo "=== verifying and collecting all eight"
bash tools/collect_movies.sh
