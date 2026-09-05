#!/usr/bin/env bash
set -u
cd /home/arash/Projects/2D-DWE
for r in out_c6_1m_rain out_c6_1m_berm_6hr; do
  for fld in depth speed; do
    echo "=== $r / $fld"
    if [ "$fld" = "depth" ]; then
      python3 tools/make_movie.py "$r" 1 2>&1 | grep -viE "memory leak" | tail -3
    else
      python3 tools/make_movie.py "$r" 1 --field speed 2>&1 | grep -viE "memory leak" | tail -3
    fi
  done
done
echo "=== rendering done; verifying and collecting"
bash tools/collect_movies.sh
