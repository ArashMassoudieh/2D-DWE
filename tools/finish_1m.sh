#!/bin/bash
# Wait for the in-flight 1 m run to exit, then regenerate every 1 m product.
# Writes ONLY into out_c6_1m/ -- the report's Figures/ are left alone, because
# the choice between the 1 m and 3 m result has not been made yet.
cd /home/arash/Projects/2D-DWE
until ! ps -p "$1" > /dev/null 2>&1; do sleep 60; done
echo "=== run finished $(date) ==="
tail -12 out_c6_1m.log
echo "=== envelope + barrier profile ==="
python3 tools/plot_c6.py         out_c6_1m 1 2>&1 | grep -vE "Warning|memory leak"
echo "=== hazard bands ==="
python3 tools/plot_depth_bands.py out_c6_1m 1 2>&1 | grep -vE "Warning|memory leak"
echo "=== animations ==="
python3 tools/make_movie.py out_c6_1m 1 --dt-frame 3600 --fps 2 2>&1 | tail -2
python3 tools/make_movie.py out_c6_1m 1 --dt-frame 3600 --fps 2 --field speed 2>&1 | tail -2
echo "=== done $(date) ==="
ls -la out_c6_1m/*.png out_c6_1m/*.mp4 out_c6_1m/*.tif 2>/dev/null
