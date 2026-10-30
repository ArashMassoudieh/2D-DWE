#!/usr/bin/env bash
set -euo pipefail

solver=${1:-bin/2D-DWE}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
solver=$(realpath "$solver")
case_dir="$repo_root/cases/ea_test4"
workdir=$(mktemp -d)
trap 'rm -rf "$workdir"' EXIT

cp "$case_dir/config.txt" "$case_dir/bc_test4.csv" "$case_dir/make_case.py" "$workdir/"
(
    cd "$workdir"
    python3 make_case.py
    OMP_NUM_THREADS=2 "$solver" --config config.txt --out out >run.log
)

grep -q 'mass balance .* OK' "$workdir/run.log"
grep -q 'finished' "$workdir/run.log"

python3 - "$workdir/out/gauges.csv" <<'PY'
import csv
import sys

path = sys.argv[1]
with open(path, newline="") as handle:
    rows = list(csv.DictReader(line for line in handle if not line.startswith("#")))

def nearest(point, hour):
    target = hour * 3600.0
    row = min(rows, key=lambda r: abs(float(r["t_s"]) - target))
    return float(row[f"{point}_depth_m"])

def arrival(point, threshold=0.002):
    for row in rows:
        if float(row[f"{point}_depth_m"]) > threshold:
            return float(row["t_s"]) / 60.0
    raise AssertionError(f"{point} never became wet")

arrival_ranges = {
    "P1": (8.0, 11.0),
    "P3": (25.0, 31.0),
    "P5": (52.0, 63.0),
    "P6": (55.0, 65.0),
}
depth_ranges = {
    "P1": {2: (0.30, 0.33), 3: (0.30, 0.33), 4: (0.30, 0.33)},
    "P3": {2: (0.20, 0.22), 4: (0.21, 0.23), 5: (0.147, 0.158)},
    "P5": {2: (0.14, 0.155), 3: (0.16, 0.172), 4: (0.167, 0.18), 5: (0.146, 0.155)},
    "P6": {2: (0.137, 0.151), 3: (0.152, 0.168), 4: (0.162, 0.176), 5: (0.145, 0.154)},
}

for point, limits in arrival_ranges.items():
    value = arrival(point)
    assert limits[0] <= value <= limits[1], (point, "arrival_min", value, limits)

for point, samples in depth_ranges.items():
    for hour, limits in samples.items():
        value = nearest(point, hour)
        assert limits[0] <= value <= limits[1], (point, hour, value, limits)

print("EA Test 4 arrival times and depths are within the published benchmark spread.")
PY
