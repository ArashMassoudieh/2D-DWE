#!/usr/bin/env bash
set -euo pipefail

solver=${1:-bin/2D-DWE}
workdir=$(mktemp -d)
trap 'rm -rf "$workdir"' EXIT

run_and_check() {
    local case_name=$1
    local log="$workdir/${case_name}.log"
    OMP_NUM_THREADS=2 "$solver" "$case_name" --out "$workdir/${case_name}" >"$log"

    if grep -q '\*\*\*' "$log"; then
        cat "$log"
        echo "$case_name verification reported a failure" >&2
        return 1
    fi

    grep -q 'mass balance .* OK' "$log"
    grep -q 'finished' "$log"
}

run_and_check plane
grep -q 'normal depth .* OK' "$workdir/plane.log"

run_and_check pond
grep -q 'draining to the notch' "$workdir/pond.log"
grep -q 'rim held' "$workdir/pond.log"

echo 'All analytical verification cases passed.'
