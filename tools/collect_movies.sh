#!/usr/bin/env bash
# Copy the eight report animations into one folder under self-describing names.
# Source names are identical across runs, so they cannot share a folder as-is.
cd /home/arash/Projects/2D-DWE
DEST=/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097/Report_for_C6/Animations
mkdir -p "$DEST"
declare -A SC=( [out_c6_1m_rain]=existing [out_c6_1m_berm_6hr]=berm
                [out_c6_1m_24hr_v2]=existing [out_c6_1m_berm_24hr]=berm )
declare -A ST=( [out_c6_1m_rain]=6hr [out_c6_1m_berm_6hr]=6hr
                [out_c6_1m_24hr_v2]=24hr [out_c6_1m_berm_24hr]=24hr )
declare -A FL=( [depth]=depth [speed]=velocity )
ok=0; bad=0
for r in "${!SC[@]}"; do
  n=$(ls $r/c6_*.vti 2>/dev/null | wc -l); want=$(python3 -c "print($n/3.0)")
  for fld in depth speed; do
    src=$r/c6_coyote_canyon_$fld.mp4
    dst="$DEST/NM15-097_C6_2D_${SC[$r]}_${ST[$r]}_${FL[$fld]}.mp4"
    if [ ! -s "$src" ]; then echo "MISSING  $src"; bad=$((bad+1)); continue; fi
    got=$(ffprobe -v error -show_entries format=duration -of csv=p=0 "$src" 2>/dev/null || echo 0)
    if [ "$(python3 -c "print(1 if abs($got-$want)<0.5 else 0)")" != "1" ]; then
      echo "TRUNCATED $src (${got}s, want ${want}s) -- not copied"; bad=$((bad+1)); continue
    fi
    cp "$src" "$dst"; ok=$((ok+1))
    printf "  %-52s %5.1fs\n" "$(basename "$dst")" "$got"
  done
done
echo "copied $ok, problems $bad  ->  $DEST"
