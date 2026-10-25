#!/usr/bin/env bash
# vsb_survey.sh -- one second of each UHF ATSC 1.0 channel on the yagi, then
# rank them by pilot SNR. Waits for the radio lock to be free first.
cd "$(dirname "$0")/.." || exit 1
PY=/c/Users/emane/radioconda/python.exe
until python -c "import sys; sys.path.insert(0, r'Z:\src\gr-radiotuna\tools'); import radio_lock; sys.exit(0 if radio_lock.status() is None else 1)"; do
  sleep 60
done
echo "radio free at $(date -u +%H:%M:%SZ)"
# UHF 8-VSB channels per the FCC TV query, 2026-10-24 (post-repack)
for rf in 15 21 22 26 27 31 34 35 36; do
  fc=$(( 470 + (rf-14)*6 + 3 ))e6
  $PY tools/iq_capture.py --fc $fc --rate 8e6 --secs 1 --ant "Antenna B" --rfgain 4 --ifgr 32 \
      --out lab/vsb_survey_rf$rf.cs16 2>&1 | grep -E "lock|headroom|VERDICT" | tr '\n' ' ' | sed "s/^/rf$rf: /"
  echo
done
python tools/vsb_ref.py survey lab/vsb_survey_rf*.cs16
echo SURVEY DONE
