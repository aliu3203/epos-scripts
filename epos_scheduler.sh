#!/bin/bash
# Run cycle_epos.sh every 7 hours, forever.
#
# Start it (from a terminal where ROOT is set up):
#   nohup ./epos_scheduler.sh auau200-7 > scheduler.log 2>&1 &
# Watch it:   tail -f scheduler.log
# Stop it:    pkill -f epos_scheduler.sh   (running simulations keep going)

DEST=${1:-auau200-7}
INTERVAL=${2:-7h}

while true; do
    echo "=== cycle started $(date '+%F %T') ==="
    /home/aliu/cycle_epos.sh "$DEST"
    echo "=== sleeping $INTERVAL ==="
    sleep "$INTERVAL"
done
