#!/bin/bash
# Run all 5 EPOS simulations in parallel, each doing NBATCH consecutive runs.
# After each run completes, its root file is renamed to the next free index
# and moved into the destination folder, then the next run starts immediately.
#
# Usage:  ./run_all.sh [runs-per-simulation] [destination-folder]
# Example: ./run_all.sh 3 auau200-7   -> 15 root files total in auau200-7

NBATCH=${1:-1}
DEST=${2:-auau200-7}

JINBASE=/media/Disk_Jin/aliu
HOMEBASE=/home/aliu

mkdir -p "$JINBASE/$DEST"

run_sim() {
    local n=$1
    for ((b = 1; b <= NBATCH; b++)); do
        echo "[epos$n] starting run $b/$NBATCH  ($(date '+%F %T'))"

        # Subshell isolates this simulation's environment (startup.sh exports)
        (
            cd "$HOMEBASE/epos$n" || exit 1
            source startup.sh
            cd epos4.0.3/auau200 || exit 1
            "$BIN_DIR/bin/epos" -root "auau_run_$n" > "$JIN/epos$n/output.log" 2>&1
        )

        local src="$JINBASE/epos$n/z-auau_run_$n.root"
        if [[ -f $src ]]; then
            # flock: only one simulation picks an index at a time, so two
            # finishing simultaneously can't grab the same number
            (
                flock 9
                i=1
                while [[ -e "$JINBASE/$DEST/z-auau_run_$i.root" ]]; do ((i++)); done
                mv "$src" "$JINBASE/$DEST/z-auau_run_$i.root"
                echo "[epos$n] run $b done -> $DEST/z-auau_run_$i.root  ($(date '+%F %T'))"
            ) 9>"$JINBASE/$DEST/.movelock"
        else
            echo "[epos$n] run $b: no root file produced, see $JINBASE/epos$n/output.log" >&2
        fi
    done
}

for n in 1 2 3 4 5; do
    run_sim "$n" &
done
wait
echo "All simulations finished."
