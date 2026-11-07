#!/bin/bash
# One cycle over simulations 1-5:
#   1. skip any simulation whose epos process is still running
#   2. move its finished z-auau_run_{n}.root into DEST with the next free index
#   3. start its next run (unless RESTART=0)
#
# Usage:  ./cycle_epos.sh [destination-folder]
#         RESTART=0 ./cycle_epos.sh auau200-7   # move files only, no new runs

DEST=${1:-auau200-7}
RESTART=${RESTART:-1}

JINBASE=/media/Students/aliu
HOMEBASE=/home/aliu

mkdir -p "$JINBASE/$DEST"

for n in 1 2 3 4 5; do
    # Safety: never touch a root file while its epos process is alive
    if pgrep -f "epos -root auau_run_$n$" > /dev/null; then
        echo "[epos$n] still running - skipped this cycle"
        continue
    fi

    src="$JINBASE/epos$n/z-auau_run_$n.root"
    if [[ -f $src ]]; then
        i=1
        while [[ -e "$JINBASE/$DEST/z-auau_run_$i.root" ]]; do ((i++)); done
        mv "$src" "$JINBASE/$DEST/z-auau_run_$i.root"
        echo "[epos$n] moved root file -> $DEST/z-auau_run_$i.root"
    else
        echo "[epos$n] no root file to move"
    fi

    if [[ $RESTART == 1 ]]; then
        # Subshell isolates this simulation's startup.sh environment
        (
            cd "$HOMEBASE/epos$n" || exit 1
            source startup.sh
            cd epos4.0.3/auau200 || exit 1
            nohup "$BIN_DIR/bin/epos" -root "auau_run_$n" > "$JIN/epos$n/output.log" 2>&1 &
        )
        echo "[epos$n] new run started"
    fi
done
