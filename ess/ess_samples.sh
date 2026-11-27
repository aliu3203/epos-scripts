#!/bin/bash
# ess_samples.sh -- run ess_all.sh on every EPOS sample on the data disk, results organized as $ROOT/<energy>/<centrality>/.
#
# A sample is the set of folders $SRC/cen<lo><hi>_<E>GeV_dayN; e.g. cen3040_27GeV_day1 -> $ROOT/27GeV/30-40.
# Every finished day is used. A day counts as finished once its newest file is older than $MINAGE minutes, so a day
# still being produced is left for the next invocation (ess_all.sh freezes a day's file list when it first chunks it).
# The impact-parameter window of each centrality is the production window, from the table below (checked against
# the files: 30-40% 8.4-9.2 fm at 200 and 27 GeV, 70-80% 12.5-13.4 fm).
# Samples run one after another, smallest first: all are on one disk and every pass is disk-bound. Run it again
# after new days appear: per sample, pass 1 runs only on the new chunks, then pass 2 is redone with the updated
# constants (see ess_all.sh).
#
# Usage: mkdir -p /media/Students/aliu/ess_results && cd ~ && setsid nohup /media/Students/epos-scripts/ess/ess_samples.sh \
#            >> /media/Students/aliu/ess_results/samples.log 2>&1 < /dev/null &
#        DRYRUN=1 /media/Students/epos-scripts/ess/ess_samples.sh      # only print what would run

ANA=${ANA:-$(cd $(dirname $0) && pwd)}
ROOT=${ROOT:-/media/Students/aliu/ess_results}
SRC=${SRC:-/media/Disk_Alpha/EPOS}
MINAGE=${MINAGE:-360}
export NPAR=${NPAR:-4}

declare -A BWIN=( [3040]="8.4 9.2" [7080]="12.5 13.4" )     # centrality -> impact-parameter window (fm) of the production

log() { echo "[$(date '+%F %T')] $*"; }

# samples, smallest (total bytes) first
for s in $(cd $SRC && ls -d cen*_*GeV_day* | sed 's/_day[0-9]*$//' | sort -u); do
    echo "$(du -s --apparent-size $SRC/${s}_day* | awk '{t += $1} END {print t}') $s"
done | sort -n | while read size s; do
    c=${s#cen}; c=${c%%_*}; e=${s#*_}                          # e.g. 3040, 27GeV
    cent=${c:0:2}-${c:2:2}
    [ -n "${BWIN[$c]}" ] || { log "$s: no impact-parameter window for centrality $cent in ess_samples.sh -- skipped"; continue; }
    days=""
    for d in $(cd $SRC && ls -d ${s}_day* | sort -V); do
        if [ -n "$(find $SRC/$d -name '*.root' -mmin -$MINAGE -print -quit)" ]; then
            log "$d: still being written -- left for the next run"; continue
        fi
        days="$days ${d##*_}"
    done
    [ -n "$days" ] || continue
    out=$ROOT/$e/$cent
    set -- ${BWIN[$c]}
    log "$s -> $out (days:$days; b = $1-$2 fm)"
    [ -n "$DRYRUN" ] && continue
    mkdir -p $out
    OUT=$out SAMPLE=$s DAYS="${days# }" BMIN=$1 BMAX=$2 CENT="$cent%" ENERGY="${e%GeV} GeV" ANA=$ANA \
        $ANA/ess_all.sh >> $out/daemon.log 2>&1 < /dev/null || log "$s: ess_all.sh stopped early, see $out/daemon.log"
    log "$s: $(tail -1 $out/daemon.log)"
done
log "all samples done"
