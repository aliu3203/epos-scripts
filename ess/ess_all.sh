#!/bin/bash
# ess_all.sh -- event-shape-selection analysis (arXiv:2307.14997) over the four Disk_Alpha EPOS days.
#
#   pass 1: V2Const.C on every chunk -> hadd -> $OUT/v2const.root   (v2{RP}, v2pair{RP}, v2B{RP} for the q2 normalization)
#   pass 2: ESS.C on every chunk (reads v2const.root)              -> MakeFigure_ESS.C (combine, fits, figures, summary)
#
# Files are chained as in correlators/correlators_all.sh: each day is split into chunks of whole files (~50k events),
# symlinked as z-auau_run_*.root so the macros' TChain glob picks them up. Chunks are independent subsamples, which
# MakeFigure_ESS.C uses for the errors. Per chunk, one sequential `cat` pulls the files into page cache, then the
# macro runs on it in the background (at most $NPAR at a time) while the next chunk is read. All days are on one
# disk, so they go through one worker pool in sequence.
#
# Resumable: finished chunks are marked .done.<macro> and skipped on restart.
#
# Usage: mkdir -p /media/Students/aliu/ess_all && cd ~ && setsid nohup /media/Students/epos-scripts/ess/ess_all.sh \
#            > /media/Students/aliu/ess_all/daemon.log 2>&1 < /dev/null &

ANA=${ANA:-/media/Students/epos-scripts/ess}
OUT=${OUT:-/media/Students/aliu/ess_all}
SRC=${SRC:-/media/Disk_Alpha/EPOS}
DAYS=${DAYS:-"day1 day2 day3 day4"}
NPER=${NPER:-250}          # files per chunk (200 events each)
NPAR=${NPAR:-3}
BUILD=$OUT/build

log() { echo "[$(date '+%F %T')] $*"; }
rootb() { root -l -b -q -e "gSystem->SetBuildDir(\"$BUILD\",true)" "$@"; }

# split day $1 into chunks of symlinks $OUT/$1/chunks/cNNN/z-auau_run_*.root (frozen after the first run)
make_chunks() {
    local run=$OUT/$1 src=$SRC/cen3040_200GeV_$1
    [ -d $run/chunks ] && return
    local files=( $(cd $src && ls auau_*.root | sort -V) )
    local n=${#files[@]} K=$(( (${#files[@]} + NPER - 1) / NPER ))
    for ((c = 0; c < K; c++)); do
        local d=$(printf "%s/chunks/c%03d" $run $c)
        mkdir -p $d
        for ((i = c*n/K; i < (c+1)*n/K; i++)); do ln -s $src/${files[$i]} $d/z-auau_run_${files[$i]#auau_}; done
    done
    log "$1: $n files -> $K chunks"
}

# run macro $1 with arguments $2 in chunk directory $3
run_chunk() {
    local m=$1 args=$2 d=$3
    (cd $d && stdbuf -oL -eL root -l -b -q -e "gSystem->SetBuildDir(\"$BUILD\",true)" "$ANA/$m.C+($args)" > $m.log 2>&1)
    if grep -q '^events used' $d/$m.log; then touch $d/.done.$m; else log "chunk $d FAILED, see $d/$m.log"; fi
}

# run macro $1 (arguments $2) over every unfinished chunk of every day
pass() {
    local m=$1 args=$2 t0=$(date +%s)
    for day in $DAYS; do
        for d in $OUT/$day/chunks/c*; do
            [ -e $d/.done.$m ] && continue
            cat $d/z-auau_run_*.root > /dev/null 2>&1                    # prefetch into page cache
            while [ $(jobs -rp | wc -l) -ge $NPAR ]; do sleep 10; done
            run_chunk $m "$args" $d &
            log "$m: $day/$(basename $d) started ($(( $(date +%s) - t0 )) s)"
        done
    done
    wait
    local ntot=$(ls -d $OUT/*/chunks/c* | wc -l) ndone=$(ls $OUT/*/chunks/c*/.done.$m 2>/dev/null | wc -l)
    log "$m: $ndone/$ntot chunks done, $(( $(date +%s) - t0 )) s"
    [ $ndone -eq $ntot ] || { log "$m: not all chunks finished -- rerun to resume"; exit 1; }
}

mkdir -p $OUT $BUILD
exec 9> $OUT/.lock
flock -n 9 || { log "another ess_all.sh is running on $OUT"; exit 1; }
log "start: out=$OUT days=$DAYS (ess at git $(git -C $ANA rev-parse --short HEAD 2>/dev/null))"

for m in V2Const ESS MakeFigure_ESS; do
    rootb -e ".L $ANA/$m.C+" > $BUILD/$m.build.log 2>&1 || { log "build of $m failed, see $BUILD/$m.build.log"; exit 1; }
done
for day in $DAYS; do make_chunks $day; done

# ---- pass 1: normalization constants ----
pass V2Const 0
hadd -f $OUT/v2const.root $OUT/*/chunks/c*/cen5.V2Const_job0.root > $OUT/hadd_v2const.log 2>&1
rootb -e "TFile f(\"$OUT/v2const.root\"); TProfile* p=(TProfile*)f.Get(\"pV2\"); \
          printf(\"v2{RP} = %.5f  v2pair{RP} = %.5f  v2B{RP} = %.5f\n\", p->GetBinContent(1), p->GetBinContent(2), p->GetBinContent(3));" \
    | tee $OUT/v2const.txt
# constants changed (e.g. new days added) -> every chunk must be redone with them, not only the new ones
if ! cmp -s $OUT/v2const.txt $OUT/v2const.used; then
    rm -f $OUT/*/chunks/c*/.done.ESS
    cp $OUT/v2const.txt $OUT/v2const.used
fi

# ---- pass 2: event-shape selection ----
pass ESS "0,\"z-auau_run_*.root\",\"$OUT/v2const.root\""
rootb "$ANA/MakeFigure_ESS.C+(\"$OUT\",\"$DAYS\")" > $OUT/MakeFigure_ESS.log 2>&1
log "done: $OUT/ESS_summary.txt, $OUT/ESS_fig6.png, $OUT/ESS_fig5_9.png"
