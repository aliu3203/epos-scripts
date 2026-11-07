#!/bin/bash
# correlators_all.sh -- run auau200-7/Correlators.C (v2, gamma, delta, kappa112 for charged pions vs the
# true RP, 30-40%) over every event of both EPOS 30-40% samples and combine them with
# auau200-7/Combine_Correlators.C.
#
# Correlators.C runs unmodified on chunks of whole files (~50k events), symlinked as z-auau_run_*.root
# to match the macro's glob. Chunks are independent subsamples, and the combine step takes all errors
# from the chunk-to-chunk spread. Correlators.C's own sub-sampling error is not usable here: it deals
# events round-robin into 20 subgroups, so with nfreeze 20 the 20 events of one hydro event land in
# 20 different subgroups and the subgroups are not independent.
#
# Per dataset: one sequential `cat` pulls each chunk into page cache (the fastest read pattern on the
# USB disk), then Correlators.C runs on it in the background, at most $NPAR at a time, while the next
# chunk is read. The two datasets are on different disks and run concurrently.
#
# Resumable: finished chunks are marked .done and skipped on restart; chunk lists are frozen at the
# first run (delete $OUT/<name>/chunks to rebuild them).
#
# Usage (detached from the terminal, survives logout):
#   mkdir -p /media/Students/aliu/correlators_all && cd ~ && setsid nohup /media/Students/epos-scripts/correlators_all.sh \
#       > /media/Students/aliu/correlators_all/daemon.log 2>&1 < /dev/null &

ANA=${ANA:-/media/Students/epos-scripts}
OUT=${OUT:-/media/Students/aliu/correlators_all}
BUILD=$OUT/build
NPAR=${NPAR:-2}

# name | source dir | file glob | files per chunk
DATASETS=(
    "auau200-7|/media/Students/aliu/auau200-7|z-*.root|50"
    "alpha_day1|/media/Disk_Alpha/EPOS/cen3040_200GeV_day1|auau_*.root|250"
)
MACROS="Correlators Combine_Correlators MakeFigure_Correlators MakeFigure_Kappa_STAR"

log() { echo "[$(date '+%F %T')] $*"; }

# root with the private build dir; usage: rootb <macro+args>
rootb() { root -l -b -q -e "gSystem->SetBuildDir(\"$BUILD\",true)" "$@"; }

build() {
    mkdir -p $BUILD
    for m in $MACROS; do
        rootb -e ".L $ANA/$m.C+" > $BUILD/$m.build.log 2>&1 || { log "build of $m failed, see $BUILD/$m.build.log"; exit 1; }
    done
    log "built macros into $BUILD (auau200-7 at git $(git -C $ANA rev-parse --short HEAD 2>/dev/null))"
}

# split the sorted file list into balanced chunks of symlinks named z-auau_run_*.root (what Correlators.C globs)
make_chunks() {
    local run=$1 src=$2 glob=$3 nper=$4
    [ -d $run/chunks ] && return
    local files=( $(cd $src && ls $glob | sort -V) )
    local n=${#files[@]} K=$(( (${#files[@]} + nper - 1) / nper ))
    for ((c = 0; c < K; c++)); do
        local d=$(printf "%s/chunks/c%03d" $run $c)
        mkdir -p $d
        for ((i = c*n/K; i < (c+1)*n/K; i++)); do
            local f=${files[$i]} b=${files[$i]#z-}
            b=${b#auau_run_}; b=${b#auau_}
            ln -s $src/$f $d/z-auau_run_$b
        done
    done
    log "$run: $n files -> $K chunks"
}

prefetch() { cat $1/z-auau_run_*.root > /dev/null 2>&1; }

run_chunk() {
    local d=$1
    (cd $d && stdbuf -oL -eL root -l -b -q -e "gSystem->SetBuildDir(\"$BUILD\",true)" "$ANA/Correlators.C+(5,0)" > Correlators.log 2>&1)
    if [ -s $d/cen5.correlators_job0.root ] && grep -q '^events used' $d/Correlators.log; then
        touch $d/.done
    else
        log "chunk $d FAILED, see $d/Correlators.log"
    fi
}

pipeline() {
    local name=$1 src=$2 glob=$3 nper=$4
    local run=$OUT/$name t0=$(date +%s) pids=()
    mkdir -p $run
    make_chunks $run $src "$glob" $nper
    local pending=()
    for d in $run/chunks/c*; do [ -e $d/.done ] || pending+=($d); done
    log "$name: ${#pending[@]} chunks to process"
    running() { local n=0 p; for p in ${pids[@]}; do kill -0 $p 2>/dev/null && n=$((n+1)); done; echo $n; }
    for ((i = 0; i < ${#pending[@]}; i++)); do
        prefetch ${pending[$i]}
        while [ $(running) -ge $NPAR ]; do sleep 10; done
        run_chunk ${pending[$i]} &
        pids+=($!)
        log "$name: $(basename ${pending[$i]}) cached and started ($((i+1))/${#pending[@]}, $(( $(date +%s) - t0 )) s)"
    done
    wait
    log "$name: finished, $(ls $run/chunks/c*/.done 2>/dev/null | wc -l)/$(ls -d $run/chunks/c* | wc -l) chunks ok, $(( $(date +%s) - t0 )) s"
}

mkdir -p $OUT
exec 9> $OUT/.lock
flock -n 9 || { log "another correlators_all.sh is running on $OUT"; exit 1; }

log "start: out=$OUT"
build
pids=()
for ds in "${DATASETS[@]}"; do
    IFS='|' read -r name src glob nper <<< "$ds"
    pipeline $name $src "$glob" $nper > $OUT/$name.log 2>&1 &
    pids+=($!)
done
wait ${pids[@]}

names=$(for ds in "${DATASETS[@]}"; do echo -n "${ds%%|*} "; done)
names=${names% }
rootb "$ANA/Combine_Correlators.C+(\"$OUT\",\"$names\")" > $OUT/combine.log 2>&1
# standard figures from the unmodified figure macros, on each merged file
for n in combined $names; do
    (cd $OUT/$n && rootb "$ANA/MakeFigure_Correlators.C+(5,0)" > MakeFigure_Correlators.log 2>&1
                   rootb "$ANA/MakeFigure_Kappa_STAR.C+(5,0)" > MakeFigure_Kappa_STAR.log 2>&1)
done
log "done: $OUT/summary.txt; merged files and figures in $OUT/{combined,${names// /,}}"
