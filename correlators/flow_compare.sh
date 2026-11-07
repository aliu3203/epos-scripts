#!/bin/bash
# flow_compare.sh -- run the auau200-7 flow analyses (Elliptic, Elliptic_KP, Elliptic_EP, Elliptic_ESS)
# on two EPOS datasets and compare them with auau200-7/Compare_Datasets.C.
#
# Each dataset is split into chunks of whole files (~50k events). Files never share a hydro event,
# so chunks are independent subsamples; Compare_Datasets.C takes its errors from the chunk spread
# (nfreeze 20: every hydro event is hadronised 20x, so per-particle TProfile errors are too small).
#
# Per chunk: one sequential `cat` pulls the files into page cache (the fastest read pattern on the
# USB disks), then the four unmodified macros run on it in parallel from cache while the next chunk
# prefetches. Both datasets are on different physical disks and run concurrently.
# The macros are compiled into $OUT/build, leaving the auau200-7 ACLiC artifacts untouched.
#
# Resumable: finished chunks are marked .done and skipped on restart; chunk lists are frozen at the
# first run (delete $OUT/<name>/chunks to rebuild them).
#
# Usage:  cd ~ && nohup /media/Students/epos-scripts/flow_compare.sh > /media/Students/aliu/dataset_compare/daemon.log 2>&1 &

ANA=${ANA:-/media/Students/epos-scripts}
OUT=${OUT:-/media/Students/aliu/dataset_compare}
BUILD=$OUT/build

# name | source dir | file glob | files per chunk
REF_NAME=${REF_NAME:-auau200-7};  REF_SRC=${REF_SRC:-/media/Students/aliu/auau200-7};          REF_GLOB=${REF_GLOB:-'z-*.root'};    REF_CHUNK=${REF_CHUNK:-50}
NEW_NAME=${NEW_NAME:-alpha_day1}; NEW_SRC=${NEW_SRC:-/media/Disk_Alpha/EPOS/cen3040_200GeV_day1}; NEW_GLOB=${NEW_GLOB:-'auau_*.root'}; NEW_CHUNK=${NEW_CHUNK:-250}
REF_LABEL=${REF_LABEL:-"auau200-7 (local)"}
NEW_LABEL=${NEW_LABEL:-"Disk_Alpha day1"}

ANALYSES="Elliptic Elliptic_KP Elliptic_EP Elliptic_ESS"
OUTPUTS="v2pT_pion v2pT_kaon_proton v2pT_EP v2q2_ESS"
FIGURES="MakeFigure_Elliptic MakeFigure_Elliptic_KP MakeFigure_Elliptic_EP MakeFigure_Elliptic_ESS"

log() { echo "[$(date '+%F %T')] $*"; }

# root with the private build dir; usage: rootb <macro+args>
rootb() { root -l -b -q -e "gSystem->SetBuildDir(\"$BUILD\",true)" "$@"; }

build() {
    mkdir -p $BUILD
    for m in $ANALYSES $FIGURES Compare_Datasets; do
        rootb -e ".L $ANA/$m.C+" > $BUILD/$m.build.log 2>&1 || { log "build of $m failed, see $BUILD/$m.build.log"; exit 1; }
    done
    log "built macros into $BUILD (auau200-7 at git $(git -C $ANA rev-parse --short HEAD 2>/dev/null))"
}

# split the sorted file list into balanced chunks of symlinks named z-*.root (what the macros glob)
make_chunks() {
    local run=$1 src=$2 glob=$3 nper=$4
    [ -d $run/chunks ] && return
    local files=( $(cd $src && ls $glob | sort -V) )
    local n=${#files[@]} K=$(( (${#files[@]} + nper - 1) / nper ))
    for ((c = 0; c < K; c++)); do
        local d=$(printf "%s/chunks/c%03d" $run $c)
        mkdir -p $d
        for ((i = c*n/K; i < (c+1)*n/K; i++)); do
            local f=${files[$i]}
            ln -s $src/$f $d/z-${f#z-}
        done
    done
    log "$run: $n files -> $K chunks"
}

prefetch() { cat $1/z-*.root > /dev/null 2>&1; }

run_chunk() {
    local d=$1 pids=() ok=1
    for m in $ANALYSES; do
        (cd $d && stdbuf -oL -eL root -l -b -q -e "gSystem->SetBuildDir(\"$BUILD\",true)" "$ANA/$m.C+(5,0)" > $m.log 2>&1) &
        pids+=($!)
    done
    for p in ${pids[@]}; do wait $p || ok=0; done
    for o in $OUTPUTS; do [ -s $d/cen5.${o}_job0.root ] || ok=0; done
    if [ $ok = 1 ]; then touch $d/.done; else log "chunk $d FAILED, see $d/*.log"; fi
}

merge() {
    local run=$1
    local nchunk=$(ls -d $run/chunks/c*/.done 2>/dev/null | wc -l)
    for o in $OUTPUTS; do
        hadd -f $run/cen5.${o}_job0.root $(for f in $run/chunks/c*/.done; do echo ${f%/.done}/cen5.${o}_job0.root; done) > $run/hadd_$o.log 2>&1
    done
    # hadd sums hConst (v2{2}, v2pair{2}) over chunks; store the chunk average instead
    root -l -b -q -e "TFile f(\"$run/cen5.v2q2_ESS_job0.root\",\"UPDATE\"); TH1* h=(TH1*)f.Get(\"hConst\"); h->Scale(1./$nchunk); h->Write(\"\",TObject::kOverwrite);" > /dev/null 2>&1
    # per-dataset figures from the standard figure macros
    (cd $run && for m in $FIGURES; do rootb "$ANA/$m.C+(5,0)" > $m.log 2>&1; done)
    log "$run: merged $nchunk chunks, figures done"
}

pipeline() {
    local name=$1 src=$2 glob=$3 nper=$4
    local run=$OUT/$name
    mkdir -p $run
    make_chunks $run $src "$glob" $nper
    local pending=()
    for d in $run/chunks/c*; do [ -e $d/.done ] || pending+=($d); done
    log "$name: ${#pending[@]} chunks to process"
    [ ${#pending[@]} -gt 0 ] && prefetch ${pending[0]}
    for ((i = 0; i < ${#pending[@]}; i++)); do
        local next=${pending[$((i+1))]:-} pf=
        if [ -n "$next" ]; then prefetch $next & pf=$!; fi
        local t0=$(date +%s)
        run_chunk ${pending[$i]}
        [ -n "$pf" ] && wait $pf
        log "$name: $(basename ${pending[$i]}) done ($((i+1))/${#pending[@]}, $(( $(date +%s) - t0 )) s incl. prefetch of next)"
    done
    merge $run
}

mkdir -p $OUT
exec 9> $OUT/.lock
flock -n 9 || { log "another flow_compare.sh is running on $OUT"; exit 1; }

log "start: ref=$REF_SRC new=$NEW_SRC out=$OUT"
build
pipeline $REF_NAME $REF_SRC "$REF_GLOB" $REF_CHUNK > $OUT/$REF_NAME.log 2>&1 &
p1=$!
pipeline $NEW_NAME $NEW_SRC "$NEW_GLOB" $NEW_CHUNK > $OUT/$NEW_NAME.log 2>&1 &
p2=$!
wait $p1 $p2

mkdir -p $OUT/comparison
rootb "$ANA/Compare_Datasets.C+(\"$OUT/$REF_NAME\",\"$OUT/$NEW_NAME\",\"$OUT/comparison\",\"$REF_LABEL\",\"$NEW_LABEL\")" > $OUT/comparison/compare.log 2>&1
log "comparison written to $OUT/comparison (compare_summary.txt, compare_*.png)"
