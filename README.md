# EPOS N-job production daemon — setup guide

How to reproduce, on a fresh machine, the automated EPOS 4.0.3 Au+Au production
that fills `auau200-7/`. Written 2026-08-14 from the running setup on this box.

The daemon keeps **N independent EPOS instances** generating events forever,
harvesting each finished ROOT file into one dataset directory with a
monotonically increasing index. On this machine N = 5 and the dataset is at
682 files / 151 GB.

---

## 1. What the daemon actually does

```
epos_scheduler.sh   (the daemon — one process, sleeps 7 h between cycles)
   └── cycle_epos.sh   (one cycle: for each instance n = 1..N)
         ├── is `epos -root auau_run_n` still alive?  →  skip this instance
         ├── move  $JIN/epos<n>/z-auau_run_<n>.root
                →  $JIN/<DEST>/z-auau_run_<i>.root     (i = next free index)
         └── start the next run for instance n, detached
```

Key design points, all deliberate:

- **Fixed working filename, indexed archive name.** Every instance always writes
  the same file (`z-auau_run_<n>.root`); the indexing happens at harvest time, so
  the generator never needs to know how many files came before.
- **`pgrep` guard.** A ROOT file is never touched while its producer is alive.
  This is the entire safety mechanism — do not remove it.
- **Per-instance subshell.** Each instance's `startup.sh` exports `HTO`/`CHK`
  pointing at *its own* output directory. Sourcing them in one shell would make
  every instance write to the same place, so `cycle_epos.sh` sources each inside
  `( … )`. **The epos online documentation explains what these environment variables are, some can be left as the epos installation directory**
- **Polling, not event-driven.** A run takes ~10 h; the cycle interval is 7 h.
  An instance that finishes just after a tick sits idle until the next one, so
  expect some idle capacity. Shortening the interval costs nothing but log noise.

---

## 2. Sizing: how to choose N

Measured on the running instances:

| resource | per EPOS instance |
|---|---|
| CPU | ~0.8 of one core (single-threaded) |
| **RAM** | **~17 GB resident** | (ps, i think this number could actually be a bit higher for some jobs)
| wall clock | ~10 h for 1000 events (`nfull 50` × `nfreeze 20`) |
| disk | ~240 MB per output ROOT file |
| disk (install) | ~4.4 GB per instance (each is a full private copy) |

**RAM is the binding constraint, not cores.** This box has 24 cores and 94 GB,
so N = 5 (≈85 GB) is already at the edge — cores are not the limit, memory is.

```
N  =  floor( (usable_RAM_GB - 8) / 17 )
```

Leave headroom: analysis jobs (`root -l -b -q 'Elliptic.C+'`) run on the same
machine and want several GB. Also budget install space: N × 4.4 GB on the
instance disk, plus ~240 MB × (runs you intend to keep) on the data disk.

---

## 3. Dependencies

Install these in order — EPOS's CMake requires all of them at configure time
when built with `-DCOMPILE_OPTION=BASIC` (which turns on ROOT, FastJet, HepMC3,
hydro and XA).

### 3.1 System packages

(note this is likely already installed on your machine if it has ROOT installed)

```bash
sudo apt install build-essential gfortran cmake git zlib1g-dev \
                 python3-dev libx11-dev libxpm-dev libxft-dev libxext-dev
```

`gfortran` is not optional — EPOS is `project(EPOS LANGUAGES Fortran C CXX)` and
the bulk of the physics is Fortran.

### 3.2 ROOT ≥ 6.16 (this machine: 6.30.06)

Either install a binary release or build from source. EPOS needs the `Tree`,
`Hist` and `Physics` components and a CMake config package (`find_package(ROOT
… CONFIG)`), so a distro package without `ROOTConfig.cmake` will not work.

ROOT must be in the environment **before** anything else is configured, because
`startup.sh` interpolates `$ROOTSYS` into `LD_LIBRARY_PATH`:

```bash
source /path/to/root/bin/thisroot.sh
```

On this machine ROOT is a shared group install and `/etc/profile` line 29 sources
`/home/zji/programs/root/install/bin/thisroot.sh` for every login shell. On your
own machine, put the `source` line in `~/.bashrc` — and see §8 for why a daemon
started from a shell without it will fail.

### 3.3 FastJet 3.5.1

```bash
tar xzf fastjet-3.5.1.tar.gz && cd fastjet-3.5.1
./configure --prefix=$HOME/fastjet-install
make -j8 && make install
```

`startup.sh` locates it via `fastjet-config --prefix`, so the binary at
`$HOME/fastjet-install/bin/fastjet-config` must exist and be executable.

### 3.4 HepMC3 3.2.6

The exact configure line used here (from `~/compile.txt`) — ROOT I/O and
protobuf off, Python on:

```bash
mkdir hepmc3-build && cd hepmc3-build
cmake -DHEPMC3_ENABLE_ROOTIO=OFF \
      -DHEPMC3_ENABLE_PROTOBUFIO:BOOL=OFF \
      -DHEPMC3_ENABLE_PYTHON:BOOL=ON \
      -DHEPMC3_PYTHON_VERSIONS=3.10 \
      -DHEPMC3_PYTHON_SITEARCH310=../hepmc3-install/lib/python3.10/site-packages \
      -DCMAKE_INSTALL_PREFIX=../hepmc3-install \
      ../HepMC3-3.2.6
make -j8 && make install
```

Adjust `3.10` to your system Python. EPOS finds it through `HepMC3_DIR`, which
`startup.sh` sets to `<prefix>/share/HepMC3/cmake`.

### 3.5 EPOS 4.0.3 source

Not publicly downloadable — request the tarball from the EPOS4 project page
(<https://klaus.pages.in2p3.fr/epos4/>) and cite the papers listed under "Cite"
when you publish. The archive used here is `~/epos4.0.3.tar` (1.9 GB; most of
that is the `z-eos4f.eos` equation-of-state table, 419 MB, which the hydro needs
at runtime and must sit next to the run card).

---

## 4. Build one instance

Every instance is a **complete, independent copy** of the source tree and its
build — they are not sharing a binary. That is what lets each one keep its own
`auau200/` working directory and write a fixed-name output file without
collisions.

```bash
mkdir ~/epos1 && cd ~/epos1
tar xf ~/epos4.0.3.tar          # → ~/epos1/epos4.0.3/
```

Create `~/epos1/startup.sh` (see §5 for the real, symlinked version) and build:

```bash
cd ~/epos1
source startup.sh
cmake -S $EPO -B $BUILD_DIR -DCMAKE_INSTALL_PREFIX=$BIN_DIR \
      -DCOMPILE_OPTION=BASIC -DCMAKE_BUILD_TYPE=Release \
      -DFASTSYS=$FASTJET_DIR -DCMAKE_INSTALL_MESSAGE=LAZY
cmake --build $BUILD_DIR -j8
cmake --install $BUILD_DIR
```

This is verbatim `~/epos1/compile.sh`. It produces `$EPO/bin/epos` (a shell
wrapper) and `$EPO/bin/Xepos` (the actual executable that shows up in `ps` and
owns the 17 GB).

Smoke-test before replicating:

```bash
cd ~/epos1 && source startup.sh && cd epos4.0.3/auau200
$BIN_DIR/bin/epos -root auau_run_1
```

Fix any missing-library errors here. Debugging five broken copies is five times
the work.

### VERY VERY VERY IMPORTANT ###

to obtain random events, you actually have to go to `epos1/epos4.0.3/bin/epos` and modify somewhere on lines 209 and 210 seedi and seedj to be `date  '+%N'`. For example , right now on the current installation

`
seedj=`date  '+%N'`
seedi=`date  '+%N'`
`


---

## 5. Environment: `startup.sh` per instance

Each instance gets its own, differing **only in the instance number**. Instance
1's (`scripts/epos1/startup.sh`):

```bash
export EPOVSN=4.0.3
export MYDIR=$PWD/
export EPO=${MYDIR}epos${EPOVSN}/
export BUILD_DIR=${MYDIR}epos-build
export BIN_DIR=${EPO}
export JIN=/media/Students/aliu/
export OPT=./
export HTO=${JIN}epos1/          # ← instance number here
export CHK=${JIN}epos1/          # ← and here
export HepMC3_DIR=/home/aliu/hepmc3-install/share/HepMC3/cmake
export FASTJETSYS=`/home/aliu/fastjet-install/bin/fastjet-config --prefix`
export FASTJET_DIR=${FASTJETSYS}
export PATH=$PATH:$FASTSYS/bin
export LD_LIBRARY_PATH=$ROOTSYS/lib:$ROOTSYS/lib/root:/usr/local/lib:/home/aliu/hepmc3-install/lib:/home/aliu/fastjet-install/lib:$LD_LIBRARY_PATH
```

Three things to change for your machine: `JIN` (the data disk), and the two
absolute paths to your HepMC3 and FastJet installs.

Two things that must **not** change:

- `MYDIR=$PWD/` — this file is only correct when sourced from inside
  `~/epos<n>/`. `cycle_epos.sh` `cd`s there first; do the same by hand.
- `HTO`/`CHK` pointing at a **per-instance** directory. If two instances share
  one, they overwrite each other's output mid-run.

Generate the whole set:

```bash
N=5
SCRIPTS=/media/Students/epos-scripts    # where the scripts live (data disk)
JIN=/media/Students/aliu                # where output goes
for n in $(seq 1 $N); do
  mkdir -p "$SCRIPTS/epos$n" "$JIN/epos$n"
  sed "s/epos1\//epos$n\//g" "$SCRIPTS/epos1/startup.sh" > "$SCRIPTS/epos$n/startup.sh"
  chmod +x "$SCRIPTS/epos$n/startup.sh"
  ln -sf "$SCRIPTS/epos$n/startup.sh" "$HOME/epos$n/startup.sh"
done
```

### Why the symlinks

The scripts physically live in the versioned `scripts/` directory on the data
disk; everything under `$HOME` is a symlink into it:

```
~/epos_scheduler.sh   → <project>/scripts/epos_scheduler.sh
~/cycle_epos.sh       → <project>/scripts/cycle_epos.sh
~/run_all.sh          → <project>/scripts/run_all.sh
~/run.sh              → <project>/scripts/run.sh
~/startup.sh          → <project>/scripts/startup.sh
~/epos<n>/startup.sh  → <project>/scripts/epos<n>/startup.sh
```

Edit either path — same file, and edits are captured by git. The cost is that
the daemon **cannot start until the data disk is mounted** (§8).

---

## 6. The run card

Each instance needs a config file. `~/epos<n>/epos4.0.3/auau200/auau_run_<n>.optns` is the config used here (30–40 % central Au+Au
at 200 GeV, 1000 events per file):

```
application hadron
set laproj 79       ! projectile Z          (Au)
set maproj 197      ! projectile A
set latarg 79       ! target Z              (Au)
set matarg 197      ! target A
set ecms 200        ! sqrt(s)_NN [GeV]

set istmax 25
set iranphi 1       ! random reaction-plane rotation per event  ← branch `phi`
ftime on

nodecays
 120 -120 130 -130  ! keep π±, K± stable (EPOS ids, not PDG)
end

set ninicon 5       ! initial conditions per hydro evolution
core full           ! core/corona
hydro hlle
eos x3ff
hacas full          ! hadronic cascade
set nfull 50        ! hydro events
set nfreeze 20      ! freeze-out events per hydro event → 50×20 = 1000 events
set centrality 0    ! 0 = min bias, then cut on b:
set bminim 8.4      ! impact parameter window [fm] = 30–40 % centrality
set bmaxim 9.2

fillTree4(C1)       ! ROOT output → z-<basename>.root
```

`iranphi 1` is what makes the `phi` branch a genuine random reaction-plane angle;
the flow analysis depends on it. `nfull`/`nfreeze` set the ~10 h run length —
halve them if you want a faster harvest cadence at the cost of more, smaller
files. 

Copy it per instance, renaming to match the `-root` argument:

```bash
for n in $(seq 2 $N); do
  sed "s/auau_run_1/auau_run_$n/g" \
      ~/epos1/epos4.0.3/auau200/auau_run_1.optns \
    > ~/epos$n/epos4.0.3/auau200/auau_run_$n.optns
done
```

(The `z-auau_run_<n>.optns` and `.clinput` files that appear next to it are
per-run dumps written by EPOS, not inputs. Don't edit them.)

---

## 7. Adapting the scripts

### Paths hardcoded in the scripts

These are absolute and written for this machine. Change all four before the
first run:

| file | line | change to |
|---|---|---|
| `cycle_epos.sh` | `JINBASE=/media/Students/aliu` | your data disk |
| `cycle_epos.sh` | `HOMEBASE=/home/aliu` | wherever `epos<n>/` trees live |
| `epos_scheduler.sh` | `/home/aliu/cycle_epos.sh "$DEST"` | your path to `cycle_epos.sh` |
| `run_all.sh` | `JINBASE=/media/Disk_Jin/aliu` | your data disk (**this one is stale even here** — see §8) |

`HOMEBASE` and `JINBASE` are intentionally separate: the instance trees live on
the fast local NVMe, the output goes to the big shared disk.

### Making it an N-job daemon

Both scripts currently hardcode the instance list:

```bash
for n in 1 2 3 4 5; do
```

For a different N, change that line in `cycle_epos.sh` (and `run_all.sh` if you
use it) to read from the environment:

```bash
NJOBS=${NJOBS:-5}
for n in $(seq 1 "$NJOBS"); do
```

then `NJOBS=8 ./epos_scheduler.sh auau200-7` — but note `epos_scheduler.sh`
would need to export it, so add `export NJOBS` near the top. Everything else in
both scripts is already N-agnostic; the loop bound is the only assumption.

Nothing anywhere assumes the instances are *contiguous* either — if instance 3's
build is broken, dropping it from the list is safe, and the harvest indices stay
gapless because they're allocated at move time.

---

## 8. Running it

```bash
cd ~
nohup ./epos_scheduler.sh auau200-7 > scheduler.log 2>&1 &
```

Arguments: `epos_scheduler.sh [dataset-dir] [interval]`, defaults `auau200-7`
and `7h`. The dataset dir is created under `$JINBASE` if missing.

```bash
tail -f ~/scheduler.log                  # watch the daemon
tail -f $JIN/epos1/output.log            # watch one instance
ps -eo pid,etime,rss,cmd | grep Xepos    # what's actually running
pkill -f epos_scheduler.sh               # stop the daemon
```

Stopping the daemon **does not stop running simulations** — they're detached
`nohup` children and will finish their current run, then sit idle. That is the
intended way to drain the machine: kill the daemon, wait for the last runs, then
`RESTART=0 ./cycle_epos.sh auau200-7` to harvest what they produced without
starting anything new.

**Three ways this bites you:**

1. **ROOT environment.** The daemon inherits its environment from the shell that
   launched it. `startup.sh` uses `$ROOTSYS` but never sets it — start the daemon
   from a login shell where `thisroot.sh` has been sourced, or every run fails at
   link time.

2. **Mount order after a reboot.** The data disk must be mounted *before* the
   daemon starts, or every `~/…startup.sh` symlink dangles and `source startup.sh`
   silently does nothing. Check first:

   ```bash
   readlink -e ~/cycle_epos.sh    # empty output = not mounted, do not start
   ```

3. **Harvest lag.** `cycle_epos.sh` only acts on ticks. A run that ends one
   minute after a cycle waits the full interval before its file moves and its
   successor starts. Normal — an idle instance in `ps` is not a failure.

### `run_all.sh` — the other mode

Batch, not daemon: runs each instance's runs back-to-back `NBATCH` times, in
parallel, then exits. `./run_all.sh 3 auau200-7` → 15 files. It restarts each
instance the instant its run ends (no polling lag) and uses `flock` to allocate
indices, so it's strictly better for a *finite* production. The daemon is for
open-ended running.

**Known bug, unfixed:** `run_all.sh` still points `JINBASE` at the dead
`/media/Disk_Jin/aliu` path from before the 2026-07-02 move. Fix that line before
using it. `cycle_epos.sh` — what the daemon actually calls — is correct.

---

## 9. File map

| path | role |
|---|---|
| `scripts/epos_scheduler.sh` | the daemon: `cycle_epos.sh` in a `while true; sleep` loop |
| `scripts/cycle_epos.sh` | one harvest+restart pass over all instances |
| `scripts/run_all.sh` | finite batch alternative (see bug above) |
| `scripts/startup.sh` | generic env template, no `JIN`/`HTO`/`CHK` — not used by the daemon |
| `scripts/epos<n>/startup.sh` | the real per-instance environment |
| `scripts/run.sh` | one-liner to launch instance 1 by hand |
| `~/epos<n>/` | private source tree + build, 4.4 GB each |
| `~/epos<n>/compile.sh` | the CMake invocation from §4 |
| `~/epos<n>/epos4.0.3/auau200/` | run card, `z-eos4f.eos`, per-run dumps |
| `$JIN/epos<n>/` | live output: `output.log`, `z-auau_run_<n>.root` in progress |
| `$JIN/auau200-7/` | harvested dataset, `z-auau_run_<i>.root` |

For what's inside the ROOT files — and the several misleadingly-named branches
you must not take at face value — see `CLAUDE.md` at the project root.
