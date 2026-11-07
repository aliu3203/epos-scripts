export EPOVSN=4.0.3
export MYDIR=$PWD/
export EPO=${MYDIR}epos${EPOVSN}/
export BUILD_DIR=${MYDIR}epos-build
export BIN_DIR=${EPO}
export OPT=./
export HTO=./
export CHK=./
export FASTJETSYS=`/home/aliu/fastjet-install/bin/fastjet-config --prefix`
export FASTJET_DIR=${FASTJETSYS}
export PATH=$PATH:$FASTSYS/bin
export LD_LIBRARY_PATH=$ROOTSYS/lib:$ROOTSYS/lib/root:/usr/local/lib:$LD_LIBRARY_PATH


# THIS SCRIPT IS OUTDATED, startup.sh IS LOCATED UNDER INDIVIDUAL EPOS FOLDERS (epos1, epos2, ...)
