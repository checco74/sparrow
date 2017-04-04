#! /usr/bin/tcsh

set ANN = /scratch/scratch/rasinski/ann.prg
set WDIR = /public/scratch/rasinski/
set PROG = ./ann

set DIRNAME = ann
if ("$1" != "") then
  set DIRNAME = $1
endif


echo "copying files..."
mkdir -p ${WDIR}/${DIRNAME}
cd ${WDIR}/${DIRNAME}
echo "building ANN..."
make fullclean
make ann
echo "starting program..."
${PROG} ${1}

mkdir -p ${OPT}/${DIRNAME}
cp *.log ${OPT}/${DIRNAME}
cp *.csv ${OPT}/${DIRNAME}
cp *.prm ${OPT}/${DIRNAME}
cp ann.* ${OPT}/${DIRNAME}
cp *.raw ${OPT}/${DIRNAME}
cp ann.cfg ${OPT}/${DIRNAME}
cd ${WDIR}
rm -rf ${DIRNAME}
echo "script ended..."
