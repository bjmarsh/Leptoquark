#!/bin/bash

make -j 8

OUTDIR=output/test

LOGDIR=logs

mkdir -p ${OUTDIR}
mkdir -p ${LOGDIR}

#most recent MC
INDIR=/nfs-6/userdata/bemarsh/leptoquark/babies/v2

# declare -a Samples=(singletop ttdl ttsl dyjetsll wjets ttz ttw ttg tth ww_2l2nu)
# declare -a Samples=(data_Run2016B data_Run2016C data_Run2016D data_Run2016E data_Run2016F data_Run2016G data_Run2016H)
declare -a Samples=(T2tt_RPV_700 T2tt_RPV_900 T2tt_RPV_1100)
# declare -a Samples=(singletop)

for SAMPLE in ${Samples[@]};
  do echo ./runLooper ${INDIR} ${SAMPLE} ${OUTDIR}
  nohup ./runLooper ${INDIR} ${SAMPLE} ${OUTDIR} >& ${LOGDIR}/log_${SAMPLE}.txt &
done
