#!/bin/bash

make -j 8

OUTDIR=output/test2

LOGDIR=logs

mkdir -p ${OUTDIR}
mkdir -p ${LOGDIR}

#most recent MC
INDIR=/nfs-6/userdata/bemarsh/leptoquark/babies/v2

declare -a Samples=(singletop ttdl ttsl dyjetsll wjets ttz ttw ttg tth ww_2l2nu)
# declare -a Samples=(T2tt_RPV_700 T2tt_RPV_900 T2tt_RPV_1100)

for SAMPLE in ${Samples[@]};
  do echo ./runLooper ${INDIR} ${SAMPLE} ${OUTDIR}
  nohup ./runLooper ${INDIR} ${SAMPLE} ${OUTDIR} >& ${LOGDIR}/log_${SAMPLE}.txt &
  # ./runLooper ${INDIR} ${SAMPLE} ${OUTDIR} #>& ${LOGDIR}/log_${SAMPLE}.txt &
done
