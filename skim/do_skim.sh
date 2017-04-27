#! /bin/bash

LOGDIR=logs/

# INPATH=/nfs-6/userdata/mt2/V00-08-14_moriondmc
# INPATH=/nfs-6/userdata/mt2/V00-08-12_json_271036-284044_23Sep2016ReReco_36p26fb
INPATH=/nfs-6/userdata/mt2/T2tt_RPV
OUTPATH=/nfs-6/userdata/bemarsh/leptoquark/babies/v2/

# declare -a Samples=(wjets ttsl ttdl singletop dyjetsll ttw ttz ttg tth ww_2l2nu)
# declare -a Samples=(data_Run2016B_SingleMuon data_Run2016B_DoubleMuon data_Run2016C_SingleMuon data_Run2016C_DoubleMuon data_Run2016D_SingleMuon data_Run2016D_DoubleMuon data_Run2016E_SingleMuon data_Run2016E_DoubleMuon data_Run2016F_SingleMuon data_Run2016F_DoubleMuon data_Run2016G_SingleMuon data_Run2016G_DoubleMuon data_Run2016H_SingleMuon data_Run2016H_DoubleMuon )
declare -a Samples=(T2tt_RPV_700 T2tt_RPV_900 T2tt_RPV_1100)

mkdir -p $OUTPATH
mkdir -p $LOGDIR

for SAMPLE in ${Samples[@]}

do nohup root -b -q -l skim_mt2_babies.C\(\"$INPATH\",\"$OUTPATH\",\"$SAMPLE\"\) >& $LOGDIR/log_${SAMPLE}.txt &

done