#!/bin/bash
#ref of command awk: https://ubuntuforums.org/showthread.php?t=840054
set -e

#INPUTDIR="ntuples_skimmed"
#INPUTDIR="/wk_cms/ykao/tqHGG/ntuples_skimmed_trigger_pustudy"
INPUTDIR="/wk_cms/ykao/tqHGG/ntuples_skimmed"
OUTPUTDIR="plots"
EXECUTABLE=./bin/selection
EXECUTABLE_npu=./bin/selection_npu_float

#--------------- Function ---------------#
function ExeAnalysis(){
    file=$1
    channel=$2
    log=log/log_skimmed_ntuple_${file}.txt
    echo "[MESSAGE] Start to analyze ${file}" | tee ${log} 
    echo $file | awk -F "." '{print "'$INPUTDIR'""/ntuple_"$1".root", "'$OUTPUTDIR'""/"$1"/hist_"$1".root", $1, "'$OUTPUTDIR'""/"$1, "'$channel'"}' |\
    xargs -n5 ${EXECUTABLE} | grep -v ientry | tee -a ${log}
    # REMARK: input_file / output_file / datasets / output_dir
}

function ExeAnalysis_npu_float(){
    file=$1
    channel=$2
    log=log/log_skimmed_ntuple_${file}.txt
    echo "[MESSAGE_npu] Start to analyze ${file}" | tee ${log} 
    echo $file | awk -F "." '{print "'$INPUTDIR'""/ntuple_"$1".root", "'$OUTPUTDIR'""/"$1"/hist_"$1".root", $1, "'$OUTPUTDIR'""/"$1, "'$channel'"}' |\
    xargs -n5 ${EXECUTABLE_npu} | grep -v ientry | tee -a ${log}
    # REMARK: input_file / output_file / datasets / output_dir
}
#--------------- DryRun ---------------#
if [[ $1 == '-d' || $1 == '--dryRun' ]]; then
    #ExeAnalysis ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root leptonic 
    ExeAnalysis ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root hadronic 
    echo "[MESSAGE] This is the end of dryRun execution." && exit 0
fi

#--------------- Execution ---------------#
FILE=$1; CHANNEL=$2
if [[
      $FILE == 'TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8' ||
      $FILE == 'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8' ||
      $FILE == 'WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8' ||
      $FILE == 'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8' ||
      $FILE == 'WW_TuneCP5_13TeV-pythia8' ||
      $FILE == 'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8' ||
      $FILE == 'ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8' ||
      $FILE == 'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'
   ]]; then
    ExeAnalysis_npu_float $FILE $CHANNEL;
    echo "[MESSAGE] This is the end of $FILE" && exit 0
else
    ExeAnalysis $FILE $CHANNEL;
    echo "[MESSAGE] This is the end of $FILE" && exit 0
fi
