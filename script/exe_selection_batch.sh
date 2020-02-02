#!/bin/bash
#ref of command awk: https://ubuntuforums.org/showthread.php?t=840054
set -e

#INPUTDIR="/wk_cms/ykao/tqHGG/2017_oldNtuples/ntuples_skimmed"
INPUTDIR="/wk_cms/ykao/tqHGG/2017_oldNtuples/ntuples_skimmed_mgg_tighter"
OUTPUTDIR="plots"
EXECUTABLE=./bin/selection
EXECUTABLE_npu=./bin/selection_npu_float

# functions{{{
#--------------- Function ---------------#
function ExeAnalysis() {
    file=$1; channel=$2; tag=$3;
    if   [[ ${tag} == 'npu_is_int' ]];   then exe=${EXECUTABLE};     message="[MESSAGE] Start to analyze ${file}";
    elif [[ ${tag} == 'npu_is_float' ]]; then exe=${EXECUTABLE_npu}; message="[MESSAGE_npu] Start to analyze ${file}";
    else echo "[WARNING] something went wrong in the file: script/exe_selection_batch.sh"; exit 1;
    fi
    #echo "manual check: "${exe}
    #--------------------------------------------------#
    log=log/log_skimmed_ntuple_${file}_${channel}.txt
    echo ${message} | tee ${log} 
    echo $file | awk -F "." '{print "'$INPUTDIR'""/ntuple_"$1".root", "'$OUTPUTDIR'""/"$1"/hist_"$1".root", "'$OUTPUTDIR'""/mva/tree_"$1".root", $1, "'$OUTPUTDIR'""/"$1, "'$channel'"}' |\
    xargs -n6 ${exe} | tee -a ${log} # REMARK: input_file / output_file / output_tree / datasets / output_dir / channel
    cp -p ${log} ${OUTPUTDIR}/log
}

function Get_type_based_on_file_name() {
    FILE=$1
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
        local function_result="npu_is_float"
        echo ${function_result}
    else
        local function_result="npu_is_int"
        echo ${function_result}
    fi
}
#}}}

#--------------- Execution ---------------#
FILE=$1; CHANNEL=$2; TAG=`Get_type_based_on_file_name ${FILE}`
ExeAnalysis ${FILE} ${CHANNEL} ${TAG}
echo "[MESSAGE] This is the end of $FILE" && exit 0
