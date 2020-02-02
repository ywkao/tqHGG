#!/bin/bash
#ref of command awk: https://ubuntuforums.org/showthread.php?t=840054
set -e

INPUTDIR="/wk_cms2/youying/public/2017_94X_3_1_X_and_3_2_0"
INPUTDIR_ywk="/wk_cms2/ykao/public/2017_94X_3_1_X_and_3_2_0"
OUTPUTDIR="/wk_cms/ykao/tqHGG/2017_oldNtuples/ntuples_skimmed"
EXECUTABLE=./bin/preselection
EXECUTABLE_npu=./bin/preselection_npu_float

# functions{{{
function ExeAnalysis() {
    file=$1; tag=$2;
    if   [[ ${tag} == 'is_regular'           ]]; then exe=${EXECUTABLE}    ; input_directory=${INPUTDIR}    ; root_extension=".root";
    elif [[ ${tag} == 'npu_is_float'         ]]; then exe=${EXECUTABLE_npu}; input_directory=${INPUTDIR_ywk}; root_extension=""     ;
    elif [[ ${tag} == 'is_single_top'        ]]; then exe=${EXECUTABLE}    ; input_directory=${INPUTDIR_ywk}; root_extension=".root";
    elif [[ ${tag} == 'has_multi_root_files' ]]; then exe=${EXECUTABLE}    ; input_directory=${INPUTDIR_ywk}; root_extension=""     ;
    else echo "[WARNING] something went wrong in the file: script/exe_preselection_batch.sh"; exit 1;
    fi
    #echo "manual check: ${exe} ${input_directory} ${root_extension}"
    #echo $file | awk -F "." '{print "'$input_directory'""/"$1"'$root_extension'", "'$OUTPUTDIR'""/ntuple_"$1".root", $1}'
    #--------------------------------------------------
    log=log/log_${file}.txt
    echo "[MESSAGE] Start to analyze ${file}" | tee ${log} 
    echo $file | awk -F "." '{print "'$input_directory'""/"$1"'$root_extension'", "'$OUTPUTDIR'""/ntuple_"$1".root", $1}' |\
    xargs -n3 ${exe} | grep -v ientry | tee -a ${log} # REMARK: input_file / output_file / datasets
    cp -p ${log} ${OUTPUTDIR}/log/
}

function Get_type_based_on_file_name() {
    if [[
          $FILE == 'DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa' ||
          $FILE == 'GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8'
       ]]; then
        local function_result="has_multi_root_files"
        echo ${function_result}
    elif [[
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
    elif [[
            $FILE == 'ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root' ||\
            $FILE == 'ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root' ||\
            $FILE == 'ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root' ||\
            $FILE == 'ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root'
         ]]; then
        local function_result="is_single_top"
        echo ${function_result}
    else
        local function_result="is_regular"
        echo ${function_result}
    fi 
}
#}}}

#--------------- Execution ---------------#
FILE=$1; TAG=`Get_type_based_on_file_name ${FILE}`;
ExeAnalysis ${FILE} ${TAG}
echo "[MESSAGE] This is the end of $FILE" && exit 0
