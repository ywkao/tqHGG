#!/bin/bash
#ref of command awk: https://ubuntuforums.org/showthread.php?t=840054
set -e; bool_exe=true; #bool_exe=false;
# I/O{{{
INPUTDIR="/wk_cms2/youying/public/2017_94X_3_1_X_and_3_2_0"
INPUTDIR_ywk="/wk_cms2/ykao/public/2017_94X_3_1_X_and_3_2_0"
INPUTDIR_2016="/wk_cms2/ykao/public/thFCNC/flashgg_105X/2016"
INPUTDIR_2017="/wk_cms2/youying/public/thFCNC/flashgg_105X/2017"
INPUTDIR_2018="/wk_cms2/youying/public/thFCNC/flashgg_105X/2018"

BIN=bin
OUTPUTDIR="/wk_cms/ykao/tqHGG/2017old/output_preselection"
OUTPUTDIR_2016="/wk_cms/ykao/tqHGG/2016/output_preselection"
OUTPUTDIR_2017="/wk_cms/ykao/tqHGG/2017/output_preselection"
OUTPUTDIR_2018="/wk_cms/ykao/tqHGG/2018/output_preselection"

#BIN=bin_2
#OUTPUTDIR_2016="/wk_cms/ykao/tqHGG/2016/output_preselection_tighterPhotonPT"
#OUTPUTDIR_2017="/wk_cms/ykao/tqHGG/2017/output_preselection_tighterPhotonPT"
#OUTPUTDIR_2018="/wk_cms/ykao/tqHGG/2018/output_preselection_tighterPhotonPT"

#BIN=bin_3
#OUTPUTDIR_2016="/wk_cms/ykao/tqHGG/2016/output_preselection_noPhotonPT"
#OUTPUTDIR_2017="/wk_cms/ykao/tqHGG/2017/output_preselection_noPhotonPT"
#OUTPUTDIR_2018="/wk_cms/ykao/tqHGG/2018/output_preselection_noPhotonPT"

#OUTPUTDIR="/wk_cms/ykao/tqHGG/2017_oldNtuples/ntuples_skimmed"
#}}}
EXECUTABLE=./${BIN}/preselection
EXECUTABLE_npu=./${BIN}/preselection_npu_float

# function KeepTrack(){{{
function KeepTrack(){
    log_dir=$1; dir=$2; code=$3;
    source_code=./${dir}/${code};
    target=${log_dir}/${code};
    if [ ! -f ${target} ]; then echo "cp -p ${source_code} ${target}" ; cp -p ${source_code} ${target}; fi
}
#}}}
# function ExeAnalysis() {{{
function ExeAnalysis() {
    file=$1; year=$2; tag=$3;
    if   [[ ${tag} == 'is_regular'                           ]]; then exe=${EXECUTABLE}    ; input_directory=${INPUTDIR}     ; output_directory=${OUTPUTDIR}     ; root_extension=".root";
    elif [[ ${tag} == 'npu_is_float'                         ]]; then exe=${EXECUTABLE_npu}; input_directory=${INPUTDIR_ywk} ; output_directory=${OUTPUTDIR}     ; root_extension=""     ;
    elif [[ ${tag} == 'is_single_top'                        ]]; then exe=${EXECUTABLE}    ; input_directory=${INPUTDIR_ywk} ; output_directory=${OUTPUTDIR}     ; root_extension=".root";
    elif [[ ${tag} == 'has_multi_root_files'                 ]]; then exe=${EXECUTABLE}    ; input_directory=${INPUTDIR_ywk} ; output_directory=${OUTPUTDIR}     ; root_extension=""     ;
    elif [[ ${year} == '2016' ]] && [[ ${tag} == 'rootfile'  ]]; then exe=${EXECUTABLE_npu}; input_directory=${INPUTDIR_2016}; output_directory=${OUTPUTDIR_2016}; root_extension=".root";
    elif [[ ${year} == '2016' ]] && [[ ${tag} == 'directory' ]]; then exe=${EXECUTABLE_npu}; input_directory=${INPUTDIR_2016}; output_directory=${OUTPUTDIR_2016}; root_extension=""     ;
    elif [[ ${year} == '2017' ]] && [[ ${tag} == 'rootfile'  ]]; then exe=${EXECUTABLE_npu}; input_directory=${INPUTDIR_2017}; output_directory=${OUTPUTDIR_2017}; root_extension=".root";
    elif [[ ${year} == '2017' ]] && [[ ${tag} == 'directory' ]]; then exe=${EXECUTABLE_npu}; input_directory=${INPUTDIR_2017}; output_directory=${OUTPUTDIR_2017}; root_extension=""     ;
    elif [[ ${year} == '2018' ]] && [[ ${tag} == 'rootfile'  ]]; then exe=${EXECUTABLE_npu}; input_directory=${INPUTDIR_2018}; output_directory=${OUTPUTDIR_2018}; root_extension=".root";
    elif [[ ${year} == '2018' ]] && [[ ${tag} == 'directory' ]]; then exe=${EXECUTABLE_npu}; input_directory=${INPUTDIR_2018}; output_directory=${OUTPUTDIR_2018}; root_extension=""     ;
    else echo "[WARNING] something went wrong in the file: script/exe_preselection_batch.sh"; exit 1;
    fi
    mkdir -p ${output_directory} # create output directory if it is new
    if [[ ${root_extension} == '' ]]; then isDirecotry="directory"; else isDirecotry="rootfile"; fi; message="[MESSAGE] Start to analyze ${file} [check] file type: ${isDirecotry}"
    #--------------------------------------------------
    if [ ${bool_exe} = true ]; then
        dataset=`echo $file | awk -F "." '{print $1}'`
        log=log/log_preselection_${year}_${dataset}.txt
        input_file="${input_directory}/${dataset}${root_extension}"
        output_file="${output_directory}/ntuple_${dataset}.root"
        echo "${message}" | tee ${log} 
        ${exe} ${input_file} ${output_file} ${dataset} ${isDirecotry} ${year} | tee -a ${log}
        # store log messages
        log_dir="${output_directory}/log"; mkdir -p ${log_dir}; cp -p ${log} ${log_dir};
        KeepTrack ${log_dir} "src" "preselection.cpp"
        KeepTrack ${log_dir} "include" "preselection_criteria.h"
    else
        dataset=`echo $file | awk -F "." '{print $1}'`
        input_file="${input_directory}/${dataset}${root_extension}"
        output_file="${output_directory}/ntuple_${dataset}.root"
        exe=echo
        ${exe} ${message}
        ${exe} ${input_file}
        ${exe} ${output_file}
        #${exe} ${input_file} ${output_file} ${dataset} ${isDirecotry} ${year}

        log_dir="${output_directory}/log"
        KeepTrack ${log_dir} "src" "preselection.cpp"
        KeepTrack ${log_dir} "include" "preselection_criteria.h"
    fi
}
#}}}
#function Get_type_based_on_file_name() {{{
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
FILE=$1; YEAR=$2;
echo "check year=$YEAR"
if [[ ${YEAR} == "2017old" ]]; then TAG=`Get_type_based_on_file_name ${FILE}`; else TAG=$3; fi
ExeAnalysis ${FILE} ${YEAR} ${TAG}
echo "[MESSAGE] This is the end of $FILE" && exit 0
