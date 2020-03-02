#!/bin/bash
#ref of command awk: https://ubuntuforums.org/showthread.php?t=840054
set -e; bool_exe=true; #bool_exe=false;
# I/O {{{
INPUTDIR="/wk_cms/ykao/tqHGG/2017old/output_preselection"
INPUTDIR_2016="/wk_cms/ykao/tqHGG/2016/input_selection"
INPUTDIR_2017="/wk_cms/ykao/tqHGG/2017/input_selection"
INPUTDIR_2018="/wk_cms/ykao/tqHGG/2018/input_selection"
INPUTDIR_161718="/wk_cms/ykao/tqHGG/161718/input_selection"

OUTPUTDIR="plots/2017old"
OUTPUTDIR_2016="plots/2016"
OUTPUTDIR_2017="plots/2017"
OUTPUTDIR_2018="plots/2018"
OUTPUTDIR_161718="plots/161718"
#}}}
EXECUTABLE=./bin/selection
EXECUTABLE_npu=./bin/selection_npu_float

#--------------- Function ---------------#
# function ExeAnalysis() {{{
function ExeAnalysis() {
    file=$1; year=$2; channel=$3; tag=$4;
    if   [[ ${tag} == 'npu_is_int' ]];   then exe=${EXECUTABLE};     input_directory=${INPUTDIR}; output_directory=${OUTPUTDIR}; message="[MESSAGE] Start to analyze ${file}";
    elif [[ ${tag} == 'npu_is_float' ]]; then exe=${EXECUTABLE_npu}; input_directory=${INPUTDIR}; output_directory=${OUTPUTDIR}; message="[MESSAGE_npu] Start to analyze ${file}";
    elif [[ ${year} == '2016' ]]; then exe=${EXECUTABLE_npu}; input_directory=${INPUTDIR_2016}; output_directory=${OUTPUTDIR_2016}; message="[MESSAGE] Start to anlayze ${year} ${file}";
    elif [[ ${year} == '2017' ]]; then exe=${EXECUTABLE_npu}; input_directory=${INPUTDIR_2017}; output_directory=${OUTPUTDIR_2017}; message="[MESSAGE] Start to anlayze ${year} ${file}";
    elif [[ ${year} == '2018' ]]; then exe=${EXECUTABLE_npu}; input_directory=${INPUTDIR_2018}; output_directory=${OUTPUTDIR_2018}; message="[MESSAGE] Start to anlayze ${year} ${file}";
    elif [[ ${year} == '161718' ]]; then exe=${EXECUTABLE_npu}; input_directory=${INPUTDIR_161718}; output_directory=${OUTPUTDIR_161718}; message="[MESSAGE] Start to anlayze ${year} ${file}";
    else echo "[WARNING] something went wrong in the file: script/exe_selection_batch.sh"; exit 1;
    fi
    #echo "manual check: "${exe}
    #--------------------------------------------------#
    if [ ${bool_exe} = true ]; then
        dataset=`echo $file | awk -F "." '{print $1}'`
        # log messages
        log=log/log_selection_${year}_${dataset}_${channel}.txt
        log_dir=${output_directory}/log
        echo ${message} | tee ${log} 
        # I/O
        input_file=${input_directory}/ntuple_${dataset}.root
        output_sub_dir=${output_directory}/${dataset}
        output_mva_dir=${output_directory}/mva
        output_logscale=${output_directory}/log_scale
        output_file=${output_sub_dir}/hist_${dataset}.root
        output_tree=${output_mva_dir}/tree_${dataset}.root

        mkdir -p ${output_sub_dir}
        mkdir -p ${output_mva_dir}
        mkdir -p ${output_logscale}
        mkdir -p ${log_dir}

        # execution
        ${exe} ${input_file} ${output_file} ${output_tree} ${dataset} ${output_sub_dir} ${channel} | tee -a ${log}
        # backup individual log messages
        cp -p ${log} ${log_dir}
    else
        echo ${message}
        dataset=`echo $file | awk -F "." '{print $1}'`
        input_file=${input_directory}/ntuple_${dataset}.root
        output_sub_dir=${output_directory}/${dataset}
        output_mva_dir=${output_directory}/mva
        output_file=${output_sub_dir}/hist_${dataset}.root
        output_tree=${output_mva_dir}/tree_${dataset}.root

        exe=echo
        ${exe} "test\n" ${input_file} ${output_file} ${output_tree} ${dataset} ${output_sub_dir} ${channel}
        echo "done! $?"
    fi
}
#}}}
# function Get_type_based_on_file_name() {{{
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
FILE=$1; YEAR=$2; CHANNEL=$4;
echo "check year=$YEAR"
if [[ ${YEAR} == "2017old" ]]; then TAG=`Get_type_based_on_file_name ${FILE}`; else TAG=$3; fi
ExeAnalysis ${FILE} ${YEAR} ${CHANNEL} ${TAG}
echo "[MESSAGE] This is the end of $FILE" && exit 0
