#!/bin/bash
set -e; bool_exe=true; #bool_exe=false;
dir="/wk_cms/ykao/tqHGG"
tag=""
#tag="_noPhotonPT"
#tag="_test"
#tag="_tighterPhotonPT"
output_preselection="output_preselection${tag}"
input_selection="input_selection${tag}"

# function combine_dataset(){{{
function combine_dataset(){
    year="$1"
    source_dir="${dir}/${year}/${output_preselection}"
    target_dir="${dir}/${year}/${input_selection}"; mkdir -p ${target_dir};
    if [ $bool_exe = "true" ];
    then
        exe=./script/hadd_${year}.sh
    else
        exe=./script/test_${year}.sh
    fi
    ${exe} ${target_dir} ${source_dir}
}
#}}}
# function mergeThreeYears(){{{
function mergeThreeYears(){
    rootfile=$1
    if [ $bool_exe = "true" ];
    then
        hadd ${target_dir}/${rootfile} \
             ${source_2016}/${rootfile} \
             ${source_2017}/${rootfile} \
             ${source_2018}/${rootfile}
    else
        echo ${target_dir}/${rootfile}
        echo ${source_2016}/${rootfile}
        echo ${source_2017}/${rootfile}
        echo ${source_2018}/${rootfile}
    fi
}
#}}}
# function mergeTwoYears(){{{
function mergeTwoYears(){
    # 2016 lacks VH samples......
    rootfile=$1
    if [ $bool_exe = "true" ];
    then
        hadd ${target_dir}/${rootfile} \
             ${source_2017}/${rootfile} \
             ${source_2018}/${rootfile}
    else
        echo ${target_dir}/${rootfile}
        echo ${source_2017}/${rootfile}
        echo ${source_2018}/${rootfile}
    fi
}
#}}}

#------------------------ Execution Section --------------------------#

# combine similar processes
combine_dataset 2016
combine_dataset 2017
combine_dataset 2018

# merge three years together
source_2016="${dir}/2016/${input_selection}"
source_2017="${dir}/2017/${input_selection}"
source_2018="${dir}/2018/${input_selection}"
target_dir="${dir}/161718/${input_selection}"; mkdir -p ${target_dir};

mergeThreeYears ntuple_Data.root
mergeThreeYears ntuple_DYJetsToLL.root
mergeThreeYears ntuple_DiPhotonJetsBox.root
mergeThreeYears ntuple_GJet.root
mergeThreeYears ntuple_GluGluHToGG.root
mergeThreeYears ntuple_QCD.root
mergeThreeYears ntuple_TGJets.root
mergeThreeYears ntuple_TTGG.root
mergeThreeYears ntuple_TTGJets.root
mergeThreeYears ntuple_TTJets.root
mergeThreeYears ntuple_VBFHToGG.root
mergeThreeYears ntuple_ttHJetToGG.root
mergeThreeYears ntuple_WGToLNuG.root
mergeThreeYears ntuple_WW.root
mergeThreeYears ntuple_WZTo2L2Q.root
mergeThreeYears ntuple_ZGToLLG.root
mergeThreeYears ntuple_ZZTo2L2Q.root

mergeTwoYears ntuple_VHToGG.root

### use 2017 signal samples as representatives 
### WARNING: should be done after src/selection.cpp because of int/float issue!
##source_2017old="${dir}/2017old/${output_preselection}"
##for fcnc_signal in `ls ${source_2017old} | grep -i fcnc`
##do
##    file=${source_2017old}/${fcnc_signal}
##
##done
##
##ls ${source_2016}
##ls ${source_2017}
##ls ${source_2018}
##ls ${target_dir}

#--------------------------------------------------#
## From 2017 old ntuples (DiPhotonJetBox_M40_80 is given up)
#path_2017old="${dir}/keep_temporarily_2017old"
#source_2017old="${dir}/2017old/${input_selection}"
#hadd ${path_161718}/ntuple_DiPhotonJetsBox_M40_80-Sherpa.root ${path_2017old}/ntuple_DiPhotonJetsBox_M40_80-Sherpa.root
#hadd ${target_dir}/ntuple_WJetsToLNu.root ${path_2017old}/ntuple_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root
