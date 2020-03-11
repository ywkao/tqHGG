#!/bin/bash
# vim: set fdm=marker:
set -e
#===== Functions =====#
# function Preselection(){{{
function Preselection(){
    year="$1"
    log=log/stdout_preselection_${year}.txt
    err=log/stderr_preselection_${year}.txt
    exe=./script/qsubMultiJob_preselection.py; opt=""; #exe=./script/fireBatchJobs.sh; opt="-p";
    time ${exe} ${opt} ${year} > >(tee ${log}) 2> >(tee ${err} >&2)
}
#}}}
# function Intermission(){{{
function Intermission(){
    year=$1
    #./script/check_submit_status.sh -p | tee -a log/stdout.log
    echo "Rename log/stdout.log to log/stdout_pre.log";
    #mv log/stdout.log log/stdout_pre.log
    #mv log/stderr.log log/stderr_pre.log
    #backup condition of preselection
}
#}}}
# function Selection(){{{
function Selection(){
    year="$1"; channel=$2;
    log=log/stdout_selection_${year}_${channel}.txt; touch ${log}
    err=log/stderr_selection_${year}_${channel}.txt; touch ${err}
    exe=./script/qsubMultiJob_selection.py; opt=""; #exe=./script/fireBatchJobs.sh; opt="-s";
    time ${exe} ${opt} ${year} ${channel} > >(tee ${log}) 2> >(tee ${err} >&2)
    sleep 30s; # Depends on the condition of ntugrid5 computers...

    # get sample list
    if   [[ $year == '2016' ]]; then list="lists/list_selection_2016";
    elif [[ $year == '2017' ]]; then list="lists/list_selection";
    elif [[ $year == '2018' ]]; then list="lists/list_selection";
    elif [[ $year == '161718' ]]; then list="lists/list_selection";
    else list="lists/ListRootFiles";
    fi
    echo ${list}

    # monitor processing status
    ./script/monitor.sh ${year} ${channel} ${list}

    # backup log messages
    log_dir="plots/${year}/log" # log dir is created by script/exe_selection_batch.sh
    cp -p src/selection.cpp ${log_dir}
    cp -p ${log} ${log_dir}
    cp -p ${err} ${log_dir}
}
#}}}
# function AfterSelection(){{{
function AfterSelection(){
    year=$1; channel=$2; local log=log/info_stack_plots_${channel}.txt;
    ./script/run_macro_stackPlots.sh ${year} ${channel}  | tee ${log}
    log_dir="plots/${year}/log"
    cp -p src/stackHist*.C ${log_dir}
    cp ${log} ${log_dir}
}
#}}}
# function ReRunStackPlotsOnly(){{{
function ReRunStackPlotsOnly(){
    year="$1"; channel=$2;
    local log=log/info_stack_plots_${channel}.txt
    if [ ! -d plots ]; then echo "[INFO] mv plots_${channel} plots"; mv plots_${channel} plots;
    else echo "[WARNNING] dir plots exists! (abort)"; exit 1;
    fi
    ./script/run_macro_stackPlots.sh ${channel}  | tee ${log}
    cp ${log} plots/${year}/log
    ./script/resetPlotsChannels.sh plots_${channel}
}
#}}}
#------------------------------ Test Section ------------------------------#
#./script/setup_before_fireBatchJobs.sh 
# preselection{{{
#./script/prepareExeForNewMC_npu_float.sh "preselection_npustudy"
#./script/run_macro_stackPlots.sh "hadronic"
#}}}
# selection{{{
#./script/prepareExeForNewMC_npu_float.sh "selection"
#./script/exe_selection_batch.sh TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root "leptonic"
#time ./script/exe_selection_batch.sh "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8" "leptonic"
#time ./script/exe_selection_batch.sh "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8" "hadronic"
#}}}
# top reconstruction study hadronic{{{
#export LD_LIBRARY_PATH=../TopKinFit/:$LD_LIBRARY_PATH ###!!! For topKinFit method purpose.

## TT signal
#make && time ./script/exe_covarianceMatrixStudy.sh "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8.root"
#make && time ./script/exe_generalChiSquareStudy.sh "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8.root" "hadronic"
#
#make && time ./script/exe_covarianceMatrixStudy.sh "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8.root"
#make && time ./script/exe_generalChiSquareStudy.sh "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8.root" "hadronic"
#
#make && time ./script/exe_covarianceMatrixStudy.sh "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root"
#make && time ./script/exe_generalChiSquareStudy.sh "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root" "hadronic"
#
#make && time ./script/exe_covarianceMatrixStudy.sh "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root"
#make && time ./script/exe_generalChiSquareStudy.sh "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root" "hadronic"
#
#
## ST signal
#make && time ./script/exe_covarianceMatrixStudy.sh "ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root"
#make && time ./script/exe_generalChiSquareStudy.sh "ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root" "hadronic"
#
#make && time ./script/exe_covarianceMatrixStudy.sh "ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root"
#make && time ./script/exe_generalChiSquareStudy.sh "ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root" "hadronic"
#}}}
# top reconstruction study leptonic{{{
#export LD_LIBRARY_PATH=../TopKinFit/:$LD_LIBRARY_PATH ###!!! For topKinFit method purpose.
#make && time ./script/exe_generalChiSquareStudy.sh "TT_FCNC-T2HJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8.root" "leptonic"
#make && time ./script/exe_generalChiSquareStudy.sh "TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8.root" "leptonic"
#make && time ./script/exe_generalChiSquareStudy.sh "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root" "leptonic"
#make && time ./script/exe_generalChiSquareStudy.sh "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root" "leptonic"
#make && time ./script/exe_generalChiSquareStudy.sh "ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root" "leptonic"
#make && time ./script/exe_generalChiSquareStudy.sh "ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root" "leptonic"
#}}}
#------------------------- Main Exe Section -------------------------#
#./script/setup_before_fireBatchJobs.sh 
#./script/prepareExeForNewMC_npu_float.sh "preselection"
#for year in "2017old" "2018" "2017" "2016"
#do
#    Preselection "${year}"
#done
#Intermission # not yet
#time ./script/convert_outputPreselection_to_inputSelection.sh

#--------------------------------------------------#
#export LD_LIBRARY_PATH=../TopKinFit/:$LD_LIBRARY_PATH #!For topKinFit method
#./script/prepareExeForNewMC_npu_float.sh "selection"

tag="_latest"
#tag="_update"
for channel in "leptonic" "hadronic"
do
    echo "mv plots_${channel}${tag} plots"; mv plots_${channel}${tag} plots; # re-stack
    #echo "mv plots_${channel} plots"; mv plots_${channel} plots; # re-stack
    #for year in "2017old" "161718" "2018" "2017" "2016"
    for year in "161718"
    do
        #Selection ${year} ${channel}
        #if [ $year == "2016" ]; then
        #    ./script/warning_use2017VHToGG_as2016VHToGG.sh
        #    ls "plots/2016"
        #fi
        AfterSelection ${year} ${channel}
    done
    #./script/resetPlotsChannels.sh plots_${channel}
    ./script/resetPlotsChannels.sh plots_${channel}${tag}
    #./script/resetPlotsChannels.sh ${target}
done
./script/tableMaker.sh "${tag}"
#--------------------------------------------------
#target_dir="plots_leptonic_latest/161718/mva"

### MVA training ###
### please edit before training
### 1. src/mva_TMVAClassification_leptonic.C
### 2. script/train_MVA_method.sh
#time ./script/train_MVA_method.sh
#cp src/mva_TMVAClassification_leptonic.C ${target_dir}

### MVA application (I) ###
### please check following steps before submission
### 1. put the dataset/weights in a proper directory
### 2. src/app_TMVAClassificationApplication.C # Directory, tag, var
#cp -r dataset_testALL_24 ${target_dir}
#time ./script/qsubMultiJob_mvaApplication.py

### stack plots (check significance) ###
### please check following steps before stacking
### 1. include/stack_mva_output_score.h # TAG, DIR
### 2. isCoraser=false/true
### 3. mkdir -p plots_leptonic_latest/161718/mva/plots
### 4. script/run_macro_stackPlots_mvaScore.sh # tag
#time ./script/run_macro_stackPlots_mvaScore.sh

### MVA application (II) ###
### please check following steps before submission
### 1. src/app_TMVAClassificationApplication.C # FINAL_SELECTION
#time ./script/qsubMultiJob_mvaApplication.py
#cp src/app_TMVAClassificationApplication.C ${target_dir}

### stack plots (check final yields) ###
#time ./script/run_macro_stackPlots_mvaScore.sh

### upperlimit ###
### /wk_cms2/ykao/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/
### 1. quick.sh # check tag, str, signalRater
### 2. upperlimit/macro.C # check xtitle
#time ./quick.sh

#--------------------------------------------------
#Selection 161718 "leptonic" #selection and make plots for leptonic channel
#./script/tableMaker.sh
#ReRunStackPlotsOnly 161718 "hadronic"
#ReRunStackPlotsOnly 161718 "leptonic"
#./script/tableMaker.sh
#./script/check_errors.sh

