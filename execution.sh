#!/bin/bash
set -e

#===== Preselection =====#
function Preselection_npustudy(){
    time ./fireBatchJobs.sh -npustudy > >(tee log/stdout.log) 2> >(tee log/stderr.log >&2)
}
function Preselection(){
    time ./fireBatchJobs.sh -p > >(tee log/stdout.log) 2> >(tee log/stderr.log >&2)
}
function Intermission(){
    ./check_submit_status.sh | tee -a log/stdout.log
    echo "Rename log/stdout.log to log/stdout_pre.log";
    mv log/stdout.log log/stdout_pre.log
    mv log/stderr.log log/stderr_pre.log
}

#===== Selection =====#
function Selection(){
    CHANNEL=$1
    time ./fireBatchJobs.sh -s ${CHANNEL} > >(tee log/stdout_selection.log) 2> >(tee log/stderr_selection.log >&2)
    sleep 1m; # Depends on the condition of ntugrid5 computers...
    NUM=`cat ListRootFiles | grep -v "#" | nl | cut -f 1 | tail -n1 | tr -d " "`
    while [ `grep -ic end log/stdout_selection.log` != ${NUM} ]
    do
        echo "[INFO-execution] Wait another 10 sec.."
        sleep 10s;
    done
    ./check_submit_status.sh | tee -a log/stdout_selection.log
}
function AfterSelection(){
    CHANNEL=$1
    ./script/run_macro_stackPlots.sh ${CHANNEL}
    ./script/resetPlotsChannels.sh plots_${CHANNEL}
    cp -p log/stdout_selection.log plots_${CHANNEL}/log/stdout_${CHANNEL}.log
    cp -p log/stderr_selection.log plots_${CHANNEL}/log/stderr_${CHANNEL}.log
    cp -p log/info_stack_plots_${CHANNEL} plots_${CHANNEL}/log
    mv log/stdout_selection.log log/stdout_${CHANNEL}.log
    mv log/stderr_selection.log log/stderr_${CHANNEL}.log
}
function ReRunStackPlotsOnly(){
    CHANNEL=$1
    if [ ! -d plots ]; then echo "[INFO] mv plots_${CHANNEL} plots"; mv plots_${CHANNEL} plots;
    else echo "[WARNNING] dir plots exists! (abort)"; exit 1;
    fi
    ./script/run_macro_stackPlots.sh ${CHANNEL}
    ./script/resetPlotsChannels.sh plots_${CHANNEL}
    cp -p log/info_stack_plots_${CHANNEL} plots_${CHANNEL}/log
}


#------------------------------ Test Section ------------------------------#
#./script/prepareExeForNewMC_npu_float.sh "preselection_npustudy"
#Preselection_npustudy
#./script/run_macro_stackPlots.sh "hadronic"

export LD_LIBRARY_PATH=../TopKinFit/:$LD_LIBRARY_PATH ###!!! For topKinFit method purpose.
#make && time ./script/exe_covarianceMatrixStudy.sh "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8.root"
#make && time ./script/exe_generalChiSquareStudy.sh "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8.root" "hadronic"
make && time ./script/exe_generalChiSquareStudy.sh "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root" "leptonic"

#./script/prepareExeForNewMC_npu_float.sh "selection"
#./script/exe_selection_batch.sh TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root "leptonic"
#time ./script/exe_selection_batch.sh "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8" "leptonic"
#time ./script/exe_selection_batch.sh "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8" "hadronic"


#------------------------- Main Exe Section -------------------------#
#./script/prepareExeForNewMC_npu_float.sh "preselection"
#Preselection
#Intermission

#./script/prepareExeForNewMC_npu_float.sh "selection"
#Selection "hadronic" #selection and make plots for hadronic channel
#AfterSelection "hadronic" 
#Selection "leptonic" #selection and make plots for leptonic channel
#AfterSelection "leptonic" 
#./script/check_errors.sh
#./script/tableMaker.sh

#ReRunStackPlotsOnly "hadronic"
#ReRunStackPlotsOnly "leptonic"
