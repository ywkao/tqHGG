#!/bin/bash
#ref of command awk: https://ubuntuforums.org/showthread.php?t=840054
set -e

INPUTDIR="/wk_cms2/youying/public/tH_FCNC/Era2017_RR-31Mar2018_v2"
#INPUTDIR="/wk_cms2/youying/public/2017_94X_3_1_X_and_3_2_0"
INPUTDIR_ywk="/wk_cms2/ykao/public/2017_94X_3_1_X_and_3_2_0"
OUTPUTDIR="ntuples_skimmed"
EXECUTABLE=./bin/covarianceMatrixStudy

ExeAnalysis(){
    file=$1
    log=log/log_${file}.txt
    echo "[MESSAGE] Start to analyze ${file}" | tee ${log} 
    echo $file | awk -F "." '{print "'$INPUTDIR'""/"$1".root", "'$OUTPUTDIR'""/covMatrix_ntuple_"$1".root", $1}' |\
    xargs -n3 ${EXECUTABLE} | tee -a ${log}
    #xargs -n3 ${EXECUTABLE} | grep -v ientry | tee -a ${log}
    # REMARK: input_file / output_file / datasets / output_dir
}
ExeAnalysis_multi(){
    file=$1
    log=log/log_${file}.txt
    echo "[MESSAGE] Start to analyze ${file}" | tee ${log} 
    echo $file | awk -F "." '{print "'$INPUTDIR_ywk'""/"$1, "'$OUTPUTDIR'""/covMatrix_ntuple_"$1".root", $1}' |\
    xargs -n3 ${EXECUTABLE} | tee -a ${log}
    # REMARK: input_file / output_file / datasets / output_dir
}
ExeAnalysis_singletop(){
    file=$1
    log=log/log_${file}.txt
    echo "[MESSAGE] Start to analyze ${file}" | tee ${log} 
    echo $file | awk -F "." '{print "'$INPUTDIR_ywk'""/"$1".root", "'$OUTPUTDIR'""/covMatrix_ntuple_"$1".root", $1}' |\
    xargs -n3 ${EXECUTABLE} | tee -a ${log}
    # REMARK: input_file / output_file / datasets / output_dir
}

#--------------- DryRun ---------------#
if [[ $1 == '-d' || $1 == '--dryRun' ]]; then
    #ExeAnalysis ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root
    ExeAnalysis GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8.root
    echo "[MESSAGE] This is the end of dryRun execution." && exit 0
fi

#--------------- Execution ---------------#
FILE=$1
if [[ $FILE == 'DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa' || $FILE == 'GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8' ]]; then
    ExeAnalysis_multi $FILE;
    echo "[MESSAGE] This is the end of $FILE" && exit 0
elif [[
        $FILE == 'ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root' ||\
        $FILE == 'ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root' ||\
        $FILE == 'ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root' ||\
        $FILE == 'ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root'
     ]]; then
    ExeAnalysis_singletop $FILE;
    echo "[MESSAGE] This is the end of $FILE" && exit 0
else
    ExeAnalysis $FILE;
    echo "[MESSAGE] This is the end of $FILE" && exit 0
fi 
