#!/bin/bash
set -e
#channel="hadronic"
channel="leptonic"
src_file="src/generalChiSquareStudy_${channel}.cpp"
exe_file="src/generalChiSquareStudy_${channel}_exe.cpp"
rootfile="null.root"

#function evaluate_cov(){{{
function evaluate_cov(){
    log=$1
    make && time ./script/exe_covarianceMatrixStudy.sh ${rootfile} ${channel} | tee ${log}
}
#}}}
#function execution(){{{
function execution(){
    log=$1; ext=$2
    export LD_LIBRARY_PATH=../TopKinFit/:$LD_LIBRARY_PATH ###!!! For topKinFit method purpose.
    make && time ./script/exe_generalChiSquareStudy.sh ${rootfile} ${channel} ${ext} | tee ${log}
}
#}}}
#function EditFile(){{{
function EditFile(){
    # example{{{
    #bool bool_bjet_is_loose  = true;
    #bool bool_bjet_is_medium = true;
    #bool bool_bjet_is_tight  = false;
    #bool bool_num_bjets_is_exactly_one = true;
    #}}}

    loose_condition="bool bool_bjet_is_loose  = $1;"
    medium_condition="bool bool_bjet_is_medium = $2;"
    tight_condition="bool bool_bjet_is_tight  = $3;"
    num_condition="bool bool_num_bjets_is_exactly_one = $4;"
    
    sed -i "20c ${loose_condition}"  ${exe_file}
    sed -i "21c ${medium_condition}" ${exe_file}
    sed -i "22c ${tight_condition}"  ${exe_file}
    sed -i "23c ${num_condition}"    ${exe_file}

    execution "result_top_reco_study/log_chi_${channel}_$1_$2_$3_$4" "$1_$2_$3_$4"
}
#}}}

cp -p ${src_file} ${exe_file}

#--------------- Execution ---------------#
#dataset=`echo ${rootfile} | awk -F "." '{print $1}'`
#evaluate_cov "log/covariance_matrix_${dataset}.txt"

#EditFile false false true true
#EditFile false true false true
#EditFile true false false true
#EditFile false false true false
#EditFile false true false false
EditFile true false false false

#echo "[MESSAGE] python script/chi_examine.py"; python script/chi_examine.py;


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
