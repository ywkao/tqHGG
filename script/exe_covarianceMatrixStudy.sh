#!/bin/bash
./bin/covarianceMatrixStudy



## previous scripts{{{
##ref of command awk: https://ubuntuforums.org/showthread.php?t=840054
##INPUTDIR="/wk_cms2/youying/public/tH_FCNC/Era2017_RR-31Mar2018_v2/"
##OUTPUTDIR="result_top_reco_study"
#EXECUTABLE=./bin/covarianceMatrixStudy
#
#ExeAnalysis(){
#    file=$1; channel=$2;
#    exe_cov=${EXECUTABLE}
#    message_cov="[INFO] Start to calculate covariance matrix..."
#    dataset=`echo $file | awk -F "." '{print $1}'`
#    #outputfile="${OUTPUTDIR}/covariance_${channel}_${file}"
#
#    mkdir -p ${OUTPUTDIR}
#
#    # execution
#    #echo ${message_cov}; ${exe_cov} ${outputfile} ${dataset};
#    echo ${message_cov}; ${exe_cov} ${dataset};
#}
##--------------- Execution ---------------#
#FILE="null.root"; CHANNEL="hadronic";
#ExeAnalysis $FILE $CHANNEL;
##}}}
