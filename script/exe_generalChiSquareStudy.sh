#!/bin/bash
#ref of command awk: https://ubuntuforums.org/showthread.php?t=840054

INPUTDIR="/wk_cms2/youying/public/tH_FCNC/Era2017_RR-31Mar2018_v2/"
OUTPUTDIR="result_top_reco_study"
EXECUTABLE_hadronic=./bin/generalChiSquareStudy_hadronic_exe
EXECUTABLE_leptonic=./bin/generalChiSquareStudy_leptonic_exe

ExeAnalysis(){
    file=$1; channel=$2; ext=$3;
    exe=""; if [[ ${channel} == "hadronic" ]]; then exe=${EXECUTABLE_hadronic}; else exe=${EXECUTABLE_leptonic}; fi
    message="[INFO] Start analyze ${file} with ${exe}"
    message_chi="[INFO] Start to estimate reco performance..."
    dataset=`echo $file | awk -F "." '{print $1}'`
    inputfile="${INPUTDIR}/${file}"
    outputfile="${OUTPUTDIR}/chi2study_${channel}_${dataset}_${ext}.root"

    mkdir -p "${OUTPUTDIR}"
    mkdir -p "${OUTPUTDIR}/hist_factory_hadronic"
    mkdir -p "${OUTPUTDIR}/hist_factory_leptonic"

    # execution
    echo ${message}
    echo ${message_chi}; ${exe} ${inputfile} ${outputfile} ${dataset};

    #cp -rp "${OUTPUTDIR}/hist_factory_${channel}" "${OUTPUTDIR}/hist_factory_${channel}_${ext}"
}
#--------------- Execution ---------------#
FILE=$1; CHANNEL=$2; EXT=$3;
ExeAnalysis $FILE $CHANNEL $EXT;
echo "[MESSAGE] This is the end of $FILE $CHANNEL ${EXT}" && exit 0
