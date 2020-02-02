#!/bin/bash
set -e
#function GetNumber(){{{
#--------------- Function ---------------#
function GetNumber(){
    number=0
    #while [ $number -eq 0 ] # make sure no node00 occurs
    #while [[ $number -eq 0 || $number -eq 2 || $number -eq 3 || $number -eq 4 || $number -eq 6 || $number -eq 9 || $number -eq 10 || $number -eq 15 || $number -eq 17 ]] # make sure no node00 occurs
    while [[ $number -eq 0 || $number -eq 2 || $number -eq 3 || $number -eq 4 || $number -eq 9 ]] # make sure no node00 occurs
    do 
        number=$((RANDOM%20))
    done

    if [ $number -lt 10 ]
    then
        echo 0$number
    else
        echo $number
    fi
}
if [[ $1 == '-t' ]]; then # test purpose only
    for i in {1..100}; do echo $(GetNumber); done
    exit 0
fi
#}}}
# funciont submit(){{{
function submit() {
    option=$1; channel=$2;
    if   [[ $option == '-p' ]]; then message="[MESSAGE] Start submit preselection!"; exe=./script/exe_preselection_batch.sh;
    elif [[ $option == '-s' ]]; then message="[MESSAGE] Start submit selection!"   ; exe=./script/exe_selection_batch.sh   ;
    elif [[ $option == '-c' ]]; then message="[MESSAGE] Start submit check work!"  ; exe=./script/doCheckYields            ;
    else echo "[WARNING] something went wrong in the file: ./script/fireBatchJobs.sh"; exit 1;
    fi
    #--------------------------------------------------
    echo ${message}
    for dataset in `cat ListRootFiles | grep -v "#"`
    do
        name="test_${dataset}"
        number=$(GetNumber); echo "[MESSAGE] ssh node${number} dataset=${dataset}"
        command="source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd; eval $(scramv1 runtime -sh) time ${exe} ${dataset} ${channel}"
        echo ${command}
        ssh node${number} ${command} & # self submission
        # qsub (not work yet){{{
        #command='"cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd; time ${exe} ${dataset} ${channel}"' # do not work
        #if [[ $option == '-s' ]];
        #then ./script/submitJOB.py --command=${command} --name=${name} # qsub
        #else ssh node${number} ${command} & # self submission
        #fi
        #}}}
    done
}
#}}}
    
#--------------- Setup ---------------#
if [ $# -eq 0 ]; then
    echo "[WARNING] Please input an option."
    echo "./script/fireBatchJobs -c  # Check Entries/Yields"
    echo "./script/fireBatchJobs -p  # Preselection"
    echo "./script/fireBatchJobs -s  # Selection"
    exit 0
fi
echo "[MESSAGE] Check directories for plots..." && ./script/mkplotdir.sh
echo "[MESSAGE] Check executable..." && make
echo "[MESSAGE] Ready!"

#--------------- Submit ---------------#
option=$1; channel=$2; submit ${option} ${channel}
