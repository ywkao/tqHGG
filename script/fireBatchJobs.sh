#!/bin/bash
set -e; bool_exe=true; #bool_exe=false;
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
    option=$1; year="$2"; channel=$3;
    if   [[ $option == '-p' ]]; then message="[MESSAGE] Start submit preselection!"; exe=./script/exe_preselection_batch.sh;
    elif [[ $option == '-s' ]]; then message="[MESSAGE] Start submit selection!"   ; exe=./script/exe_selection_batch.sh   ;
    elif [[ $option == '-c' ]]; then message="[MESSAGE] Start submit check work!"  ; exe=./script/doCheckYields            ;
    else echo "[WARNING] something went wrong in the file: ./script/fireBatchJobs.sh"; exit 1;
    fi

    if   [[ $option == '-p' ]] && [[ $year == '2016' ]]; then list="./lists/list_samples_2016";
    elif [[ $option == '-p' ]] && [[ $year == '2017' ]]; then list="./lists/list_samples_2017";
    elif [[ $option == '-p' ]] && [[ $year == '2018' ]]; then list="./lists/list_samples_2018";
    elif [[ $option == '-s' ]] && [[ $year == '2016' ]]; then list="./lists/list_selection_2016";
    elif [[ $option == '-s' ]] && [[ $year == '2017' ]]; then list="./lists/list_selection";
    elif [[ $option == '-s' ]] && [[ $year == '2018' ]]; then list="./lists/list_selection";
    elif [[ $option == '-s' ]] && [[ $year == '161718' ]]; then list="./lists/list_selection"; # selection stage only
    else list="./lists/ListRootFiles";
    fi
    echo ${list}
    #--------------------------------------------------
    if [ ${bool_exe} = true ]; then echo ${message}; fi
    for dataset in `cat ${list} | grep -v "#"`
    do
        if [[ $dataset == 'DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa' ]]; then tag="directory"; else tag="rootfile"; fi # Note: the tag for 2017old is not correct here though, modifies in exe_pre*.sh
        if [[ $option  == '-s' ]]; then tag="rootfile"; fi # Note: in selection stage, all files are rootfiles
        number=$(GetNumber);
        if [ ${bool_exe} = true ]; then
            echo "[MESSAGE] ssh node${number} dataset=${dataset}"
            command="source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd; eval $(scramv1 runtime -sh) time ${exe} ${dataset} ${year} ${tag} ${channel}"
            ssh node${number} ${command} & # self submission
        fi
        # test{{{
        #--- for test only ---#
        if [ ${bool_exe} = false ]; then
            command="cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; time ${exe} ${dataset} ${year} ${tag} ${channel}"
            echo ${command}
        fi
        #}}}
    done
}
#}}}
# warning message{{{
#--------------- Warning Message ---------------#
if [ $# -eq 0 ]; then
    echo "[WARNING] Please input an option."
    echo "./script/fireBatchJobs -c  # Check Entries/Yields"
    echo "./script/fireBatchJobs -p  # Preselection"
    echo "./script/fireBatchJobs -s  # Selection"
    exit 0
fi
#}}}

#--------------- Submit ---------------#
option=$1; year=$2; channel=$3; submit ${option} ${year} ${channel}
