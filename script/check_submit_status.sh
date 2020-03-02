#!/bin/bash
set -e
dir='/wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG'
log=""

function GetResult(){
    log=$1; list=$2;
    grep submit ${log}
    NUM=`cat ${list} | grep -v "#" | nl | cut -f 1 | tail -n1 | tr -d " "`
    echo "[CHECK] Finished files: `grep -c end ${log}` / ${NUM}"
}

WARNING_MESSAGE="[WARNNING] please add either -p or -s for the status of preselection and selection respectively"
if [ $# -ne 1 ]; then echo ${WARNING_MESSAGE}; exit 1; fi

if   [ $1 == "--preselection" ] || [ $1 == "-p" ]; then log=$dir/log/stdout.log;
elif [ $1 == "--selection" ]    || [ $1 == "-s" ]; then log=$dir/log/stdout_selection.log;
else echo ${WARNING_MESSAGE}; exit 1;
fi

for year in 2016 2017 2018 "2017old";
do
    log=$dir/log/stdout_${year}.log
    if   [[ $year == '2016' ]]; then list="x_check_new_ntuples/list_samples_2016";
    elif [[ $year == '2017' ]]; then list="x_check_new_ntuples/list_samples_2017";
    elif [[ $year == '2018' ]]; then list="x_check_new_ntuples/list_samples_2018";
    else list="ListRootFiles";
    fi

    GetResult ${log} ${list}

done

#grep submit ${log}
#NUM=`cat ListRootFiles | grep -v "#" | nl | cut -f 1 | tail -n1 | tr -d " "`
#echo "[CHECK] Finished files: `grep -c end ${log}` / ${NUM}"
