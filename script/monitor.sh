#!/bin/bash
#set -e # should be commented out for the use of grep "str" file = nothing

function monitor(){
    local year="$1"; local channel=$2; local list=$3;
    local NUM=`cat ${list} | grep -v "#" | nl | cut -f 1 | tail -n1 | tr -d " "`
    local log_monitor="log/monitor.txt"; cat ${list} | grep -v "#"> ${log_monitor};
    echo "[INFO] script/monitor.sh: Expected = ${NUM}"

    #--------------------------------------------------#
    local instant_num=0
    while [ ${instant_num} -ne ${NUM} ]
    do
        # scan each dataset
        for file in `cat ${list} | grep -v "#"`;
        do
            _dataset_=`echo ${file} | awk -F "." '{print $1}'`
            _log_=log/log_selection_${year}_${_dataset_}_${channel}.txt

            # check existence of a file and a keyword
            if [ -f ${_log_} ] && [ `grep -ic yield ${_log_}` ]; then
                _num_=`grep -ic yield ${_log_}`
            else
                _num_=0;
            fi #; echo "${_log_}: ${_num_}"

            # update the log/monitor.txt
            if [ ${_num_} -eq 1 ]; then sed -i "s/^${file}/[V] ${file}/g" ${log_monitor}; fi 
        done
        # record num of finished jobs
        instant_num=`grep -ic "\[V\]" ${log_monitor}`
        if [ ${instant_num} -ne ${NUM} ]; then
            echo "[INFO] script/monitor.sh: instant report && wait for another 10 sec.."
            cat ${log_monitor}; sleep 10s;
        fi
    done
    #--------------------------------------------------#
    echo "[INFO] script/monitor.sh: Finished!"
    cat ${log_monitor}
}

monitor $1 $2 $3

#counter=0
#while [ $counter -ne 3 ]
#do
#    counter=$((counter+1))
#    echo "[INFO] ${counter} check..."
#    cat log/log_selection_161718_WW_hadronic.txt
#    sleep 2s;
#    grep -c yield log/log_selection_161718_WW_hadronic.txt
#done
