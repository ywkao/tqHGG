#!/bin/bash
#for code in {01..09} {11..14} 16 {18..20}
#command="ps aux | grep ykao | grep pre | grep -v source"
command="ps aux | grep ykao"
for code in {01..09} {10..20}
do
    echo "[INFO] in node"$code"..."
    num_running=`ssh node${code} "${command} | grep -c \" R \""`
    if [[ ${num_running} -eq 0 ]];
    then
        echo "[INFO] no exe running.";
    else
        echo "[KILLING]";
        ssh node${code} "${command} | grep \" R \""
        #ssh node${code} "ps aux | grep ykao | grep pre | grep -v source | grep \" R \" | awk '{print \$2}' | xargs -n1 kill -9"
    fi
done
