#!/bin/bash
set -e
dir='/wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG'
log=""

WARNING_MESSAGE="[WARNNING] please add either -p or -s for the status of preselection and selection respectively"

if [ $# -ne 1 ]; then echo ${WARNING_MESSAGE}; exit 1; fi

if   [ $1 == "--preselection" ] || [ $1 == "-p" ]; then log=$dir/log/stdout.log;
elif [ $1 == "--selection" ]    || [ $1 == "-s" ]; then log=$dir/log/stdout_selection.log;
else echo ${WARNING_MESSAGE}; exit 1;
fi

grep submit ${log}
NUM=`cat ListRootFiles | grep -v "#" | nl | cut -f 1 | tail -n1 | tr -d " "`
echo "[CHECK] Finished files: `grep -c end ${log}` / ${NUM}"
