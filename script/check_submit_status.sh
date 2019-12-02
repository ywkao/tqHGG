#!/bin/bash
set -e
dir='/wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG'
#log=$dir/log/stdout.log
log=$dir/log/stdout_selection.log
grep submit ${log}
NUM=`cat ListRootFiles | grep -v "#" | nl | cut -f 1 | tail -n1 | tr -d " "`
echo "[CHECK] Finished files: `grep -c end ${log}` / ${NUM}"


