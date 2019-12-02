#!/bin/bash
set -e
CHANNEL=$1
if [[ $CHANNEL == "" ]]; then
    echo "[WARNNING] Please input valid option (hadronic/leptonic)"
    exit 0
fi
root -l -b -q src/stackHist.C\(\"${CHANNEL}\"\) | tee log/info_stack_plots_${CHANNEL}
