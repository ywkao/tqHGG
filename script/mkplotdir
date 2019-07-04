#!/bin/bash
set -e
if [ ! -d bin ]; then
    mkdir bin;
fi

if [ ! -d ntuples_skimmed ]; then
    mkdir ntuples_skimmed;
fi

if [ ! -d plots ]; then
    mkdir plots;
fi

if [ ! -d plots/log_scale ]; then
    mkdir plots/log_scale;
fi

for file in `cat ListRootFiles | grep -v "#"`
do
    DIRECTORY=plots/"`echo $file | awk -F "." '{print $1}'`"
    if [ ! -d "${DIRECTORY}" ]; then
        mkdir ${DIRECTORY}
    fi
done

if [ ! -d log ]; then
    mkdir log;
fi

echo "[MESSAGE] Directory are checked and prepared!"
