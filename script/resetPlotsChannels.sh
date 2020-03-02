#!/bin/bash
set -e

targetDir=$1
if [ ! -d plots ]; then
    echo "[WARNNING] directory ./plots does not exist."
    exit 1;
elif [ -d ${targetDir} ]; then
    echo "[WARNNING] ${targetDir} exists. Abort."
    #echo "[WARNNING] rm -r ${targetDir}; mv plots ${targetDir}"
    #rm -r ${targetDir}; mv plots ${targetDir}
    exit 1;
else
    echo "[MESSAGE] mv plots ${targetDir}"
    mv plots ${targetDir}
    exit 0;
fi

