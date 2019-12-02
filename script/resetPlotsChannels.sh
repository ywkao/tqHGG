#!/bin/bash
set -e

targetDir=$1
if [ -d ${targetDir} ]; then
    echo "[WARNNING] rm -r ${targetDir}; mv plots ${targetDir}"
    rm -r ${targetDir}; mv plots ${targetDir}
    exit 0;
else
    echo "[MESSAGE] mv plots ${targetDir}"
    mv plots ${targetDir}
    exit 0;
fi

