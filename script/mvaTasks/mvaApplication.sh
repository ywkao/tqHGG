#!/bin/bash
set -e
dataset=$1; tmpDir=src/tmp/app_${dataset};

# function{{{
function train(){
    channel=$1; dataset=$2; log=log/log_tmva_application_${channel}_${dataset};
    # execution
    time root -l -b -q ${tmpDir}/TMVAClassificationApplication.C | tee ${log}
}
#}}}

train "leptonic" ${dataset}
