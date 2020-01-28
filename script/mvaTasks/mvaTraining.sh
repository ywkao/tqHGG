#!/bin/bash
set -e
tag=$1
tmpDir=src/tmp/task_${tag}

# function{{{
function train(){
    channel=$1; tag=$2;
    dir=plots_${channel}/mva/
    log=log/log_tmva_${channel}_${tag}

    # execution
    time root -l -b -q ${tmpDir}/TMVAClassification.C | tee ${log}

    # backup log and dataset
    cp -p ${log} ${dir}
    cp -rp dataset_${tag} ${dir}
}
#}}}

#train "hadronic"
train "leptonic" ${tag}
