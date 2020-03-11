#!/bin/bash
set -e
dataset=$1
dirtype=$2
tmpDir=src/tmp/app_${dataset}
macro="TMVAClassificationApplication.C"
file="src/app_${macro}"
target=${tmpDir}/${macro}

# function{{{
function Edit(){
    line=$1; pre=$2; content=$3;
    sed -i -e "${line}c ${content}" ${target}
    sed -i -e ''"$line"'s/^'"$pre"'/    '"$pre"'/g' ${target}
}
#}}}

mkdir -p ${tmpDir}; cp -p ${file} ${target}; ls ${target}

content_01="const char dataset[256] = \"${dataset}\";"
Edit 133 "const" "${content_01}"; sed -n "133p" ${target} # check if it is correctly modified

content_02="TString fname = Form(\"%s/tree_%s.root\", ${dirtype}, dataset);"
Edit 137 "TString" "${content_02}"; sed -n "137p" ${target} # check if it is correctly modified
