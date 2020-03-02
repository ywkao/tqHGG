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
    sed -i -e "${line}c ${content}" $file
    sed -i -e ''"$line"'s/^'"$pre"'/   '"$pre"'/g' $file
}
#}}}

content_01="const char dataset[256] = \"${dataset}\";"
Edit 133 "const" "${content_01}"; sed -n "133p" ${file} # check if it is correctly modified

content_02="TString fname = Form(\"%s/tree_%s.root\", ${dirtype}, dataset);"
Edit 137 "TString" "${content_02}"; sed -n "137p" ${file} # check if it is correctly modified

mkdir -p ${tmpDir}; cp -p ${file} ${target}; ls ${target}
