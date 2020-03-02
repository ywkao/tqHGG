#!/bin/bash
set -e
year=$1; channel=$2;

if [[ $channel == "" ]]; then
    echo "[WARNNING] Please input valid option (hadronic/leptonic)"
    exit 0
fi

#------------------------- setup -------------------------#

string_ori="const char TARGET_DIR[128] = \"plots/161718\";"
string_mod="const char TARGET_DIR[128] = \"plots/${year}\";"
string_year_ori="const char year[32] = \"161718\";"
string_year_mod="const char year[32] = \"${year}\";"
header_file="include/stack.h"

if [ $year == "2017old" ]; then
    command="src/stackHist_2017old.C(\"${channel}\")"
else
    command="src/stackHist.C(\"${channel}\")"
fi

#------------------------- execution -------------------------#
sed -i "7c${string_mod}" ${header_file}
sed -i "10c${string_year_mod}" ${header_file}
root -l -b -q ${command}
sed -i "7c${string_ori}" ${header_file}
sed -i "10c${string_year_ori}" ${header_file}
