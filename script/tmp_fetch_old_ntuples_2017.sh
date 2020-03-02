#!/bin/bash
set -e
channel=$1
for dir in `ls plots_${channel}_latest/ | grep -v pdf | grep -E -i "fcnc|di" | grep -v Inf`
do
    #echo plots_${channel}_latest/${dir}
    mkdir plots/161718/${dir};
    cp -p plots_${channel}_latest/${dir}/*root plots/161718/${dir};
done
