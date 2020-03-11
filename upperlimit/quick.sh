#!/bin/bash
set -e
#tag="hct_tt"
#tag="24_st_hut"
#tag="24_st_hct"
#tag="25_st_hct"
tag="25_st_hut"
channel="lep"

dir="./upperlimit/${tag}_${channel}"
log="${dir}/record.log"
str="fcnc_ST_leptonic_counting"
# function exe(){{{
function exe(){
    txt=$1
    file=${dir}/${txt}; echo ${file} | tee -a ${log}
    sed -n "11p" ${file} | tee -a ${log}; ### double check
    combine -M AsymptoticLimits ${file} | tee -a ${log}
}
#}}}
#function get2digit(){{{
function get2digit(){
    num=$1
    if [ $num -lt 10 ]; then echo "0${num}";
    else echo "${num}";
    fi
}
#}}}
#function createLine(){{{
function createLine(){
    content=$1
    echo ${content} | tee -a ${headerFile}
}
#}}}
# prepare for scanning (once) (skipped) {{{
src="${dir}/${str}.txt"
#signalRate="2.000"
#signalRate="0.967"
#signalRate="0.958"
#signalRate="0.871"
signalRate="1.104"
for i in {1..10}
do
    ans=`get2digit $i`
    target="${dir}/${str}_br0.${ans}.txt"; echo ${target};
    cp ${src} ${target}
    num=`printf "%.3f\n" $(echo "scale=7; ${signalRate}*(${i})" | bc)`; #echo $num;
    sed -i "s/${signalRate}/${num}/g" ${target}
    sed -n "11p" ${target}; ### double check
done
#}}}
# scanning{{{
echo "quick.sh::start!" | tee ${log}
exe fcnc_ST_leptonic_counting_br0.01.txt
exe fcnc_ST_leptonic_counting_br0.02.txt
exe fcnc_ST_leptonic_counting_br0.03.txt
exe fcnc_ST_leptonic_counting_br0.04.txt
exe fcnc_ST_leptonic_counting_br0.05.txt
exe fcnc_ST_leptonic_counting_br0.06.txt
exe fcnc_ST_leptonic_counting_br0.07.txt
exe fcnc_ST_leptonic_counting_br0.08.txt
exe fcnc_ST_leptonic_counting_br0.09.txt
exe fcnc_ST_leptonic_counting_br0.10.txt
#}}}
# prepare signalStrength.h{{{
headerFile="${dir}/signalStrength.h"
echo "//--- Automatically created by quick.sh ---//" | tee ${headerFile}
createLine "const int NUM = 10;"

#content=`grep "upper" ${log} | awk -F "br" '{print $2}'\
#| awk -F ".txt" '(FNR==1){printf("double BF[NUM] = { %.2f, ", $1)}     (FNR==10){print $1" };"}     (FNR!=1 && FNR!=10){print $1","}'\
#| xargs -n 10 echo`
content="double BF[NUM] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };"
createLine "${content}"

content=`grep "Expected  2" ${log}\
| awk '(FNR==1){print "double Expected_025[NUM] = { "$5", "}     (FNR==10){print $5" };"}     (FNR!=1 && FNR!=10){print $5","}'\
| xargs -n 10 echo`
createLine "${content}"

content=`grep "Expected 16" ${log}\
| awk '(FNR==1){print "double Expected_160[NUM] = { "$5", "}     (FNR==10){print $5" };"}     (FNR!=1 && FNR!=10){print $5","}'\
| xargs -n 10 echo`
createLine "${content}"

content=`grep "Expected 50" ${log}\
| awk '(FNR==1){print "double Expected_500[NUM] = { "$5", "}     (FNR==10){print $5" };"}     (FNR!=1 && FNR!=10){print $5","}'\
| xargs -n 10 echo`
createLine "${content}"

content=`grep "Expected 84" ${log}\
| awk '(FNR==1){print "double Expected_840[NUM] = { "$5", "}     (FNR==10){print $5" };"}     (FNR!=1 && FNR!=10){print $5","}'\
| xargs -n 10 echo`
createLine "${content}"

content=`grep "Expected 97" ${log}\
| awk '(FNR==1){print "double Expected_975[NUM] = { "$5", "}     (FNR==10){print $5" };"}     (FNR!=1 && FNR!=10){print $5","}'\
| xargs -n 10 echo`
createLine "${content}"
#}}}

template_macro="./upperlimit/macro.C"
cp ${template_macro} ${dir}

command="${dir}/macro.C(\"${tag}\")"
root -l -b -q "${command}"


#test{{{
#combine -M AsymptoticLimits realistic-counting-experiment.txt
#combine -M AsymptoticLimits ${dir}/fcnc_TT_leptonic_counting.txt
#}}}
