#!/bin/bash
set -e
channel="leptonic";

if [[ $channel == "" ]]; then
    echo "[WARNNING] Please input valid option (hadronic/leptonic)"
    exit 0
fi

#------------------------- setup -------------------------#
command="src/output_mvaScore_stackHist.C(\"${channel}\")"

#------------------------- execution -------------------------#
log="log/report_output_mvaScore_stackHist.txt"
macro="src/output_mvaScore_stackHist.C"
string_true="bool isCoarser = true;"
string_false="bool isCoarser = false;"

echo "Log from output_mvaScore_stackHist.C" > ${log}
sed -i "17c${string_true}" ${macro}
root -l -b -q ${command} # print info but do not save to the log file
sed -i "17c${string_false}" ${macro}
root -l -b -q ${command} >> ${log}
#
#output_dir=`grep OUTPUT include/stack_mva_output_score.h | awk '{print $5}'`
#echo "src/output_mvaScore_stackHist.C::output ${output_dir}"
#echo "script/run_macro_stackPlots_mvaScore.sh::log ${log}"


###------------------------- init upper limit -------------------------#
tag="25_st_hut_lep"
limit_dir="/wk_cms2/ykao/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/upperlimit"
tag_dir=${limit_dir}/${tag}
mkdir -p ${tag_dir}

tag_file=${tag_dir}/"fcnc_ST_leptonic_counting.txt"
cp ${log} ${tag_file}
sed -i '1,4d' ${tag_file}

ls -lhrt ${tag_dir}
cat ${tag_file}
