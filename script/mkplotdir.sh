#!/bin/bash
# vim: set fdm=marker:
set -e
OUTPUTDIR="/wk_cms/ykao/tqHGG/2017_oldNtuples/ntuples_skimmed"
OUTPUTDIR_2016="/wk_cms/ykao/tqHGG/2016/ntuples_skimmed"
OUTPUTDIR_2017="/wk_cms/ykao/tqHGG/2017/ntuples_skimmed"
OUTPUTDIR_2018="/wk_cms/ykao/tqHGG/2018/ntuples_skimmed"

mkdir -p bin;
mkdir -p log;
# direcotry ntuples_skimmed{{{
#mkdir -p ntuples_skimmed/chi2_study_leptonic_1D_plots; # for output from chi-2 study
#}}}
# direcotry plots {{{
#mkdir -p plots/log_scale;
#mkdir -p plots/mva;
#mkdir -p plots/log;

#for file in `cat ListRootFiles | grep -v "#"`
#do
#    DIRECTORY=plots/"`echo $file | awk -F "." '{print $1}'`"
#    mkdir -p ${DIRECTORY}
#done
#}}}

echo "[MESSAGE] Directory are checked and prepared!"
