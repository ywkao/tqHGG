#!/bin/bash
# vim: set fdm=marker:
set -e
OUTPUTDIR="/wk_cms/ykao/tqHGG/ntuples_skimmed"

if [ ! -d bin ]; then
    mkdir bin;
fi

# direcotry ntuples_skimmed{{{
# for output from preselection
mkdir -p ${OUTPUTDIR}/log;
# for output from chi-2 study
mkdir -p ntuples_skimmed/chi2_study_leptonic_1D_plots;
#}}}

# direcotry plots {{{
if [ ! -d plots ]; then
    mkdir plots;
fi

if [ ! -d plots/log_scale ]; then
    mkdir plots/log_scale;
fi

if [ ! -d plots/mva ]; then
    mkdir plots/mva;
fi

if [ ! -d plots/log ]; then
    mkdir plots/log;
fi
#}}}

for file in `cat ListRootFiles | grep -v "#"`
do
    DIRECTORY=plots/"`echo $file | awk -F "." '{print $1}'`"
    if [ ! -d "${DIRECTORY}" ]; then
        mkdir ${DIRECTORY}
    fi
done

if [ ! -d log ]; then
    mkdir log;
fi

echo "[MESSAGE] Directory are checked and prepared!"
