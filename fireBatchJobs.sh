#!/bin/bash
set -e
#function GetNumber(){{{
#--------------- Function ---------------#
function GetNumber(){
    number=0
    #while [ $number -eq 0 ] # make sure no node00 occurs
    #while [[ $number -eq 0 || $number -eq 2 || $number -eq 3 || $number -eq 4 || $number -eq 6 || $number -eq 9 || $number -eq 10 || $number -eq 15 || $number -eq 17 ]] # make sure no node00 occurs
    while [[ $number -eq 0 || $number -eq 2 || $number -eq 3 || $number -eq 4 || $number -eq 9 ]] # make sure no node00 occurs
    do 
        number=$((RANDOM%20))
    done

    if [ $number -lt 10 ]
    then
        echo 0$number
    else
        echo $number
    fi
}
if [[ $1 == '-t' ]]; then
    for i in {1..100}; do echo $(GetNumber); done
    exit 0
fi
#}}}
# Setup {{{
#--------------- Setup ---------------#
echo "[MESSAGE] Check directories for plots..." && ./script/mkplotdir.sh
echo "[MESSAGE] Check executable..." && make
echo "[MESSAGE] Ready!"
#}}}
# Option needed {{{
#--------------- Option needed ---------------#
if [[ $1 == '' ]]; then
    echo "[WARNNING] Please input an option."
    echo "./fireBatchJobs -c  # Check Entries/Yields"
    echo "./fireBatchJobs -d  # Dryrun"
    echo "./fireBatchJobs -p  # Preselection"
    echo "./fireBatchJobs -s  # Selection"
    echo "./fireBatchJobs -pu # Pileup study"
    exit 0
fi
#}}}
#--------------- DryRun ---------------#
if [[ $1 == '-d' || $1 == '--dryRun' ]]; then
    number=$(GetNumber)
    echo "[MESSAGE] ssh node${number}"
    #ssh node${number} "source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd; eval `scramv1 runtime -sh` time echo \"Hello World!\"" &
    ssh node${number} "source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd; eval `scramv1 runtime -sh` time ./script/exe_preselection_batch.sh -d" &
    #ssh node${number} "source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd; eval `scramv1 runtime -sh` time ./script/exe_selection_batch.sh -d" &
    #ssh node${number} "source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd; eval `scramv1 runtime -sh` time ./script/exe_preselection_npustudy_batch -d" &
    #ssh node${number} "source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd; eval `scramv1 runtime -sh` time ./script/doCheckYields \"ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root\"" &
    #ssh node${number} "source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd; eval $(scramv1 runtime -sh) time ./script/doCheckYields \"ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root\"" &
    #ssh node${number} "source /wk_cms/ykao/root/bin/thisroot.sh; cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd; time ./script/doCheckYields \"ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root\"" &
    echo "[MESSAGE] This is the end of dryRun execution." && exit 0
#--------------- Submit ---------------#
#elif [[ $1 == '-pu' ]]; then
#    echo "[MESSAGE] Start submit pustudy!"
#
#    if [ ! -d ntuples_pustudy ]; then mkdir ntuples_pustudy; fi
#
#    for dataset in `cat ListRootFiles | grep -v "#"`
#    do
#        echo "[MESSAGE] ssh node$(GetNumber) dataset=${dataset}"
#        ssh node$(GetNumber) "source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd;\
#        eval `scramv1 runtime -sh` time ./script/exe_pustudy_batch.sh ${dataset}" &
#    done
#--------------------------------------#
elif [[ $1 == '-c' ]]; then
    echo "[MESSAGE] Start submit check work!"
    for dataset in `cat ListRootFiles | grep -v "#"`
    do
        number=$(GetNumber)
        echo "[MESSAGE] ssh node${number} dataset=${dataset}"
        #ssh node$(GetNumber) "cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd; time ./script/doCheckYields ${dataset}" &
        ssh node${number} "source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd;\
        eval $(scramv1 runtime -sh) time ./script/doCheckYields ${dataset}" &
    done
#--------------------------------------#
elif [[ $1 == '-npustudy' ]]; then
    echo "[MESSAGE] Start submit preselection_npustudy!"
    for dataset in `cat ListRootFiles | grep -v "#"`
    do
        number=$(GetNumber)
        echo "[MESSAGE] ssh node${number} dataset=${dataset}"
        ssh node${number} "source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd;\
        eval `scramv1 runtime -sh` time ./script/exe_preselection_npustudy_batch ${dataset}" &
    done
#--------------------------------------#
elif [[ $1 == '-p' ]]; then
    echo "[MESSAGE] Start submit preselection!"
    for dataset in `cat ListRootFiles | grep -v "#"`
    do
        number=$(GetNumber)
        echo "[MESSAGE] ssh node${number} dataset=${dataset}"
        ssh node${number} "source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd;\
        eval `scramv1 runtime -sh` time ./script/exe_preselection_batch.sh ${dataset}" &
    done
#--------------------------------------#
elif [[ $1 == '-s' ]]; then
    CHANNEL=$2
    echo "[MESSAGE] Start submit selection! (${CHANNEL})"
    for dataset in `cat ListRootFiles | grep -v "#"`
    do
        number=$(GetNumber)
        echo "[MESSAGE] ssh node${number} dataset=${dataset}"
        ssh node${number} "source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd;\
        eval `scramv1 runtime -sh` time ./script/exe_selection_batch.sh ${dataset} ${CHANNEL}" &
    done
else
    echo "[WARNNING] Unknown parameters."
fi

#ssh node01 "source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; pwd; eval `scramv1 runtime -sh` time ./script/exe_pustudy_batch.sh -d" &
##seq 1 4 | xargs -I{} -n1 -P4 ssh node01 "cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG && time ./script/exePreselection_batch_0{} -d" #problematic
