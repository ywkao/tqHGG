# tqHGG

Searches for top FCNH with H decays to diphoton analysis. (2016)

Although the analysis code is performed under CMSSW_9_4_10,
it should be independent of the version of CMSSW because the only used tool is ROOT.

Usage:
1. Setup the cmsenv #In order to use ROOT. <br />
   $ source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv
2. Check the directory of input root files in the shell script, script/exe_preselection_batch.
3. Perform preselection with parallel computing. (currently, on ntugrid5 only.) <br />
    <br />
   Step 3-1: in ~/.bashrc, comment out the command(s) related to source any cvmfs. 
   This is meant to prevent unexpected errors when using nodes for computation.<br />
   e.g. #source /cvmfs/cms.cern.ch/cmsset_default.sh <br />
    <br />
   Step 3-2: dryRun test <br />
   $ ./fireBatchJobs -d <br />
   NOTE: if there is any error message, please use another new terminal to login ntugrid5, set up the environment, and repeat this step again.
   If there is any further problems, please contact me, ykao@cern.ch . <br />
    <br />
   Step 3-3: check the option in the script, execution. Then, execute it. <br />
   $ ./execution <br />
   NOTE: execution is kind of for lazy operation. It will also output log message (log/stdout.log), such as on which node the datasets are processed.
    <br />
4. Perform selection with parallel computing. (currently, on ntugrid5 only.) <br />
   Step 4-1: check the option in the script, execution. Then, execute it. <br />
   $ ./execution <br />
5. root -l -b -q src/stackHist.C #Make stack plots.
6. ./makePUweight.py -x 69200 #Produce rootfile of pileup reweighting factor ~5min <br />
    - Step 6 is optional. Need not execute once having its output file, data/MCPileUp.root<br />
    - [Python module](https://github.com/cms-sw/cmssw/tree/master/SimGeneral/MixingModule/python?fbclid=IwAR2ehfE0hR8uaewPro4vQXos5I_IU6O7cyrtefQxTT6bMpyMETCTzpSuK58): 'mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi' <br />
    - [Golden json](https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2017Analysis#13_TeV_pp_runs_ReReco): 'data/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt' <br />
    - [Pileup json](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Pileup_JSON_Files_For_Run_II): '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt'
