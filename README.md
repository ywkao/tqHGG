# tqHGG
t to cH analysis

Searches for top FCNH with H decays to diphoton analysis.

Although the analysis code is performed under CMSSW_9_4_10,
it should be independent of the version of CMSSW because the only used tool is ROOT.

Usage:
1. cmsenv #In order to use ROOT. 
2. Modify the directory of input root files in the shell script, exePreselection.
3. ./exePreselection #Perform preselection ~7hrs. Additional option: -d or --dryRun for testing purpose. 
4. ./doSelection #Perform selection ~4min.
5. root -l -b -q src/stackHist.C #Make stack plots.

(optional)
6 ./makePUweight.py -x 69200 #Produce rootfile of pileup reweighting factor ~5min
    - [Python module](https://github.com/cms-sw/cmssw/tree/master/SimGeneral/MixingModule/python?fbclid=IwAR2ehfE0hR8uaewPro4vQXos5I_IU6O7cyrtefQxTT6bMpyMETCTzpSuK58): 'mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi'
    - [Golden json](https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2017Analysis#13_TeV_pp_runs_ReReco): 'data/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
    - [Pileup json](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Pileup_JSON_Files_For_Run_II): '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt',
    #Need not execute when having data/MCPileUp.root
