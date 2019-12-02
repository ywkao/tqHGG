# tqHGG

Searches for top FCNH with H decays to diphoton analysis. (2017)
Analysis code for 2016 is in branch named after "2016".

Although the analysis code is performed under CMSSW_9_4_10,
it should be independent of the version of CMSSW because the only used tool is ROOT.

The main computing resouces are the nodes on ntugrid5.

## Usage
All the actions are written in execution.sh :
```javascript=69
#------------------------- Main Exe Section -------------------------#
#./script/prepareExeForNewMC_npu_float.sh "preselection"
#Preselection
#Intermission

#./script/prepareExeForNewMC_npu_float.sh "selection"
#Selection "hadronic" #selection and make plots for hadronic channel
#AfterSelection "hadronic" 
#Selection "leptonic" #selection and make plots for leptonic channel
#AfterSelection "leptonic" 
#./check_errors.sh

ReRunStackPlotsOnly "hadronic"
#ReRunStackPlotsOnly "leptonic"
```

One can uncomment the commands and execute it.
```javascript=
$ cmsenv # in order to use ROOT
$ ./execution.sh
```

<p><b>Steps</b></p>
1. Preselection stage (~ 3 hrs)
    * Purpose: reduce information of the ntuples from the flashgg package
    * Action: uncomment the lines in exection.sh:
```javascript=70
./script/prepareExeForNewMC_npu_float.sh "preselection"
Preselection
Intermission
```

2. Selection stage (~ 5 mins)
    * Purpose: apply selection condition & top reconstruction & stack plots
    * Action: uncomment the lines in exection.sh:
```javascript=74
./script/prepareExeForNewMC_npu_float.sh "selection"
Selection "hadronic" #selection and make plots for hadronic channel
AfterSelection "hadronic" 
Selection "leptonic" #selection and make plots for leptonic channel
AfterSelection "leptonic" 
./check_errors.sh
```

3. Re-make the stack plots (~ 2 mins)
    * Purpose: performed only when modifying stack plots ( src/stackHist.C )
    * Action: uncomment the lines in exection.sh:
```javascript=81
ReRunStackPlotsOnly "hadronic"
ReRunStackPlotsOnly "leptonic"
```


## Pileup related
The pu reweighting factor can be calculated by a python script.
```javascript
$ ./makePUweight.py -x 69200 #Produce rootfile of pileup reweighting factor (~ 5 mins)
```

This step is optional. Need not execute once having its output file, data/MCPileUp.root<br />
[Python module](https://github.com/cms-sw/cmssw/tree/master/SimGeneral/MixingModule/python?fbclid=IwAR2ehfE0hR8uaewPro4vQXos5I_IU6O7cyrtefQxTT6bMpyMETCTzpSuK58): 'mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi' <br />
[Golden json](https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2017Analysis#13_TeV_pp_runs_ReReco): 'data/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt' <br />
[Pileup json](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Pileup_JSON_Files_For_Run_II): '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt' <br />


