#!/usr/bin/env python2
command='"cd /wk_cms2/ykao/CMSSW_9_4_10/src/2017/tqHGG; ./script/mvaTasks/mvaApplication.sh {}"'
import os

#signal = "TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8"
#os.system('./script/mvaTasks/mvaApp_edition.sh {0} {1}'.format(signal, "signal_dir"))
#os.system('./script/submitJOB.py --command={} --name={}'.format(command.format(signal), "mvapp_signal"))

#list = [ "DYJetsToLL" ]

# dataset{{{
list_signal = [
"ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8", \
"ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8", \
"ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8", \
"ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8", \
"TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8", \
"TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8", \
"TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8", \
"TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8", \
"TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8", \
"TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8", \
"TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8", \
"TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8"
]

list = [
"DYJetsToLL", \
"Data", \
"DiPhotonJetsBox", \
"GJet", \
"GluGluHToGG", \
"QCD", \
"TGJets", \
"TTGG", \
"TTGJets", \
"TTJets", \
"VBFHToGG", \
"VHToGG", \
"WGToLNuG", \
"WW", \
"WZTo2L2Q", \
"ZGToLLG", \
"ZZTo2L2Q", \
"ttHJetToGG"
]
#}}}

for dataset in list:
    tagName = 'mvapp_{}'.format(dataset)
    print tagName
    #--------------------------------------------------
    # 1) edit src/TMVAClassificationApplication.C
    # 2) cp src/TMVAClassificationApplication.C src/tmp/app_${dataset}/TMVAClassificationApplication.C
    os.system('./script/mvaTasks/mvaApp_edition.sh {0} {1}'.format(dataset, "Directory"))
    #--------------------------------------------------
    # submit jobs; root -l -b -q src/tmp/app_${dataset}/TMVAClassificationApplication.C
    os.system('./script/submitJOB.py --command={} --name={}'.format(command.format(dataset), tagName))

counter = 0
for dataset in list_signal:
    counter = counter + 1
    tagName = 'mvapp_sig_{}'.format(counter)
    print tagName
    #--------------------------------------------------
    # 1) edit src/TMVAClassificationApplication.C
    # 2) cp src/TMVAClassificationApplication.C src/tmp/app_${dataset}/TMVAClassificationApplication.C
    os.system('./script/mvaTasks/mvaApp_edition.sh {0} {1}'.format(dataset, "signal_dir"))
    #--------------------------------------------------
    # submit jobs; root -l -b -q src/tmp/app_${dataset}/TMVAClassificationApplication.C
    os.system('./script/submitJOB.py --command={} --name={}'.format(command.format(dataset), tagName))
