#ifndef __CROSS_SECTION_H__
#define __CROSS_SECTION_H__
#include <string>
using namespace std;

double GetXsec(char* dataset){
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016B-07Aug17_ver2-v2") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016C-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016D-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016E-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016F-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016G-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016H-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016B-07Aug17_ver2-v2") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016C-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016D-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016E-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016F-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016G-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016H-07Aug17-v1") return 1.;
    if((string)dataset == "GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8") return 31.83;
    if((string)dataset == "VBFHToGG_M125_13TeV_amcatnlo_pythia8") return 3.858;
    if((string)dataset == "VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 2.162;
    if((string)dataset == "ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 0.4793;
    if((string)dataset == "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") return 5941.0;
    if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return 82.51;
    if((string)dataset == "GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8") return 219.2;
    if((string)dataset == "GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8") return 3255.0;
    if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8") return 862.4;
    if((string)dataset == "QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8") return 22180.0;
    if((string)dataset == "QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8") return 247000.0;
    if((string)dataset == "QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8") return 113100.0;
    if((string)dataset == "TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8") return 2.967;
    if((string)dataset == "TTGG_0Jets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8") return 0.01731;
    if((string)dataset == "TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8") return 3.795;
    if((string)dataset == "TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8") return 1.;//Not used in Yi-Jhih's analysis
    if((string)dataset == "WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") return 1.;//Not used in Yi-Jhih's analysis
    if((string)dataset == "WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8") return 1.;//Not used in Yi-Jhih's analysis
    if((string)dataset == "WW_TuneCUETP8M1_13TeV-pythia8") return 1.;//Not used in Yi-Jhih's analysis
    if((string)dataset == "WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8") return 1.;//Not used in Yi-Jhih's analysis
    if((string)dataset == "ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") return 1.;//Not used in Yi-Jhih's analysis
    if((string)dataset == "ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8") return 1.;//Not used in Yi-Jhih's analysis

    /*
    if((string)dataset == "DiPhotonJetsBox_M40_80-Sherpa") return 299.3;
    if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return 84.0;
    if((string)dataset == "DoubleEG_B") return 1.;
    if((string)dataset == "DoubleEG_C") return 1.;
    if((string)dataset == "DoubleEG_D") return 1.;
    if((string)dataset == "DoubleEG_E") return 1.;
    if((string)dataset == "DoubleEG_F") return 1.;
    if((string)dataset == "GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 218.61;
    if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 863;
    if((string)dataset == "GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8") return 3329.475;
    if((string)dataset == "QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 24810;
    if((string)dataset == "QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8") return 241400;
    if((string)dataset == "QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 118100;
    if((string)dataset == "TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8") return 2.967; //TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8
    if((string)dataset == "TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8") return 0.01687;
    if((string)dataset == "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8") return 3.697; //TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8
    if((string)dataset == "ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8") return 831.76;
    if((string)dataset == "ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8") return 831.76;
    if((string)dataset == "ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8") return 831.76;
    if((string)dataset == "ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8") return 831.76;
    if((string)dataset == "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8") return 831.76;
    if((string)dataset == "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8") return 831.76;
    if((string)dataset == "TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8") return 831.76;
    if((string)dataset == "TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8") return 831.76;
    if((string)dataset == "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8") return 831.76;
    if((string)dataset == "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8") return 831.76;
    if((string)dataset == "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8") return 831.76;
    if((string)dataset == "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8") return 831.76;
    if((string)dataset == "GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8") return 48.58;
    if((string)dataset == "VBFHToGG_M125_13TeV_amcatnlo_pythia8") return 3.7820;
    if((string)dataset == "VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 2.2569;
    if((string)dataset == "ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 0.5071;
    */
}

double GetBranchingFraction(char* dataset){
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016B-07Aug17_ver2-v2") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016C-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016D-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016E-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016F-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016G-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016H-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016B-07Aug17_ver2-v2") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016C-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016D-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016E-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016F-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016G-07Aug17-v1") return 1.;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016H-07Aug17-v1") return 1.;
    if((string)dataset == "GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8") return 0.00227;
    if((string)dataset == "VBFHToGG_M125_13TeV_amcatnlo_pythia8") return 0.00227;
    if((string)dataset == "VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 0.00227;
    if((string)dataset == "ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 0.00227;
    if((string)dataset == "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") return 1.;
    if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return 1.;
    if((string)dataset == "GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8") return 1.;
    if((string)dataset == "GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8") return 1.;
    if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8") return 1.;
    if((string)dataset == "QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8") return 1.;
    if((string)dataset == "QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8") return 1.;
    if((string)dataset == "QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8") return 1.;
    if((string)dataset == "TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8") return 1.;
    if((string)dataset == "TTGG_0Jets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8") return 1.;
    if((string)dataset == "TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8") return 1.;
    if((string)dataset == "TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8") return 1.;//Not used in Yi-Jhih's analysis
    if((string)dataset == "WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") return 1.;//Not used in Yi-Jhih's analysis
    if((string)dataset == "WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8") return 1.;//Not used in Yi-Jhih's analysis
    if((string)dataset == "WW_TuneCUETP8M1_13TeV-pythia8") return 1.;//Not used in Yi-Jhih's analysis
    if((string)dataset == "WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8") return 1.;//Not used in Yi-Jhih's analysis
    if((string)dataset == "ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") return 1.;//Not used in Yi-Jhih's analysis
    if((string)dataset == "ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8") return 1.;//Not used in Yi-Jhih's analysis


    /*
    if((string)dataset == "DiPhotonJetsBox_M40_80-Sherpa") return 1.;
    if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return 1.;
    if((string)dataset == "DoubleEG_B") return 1.;
    if((string)dataset == "DoubleEG_C") return 1.;
    if((string)dataset == "DoubleEG_D") return 1.;
    if((string)dataset == "DoubleEG_E") return 1.;
    if((string)dataset == "DoubleEG_F") return 1.;
    if((string)dataset == "GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 1.;
    if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 1.;
    if((string)dataset == "GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8") return 1.;
    if((string)dataset == "QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 1.;
    if((string)dataset == "QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8") return 1.;
    if((string)dataset == "QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 1.;
    if((string)dataset == "TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8") return 1.; //TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8
    if((string)dataset == "TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8") return 1.;
    if((string)dataset == "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8") return 1.; //TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8
    if((string)dataset == "ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8") return 0.00227;
    if((string)dataset == "ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8") return 0.00227;
    if((string)dataset == "ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8") return 0.00227;
    if((string)dataset == "ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8") return 0.00227;
    if((string)dataset == "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8") return 0.00227;
    if((string)dataset == "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8") return 0.00227;
    if((string)dataset == "TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8") return 0.00227;
    if((string)dataset == "TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8") return 0.00227;
    if((string)dataset == "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8") return 0.00227;
    if((string)dataset == "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8") return 0.00227;
    if((string)dataset == "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8") return 0.00227;
    if((string)dataset == "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8") return 0.00227;
    if((string)dataset == "GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8") return 0.00227;
    if((string)dataset == "VBFHToGG_M125_13TeV_amcatnlo_pythia8") return 0.00227;
    if((string)dataset == "VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 0.00227;
    if((string)dataset == "ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 0.00227;
    */
}
#endif

//ref: https://github.com/cms-analysis/flashgg/blob/master/MetaData/data/cross_sections.json
