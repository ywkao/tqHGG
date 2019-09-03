#ifndef __CROSS_SECTION_H__
#define __CROSS_SECTION_H__
#include <string>
using namespace std;
//NOTE: these are final cross section values with the branching fraction(s) taken into account.

double GetXsec(char* dataset){
    if((string)dataset == "DiPhotonJetsBox_M40_80-Sherpa") return 299.3;
    else if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return 84.0;
    else if((string)dataset == "DoubleEG_B") return 1.;
    else if((string)dataset == "DoubleEG_C") return 1.;
    else if((string)dataset == "DoubleEG_D") return 1.;
    else if((string)dataset == "DoubleEG_E") return 1.;
    else if((string)dataset == "DoubleEG_F") return 1.;
    else if((string)dataset == "GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 232.9;
    else if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 878.1;
    else if((string)dataset == "GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8") return 3186.0;
    else if((string)dataset == "QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 24810.0;
    else if((string)dataset == "QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8") return 241400.0;
    else if((string)dataset == "QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 118100.0;
    else if((string)dataset == "TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8") return 2.967; //TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8
    else if((string)dataset == "TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8") return 0.01687;
    else if((string)dataset == "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8") return 3.819; //TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8
    else if((string)dataset == "ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8") return 0.067;
    else if((string)dataset == "ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8") return 0.139;
    else if((string)dataset == "ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8") return 0.009;
    else if((string)dataset == "ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8") return 0.019;
    else if((string)dataset == "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8") return 0.230;
    else if((string)dataset == "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8") return 0.230;
    else if((string)dataset == "TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8") return 0.110;
    else if((string)dataset == "TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8") return 0.110;
    else if((string)dataset == "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8") return 0.230;
    else if((string)dataset == "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8") return 0.230;
    else if((string)dataset == "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8") return 0.110;
    else if((string)dataset == "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8") return 0.110;
    else if((string)dataset == "GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8") return 0.11;
    else if((string)dataset == "VBFHToGG_M125_13TeV_amcatnlo_pythia8") return 0.0086;
    else if((string)dataset == "VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 0.00512; 
    else if((string)dataset == "ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 0.0016;

    //NOTE: not fully checked!!
    else if((string)dataset == "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 830;
    else if((string)dataset == "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 5765.4;
    else if((string)dataset == "WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 489.0;//Not found yet //WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8
    //else if((string)dataset == "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8") return 52940.;//DAS
    else if((string)dataset == "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8") return 55690.;//GenXSecAnalyzer
    else if((string)dataset == "WW_TuneCP5_13TeV-pythia8") return 75.8;//DAS
    else if((string)dataset == "WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8") return 5.52;
    else if((string)dataset == "ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8") return 3.38;
    else if((string)dataset == "ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 50.43;

    else return 1.;
}

double GetBranchingFraction(char* dataset){
    if((string)dataset == "DiPhotonJetsBox_M40_80-Sherpa") return 1.;
    else if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return 1.;
    else if((string)dataset == "DoubleEG_B") return 1.;
    else if((string)dataset == "DoubleEG_C") return 1.;
    else if((string)dataset == "DoubleEG_D") return 1.;
    else if((string)dataset == "DoubleEG_E") return 1.;
    else if((string)dataset == "DoubleEG_F") return 1.;
    else if((string)dataset == "GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 1.;
    else if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 1.;
    else if((string)dataset == "GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8") return 1.;
    else if((string)dataset == "QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 1.;
    else if((string)dataset == "QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8") return 1.;
    else if((string)dataset == "QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 1.;
    else if((string)dataset == "TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8") return 1.; //TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8
    else if((string)dataset == "TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8") return 1.;
    else if((string)dataset == "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8") return 1.; //TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8
    else if((string)dataset == "ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8") return 1.0;
    else if((string)dataset == "ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8") return 1.0;
    else if((string)dataset == "ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8") return 1.0;
    else if((string)dataset == "ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8") return 1.0;
    else if((string)dataset == "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8") return 1.0;
    else if((string)dataset == "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8") return 1.0;
    else if((string)dataset == "TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8") return 1.0;
    else if((string)dataset == "TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8") return 1.0;
    else if((string)dataset == "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8") return 1.0;
    else if((string)dataset == "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8") return 1.0;
    else if((string)dataset == "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8") return 1.0;
    else if((string)dataset == "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8") return 1.0;
    else if((string)dataset == "GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8") return 1.0;
    else if((string)dataset == "VBFHToGG_M125_13TeV_amcatnlo_pythia8") return 1.0;
    else if((string)dataset == "VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 1.0;
    else if((string)dataset == "ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 1.0;

    //NOTE: not fully checked!!
    else if((string)dataset == "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 1.;
    else if((string)dataset == "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 1.;
    else if((string)dataset == "WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 1.;
    else if((string)dataset == "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8") return 1.;
    else if((string)dataset == "WW_TuneCP5_13TeV-pythia8") return 1.;
    else if((string)dataset == "WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8") return 1.;
    else if((string)dataset == "ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8") return 1.;
    else if((string)dataset == "ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 1.;

    else return 1.;
}
#endif

//ref: https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/MetaData/data/cross_sections.json
