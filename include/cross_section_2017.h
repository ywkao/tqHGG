#ifndef __CROSS_SECTION_2017_H__
#define __CROSS_SECTION_2017_H__
#include <string>
using namespace std;
//NOTE: these are final cross section values with the branching fraction(s) taken into account.

float GetLuminosity_2017(){
    return 41.53;
}

float GetXsec_2017(char* dataset){
    if((string)dataset == "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 5765.4;
    else if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return 84.0;
    else if((string)dataset == "DoubleEG") return 1.;
    else if((string)dataset == "GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 232.9;
    else if((string)dataset == "GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8") return 3186.0;
    else if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 878.1;
    else if((string)dataset == "GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8") return 0.11;
    else if((string)dataset == "QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 24810;
    else if((string)dataset == "QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8") return 241400;
    else if((string)dataset == "QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 118100;

    else if((string)dataset == "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8") return 34.97;// to be updated
    else if((string)dataset == "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8") return 34.91;// to be updated
    else if((string)dataset == "TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8") return 2.967;// to be updated
    else if((string)dataset == "THQ_ctcvcp_HToGG_M125_13TeV-madgraph-pythia8_TuneCP5") return 1.;// to be updated
    else if((string)dataset == "THW_ctcvcp_HToGG_M125_13TeV-madgraph-pythia8_TuneCP5") return 1.;// to be updated

    else if((string)dataset == "TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8") return 0.01687;
    else if((string)dataset == "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8") return 3.819;
    else if((string)dataset == "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 830;

    else if((string)dataset == "TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8") return 0.2149;// to be updated
    else if((string)dataset == "TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8") return 0.2432;// to be updated

    else if((string)dataset == "VBFHToGG_M125_13TeV_amcatnlo_pythia8") return 0.0086;
    else if((string)dataset == "VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 0.00512;
    else if((string)dataset == "WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 489;
    else if((string)dataset == "WW_TuneCP5_13TeV-pythia8") return 75.8;
    else if((string)dataset == "WZ_TuneCP5_13TeV-pythia8") return 5.52;
    else if((string)dataset == "ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 50.43;
    else if((string)dataset == "ZZ_TuneCP5_13TeV-pythia8") return 3.38;
    else if((string)dataset == "ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 0.0016;

    else return 1.;
}

float GetBranchingFraction_2017(char* dataset){
    return 1.;
}

float GetTotalGenweight_2017(char* dataset){
    if((string)dataset == "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 3156349621435.696777;
    else if((string)dataset == "DoubleEG") return 1.000000;
    else if((string)dataset == "GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 35284140.000000;
    else if((string)dataset == "GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8") return 40198250.000000;
    else if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 158315792.000000;
    else if((string)dataset == "GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8") return 219403548.242235;
    else if((string)dataset == "QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 14594575.000000;
    else if((string)dataset == "QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8") return 41974430.000000;
    else if((string)dataset == "QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 20277088.000000;
    else if((string)dataset == "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8") return 272182819.816308;
    else if((string)dataset == "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8") return 222419374.928597;
    else if((string)dataset == "TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8") return 5919452.569381;
    else if((string)dataset == "THQ_ctcvcp_HToGG_M125_13TeV-madgraph-pythia8_TuneCP5") return 1887041.000000;
    else if((string)dataset == "THW_ctcvcp_HToGG_M125_13TeV-madgraph-pythia8_TuneCP5") return 866503.000000;
    else if((string)dataset == "TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8") return 25145.499084;
    else if((string)dataset == "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8") return 66911664.249165;
    else if((string)dataset == "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 315906974408.811157;
    else if((string)dataset == "TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8") return 1644338.821542;
    else if((string)dataset == "TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8") return 1512335.496318;
    else if((string)dataset == "VBFHToGG_M125_13TeV_amcatnlo_pythia8") return 7715490.292197;
    else if((string)dataset == "VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 4100171.411237;
    else if((string)dataset == "WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 9250613153.254873;
    else if((string)dataset == "WW_TuneCP5_13TeV-pythia8") return 7791560.886900;
    else if((string)dataset == "WZ_TuneCP5_13TeV-pythia8") return 3928630.000000;
    else if((string)dataset == "ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 3305019647.078636;
    else if((string)dataset == "ZZ_TuneCP5_13TeV-pythia8") return 1925931.000000;
    else if((string)dataset == "ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 348168.008315;
    else if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return 21984000.000000;
    else return -1.;
}
#endif

//ref: https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/MetaData/data/cross_sections.json
