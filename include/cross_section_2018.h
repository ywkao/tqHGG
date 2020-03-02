#ifndef __CROSS_SECTION_2018_H__
#define __CROSS_SECTION_2018_H__
#include <string>
using namespace std;
//NOTE: these are final cross section values with the branching fraction(s) taken into account.

float GetLuminosity_2018(){
    return 59.74;
}

float GetXsec_2018(char* dataset){
    if((string)dataset == "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 5765.4;
    else if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return 84;
    else if((string)dataset == "EGamma") return 1.;
    else if((string)dataset == "GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 232.9;
    else if((string)dataset == "GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8") return 3186;
    else if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 878.1;
    else if((string)dataset == "GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 0.11;
    else if((string)dataset == "QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 118100;

    else if((string)dataset == "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8") return 34.97;// to be updated
    else if((string)dataset == "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8") return 34.91;// to be updated

    else if((string)dataset == "TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8") return 2.967;
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

float GetBranchingFraction_2018(char* dataset){
    return 1.;
}

float GetTotalGenweight_2018(char* dataset){
    if((string)dataset == "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 3167968837617.244141;
    else if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return 6283025.499998;
    else if((string)dataset == "EGamma") return 1.000000;
    else if((string)dataset == "GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 13990283.000000;
    else if((string)dataset == "GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8") return 6072002.000000;
    else if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 9959511.000000;
    else if((string)dataset == "GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 213408144.617042;
    else if((string)dataset == "QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return 10766056.505585;
    else if((string)dataset == "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8") return 266470421.965493;
    else if((string)dataset == "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8") return 333129980.106007;
    else if((string)dataset == "TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8") return 5950566.862400;
    else if((string)dataset == "TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8") return 23938.678612;
    else if((string)dataset == "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8") return 32899050.301479;
    else if((string)dataset == "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 220589337894.843170;
    else if((string)dataset == "TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8") return 1690091.782056;
    else if((string)dataset == "TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8") return 3226764.905734;
    else if((string)dataset == "VBFHToGG_M125_13TeV_amcatnlo_pythia8") return 7681928.604484;
    else if((string)dataset == "VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 3800454.188507;
    else if((string)dataset == "WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 9954974046.222622;
    else if((string)dataset == "WW_TuneCP5_13TeV-pythia8") return 7846135.953322;
    else if((string)dataset == "WZ_TuneCP5_13TeV-pythia8") return 3884167.004738;
    else if((string)dataset == "ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8") return 1514115614.242440;
    else if((string)dataset == "ZZ_TuneCP5_13TeV-pythia8") return 1978776.751737;
    else if((string)dataset == "ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 263120.046941;
    else return -1.;
}
#endif

//ref: https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/MetaData/data/cross_sections.json
