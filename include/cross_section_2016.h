#ifndef __CROSS_SECTION_2016_H__
#define __CROSS_SECTION_2016_H__
#include <string>
using namespace std;
//NOTE: these are final cross section values with the branching fraction(s) taken into account.

float GetLuminosity_2016(){
    return 35.92;
}

float GetXsec_2016(char* dataset){
    if((string)dataset == "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") return 5765.4;
    else if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return 84.0;
    else if((string)dataset == "DoubleEG_B") return 1.;
    else if((string)dataset == "DoubleEG_C") return 1.;
    else if((string)dataset == "DoubleEG_D") return 1.;
    else if((string)dataset == "DoubleEG_E") return 1.;
    else if((string)dataset == "DoubleEG_F") return 1.;
    else if((string)dataset == "DoubleEG_G") return 1.;
    else if((string)dataset == "DoubleEG_H") return 1.;
    else if((string)dataset == "GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8") return 232.9;
    else if((string)dataset == "GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8") return 3186.0;
    else if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8") return 878.1;
    else if((string)dataset == "GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8") return 0.11;
    else if((string)dataset == "QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8") return 24810.0;
    else if((string)dataset == "QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8") return 118100.0;

    else if((string)dataset == "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8") return 34.97;// to be updated
    else if((string)dataset == "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8") return 34.91;// to be updated
    else if((string)dataset == "TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8") return 2.967;// to be updated
    else if((string)dataset == "THQ_ctcvcp_HToGG_M125_13TeV-madgraph-pythia8") return 1.;// to be updated
    else if((string)dataset == "THW_ctcvcp_HToGG_M125_13TeV-madgraph-pythia8") return 1.;// to be updated

    else if((string)dataset == "TTGG_0Jets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8") return 0.01687;
    else if((string)dataset == "TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8") return 3.819;
    else if((string)dataset == "TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8") return 830;

    else if((string)dataset == "TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8") return 0.2149;// to be updated
    else if((string)dataset == "TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8") return 0.2432;// to be updated

    else if((string)dataset == "VBFHToGG_M125_13TeV_amcatnlo_pythia8") return 0.0086;
    else if((string)dataset == "WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") return 489.0;
    else if((string)dataset == "WW_TuneCUETP8M1_13TeV-pythia8") return 75.8;
    else if((string)dataset == "WZ_TuneCUETP8M1_13TeV-pythia8") return 5.52;
    else if((string)dataset == "ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") return 50.43;
    else if((string)dataset == "ZZ_TuneCUETP8M1_13TeV-pythia8") return 3.38;
    else if((string)dataset == "ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 0.0016;

    else return 1.;
}

float GetBranchingFraction_2016(char* dataset){
    return 1.;
}

float GetTotalGenweight_2016(char* dataset){
    if((string)dataset == "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") return 1400254241336.664062;
    else if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return 5.55474e+07;
    else if((string)dataset == "DoubleEG_B") return 1.000000;
    else if((string)dataset == "DoubleEG_C") return 1.000000;
    else if((string)dataset == "DoubleEG_D") return 1.000000;
    else if((string)dataset == "DoubleEG_E") return 1.000000;
    else if((string)dataset == "DoubleEG_F") return 1.000000;
    else if((string)dataset == "DoubleEG_G") return 1.000000;
    else if((string)dataset == "DoubleEG_H") return 1.000000;
    else if((string)dataset == "GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8") return 25030146.000000;
    else if((string)dataset == "GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8") return 1520381.000000;
    else if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8") return 7149943.000000;
    else if((string)dataset == "GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8") return 89212192.123418;
    else if((string)dataset == "QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8") return 16739745.504071;
    else if((string)dataset == "QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8") return 20816188.000000;
    else if((string)dataset == "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8") return 6858196.000000;
    else if((string)dataset == "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8") return 992024.000000;
    else if((string)dataset == "TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8") return 796307.473350;
    else if((string)dataset == "THQ_ctcvcp_HToGG_M125_13TeV-madgraph-pythia8") return 1980260.000000;
    else if((string)dataset == "THW_ctcvcp_HToGG_M125_13TeV-madgraph-pythia8") return 799260.000000;
    else if((string)dataset == "TTGG_0Jets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8") return 25098.827878;
    else if((string)dataset == "TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8") return 25990976.734347;
    else if((string)dataset == "TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8") return 54866099784.715492;
    else if((string)dataset == "TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8") return 1054697.948001;
    else if((string)dataset == "TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8") return 502198.859443;
    else if((string)dataset == "VBFHToGG_M125_13TeV_amcatnlo_pythia8") return 7355317.204163;
    else if((string)dataset == "WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") return 4697824264.060372;
    else if((string)dataset == "WW_TuneCUETP8M1_13TeV-pythia8") return 994032.181008;
    else if((string)dataset == "WZ_TuneCUETP8M1_13TeV-pythia8") return 1000000.000000;
    else if((string)dataset == "ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") return 3414602356.741301;
    else if((string)dataset == "ZZ_TuneCUETP8M1_13TeV-pythia8") return 965352.000000;
    else if((string)dataset == "ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8") return 198125.438385;
    else return -1.;
}
#endif

//ref: https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/MetaData/data/cross_sections.json

