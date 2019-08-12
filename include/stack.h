#ifndef __STACK_H__
#define __STACK_H__
#include <string>
using namespace std;


const char TARGET_DIR[128] = "plots";
//const char TARGET_DIR[128] = "plots_leptonic/plots_noBtag";
//const char TARGET_DIR[128] = "plots_leptonic/plots_Btag_pu";
//const char TARGET_DIR[128] = "plots_leptonic/plots_Btag";
//
//const char TARGET_DIR[128] = "plots_hadronic/plots_noBtag";
//const char TARGET_DIR[128] = "plots_hadronic/plots_Btag_pu";
//const char TARGET_DIR[128] = "plots_hadronic/plots_Btag";

const int NUM_sig = 12;
const int NUM_resbkg = 4;
const int NUM_nonresbkg = 19;
const int NUM_data = 5;
const int NUM = NUM_sig + NUM_resbkg + NUM_nonresbkg + NUM_data;

void MakeStackHist(const char* histName);
void RegisterHistogram(TFile *&file, const char* fileName, TH1D* &hist, const char* histName, int color, bool isSigMC = true, bool isData = false);
void CalculateHistYields(const char *process, TH1D* hist);
double SumErrors(TH1D* hist);
string GetXtitleAccordingToHistName(const char* histName);
string GetYtitleAccordingToHistName(const char* histName, double BinWidth);
bool isThisIDMVA(const char* histName);
bool isThisNum(const char* histName);
bool isThisEtaPhi(const char* histName);
bool isThisMassSpectrum(const char* histName);
bool isThisDijetSpectrum(const char* histName);
bool isThisDiPhotonSpectrum(const char* histName);
bool isThisTopSpectrum(const char* histName);

//#--------------- Signal ---------------#
string fileNames_sig[NUM_sig] = {
Form("%s/TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8/hist_TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root", TARGET_DIR),
Form("%s/TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8/hist_TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8.root", TARGET_DIR),
Form("%s/TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8/hist_TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root", TARGET_DIR),
Form("%s/TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8/hist_TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8.root", TARGET_DIR),
//-----
Form("%s/TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8/hist_TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root", TARGET_DIR),
Form("%s/TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8/hist_TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8.root", TARGET_DIR),
Form("%s/TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8/hist_TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root", TARGET_DIR),
Form("%s/TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8/hist_TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8.root", TARGET_DIR),
//-----
Form("%s/ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8/hist_ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root", TARGET_DIR),
Form("%s/ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8/hist_ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root", TARGET_DIR),
Form("%s/ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8/hist_ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root", TARGET_DIR),
Form("%s/ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8/hist_ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root", TARGET_DIR)
};

//#--------------- Resonant bkg ---------------#
string fileNames_resbkg[NUM_resbkg] = {
Form("%s/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/hist_GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8.root", TARGET_DIR),
Form("%s/VBFHToGG_M125_13TeV_amcatnlo_pythia8/hist_VBFHToGG_M125_13TeV_amcatnlo_pythia8.root", TARGET_DIR),
Form("%s/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/hist_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root", TARGET_DIR),
Form("%s/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/hist_ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root", TARGET_DIR)
};

//#--------------- non-Resonant bkg ---------------#
string fileNames_nonresbkg[NUM_nonresbkg] = {
Form("%s/DiPhotonJetsBox_M40_80-Sherpa/hist_DiPhotonJetsBox_M40_80-Sherpa.root", TARGET_DIR),
Form("%s/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa/hist_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa.root", TARGET_DIR),
Form("%s/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/hist_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root", TARGET_DIR),
Form("%s/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8/hist_GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8.root", TARGET_DIR),
Form("%s/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/hist_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root", TARGET_DIR),
Form("%s/TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8/hist_TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8.root", TARGET_DIR),
Form("%s/TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8/hist_TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8.root", TARGET_DIR),
Form("%s/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/hist_TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root", TARGET_DIR),
Form("%s/QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/hist_QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root", TARGET_DIR),
Form("%s/QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8/hist_QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8.root", TARGET_DIR),
Form("%s/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/hist_QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root", TARGET_DIR),
Form("%s/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/hist_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8.root", TARGET_DIR),
Form("%s/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/hist_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8.root", TARGET_DIR),
Form("%s/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/hist_ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8.root", TARGET_DIR),
Form("%s/WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/hist_WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8.root", TARGET_DIR),
Form("%s/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/hist_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root", TARGET_DIR),
Form("%s/WW_TuneCP5_13TeV-pythia8/hist_WW_TuneCP5_13TeV-pythia8.root", TARGET_DIR),
Form("%s/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/hist_WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root", TARGET_DIR),
Form("%s/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/hist_ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root", TARGET_DIR)
};

//#--------------- Data Samples ---------------#
string fileNames_data[NUM_data] = {
Form("%s/DoubleEG_B/hist_DoubleEG_B.root", TARGET_DIR),
Form("%s/DoubleEG_C/hist_DoubleEG_C.root", TARGET_DIR),
Form("%s/DoubleEG_D/hist_DoubleEG_D.root", TARGET_DIR),
Form("%s/DoubleEG_E/hist_DoubleEG_E.root", TARGET_DIR),
Form("%s/DoubleEG_F/hist_DoubleEG_F.root", TARGET_DIR)
};
#endif
