#ifndef __STACK_H__
#define __STACK_H__
#include <string>
using namespace std;

const int NUM_sig = 12;
const int NUM_resbkg = 4;
const int NUM_nonresbkg = 11;
const int NUM_data = 5;
const int NUM = NUM_sig + NUM_resbkg + NUM_nonresbkg + NUM_data;

void MakeStackHist(const char* histName);
void RegisterHistogram(const char* fileName, TH1D* &hist, const char* histName, int color, bool isSigMC = true, bool isData = false);
string GetXtitleAccordingToHistName(const char* histName);
string GetYtitleAccordingToHistName(const char* histName, double BinWidth);
bool isThisIDMVA(const char* histName);
bool isThisNumEtaPhi(const char* histName);
bool isThisMassSpectrum(const char* histName);
bool isThisDijetSpectrum(const char* histName);
bool isThisDiPhotonSpectrum(const char* histName);
bool isThisTopSpectrum(const char* histName);

//#--------------- Signal ---------------#
string fileNames_sig[NUM_sig] = {
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8/hist_TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8/hist_TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8/hist_TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8/hist_TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8/hist_TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8/hist_TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8/hist_TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8/hist_TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8/hist_ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8/hist_ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8/hist_ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8/hist_ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root"
};

//#--------------- Resonant bkg ---------------#
string fileNames_resbkg[NUM_resbkg] = {
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/hist_GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/VBFHToGG_M125_13TeV_amcatnlo_pythia8/hist_VBFHToGG_M125_13TeV_amcatnlo_pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/hist_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/hist_ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root"
};

//#--------------- non-Resonant bkg ---------------#
string fileNames_nonresbkg[NUM_nonresbkg] = {
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/DiPhotonJetsBox_M40_80-Sherpa/hist_DiPhotonJetsBox_M40_80-Sherpa.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa/hist_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/hist_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8/hist_GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/hist_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8/hist_TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8/hist_TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/hist_TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/hist_QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8/hist_QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/hist_QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root"
};

//#--------------- Data Samples ---------------#
string fileNames_data[NUM_data] = {
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/DoubleEG_B/hist_DoubleEG_B.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/DoubleEG_C/hist_DoubleEG_C.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/DoubleEG_D/hist_DoubleEG_D.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/DoubleEG_E/hist_DoubleEG_E.root",
"/home/ykao/legacy/CMSSW_9_4_10/src/tqHGG/plots/DoubleEG_F/hist_DoubleEG_F.root"
};
#endif
