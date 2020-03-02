#ifndef __STACK_H__
#define __STACK_H__
#include <string>
using namespace std;


const char TARGET_DIR[128] = "plots/161718";
const char SIGNAL_DIR[128] = "plots/2017old";

const char year[32] = "161718";

const int NUM_sig = 12;
const int NUM_resbkg = 4;
const int NUM_nonresbkg = 13;
const int NUM_data = 1;
const int NUM = NUM_sig + NUM_resbkg + NUM_nonresbkg + NUM_data;

void MakeStackHist(const char* histName);
void RegisterHistogram(TFile *&file, const char* fileName, TH1D* &hist, const char* histName, int color, bool isSigMC = true, bool isData = false);
void CalculateHistYields_signalRegion(const char *process, TH1D* hist);
void CalculateHistYields_sidebandRegion(const char *process, TH1D* hist);
void CalculateHistYields_fullRegion(const char *process, TH1D* hist);
double SumErrors(TH1D* hist);
double GetMaxScope(THStack* h1, std::vector<TH1D*> vec_hists);
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
Form("%s/TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8/hist_TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8/hist_TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8/hist_TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8/hist_TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8.root", SIGNAL_DIR),
//-----
Form("%s/TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8/hist_TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8/hist_TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8/hist_TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8/hist_TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8.root", SIGNAL_DIR),
//-----
Form("%s/ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8/hist_ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8/hist_ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8/hist_ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8/hist_ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root", SIGNAL_DIR)
};

//#--------------- Resonant bkg ---------------#
string fileNames_resbkg[NUM_resbkg] = {
Form("%s/GluGluHToGG/hist_GluGluHToGG.root", TARGET_DIR),
Form("%s/VBFHToGG/hist_VBFHToGG.root", TARGET_DIR),
Form("%s/VHToGG/hist_VHToGG.root", TARGET_DIR),
Form("%s/ttHJetToGG/hist_ttHJetToGG.root", TARGET_DIR)
};

//#--------------- non-Resonant bkg ---------------#
string fileNames_nonresbkg[NUM_nonresbkg] = {
Form("%s/DiPhotonJetsBox/hist_DiPhotonJetsBox.root", TARGET_DIR),
Form("%s/GJet/hist_GJet.root", TARGET_DIR),
Form("%s/TGJets/hist_TGJets.root", TARGET_DIR),
Form("%s/TTGG/hist_TTGG.root", TARGET_DIR),
Form("%s/TTGJets/hist_TTGJets.root", TARGET_DIR),
Form("%s/QCD/hist_QCD.root", TARGET_DIR),
Form("%s/TTJets/hist_TTJets.root", TARGET_DIR),
Form("%s/DYJetsToLL/hist_DYJetsToLL.root", TARGET_DIR),
Form("%s/ZGToLLG/hist_ZGToLLG.root", TARGET_DIR),
Form("%s/WGToLNuG/hist_WGToLNuG.root", TARGET_DIR),
Form("%s/WW/hist_WW.root", TARGET_DIR),
Form("%s/WZTo2L2Q/hist_WZTo2L2Q.root", TARGET_DIR),
Form("%s/ZZTo2L2Q/hist_ZZTo2L2Q.root", TARGET_DIR)
};

//#--------------- Data Samples ---------------#
string fileNames_data[NUM_data] = {
Form("%s/Data/hist_Data.root", TARGET_DIR)
};
#endif
