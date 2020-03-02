#ifndef __STACK_H__
#define __STACK_H__
#include <string>
using namespace std;


//const char TARGET_DIR[128] = "plots_leptonic_test_star/161718/mva";
//const char SIGNAL_DIR[128] = "plots_leptonic_test_star/161718/mva";
const char TARGET_DIR[128] = "plots_leptonic_latest/161718/mva";
const char SIGNAL_DIR[128] = "plots_leptonic_latest/161718/mva";
//const char TARGET_DIR[128] = "plots_leptonic_update/161718/mva";
//const char SIGNAL_DIR[128] = "plots_leptonic_update/161718/mva";
//const char TARGET_DIR[128] = "collections/plots_leptonic_update/161718/app_stack";
//const char SIGNAL_DIR[128] = "collections/plots_leptonic_update/161718/app_stack";
//const char TARGET_DIR[128] = "collections/plots_leptonic_latest/161718/app_stack";
//const char SIGNAL_DIR[128] = "collections/plots_leptonic_latest/161718/app_stack";

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
void CalculateHistYields_mva_score(const char *process, TH1D* hist);
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

//TH1D* GetSignificanceHist(TH1D* sig, TH1D* bkg);
//void GetSignificanceHist(TH1D *& hist_significance, TH1D* sig, TH1D* bkg);
void GetSignificanceHist(TH1D *& hist_significance, TH1D *& hist_yield_s, TH1D *& hist_yield_b, TH1D* sig, TH1D* bkg);
double calculate_ZA(float s, float b);
double unc_ZA(float s, float b, float ds, float db);

//#--------------- Signal ---------------#
string fileNames_sig[NUM_sig] = {
Form("%s/output_score_TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/output_score_TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/output_score_TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/output_score_TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8.root", SIGNAL_DIR),
//-----
Form("%s/output_score_TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/output_score_TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/output_score_TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/output_score_TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8.root", SIGNAL_DIR),
//-----
Form("%s/output_score_ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/output_score_ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/output_score_ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root", SIGNAL_DIR),
Form("%s/output_score_ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root", SIGNAL_DIR)
};

//#--------------- Resonant bkg ---------------#
string fileNames_resbkg[NUM_resbkg] = {
Form("%s/output_score_GluGluHToGG.root", TARGET_DIR),
Form("%s/output_score_VBFHToGG.root", TARGET_DIR),
Form("%s/output_score_VHToGG.root", TARGET_DIR),
Form("%s/output_score_ttHJetToGG.root", TARGET_DIR)
};

//#--------------- non-Resonant bkg ---------------#
string fileNames_nonresbkg[NUM_nonresbkg] = {
Form("%s/output_score_DiPhotonJetsBox.root", TARGET_DIR),
Form("%s/output_score_GJet.root", TARGET_DIR),
Form("%s/output_score_TGJets.root", TARGET_DIR),
Form("%s/output_score_TTGG.root", TARGET_DIR),
Form("%s/output_score_TTGJets.root", TARGET_DIR),
Form("%s/output_score_QCD.root", TARGET_DIR),
Form("%s/output_score_TTJets.root", TARGET_DIR),
Form("%s/output_score_DYJetsToLL.root", TARGET_DIR),
Form("%s/output_score_ZGToLLG.root", TARGET_DIR),
Form("%s/output_score_WGToLNuG.root", TARGET_DIR),
Form("%s/output_score_WW.root", TARGET_DIR),
Form("%s/output_score_WZTo2L2Q.root", TARGET_DIR),
Form("%s/output_score_ZZTo2L2Q.root", TARGET_DIR)
};

//#--------------- Data Samples ---------------#
string fileNames_data[NUM_data] = {
Form("%s/output_score_Data.root", TARGET_DIR)
};
#endif
