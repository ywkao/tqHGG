#ifndef __STACK_H__
#define __STACK_H__
#include <string>
using namespace std;

const int NUM_sig = 12;
const int NUM_resbkg = 4;
const int NUM_nonresbkg = 10;//lack of DiPhotonJetsBox_M40to80_13TeV-Sherpa
const int NUM_data = 14;
const int NUM = NUM_sig + NUM_resbkg + NUM_nonresbkg + NUM_data;

void MakeStackHist(const char* histName);
//void RegisterHistogram(TFile *&file, const char* fileName, TH1D* &hist, const char* histName, int color, bool isSigMC = true, bool isData = false);
void RegisterHistogram(TFile *&file, const char* dataset, TH1D* &hist, const char* histName, int color, bool isSigMC = true, bool isData = false);
void CalculateHistYields(const char *process, TH1D* hist);
double SumErrors(TH1D* hist);

string GetXtitleAccordingToHistName(const char* histName);
string GetYtitleAccordingToHistName(const char* histName, double BinWidth);
bool isThisNVtxRho(const char* histName);
bool isThisIDMVA(const char* histName);
bool isThisNumEtaPhi(const char* histName);
bool isThisRho(const char* histName);
bool isThisMassSpectrum(const char* histName);
bool isThisDijetSpectrum(const char* histName);
bool isThisDiPhotonSpectrum(const char* histName);
bool isThisTopSpectrum(const char* histName);

//#--------------- Signal ---------------#
string dataset_sig[NUM_sig] = {
"TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8",
"TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8",
"TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8",
"TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8",
"TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8",
"TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8",
"TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8",
"TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8",
"ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8",
"ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8",
"ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8",
"ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8"
};

//#--------------- Resonant bkg ---------------#
string dataset_resbkg[NUM_resbkg] = {
"GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8",
"VBFHToGG_M125_13TeV_amcatnlo_pythia8",
"VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8",
"ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8"
};

//#--------------- non-Resonant bkg ---------------#
string dataset_nonresbkg[NUM_nonresbkg] = {
"DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa",
"GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8",
"GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8",
"GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8",
"QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8",
"QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8",
"QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8",
"TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8",
"TTGG_0Jets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8",
"TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8",
};

//#--------------- Data Samples ---------------#
string dataset_data[NUM_data] = {
"DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016B-07Aug17_ver2-v2",
"DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016C-07Aug17-v1",
"DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016D-07Aug17-v1",
"DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016E-07Aug17-v1",
"DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016F-07Aug17-v1",
"DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016G-07Aug17-v1",
"DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016H-07Aug17-v1",
"DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016B-07Aug17_ver2-v2",
"DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016C-07Aug17-v1",
"DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016D-07Aug17-v1",
"DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016E-07Aug17-v1",
"DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016F-07Aug17-v1",
"DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016G-07Aug17-v1",
"DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016H-07Aug17-v1"
};
#endif
