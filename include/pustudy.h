#ifndef __PUSTUDY_H__
#define __PUSTUDY_H__
#include <TFile.h>
#include <TH1D.h>
#include <string>
using namespace std;

const int NUM_data = 14;

void RegisterHistogram(TFile *&file, const char* fileName, TH1D* &hist, const char* histName, int color);

string fileNames_data[NUM_data] = {
"ntuples_pustudy/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016B-07Aug17_ver2-v2.root",
"ntuples_pustudy/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016C-07Aug17-v1.root",
"ntuples_pustudy/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016D-07Aug17-v1.root",
"ntuples_pustudy/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016E-07Aug17-v1.root",
"ntuples_pustudy/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016F-07Aug17-v1.root",
"ntuples_pustudy/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016G-07Aug17-v1.root",
"ntuples_pustudy/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016H-07Aug17-v1.root",
"ntuples_pustudy/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016B-07Aug17_ver2-v2.root",
"ntuples_pustudy/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016C-07Aug17-v1.root",
"ntuples_pustudy/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016D-07Aug17-v1.root",
"ntuples_pustudy/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016E-07Aug17-v1.root",
"ntuples_pustudy/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016F-07Aug17-v1.root",
"ntuples_pustudy/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016G-07Aug17-v1.root",
"ntuples_pustudy/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016H-07Aug17-v1.root"
};
#endif
