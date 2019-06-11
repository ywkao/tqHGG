#ifndef __PUSTUDY_H__
#define __PUSTUDY_H__
#include <TFile.h>
#include <TH1D.h>
#include <string>
using namespace std;

const int NUM_data = 14;
const char* dir = "collections/ntuples_pustudy";

void LoadFile(TFile *&file, const char* fileName);
void RegisterHistogram(TFile *file, TH1D* &hist, const char* histName, int color);

string fileNames_data[NUM_data] = {
Form("%s/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016B-07Aug17_ver2-v2.root", dir),
Form("%s/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016C-07Aug17-v1.root", dir),
Form("%s/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016D-07Aug17-v1.root", dir),
Form("%s/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016E-07Aug17-v1.root", dir),
Form("%s/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016F-07Aug17-v1.root", dir),
Form("%s/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016G-07Aug17-v1.root", dir),
Form("%s/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016H-07Aug17-v1.root", dir),
Form("%s/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016B-07Aug17_ver2-v2.root", dir),
Form("%s/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016C-07Aug17-v1.root", dir),
Form("%s/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016D-07Aug17-v1.root", dir),
Form("%s/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016E-07Aug17-v1.root", dir),
Form("%s/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016F-07Aug17-v1.root", dir),
Form("%s/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016G-07Aug17-v1.root", dir),
Form("%s/ntuple_DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016H-07Aug17-v1.root", dir)
};
const char* GetXtitle(const char* histName){
    if((string)histName == "hist_Rho")     return "Rho";
    if((string)histName == "hist_Rho_pu")  return "Rho";
    if((string)histName == "hist_NVtx")    return "Number of vertices";
    if((string)histName == "hist_NVtx_pu") return "Number of vertices";
    return "";
}
#endif
