#ifndef __MAIN_H__
#define __MAIN_H__

void MakePlots(TCanvas *c1, TH1D* hist, const char* title, const char* outputFile);
bool isThisDataOrNot(char* dataset);
bool isThisMultiFile(char* dataset);
bool isThisMCsignal(char* dataset);

#endif
