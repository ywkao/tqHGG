#ifndef __MAIN_H__
#define __MAIN_H__

double Chi2_calculator(double w_mass, double t_mass);
double Chi2_calculator_w_only(double w_mass);
void MakePlots(TCanvas *c1, TH1D* hist, const char* title, const char* outputFile);
bool isThisDataOrNot(char* dataset);
bool isThisMultiFile(char* dataset);
bool isThisMCsignal(char* dataset);

#endif
