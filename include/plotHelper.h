#ifndef _plotHelper_
#define _plotHelper_

#include <string>
#include <vector>
#include <TCanvas.h>
#include <TColor.h>
#include <TLegend.h>
#include <TPad.h>
#include <TH1D.h>
using namespace std;

double ph_getMaximumScope(vector<TH1D*> hists);
void ph_setMaximumScope(vector<TH1D*> hists);
void ph_setMaximumScope(vector<TH1D*> hists, double max);
void ph_setLegendPosition(TLegend *legend, double x1, double y1, double x2, double y2);
void ph_setPlot(TH1D* hist, string style, int color);
void ph_setBarPlot(TH1D* hist, string style, int color, double offset);
void ph_makePlot(TCanvas *c1, TH1D *hist, string name);
void ph_bjet_study(TCanvas *c1, TH1D *h_ref, TH1D *h_check, string name);
void ph_makePlots(TCanvas *c1, vector<TH1D*> hists, vector<string> styles, vector<int> colors, vector<int> orders, vector<string> legends, string xtitle, string name);

#endif
