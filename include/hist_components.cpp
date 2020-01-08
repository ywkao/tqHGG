#include "hist_components.h"

component_plots::component_plots()
{
}
component_plots::~component_plots()
{
    delete c1;
    delete legend;
    delete hist;
    delete hist_positive;
    delete hist_negative;
}

void component_plots::Init(const char* dirName, const char* histName, const char* ext_01, const char* ext_02, int nbins, double lower_bound, double upper_bound)
{
    sprintf(dirName_, "%s", dirName);
    sprintf(histName_, "%s", histName);
    c1 = new TCanvas("c1", "c1", 800, 600);
    legend = new TLegend(0.65,0.53,0.85,0.73);
    hist = new TH1D(histName, "", nbins, lower_bound, upper_bound);
    hist_positive = new TH1D(Form("%s_%s", histName, ext_01), "", nbins, lower_bound, upper_bound);
    hist_negative = new TH1D(Form("%s_%s", histName, ext_02), "", nbins, lower_bound, upper_bound);
}

void component_plots::Fill_hist(bool canFindSolution, double value)
{
    hist -> Fill(value);
    if(canFindSolution) hist_positive -> Fill(value);
    else                hist_negative -> Fill(value);
}

void component_plots::Draw_hist(const char* label1, const char* label2)
{
    hist -> Draw("hist");
    hist -> SetMinimum(0.);
    hist -> SetLineWidth(2);
    hist_positive -> Draw("hist;same");
    hist_positive -> SetLineColor(kRed);
    hist_positive -> SetLineWidth(2);
    hist_positive -> SetLineStyle(2);
    hist_negative -> Draw("hist;same");
    hist_negative -> SetLineColor(kGreen);
    hist_negative -> SetLineWidth(2);
    legend->Clear();
    legend->AddEntry(hist, "All", "l");
    legend->AddEntry(hist_positive, label1, "l");
    legend->AddEntry(hist_negative, label2, "l");
    legend->SetLineColor(0);
    legend->Draw("same");
    c1->SaveAs(Form("%s/%s.png", dirName_, histName_));
    //c1->SaveAs(Form("ntuples_skimmed/%s.png", histName));
}
void component_plots::Report_hist(const char* label1, const char* label2)
{
    int _n1, _n2;
    _n1 = hist_positive->GetEntries(); _n2 = hist->GetEntries();
    printf("[INFO] %s / all = %d / %d (%6.2f%%)\n", label1, _n1, _n2, 100. * (double)_n1 / (double)_n2 );
    _n1 = hist_negative->GetEntries(); _n2 = hist->GetEntries();
    printf("[INFO] %s / all = %d / %d (%6.2f%%)\n", label2, _n1, _n2, 100. * (double)_n1 / (double)_n2 );
}
