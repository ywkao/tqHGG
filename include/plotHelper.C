#include "plotHelper.h"
//void ph_setLegendPosition(TLegend *legend, double x1, double y1, double x2, double y2){{{
void ph_setLegendPosition(TLegend *legend, double x1, double y1, double x2, double y2)
{
    legend->SetX1NDC(x1);
    legend->SetX2NDC(x2);
    legend->SetY1NDC(y1);
    legend->SetY2NDC(y2);
    return;
}
//}}}

//double ph_getMaximumScope(vector<TH1D*> hists){{{
double ph_getMaximumScope(vector<TH1D*> hists)
{
    vector<double> vec;
    size_t total = hists.size();
    for(size_t i=0; i!=total; ++i)
        vec.push_back(hists[i]->GetMaximum());

    int max_index =  std::max_element(vec.begin(), vec.end()) - vec.begin();
    double max    = *std::max_element(vec.begin(), vec.end());
    return max;
}
//}}}
//void ph_setMaximumScope(vector<TH1D*> hists){{{
void ph_setMaximumScope(vector<TH1D*> hists)
{
    double scale = 1.2;
    double max = ph_getMaximumScope(hists);
    hists[0]->SetMaximum(max*scale);
    return;
}
//}}}
//void ph_setMaximumScope(vector<TH1D*> hists, double max){{{
void ph_setMaximumScope(vector<TH1D*> hists, double max)
{
    hists[0]->SetMaximum(max);
    return;
}
//}}}

//void ph_setBarPlot(TH1D* hist, string style, int color, double offset){{{
void ph_setBarPlot(TH1D* hist, string style, int color, double offset)
{
    hist->SetStats(0);
    hist->SetMinimum(0);
    hist->GetXaxis()->SetLabelSize(25);
    hist->GetXaxis()->SetLabelFont(43);
    hist->GetYaxis()->SetTitleSize(25);
    hist->GetYaxis()->SetTitleFont(43);
    hist->GetYaxis()->SetTitleOffset(1.);

    printf("ph_setBarPlot::color %d\n", color);

    if(style=="l"){
        hist->SetLineWidth(2);
        hist->SetLineColor(color);
    }

    if(style=="f"){
        hist->SetLineWidth(0);
        hist->SetFillColor(color);
        hist->SetFillStyle(3001);
    }

    hist->SetBarWidth(0.3);
    hist->SetBarOffset(offset);

    return;
}
//}}}
//void ph_setPlot(TH1D* hist, string style, int color){{{
void ph_setPlot(TH1D* hist, string style, int color)
{
    hist->SetStats(0);

    if(style=="l"){
        hist->SetLineWidth(2);
        hist->SetLineColor(color);
    }

    if(style=="f"){
        hist->SetLineWidth(0);
        hist->SetFillColor(color);
        hist->SetFillStyle(3001);
    }

    return;
}
//}}}

//void ph_makePlot(TCanvas *c1, TH1D *hist, string name){{{
void ph_makePlot(TCanvas *c1, TH1D *hist, string name)
{
    hist->SetLineWidth(2);
    hist->Draw();
    c1->SaveAs(name.c_str());
}
//}}}
//void ph_bjet_study(TCanvas *c1, TH1D *h_ref, TH1D *h_check, string name){{{
void ph_bjet_study(TCanvas *c1, TH1D *h_ref, TH1D *h_check, string name)
{
    gPad->SetTicks(0,1);
    h_ref->SetStats(0);
    h_ref->SetMaximum(1.);
    h_ref->SetMinimum(0.);
    h_ref->SetLineWidth(2);
    h_ref->SetLineColor(kBlue);
    h_ref->GetXaxis()->SetLabelSize(25);
    h_ref->GetXaxis()->SetLabelFont(43);
    h_ref->GetYaxis()->SetTitleSize(25);
    h_ref->GetYaxis()->SetTitleFont(43);
    h_ref->Draw();

    h_check->SetLineWidth(2);
    h_check->SetLineColor(2);
    h_check->Draw("same");

    TLegend *legend = new TLegend(0.50,0.20,0.88,0.45);
    legend->SetTextFont(43);
    legend->AddEntry(h_ref,    "the bjet is included", "l");
    legend->AddEntry(h_check,  "the bjet is the leading one", "l");
    legend->SetLineColor(0);
    legend->Draw("same");

    c1->SaveAs(name.c_str());
}
//}}}
//void ph_makePlots(TCanvas *c1, vector<TH1D*> hists, vector<string> styles, vector<int> colors, vector<int> orders, vector<string> legends, string xtitle, string name){{{
void ph_makePlots(TCanvas *c1, vector<TH1D*> hists, vector<string> styles, vector<int> colors, vector<int> orders, vector<string> legends, string xtitle, string name)
{
    size_t total = hists.size();

    // set plots
    double max = ph_getMaximumScope(hists);
    hists[0] -> SetMaximum(max*1.2);
    hists[0] -> GetXaxis() -> SetTitle(xtitle.c_str());
    for(size_t i=0; i!=total; ++i)
        ph_setPlot(hists[i], styles[i], colors[i]);

    // draw plots
    for(size_t i=0; i!=total; ++i)
    {
        if(i==0) hists[orders[0]] -> Draw();
        else     hists[orders[i]] -> Draw("same");
    }
    hists[orders[0]] -> Draw("same"); // on top of others

    // set legends
    TLegend *legend = new TLegend(0.50,0.55,0.85,0.85);
    legend->SetLineColor(0);
    //legend->SetTextSize(20);
    for(size_t i=0; i!=total; ++i)
        legend->AddEntry(hists[i], legends[i].c_str(), styles[i].c_str());
    legend->Draw("same");

    c1->SaveAs(name.c_str());
}
//}}}
