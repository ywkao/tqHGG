//# include{{{
#include <stdio.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPolyLine.h>
#include <TPad.h>
#include "signalStrength.h"
using namespace std;
//}}}
const double UPPER_BOUND=1e+03, LOWER_BOUND=5e-02;
const double BF_min=0., BF_max=1.2;
const char xtitle[128] = "Branching fraction (t #rightarrow uH) [%]";
const char ytitle[128] = "95\% CL upper limit on signal strength";

void macro(const char* tag){
    //# canvas{{{
	TCanvas *can = new TCanvas("can","",800,600);
	can -> SetLogy();
	can -> SetLeftMargin(0.12);
	can -> SetBottomMargin(0.12);
    gPad->SetTicks(1,1);
    //gPad->SetGrid();
    //}}}
    //# plots{{{
    TH2D *frame = new TH2D("frame", Form(";%s;%s", xtitle, ytitle), 12, BF_min, BF_max, 50, LOWER_BOUND, UPPER_BOUND);
    TGraph *graph     = new TGraph(NUM, BF, Expected_500);
    //TGraph *graph_obs = new TGraph(NUM, BF, obs);
    TPolyLine *P1     = new TPolyLine(NUM*2);
    TPolyLine *P2     = new TPolyLine(NUM*2);

    for(int i=0; i<NUM*2; i++){
        if(i<NUM){
            P1->SetNextPoint(BF[i], Expected_160[i]);
            P2->SetNextPoint(BF[i], Expected_025[i]);
        } else{
            P1->SetNextPoint(BF[2*NUM-1-i], Expected_840[2*NUM-1-i]);
            P2->SetNextPoint(BF[2*NUM-1-i], Expected_975[2*NUM-1-i]);
        }
    }

	frame -> SetStats(0);
    frame -> GetXaxis() -> SetTitleOffset(1.1);
    frame -> GetXaxis() -> SetTitleFont(43);
    frame -> GetXaxis() -> SetTitleSize(25);
    frame -> GetYaxis() -> SetTitleOffset(1.1);
    frame -> GetYaxis() -> SetTitleFont(43);
    frame -> GetYaxis() -> SetTitleSize(25);
    P1->SetFillColor(kGreen);
    P2->SetFillColor(kYellow);
    P1->SetFillStyle(1001);
    P2->SetFillStyle(1001);
	graph -> SetLineColor(2);
	graph -> SetLineStyle(2);
	graph -> SetLineWidth(2);
	//graph_obs -> SetLineStyle(1);
	//graph_obs -> SetLineWidth(2);

	frame -> Draw();
    P2->Draw("f,same");
    P1->Draw("f,same");
    graph->Draw("pl,same");
    //graph_obs->Draw("pl,same");
    //}}}
    //# legend, text, line{{{
	TLegend *legend = new TLegend(0.55,0.60,0.80,0.85);
	//legend-> AddEntry(graph_obs," Observed Limit","l");
	legend->AddEntry(graph," Mediam expected","l");
	legend->AddEntry(P1," 68% expected","f");
	legend->AddEntry(P2," 95% expected","f");
	//legend->AddEntry(graph," Expected Limit","l");
	//legend->AddEntry(P1," #pm 1 #sigma Expected","f");
	//legend->AddEntry(P2," #pm 2 #sigma Expected","f");
	legend->SetLineColor(0);
	legend->SetFillStyle(0);
    legend->SetTextFont(45);
    legend->SetTextSize(25);
	legend->Draw("same");

	TLatex *text = new TLatex(0,0,"");
    text->SetNDC(kTRUE);
    text->SetTextFont(43);
    text->SetTextSize(25);
    text->SetTextAlign(11);
    text->DrawLatex(0.13, 0.92, "#bf{CMS} #it{Preliminary}");
    text->SetTextAlign(31);
	text->DrawLatex(0.90, 0.92,"137.2 fb^{-1} (13 TeV)");
    
    TLine *line = new TLine(0,0,0,0);
    line->SetLineColor(1);
    line->SetLineStyle(1);
    line->SetLineWidth(2);
    line->DrawLine(BF_min, 1, BF_max, 1);
    //}}}
    //# interpolation{{{
    int index=-999, counter=0;
    for(int i=0; i<NUM; ++i){
        if(Expected_500[i]<1.0){index=i; counter+=1;}
        if(counter!=0) break;
    }
    bool isExtrapolation = (index==-999);
    if(isExtrapolation) index = NUM - 1; // use the last two bins for extrapolation
    double branchingFraction_a = BF[index-1];
    double branchingFraction_b = BF[index];
    double signalStrength_a = Expected_500[index-1];
    double signalStrength_b = Expected_500[index];
    double coeff_a = (1.0000 - signalStrength_b) / (signalStrength_a - signalStrength_b);
    double coeff_b = (signalStrength_a - 1.0000) / (signalStrength_a - signalStrength_b);
    double signalStrength_upperlimit = coeff_a * branchingFraction_a + coeff_b * branchingFraction_b;
    if(!isExtrapolation) printf("[INFO] upperlimit on signal strength (interpolation) = %f\n", signalStrength_upperlimit);
    else                 printf("[INFO] upperlimit on signal strength (extrapolation) = %f\n", signalStrength_upperlimit);
    //------------------------------
    text->SetTextColor(kBlue);
    text->SetTextSize(30);
    text->SetTextAlign(22);
    text->DrawLatex(0.50, 0.20, Form("Upper limit on BF = %.3f %%", signalStrength_upperlimit) );
    //}}}
    can->SaveAs(Form("upperlimit_signalStrength_%s_lep.png", tag));
}
//# Usage::SetTextAlign{{{
//- align = 10*HorizontalAlign + VerticalAlign
//- 1=left adjusted, 2=centered, 3=right adjusted
//- 1=bottom adjusted, 2=centered, 3=top adjusted
//}}}
