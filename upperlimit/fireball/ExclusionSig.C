#include <stdio.h>
#include <TCanvas.h>
#include "TGaxis.h"
#include "TGraph.h"
#include <TH2D.h>
#include "TLatex.h"
#include "TLegend.h"
#include <TPolyLine.h>
#include <TPad.h>
#include "signalStrength.h"
using namespace std;

const int UP=1, LOW=0;
//const double UPPER_BOUND=2e-02, LOWER_BOUND=-1e-03;
const double UPPER_BOUND=1e+03, LOWER_BOUND=5e-02;

void ExclusionSig(){
	TCanvas *can = new TCanvas("can","",800,600);
	can -> SetLogy();
	can -> SetLeftMargin(0.12);
	can -> SetBottomMargin(0.12);

    //MASS, Xsec, Err_Xsec, signal strength
    //#include "signalStrength.h"

    TH2D *frame = new TH2D("frame",";M_{Q} [TeV];95\%  CL  signal strength (Q#bar{Q} #rightarrow nV)",10, MASS_min, MASS_max, 50, LOWER_BOUND, UPPER_BOUND);
    TGraph *graph     = new TGraph(NUM, MASS, mean);
    TGraph *graph_obs = new TGraph(NUM, MASS, obs);
    TPolyLine *P1     = new TPolyLine(NUM*2);
    TPolyLine *P2     = new TPolyLine(NUM*2);

    for(int i=0; i<NUM*2; i++){
        if(i<NUM){
            P1->SetNextPoint(MASS[i], s1[LOW][i]);
            P2->SetNextPoint(MASS[i], s2[LOW][i]);
        } else{
            P1->SetNextPoint(MASS[2*NUM-1-i], s1[UP][2*NUM-1-i]);
            P2->SetNextPoint(MASS[2*NUM-1-i], s2[UP][2*NUM-1-i]);
        }
    }

	frame -> SetStats(0);
    frame -> GetYaxis() -> SetTitleOffset(1.5);
    P1->SetFillColor(kGreen);
    P2->SetFillColor(kYellow);
    P1->SetFillStyle(1001);
    P2->SetFillStyle(1001);
	graph -> SetLineStyle(2);
	graph -> SetLineWidth(2);
	graph_obs -> SetLineStyle(1);
	graph_obs -> SetLineWidth(2);
	//graph -> SetMarkerStyle(20);
	//graph -> SetMarkerColor(1);

	frame -> Draw();
    P2->Draw("f,same");
    P1->Draw("f,same");
    graph->Draw("pl,same");
    graph_obs->Draw("pl,same");

	TLegend *legend = new TLegend(0.20,0.66,0.45,0.87);
	legend-> AddEntry(graph_obs," Observed Limit","l");
	legend-> AddEntry(graph," Expected Limit","l");
	legend-> AddEntry(P1," #pm 1 #sigma Expected","f");
	legend-> AddEntry(P2," #pm 2 #sigma Expected","f");
	legend-> SetLineColor(0);
	legend-> SetFillStyle(0);
	legend-> Draw("same");

	TLatex *text = new TLatex(0,0,"");
    text->SetTextFont(43);
    text->SetTextSize(18);
	text->DrawLatex(0.425,1.2e+03,"MC Simulation from MadGraph5");
	text->DrawLatex(1.750,1.2e+03,"L = 100 fb^{-1} (13TeV)");
    //text->SetTextSize(25);
	//text->DrawLatex(0.9,5e+03,"#int Ldt = 100 fb^{-1}");
	//text->DrawLatex(0.9,5e+02,"#sqrt{s} = 13 TeV");
    
    TLine *line = new TLine(0,0,0,0);
    line->SetLineStyle(1);
    line->SetLineColor(1);
    line->DrawLine(MASS_min,1,MASS_max,1);

    can->SaveAs("UpperLimitSignalStrength.png");
}
