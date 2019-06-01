//**************************************************
//
// Study of pile up reweighting.
// Compare shape of NVtx and Rho between data & MC
// Input: rootfiles coming from study_pileup_reweighting.cpp
// Output: png files (NVtx, Rho)
//
// ### NOTE: The plots are buggy due to the way processing histograms (pointers reset issue...?
//
//**************************************************
#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include "../include/putools.h"
#include "../include/pustudy.h"
void process(TH1D *hist, TH1D *h_data, const char *output);
void interface(const char* histName, const char* output);

TCanvas *c1 = new TCanvas("c1", "c1", 700, 800);

void study_pileup_reweighting_plots(){
    //interface("hist_Rho","plots/puReweighting_Rho.png");
    //interface("hist_Rho_pu","plots/puReweighting_Rho_pu.png");
    //interface("hist_NVtx","plots/puReweighting_vertices.png");
    interface("hist_NVtx_pu","plots/puReweighting_vertices_pu.png");
}

void interface(const char* histName, const char* output){
    TFile *f_data[NUM_data];
    TH1D  *h_data[NUM_data+1], *h;
    //--------------------------------------------------------------------------------
    //prepare data
    //--------------------------------------------------------------------------------
    for(int i=0; i<NUM_data; i++){
        RegisterHistogram(f_data[i], fileNames_data[i].c_str(), h_data[i], histName, kBlack);
    }

    h_data[NUM_data] = (TH1D*) h_data[0]->Clone();
    for(int i=1; i<NUM_data; i++) h_data[NUM_data]->Add(h_data[i]);
    //--------------------------------------------------------------------------------
    //prepare mc
    //--------------------------------------------------------------------------------
    TFile *f_mc = TFile::Open("ntuples_pustudy/ntuple_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa.root");
    h = (TH1D*) f_mc->Get(histName);
    //--------------------------------------------------------------------------------
    //make plot
    //--------------------------------------------------------------------------------
    process(h, h_data[NUM_data], output);
    //--------------------------------------------------------------------------------
    //clear ## Note: Do not close files, otherwise, no subsequent plots...
    //--------------------------------------------------------------------------------
    //f_mc->Close();
    //for(int i=0; i<NUM_data; i++) f_data[i]->Close();
}


void process(TH1D *hist, TH1D *h_data, const char *output){
    //TH1D *h_mc = new TH1D("h_mc",";;", 75, 0, 75); h_mc->Sumw2();
    TH1D *h_mc = (TH1D*)hist->Clone();

    //Rebin(hist, h_mc);

    double n_mc, n_data;
    n_mc = SumContents(h_mc);
    n_data = SumContents(h_data);

    Normalization(h_mc, n_mc);
    Normalization(h_data, n_data);

    TH1D *h_ratio = (TH1D*)h_data->Clone();
    h_ratio->Divide(h_mc);

    for(int i=0; i<100; ++i){
        double mc = h_mc->GetBinContent(i+1);
        double data = h_data->GetBinContent(i+1);
        double ratio = data/mc;
        double content = h_ratio->GetBinContent(i+1);
        printf("bin = %d, mc = %.5f, data = %.5f, ratio = %.5f, content = %5f\n", i+1, mc, data, ratio, content);
    }

    /*
    h_data->SetLineColor(kRed);
    h_data->SetMarkerStyle(20);
    h_data->Draw("p, E1");
    h_mc->SetLineColor(kBlue);
    h_mc->Draw("hist,same");
    c1->SaveAs(output);
    */

    //=======================//
    //=====  Make Plot  =====//
    //=======================//
    //--------------------------------------------------------------------------------
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
    pad1->SetBottomMargin(0); //Upper and lower pads are joined
    pad1->SetAttLinePS(kBlack,1,2);
    pad1->Draw();
    pad1->cd(); //pad1 becomes current pad
    //--------------------------------------------------------------------------------
    h_data->SetTitle("");
    h_data->SetStats(0);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(1.2);
    h_data->SetLineColor(kBlack);
    h_data->Draw("p,E1");
    h_mc->SetFillColor(kRed);
    h_mc->SetFillStyle(3001);
    h_mc->SetLineColor(kBlack);
    h_mc->SetLineWidth(2);
    h_mc->Draw("hist,same");
    h_data->Draw("p,E1,same");
    //--------------------------------------------------------------------------------
    TLegend *legend = new TLegend(0.60,0.65,0.85,0.85);
    legend->AddEntry(h_data, "data", "ep");
    legend->AddEntry(h_mc, "#gamma#gamma+jets", "f");
    legend->SetLineColor(0);
    legend->Draw("same");
    //--------------------------------------------------------------------------------
    TLatex latex;
    latex.SetNDC(kTRUE);
    latex.SetTextFont(43);
    latex.SetTextSize(22);
    latex.SetTextAlign(0);
    latex.DrawLatex(0.10, 0.91, "#bf{CMS} #it{Preliminary}");
    latex.DrawLatex(0.58, 0.91, "35.9 fb^{-1} (2016, 13 TeV)");
    //--------------------------------------------------------------------------------
    c1->cd(); //Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
    pad2->SetTopMargin(0.01); 
    pad2->SetBottomMargin(0.4);
    pad2->SetAttLinePS(kBlack,1,2);
    pad2->SetGridx(1); //Upper and lower pads are joined
    pad2->Draw();
    pad2->cd(); //pad2 becomes current pad
    gPad->SetTicks(1,1);
    //--------------------
    h_ratio->SetMaximum(2.0);
    h_ratio->SetMinimum(0.0);
    h_ratio->SetTitle("");
    h_ratio->SetStats(0); //No statistics on lower plot
    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerSize(1.2);
    h_ratio->SetLineColor(kBlack);
    h_ratio->Draw("p,E1");
    //--------------------
    h_ratio->GetXaxis()->SetTitle("Number of vertices");
    h_ratio->GetXaxis()->SetTitleSize(20);
    h_ratio->GetXaxis()->SetTitleFont(43);
    h_ratio->GetXaxis()->SetTitleOffset(4.);
    h_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h_ratio->GetXaxis()->SetLabelSize(15);
    //--------------------
    h_ratio->GetYaxis()->SetTitle("Data/MC");
    h_ratio->GetYaxis()->SetNdivisions(5);
    h_ratio->GetYaxis()->SetTitleSize(20);
    h_ratio->GetYaxis()->SetTitleFont(43);
    h_ratio->GetYaxis()->SetTitleOffset(1.2);
    h_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h_ratio->GetYaxis()->SetLabelSize(15);
    //--------------------
    c1->Update();//update the value of pad2->GetUxmax().
    TLine line;
    line.SetLineStyle(2);
    line.DrawLine(pad2->GetUxmin(),0.5,pad2->GetUxmax(),0.5);
    line.DrawLine(pad2->GetUxmin(),1.0,pad2->GetUxmax(),1.0);
    line.DrawLine(pad2->GetUxmin(),1.5,pad2->GetUxmax(),1.5);
    line.DrawLine(pad2->GetUxmin(),2.0,pad2->GetUxmax(),2.0);
    //--------------------
    c1->SaveAs(output);
}
void RegisterHistogram(TFile *&file, const char* fileName, TH1D* &hist, const char* histName, int color){
    //printf("Registering histogram of %s\n", fileName);
    file = TFile::Open(fileName);
    hist = (TH1D*)file->Get(histName);
    hist->SetLineColor(color);//problematic in second run
    hist->SetLineWidth(2);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1.2);
}
