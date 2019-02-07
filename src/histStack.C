#include <string>
#include <TCanvas.h>
#include <TFile.h>
#include <THStack.h>
#include <TH1D.h>
#include <TPad.h>
#include "/home/ykao/legacy/CMSSW_9_4_10/src/t2cH/include/stack.h"

const double TunableSigBranchingFraction = 0.0003; //The branching fraction of signal MC = 0.03%

void histStack(){
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    THStack *stackHist = new THStack("stackHist", ";Invariant mass of diphoton + jet [GeV/c^{2}];Entries");

    TH1D  *hist_inv_mass_tch_sig_hadronic;
    TH1D  *hist_inv_mass_tch_sig_leptonic;
    TH1D  *hist_inv_mass_tch_sig[NUM_sig];
    TH1D  *hist_inv_mass_tch_resbkg[NUM_resbkg+1];
    TH1D  *hist_inv_mass_tch_nonresbkg[NUM_nonresbkg+1];
    TH1D  *hist_inv_mass_tch_data[NUM_data+1];
    TH1D  *hist_inv_mass_tch_mc_wosig;
    TH1D  *hist_inv_mass_tch_ratio;
    //===== Register histograms =====//
    for(int i=0; i<NUM_sig; i++) RegisterHistogram(fileNames_sig[i].c_str(), hist_inv_mass_tch_sig[i], "hist_inv_mass_tch", kRed-i, true, false);
    for(int i=0; i<NUM_resbkg; i++) RegisterHistogram(fileNames_resbkg[i].c_str(), hist_inv_mass_tch_resbkg[i], "hist_inv_mass_tch", kBlue-i-2, false, false);
    for(int i=0; i<NUM_nonresbkg; i++) RegisterHistogram(fileNames_nonresbkg[i].c_str(), hist_inv_mass_tch_nonresbkg[i], "hist_inv_mass_tch", kGreen+2, false, false);
    for(int i=0; i<NUM_data; i++) RegisterHistogram(fileNames_data[i].c_str(), hist_inv_mass_tch_data[i], "hist_inv_mass_tch", kBlack, false, true);
    //===== Combine had/lep signal =====//
    hist_inv_mass_tch_sig_hadronic = (TH1D*) hist_inv_mass_tch_sig[0]->Clone();
    for(int i=1; i<NUM_sig-4; i++) hist_inv_mass_tch_sig_hadronic->Add(hist_inv_mass_tch_sig[i]);
    hist_inv_mass_tch_sig_leptonic = (TH1D*) hist_inv_mass_tch_sig[4]->Clone();
    for(int i=5; i<NUM_sig; i++) hist_inv_mass_tch_sig_leptonic->Add(hist_inv_mass_tch_sig[i]);
    //===== Combine backgournds =====//
    hist_inv_mass_tch_resbkg[NUM_resbkg] = (TH1D*) hist_inv_mass_tch_resbkg[0]->Clone();
    for(int i=1; i<NUM_resbkg; i++) hist_inv_mass_tch_resbkg[NUM_resbkg]->Add(hist_inv_mass_tch_resbkg[i]);
    hist_inv_mass_tch_nonresbkg[NUM_nonresbkg] = (TH1D*) hist_inv_mass_tch_nonresbkg[0]->Clone();
    for(int i=1; i<NUM_nonresbkg; i++) hist_inv_mass_tch_nonresbkg[NUM_nonresbkg]->Add(hist_inv_mass_tch_nonresbkg[i]);
    //===== Combine data =====//
    hist_inv_mass_tch_data[NUM_data] = (TH1D*) hist_inv_mass_tch_data[0]->Clone();
    for(int i=1; i<NUM_data; i++) hist_inv_mass_tch_data[NUM_data]->Add(hist_inv_mass_tch_data[i]);
    //===== Combine mc w/o sig =====//
    hist_inv_mass_tch_mc_wosig = (TH1D*) hist_inv_mass_tch_resbkg[NUM_resbkg]->Clone();
    hist_inv_mass_tch_mc_wosig->Add(hist_inv_mass_tch_nonresbkg[NUM_nonresbkg]);
    //===== Make data-mc ratio histogram =====//
    hist_inv_mass_tch_ratio = (TH1D*) hist_inv_mass_tch_data[NUM_data]->Clone();
    hist_inv_mass_tch_ratio->Divide(hist_inv_mass_tch_mc_wosig);
    //===== Make stack histograms =====//
    stackHist->Add(hist_inv_mass_tch_resbkg[NUM_resbkg]);
    stackHist->Add(hist_inv_mass_tch_nonresbkg[NUM_nonresbkg]);
    //===== Draw upper plots =====//
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
    pad1->SetBottomMargin(0); //Upper and lower pads are joined
    pad1->Draw();
    pad1->cd(); //pad1 becomes current pad
    gPad->SetTicks(1,1);
    gPad->SetLogy(1);
    //--------------------
    //stackHist->SetMaximum(1e+3);
    stackHist->SetMaximum(1e+7);
    stackHist->SetMinimum(1e-4);
    stackHist->Draw("hist");
    hist_inv_mass_tch_sig_leptonic->Draw("h,same");
    hist_inv_mass_tch_sig_hadronic->Draw("h,same");
    hist_inv_mass_tch_data[NUM_data]->Draw("p,E1,same");
    //--------------------
    stackHist->GetYaxis()->SetTitleSize(20);
    stackHist->GetYaxis()->SetTitleFont(43);
    stackHist->GetYaxis()->SetTitleOffset(1.2);
    stackHist->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    stackHist->GetYaxis()->SetLabelSize(15);
    //--------------------
    TLegend *legend = new TLegend(0.60,0.62,0.87,0.87);
    legend->AddEntry(hist_inv_mass_tch_resbkg[NUM_resbkg], "Resonant bkg", "f");
    legend->AddEntry(hist_inv_mass_tch_nonresbkg[NUM_nonresbkg], "Non-resonant bkg", "f");
    legend->AddEntry(hist_inv_mass_tch_sig_hadronic, Form("Hadronic t#bar{t} (tqH, BF=%.2f%%)", TunableSigBranchingFraction*100), "f");
    legend->AddEntry(hist_inv_mass_tch_sig_leptonic, Form("Leptonic t#bar{t} (tqH, BF=%.2f%%)", TunableSigBranchingFraction*100), "f");
    legend->AddEntry(hist_inv_mass_tch_data[NUM_data], "DoubleEG B-F", "lep");
    legend->SetLineColor(0);
    legend->Draw("same");
    //--------------------
    TLatex latex, latex_lumi;
    latex.SetNDC(kTRUE);
    latex.SetTextFont(43);
    latex.SetTextSize(22);
    latex.SetTextAlign(13);
    latex.DrawLatex(0.14, 0.84, "#bf{CMS}");
    latex.DrawLatex(0.14, 0.78, "#it{Preliminary}");
    //--------------------
    latex_lumi.SetNDC(kTRUE);
    latex_lumi.SetTextFont(43);
    latex_lumi.SetTextSize(22);
    latex_lumi.SetTextAlign(31);
    latex_lumi.DrawLatex(0.89, 0.92, "42 fb^{-1} (2017, 13 TeV)");
    //===== Draw lower plots =====//
    c1->cd(); //Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
    pad2->SetTopMargin(0.01); 
    pad2->SetBottomMargin(0.4);
    pad2->SetGridx(1); //Upper and lower pads are joined
    pad2->Draw();
    pad2->cd(); //pad2 becomes current pad
    gPad->SetTicks(1,1);
    //--------------------
    //hist_inv_mass_tch_ratio->SetMaximum(2);
    //hist_inv_mass_tch_ratio->SetMinimum(0);
    hist_inv_mass_tch_ratio->SetTitle("");
    hist_inv_mass_tch_ratio->SetStats(0); //No statistics on lower plot
    hist_inv_mass_tch_ratio->Draw("p,E1");
    //--------------------
    hist_inv_mass_tch_ratio->GetYaxis()->SetTitle("Data/MC");
    hist_inv_mass_tch_ratio->GetYaxis()->SetNdivisions(505);
    hist_inv_mass_tch_ratio->GetYaxis()->SetTitleSize(20);
    hist_inv_mass_tch_ratio->GetYaxis()->SetTitleFont(43);
    hist_inv_mass_tch_ratio->GetYaxis()->SetTitleOffset(1.2);
    hist_inv_mass_tch_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hist_inv_mass_tch_ratio->GetYaxis()->SetLabelSize(15);
    //--------------------
    hist_inv_mass_tch_ratio->GetXaxis()->SetTitleSize(20);
    hist_inv_mass_tch_ratio->GetXaxis()->SetTitleFont(43);
    hist_inv_mass_tch_ratio->GetXaxis()->SetTitleOffset(4.);
    hist_inv_mass_tch_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hist_inv_mass_tch_ratio->GetXaxis()->SetLabelSize(15);
    c1->SaveAs("plots/stack_hist_inv_mass_tqh.png");
}

void RegisterHistogram(const char* fileName, TH1D* &hist, const char* histName, int color, bool isSigMC = true, bool isData = false){
    printf("Registering histogram of %s\n", fileName);
    TFile *file = TFile::Open(fileName);
    hist = (TH1D*)file->Get(histName);
    if(isSigMC) hist->Scale(TunableSigBranchingFraction);
    if(!isSigMC) hist->SetFillColor(color);
    hist->SetLineColor(color);
    hist->SetLineWidth(2);
    if(isData) hist->SetMarkerStyle(20);
    if(isData) hist->SetMarkerSize(1.2);
}
