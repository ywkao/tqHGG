#include <string>
#include <TCanvas.h>
#include <TFile.h>
#include <THStack.h>
#include <TH1D.h>
#include <TLine.h>
#include <TPad.h>
#include <string>
#include "/home/ykao/legacy/CMSSW_9_4_10/src/t2cH/include/stack.h"
using namespace std;

const double TunableSigBranchingFraction = 0.0001; //The branching fraction of signal MC = 0.01%

void stackHist(){
    MakeStackHist("hist_num_jets");
    MakeStackHist("hist_num_btagged_jets");
    MakeStackHist("hist_num_nonbtagged_jets");
    MakeStackHist("hist_wboson_mass_spectrum");
    MakeStackHist("hist_diphoton_mass_spectrum");
    MakeStackHist("hist_bjet_pt");
    MakeStackHist("hist_bjet_eta");
    MakeStackHist("hist_bjet_phi");
    //MakeStackHist("hist_cjet_pt");
    //MakeStackHist("hist_cjet_eta");
    //MakeStackHist("hist_cjet_phi");
    MakeStackHist("hist_jet1_pt");
    MakeStackHist("hist_jet1_eta");
    MakeStackHist("hist_jet1_phi");
    MakeStackHist("hist_jet2_pt");
    MakeStackHist("hist_jet2_eta");
    MakeStackHist("hist_jet2_phi");
    MakeStackHist("hist_inv_mass_tqh");
    MakeStackHist("hist_inv_mass_tbw");
}
void MakeStackHist(const char* histName){
    TCanvas *c1 = new TCanvas("c1", "c1", 700, 800);
    THStack *stackHist = new THStack("stackHist", "");//If setting titles here, the x(y)-title will NOT be able to set later.
    bool isNumEtaPhi = isThisNumEtaPhi(histName);
    bool isMassSpectrum = isThisMassSpectrum(histName);

    TH1D  *hist_tqh_sig_ttpair;
    TH1D  *hist_tqh_sig_singletop;
    TH1D  *hist_tqh_sig_tt_hadronic;
    TH1D  *hist_tqh_sig_tt_leptonic;
    TH1D  *hist_tqh_sig[NUM_sig];
    TH1D  *hist_tqh_resbkg[NUM_resbkg+1];
    TH1D  *hist_tqh_nonresbkg[NUM_nonresbkg+1];
    //--------------------
    TH1D  *hist_tqh_GluGluHToGG;
    TH1D  *hist_tqh_DiPhotonJetsBox;
    TH1D  *hist_tqh_GJet;
    TH1D  *hist_tqh_QCD;
    TH1D  *hist_tqh_TGJets;
    TH1D  *hist_tqh_TTGG;
    TH1D  *hist_tqh_TTGJets;
    //--------------------
    TH1D  *hist_tqh_data[NUM_data+1];
    TH1D  *hist_tqh_mc_wosig; //1.Calculate data/bkg ratio. 2.Statistical bkg uncertainty.
    TH1D  *hist_tqh_mc_uncertainty; //1.Statistical MC uncertainty.
    TH1D  *hist_tqh_ratio;
    //===== Register histograms =====//
    for(int i=0; i<NUM_sig; i++) RegisterHistogram(fileNames_sig[i].c_str(), hist_tqh_sig[i], histName, kRed, true, false);
    for(int i=0; i<NUM_resbkg; i++) RegisterHistogram(fileNames_resbkg[i].c_str(), hist_tqh_resbkg[i], histName, kBlue-i-2, false, false);
    for(int i=0; i<NUM_nonresbkg; i++) RegisterHistogram(fileNames_nonresbkg[i].c_str(), hist_tqh_nonresbkg[i], histName, kGreen+2, false, false);
    for(int i=0; i<NUM_data; i++) RegisterHistogram(fileNames_data[i].c_str(), hist_tqh_data[i], histName, kBlack, false, true);
    //===== Combine had/lep signal =====//
    hist_tqh_sig_tt_hadronic = (TH1D*) hist_tqh_sig[0]->Clone();
    for(int i=1; i<4; i++) hist_tqh_sig_tt_hadronic->Add(hist_tqh_sig[i]);
    hist_tqh_sig_tt_leptonic = (TH1D*) hist_tqh_sig[4]->Clone();
    for(int i=5; i<8; i++) hist_tqh_sig_tt_leptonic->Add(hist_tqh_sig[i]);
    hist_tqh_sig_singletop = (TH1D*) hist_tqh_sig[8]->Clone();
    for(int i=9; i<NUM_sig; i++) hist_tqh_sig_singletop->Add(hist_tqh_sig[i]);
    hist_tqh_sig_singletop->SetLineColor(kPink);
    hist_tqh_sig_singletop->SetLineStyle(2);
    hist_tqh_sig_ttpair = (TH1D*) hist_tqh_sig_tt_hadronic->Clone();
    hist_tqh_sig_ttpair->Add(hist_tqh_sig_tt_leptonic);
    //===== Combine backgournds =====//
    hist_tqh_resbkg[NUM_resbkg] = (TH1D*) hist_tqh_resbkg[0]->Clone();
    for(int i=1; i<NUM_resbkg; i++) hist_tqh_resbkg[NUM_resbkg]->Add(hist_tqh_resbkg[i]);
    hist_tqh_nonresbkg[NUM_nonresbkg] = (TH1D*) hist_tqh_nonresbkg[0]->Clone();
    for(int i=1; i<NUM_nonresbkg; i++) hist_tqh_nonresbkg[NUM_nonresbkg]->Add(hist_tqh_nonresbkg[i]);
    //===== Combine data =====//
    hist_tqh_data[NUM_data] = (TH1D*) hist_tqh_data[0]->Clone();
    for(int i=1; i<NUM_data; i++) hist_tqh_data[NUM_data]->Add(hist_tqh_data[i]);
    //===== Combine mc w/o sig =====//
    hist_tqh_mc_wosig = (TH1D*) hist_tqh_resbkg[NUM_resbkg]->Clone();
    hist_tqh_mc_wosig->Add(hist_tqh_nonresbkg[NUM_nonresbkg]);
    //===== Calculate relative MC uncertainty =====//
    hist_tqh_mc_uncertainty = (TH1D*)hist_tqh_mc_wosig->Clone();
    int nbins = hist_tqh_mc_wosig->GetNbinsX();
    //printf("nbins=%d\n", nbins);
    for(int i=0; i<nbins; i++){
        double mean = hist_tqh_mc_wosig->GetBinContent(i+1);
        double error = hist_tqh_mc_wosig->GetBinError(i+1);
        double upper_rel_error = (mean==0.) ? 0. : (mean+error)/mean;
        double lower_rel_error = (mean==0.) ? 0. : (mean-error)/mean;
        double error_mean = (upper_rel_error + lower_rel_error)/2.;
        double error_error = (upper_rel_error - lower_rel_error)/2.;
        hist_tqh_mc_uncertainty->SetBinContent(i+1, error_mean);
        hist_tqh_mc_uncertainty->SetBinError(i+1, error_error);
        //printf("bin = %d, mean=%f, error=%f, error_error=%f\n", i+1, mean, error, error_error);
    }
    //===== Make data-mc ratio histogram =====//
    hist_tqh_ratio = (TH1D*) hist_tqh_data[NUM_data]->Clone();
    hist_tqh_ratio->Divide(hist_tqh_mc_wosig);
    //===== Make stack histograms =====//
    stackHist->Add(hist_tqh_resbkg[NUM_resbkg]);
    stackHist->Add(hist_tqh_nonresbkg[NUM_nonresbkg]);
    //===== Draw upper plots =====//
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
    pad1->SetBottomMargin(0); //Upper and lower pads are joined
    pad1->SetAttLinePS(kBlack,1,2);
    pad1->Draw();
    pad1->cd(); //pad1 becomes current pad
    gPad->SetTicks(1,1);
    //--------------------
    gPad->SetLogy(1);
    stackHist->SetMaximum(isNumEtaPhi ? 1e+9 : 1e+6);
    stackHist->SetMinimum(5e-1);
    //if(isMassSpectrum){
    //    stackHist->SetMaximum(100);
    //    stackHist->SetMinimum(0);
    //} else{
    //    gPad->SetLogy(1);
    //    stackHist->SetMaximum(isNumEtaPhi ? 1e+9 : 1e+6);
    //    stackHist->SetMinimum(5e-1);
    //}
    stackHist->Draw("hist");
    hist_tqh_sig_ttpair->Draw("hist,same");
    hist_tqh_sig_singletop->Draw("hist,same");
    hist_tqh_mc_wosig->SetLineWidth(0);
    hist_tqh_mc_wosig->SetMarkerColor(kGray+1);
    hist_tqh_mc_wosig->SetFillColor(kGray+1);
    hist_tqh_mc_wosig->SetFillStyle(3002);
    hist_tqh_mc_wosig->Draw("E2,same");
    hist_tqh_data[NUM_data]->Draw("p,E1,same");
    //--------------------
    double BinWidth = hist_tqh_data[NUM_data]->GetXaxis()->GetBinWidth(1);//Take the width of the first bin as a representative.
    string yTitle = GetYtitleAccordingToHistName(histName, BinWidth);
    stackHist->GetYaxis()->SetTitle(yTitle.c_str());
    stackHist->GetYaxis()->SetTitleSize(20);
    stackHist->GetYaxis()->SetTitleFont(43);
    stackHist->GetYaxis()->SetTitleOffset(1.2);
    stackHist->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    stackHist->GetYaxis()->SetLabelSize(15);
    //--------------------
    TLegend *legend = new TLegend(0.30,0.55,0.85,0.85);
    legend->SetNColumns(2);
    //legend->AddEntry(hist_tqh_data[NUM_data], "Data (DoubleEG B-F)", "lep");
    legend->AddEntry(hist_tqh_data[NUM_data], "Observed", "lep");
    legend->AddEntry(hist_tqh_sig_ttpair, Form("t#bar{t} (tqH, BF=%.2f%%)", TunableSigBranchingFraction*100), "f");
    legend->AddEntry(hist_tqh_resbkg[NUM_resbkg], "Resonant bkg", "f");
    legend->AddEntry(hist_tqh_sig_singletop, Form("single top (tqH, BF=%.2f%%)", TunableSigBranchingFraction*100), "f");
    legend->AddEntry(hist_tqh_nonresbkg[NUM_nonresbkg], "Non-resonant bkg", "f");
    legend->AddEntry(hist_tqh_mc_wosig, "Bkg uncertainty (stat. only)", "f");
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
    pad2->SetAttLinePS(kBlack,1,2);
    pad2->SetGridx(1); //Upper and lower pads are joined
    pad2->Draw();
    pad2->cd(); //pad2 becomes current pad
    gPad->SetTicks(1,1);
    //--------------------
    hist_tqh_ratio->SetMaximum(2.25);
    hist_tqh_ratio->SetMinimum(0);
    hist_tqh_ratio->SetTitle("");
    hist_tqh_ratio->SetStats(0); //No statistics on lower plot
    hist_tqh_ratio->Draw("p,E1");
    hist_tqh_mc_uncertainty->SetFillColor(kGray);
    hist_tqh_mc_uncertainty->SetLineColor(kGray);
    hist_tqh_mc_uncertainty->Draw("E2,same");
    hist_tqh_ratio->Draw("p,E1,same");
    //--------------------
    hist_tqh_ratio->GetYaxis()->SetTitle("Obs/Exp");
    hist_tqh_ratio->GetYaxis()->SetNdivisions(5);
    hist_tqh_ratio->GetYaxis()->SetTitleSize(20);
    hist_tqh_ratio->GetYaxis()->SetTitleFont(43);
    hist_tqh_ratio->GetYaxis()->SetTitleOffset(1.2);
    hist_tqh_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hist_tqh_ratio->GetYaxis()->SetLabelSize(15);
    //--------------------
    //hist_tqh_ratio->GetXaxis()->SetTitle("Invariant mass of diphoton + jet [GeV/c^{2}]");
    string xTitle = GetXtitleAccordingToHistName(histName);
    hist_tqh_ratio->GetXaxis()->SetTitle(xTitle.c_str());
    hist_tqh_ratio->GetXaxis()->SetTitleSize(20);
    hist_tqh_ratio->GetXaxis()->SetTitleFont(43);
    hist_tqh_ratio->GetXaxis()->SetTitleOffset(4.);
    hist_tqh_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hist_tqh_ratio->GetXaxis()->SetLabelSize(15);
    //--------------------
    c1->Update();//update the value of pad2->GetUxmax().
    TLine line;
    line.SetLineStyle(2);
    line.DrawLine(pad2->GetUxmin(),0.5,pad2->GetUxmax(),0.5);
    line.DrawLine(pad2->GetUxmin(),1.0,pad2->GetUxmax(),1.0);
    line.DrawLine(pad2->GetUxmin(),1.5,pad2->GetUxmax(),1.5);
    line.DrawLine(pad2->GetUxmin(),2.0,pad2->GetUxmax(),2.0);
    c1->SaveAs(Form("plots/stack_%s.png", histName));
}

bool isThisMassSpectrum(const char* histName){
    if((string)histName == "hist_num_jets") return false;
    if((string)histName == "hist_num_btagged_jets") return false;
    if((string)histName == "hist_num_nonbtagged_jets") return false;
    if((string)histName == "hist_wboson_mass_spectrum") return true;
    if((string)histName == "hist_diphoton_mass_spectrum") return true;
    if((string)histName == "hist_bjet_pt") return false;
    if((string)histName == "hist_jet1_pt") return false;
    if((string)histName == "hist_jet2_pt") return false;
    if((string)histName == "hist_cjet_pt") return false;
    if((string)histName == "hist_bjet_eta") return false;
    if((string)histName == "hist_jet1_eta") return false;
    if((string)histName == "hist_jet2_eta") return false;
    if((string)histName == "hist_cjet_eta") return false;
    if((string)histName == "hist_bjet_phi") return false;
    if((string)histName == "hist_jet1_phi") return false;
    if((string)histName == "hist_jet2_phi") return false;
    if((string)histName == "hist_cjet_phi") return false;
    if((string)histName == "hist_inv_mass_tch") return true;
    if((string)histName == "hist_inv_mass_tbw") return true;
    return false;
}

bool isThisNumEtaPhi(const char* histName){
    if((string)histName == "hist_num_jets") return true;
    if((string)histName == "hist_num_btagged_jets") return true;
    if((string)histName == "hist_num_nonbtagged_jets") return true;
    if((string)histName == "hist_wboson_mass_spectrum") return false;
    if((string)histName == "hist_diphoton_mass_spectrum") return false;
    if((string)histName == "hist_bjet_pt") return false;
    if((string)histName == "hist_jet1_pt") return false;
    if((string)histName == "hist_jet2_pt") return false;
    if((string)histName == "hist_cjet_pt") return false;
    if((string)histName == "hist_bjet_eta") return true;
    if((string)histName == "hist_jet1_eta") return true;
    if((string)histName == "hist_jet2_eta") return true;
    if((string)histName == "hist_cjet_eta") return true;
    if((string)histName == "hist_bjet_phi") return true;
    if((string)histName == "hist_jet1_phi") return true;
    if((string)histName == "hist_jet2_phi") return true;
    if((string)histName == "hist_cjet_phi") return true;
    if((string)histName == "hist_inv_mass_tch") return false;
    if((string)histName == "hist_inv_mass_tbw") return false;
    return false;
}

string GetXtitleAccordingToHistName(const char* histName){
    if((string)histName == "hist_num_jets") return "Number of jets";
    if((string)histName == "hist_num_btagged_jets") return "Number of btagged jets";
    if((string)histName == "hist_num_nonbtagged_jets") return "Number of nonbtagged jets";
    if((string)histName == "hist_wboson_mass_spectrum") return "Invariant mass of di-jets [GeV/c^{2}]";
    if((string)histName == "hist_diphoton_mass_spectrum") return "Invariant mass of di-photon [GeV/c^{2}]";
    if((string)histName == "hist_bjet_pt") return "Pt of bjet [GeV/c]";
    if((string)histName == "hist_jet1_pt") return "Pt of jet1 [GeV/c]";
    if((string)histName == "hist_jet2_pt") return "Pt of jet2 [GeV/c]";
    if((string)histName == "hist_cjet_pt") return "Pt of cjet [GeV/c]";
    if((string)histName == "hist_bjet_eta") return "Eta of bjet";
    if((string)histName == "hist_jet1_eta") return "Eta of jet1";
    if((string)histName == "hist_jet2_eta") return "Eta of jet2";
    if((string)histName == "hist_cjet_eta") return "Eta of cjet";
    if((string)histName == "hist_bjet_phi") return "Phi of bjet";
    if((string)histName == "hist_jet1_phi") return "Phi of jet1";
    if((string)histName == "hist_jet2_phi") return "Phi of jet2";
    if((string)histName == "hist_cjet_phi") return "Phi of cjet";
    if((string)histName == "hist_inv_mass_tch") return "Invariant mass of diphoton + jet [GeV/c^{2}]";
    if((string)histName == "hist_inv_mass_tbw") return "Invariant mass of w boson + b-jet [GeV/c^{2}]";
    return "";
}

string GetYtitleAccordingToHistName(const char* histName, double BinWidth){
    string str_ytitle_1("Entries");
    string str_ytitle_2(Form("Entries / %.2f [GeV]", BinWidth));
    if((string)histName == "hist_num_jets") return str_ytitle_1;
    if((string)histName == "hist_num_btagged_jets") return str_ytitle_1;
    if((string)histName == "hist_num_nonbtagged_jets") return str_ytitle_1;
    if((string)histName == "hist_bjet_eta") return str_ytitle_1;
    if((string)histName == "hist_jet1_eta") return str_ytitle_1;
    if((string)histName == "hist_jet2_eta") return str_ytitle_1;
    if((string)histName == "hist_cjet_eta") return str_ytitle_1;
    if((string)histName == "hist_bjet_phi") return str_ytitle_1;
    if((string)histName == "hist_jet1_phi") return str_ytitle_1;
    if((string)histName == "hist_jet2_phi") return str_ytitle_1;
    if((string)histName == "hist_cjet_phi") return str_ytitle_1;
    if((string)histName == "hist_bjet_pt") return str_ytitle_2;
    if((string)histName == "hist_jet1_pt") return str_ytitle_2;
    if((string)histName == "hist_jet2_pt") return str_ytitle_2;
    if((string)histName == "hist_cjet_pt") return str_ytitle_2;
    if((string)histName == "hist_wboson_mass_spectrum") return str_ytitle_2;
    if((string)histName == "hist_diphoton_mass_spectrum") return str_ytitle_2;
    if((string)histName == "hist_inv_mass_tch") return str_ytitle_2;
    if((string)histName == "hist_inv_mass_tbw") return str_ytitle_2;
    return "";
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
