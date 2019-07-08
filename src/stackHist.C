#include <string>
#include <TCanvas.h>
#include <TFile.h>
#include <THStack.h>
#include <TH1D.h>
#include <TLine.h>
#include <TPad.h>
#include <string>
#include "../include/stack.h"
using namespace std;

bool bool_isHadronic;
bool bool_isHut = true;
const double TunableSigBranchingFraction = 2.0; //The branching fraction of signal MC = 0.01%
//const double TunableSigBranchingFraction = 0.001; //The branching fraction of signal MC = 0.01%
bool PrintHistInfo = false;
bool PrintTexStyle;

void stackHist(const char* channel){
    if((string)channel == "hadronic") bool_isHadronic = true; else bool_isHadronic = false;
    if(bool_isHadronic) printf("[CHECK] from macro src/stackHist.C: isHadronic!\n");
    if(!bool_isHadronic) printf("[CHECK] from macro src/stackHist.C: isLeptonic!\n");
    //------------------------
    MakeStackHist("hist_EvtInfo_NPu");
    MakeStackHist("hist_EvtInfo_Rho");
    MakeStackHist("hist_EvtInfo_NVtx");
    //------------------------
    MakeStackHist("hist_DiPhoInfo_mass");
    MakeStackHist("hist_DiPhoInfo_pt");
    MakeStackHist("hist_DiPhoInfo_eta");
    MakeStackHist("hist_DiPhoInfo_phi");
    MakeStackHist("hist_DiPhoInfo_energy");
    MakeStackHist("hist_DiPhoInfo_leadPt");
    MakeStackHist("hist_DiPhoInfo_leadEta");
    MakeStackHist("hist_DiPhoInfo_leadPhi");
    MakeStackHist("hist_DiPhoInfo_leadE");
    MakeStackHist("hist_DiPhoInfo_leadhoe");
    MakeStackHist("hist_DiPhoInfo_leadIDMVA");
    MakeStackHist("hist_DiPhoInfo_subleadPt");
    MakeStackHist("hist_DiPhoInfo_subleadEta");
    MakeStackHist("hist_DiPhoInfo_subleadPhi");
    MakeStackHist("hist_DiPhoInfo_subleadE");
    MakeStackHist("hist_DiPhoInfo_subleadhoe");
    MakeStackHist("hist_DiPhoInfo_subleadIDMVA");
    //------------------------
    MakeStackHist("hist_ElecInfo_Size");
    MakeStackHist("hist_MuonInfo_Size");
    MakeStackHist("hist_num_leptons");// # of selected objects.
    MakeStackHist("hist_num_electrons");// # of selected objects.
    MakeStackHist("hist_num_muons");// # of selected objects.
    MakeStackHist("hist_ElecInfo_electron_pt");
    MakeStackHist("hist_ElecInfo_electron_eta");
    MakeStackHist("hist_ElecInfo_electron_phi");
    MakeStackHist("hist_ElecInfo_electron_energy");
    MakeStackHist("hist_ElecInfo_electron_diphoton_deltaR");
    MakeStackHist("hist_MuonInfo_muon_pt");
    MakeStackHist("hist_MuonInfo_muon_eta");
    MakeStackHist("hist_MuonInfo_muon_phi");
    MakeStackHist("hist_MuonInfo_muon_energy");
    MakeStackHist("hist_MuonInfo_muon_diphoton_deltaR");
    //------------------------
    MakeStackHist("hist_jets_size");
    MakeStackHist("hist_num_jets");
    MakeStackHist("hist_JetInfo_jet_pt");
    MakeStackHist("hist_JetInfo_jet_eta");
    MakeStackHist("hist_JetInfo_jet_phi");
    MakeStackHist("hist_JetInfo_jet_energy");
    MakeStackHist("hist_JetInfo_jet_diphoton_deltaR");
    //------------------------
    MakeStackHist("hist_lepton_pt");
    MakeStackHist("hist_lepton_eta");
    MakeStackHist("hist_lepton_phi");
    MakeStackHist("hist_lepton_energy");
    MakeStackHist("hist_lepton_diphoton_deltaR");
    //------------------------
    MakeStackHist("hist_jet1_diphoton_deltaR");
    MakeStackHist("hist_jet2_diphoton_deltaR");
    MakeStackHist("hist_jet1_lepton_deltaR");
    MakeStackHist("hist_jet2_lepton_deltaR");
    //------------------------
    MakeStackHist("hist_jet1_pt");
    MakeStackHist("hist_jet1_eta");
    MakeStackHist("hist_jet1_phi");
    MakeStackHist("hist_jet1_energy");
    MakeStackHist("hist_jet2_pt");
    MakeStackHist("hist_jet2_eta");
    MakeStackHist("hist_jet2_phi");
    MakeStackHist("hist_jet2_energy");
    //------------------------
    MakeStackHist("hist_leading_bjet_pt");
    MakeStackHist("hist_leading_bjet_eta");
    MakeStackHist("hist_leading_bjet_phi");
    MakeStackHist("hist_leading_bjet_energy");
    MakeStackHist("hist_chosen_bjet_pt");
    MakeStackHist("hist_chosen_bjet_eta");
    MakeStackHist("hist_chosen_bjet_phi");
    MakeStackHist("hist_chosen_bjet_energy");
    //------------------------
    MakeStackHist("hist_inv_mass_dijet");
    MakeStackHist("hist_inv_mass_diphoton");
    MakeStackHist("hist_inv_mass_tbw");
}
void MakeStackHist(const char* histName){
    if((string)histName == "hist_EvtInfo_NVtx" || (string)histName == "hist_EvtInfo_NVtx_wopu") PrintTexStyle = true; else PrintTexStyle = false;
    TCanvas *c1 = new TCanvas("c1", "c1", 700, 800);
    THStack *stackHist = new THStack("stackHist", "");//If setting titles here, the x(y)-title will NOT be able to set later.
    bool considerQCD = true;
    bool isIDMVA = isThisIDMVA(histName);
    bool isNum = isThisNum(histName);
    bool isEtaPhi = isThisEtaPhi(histName);
    bool isMassSpectrum = isThisMassSpectrum(histName);
    bool isDijetSpectrum = isThisDijetSpectrum(histName);
    bool isDiPhotonSpectrum = isThisDiPhotonSpectrum(histName);
    bool isTopSpectrum = isThisTopSpectrum(histName);
    //--------------------
    TH1D  *hist_tqh_sig_ttpair;
    TH1D  *hist_tqh_sig_singletop;
    TH1D  *hist_tqh_sig_st_hadronic;
    TH1D  *hist_tqh_sig_st_hadronic_hut;
    TH1D  *hist_tqh_sig_st_hadronic_hct;
    TH1D  *hist_tqh_sig_st_leptonic;
    TH1D  *hist_tqh_sig_st_leptonic_hut;
    TH1D  *hist_tqh_sig_st_leptonic_hct;
    //---
    TH1D  *hist_tqh_sig_tt_hadronic;
    TH1D  *hist_tqh_sig_tt_hadronic_hut;
    TH1D  *hist_tqh_sig_tt_hadronic_hct;
    TH1D  *hist_tqh_sig_tt_leptonic;
    TH1D  *hist_tqh_sig_tt_leptonic_hut;
    TH1D  *hist_tqh_sig_tt_leptonic_hct;
    TH1D  *hist_tqh_sig[NUM_sig];
    TH1D  *hist_tqh_resbkg[NUM_resbkg+1];
    TH1D  *hist_tqh_nonresbkg[NUM_nonresbkg+1];
    //--------------------
    TH1D  *hist_tqh_ggH;
    TH1D  *hist_tqh_VBF;
    TH1D  *hist_tqh_VH;
    TH1D  *hist_tqh_ttH;
    //--------------------
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
    //===============================//
    //===== Register histograms =====//
    //===============================//
    TFile *file[NUM];
    for(int i=0; i<NUM_sig; i++)
        RegisterHistogram(file[i], fileNames_sig[i].c_str(), hist_tqh_sig[i], histName, kRed, true, false);
    for(int i=0; i<NUM_resbkg; i++)
        RegisterHistogram(file[i+NUM_sig], fileNames_resbkg[i].c_str(), hist_tqh_resbkg[i], histName, kBlue-i-2, false, false);
    for(int i=0; i<NUM_nonresbkg; i++)
        RegisterHistogram(file[i+NUM_sig+NUM_resbkg], fileNames_nonresbkg[i].c_str(), hist_tqh_nonresbkg[i], histName, kGreen+2, false, false);
    for(int i=0; i<NUM_data; i++)
        RegisterHistogram(file[i+NUM_sig+NUM_resbkg+NUM_nonresbkg], fileNames_data[i].c_str(), hist_tqh_data[i], histName, kBlack, false, true);
    //==================================//
    //===== Combine had/lep signal =====//
    //==================================//
    hist_tqh_sig_tt_hadronic = (TH1D*) hist_tqh_sig[0]->Clone(); for(int i=1; i<4; i++) hist_tqh_sig_tt_hadronic->Add(hist_tqh_sig[i]);
    hist_tqh_sig_tt_hadronic_hut = (TH1D*) hist_tqh_sig[0]->Clone();
    hist_tqh_sig_tt_hadronic_hut -> Add(hist_tqh_sig[1]);
    hist_tqh_sig_tt_hadronic_hct = (TH1D*) hist_tqh_sig[2]->Clone();
    hist_tqh_sig_tt_hadronic_hct -> Add(hist_tqh_sig[3]);
    //---
    hist_tqh_sig_tt_leptonic = (TH1D*) hist_tqh_sig[4]->Clone(); for(int i=5; i<8; i++) hist_tqh_sig_tt_leptonic->Add(hist_tqh_sig[i]);
    hist_tqh_sig_tt_leptonic_hut = (TH1D*) hist_tqh_sig[4]->Clone();
    hist_tqh_sig_tt_leptonic_hut -> Add(hist_tqh_sig[5]);
    hist_tqh_sig_tt_leptonic_hct = (TH1D*) hist_tqh_sig[6]->Clone();
    hist_tqh_sig_tt_leptonic_hct -> Add(hist_tqh_sig[7]);
    //--------------------
    hist_tqh_sig_st_hadronic = (TH1D*) hist_tqh_sig[8]->Clone(); hist_tqh_sig_st_hadronic -> Add(hist_tqh_sig[9]);
    hist_tqh_sig_st_hadronic_hut = (TH1D*) hist_tqh_sig[8]->Clone();
    hist_tqh_sig_st_hadronic_hct = (TH1D*) hist_tqh_sig[9]->Clone();
    //---
    hist_tqh_sig_st_leptonic = (TH1D*) hist_tqh_sig[10]->Clone(); hist_tqh_sig_st_leptonic -> Add(hist_tqh_sig[11]);
    hist_tqh_sig_st_leptonic_hut = (TH1D*) hist_tqh_sig[10]->Clone();
    hist_tqh_sig_st_leptonic_hct = (TH1D*) hist_tqh_sig[11]->Clone();
    //--------------------
    hist_tqh_sig_singletop = (TH1D*) hist_tqh_sig_st_hadronic->Clone();
    hist_tqh_sig_singletop->Add(hist_tqh_sig_st_leptonic);
    hist_tqh_sig_ttpair = (TH1D*) hist_tqh_sig_tt_hadronic->Clone();
    hist_tqh_sig_ttpair->Add(hist_tqh_sig_tt_leptonic);
    //===============================//
    //===== Combine backgournds =====//
    //===============================//
    hist_tqh_resbkg[NUM_resbkg] = (TH1D*) hist_tqh_resbkg[0]->Clone();
    for(int i=1; i<NUM_resbkg; i++) hist_tqh_resbkg[NUM_resbkg]->Add(hist_tqh_resbkg[i]);
    hist_tqh_nonresbkg[NUM_nonresbkg] = (TH1D*) hist_tqh_nonresbkg[0]->Clone();
    for(int i=1; i<NUM_nonresbkg; i++) hist_tqh_nonresbkg[NUM_nonresbkg]->Add(hist_tqh_nonresbkg[i]);
    //--------------------
    hist_tqh_ggH = (TH1D*) hist_tqh_resbkg[0]->Clone();
    hist_tqh_VBF = (TH1D*) hist_tqh_resbkg[1]->Clone();
    hist_tqh_VH  = (TH1D*) hist_tqh_resbkg[2]->Clone();
    hist_tqh_ttH = (TH1D*) hist_tqh_resbkg[3]->Clone();
    hist_tqh_ggH -> SetFillColor(kViolet);
    hist_tqh_ggH -> SetLineColor(kViolet);
    hist_tqh_VBF -> SetFillColor(kBlue-7);
    hist_tqh_VBF -> SetLineColor(kBlue-7);
    hist_tqh_VH  -> SetFillColor(kBlue-9);
    hist_tqh_VH  -> SetLineColor(kBlue-9);
    hist_tqh_ttH -> SetFillColor(kBlue-10);
    hist_tqh_ttH -> SetLineColor(kBlue-10);
    //--------------------
    hist_tqh_DiPhotonJetsBox = (TH1D*) hist_tqh_nonresbkg[0]->Clone();
    hist_tqh_DiPhotonJetsBox->Add(hist_tqh_nonresbkg[1]);
    hist_tqh_GJet = (TH1D*) hist_tqh_nonresbkg[2]->Clone();
    hist_tqh_GJet->Add(hist_tqh_nonresbkg[3]);
    hist_tqh_GJet->Add(hist_tqh_nonresbkg[4]);
    hist_tqh_TGJets = (TH1D*) hist_tqh_nonresbkg[5]->Clone();
    hist_tqh_TTGG = (TH1D*) hist_tqh_nonresbkg[6]->Clone();
    hist_tqh_TTGJets = (TH1D*) hist_tqh_nonresbkg[7]->Clone();
    hist_tqh_QCD = (TH1D*) hist_tqh_nonresbkg[8]->Clone();
    hist_tqh_QCD->Add(hist_tqh_nonresbkg[9]);
    hist_tqh_QCD->Add(hist_tqh_nonresbkg[10]);
    //--------------------
    hist_tqh_DiPhotonJetsBox->SetFillColor(kOrange+1);
    hist_tqh_DiPhotonJetsBox->SetLineColor(kOrange+1);
    hist_tqh_GJet->SetFillColor(kOrange+2);
    hist_tqh_GJet->SetLineColor(kOrange+2);
    hist_tqh_TTGJets->SetFillColor(kOrange+3);
    hist_tqh_TTGJets->SetLineColor(kOrange+3);
    hist_tqh_TGJets->SetFillColor(kGreen+1);
    hist_tqh_TGJets->SetLineColor(kGreen+1);
    hist_tqh_TTGG->SetFillColor(kGreen+4);
    hist_tqh_TTGG->SetLineColor(kGreen+4);
    hist_tqh_QCD->SetFillColor(kOrange);
    hist_tqh_QCD->SetLineColor(kOrange);
    //========================//
    //===== Combine data =====//
    //========================//
    hist_tqh_data[NUM_data] = (TH1D*) hist_tqh_data[0]->Clone();
    for(int i=1; i<NUM_data; i++) hist_tqh_data[NUM_data]->Add(hist_tqh_data[i]);
    //==============================//
    //===== Combine mc w/o sig =====//
    //==============================//
    hist_tqh_mc_wosig = (TH1D*) hist_tqh_resbkg[NUM_resbkg]->Clone();
    if(considerQCD) hist_tqh_mc_wosig->Add(hist_tqh_nonresbkg[NUM_nonresbkg]);
    else{for(int i=1; i<NUM_nonresbkg-3; i++) hist_tqh_mc_wosig->Add(hist_tqh_nonresbkg[i]);}
    //=====================================================//
    //===== Calculate yields from combined histograms =====//
    //=====================================================//
    if(PrintHistInfo){
        printf("===========================================================================\n");
        printf("Process \t\t\t Yields \t\t   Error \n");
        printf("===========================================================================\n");
        CalculateHistYields("QCD\t\t", hist_tqh_QCD);
        CalculateHistYields("GJet\t\t", hist_tqh_GJet);
        CalculateHistYields("DiPhotonJetsBox", hist_tqh_DiPhotonJetsBox);
        CalculateHistYields("TGJets\t\t", hist_tqh_TGJets);
        CalculateHistYields("TTGJets\t\t", hist_tqh_TTGJets);
        CalculateHistYields("TTGG\t\t", hist_tqh_TTGG);
        CalculateHistYields("ggH\t\t", hist_tqh_ggH);
        CalculateHistYields("VBF\t\t", hist_tqh_VBF);
        CalculateHistYields("VH\t\t", hist_tqh_VH);
        CalculateHistYields("ttH\t\t", hist_tqh_ttH);
        printf("---------------------------------------------------------------------------\n");
        CalculateHistYields("MC bkg\t\t", hist_tqh_mc_wosig);
        CalculateHistYields("Data\t\t", hist_tqh_data[NUM_data]);
        printf("\n\n");
    }
    if(PrintTexStyle){
        printf("Processes & \\multicolumn{2}{c}{Yields}\\\\ \n");
        printf("\\hline\\hline\n");
        CalculateHistYields("QCD\t\t", hist_tqh_QCD);
        CalculateHistYields("GJet\t\t", hist_tqh_GJet);
        CalculateHistYields("DiPhotonJetsBox", hist_tqh_DiPhotonJetsBox);
        CalculateHistYields("TGJets\t\t", hist_tqh_TGJets);
        CalculateHistYields("TTGJets\t\t", hist_tqh_TTGJets);
        CalculateHistYields("TTGG\t\t", hist_tqh_TTGG);
        CalculateHistYields("ggH\t\t", hist_tqh_ggH);
        CalculateHistYields("VBF\t\t", hist_tqh_VBF);
        CalculateHistYields("VH\t\t", hist_tqh_VH);
        CalculateHistYields("ttH\t\t", hist_tqh_ttH);
        printf("\\hline\n");
        CalculateHistYields("MC background\t\t", hist_tqh_mc_wosig);
        CalculateHistYields("Data\t\t", hist_tqh_data[NUM_data]);
        printf("\n\n");
    }
    //=============================================//
    //===== Calculate relative MC uncertainty =====//
    //=============================================//
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
    //========================================//
    //===== Make data-mc ratio histogram =====//
    //========================================//
    hist_tqh_ratio = (TH1D*) hist_tqh_data[NUM_data]->Clone();
    hist_tqh_ratio->Divide(hist_tqh_mc_wosig);
    //=================================//
    //===== Make stack histograms =====//
    //=================================//
    //stackHist->Add(hist_tqh_resbkg[NUM_resbkg]);
    stackHist->Add(hist_tqh_ggH);
    stackHist->Add(hist_tqh_VBF);
    stackHist->Add(hist_tqh_VH); 
    stackHist->Add(hist_tqh_ttH);
    //stackHist->Add(hist_tqh_nonresbkg[NUM_nonresbkg]);
    stackHist->Add(hist_tqh_DiPhotonJetsBox);
    stackHist->Add(hist_tqh_GJet);
    stackHist->Add(hist_tqh_TTGJets);
    stackHist->Add(hist_tqh_TGJets);
    stackHist->Add(hist_tqh_TTGG);
    if(considerQCD) stackHist->Add(hist_tqh_QCD);
    //================================================//
    //===== Choose proper signal hist to present =====//
    //================================================//
    TH1D *hist_sig_tt, *hist_sig_st;
    if(bool_isHadronic && bool_isHut){
        hist_sig_tt = (TH1D*) hist_tqh_sig_tt_hadronic_hut->Clone(); 
        hist_sig_st = (TH1D*) hist_tqh_sig_st_hadronic_hut->Clone(); 
    } else if(bool_isHadronic && !bool_isHut){
        hist_sig_tt = (TH1D*) hist_tqh_sig_tt_hadronic_hct->Clone(); 
        hist_sig_st = (TH1D*) hist_tqh_sig_st_hadronic_hct->Clone(); 
    } else if(!bool_isHadronic && bool_isHut){
        hist_sig_tt = (TH1D*) hist_tqh_sig_tt_leptonic_hut->Clone(); 
        hist_sig_st = (TH1D*) hist_tqh_sig_st_leptonic_hut->Clone(); 
    } else{
        hist_sig_tt = (TH1D*) hist_tqh_sig_tt_leptonic_hct->Clone(); 
        hist_sig_st = (TH1D*) hist_tqh_sig_st_leptonic_hct->Clone(); 
    }
    //--------------------
    hist_sig_tt->SetLineColor(kRed);
    hist_sig_st->SetLineColor(kPink);
    hist_sig_st->SetLineStyle(2);
    //============================//
    //===== Draw upper plots =====//
    //============================//
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
    pad1->SetBottomMargin(0); //Upper and lower pads are joined
    pad1->SetAttLinePS(kBlack,1,2);
    pad1->Draw();
    pad1->cd(); //pad1 becomes current pad
    gPad->SetTicks(1,1);
    //--------------------
    hist_tqh_mc_wosig->SetLineWidth(0);
    hist_tqh_mc_wosig->SetMarkerColor(kGray+2);
    hist_tqh_mc_wosig->SetFillColor(kGray+2);
    hist_tqh_mc_wosig->SetFillStyle(3001);
    //--------------------
    stackHist->SetMinimum(0);
    double scale = 2.5;
    double max_stack = stackHist->GetMaximum();
    double max_data = hist_tqh_data[NUM_data]->GetMaximum();
    if(max_stack > max_data){
        stackHist->Draw("hist");
        hist_sig_tt->Draw("hist,same");
        hist_sig_st->Draw("hist,same");
        hist_tqh_mc_wosig->Draw("E2,same");
        hist_tqh_data[NUM_data]->Draw("p,E1,same");
        stackHist->SetMaximum(max_stack*scale);
    } else{
        hist_tqh_data[NUM_data]->SetStats(0);
        hist_tqh_data[NUM_data]->Draw("p,E1");
        stackHist->Draw("hist,same");
        hist_sig_tt->Draw("hist,same");
        hist_sig_st->Draw("hist,same");
        hist_tqh_mc_wosig->Draw("E2,same");
        hist_tqh_data[NUM_data]->Draw("p,E1,same");
        hist_tqh_data[NUM_data]->SetMaximum(max_data*scale);
    }
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
    legend->AddEntry(hist_tqh_data[NUM_data], "Observed", "lep");
    if(bool_isHadronic && bool_isHut)       legend->AddEntry(hist_sig_tt, Form("TT Hut (Had., BF=%.2f%%)", TunableSigBranchingFraction*100), "f");
    else if(!bool_isHadronic && bool_isHut) legend->AddEntry(hist_sig_tt, Form("TT Hut (Lep., BF=%.2f%%)", TunableSigBranchingFraction*100), "f");
    else if(bool_isHadronic && !bool_isHut) legend->AddEntry(hist_sig_tt, Form("TT Hct (Had., BF=%.2f%%)", TunableSigBranchingFraction*100), "f");
    else                                  legend->AddEntry(hist_sig_tt, Form("TT Hct (Lep., BF=%.2f%%)", TunableSigBranchingFraction*100), "f");
    legend->AddEntry(hist_tqh_mc_wosig, "Bkg uncertainty (stat. only)", "f");
    if(bool_isHadronic && bool_isHut)       legend->AddEntry(hist_sig_st, Form("ST Hut (Had., BF=%.2f%%)", TunableSigBranchingFraction*100), "f");
    else if(!bool_isHadronic && bool_isHut) legend->AddEntry(hist_sig_st, Form("ST Hut (Lep., BF=%.2f%%)", TunableSigBranchingFraction*100), "f");
    else if(bool_isHadronic && !bool_isHut) legend->AddEntry(hist_sig_st, Form("ST Hct (Had., BF=%.2f%%)", TunableSigBranchingFraction*100), "f");
    else                                  legend->AddEntry(hist_sig_st, Form("ST Hct (Lep., BF=%.2f%%)", TunableSigBranchingFraction*100), "f");
    legend->AddEntry(hist_tqh_ggH, "ggH", "f");
    legend->AddEntry(hist_tqh_VBF, "VBF", "f");
    legend->AddEntry(hist_tqh_VH , "VH" , "f"); 
    legend->AddEntry(hist_tqh_ttH, "ttH", "f");
    legend->AddEntry(hist_tqh_DiPhotonJetsBox, "DiPhotonJetsBox", "f");
    legend->AddEntry(hist_tqh_TGJets, "TGJets", "f");
    legend->AddEntry(hist_tqh_GJet, "GJet", "f");
    legend->AddEntry(hist_tqh_TTGG, "TTGG", "f");
    legend->AddEntry(hist_tqh_TTGJets, "TTGJets", "f");
    if(considerQCD) legend->AddEntry(hist_tqh_QCD, "QCD", "f");
    //legend->AddEntry(hist_tqh_sig_ttpair, Form("t#bar{t} (tqH, BF=%.2f%%)", TunableSigBranchingFraction*100), "f");
    //legend->AddEntry(hist_tqh_resbkg[NUM_resbkg], "Resonant bkg", "f");
    //legend->AddEntry(hist_tqh_sig_singletop, Form("single top (tqH, BF=%.2f%%)", TunableSigBranchingFraction*100), "f");
    //legend->AddEntry(hist_tqh_nonresbkg[NUM_nonresbkg], "Non-resonant bkg", "f");
    //legend->AddEntry(hist_tqh_mc_wosig, "Bkg uncertainty (stat. only)", "f");
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
    latex_lumi.DrawLatex(0.89, 0.92, "41.5 fb^{-1} (2017, 13 TeV)");
    //============================//
    //===== Draw lower plots =====//
    //============================//
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
    hist_tqh_ratio->SetMaximum(1.75);
    hist_tqh_ratio->SetMinimum(0.5);
    //hist_tqh_ratio->SetMaximum(2.25);
    //hist_tqh_ratio->SetMinimum(0);
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
    //============================//
    //====== Draw log scale ======//
    //============================//
    c1->cd(); pad1->cd();
    gPad->SetLogy(1);
    if(isMassSpectrum){
        //stackHist->SetMaximum(600);
        //c1->SaveAs(Form("plots/stack_%s_zoomin.png", histName));
        if(isDiPhotonSpectrum) stackHist->SetMaximum(1e+11);
        if(isDijetSpectrum)    stackHist->SetMaximum(1e+8);
        if(isTopSpectrum)      stackHist->SetMaximum(1e+8);
        stackHist->SetMinimum(5e-1);
        c1->SaveAs(Form("plots/log_scale/stack_%s_log.png", histName));
    }
    else stackHist->SetMaximum(isNum ? 1e+9 : 1e+9);

    for(int i=0; i<NUM; i++) file[i]->Close();
}

bool isThisIDMVA(const char* histName){
    if((string)histName == "hist_DiPhoInfo_leadIDMVA") return true;
    if((string)histName == "hist_DiPhoInfo_subleadIDMVA") return true;
    if((string)histName == "hist_DiPhoInfo_leadIDMVA_ori") return true;
    if((string)histName == "hist_DiPhoInfo_subleadIDMVA_ori") return true;
    return false;
}
bool isThisDijetSpectrum(const char* histName){
    if((string)histName == "hist_inv_mass_dijet") return true;
    return false;
}
bool isThisDiPhotonSpectrum(const char* histName){
    if((string)histName == "hist_inv_mass_diphoton") return true;
    if((string)histName == "hist_inv_mass_diphoton_ori") return true;
    return false;
}
bool isThisTopSpectrum(const char* histName){
    if((string)histName == "hist_inv_mass_tbw") return true;
    return false;
}
bool isThisMassSpectrum(const char* histName){
    if((string)histName == "hist_inv_mass_dijet") return true;
    if((string)histName == "hist_inv_mass_diphoton") return true;
    if((string)histName == "hist_inv_mass_diphoton_ori") return true;
    if((string)histName == "hist_inv_mass_tbw") return true;
    return false;
}
bool isThisNum(const char* histName){
    if((string)histName == "hist_ElecInfo_Size") return true;
    if((string)histName == "hist_MuonInfo_Size") return true;
    if((string)histName == "hist_jets_size") return true;
    if((string)histName == "hist_num_leptons") return true;// # of selected objects.
    if((string)histName == "hist_num_electrons") return true;// # of selected objects.
    if((string)histName == "hist_num_muons") return true;// # of selected objects.
    if((string)histName == "hist_num_jets") return true;
    return false;
}
bool isThisEtaPhi(const char* histName){
    if((string)histName == "hist_DiPhoInfo_leadPhi") return true;
    if((string)histName == "hist_DiPhoInfo_subleadPhi") return true;
    if((string)histName == "hist_DiPhoInfo_phi") return true;
    if((string)histName == "hist_ElecInfo_electron_phi") return true;
    if((string)histName == "hist_MuonInfo_muon_phi") return true;
    if((string)histName == "hist_JetInfo_jet_phi") return true;
    if((string)histName == "hist_lepton_phi") return true;
    if((string)histName == "hist_jet1_phi") return true;
    if((string)histName == "hist_jet2_phi") return true;
    if((string)histName == "hist_leading_bjet_phi") return true;
    if((string)histName == "hist_chosen_bjet_phi") return true;
    //--------------------
    if((string)histName == "hist_DiPhoInfo_leadEta") return true;
    if((string)histName == "hist_DiPhoInfo_subleadEta") return true;
    if((string)histName == "hist_DiPhoInfo_eta") return true;
    if((string)histName == "hist_ElecInfo_electron_eta") return true;
    if((string)histName == "hist_MuonInfo_muon_eta") return true;
    if((string)histName == "hist_JetInfo_jet_eta") return true;
    if((string)histName == "hist_lepton_eta") return true;
    if((string)histName == "hist_jet1_eta") return true;
    if((string)histName == "hist_jet2_eta") return true;
    if((string)histName == "hist_leading_bjet_eta") return true;
    if((string)histName == "hist_chosen_bjet_eta") return true;
    //--------------------
    return false;
}
string GetXtitleAccordingToHistName(const char* histName){
    if((string)histName == "hist_EvtInfo_NPu") return "Number of pile up";
    if((string)histName == "hist_EvtInfo_Rho") return "Rho";
    if((string)histName == "hist_EvtInfo_NVtx") return "Number of vertices";
    if((string)histName == "hist_EvtInfo_NVtx_wopu") return "Number of vertices (w/o pu)";
    //----------------------------------------------------------------------
    if((string)histName == "hist_ElecInfo_Size") return "Number of electrons (before any selection)";
    if((string)histName == "hist_MuonInfo_Size") return "Number of muons (before any selection)";
    if((string)histName == "hist_num_leptons") return "Number of leptons";// # of selected objects.
    if((string)histName == "hist_num_electrons") return "Number of electrons";// # of selected objects.
    if((string)histName == "hist_num_muons") return "Number of muons";// # of selected objects.
    //--------------------
    if((string)histName == "hist_ElecInfo_electron_pt") return "Pt of electron [GeV/c]";
    if((string)histName == "hist_ElecInfo_electron_eta") return "#eta of electron";
    if((string)histName == "hist_ElecInfo_electron_phi") return "#phi of electron";
    if((string)histName == "hist_ElecInfo_electron_energy") return "Energy of electron [GeV]";
    if((string)histName == "hist_ElecInfo_electron_diphoton_deltaR") return "deltaR between electron and diphoton";
    //--------------------
    if((string)histName == "hist_MuonInfo_muon_pt") return "Pt of muon [GeV/c]";
    if((string)histName == "hist_MuonInfo_muon_eta") return "#eta of muon";
    if((string)histName == "hist_MuonInfo_muon_phi") return "#phi of muon";
    if((string)histName == "hist_MuonInfo_muon_energy") return "Energy of muon [GeV]";
    if((string)histName == "hist_MuonInfo_muon_diphoton_deltaR") return "deltaR between muon and diphoton";
    //--------------------
    if((string)histName == "hist_lepton_pt") return "Pt of lepton [GeV/c]";
    if((string)histName == "hist_lepton_eta") return "#eta of lepton";
    if((string)histName == "hist_lepton_phi") return "#phi of lepton";
    if((string)histName == "hist_lepton_energy") return "Energy of lepton [GeV]";
    if((string)histName == "hist_lepton_diphoton_deltaR") return "deltaR between lepton and diphoton";
    //----------------------------------------------------------------------
    if((string)histName == "hist_leading_bjet_pt") return "Pt of bjet [GeV/c]";
    if((string)histName == "hist_leading_bjet_eta") return "#eta of bjet";
    if((string)histName == "hist_leading_bjet_phi") return "#phi of bjet";
    if((string)histName == "hist_leading_bjet_energy") return "Energy of bjet [GeV]";
    if((string)histName == "hist_chosen_bjet_pt") return "Pt of bjet [GeV/c]";
    if((string)histName == "hist_chosen_bjet_eta") return "#eta of bjet";
    if((string)histName == "hist_chosen_bjet_phi") return "#phi of bjet";
    if((string)histName == "hist_chosen_bjet_energy") return "Energy of bjet [GeV]";
    //----------------------------------------------------------------------
    if((string)histName == "hist_DiPhoInfo_mass") return "Diphoton invariant mass [GeV/c^2]";
    if((string)histName == "hist_DiPhoInfo_pt") return "Pt of diphoton [GeV/c]";
    if((string)histName == "hist_DiPhoInfo_eta") return "#eta of diphoton";
    if((string)histName == "hist_DiPhoInfo_phi") return "#phi of diphoton";
    if((string)histName == "hist_DiPhoInfo_energy") return "Energy of diphoton [GeV]";
    //--------------------
    if((string)histName == "hist_DiPhoInfo_leadPt") return "Pt of leading photon [GeV/c]";
    if((string)histName == "hist_DiPhoInfo_leadEta") return "#eta of leading photon";
    if((string)histName == "hist_DiPhoInfo_leadPhi") return "#phi of leading photon";
    if((string)histName == "hist_DiPhoInfo_leadE") return "Energy of leading photon [GeV]";
    if((string)histName == "hist_DiPhoInfo_leadhoe") return "Ratio of Hadronic energy over EM energy";
    if((string)histName == "hist_DiPhoInfo_leadIDMVA") return "IDMVA of leading photon";
    //--------------------
    if((string)histName == "hist_DiPhoInfo_subleadPt") return "Pt of SubleaddiPhoInfo [GeV/c]";
    if((string)histName == "hist_DiPhoInfo_subleadEta") return "#eta of subleading photon";
    if((string)histName == "hist_DiPhoInfo_subleadPhi") return "#phi of subleading photon";
    if((string)histName == "hist_DiPhoInfo_subleadE") return "Energy of subleading photon [GeV]";
    if((string)histName == "hist_DiPhoInfo_subleadhoe") return "Ratio of Hadronic energy over EM energy";
    if((string)histName == "hist_DiPhoInfo_subleadIDMVA") return "IDMVA of subleading photon";
    //----------------------------------------------------------------------
    if((string)histName == "hist_jets_size") return "Number of jets (before any selection)";
    if((string)histName == "hist_num_jets") return "Number of jets";
    if((string)histName == "hist_JetInfo_jet_pt") return "Pt of jet [GeV/c]";
    if((string)histName == "hist_JetInfo_jet_eta") return "#eta of jet";
    if((string)histName == "hist_JetInfo_jet_phi") return "#phi of jet";
    if((string)histName == "hist_JetInfo_jet_energy") return "Energy of jet [GeV]";
    if((string)histName == "hist_JetInfo_jet_diphoton_deltaR") return "deltaR between jet and diphoton";
    //------------------------
    if((string)histName == "hist_jet1_pt") return "Pt of jet1 [GeV/c]";
    if((string)histName == "hist_jet1_eta") return "#eta of jet1";
    if((string)histName == "hist_jet1_phi") return "#phi of jet1";
    if((string)histName == "hist_jet1_energy") return "Energy of jet1 [GeV]";
    if((string)histName == "hist_jet1_diphoton_deltaR") return "deltaR between jet1 and diphoton";
    if((string)histName == "hist_jet1_lepton_deltaR") return "deltaR between jet1 and lepton";
    //------------------------
    if((string)histName == "hist_jet2_pt") return "Pt of jet2 [GeV/c]";
    if((string)histName == "hist_jet2_eta") return "#eta of jet2";
    if((string)histName == "hist_jet2_phi") return "#phi of jet2";
    if((string)histName == "hist_jet2_energy") return "Energy of jet2 [GeV]";
    if((string)histName == "hist_jet2_diphoton_deltaR") return "deltaR between jet2 and diphoton";
    if((string)histName == "hist_jet2_lepton_deltaR") return "deltaR between jet2 and lepton";
    //----------------------------------------------------------------------
    if((string)histName == "hist_inv_mass_dijet") return "Invariant mass of dijet [GeV/c^{2}]";
    if((string)histName == "hist_inv_mass_diphoton") return "Invariant mass of diphoton [GeV/c^{2}]";
    if((string)histName == "hist_inv_mass_tbw") return "Invariant mass of dijet+bjet [GeV/c^{2}]";
    return "";
}

string GetYtitleAccordingToHistName(const char* histName, double BinWidth){
    string str_ytitle_1("Entries");
    string str_ytitle_2(Form("Entries / %.2f [GeV]", BinWidth));
    if((string)histName == "hist_EvtInfo_NPu") return str_ytitle_1;
    if((string)histName == "hist_EvtInfo_Rho") return str_ytitle_1;
    if((string)histName == "hist_EvtInfo_NVtx") return str_ytitle_1;
    if((string)histName == "hist_EvtInfo_NVtx_wopu") return str_ytitle_1;
    //----------------------------------------------------------------------
    if((string)histName == "hist_ElecInfo_Size") return str_ytitle_1;
    if((string)histName == "hist_MuonInfo_Size") return str_ytitle_1;
    if((string)histName == "hist_num_leptons") return str_ytitle_1;// # of selected objects.
    if((string)histName == "hist_num_electrons") return str_ytitle_1;// # of selected objects.
    if((string)histName == "hist_num_muons") return str_ytitle_1;// # of selected objects.
    //--------------------
    if((string)histName == "hist_ElecInfo_electron_pt") return str_ytitle_2;
    if((string)histName == "hist_ElecInfo_electron_eta") return str_ytitle_1;
    if((string)histName == "hist_ElecInfo_electron_phi") return str_ytitle_1;
    if((string)histName == "hist_ElecInfo_electron_energy") return str_ytitle_2;
    if((string)histName == "hist_ElecInfo_electron_diphoton_deltaR") return str_ytitle_1;
    //--------------------
    if((string)histName == "hist_MuonInfo_muon_pt") return str_ytitle_2;
    if((string)histName == "hist_MuonInfo_muon_eta") return str_ytitle_1;
    if((string)histName == "hist_MuonInfo_muon_phi") return str_ytitle_1;
    if((string)histName == "hist_MuonInfo_muon_energy") return str_ytitle_2;
    if((string)histName == "hist_MuonInfo_muon_diphoton_deltaR") return str_ytitle_1;
    //--------------------
    if((string)histName == "hist_lepton_pt") return str_ytitle_2;
    if((string)histName == "hist_lepton_eta") return str_ytitle_1;
    if((string)histName == "hist_lepton_phi") return str_ytitle_1;
    if((string)histName == "hist_lepton_energy") return str_ytitle_2;
    if((string)histName == "hist_lepton_diphoton_deltaR") return str_ytitle_1;
    //----------------------------------------------------------------------
    if((string)histName == "hist_leading_bjet_pt") return str_ytitle_2;
    if((string)histName == "hist_leading_bjet_eta") return str_ytitle_1;
    if((string)histName == "hist_leading_bjet_phi") return str_ytitle_1;
    if((string)histName == "hist_leading_bjet_energy") return str_ytitle_1;
    if((string)histName == "hist_chosen_bjet_pt") return str_ytitle_2;
    if((string)histName == "hist_chosen_bjet_eta") return str_ytitle_1;
    if((string)histName == "hist_chosen_bjet_phi") return str_ytitle_1;
    if((string)histName == "hist_chosen_bjet_energy") return str_ytitle_1;
    //----------------------------------------------------------------------
    if((string)histName == "hist_DiPhoInfo_mass") return str_ytitle_2;
    if((string)histName == "hist_DiPhoInfo_pt") return str_ytitle_2;
    if((string)histName == "hist_DiPhoInfo_eta") return str_ytitle_1;
    if((string)histName == "hist_DiPhoInfo_phi") return str_ytitle_1;
    if((string)histName == "hist_DiPhoInfo_energy") return str_ytitle_1;
    //--------------------
    if((string)histName == "hist_DiPhoInfo_leadPt") return str_ytitle_2;
    if((string)histName == "hist_DiPhoInfo_leadEta") return str_ytitle_1;
    if((string)histName == "hist_DiPhoInfo_leadPhi") return str_ytitle_1;
    if((string)histName == "hist_DiPhoInfo_leadE") return str_ytitle_2;
    if((string)histName == "hist_DiPhoInfo_leadhoe") return str_ytitle_1;
    if((string)histName == "hist_DiPhoInfo_leadIDMVA") return str_ytitle_1;
    //--------------------
    if((string)histName == "hist_DiPhoInfo_subleadPt") return str_ytitle_2;
    if((string)histName == "hist_DiPhoInfo_subleadEta") return str_ytitle_1;
    if((string)histName == "hist_DiPhoInfo_subleadPhi") return str_ytitle_1;
    if((string)histName == "hist_DiPhoInfo_subleadE") return str_ytitle_2;
    if((string)histName == "hist_DiPhoInfo_subleadhoe") return str_ytitle_1;
    if((string)histName == "hist_DiPhoInfo_subleadIDMVA") return str_ytitle_1;
    //----------------------------------------------------------------------
    if((string)histName == "hist_jets_size") return str_ytitle_1;
    if((string)histName == "hist_num_jets") return str_ytitle_1;
    if((string)histName == "hist_JetInfo_jet_pt") return str_ytitle_2;
    if((string)histName == "hist_JetInfo_jet_eta") return str_ytitle_1;
    if((string)histName == "hist_JetInfo_jet_phi") return str_ytitle_1;
    if((string)histName == "hist_JetInfo_jet_energy") return str_ytitle_2;
    if((string)histName == "hist_JetInfo_jet_diphoton_deltaR") return str_ytitle_1;
    //------------------------
    if((string)histName == "hist_jet1_pt") return str_ytitle_2;
    if((string)histName == "hist_jet1_eta") return str_ytitle_1;
    if((string)histName == "hist_jet1_phi") return str_ytitle_1;
    if((string)histName == "hist_jet1_energy") return str_ytitle_2;
    if((string)histName == "hist_jet1_diphoton_deltaR") return str_ytitle_1;
    if((string)histName == "hist_jet1_lepton_deltaR") return str_ytitle_1;
    //------------------------
    if((string)histName == "hist_jet2_pt") return str_ytitle_2;
    if((string)histName == "hist_jet2_eta") return str_ytitle_1;
    if((string)histName == "hist_jet2_phi") return str_ytitle_1;
    if((string)histName == "hist_jet2_energy") return str_ytitle_2;
    if((string)histName == "hist_jet2_diphoton_deltaR") return str_ytitle_1;
    if((string)histName == "hist_jet2_lepton_deltaR") return str_ytitle_2;
    //----------------------------------------------------------------------
    if((string)histName == "hist_inv_mass_dijet") return str_ytitle_2;
    if((string)histName == "hist_inv_mass_diphoton") return str_ytitle_2;
    if((string)histName == "hist_inv_mass_tbw") return str_ytitle_2;
    return "";
}

void RegisterHistogram(TFile *&file, const char* fileName, TH1D* &hist, const char* histName, int color, bool isSigMC = true, bool isData = false){
    //printf("Registering histogram of %s\n", fileName);
    //TFile *file = TFile::Open(fileName);// Will lead to problem of opening too many file.
    file = TFile::Open(fileName);
    hist = (TH1D*)file->Get(histName);
    if(isSigMC) hist->Scale(TunableSigBranchingFraction);
    if(!isSigMC) hist->SetFillColor(color);
    hist->SetLineColor(color);
    hist->SetLineWidth(2);
    if(isData) hist->SetMarkerStyle(20);
    if(isData) hist->SetMarkerSize(1.2);
    hist->SetTitle("");
    hist->SetStats(0);
}

void CalculateHistYields(const char *process, TH1D* hist){
    double totalYields = hist->Integral();
    double totalError = SumErrors(hist);
    if(PrintTexStyle) printf("%s & \t %15.2f & $\\pm$ \t %10.2f\\\\\n", process, totalYields, totalError);
    else printf("%s \t %15.2f \t %15.2f\n", process, totalYields, totalError);
}

double SumErrors(TH1D* hist){
    double totalError = 0.;
    for(int i=0; i<hist->GetNbinsX(); ++i){
        double error = hist->GetBinError(i+1);
        totalError += pow(error, 2);
    }
    totalError = sqrt(totalError);
    return totalError;
}
/*
    //------------------------
    hist_EvtInfo_NPu,
    hist_EvtInfo_Rho,
    hist_EvtInfo_NVtx,
    hist_EvtInfo_genweight,
    //------------------------
    hist_DiPhoInfo_mass,
    hist_DiPhoInfo_pt,
    hist_DiPhoInfo_energy,
    hist_DiPhoInfo_leadPt,
    hist_DiPhoInfo_leadE,
    hist_DiPhoInfo_leadhoe,
    hist_DiPhoInfo_leadIDMVA,
    hist_DiPhoInfo_subleadPt,
    hist_DiPhoInfo_subleadE,
    hist_DiPhoInfo_subleadhoe,
    hist_DiPhoInfo_subleadIDMVA,
    //------------------------
    hist_ElecInfo_electron_pt,
    hist_ElecInfo_electron_energy,
    hist_ElecInfo_electron_diphoton_deltaR,
    hist_MuonInfo_muon_pt,
    hist_MuonInfo_muon_energy,
    hist_MuonInfo_muon_diphoton_deltaR,
    //------------------------
    hist_JetInfo_jet_pt,
    hist_JetInfo_jet_energy,
    hist_JetInfo_jet_diphoton_deltaR,
    //------------------------
    //Define in selection stage
    //------------------------
    hist_lepton_pt,
    hist_lepton_energy,
    hist_lepton_diphoton_deltaR,
    //------------------------
    hist_jet1_diphoton_deltaR,
    hist_jet2_diphoton_deltaR,
    hist_jet1_lepton_deltaR,
    hist_jet2_lepton_deltaR,
    //------------------------
    hist_jet1_pt,
    hist_jet1_energy,
    hist_jet2_pt,
    hist_jet2_energy,
    //------------------------
    hist_leading_bjet_pt,
    hist_leading_bjet_energy,
    hist_chosen_bjet_pt,
    hist_chosen_bjet_energy,
    //------------------------
*/
