// vim: set fdm=marker:
#include <string>
#include <TCanvas.h>
#include <TFile.h>
#include <THStack.h>
#include <TH1D.h>
#include <TLine.h>
#include <TPad.h>
#include <string>
#include "../include/stack_2017old.h"
using namespace std;

bool considerQCD = false;
bool normalizeSigEntryToData = false;
bool bool_isHadronic;
bool bool_isLeptonic;
const double TunableSigBranchingFraction = 0.005; //The branching fraction of signal MC = 0.5%
//const double TunableSigBranchingFraction = 2.0; //The branching fraction of signal MC = 2001%
//const double TunableSigBranchingFraction = 0.001; //The branching fraction of signal MC = 0.1%
bool PrintTexStyle;

void stackHist_2017old(const char* channel){
    //### bool hadronic/leptonic channel{{{
    if((string)channel == "hadronic") bool_isHadronic = true; else bool_isHadronic = false;
    bool_isLeptonic = !bool_isHadronic;
    if(bool_isHadronic) printf("[CHECK] from macro src/stackHist.C: isHadronic!\n");
    if(bool_isLeptonic) printf("[CHECK] from macro src/stackHist.C: isLeptonic!\n");
    //}}}
    MakeStackHist("hist_DiPhoInfo_mass");
//    MakeStackHist("hist_num_leptons");
//    MakeStackHist("hist_num_jets");
//    MakeStackHist("hist_num_bjets");
//    MakeStackHist("hist_num_jets_leptonicChannel");
//    MakeStackHist("hist_num_bjets_leptonicChannel");
//    MakeStackHist("hist_num_jets_hadronicChannel");
//    MakeStackHist("hist_num_bjets_hadronicChannel");
    /*
    //### MakeStackHist{{{
    MakeStackHist("hist_EvtInfo_NPu");
    MakeStackHist("hist_EvtInfo_Rho");
    MakeStackHist("hist_EvtInfo_Rho_wopu");
    MakeStackHist("hist_EvtInfo_NVtx");
    MakeStackHist("hist_EvtInfo_NVtx_wopu");
    MakeStackHist("hist_EvtInfo_genweight");
    MakeStackHist("hist_DiPhoInfo_mass");
    MakeStackHist("hist_DiPhoInfo_pt");
    MakeStackHist("hist_DiPhoInfo_pt_overM");
    MakeStackHist("hist_DiPhoInfo_eta");
    MakeStackHist("hist_DiPhoInfo_phi");
    MakeStackHist("hist_DiPhoInfo_energy");
    MakeStackHist("hist_DiPhoInfo_leadPt");
    MakeStackHist("hist_DiPhoInfo_leadPt_overM");
    MakeStackHist("hist_DiPhoInfo_leadEta");
    MakeStackHist("hist_DiPhoInfo_leadPhi");
    MakeStackHist("hist_DiPhoInfo_leadE");
    MakeStackHist("hist_DiPhoInfo_leadhoe");
    MakeStackHist("hist_DiPhoInfo_leadIDMVA_beforeCut");
    MakeStackHist("hist_DiPhoInfo_leadhasPixelSeed_beforeCut");
    MakeStackHist("hist_DiPhoInfo_leadIDMVA");
    MakeStackHist("hist_DiPhoInfo_leadhasPixelSeed");
    MakeStackHist("hist_DiPhoInfo_subleadPt");
    MakeStackHist("hist_DiPhoInfo_subleadPt_overM");
    MakeStackHist("hist_DiPhoInfo_subleadEta");
    MakeStackHist("hist_DiPhoInfo_subleadPhi");
    MakeStackHist("hist_DiPhoInfo_subleadE");
    MakeStackHist("hist_DiPhoInfo_subleadhoe");
    MakeStackHist("hist_DiPhoInfo_subleadIDMVA_beforeCut");
    MakeStackHist("hist_DiPhoInfo_subleadhasPixelSeed_beforeCut");
    MakeStackHist("hist_DiPhoInfo_subleadIDMVA");
    MakeStackHist("hist_DiPhoInfo_subleadhasPixelSeed");
    MakeStackHist("hist_ElecInfo_Size");
    MakeStackHist("hist_MuonInfo_Size");
    MakeStackHist("hist_num_leptons");
    MakeStackHist("hist_num_electrons");
    MakeStackHist("hist_num_muons");
    MakeStackHist("hist_ElecInfo_electron_charge");
    MakeStackHist("hist_ElecInfo_electron_pt");
    MakeStackHist("hist_ElecInfo_electron_eta");
    MakeStackHist("hist_ElecInfo_electron_phi");
    MakeStackHist("hist_ElecInfo_electron_energy");
    MakeStackHist("hist_ElecInfo_electron_diphoton_deltaR");
    MakeStackHist("hist_MuonInfo_muon_charge");
    MakeStackHist("hist_MuonInfo_muon_pt");
    MakeStackHist("hist_MuonInfo_muon_eta");
    MakeStackHist("hist_MuonInfo_muon_phi");
    MakeStackHist("hist_MuonInfo_muon_energy");
    MakeStackHist("hist_MuonInfo_muon_diphoton_deltaR");
    MakeStackHist("hist_jets_size");
    MakeStackHist("hist_num_jets");
    MakeStackHist("hist_num_bjets");
    MakeStackHist("hist_num_jets_leptonicChannel");
    MakeStackHist("hist_num_bjets_leptonicChannel");
    MakeStackHist("hist_num_jets_hadronicChannel");
    MakeStackHist("hist_num_bjets_hadronicChannel");
    MakeStackHist("hist_JetInfo_jet_pt");
    MakeStackHist("hist_JetInfo_jet_eta");
    MakeStackHist("hist_JetInfo_jet_phi");
    MakeStackHist("hist_JetInfo_jet_energy");
    MakeStackHist("hist_JetInfo_jet_diphoton_deltaR");
    MakeStackHist("hist_MetInfo_Pt");
    MakeStackHist("hist_MetInfo_Phi");
    MakeStackHist("hist_MetInfo_Px");
    MakeStackHist("hist_MetInfo_Py");
    MakeStackHist("hist_MetInfo_SumET");
    MakeStackHist("hist_MetInfo_Pz_solution_1");
    MakeStackHist("hist_MetInfo_Pz_solution_2");
    MakeStackHist("hist_MetInfo_coeff_D");
    MakeStackHist("hist_MetInfo_coeff_A");
    MakeStackHist("hist_MetInfo_coeff_B2A");
    MakeStackHist("hist_MetInfo_coeff_D2A");
    MakeStackHist("hist_lepton_charge");
    MakeStackHist("hist_lepton_pt");
    MakeStackHist("hist_lepton_eta");
    MakeStackHist("hist_lepton_phi");
    MakeStackHist("hist_lepton_energy");
    MakeStackHist("hist_lepton_diphoton_deltaR");
    MakeStackHist("hist_jet1_pt");
    MakeStackHist("hist_jet1_eta");
    MakeStackHist("hist_jet1_phi");
    MakeStackHist("hist_jet1_energy");
    MakeStackHist("hist_jet1_btag_score");
    MakeStackHist("hist_jet1_CvsL_score");
    MakeStackHist("hist_jet1_CvsB_score");
    MakeStackHist("hist_jet1_diphoton_deltaR");
    MakeStackHist("hist_jet1_lepton_deltaR");
    MakeStackHist("hist_jet2_pt");
    MakeStackHist("hist_jet2_eta");
    MakeStackHist("hist_jet2_phi");
    MakeStackHist("hist_jet2_energy");
    MakeStackHist("hist_jet2_btag_score");
    MakeStackHist("hist_jet2_CvsL_score");
    MakeStackHist("hist_jet2_CvsB_score");
    MakeStackHist("hist_jet2_diphoton_deltaR");
    MakeStackHist("hist_jet2_lepton_deltaR");
    MakeStackHist("hist_wjet1_pt");
    MakeStackHist("hist_wjet1_eta");
    MakeStackHist("hist_wjet1_phi");
    MakeStackHist("hist_wjet1_energy");
    MakeStackHist("hist_wjet1_btag_score");
    MakeStackHist("hist_wjet1_CvsL_score");
    MakeStackHist("hist_wjet1_CvsB_score");
    MakeStackHist("hist_wjet1_diphoton_deltaR");
    MakeStackHist("hist_wjet1_lepton_deltaR");
    MakeStackHist("hist_wjet2_pt");
    MakeStackHist("hist_wjet2_eta");
    MakeStackHist("hist_wjet2_phi");
    MakeStackHist("hist_wjet2_energy");
    MakeStackHist("hist_wjet2_btag_score");
    MakeStackHist("hist_wjet2_CvsL_score");
    MakeStackHist("hist_wjet2_CvsB_score");
    MakeStackHist("hist_wjet2_diphoton_deltaR");
    MakeStackHist("hist_wjet2_lepton_deltaR");
    MakeStackHist("hist_jetq_pt");
    MakeStackHist("hist_jetq_eta");
    MakeStackHist("hist_jetq_phi");
    MakeStackHist("hist_jetq_energy");
    MakeStackHist("hist_jetq_btag_score");
    MakeStackHist("hist_jetq_CvsL_score");
    MakeStackHist("hist_jetq_CvsB_score");
    MakeStackHist("hist_jetq_diphoton_deltaR");
    MakeStackHist("hist_jetq_lepton_deltaR");
    MakeStackHist("hist_leading_bjet_pt");
    MakeStackHist("hist_leading_bjet_eta");
    MakeStackHist("hist_leading_bjet_phi");
    MakeStackHist("hist_leading_bjet_energy");
    MakeStackHist("hist_leading_bjet_btag_score");
    MakeStackHist("hist_leading_bjet_CvsL_score");
    MakeStackHist("hist_leading_bjet_CvsB_score");
    MakeStackHist("hist_deltaR_top_top");
    MakeStackHist("hist_deltaR_qH");
    MakeStackHist("hist_deltaR_photon_photon");
    MakeStackHist("hist_deltaR_bW");
    MakeStackHist("hist_deltaR_HW");
    MakeStackHist("hist_deltaR_tH");
    MakeStackHist("hist_deltaR_lep_met");
    MakeStackHist("hist_deltaR_jet1_jet2");
    MakeStackHist("hist_deltaR_wjet1_wjet2");
    MakeStackHist("hist_top_tqh_pt");
    MakeStackHist("hist_top_tqh_eta");
    MakeStackHist("hist_top_tqh_mass");
    MakeStackHist("hist_hadronic_w_candidate_pt");
    MakeStackHist("hist_hadronic_w_candidate_eta");
    MakeStackHist("hist_hadronic_w_candidate_mass");
    MakeStackHist("hist_hadronic_top_tbw_pt");
    MakeStackHist("hist_hadronic_top_tbw_eta");
    MakeStackHist("hist_hadronic_top_tbw_mass");
    MakeStackHist("hist_leptonic_w_candidate_solution1_pt");
    MakeStackHist("hist_leptonic_w_candidate_solution1_eta");
    MakeStackHist("hist_leptonic_w_candidate_solution1_mass");
    MakeStackHist("hist_leptonic_top_tbw_solution1_pt");
    MakeStackHist("hist_leptonic_top_tbw_solution1_eta");
    MakeStackHist("hist_leptonic_top_tbw_solution1_mass");
    MakeStackHist("hist_leptonic_w_candidate_solution2_pt");
    MakeStackHist("hist_leptonic_w_candidate_solution2_eta");
    MakeStackHist("hist_leptonic_w_candidate_solution2_mass");
    MakeStackHist("hist_leptonic_top_tbw_solution2_pt");
    MakeStackHist("hist_leptonic_top_tbw_solution2_eta");
    MakeStackHist("hist_leptonic_top_tbw_solution2_mass");
    MakeStackHist("hist_leptonic_w_candidate_topKinFit_pt");
    MakeStackHist("hist_leptonic_w_candidate_topKinFit_eta");
    MakeStackHist("hist_leptonic_w_candidate_topKinFit_mass");
    MakeStackHist("hist_leptonic_top_tbw_topKinFit_pt");
    MakeStackHist("hist_leptonic_top_tbw_topKinFit_eta");
    MakeStackHist("hist_leptonic_top_tbw_topKinFit_mass");
    //end of MakeStackHist}}}
    */
}
void MakeStackHist(const char* histName){
    //if((string)histName == "hist_EvtInfo_NVtx") PrintTexStyle = true; else PrintTexStyle = false;
    if((string)histName == "hist_DiPhoInfo_mass") PrintTexStyle = true; else PrintTexStyle = false;
    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 800);
    THStack *stackHist = new THStack("stackHist", "");//If setting titles here, the x(y)-title will NOT be able to set later.
    //### Set up boolean{{{
    bool isIDMVA = isThisIDMVA(histName);
    bool isNum = isThisNum(histName);
    bool isEtaPhi = isThisEtaPhi(histName);
    bool isMassSpectrum = isThisMassSpectrum(histName);
    bool isDijetSpectrum = isThisDijetSpectrum(histName);
    bool isDiPhotonSpectrum = isThisDiPhotonSpectrum(histName);
    bool isTopSpectrum = isThisTopSpectrum(histName);
    //}}}
    //### Hists{{{
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
    TH1D  *hist_tqh_Higgs;
    //--------------------
    TH1D  *hist_tqh_DiPhotonJetsBox;
    TH1D  *hist_tqh_GJet;
    TH1D  *hist_tqh_QCD;
    TH1D  *hist_tqh_TGJets;
    TH1D  *hist_tqh_TTGG;
    TH1D  *hist_tqh_TTGJets;
    TH1D  *hist_tqh_TTJets;
    TH1D  *hist_tqh_DYJetsToLL;
    TH1D  *hist_tqh_ZGToLLG;
    TH1D  *hist_tqh_WGToLNuG;
    TH1D  *hist_tqh_WJetsToLNu;
    TH1D  *hist_tqh_WW;
    TH1D  *hist_tqh_WZTo2L2Q;
    TH1D  *hist_tqh_ZZTo2L2Q;
    //--------------------
    TH1D  *hist_tqh_VG;//ZGToLLG, WGToLNuG
    TH1D  *hist_tqh_VV;//WZTo2L2Q, ZZTo2L2Q, WW
    //--------------------
    TH1D  *hist_tqh_data[NUM_data+1];
    TH1D  *hist_tqh_mc_wosig; //1.Calculate data/bkg ratio. 2.Statistical bkg uncertainty.
    TH1D  *hist_tqh_mc_uncertainty; //1.Statistical MC uncertainty.
    TH1D  *hist_tqh_ratio;
    //}}}
    //### Register histograms{{{
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
    //}}}
    if((string)histName == "hist_DiPhoInfo_mass") printf("[check] [%6.2f, %6.2f]\n", hist_tqh_sig[0]->GetBinLowEdge(5), hist_tqh_sig[0]->GetBinLowEdge(7));
    //if((string)histName == "hist_DiPhoInfo_mass") printf("[check] [%6.2f, %6.2f]\n", hist_tqh_sig[0]->GetBinLowEdge(11), hist_tqh_sig[0]->GetBinLowEdge(16));
    //### Combine had/lep signal{{{ 
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
    //}}}
    //### Combine backgournds {{{
    //===============================//
    //===== Combine backgournds =====//
    //===============================//
    hist_tqh_resbkg[NUM_resbkg] = (TH1D*) hist_tqh_resbkg[0]->Clone();
    for(int i=1; i<NUM_resbkg; i++) hist_tqh_resbkg[NUM_resbkg]->Add(hist_tqh_resbkg[i]);
    hist_tqh_nonresbkg[NUM_nonresbkg] = (TH1D*) hist_tqh_nonresbkg[0]->Clone();
    if(considerQCD)
    {
        for(int i=1; i<NUM_nonresbkg; i++) hist_tqh_nonresbkg[NUM_nonresbkg]->Add(hist_tqh_nonresbkg[i]);
    }
    else
    {
        for(int i=1; i<NUM_nonresbkg; i++)
        {
            if(i==8 || i==9 || i==10) continue;
            else hist_tqh_nonresbkg[NUM_nonresbkg]->Add(hist_tqh_nonresbkg[i]);
        }
    }
    //--------------------
    hist_tqh_ggH = (TH1D*) hist_tqh_resbkg[0]->Clone();
    hist_tqh_VBF = (TH1D*) hist_tqh_resbkg[1]->Clone();
    hist_tqh_VH  = (TH1D*) hist_tqh_resbkg[2]->Clone();
    hist_tqh_ttH = (TH1D*) hist_tqh_resbkg[3]->Clone();
    hist_tqh_ggH -> SetFillColor(kYellow);
    hist_tqh_ggH -> SetLineColor(kYellow);
    hist_tqh_VBF -> SetFillColor(kYellow);
    hist_tqh_VBF -> SetLineColor(kYellow);
    hist_tqh_VH  -> SetFillColor(kYellow);
    hist_tqh_VH  -> SetLineColor(kYellow);
    hist_tqh_ttH -> SetFillColor(kYellow);
    hist_tqh_ttH -> SetLineColor(kYellow);
    //--------------------
    hist_tqh_Higgs = (TH1D*) hist_tqh_ggH->Clone();
    hist_tqh_Higgs -> Add(hist_tqh_VBF);
    hist_tqh_Higgs -> Add(hist_tqh_VH);
    hist_tqh_Higgs -> Add(hist_tqh_ttH);
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
    hist_tqh_TTJets = (TH1D*) hist_tqh_nonresbkg[11];
    hist_tqh_DYJetsToLL = (TH1D*) hist_tqh_nonresbkg[12];
    hist_tqh_ZGToLLG = (TH1D*) hist_tqh_nonresbkg[13];
    hist_tqh_WGToLNuG = (TH1D*) hist_tqh_nonresbkg[14];
    hist_tqh_WJetsToLNu = (TH1D*) hist_tqh_nonresbkg[15];
    hist_tqh_WW = (TH1D*) hist_tqh_nonresbkg[16];
    hist_tqh_WZTo2L2Q = (TH1D*) hist_tqh_nonresbkg[17];
    hist_tqh_ZZTo2L2Q = (TH1D*) hist_tqh_nonresbkg[18];
    //--------------------
    hist_tqh_VG = (TH1D*) hist_tqh_ZGToLLG->Clone();
    hist_tqh_VG -> Add(hist_tqh_WGToLNuG);
    hist_tqh_VV = (TH1D*) hist_tqh_WZTo2L2Q->Clone();
    hist_tqh_VV -> Add(hist_tqh_ZZTo2L2Q);
    hist_tqh_VV -> Add(hist_tqh_WW);
    //--------------------
    hist_tqh_DiPhotonJetsBox->SetFillColor(kOrange+1);
    hist_tqh_DiPhotonJetsBox->SetLineColor(kOrange+1);
    hist_tqh_GJet->SetFillColor(kOrange+2);
    hist_tqh_GJet->SetLineColor(kOrange+2);
    hist_tqh_TGJets->SetFillColor(kOrange+3);
    hist_tqh_TGJets->SetLineColor(kOrange+3);
    hist_tqh_TTGJets->SetFillColor(kGreen+1);
    hist_tqh_TTGJets->SetLineColor(kGreen+1);
    hist_tqh_TTGG->SetFillColor(kGreen+3);
    hist_tqh_TTGG->SetLineColor(kGreen+3);
    hist_tqh_TTJets -> SetFillColor(kViolet+1);
    hist_tqh_TTJets -> SetLineColor(kViolet+1);
    hist_tqh_DYJetsToLL -> SetFillColor(kBlue);
    hist_tqh_DYJetsToLL -> SetLineColor(kBlue);
    hist_tqh_WJetsToLNu -> SetFillColor(kCyan);
    hist_tqh_WJetsToLNu -> SetLineColor(kCyan);
    hist_tqh_VG -> SetFillColor(kBlue-7);
    hist_tqh_VG -> SetLineColor(kBlue-7);
    hist_tqh_VV -> SetFillColor(kCyan-10);
    hist_tqh_VV -> SetLineColor(kCyan-10);
    hist_tqh_QCD->SetFillColor(kOrange);
    hist_tqh_QCD->SetLineColor(kOrange);
    //}}}
    //### Combine data{{{ 
    //========================//
    //===== Combine data =====//
    //========================//
    hist_tqh_data[NUM_data] = (TH1D*) hist_tqh_data[0]->Clone();
    for(int i=1; i<NUM_data; i++) hist_tqh_data[NUM_data]->Add(hist_tqh_data[i]);
    //}}}
    //### Combine mc w/o sig{{{ 
    //==============================//
    //===== Combine mc w/o sig =====//
    //==============================//
    hist_tqh_mc_wosig = (TH1D*) hist_tqh_resbkg[NUM_resbkg]->Clone();
    hist_tqh_mc_wosig->Add(hist_tqh_nonresbkg[NUM_nonresbkg]);
    //}}}
    //### Calculate relative MC uncertainty {{{
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
    //}}}
    //### Make data-mc ratio histogram {{{
    //========================================//
    //===== Make data-mc ratio histogram =====//
    //========================================//
    hist_tqh_ratio = (TH1D*) hist_tqh_data[NUM_data]->Clone();
    hist_tqh_ratio->Divide(hist_tqh_mc_wosig);
    //}}}
    //### Make stack histograms {{{
    //=================================//
    //===== Make stack histograms =====//
    //=================================//
    stackHist->Add(hist_tqh_Higgs);
    stackHist->Add(hist_tqh_DiPhotonJetsBox);
    stackHist->Add(hist_tqh_GJet);
    stackHist->Add(hist_tqh_TGJets);
    stackHist->Add(hist_tqh_TTGJets);
    stackHist->Add(hist_tqh_TTGG);
    stackHist->Add(hist_tqh_TTJets );
    stackHist->Add(hist_tqh_DYJetsToLL );
    stackHist->Add(hist_tqh_WJetsToLNu );
    stackHist->Add(hist_tqh_VG );
    stackHist->Add(hist_tqh_VV );
    if(considerQCD) stackHist->Add(hist_tqh_QCD);
    //}}}
    //### Choose proper signal hist to present{{{ 
    //================================================//
    //===== Choose proper signal hist to present =====//
    //================================================//
    TH1D *hist_sig_tt_hut, *hist_sig_st_hut;
    TH1D *hist_sig_tt_hct, *hist_sig_st_hct;
    //--------------------
    if(bool_isHadronic){
        hist_sig_tt_hut = (TH1D*) hist_tqh_sig_tt_hadronic_hut->Clone(); 
        hist_sig_tt_hct = (TH1D*) hist_tqh_sig_tt_hadronic_hct->Clone(); 
        hist_sig_st_hut = (TH1D*) hist_tqh_sig_st_hadronic_hut->Clone(); 
        hist_sig_st_hct = (TH1D*) hist_tqh_sig_st_hadronic_hct->Clone(); 
    } else{
        hist_sig_tt_hut = (TH1D*) hist_tqh_sig_tt_leptonic_hut->Clone(); 
        hist_sig_tt_hct = (TH1D*) hist_tqh_sig_tt_leptonic_hct->Clone(); 
        hist_sig_st_hut = (TH1D*) hist_tqh_sig_st_leptonic_hut->Clone(); 
        hist_sig_st_hct = (TH1D*) hist_tqh_sig_st_leptonic_hct->Clone(); 
    }
    //--------------------
    hist_sig_tt_hut->SetLineColor(kPink+1);
    hist_sig_st_hut->SetLineColor(kPink+1); hist_sig_st_hut->SetLineStyle(2);
    hist_sig_tt_hct->SetLineColor(kRed);
    hist_sig_st_hct->SetLineColor(kRed); hist_sig_st_hct->SetLineStyle(2);
    //}}}
    /*
    //### Calculate yields from combined histograms (signal region){{{
    //=====================================================//
    //===== Calculate yields from combined histograms =====//
    //=====================================================//
    if(PrintTexStyle){
        printf("Processes & Entries & \\multicolumn{2}{c}{Yields}\\\\ \n");
        printf("\\hline\\hline\n");
        CalculateHistYields_signalRegion("DiPhotonJetsBox ", hist_tqh_DiPhotonJetsBox);
        CalculateHistYields_signalRegion("GJet\t\t", hist_tqh_GJet);
        CalculateHistYields_signalRegion("TGJets\t\t", hist_tqh_TGJets);
        CalculateHistYields_signalRegion("TTGJets\t\t", hist_tqh_TTGJets);
        CalculateHistYields_signalRegion("TTGG\t\t", hist_tqh_TTGG);
        CalculateHistYields_signalRegion("TTJets\t\t", hist_tqh_TTJets );
        CalculateHistYields_signalRegion("DYJetsToLL\t", hist_tqh_DYJetsToLL );
        CalculateHistYields_signalRegion("WJetsToLNu\t", hist_tqh_WJetsToLNu );
        CalculateHistYields_signalRegion("VG\t\t", hist_tqh_VG );
        CalculateHistYields_signalRegion("VV\t\t", hist_tqh_VV );
        CalculateHistYields_signalRegion("Higgs\t\t", hist_tqh_Higgs);
        CalculateHistYields_signalRegion("ggH\t\t", hist_tqh_ggH);
        CalculateHistYields_signalRegion("VBF\t\t", hist_tqh_VBF);
        CalculateHistYields_signalRegion("VH\t\t", hist_tqh_VH);
        CalculateHistYields_signalRegion("ttH\t\t", hist_tqh_ttH);
        CalculateHistYields_signalRegion("QCD\t\t", hist_tqh_QCD);
        CalculateHistYields_signalRegion("sig_tt_hut\t", hist_sig_tt_hut);
        CalculateHistYields_signalRegion("sig_st_hut\t", hist_sig_st_hut);
        CalculateHistYields_signalRegion("sig_tt_hct\t", hist_sig_tt_hct);
        CalculateHistYields_signalRegion("sig_st_hct\t", hist_sig_st_hct);
        printf("\\hline\n");
        CalculateHistYields_signalRegion("MC background\t", hist_tqh_mc_wosig);
        CalculateHistYields_signalRegion("Data\t\t", hist_tqh_data[NUM_data]);
        printf("\n\n");
    }
    //}}}
    //### Calculate yields from combined histograms (full region){{{
    //=====================================================//
    //===== Calculate yields from combined histograms =====//
    //=====================================================//
    if(PrintTexStyle){
        printf("Processes & Entries & \\multicolumn{2}{c}{Yields}\\\\ \n");
        printf("\\hline\\hline\n");
        CalculateHistYields_fullRegion("DiPhotonJetsBox ", hist_tqh_DiPhotonJetsBox);
        CalculateHistYields_fullRegion("GJet\t\t", hist_tqh_GJet);
        CalculateHistYields_fullRegion("TGJets\t\t", hist_tqh_TGJets);
        CalculateHistYields_fullRegion("TTGJets\t\t", hist_tqh_TTGJets);
        CalculateHistYields_fullRegion("TTGG\t\t", hist_tqh_TTGG);
        CalculateHistYields_fullRegion("TTJets\t\t", hist_tqh_TTJets );
        CalculateHistYields_fullRegion("DYJetsToLL\t", hist_tqh_DYJetsToLL );
        CalculateHistYields_fullRegion("WJetsToLNu\t", hist_tqh_WJetsToLNu );
        CalculateHistYields_fullRegion("VG\t\t", hist_tqh_VG );
        CalculateHistYields_fullRegion("VV\t\t", hist_tqh_VV );
        CalculateHistYields_fullRegion("Higgs\t\t", hist_tqh_Higgs);
        CalculateHistYields_fullRegion("ggH\t\t", hist_tqh_ggH);
        CalculateHistYields_fullRegion("VBF\t\t", hist_tqh_VBF);
        CalculateHistYields_fullRegion("VH\t\t", hist_tqh_VH);
        CalculateHistYields_fullRegion("ttH\t\t", hist_tqh_ttH);
        CalculateHistYields_fullRegion("QCD\t\t", hist_tqh_QCD);
        CalculateHistYields_fullRegion("sig_tt_hut\t", hist_sig_tt_hut);
        CalculateHistYields_fullRegion("sig_st_hut\t", hist_sig_st_hut);
        CalculateHistYields_fullRegion("sig_tt_hct\t", hist_sig_tt_hct);
        CalculateHistYields_fullRegion("sig_st_hct\t", hist_sig_st_hct);
        printf("\\hline\n");
        CalculateHistYields_fullRegion("MC background\t", hist_tqh_mc_wosig);
        CalculateHistYields_fullRegion("Data\t\t", hist_tqh_data[NUM_data]);
        printf("\n\n");
    }
    //}}}
    */
    //### Calculate yields from combined histograms (sideband region){{{
    //=====================================================//
    //===== Calculate yields from combined histograms =====//
    //=====================================================//
    if(PrintTexStyle){
        printf("Processes & Entries & \\multicolumn{2}{c}{Yields}\\\\ \n");
        printf("\\hline\\hline\n");
        CalculateHistYields_sidebandRegion("DiPhotonJetsBox ", hist_tqh_DiPhotonJetsBox);
        CalculateHistYields_sidebandRegion("GJet\t\t", hist_tqh_GJet);
        CalculateHistYields_sidebandRegion("TGJets\t\t", hist_tqh_TGJets);
        CalculateHistYields_sidebandRegion("TTGJets\t\t", hist_tqh_TTGJets);
        CalculateHistYields_sidebandRegion("TTGG\t\t", hist_tqh_TTGG);
        CalculateHistYields_sidebandRegion("TTJets\t\t", hist_tqh_TTJets );
        CalculateHistYields_sidebandRegion("DYJetsToLL\t", hist_tqh_DYJetsToLL );
        CalculateHistYields_sidebandRegion("WJetsToLNu\t", hist_tqh_WJetsToLNu );
        CalculateHistYields_sidebandRegion("VG\t\t", hist_tqh_VG );
        CalculateHistYields_sidebandRegion("VV\t\t", hist_tqh_VV );
        CalculateHistYields_sidebandRegion("Higgs\t\t", hist_tqh_Higgs);
        CalculateHistYields_sidebandRegion("ggH\t\t", hist_tqh_ggH);
        CalculateHistYields_sidebandRegion("VBF\t\t", hist_tqh_VBF);
        CalculateHistYields_sidebandRegion("VH\t\t", hist_tqh_VH);
        CalculateHistYields_sidebandRegion("ttH\t\t", hist_tqh_ttH);
        CalculateHistYields_sidebandRegion("QCD\t\t", hist_tqh_QCD);
        CalculateHistYields_sidebandRegion("sig_tt_hut\t", hist_sig_tt_hut);
        CalculateHistYields_sidebandRegion("sig_st_hut\t", hist_sig_st_hut);
        CalculateHistYields_sidebandRegion("sig_tt_hct\t", hist_sig_tt_hct);
        CalculateHistYields_sidebandRegion("sig_st_hct\t", hist_sig_st_hct);
        printf("\\hline\n");
        CalculateHistYields_sidebandRegion("MC background\t", hist_tqh_mc_wosig);
        CalculateHistYields_sidebandRegion("Data\t\t", hist_tqh_data[NUM_data]);
        printf("\n\n");
    }
    //}}}
    // Scale signal{{{
    //--- normalize to total entries of data ---//
    double total_entries_data = 10.;
    double scale_sig_tt_hut = 10.;
    double scale_sig_st_hut = 10.;
    double scale_sig_tt_hct = 10.;
    double scale_sig_st_hct = 10.;
    if(normalizeSigEntryToData){
        total_entries_data = hist_tqh_data[NUM_data]->Integral();
        scale_sig_tt_hut = total_entries_data / hist_sig_tt_hut->Integral() ; hist_sig_tt_hut->Scale(scale_sig_tt_hut);
        scale_sig_st_hut = total_entries_data / hist_sig_st_hut->Integral() ; hist_sig_st_hut->Scale(scale_sig_st_hut);
        scale_sig_tt_hct = total_entries_data / hist_sig_tt_hct->Integral() ; hist_sig_tt_hct->Scale(scale_sig_tt_hct);
        scale_sig_st_hct = total_entries_data / hist_sig_st_hct->Integral() ; hist_sig_st_hct->Scale(scale_sig_st_hct);
    } else{
        hist_sig_tt_hut->Scale(scale_sig_tt_hut);
        hist_sig_st_hut->Scale(scale_sig_st_hut);
        hist_sig_tt_hct->Scale(scale_sig_tt_hct);
        hist_sig_st_hct->Scale(scale_sig_st_hct);
    }
    //}}}
    //### Draw upper plots{{{ 
    //============================//
    //===== Draw upper plots =====//
    //============================//
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
    pad1->SetRightMargin(0.22);
    pad1->SetLeftMargin(0.10);
    pad1->SetBottomMargin(0.015); //Upper and lower pads are joined
    pad1->SetAttLinePS(kBlack,1,2);
    pad1->Draw();
    pad1->cd(); //pad1 becomes current pad
    gPad->SetTicks(1,1);
    gPad->SetGrid();
    //--------------------
    hist_tqh_mc_wosig->SetLineWidth(0);
    hist_tqh_mc_wosig->SetMarkerColor(kGray+2);
    hist_tqh_mc_wosig->SetFillColor(kGray+2);
    hist_tqh_mc_wosig->SetFillStyle(3001);
    //--------------------
    stackHist->SetMinimum(0);
    double scale = 1.2;

    std::vector<TH1D*> histograms_plotted;
    histograms_plotted.push_back(hist_sig_tt_hut);
    histograms_plotted.push_back(hist_sig_st_hut);
    histograms_plotted.push_back(hist_sig_tt_hct);
    histograms_plotted.push_back(hist_sig_st_hct);
    histograms_plotted.push_back(hist_tqh_data[NUM_data]);
    histograms_plotted.push_back(hist_tqh_mc_wosig);
    double max_scope = GetMaxScope( stackHist, histograms_plotted );

    stackHist->Draw("hist");
    hist_sig_tt_hut->Draw("hist,same");
    hist_sig_st_hut->Draw("hist,same");
    hist_sig_tt_hct->Draw("hist,same");
    hist_sig_st_hct->Draw("hist,same");
    hist_tqh_mc_wosig->Draw("E2,same");
    hist_tqh_data[NUM_data]->Draw("p,E1,same");
    stackHist->SetMaximum(max_scope*scale);
    //stackHist->SetMaximum(15000);

    //##### SetTitles{{{
    //--------------------
    double BinWidth = hist_tqh_data[NUM_data]->GetXaxis()->GetBinWidth(1);//Take the width of the first bin as a representative.
    string yTitle = GetYtitleAccordingToHistName(histName, BinWidth);
    stackHist->GetYaxis()->SetTitle(yTitle.c_str());
    stackHist->GetYaxis()->SetTitleSize(32);//20
    stackHist->GetYaxis()->SetTitleFont(43);
    stackHist->GetYaxis()->SetTitleOffset(1.05);
    stackHist->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    stackHist->GetYaxis()->SetLabelSize(20);//15
    stackHist->GetXaxis()->SetLabelSize(0);// remove labels on x-axis
    //}}}
    //##### TLegend{{{
    //# normal{{{
    TLegend *legend = new TLegend(0.79,0.15,0.85,0.85);
    legend->SetTextSize(0.03);
    legend->AddEntry(hist_tqh_data[NUM_data], "Observed", "lep");
    legend->AddEntry(hist_tqh_DiPhotonJetsBox, "#gamma#gamma + jets", "f");
    legend->AddEntry(hist_tqh_GJet, "#gamma + jet", "f");
    legend->AddEntry(hist_tqh_TGJets, "t#gamma + jets", "f");
    legend->AddEntry(hist_tqh_TTGJets, "t#bar{t}#gamma + jets", "f");
    legend->AddEntry(hist_tqh_TTGG, "t#bar{t}#gamma#gamma", "f");
    legend->AddEntry(hist_tqh_TTJets, "t#bar{t} + jets", "f");
    legend->AddEntry(hist_tqh_DYJetsToLL, "DY + jets", "f");
    legend->AddEntry(hist_tqh_WJetsToLNu, "W + jets", "f");
    legend->AddEntry(hist_tqh_VG, "V#gamma", "f");
    legend->AddEntry(hist_tqh_VV, "VV", "f");
    legend->AddEntry(hist_tqh_Higgs, "Higgs", "f");
    if(considerQCD) legend->AddEntry(hist_tqh_QCD, "QCD", "f");
    legend->AddEntry(hist_tqh_mc_wosig, "Bkg uncertainty", "f");
    //already determined which channel
    legend->AddEntry(hist_sig_tt_hut, Form("TT Hut #times %.0f (BF=%.2f%%)", scale_sig_tt_hut, TunableSigBranchingFraction*100), "f");
    legend->AddEntry(hist_sig_st_hut, Form("ST Hut #times %.0f (BF=%.2f%%)", scale_sig_st_hut, TunableSigBranchingFraction*100), "f");
    legend->AddEntry(hist_sig_tt_hct, Form("TT Hct #times %.0f (BF=%.2f%%)", scale_sig_tt_hct, TunableSigBranchingFraction*100), "f");
    legend->AddEntry(hist_sig_st_hct, Form("ST Hct #times %.0f (BF=%.2f%%)", scale_sig_st_hct, TunableSigBranchingFraction*100), "f");
    legend->SetLineColor(0);
    legend->Draw("same");
    //}}}
    //--------------------
    /*
    // # poster {{{
    TLegend *legend = new TLegend(0.79,0.05,0.88,0.85);
    legend->SetTextSize(0.05);
    legend->AddEntry(hist_tqh_data[NUM_data], "Observed", "lep");
    legend->AddEntry(hist_tqh_mc_wosig, "Uncertainty", "f");
    legend->AddEntry(hist_sig_tt_hut, "TT Hut", "f");
    legend->AddEntry(hist_sig_st_hut, "ST Hut", "f");
    legend->AddEntry(hist_sig_tt_hct, "TT Hct", "f");
    legend->AddEntry(hist_sig_st_hct, "ST Hct", "f");
    legend->SetLineColor(0);
    legend->Draw("same");
    //}}}
    // # poster lep{{{
    TLegend *legend_poster_lep = new TLegend(0.40,0.43,0.73,0.75);
    legend_poster_lep->SetNColumns(2);
    legend_poster_lep->SetTextSize(0.06);//0.03
    legend_poster_lep->SetFillStyle(4000);//0.03
    legend_poster_lep->SetFillColor(4000);//0.03
    legend_poster_lep->AddEntry(hist_tqh_VG, "V#gamma", "f");
    legend_poster_lep->AddEntry(hist_tqh_DYJetsToLL, "DY + jets", "f");
    legend_poster_lep->AddEntry(hist_tqh_TTJets, "t#bar{t} + jets", "f");
    legend_poster_lep->AddEntry(hist_tqh_TTGJets, "t#bar{t}#gamma + jets", "f");
    //legend_poster_lep->AddEntry(hist_tqh_TTGG, "t#bar{t}#gamma#gamma", "f");
    //legend_poster_lep->AddEntry(hist_tqh_WJetsToLNu, "W + jets", "f");
    //legend_poster_lep->AddEntry(hist_tqh_VV, "VV", "f");
    //legend_poster_lep->AddEntry(hist_tqh_Higgs, "Higgs", "f");
    legend_poster_lep->SetLineColor(0);
    legend_poster_lep->Draw("same");
    //}}}
    // # poster had{{{
    TLegend *legend_poster_had = new TLegend(0.40,0.43,0.73,0.75);
    legend_poster_had->SetNColumns(2);
    legend_poster_had->SetTextSize(0.06);//0.03
    legend_poster_had->SetFillStyle(4000);//0.03
    legend_poster_had->SetFillColor(4000);//0.03
    if(considerQCD) legend_poster_had->AddEntry(hist_tqh_QCD, "QCD", "f");
    legend_poster_had->AddEntry(hist_tqh_WJetsToLNu, "W + jets", "f");
    legend_poster_had->AddEntry(hist_tqh_GJet, "#gamma + jet", "f");
    legend_poster_had->AddEntry(hist_tqh_DiPhotonJetsBox, "#gamma#gamma + jets", "f");
    //legend_poster_had->AddEntry(hist_tqh_TGJets, "t#gamma + jets", "f");
    //legend_poster_had->AddEntry(hist_tqh_TTGJets, "t#bar{t}#gamma + jets", "f");
    //legend_poster_had->AddEntry(hist_tqh_TTGG, "t#bar{t}#gamma#gamma", "f");
    //legend_poster_had->AddEntry(hist_tqh_TTJets, "t#bar{t} + jets", "f");
    //legend_poster_had->AddEntry(hist_tqh_DYJetsToLL, "DY + jets", "f");
    //legend_poster_had->AddEntry(hist_tqh_VG, "V#gamma", "f");
    //legend_poster_had->AddEntry(hist_tqh_VV, "VV", "f");
    //legend_poster_had->AddEntry(hist_tqh_Higgs, "Higgs", "f");
    legend_poster_had->SetLineColor(0);
    legend_poster_had->Draw("same");
    //}}}
    */

    //}}}
    //##### TLatex{{{
    //--------------------
    TLatex latex;
    latex.SetNDC(kTRUE);
    latex.SetTextFont(43);
    latex.SetTextSize(32);//22
    latex.SetTextAlign(11);
    latex.DrawLatex(0.10, 0.92, "#bf{CMS} #it{Preliminary}");
    //latex.DrawLatex(0.14, 0.84, "#bf{CMS}");
    //latex.DrawLatex(0.14, 0.78, "#it{Preliminary}");
    //--------------------
    TLatex latex_channel;
    latex_channel.SetNDC(kTRUE);
    latex_channel.SetTextFont(43);
    latex_channel.SetTextSize(40);//26
    latex_channel.SetTextAlign(31);
    if(bool_isHadronic) latex_channel.DrawLatex(0.72, 0.80, "#bf{Hadronic Channel}");
    if(bool_isLeptonic) latex_channel.DrawLatex(0.72, 0.80, "#bf{Leptonic Channel}");
    //--------------------
    TLatex latex_lumi;
    latex_lumi.SetNDC(kTRUE);
    latex_lumi.SetTextFont(43);
    latex_lumi.SetTextSize(32);//22
    latex_lumi.SetTextAlign(31);
    latex_lumi.DrawLatex(0.77, 0.92, "41.5 fb^{-1} (2017, #sqrt{s} = 13 TeV)");
    //latex_lumi.DrawLatex(0.77, 0.92, "137.2 fb^{-1} (#sqrt{s} = 13 TeV)");
    //}}}
    //}}}
    //### Draw lower plots{{{ 
    //============================//
    //===== Draw lower plots =====//
    //============================//
    c1->cd(); //Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.245);//0.25
    pad2->SetRightMargin(0.22);
    pad1->SetLeftMargin(0.10);
    pad2->SetTopMargin(0.036);//0.01 //Upper and lower pads are joined
    pad2->SetBottomMargin(0.5);
    pad2->SetAttLinePS(kBlack,1,2);
    pad2->SetGridx(1);
    pad2->Draw();
    pad2->cd(); //pad2 becomes current pad
    gPad->SetTicks(1,1);
    gPad->SetGrid();
    //--------------------
    hist_tqh_ratio->SetMaximum(2.0);
    //hist_tqh_ratio->SetMinimum(0);
    //hist_tqh_ratio->SetMaximum(1.75);
    hist_tqh_ratio->SetMinimum(0.);
    //hist_tqh_ratio->SetMaximum(2.25);
    hist_tqh_ratio->SetTitle("");
    hist_tqh_ratio->SetStats(0); //No statistics on lower plot
    hist_tqh_ratio->Draw("p,E1");
    hist_tqh_mc_uncertainty->SetFillColor(kGray);
    hist_tqh_mc_uncertainty->SetLineColor(kGray);
    hist_tqh_mc_uncertainty->Draw("E2,same");
    hist_tqh_ratio->Draw("p,E1,same");
    //--------------------
    hist_tqh_ratio->GetYaxis()->SetTitle(" Obs/Exp");
    //hist_tqh_ratio->GetYaxis()->SetNdivisions(10);
    hist_tqh_ratio->GetYaxis()->SetNdivisions(5);
    hist_tqh_ratio->GetYaxis()->SetTitleSize(32);//20
    hist_tqh_ratio->GetYaxis()->SetTitleFont(43);
    hist_tqh_ratio->GetYaxis()->SetTitleOffset(1.0);//1.2
    hist_tqh_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hist_tqh_ratio->GetYaxis()->SetLabelSize(20);
    //--------------------
    //hist_tqh_ratio->GetXaxis()->SetTitle("Invariant mass of diphoton + jet [GeV/c^{2}]");
    string xTitle = GetXtitleAccordingToHistName(histName);
    hist_tqh_ratio->GetXaxis()->SetTitle(xTitle.c_str());
    hist_tqh_ratio->GetXaxis()->SetTitleSize(40);//25
    hist_tqh_ratio->GetXaxis()->SetTitleFont(43);
    hist_tqh_ratio->GetXaxis()->SetTitleOffset(3.3);
    hist_tqh_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hist_tqh_ratio->GetXaxis()->SetLabelSize(20);//15
    hist_tqh_ratio->GetXaxis()->SetLabelOffset(0.03);
    //--------------------
    c1->Update();//update the value of pad2->GetUxmax().
    TLine line;
    line.SetLineStyle(1);
    line.DrawLine(pad2->GetUxmin(),1.0,pad2->GetUxmax(),1.0);
    //line.SetLineStyle(2);
    //line.DrawLine(pad2->GetUxmin(),0.5,pad2->GetUxmax(),0.5);
    //line.DrawLine(pad2->GetUxmin(),1.0,pad2->GetUxmax(),1.0);
    //line.DrawLine(pad2->GetUxmin(),1.5,pad2->GetUxmax(),1.5);
    //line.DrawLine(pad2->GetUxmin(),2.0,pad2->GetUxmax(),2.0);
    c1->SaveAs(Form("%s/stack_%s.pdf", TARGET_DIR, histName));
    //c1->SaveAs(Form("%s/stack_%s.root", TARGET_DIR, histName));
    //}}}
    //### Draw log scale{{{ 
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
        c1->SaveAs(Form("%s/log_scale/stack_%s_log.pdf", TARGET_DIR, histName));
    }
    else stackHist->SetMaximum(isNum ? 1e+9 : 1e+9);
    //}}}

    for(int i=0; i<NUM; i++) file[i]->Close();
}

//### bool functions{{{
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
//}}}
//### title functions{{{
string GetXtitleAccordingToHistName(const char* histName){
    if((string)histName == "hist_EvtInfo_NPu") return "Number of pile up";
    if((string)histName == "hist_EvtInfo_Rho") return "Rho";
    if((string)histName == "hist_EvtInfo_Rho_wopu") return "Rho (w/o pu)";
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
    if((string)histName == "hist_DiPhoInfo_mass") return "Diphoton invariant mass [GeV/c^{2}]";
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
//}}}
//### other functions{{{
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

//to be designed later
void CalculateHistYields_fullRegion(const char *process, TH1D* hist){
    // signal region
    int totalEntries = hist->GetEntries();
    double totalYields = 0.;
    double totalError = 0.;
    for(int i=1; i<41; ++i){
        totalYields +=  hist->GetBinContent(i);
        double error =  hist->GetBinError(i);
        totalError += pow(error, 2);
    }
    totalError = sqrt(totalError);

    if(PrintTexStyle) printf("%s & \t %15d & \t %15.2f & $\\pm$ \t %10.2f\\\\\n", process, totalEntries, totalYields, totalError);
    else printf("%s \t %15.2f \t %15.2f\n", process, totalYields, totalError);
}

void CalculateHistYields_signalRegion(const char *process, TH1D* hist){
    // signal region
    int totalEntries = hist->GetEntries();
    double totalYields = 0.;
    double totalError = 0.;
    for(int i=11; i<16; ++i){
        totalYields +=  hist->GetBinContent(i);
        double error =  hist->GetBinError(i);
        totalError += pow(error, 2);
    }
    totalError = sqrt(totalError);

    if(PrintTexStyle) printf("%s & \t %15d & \t %15.2f & $\\pm$ \t %10.2f\\\\\n", process, totalEntries, totalYields, totalError);
    else printf("%s \t %15.2f \t %15.2f\n", process, totalYields, totalError);
}

void CalculateHistYields_sidebandRegion(const char *process, TH1D* hist){
    // signal region
    int totalEntries = hist->GetEntries();
    double totalYields = 0.;
    double totalError = 0.;
    for(int i=1; i<41; ++i){
        if( i>5 && i<7) continue; // exclude [120, 130]
        //if( i>10 && i<16) continue; // exclude [120, 130]
        totalYields +=  hist->GetBinContent(i);
        double error =  hist->GetBinError(i);
        totalError += pow(error, 2);
    }
    totalError = sqrt(totalError);

    if(PrintTexStyle) printf("%s & \t %15d & \t %15.2f & $\\pm$ \t %10.2f\\\\\n", process, totalEntries, totalYields, totalError);
    else printf("%s \t %15.2f \t %15.2f\n", process, totalYields, totalError);

    //--------------------------------------------------
    // whole region
    /*
    int totalEntries = hist->GetEntries();
    double totalYields = hist->Integral();
    double totalError = SumErrors(hist);
    if(PrintTexStyle) printf("%s & \t %15d & \t %15.2f & $\\pm$ \t %10.2f\\\\\n", process, totalEntries, totalYields, totalError);
    else printf("%s \t %15.2f \t %15.2f\n", process, totalYields, totalError);
    */
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

double GetMaxScope(THStack* h1, std::vector<TH1D*> vec_hists){
    std::vector<double> values;
    values.push_back(h1->GetMaximum());
    for(int i=0; i<vec_hists.size(); ++i) values.push_back(vec_hists[i]->GetMaximum());
    TH1D* h = vec_hists[vec_hists.size()-1];//hist_tqh_mc_wosig
    values.push_back(  h->GetMaximum() + h->GetBinError(h->GetMaximumBin())  );
    double max = *std::max_element(values.begin(), values.end());

    //double values[] = {h1->GetMaximum(), h2->GetMaximum(), h3->GetMaximum(), h4->GetMaximum() + h4->GetBinError(h4->GetMaximumBin()), h5->GetMaximum()};
    //double max = *std::max_element(std::begin(values), std::end(values));
    ////for(int i=0; i<5; i++) printf("[check-scope] value_%d = %.2f\n", i, values[i]);
    ////printf("[validatation] max = %f\n\n", max);
    return max;
}
//}}}
