// vim: set fdm=marker:
//***************************************************************************
//
// FileName    : selection.cpp
// Purpose     : Develop for top FCNH with H to two photons analysis
// Description : Applying event selection & Preparing histograms for individual dataset.
// Details     : PU reweighting, leptonic/hadronic selections, hists, top reconstruction.
// Author      : Yu-Wei Kao [ykao@cern.ch]
//
//***************************************************************************
//### include & bool (channel, bjet selection){{{
#include <typeinfo>
#include <iostream>
#include <stdio.h>
#include <TFile.h>
#include <TH1D.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include "../include/main.h"
#include "../include/selection.h"
#include "../include/enumhist.h"
#include "../../TopKinFit/kinfit.h"
#include "../../TopKinFit/TopLep.h"
using namespace std;

bool bool_isHadronic;
bool bool_isLeptonic;

//--- control bjet selection ---//
bool bool_bjet_is_loose  = false;
bool bool_bjet_is_medium = false;
bool bool_bjet_is_tight  = true;
bool bool_num_bjets_is_exactly_one = true;
bool bool_num_bjets_is_atleast_one = !bool_num_bjets_is_exactly_one;
//}}}

void Selection(char* input_file, char* output_file, char* output_tree, char* dataset, char* output_dir, char* channel){
    //### bool MC/Data & hadronic/leptonic{{{
    bool isData = isThisDataOrNot(dataset);
    bool isMCsignal = isThisMCsignal(dataset);
    if((string)channel == "hadronic") bool_isHadronic = true; else bool_isHadronic = false;
    bool_isLeptonic = !bool_isHadronic;
    if(bool_isHadronic) printf("[CHECK] isHadronic!\n");
    if(bool_isLeptonic) printf("[CHECK] isLeptonic!\n");
    //}}}
    //### I/O, histograms{{{
    //============================//
    //----- Input file names -----//
    //============================//
    TFile *fin  = TFile::Open(input_file);
    TFile *f_mcpu = TFile::Open("data/MCPileUp.root");
    TFile *fout = new TFile(output_file, "RECREATE");
    TH1D  *h_pu_reweighting_factor = (TH1D*)f_mcpu->Get("puhist");
    //==================================================//
    //----------------- Read input files ---------------//
    //==================================================//
    myTreeClass treeReader;
    treeReader.InitTree("mytree");
    treeReader.AddRootFile(fin);
    treeReader.SetBranchAddresses();
    //==================================================//
    //-------------------- histograms ------------------//
    //==================================================//
    TH1D* h[totalHistNum];
    for(int i=0; i<totalHistNum; i++){
        //printf("[%d] %s, %d, %f, %f\n", i, histNames[i].c_str(), histNbins[i], histBinLow[i], histBinHigh[i]);
        h[i] = new TH1D(histNames[i].c_str(), histNames[i].c_str(), histNbins[i], histBinLow[i], histBinHigh[i]);
        h[i] -> Sumw2();
    }
    //}}}
    //### I/O, tree{{{
    TFile *f_tree = new TFile(output_tree, "RECREATE");
    TTree *myAnalysisTree = new TTree("myAnalysisTree", "myAnalysisTree");
    float tree_event_weight;
    float tree_leadingPhoton_pt;
    float tree_leadingPhoton_eta;
    float tree_leadingPhoton_IDMVA;
    float tree_subleadingPhoton_pt;
    float tree_subleadingPhoton_eta;
    float tree_subleadingPhoton_IDMVA;
    float tree_diphoton_pt;
    float tree_diphoton_eta;
    float tree_diphoton_deltaR;
    float tree_lepton_pt;
    float tree_lepton_eta;
    int   tree_lepton_charge;
    float tree_jet1_pt;
    float tree_jet1_eta;
    float tree_jet1_btag;
    float tree_jet1_CvsL;
    float tree_jet1_CvsB;
    float tree_jet2_pt;
    float tree_jet2_eta;
    float tree_jet2_btag;
    float tree_jet2_CvsL;
    float tree_jet2_CvsB;
    float tree_jet3_pt;
    float tree_jet3_eta;
    float tree_jet3_btag;
    float tree_jet3_CvsL;
    float tree_jet3_CvsB;
    float tree_jet4_pt;
    float tree_jet4_eta;
    float tree_jet4_btag;
    float tree_jet4_CvsL;
    float tree_jet4_CvsB;
    float tree_neutrino_pz;
    float tree_top_mass;
    float tree_top_pt;
    float tree_top_eta;
    float tree_tH_deltaR;
    //-----
    myAnalysisTree->Branch("tree_event_weight", &tree_event_weight, "tree_event_weight/F");
    myAnalysisTree->Branch("tree_leadingPhoton_pt", &tree_leadingPhoton_pt, "tree_leadingPhoton_pt/F");
    myAnalysisTree->Branch("tree_leadingPhoton_eta", &tree_leadingPhoton_eta, "tree_leadingPhoton_eta/F");
    myAnalysisTree->Branch("tree_leadingPhoton_IDMVA", &tree_leadingPhoton_IDMVA, "tree_leadingPhoton_IDMVA/F");
    myAnalysisTree->Branch("tree_subleadingPhoton_pt", &tree_subleadingPhoton_pt, "tree_subleadingPhoton_pt/F");
    myAnalysisTree->Branch("tree_subleadingPhoton_eta", &tree_subleadingPhoton_eta, "tree_subleadingPhoton_eta/F");
    myAnalysisTree->Branch("tree_subleadingPhoton_IDMVA", &tree_subleadingPhoton_IDMVA, "tree_subleadingPhoton_IDMVA/F");
    myAnalysisTree->Branch("tree_diphoton_pt", &tree_diphoton_pt, "tree_diphoton_pt/F");
    myAnalysisTree->Branch("tree_diphoton_eta", &tree_diphoton_eta, "tree_diphoton_eta/F");
    myAnalysisTree->Branch("tree_diphoton_deltaR", &tree_diphoton_deltaR, "tree_diphoton_deltaR/F");
    myAnalysisTree->Branch("tree_lepton_pt", &tree_lepton_pt, "tree_lepton_pt/F");
    myAnalysisTree->Branch("tree_lepton_eta", &tree_lepton_eta, "tree_lepton_eta/F");
    myAnalysisTree->Branch("tree_lepton_charge", &tree_lepton_charge, "tree_lepton_charge/I");
    myAnalysisTree->Branch("tree_jet1_pt", &tree_jet1_pt, "tree_jet1_pt/F");
    myAnalysisTree->Branch("tree_jet1_eta", &tree_jet1_eta, "tree_jet1_eta/F");
    myAnalysisTree->Branch("tree_jet1_btag", &tree_jet1_btag, "tree_jet1_btag/F");
    myAnalysisTree->Branch("tree_jet1_CvsL", &tree_jet1_CvsL, "tree_jet1_CvsL/F");
    myAnalysisTree->Branch("tree_jet1_CvsB", &tree_jet1_CvsB, "tree_jet1_CvsB/F");
    myAnalysisTree->Branch("tree_jet2_pt", &tree_jet2_pt, "tree_jet2_pt/F");
    myAnalysisTree->Branch("tree_jet2_eta", &tree_jet2_eta, "tree_jet2_eta/F");
    myAnalysisTree->Branch("tree_jet2_btag", &tree_jet2_btag, "tree_jet2_btag/F");
    myAnalysisTree->Branch("tree_jet2_CvsL", &tree_jet2_CvsL, "tree_jet2_CvsL/F");
    myAnalysisTree->Branch("tree_jet2_CvsB", &tree_jet2_CvsB, "tree_jet2_CvsB/F");
    myAnalysisTree->Branch("tree_jet3_pt", &tree_jet3_pt, "tree_jet3_pt/F");
    myAnalysisTree->Branch("tree_jet3_eta", &tree_jet3_eta, "tree_jet3_eta/F");
    myAnalysisTree->Branch("tree_jet3_btag", &tree_jet3_btag, "tree_jet3_btag/F");
    myAnalysisTree->Branch("tree_jet3_CvsL", &tree_jet3_CvsL, "tree_jet3_CvsL/F");
    myAnalysisTree->Branch("tree_jet3_CvsB", &tree_jet3_CvsB, "tree_jet3_CvsB/F");
    myAnalysisTree->Branch("tree_jet4_pt", &tree_jet4_pt, "tree_jet4_pt/F");
    myAnalysisTree->Branch("tree_jet4_eta", &tree_jet4_eta, "tree_jet4_eta/F");
    myAnalysisTree->Branch("tree_jet4_btag", &tree_jet4_btag, "tree_jet4_btag/F");
    myAnalysisTree->Branch("tree_jet4_CvsL", &tree_jet4_CvsL, "tree_jet4_CvsL/F");
    myAnalysisTree->Branch("tree_jet4_CvsB", &tree_jet4_CvsB, "tree_jet4_CvsB/F");
    myAnalysisTree->Branch("tree_neutrino_pz", &tree_neutrino_pz, "tree_neutrino_pz/F");
    myAnalysisTree->Branch("tree_top_mass", &tree_top_mass, "tree_top_mass/F");
    myAnalysisTree->Branch("tree_top_pt", &tree_top_pt, "tree_top_pt/F");
    myAnalysisTree->Branch("tree_top_eta", &tree_top_eta, "tree_top_eta/F");
    myAnalysisTree->Branch("tree_tH_deltaR", &tree_tH_deltaR, "tree_tH_deltaR/F");
    //}}}
    //### counters{{{
    int counter_M1_exists_hadronic = 0;
    int counter_M1_exists_leptonic = 0;
    int counter_selected_events = 0;
    int counter_coeff_D_isNegative = 0;
    int counter_bjet_is_bquark = 0;
    int counter_irregular_disc = 0;
    //}}}
    //==================================================//
    //-------------------- Event Loop ------------------//
    //==================================================//
    // # topKinFit method{{{
    KINFIT::kfit *kf = new KINFIT::kfit();
    kf->Init(TOPLEP); // Initialize tool for ttbar with FCNC top decay to Higgs(->bb)+u/c hypothesis
    kf->SetNToy(10); // Set number of toys for minimization
    // Define PDFs{{{
    std::string pdfFileName = "../TopKinFit/test/GenAnalysis/TopLep/pdf.root";
    //kf->SetPDF("TopWMass",pdfFileName.c_str(),"TopLepWM_Fit");
    kf->SetPDF("TopWMass",pdfFileName.c_str(),"TopWM_Fit");
    kf->SetPDF("TopMass",pdfFileName.c_str(),"TopLepRecM_Fit");
    //kf->SetPDF("HiggsMass",pdfFileName.c_str(),"HiggsRecM_Fit");
    //kf->SetPDF("TopHadMass",pdfFileName.c_str(),"TopHadRecM_Fit");
    kf->SetPDF("MetPx",pdfFileName.c_str(),"dMetPx_Gaus");
    kf->SetPDF("MetPy",pdfFileName.c_str(),"dMetPy_Gaus");
    kf->SetPDF("BJetPx",pdfFileName.c_str(),"dBJetPx_Fit");
    kf->SetPDF("BJetPy",pdfFileName.c_str(),"dBJetPy_Fit");
    kf->SetPDF("BJetPz",pdfFileName.c_str(),"dBJetPz_Fit");
    kf->SetPDF("BJetE",pdfFileName.c_str(),"dBJetE_Fit");
    kf->SetPDF("NonBJetPx",pdfFileName.c_str(),"dNonBJetPx_Fit");
    kf->SetPDF("NonBJetPy",pdfFileName.c_str(),"dNonBJetPy_Fit");
    kf->SetPDF("NonBJetPz",pdfFileName.c_str(),"dNonBJetPz_Fit");
    kf->SetPDF("NonBJetE",pdfFileName.c_str(),"dNonBJetE_Fit");
    kf->SetPDF("ElecPx",pdfFileName.c_str(),"dElecPx_Fit");
    kf->SetPDF("ElecPy",pdfFileName.c_str(),"dElecPy_Fit");
    kf->SetPDF("ElecPz",pdfFileName.c_str(),"dElecPz_Fit");
    kf->SetPDF("ElecE",pdfFileName.c_str(),"dElecE_Fit");
    kf->SetPDF("MuonPx",pdfFileName.c_str(),"dMuonPx_Fit");
    kf->SetPDF("MuonPy",pdfFileName.c_str(),"dMuonPy_Fit");
    kf->SetPDF("MuonPz",pdfFileName.c_str(),"dMuonPz_Fit");
    kf->SetPDF("MuonE",pdfFileName.c_str(),"dMuonE_Fit");
    //}}}
    //}}}
    int nentries = treeReader.GetEntries();
    for(int ientry=0; ientry<nentries; ientry++){
        TTree* tmp = treeReader.GetTTree(); tmp->GetEntry(ientry);
        //### entry reporter{{{
        if(ientry==0){
            printf("[INFO] total entry before preselection = %d\n", treeReader.EvtInfo_totalEntry_before_preselection);
            printf("[INFO] total entry after  preselection = %d\n", nentries);
        }
        //if((ientry+1)%10000==0)  printf("ientry = %d\r", ientry);
        //if((ientry+1)==nentries) printf("ientry = %d\r", ientry);
        //}}}
        //### InitTree{{{
        tree_event_weight = -999;
        tree_leadingPhoton_pt = -999;
        tree_leadingPhoton_eta = -999;
        tree_leadingPhoton_IDMVA = -999;
        tree_subleadingPhoton_pt = -999;
        tree_subleadingPhoton_eta = -999;
        tree_subleadingPhoton_IDMVA = -999;
        tree_diphoton_pt = -999;
        tree_diphoton_eta = -999;
        tree_diphoton_deltaR = -999;
        tree_lepton_pt = -999;
        tree_lepton_eta = -999;
        tree_lepton_charge = -999;
        tree_jet1_pt = -999;
        tree_jet1_eta = -999;
        tree_jet1_btag = -999;
        tree_jet1_CvsL = -999;
        tree_jet1_CvsB = -999;
        tree_jet2_pt = -999;
        tree_jet2_eta = -999;
        tree_jet2_btag = -999;
        tree_jet2_CvsL = -999;
        tree_jet2_CvsB = -999;
        tree_jet3_pt = -999;
        tree_jet3_eta = -999;
        tree_jet3_btag = -999;
        tree_jet3_CvsL = -999;
        tree_jet3_CvsB = -999;
        tree_jet4_pt = -999;
        tree_jet4_eta = -999;
        tree_jet4_btag = -999;
        tree_jet4_CvsL = -999;
        tree_jet4_CvsB = -999;
        tree_neutrino_pz = -999;
        tree_top_mass = -999;
        tree_top_pt = -999;
        tree_top_eta = -999;
        tree_tH_deltaR = -999;
        //}}}
        //### PU reweighting{{{
        //========= PU Reweighting =========//
        bool apply_PU_reweighting = true;
        double PU_reweighting_factor = apply_PU_reweighting ? h_pu_reweighting_factor->GetBinContent((int)treeReader.EvtInfo_NPu+1) : 1.;
        double NormalizationFactor = treeReader.EvtInfo_genweight * treeReader.EvtInfo_NormalizationFactor_lumi * PU_reweighting_factor;
        double NormalizationFactor_wopu = treeReader.EvtInfo_genweight * treeReader.EvtInfo_NormalizationFactor_lumi;
        tree_event_weight = NormalizationFactor;
        //Reminder: EvtInfo_NormalizationFactor_lumi = 1000. * Luminosity * CrossSection * BranchingFraction / TotalGenweight;
        //}}}
        //=== pu study only (skipped){{{
        /*
        //--- pu study only ---//
        bool pass_leadingPhotonPT = treeReader.DiPhoInfo_leadPt > 35.;
        bool pass_subleadingPhotonPT = treeReader.DiPhoInfo_subleadPt > 25.;
        bool pass_photon_criteria_pt = pass_leadingPhotonPT && pass_subleadingPhotonPT;

        bool pass_leadingPhotonEta =  (treeReader.DiPhoInfo_leadEta < 1.4442) || (treeReader.DiPhoInfo_leadEta > 1.566 && treeReader.DiPhoInfo_leadEta < 2.5);
        bool pass_subleadingPhotonEta = (treeReader.DiPhoInfo_subleadEta < 1.4442) || (treeReader.DiPhoInfo_subleadEta > 1.566 && treeReader.DiPhoInfo_subleadEta < 2.5);
        bool pass_photon_criteria_eta = pass_leadingPhotonEta && pass_subleadingPhotonEta;

        bool pass_interested_region = treeReader.DiPhoInfo_mass > 100 && treeReader.DiPhoInfo_mass < 180;
        bool pass_signal_region = treeReader.DiPhoInfo_mass>120 && treeReader.DiPhoInfo_mass<130;
        
        //require the quality of photons.
        if(!pass_photon_criteria_pt) continue;
        if(!pass_photon_criteria_eta) continue;
        ////control region
        if(!pass_interested_region) continue;
        if(!isMCsignal && pass_signal_region) continue;
        
        h[hist_EvtInfo_Rho] -> Fill(treeReader.EvtInfo_Rho, isData ? 1. : NormalizationFactor);
        h[hist_EvtInfo_Rho_wopu] -> Fill(treeReader.EvtInfo_Rho, isData ? 1. : NormalizationFactor_wopu);
        h[hist_EvtInfo_NVtx] -> Fill(treeReader.EvtInfo_NVtx, isData ? 1. : NormalizationFactor);
        h[hist_EvtInfo_NVtx_wopu] -> Fill(treeReader.EvtInfo_NVtx, isData ? 1. : NormalizationFactor_wopu);
        */
        //}}}
        //### Event selection{{{
        //blind data
        bool pass_signal_region = treeReader.DiPhoInfo_mass>120 && treeReader.DiPhoInfo_mass<130;
        if(isData && pass_signal_region) continue;
        //supress QCD
        bool pass_photon_IDMVA = treeReader.DiPhoInfo_leadIDMVA>0 && treeReader.DiPhoInfo_subleadIDMVA>0;
        if(!pass_photon_IDMVA) continue;

        bool passEvent=true;
        //--------- Hadronic Channel ---------//
        if(bool_isHadronic){
            if(treeReader.num_leptons>0) passEvent = false;
            //if(treeReader.num_jets<3) passEvent = false;
            if(treeReader.num_jets<4) passEvent = false;
        //--------- Leptonic Channel ---------//
        } else{ // bool_isLeptonic == true
            if(treeReader.num_leptons<1) passEvent = false;
            if(treeReader.num_jets<1) passEvent = false;
        }
        if(!passEvent) continue;
        //--------- check bjets ---------//
        TLorentzVector bjet, jet;
        double bjet_btag_score, bjet_CvsL_score, bjet_CvsB_score;
        int num_bjets = 0, index_bjet = -999;

        std::vector<TLorentzVector> bjets_tight, bjets_loose, bjets_medium;
        std::vector<double> btag_score_tight, btag_score_loose, btag_score_medium;
        std::vector<double> CvsL_score_tight, CvsL_score_loose, CvsL_score_medium;
        std::vector<double> CvsB_score_tight, CvsB_score_loose, CvsB_score_medium;
        int num_bjets_tight = 0, num_bjets_loose = 0, num_bjets_medium = 0;
        int index_leading_bjet_tight = -999, index_leading_bjet_loose = -999, index_leading_bjet_medium = -999;

        //categorize bjet according to deepCSV score{{{
        if(!(treeReader.num_jets<1)){
            for(int i=0; i<treeReader.num_jets; ++i){
                jet.SetPtEtaPhiE(treeReader.JetInfo_jet_pt_selection->at(i),treeReader.JetInfo_jet_eta_selection->at(i),treeReader.JetInfo_jet_phi_selection->at(i),treeReader.JetInfo_jet_energy_selection->at(i));
                double btag_score = treeReader.JetInfo_jet_pfDeepCSVJetTags_probb_selection->at(i)+treeReader.JetInfo_jet_pfDeepCSVJetTags_probbb_selection->at(i);
                //---
                double ctag_probc = treeReader.JetInfo_jet_pfDeepCSVJetTags_probc_selection->at(i);
                double CvsL_score = ctag_probc / ( ctag_probc + treeReader.JetInfo_jet_pfDeepCSVJetTags_probudsg_selection->at(i) );
                double CvsB_score = ctag_probc / ( ctag_probc + btag_score );

                if(btag_score >= pfDeepCSVJetTags_loose){
                    num_bjets_loose += 1;
                    bjets_loose.push_back(jet);
                    btag_score_loose.push_back(btag_score);
                    CvsL_score_loose.push_back(CvsL_score);
                    CvsB_score_loose.push_back(CvsB_score);
                    if(num_bjets_loose == 1) index_leading_bjet_loose = i;
                }
                if(btag_score >= pfDeepCSVJetTags_medium){
                    num_bjets_medium += 1;
                    bjets_medium.push_back(jet);
                    btag_score_medium.push_back(btag_score);
                    CvsL_score_medium.push_back(CvsL_score);
                    CvsB_score_medium.push_back(CvsB_score);
                    if(num_bjets_medium == 1) index_leading_bjet_medium = i;
                }
                if(btag_score >= pfDeepCSVJetTags_tight){
                    num_bjets_tight += 1;
                    bjets_tight.push_back(jet);
                    btag_score_tight.push_back(btag_score);
                    CvsL_score_tight.push_back(CvsL_score);
                    CvsB_score_tight.push_back(CvsB_score);
                    if(num_bjets_tight == 1) index_leading_bjet_tight = i;
                }
            }//end of looping jets
        }
        //}}}
        //determine leading bjet according to chosen WP{{{
        if(bool_bjet_is_loose){
            index_bjet = index_leading_bjet_loose;
            if(index_bjet != -999){
                bjet = bjets_loose[0];
                bjet_btag_score = btag_score_loose[0];
                bjet_CvsL_score = CvsL_score_loose[0];
                bjet_CvsB_score = CvsB_score_loose[0];
                num_bjets = num_bjets_loose;
            }
        }
        if(bool_bjet_is_medium){
            index_bjet = index_leading_bjet_medium;
            if(index_bjet != -999){
                bjet = bjets_medium[0];
                bjet_btag_score = btag_score_medium[0];
                bjet_CvsL_score = CvsL_score_medium[0];
                bjet_CvsB_score = CvsB_score_medium[0];
                num_bjets = num_bjets_medium;
            }
        }
        if(bool_bjet_is_tight){
            index_bjet = index_leading_bjet_tight;
            if(index_bjet != -999){
                bjet = bjets_tight[0];
                bjet_btag_score = btag_score_tight[0];
                bjet_CvsL_score = CvsL_score_tight[0];
                bjet_CvsB_score = CvsB_score_tight[0];
                num_bjets = num_bjets_tight;
            }
        }
        //}}}

        bool pass_bjets_multiplicity_selection;
        if(bool_num_bjets_is_exactly_one) pass_bjets_multiplicity_selection = num_bjets == 1;
        if(bool_num_bjets_is_atleast_one) pass_bjets_multiplicity_selection = num_bjets >= 1;
        if(!pass_bjets_multiplicity_selection) continue;

        //store leading bjet info{{{
        if(index_bjet != -999){//at least one bjet!
            h[hist_leading_bjet_pt] -> Fill(bjet.Pt(), isData ? 1. : NormalizationFactor);
            h[hist_leading_bjet_eta] -> Fill(bjet.Eta(), isData ? 1. : NormalizationFactor);
            h[hist_leading_bjet_phi] -> Fill(bjet.Phi(), isData ? 1. : NormalizationFactor);
            h[hist_leading_bjet_energy] -> Fill(bjet.E(), isData ? 1. : NormalizationFactor);
            h[hist_leading_bjet_btag_score] -> Fill(bjet_btag_score, isData ? 1. : NormalizationFactor);
            h[hist_leading_bjet_CvsL_score] -> Fill(bjet_CvsL_score, isData ? 1. : NormalizationFactor);
            h[hist_leading_bjet_CvsB_score] -> Fill(bjet_CvsB_score, isData ? 1. : NormalizationFactor);

            tree_jet1_pt = bjet.Pt();
            tree_jet1_eta = bjet.Eta();
            tree_jet1_btag = bjet_btag_score;
            tree_jet1_CvsL = bjet_CvsL_score;
            tree_jet1_CvsB = bjet_CvsB_score;
        }
        //}}}

        //}}}
        //### Store Info{{{
        //========= Store Info =========//
        h[hist_EvtInfo_NPu] -> Fill(treeReader.EvtInfo_NPu, isData ? 1. : NormalizationFactor);
        h[hist_EvtInfo_Rho] -> Fill(treeReader.EvtInfo_Rho, isData ? 1. : NormalizationFactor);
        h[hist_EvtInfo_NVtx] -> Fill(treeReader.EvtInfo_NVtx, isData ? 1. : NormalizationFactor);
        h[hist_EvtInfo_NVtx_wopu] -> Fill(treeReader.EvtInfo_NVtx, isData ? 1. : NormalizationFactor_wopu);
        h[hist_EvtInfo_genweight] -> Fill(treeReader.EvtInfo_genweight, isData ? 1. : NormalizationFactor);
        //}}}
        //### Diphoton{{{
        //--------- Diphoton ---------//
        TLorentzVector leading_photon, subleading_photon, diphoton;
        h[hist_DiPhoInfo_mass] -> Fill(treeReader.DiPhoInfo_mass, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_pt] -> Fill(treeReader.DiPhoInfo_pt, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_eta] -> Fill(treeReader.DiPhoInfo_eta, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_phi] -> Fill(treeReader.DiPhoInfo_phi, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_energy] -> Fill(treeReader.DiPhoInfo_energy, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_leadPt] -> Fill(treeReader.DiPhoInfo_leadPt, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_leadEta] -> Fill(treeReader.DiPhoInfo_leadEta, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_leadPhi] -> Fill(treeReader.DiPhoInfo_leadPhi, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_leadE] -> Fill(treeReader.DiPhoInfo_leadE, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_leadhoe] -> Fill(treeReader.DiPhoInfo_leadhoe, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_leadIDMVA] -> Fill(treeReader.DiPhoInfo_leadIDMVA, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_subleadPt] -> Fill(treeReader.DiPhoInfo_subleadPt, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_subleadEta] -> Fill(treeReader.DiPhoInfo_subleadEta, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_subleadPhi] -> Fill(treeReader.DiPhoInfo_subleadPhi, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_subleadE] -> Fill(treeReader.DiPhoInfo_subleadE, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_subleadhoe] -> Fill(treeReader.DiPhoInfo_subleadhoe, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_subleadIDMVA] -> Fill(treeReader.DiPhoInfo_subleadIDMVA, isData ? 1. : NormalizationFactor);
        leading_photon.SetPtEtaPhiE(treeReader.DiPhoInfo_leadPt, treeReader.DiPhoInfo_leadEta, treeReader.DiPhoInfo_leadPhi, treeReader.DiPhoInfo_leadE);
        subleading_photon.SetPtEtaPhiE(treeReader.DiPhoInfo_subleadPt, treeReader.DiPhoInfo_subleadEta, treeReader.DiPhoInfo_subleadPhi, treeReader.DiPhoInfo_subleadE);
        diphoton.SetPtEtaPhiE(treeReader.DiPhoInfo_pt, treeReader.DiPhoInfo_eta, treeReader.DiPhoInfo_phi, treeReader.DiPhoInfo_energy);

        h[hist_deltaR_photon_photon] -> Fill(leading_photon.DeltaR(subleading_photon), isData ? 1. : NormalizationFactor) ;

        tree_leadingPhoton_pt = treeReader.DiPhoInfo_leadPt;
        tree_leadingPhoton_eta = treeReader.DiPhoInfo_leadEta;
        tree_leadingPhoton_IDMVA = treeReader.DiPhoInfo_leadIDMVA;
        tree_subleadingPhoton_pt = treeReader.DiPhoInfo_subleadPt;
        tree_subleadingPhoton_eta = treeReader.DiPhoInfo_subleadEta;
        tree_subleadingPhoton_IDMVA = treeReader.DiPhoInfo_subleadIDMVA;
        tree_diphoton_pt = treeReader.DiPhoInfo_pt;
        tree_diphoton_eta = treeReader.DiPhoInfo_eta;
        tree_diphoton_deltaR = leading_photon.DeltaR(subleading_photon);

        //}}}
        //### Leptons{{{
        //--------- Leptons ---------//
        std::vector<TLorentzVector> Leptons;
        std::vector<TLorentzVector> Electrons;
        std::vector<TLorentzVector> Muons;
        std::vector<int> Leptons_charge;
        h[hist_ElecInfo_Size] -> Fill(treeReader.ElecInfo_Size, isData ? 1. : NormalizationFactor);
        h[hist_MuonInfo_Size] -> Fill(treeReader.MuonInfo_Size, isData ? 1. : NormalizationFactor);
        h[hist_num_leptons] -> Fill(treeReader.num_leptons, isData ? 1. : NormalizationFactor);// # of selected objects.
        h[hist_num_electrons] -> Fill(treeReader.num_electrons, isData ? 1. : NormalizationFactor);// # of selected objects.
        h[hist_num_muons] -> Fill(treeReader.num_muons, isData ? 1. : NormalizationFactor);// # of selected objects.
        if(treeReader.num_electrons>0){
            for(int i=0; i<treeReader.num_electrons; ++i){
                h[hist_ElecInfo_electron_charge] -> Fill(treeReader.ElecInfo_electron_charge_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_ElecInfo_electron_pt] -> Fill(treeReader.ElecInfo_electron_pt_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_ElecInfo_electron_eta] -> Fill(treeReader.ElecInfo_electron_eta_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_ElecInfo_electron_phi] -> Fill(treeReader.ElecInfo_electron_phi_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_ElecInfo_electron_energy] -> Fill(treeReader.ElecInfo_electron_energy_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_ElecInfo_electron_diphoton_deltaR] -> Fill(treeReader.ElecInfo_electron_diphoton_deltaR_selection->at(i), isData ? 1. : NormalizationFactor);
                //------------------------
                h[hist_lepton_charge] -> Fill(treeReader.ElecInfo_electron_charge_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_pt] -> Fill(treeReader.ElecInfo_electron_pt_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_eta] -> Fill(treeReader.ElecInfo_electron_eta_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_phi] -> Fill(treeReader.ElecInfo_electron_phi_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_energy] -> Fill(treeReader.ElecInfo_electron_energy_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_diphoton_deltaR] -> Fill(treeReader.ElecInfo_electron_diphoton_deltaR_selection->at(i), isData ? 1. : NormalizationFactor);
                //------------------------
                TLorentzVector electron; electron.SetPtEtaPhiE(treeReader.ElecInfo_electron_pt_selection->at(i), treeReader.ElecInfo_electron_eta_selection->at(i), treeReader.ElecInfo_electron_phi_selection->at(i), treeReader.ElecInfo_electron_energy_selection->at(i));
                Electrons.push_back(electron);
                Leptons.push_back(electron);
                Leptons_charge.push_back(treeReader.ElecInfo_electron_charge_selection->at(i));

                //printf("[Check-ele] ");
                //printf("pt = %6.2f, ", treeReader.ElecInfo_electron_pt_selection->at(i));
                //printf("eta = %6.2f, ", treeReader.ElecInfo_electron_eta_selection->at(i));
                //printf("phi = %6.2f, ", treeReader.ElecInfo_electron_phi_selection->at(i));
                //printf("energy = %6.2f\n", treeReader.ElecInfo_electron_energy_selection->at(i));

                //printf("[Check-lep] ");
                //printf("pt = %6.2f, ", Leptons.at(i).Pt());
                //printf("eta = %6.2f, ", Leptons.at(i).Eta());
                //printf("phi = %6.2f, ", Leptons.at(i).Phi());
                //printf("energy = %6.2f\n", Leptons.at(i).Energy());


                //printf("[Check-ele] ");
                //printf("px = %6.2f, ", electron.Px());
                //printf("py = %6.2f, ", electron.Py());
                //printf("pz = %6.2f, ", electron.Pz());
                //printf("energy = %6.2f\n", electron.Energy());

                //printf("[Check-lep] ");
                //printf("px = %6.2f, ", Leptons.at(i).Px());
                //printf("py = %6.2f, ", Leptons.at(i).Py());
                //printf("pz = %6.2f, ", Leptons.at(i).Pz());
                //printf("energy = %6.2f\n", Leptons.at(i).Energy());

                //printf("[Check-self] ");
                //printf("px = %6.2f, ", electron.Pt() * TMath::Cos(electron.Phi()));
                //printf("py = %6.2f, ", electron.Pt() * TMath::Sin(electron.Phi()));
                //printf("pz = %6.2f, ", electron.Pt() * TMath::SinH(electron.Eta()));
                //float momentum = electron.Pt() * TMath::CosH(electron.Eta());
                //printf("p  = %6.2f, ", momentum);
                //float mass = electron.M();
                //printf("m  = %6.2f, ", mass );
                //printf("energy = %6.2f\n", sqrt(momentum*momentum + mass*mass) );

            }
        }
        if(treeReader.num_muons>0){
            for(int i=0; i<treeReader.num_muons; ++i){
                h[hist_MuonInfo_muon_charge] -> Fill(treeReader.MuonInfo_muon_charge_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_MuonInfo_muon_pt] -> Fill(treeReader.MuonInfo_muon_pt_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_MuonInfo_muon_eta] -> Fill(treeReader.MuonInfo_muon_eta_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_MuonInfo_muon_phi] -> Fill(treeReader.MuonInfo_muon_phi_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_MuonInfo_muon_energy] -> Fill(treeReader.MuonInfo_muon_energy_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_MuonInfo_muon_diphoton_deltaR] -> Fill(treeReader.MuonInfo_muon_diphoton_deltaR_selection->at(i), isData ? 1. : NormalizationFactor);
                //------------------------
                h[hist_lepton_charge] -> Fill(treeReader.MuonInfo_muon_charge_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_pt] -> Fill(treeReader.MuonInfo_muon_pt_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_eta] -> Fill(treeReader.MuonInfo_muon_eta_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_phi] -> Fill(treeReader.MuonInfo_muon_phi_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_energy] -> Fill(treeReader.MuonInfo_muon_energy_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_diphoton_deltaR] -> Fill(treeReader.MuonInfo_muon_diphoton_deltaR_selection->at(i), isData ? 1. : NormalizationFactor);
                //------------------------
                TLorentzVector muon; muon.SetPtEtaPhiE(treeReader.MuonInfo_muon_pt_selection->at(i), treeReader.MuonInfo_muon_eta_selection->at(i), treeReader.MuonInfo_muon_phi_selection->at(i), treeReader.MuonInfo_muon_energy_selection->at(i));
                Muons.push_back(muon);
                Leptons.push_back(muon);
                Leptons_charge.push_back(treeReader.MuonInfo_muon_charge_selection->at(i));
            }
        }

        // less than 0.3% events have more than one charged lepton.
        if(bool_isLeptonic){
            tree_lepton_pt = Leptons[0].Pt();
            tree_lepton_eta = Leptons[0].Eta();
            tree_lepton_charge = Leptons_charge[0];
        }
        //}}}
        //### Jets{{{
        //--------- Jets ---------//
        std::vector<TLorentzVector> Jets;
        std::vector<double> Jets_btag_score, Jets_CvsL_score, Jets_CvsB_score;
        h[hist_jets_size] -> Fill(treeReader.jets_size, isData ? 1. : NormalizationFactor);
        h[hist_num_jets] -> Fill(treeReader.num_jets, isData ? 1. : NormalizationFactor);
        if(treeReader.num_jets>0){
            for(int i=0; i<treeReader.num_jets; ++i){
                h[hist_JetInfo_jet_pt] -> Fill(treeReader.JetInfo_jet_pt_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_JetInfo_jet_eta] -> Fill(treeReader.JetInfo_jet_eta_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_JetInfo_jet_phi] -> Fill(treeReader.JetInfo_jet_phi_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_JetInfo_jet_energy] -> Fill(treeReader.JetInfo_jet_energy_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_JetInfo_jet_diphoton_deltaR] -> Fill(treeReader.JetInfo_jet_diphoton_deltaR_selection->at(i), isData ? 1. : NormalizationFactor);
                TLorentzVector jet; jet.SetPtEtaPhiE(treeReader.JetInfo_jet_pt_selection->at(i),treeReader.JetInfo_jet_eta_selection->at(i),treeReader.JetInfo_jet_phi_selection->at(i),treeReader.JetInfo_jet_energy_selection->at(i));
                Jets.push_back(jet);
                double ctag_probc = treeReader.JetInfo_jet_pfDeepCSVJetTags_probc_selection->at(i);
                double btag_score = treeReader.JetInfo_jet_pfDeepCSVJetTags_probb_selection->at(i)+treeReader.JetInfo_jet_pfDeepCSVJetTags_probbb_selection->at(i);
                double CvsL_score = ctag_probc / ( ctag_probc + treeReader.JetInfo_jet_pfDeepCSVJetTags_probudsg_selection->at(i) );
                double CvsB_score = ctag_probc / ( ctag_probc + btag_score );
                Jets_btag_score.push_back(btag_score);
                Jets_CvsL_score.push_back(CvsL_score);
                Jets_CvsB_score.push_back(CvsB_score);
                if(i==0){//leading jet
                    h[hist_jet1_pt] -> Fill(treeReader.JetInfo_jet_pt_selection->at(i), isData ? 1. : NormalizationFactor);
                    h[hist_jet1_eta] -> Fill(treeReader.JetInfo_jet_eta_selection->at(i), isData ? 1. : NormalizationFactor);
                    h[hist_jet1_phi] -> Fill(treeReader.JetInfo_jet_phi_selection->at(i), isData ? 1. : NormalizationFactor);
                    h[hist_jet1_energy] -> Fill(treeReader.JetInfo_jet_energy_selection->at(i), isData ? 1. : NormalizationFactor);
                    h[hist_jet1_btag_score] -> Fill(btag_score, isData ? 1. : NormalizationFactor);
                    h[hist_jet1_CvsL_score] -> Fill(CvsL_score, isData ? 1. : NormalizationFactor);
                    h[hist_jet1_CvsB_score] -> Fill(CvsB_score, isData ? 1. : NormalizationFactor);
                    h[hist_jet1_diphoton_deltaR] -> Fill(treeReader.JetInfo_jet_diphoton_deltaR_selection->at(i), isData ? 1. : NormalizationFactor);
                    if(treeReader.num_leptons>0){
                        for(int i=0; i<treeReader.num_leptons; ++i){
                            double delta_R = jet.DeltaR(Leptons.at(i));
                            h[hist_jet1_lepton_deltaR] -> Fill(delta_R, isData ? 1. : NormalizationFactor);
                        }
                    }
                }
                if(i==1){//subleading jet
                    h[hist_jet2_pt] -> Fill(treeReader.JetInfo_jet_pt_selection->at(i), isData ? 1. : NormalizationFactor);
                    h[hist_jet2_eta] -> Fill(treeReader.JetInfo_jet_eta_selection->at(i), isData ? 1. : NormalizationFactor);
                    h[hist_jet2_phi] -> Fill(treeReader.JetInfo_jet_phi_selection->at(i), isData ? 1. : NormalizationFactor);
                    h[hist_jet2_energy] -> Fill(treeReader.JetInfo_jet_energy_selection->at(i), isData ? 1. : NormalizationFactor);
                    h[hist_jet2_btag_score] -> Fill(btag_score, isData ? 1. : NormalizationFactor);
                    h[hist_jet2_CvsL_score] -> Fill(CvsL_score, isData ? 1. : NormalizationFactor);
                    h[hist_jet2_CvsB_score] -> Fill(CvsB_score, isData ? 1. : NormalizationFactor);
                    h[hist_jet2_diphoton_deltaR] -> Fill(treeReader.JetInfo_jet_diphoton_deltaR_selection->at(i), isData ? 1. : NormalizationFactor);
                    if(treeReader.num_leptons>0){
                        for(int i=0; i<treeReader.num_leptons; ++i){
                            double delta_R = jet.DeltaR(Leptons.at(i));
                            h[hist_jet2_lepton_deltaR] -> Fill(delta_R, isData ? 1. : NormalizationFactor);
                        }
                    }
                }
            }//end of looping jets
        }
        //}}}
        //### MET{{{
        h[hist_MetInfo_Pt] -> Fill(treeReader.MetInfo_Pt, isData ? 1. : NormalizationFactor);
        h[hist_MetInfo_Phi] -> Fill(treeReader.MetInfo_Phi, isData ? 1. : NormalizationFactor);
        h[hist_MetInfo_Px] -> Fill(treeReader.MetInfo_Px, isData ? 1. : NormalizationFactor);
        h[hist_MetInfo_Py] -> Fill(treeReader.MetInfo_Py, isData ? 1. : NormalizationFactor);
        h[hist_MetInfo_SumET] -> Fill(treeReader.MetInfo_SumET, isData ? 1. : NormalizationFactor);
        //}}}
        //### top reconstruction
        //================================================//
        //-----------   top reconstruction     -----------//
        //================================================//
        std::vector<int> index_jet_chi2_modified(2, -999);
        std::vector<TLorentzVector> jet_chi2_modified(2);

        if(bool_isHadronic){
            //### chi-2 sorting{{{
            //----- chi-2 modified -----//
            double chi2_min = 99999;
            for(int i=0; i<treeReader.num_jets; ++i){
                if(i==index_bjet) continue;//bypass bjet
                for(int j=i+1; j<treeReader.num_jets; ++j){
                    if(j==index_bjet) continue;//bypass bjet
                    TLorentzVector w_candidate = Jets[i] + Jets[j];
                    double w_mass = w_candidate.M();
                    TLorentzVector top_candidate = w_candidate + bjet;
                    double t_mass = top_candidate.M();
                    double chi2 = Chi2_calculator_modified(w_mass, t_mass);
                    if(chi2 < chi2_min){
                        index_jet_chi2_modified[0] = i;
                        index_jet_chi2_modified[1] = j;
                        jet_chi2_modified[0] = Jets[i];
                        jet_chi2_modified[1] = Jets[j];
                        chi2_min = chi2;
                    }
                }
            }//end of looping jets
            if(chi2_min!=99999){
                h[hist_wjet1_pt] -> Fill(jet_chi2_modified[0].Pt(), isData ? 1. : NormalizationFactor);
                h[hist_wjet1_eta] -> Fill(jet_chi2_modified[0].Eta(), isData ? 1. : NormalizationFactor);
                h[hist_wjet1_phi] -> Fill(jet_chi2_modified[0].Phi(), isData ? 1. : NormalizationFactor);
                h[hist_wjet1_energy] -> Fill(jet_chi2_modified[0].E(), isData ? 1. : NormalizationFactor);
                h[hist_wjet1_btag_score] -> Fill(Jets_btag_score[index_jet_chi2_modified[0]], isData ? 1. : NormalizationFactor);
                h[hist_wjet1_CvsL_score] -> Fill(Jets_CvsL_score[index_jet_chi2_modified[0]], isData ? 1. : NormalizationFactor);
                h[hist_wjet1_CvsB_score] -> Fill(Jets_CvsB_score[index_jet_chi2_modified[0]], isData ? 1. : NormalizationFactor);
                h[hist_wjet1_diphoton_deltaR] -> Fill(jet_chi2_modified[0].DeltaR(diphoton), isData ? 1. : NormalizationFactor);
                if(treeReader.num_leptons>0){
                    for(int i=0; i<treeReader.num_leptons; ++i){
                        double delta_R = jet_chi2_modified[0].DeltaR(Leptons.at(i));
                        h[hist_wjet1_lepton_deltaR] -> Fill(delta_R, isData ? 1. : NormalizationFactor);
                    }
                }
                //---
                h[hist_wjet2_pt] -> Fill(jet_chi2_modified[1].Pt(), isData ? 1. : NormalizationFactor);
                h[hist_wjet2_eta] -> Fill(jet_chi2_modified[1].Eta(), isData ? 1. : NormalizationFactor);
                h[hist_wjet2_phi] -> Fill(jet_chi2_modified[1].Phi(), isData ? 1. : NormalizationFactor);
                h[hist_wjet2_energy] -> Fill(jet_chi2_modified[1].E(), isData ? 1. : NormalizationFactor);
                h[hist_wjet2_btag_score] -> Fill(Jets_btag_score[index_jet_chi2_modified[0]], isData ? 1. : NormalizationFactor);
                h[hist_wjet2_CvsL_score] -> Fill(Jets_CvsL_score[index_jet_chi2_modified[0]], isData ? 1. : NormalizationFactor);
                h[hist_wjet2_CvsB_score] -> Fill(Jets_CvsB_score[index_jet_chi2_modified[0]], isData ? 1. : NormalizationFactor);
                h[hist_wjet2_diphoton_deltaR] -> Fill(jet_chi2_modified[1].DeltaR(diphoton), isData ? 1. : NormalizationFactor);
                if(treeReader.num_leptons>0){
                    for(int i=0; i<treeReader.num_leptons; ++i){
                        double delta_R = jet_chi2_modified[1].DeltaR(Leptons.at(i));
                        h[hist_wjet2_lepton_deltaR] -> Fill(delta_R, isData ? 1. : NormalizationFactor);
                    }
                }
            
                tree_jet3_pt   = jet_chi2_modified[0].Pt();
                tree_jet3_eta  = jet_chi2_modified[0].Eta();
                tree_jet3_btag = Jets_btag_score[index_jet_chi2_modified[0]];
                tree_jet3_CvsL = Jets_CvsL_score[index_jet_chi2_modified[0]];
                tree_jet3_CvsB = Jets_CvsB_score[index_jet_chi2_modified[0]];

                tree_jet4_pt   = jet_chi2_modified[1].Pt();
                tree_jet4_eta  = jet_chi2_modified[1].Eta();
                tree_jet4_btag = Jets_btag_score[index_jet_chi2_modified[1]];
                tree_jet4_CvsL = Jets_CvsL_score[index_jet_chi2_modified[1]];
                tree_jet4_CvsB = Jets_CvsB_score[index_jet_chi2_modified[1]];
            }
            //}}}
            //### Reconstruct Mass (W, M2, M1){{{
            TLorentzVector w_candidate_chi2_modified = jet_chi2_modified[0] + jet_chi2_modified[1];
            TLorentzVector top_candidate_chi2_modified = w_candidate_chi2_modified + bjet;
            h[hist_hadronic_w_candidate_pt]->Fill(w_candidate_chi2_modified.Pt(), isData ? 1. : NormalizationFactor);
            h[hist_hadronic_w_candidate_eta]->Fill(w_candidate_chi2_modified.Eta(), isData ? 1. : NormalizationFactor);
            h[hist_hadronic_w_candidate_mass]->Fill(w_candidate_chi2_modified.M(), isData ? 1. : NormalizationFactor);
            h[hist_hadronic_top_tbw_pt]->Fill(top_candidate_chi2_modified.Pt(), isData ? 1. : NormalizationFactor);
            h[hist_hadronic_top_tbw_eta]->Fill(top_candidate_chi2_modified.Eta(), isData ? 1. : NormalizationFactor);
            h[hist_hadronic_top_tbw_mass]->Fill(top_candidate_chi2_modified.M(), isData ? 1. : NormalizationFactor);
            double M1 = -999;
            int index_q;
            TLorentzVector jet_q;
            TLorentzVector top_fcnh = GetBestM1(M1, treeReader.num_jets, index_bjet, index_jet_chi2_modified, diphoton, Jets, index_q, jet_q);
            if(M1 != -999){
                counter_M1_exists_hadronic += 1;
                h[hist_top_tqh_pt]->Fill(top_fcnh.Pt(), isData ? 1. : NormalizationFactor);//exclude event without suitable candidate
                h[hist_top_tqh_eta]->Fill(top_fcnh.Eta(), isData ? 1. : NormalizationFactor);//exclude event without suitable candidate
                h[hist_top_tqh_mass]->Fill(top_fcnh.M(), isData ? 1. : NormalizationFactor);//exclude event without suitable candidate
                //-----
                h[hist_jetq_pt] -> Fill(jet_q.Pt(), isData ? 1. : NormalizationFactor);
                h[hist_jetq_eta] -> Fill(jet_q.Eta(), isData ? 1. : NormalizationFactor);
                h[hist_jetq_phi] -> Fill(jet_q.Phi(), isData ? 1. : NormalizationFactor);
                h[hist_jetq_energy] -> Fill(jet_q.E(), isData ? 1. : NormalizationFactor);
                h[hist_jetq_btag_score] -> Fill(Jets_btag_score[index_q], isData ? 1. : NormalizationFactor);
                h[hist_jetq_CvsL_score] -> Fill(Jets_CvsL_score[index_q], isData ? 1. : NormalizationFactor);
                h[hist_jetq_CvsB_score] -> Fill(Jets_CvsB_score[index_q], isData ? 1. : NormalizationFactor);
                h[hist_jetq_diphoton_deltaR] -> Fill(jet_q.DeltaR(diphoton), isData ? 1. : NormalizationFactor);
                if(treeReader.num_leptons>0){
                    for(int i=0; i<treeReader.num_leptons; ++i){
                        double delta_R = jet_q.DeltaR(Leptons.at(i));
                        h[hist_jetq_lepton_deltaR] -> Fill(delta_R, isData ? 1. : NormalizationFactor);
                    }
                }

                tree_jet2_pt   = jet_q.Pt();
                tree_jet2_eta  = jet_q.Eta();
                tree_jet2_btag = Jets_btag_score[index_q];
                tree_jet2_CvsL = Jets_CvsL_score[index_q];
                tree_jet2_CvsB = Jets_CvsB_score[index_q];
            }
            //}}}
            //### deltaR{{{
            double deltaR;
            if(M1 != -999){ deltaR = top_candidate_chi2_modified.DeltaR(top_fcnh) ; h[hist_deltaR_top_top] -> Fill(deltaR, isData ? 1. : NormalizationFactor)       ; }
            if(M1 != -999){ deltaR = diphoton.DeltaR(jet_q)                       ; h[hist_deltaR_qH] -> Fill(deltaR, isData ? 1. : NormalizationFactor)            ; }
            deltaR = Jets[0].DeltaR(Jets[1])                                      ; h[hist_deltaR_jet1_jet2] -> Fill(deltaR, isData ? 1. : NormalizationFactor)     ;
            deltaR = jet_chi2_modified[0].DeltaR(jet_chi2_modified[1])            ; h[hist_deltaR_wjet1_wjet2] -> Fill(deltaR, isData ? 1. : NormalizationFactor)   ;
            deltaR = bjet.DeltaR(w_candidate_chi2_modified)                       ; h[hist_deltaR_bW] -> Fill(deltaR, isData ? 1. : NormalizationFactor)            ;
            deltaR = diphoton.DeltaR(w_candidate_chi2_modified)                   ; h[hist_deltaR_HW] -> Fill(deltaR, isData ? 1. : NormalizationFactor)            ;
            deltaR = diphoton.DeltaR(top_candidate_chi2_modified)                 ; h[hist_deltaR_tH] -> Fill(deltaR, isData ? 1. : NormalizationFactor)            ;
            //}}}

            tree_top_pt = top_candidate_chi2_modified.Pt();
            tree_top_eta = top_candidate_chi2_modified.Eta();
            tree_top_mass = top_candidate_chi2_modified.M();
            tree_tH_deltaR = diphoton.DeltaR(top_candidate_chi2_modified);

        } else{ // bool_isLeptonic == true
        //### set known info of met & leading lepton{{{
        float met_pt = treeReader.MetInfo_Pt;
        float met_phi = treeReader.MetInfo_Phi;
        float met_px = treeReader.MetInfo_Px;
        float met_py = treeReader.MetInfo_Py;
        float met_sumET = treeReader.MetInfo_SumET;
        TLorentzVector lepton = Leptons[0]; // leading lepton
        float lepton_px = lepton.Px();
        float lepton_py = lepton.Py();
        float lepton_pz = lepton.Pz();
        float lepton_energy = lepton.E();
        //}}}
        /*
        //--- solve met_pz{{{
        //float coefficient_factor = ( w_boson_mass*w_boson_mass + 2*lepton_px*met_px + 2*lepton_py*met_py ) / (2.*lepton_energy);
        float coefficient_factor = ( 80.375*80.375 + 2.*lepton_px*met_px + 2.*lepton_py*met_py ) / (2.*lepton_energy);
        //float coefficient_A = 1. - (lepton_pz/lepton_energy)*(lepton_pz/lepton_energy);
        float coefficient_A = 1. - (lepton_pz*lepton_pz)/(lepton_energy*lepton_energy);
        float coefficient_B = 2.*coefficient_factor*lepton_pz/lepton_energy;
        float coefficient_C = met_pt*met_pt - coefficient_factor*coefficient_factor;
        float coefficient_D = coefficient_B*coefficient_B - 4.*coefficient_A*coefficient_C;
        
        //printf("[CHECK-coeff] coefficient_A = %6.2f, ", coefficient_A);
        //printf("lepton_pz = %6.2f, ", lepton_pz);
        //printf("lepton_energy = %6.2f \n", lepton_energy);

        h[hist_MetInfo_coeff_A] -> Fill(coefficient_A, isData ? 1. : NormalizationFactor);
        h[hist_MetInfo_coeff_B2A] -> Fill(-coefficient_B / (2*coefficient_A), isData ? 1. : NormalizationFactor);

        float met_pz_solution_1 = 0.0;
        float met_pz_solution_2 = 0.0;
        if(coefficient_D < 0){
            counter_coeff_D_isNegative += 1;
            //printf("[check] coefficient_D = %f\n", coefficient_D);
            met_pz_solution_1 = coefficient_B / (2.*coefficient_A);
            met_pz_solution_2 = coefficient_B / (2.*coefficient_A);
            //met_pz_solution_1 = sqrt( coefficient_C / coefficient_A);
            //met_pz_solution_2 = sqrt( coefficient_C / coefficient_A);
            h[hist_MetInfo_coeff_D] -> Fill(-sqrt(-coefficient_D), isData ? 1. : NormalizationFactor);//keep tracking negative value
            h[hist_MetInfo_coeff_D2A] -> Fill(-sqrt(-coefficient_D) / (2*coefficient_A), isData ? 1. : NormalizationFactor);
        } else{
            met_pz_solution_1 = (coefficient_B + TMath::Sqrt(coefficient_D))/(2.*coefficient_A);
            met_pz_solution_2 = (coefficient_B - TMath::Sqrt(coefficient_D))/(2.*coefficient_A);
            h[hist_MetInfo_coeff_D] -> Fill(sqrt(coefficient_D), isData ? 1. : NormalizationFactor);
            h[hist_MetInfo_coeff_D2A] -> Fill(sqrt(coefficient_D) / (2*coefficient_A), isData ? 1. : NormalizationFactor);
        }
        float larger_pz  = (abs(met_pz_solution_1) > abs(met_pz_solution_2) ) ? met_pz_solution_1 : met_pz_solution_2;
        float smaller_pz = (abs(met_pz_solution_1) < abs(met_pz_solution_2) ) ? met_pz_solution_1 : met_pz_solution_2;
        met_pz_solution_1 = larger_pz;
        met_pz_solution_2 = smaller_pz;

        h[hist_MetInfo_Pz_solution_1] -> Fill(met_pz_solution_1, isData ? 1. : NormalizationFactor);
        h[hist_MetInfo_Pz_solution_2] -> Fill(met_pz_solution_2, isData ? 1. : NormalizationFactor);
        //}}}
        //--- W, M2{{{
        TLorentzVector L_b_lep = bjet;
        TLorentzVector L_met_lep[2];
        TLorentzVector L_w_lep[2];
        TLorentzVector L_bw_lep[2];

        float met_energy_solution_1 = TMath::Sqrt(met_pt*met_pt + met_pz_solution_1*met_pz_solution_1);
        float met_energy_solution_2 = TMath::Sqrt(met_pt*met_pt + met_pz_solution_2*met_pz_solution_2);
        L_met_lep[0].SetPxPyPzE( met_px, met_py, met_pz_solution_1, met_energy_solution_1 );
        L_met_lep[1].SetPxPyPzE( met_px, met_py, met_pz_solution_2, met_energy_solution_2 );

        L_w_lep[0].SetPxPyPzE(  (lepton_px + met_px), (lepton_py + met_py), (lepton_pz + met_pz_solution_1), (lepton_energy + met_energy_solution_1)  );
        L_w_lep[1].SetPxPyPzE(  (lepton_px + met_px), (lepton_py + met_py), (lepton_pz + met_pz_solution_2), (lepton_energy + met_energy_solution_2)  );

        L_bw_lep[0] = L_b_lep + L_w_lep[0];
        L_bw_lep[1] = L_b_lep + L_w_lep[1];

        h[hist_leptonic_w_candidate_solution1_pt]->Fill(L_w_lep[0].Pt(), isData ? 1. : NormalizationFactor);
        h[hist_leptonic_w_candidate_solution1_eta]->Fill(L_w_lep[0].Eta(), isData ? 1. : NormalizationFactor);
        h[hist_leptonic_w_candidate_solution1_mass]->Fill(L_w_lep[0].M(), isData ? 1. : NormalizationFactor);
        h[hist_leptonic_top_tbw_solution1_pt]->Fill(L_bw_lep[0].Pt(), isData ? 1. : NormalizationFactor);
        h[hist_leptonic_top_tbw_solution1_eta]->Fill(L_bw_lep[0].Eta(), isData ? 1. : NormalizationFactor);
        h[hist_leptonic_top_tbw_solution1_mass]->Fill(L_bw_lep[0].M(), isData ? 1. : NormalizationFactor);

        h[hist_leptonic_w_candidate_solution2_pt]->Fill(L_w_lep[1].Pt(), isData ? 1. : NormalizationFactor);
        h[hist_leptonic_w_candidate_solution2_eta]->Fill(L_w_lep[1].Eta(), isData ? 1. : NormalizationFactor);
        h[hist_leptonic_w_candidate_solution2_mass]->Fill(L_w_lep[1].M(), isData ? 1. : NormalizationFactor);
        h[hist_leptonic_top_tbw_solution2_pt]->Fill(L_bw_lep[1].Pt(), isData ? 1. : NormalizationFactor);
        h[hist_leptonic_top_tbw_solution2_eta]->Fill(L_bw_lep[1].Eta(), isData ? 1. : NormalizationFactor);
        h[hist_leptonic_top_tbw_solution2_mass]->Fill(L_bw_lep[1].M(), isData ? 1. : NormalizationFactor);
//
//        //## comments: further MC truth study is needed{{{
//        if(L_w_lep[1].M()>90){
//            printf("[check - wboson] mass = %f\n", L_w_lep[1].M());
//            printf("[check - wboson] coefficient_D = %f\n", coefficient_D);
//            printf("[check - wboson] lepton px = %6.2f, ", lepton.Px());
//            printf("py = %6.2f, ", lepton.Py());
//            printf("pz = %6.2f, ", lepton.Pz());
//            printf("energy = %6.2f\n", lepton.E());
//            printf("[check - wboson] metinf px = %6.2f, ", met_px);
//            printf("py = %6.2f, ", met_py);
//            printf("pz = %6.2f, ", met_pz_solution_2);
//            printf("energy = %6.2f\n", met_energy_solution_2);
//            printf("[check - wboson] wboson px = %6.2f, ", L_w_lep[1].Px());
//            printf("py = %6.2f, ", L_w_lep[1].Py());
//            printf("pz = %6.2f, ", L_w_lep[1].Pz());
//            printf("energy = %6.2f\n", L_w_lep[1].E());
//        }
//        //}}}
//
        //}}}
        //--- M1{{{
        double M1;
        int index_q;
        TLorentzVector jet_q;
        TLorentzVector top_fcnh = GetBestM1(M1, treeReader.num_jets, index_bjet, index_jet_chi2_modified, diphoton, Jets, index_q, jet_q);
        if(M1 != -999){
            counter_M1_exists_leptonic += 1;
            //-----
            h[hist_top_tqh_pt]->Fill(top_fcnh.Pt(), isData ? 1. : NormalizationFactor);//exclude event without suitable candidate
            h[hist_top_tqh_eta]->Fill(top_fcnh.Eta(), isData ? 1. : NormalizationFactor);//exclude event without suitable candidate
            h[hist_top_tqh_mass]->Fill(top_fcnh.M(), isData ? 1. : NormalizationFactor);//exclude event without suitable candidate
            //-----
            h[hist_jetq_pt] -> Fill(jet_q.Pt(), isData ? 1. : NormalizationFactor);
            h[hist_jetq_eta] -> Fill(jet_q.Eta(), isData ? 1. : NormalizationFactor);
            h[hist_jetq_phi] -> Fill(jet_q.Phi(), isData ? 1. : NormalizationFactor);
            h[hist_jetq_energy] -> Fill(jet_q.E(), isData ? 1. : NormalizationFactor);
            h[hist_jetq_btag_score] -> Fill(Jets_btag_score[index_q], isData ? 1. : NormalizationFactor);
            h[hist_jetq_CvsL_score] -> Fill(Jets_CvsL_score[index_q], isData ? 1. : NormalizationFactor);
            h[hist_jetq_CvsB_score] -> Fill(Jets_CvsB_score[index_q], isData ? 1. : NormalizationFactor);
            h[hist_jetq_diphoton_deltaR] -> Fill(jet_q.DeltaR(diphoton), isData ? 1. : NormalizationFactor);
            if(treeReader.num_leptons>0){
                for(int i=0; i<treeReader.num_leptons; ++i){
                    double delta_R = jet_q.DeltaR(Leptons.at(i));
                    h[hist_jetq_lepton_deltaR] -> Fill(delta_R, isData ? 1. : NormalizationFactor);
                }
            }
            tree_jet2_pt   = jet_q.Pt();
            tree_jet2_eta  = jet_q.Eta();
            tree_jet2_btag = Jets_btag_score[index_q];
            tree_jet2_CvsL = Jets_CvsL_score[index_q];
            tree_jet2_CvsB = Jets_CvsB_score[index_q];
        }

        //}}}
        //deltaR {{{
        double deltaR;
        if(M1 != -999){ deltaR = diphoton.DeltaR(jet_q)                       ; h[hist_deltaR_qH] -> Fill(deltaR, isData ? 1. : NormalizationFactor)            ; }
        deltaR = leading_photon.DeltaR(subleading_photon)                     ; h[hist_deltaR_photon_photon] -> Fill(deltaR, isData ? 1. : NormalizationFactor) ;
        if(M1 != -999){ deltaR = top_fcnh.DeltaR(L_bw_lep[1])                 ; h[hist_deltaR_top_top] -> Fill(deltaR, isData ? 1. : NormalizationFactor)       ; }
        deltaR = diphoton.DeltaR(L_bw_lep[1])                                 ; h[hist_deltaR_tH] -> Fill(deltaR, isData ? 1. : NormalizationFactor)            ;
        deltaR = diphoton.DeltaR(L_w_lep[1])                                  ; h[hist_deltaR_HW] -> Fill(deltaR, isData ? 1. : NormalizationFactor)            ;
        deltaR = bjet.DeltaR(L_w_lep[1])                                      ; h[hist_deltaR_bW] -> Fill(deltaR, isData ? 1. : NormalizationFactor)            ;
        //deltaR = Jets[0].DeltaR(Jets[1])                                      ; h[hist_deltaR_jet1_jet2] -> Fill(deltaR, isData ? 1. : NormalizationFactor)     ;
        //}}}
        // store tree parameters{{{
        tree_neutrino_pz = met_pz_solution_2;
        tree_top_pt = L_bw_lep[1].Pt();
        tree_top_eta = L_bw_lep[1].Eta();
        tree_top_mass = L_bw_lep[1].M();
        tree_tH_deltaR = diphoton.DeltaR(L_bw_lep[1]);
        //}}}
        */

        //TopKinFit Method{{{
        // set up input parameters{{{
        std::vector<float> BJetPt;
        std::vector<float> BJetEta;
        std::vector<float> BJetPhi;
        std::vector<float> BJetE;
        std::vector<float> NonBJetFilteredPt;
        std::vector<float> NonBJetFilteredEta;
        std::vector<float> NonBJetFilteredPhi;
        std::vector<float> NonBJetFilteredE;
        std::vector<float> ElectronPt;
        std::vector<float> ElectronEta;
        std::vector<float> ElectronPhi;
        std::vector<float> ElectronE;
        std::vector<float> MuonPt;
        std::vector<float> MuonEta;
        std::vector<float> MuonPhi;
        std::vector<float> MuonE;
        float MetRecPx = treeReader.MetInfo_Px;
        float MetRecPy = treeReader.MetInfo_Py;

        //exactly one tight bjet
        BJetPt.push_back(bjet.Pt());
        BJetEta.push_back(bjet.Eta());
        BJetPhi.push_back(bjet.Phi());
        BJetE.push_back(bjet.E());
        for(std::size_t i=0; i<Jets.size(); ++i){
            if(i==index_bjet) continue;
            NonBJetFilteredPt.push_back(Jets[i].Pt());
            NonBJetFilteredEta.push_back(Jets[i].Eta());
            NonBJetFilteredPhi.push_back(Jets[i].Phi());
            NonBJetFilteredE.push_back(Jets[i].E());
        }
        for(std::size_t i=0; i<Electrons.size(); ++i){
            ElectronPt.push_back(Electrons[i].Pt());
            ElectronEta.push_back(Electrons[i].Eta());
            ElectronPhi.push_back(Electrons[i].Phi());
            ElectronE.push_back(Electrons[i].E());
        }
        for(std::size_t i=0; i<Muons.size(); ++i){
            MuonPt.push_back(Muons[i].Pt());
            MuonEta.push_back(Muons[i].Eta());
            MuonPhi.push_back(Muons[i].Phi());
            MuonE.push_back(Muons[i].E());
        }

        // Pass reconstructed objects to the tool (std::vector<float> for all objects except float for Met)
        kf->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
        kf->SetNonBJet(NonBJetFilteredPt,NonBJetFilteredEta,NonBJetFilteredPhi,NonBJetFilteredE);
        kf->SetElectron(ElectronPt,ElectronEta,ElectronPhi,ElectronE);
        kf->SetMuon(MuonPt,MuonEta,MuonPhi,MuonE);
        kf->SetMet(MetRecPx,MetRecPy);
        //}}}
        //run{{{
        kf->Run(); // Run the tool
        int NPerm = kf->GetNPerm(); // Get number of permutations
        std::vector<float> NuPz, disc;

        for(int ip=0;ip<NPerm;ip++) // Loop over permutations - already sorted in likelihood value from min to max
        {
            //printf("[INFO-kinfit] disc = %f\n", disc);
            disc.push_back( kf->GetDisc(ip) ); // Get minimized likelihood value
            // Get reconstructed neutrino
            float NuPx = kf->GetNuPx(ip,0);
            float NuPy = kf->GetNuPy(ip,0);
            NuPz.push_back( kf->GetNuPz(ip,0) );
        }
        //}}}

        if(disc[0] > 100000.) counter_irregular_disc += 1;

        //reconstruct W, top
        TLorentzVector L_met_topKinFit;
        TLorentzVector L_w_topKinFit;
        TLorentzVector L_bw_topKinFit;
        float met_pz_topKinFit = NuPz[0];
        float met_energy_topKinFit = TMath::Sqrt(met_pt*met_pt + met_pz_topKinFit*met_pz_topKinFit);
        L_met_topKinFit.SetPxPyPzE( met_px, met_py, met_pz_topKinFit, met_energy_topKinFit );
        L_w_topKinFit.SetPxPyPzE( (lepton_px + met_px), (lepton_py + met_py), (lepton_pz + met_pz_topKinFit), (lepton_energy + met_energy_topKinFit) );
        L_bw_topKinFit = bjet + L_w_topKinFit;

        bool canFoundSolution_topKinFit = !(disc[0] > 100000.);
        bool is_reg_and_positive = canFoundSolution_topKinFit && met_pz_topKinFit >= 0;

        //}}}

        tree_neutrino_pz = L_met_topKinFit.Pz();
        tree_top_pt = L_bw_topKinFit.Pt();
        tree_top_eta = L_bw_topKinFit.Eta();
        tree_top_mass = L_bw_topKinFit.M();
        tree_tH_deltaR = diphoton.DeltaR(L_bw_topKinFit);

        } // end of else
        //}}}
        counter_selected_events += 1;
        myAnalysisTree->Fill();
    }//end of event loop

    //### counters, yields, plots, close{{{
    printf("[INFO] total entry after  selection    = %d\n", counter_selected_events);
    if(bool_isHadronic)
        printf("[INFO] num of events with an existing M1 (hadronic) = %d / %d \n", counter_M1_exists_hadronic, counter_selected_events);
    if(bool_isLeptonic){
        printf("[INFO] num of events with an existing M1 (leptonic) = %d / %d \n", counter_M1_exists_leptonic, counter_selected_events);
        printf("[INFO] num of negative coeff D = %d / %d \n", counter_coeff_D_isNegative, counter_selected_events);
    }
    //==================================================//
    //---------------------  Yields  -------------------//
    //==================================================//
    double yields = h[hist_EvtInfo_NPu] -> Integral();
    printf("[INFO] Expected yields = %f\n", yields);
    //==================================================//
    //--------------------- MakePlot -------------------//
    //==================================================//
    // Write/Close hist file
    fout -> cd();
    TCanvas *c1 = new TCanvas("c1", "c1", 700, 800);
    for(int i=0; i<totalHistNum; ++i) MakePlots(c1, h[i], Form("%s/%s.png", output_dir, histNames[i].c_str()));
    fout -> Close();

    // Write/Close tree file
    f_tree -> cd();
    myAnalysisTree -> Write();
    f_tree -> Close();

    // Close input file
    f_mcpu -> Close();
    fin  -> Close();
    //}}}

}

//### Main function{{{
int main(int argc, char *argv[]){
    char input_file[512]; sprintf(input_file, "%s", argv[1]); printf("[INFO] input_file  = %s\n", input_file);
    char output_file[512]; sprintf(output_file, "%s", argv[2]); printf("[INFO] output_file = %s\n", output_file);
    char output_tree[512]; sprintf(output_tree, "%s", argv[3]); printf("[INFO] output_tree = %s\n", output_tree);
    char dataset[512]; sprintf(dataset, "%s", argv[4]); printf("[INFO] dataset     = %s\n", dataset);
    char output_dir[512];  sprintf(output_dir, "%s", argv[5]); printf("[INFO] output_dir  = %s\n", output_dir);
    char channel[512];  sprintf(channel, "%s", argv[6]); printf("[INFO] channel     = %s\n", channel);
    Selection(input_file, output_file, output_tree, dataset, output_dir, channel);
    return 1;
}
//}}}
// ### GetBestM1{{{
TLorentzVector GetBestM1(double &M1, int num_jets, int index_bjet, std::vector<int> index_jet, TLorentzVector diphoton, std::vector<TLorentzVector> Jets, int &index_q, TLorentzVector &jet_q){
    //record all the possible combinations
    std::vector<int> indices;
    std::vector<double> top_fcnh_chi2;
    std::vector<TLorentzVector> top_fcnh_candidates;
    for(int i=0; i<num_jets; ++i){
        if(i==index_bjet || i==index_jet[0] || i==index_jet[1]) continue;//skip the jets for bjj
        TLorentzVector top_fcnh_tmp = diphoton + Jets[i];
        double chi2 = (top_fcnh_tmp.M() - top_quark_mass) * (top_fcnh_tmp.M() - top_quark_mass);
        indices.push_back(i);
        top_fcnh_chi2.push_back(chi2);
        top_fcnh_candidates.push_back(top_fcnh_tmp);
    }
    //choose the candidate with Mass closest to top mass
    int index_M1_min_chi2 = -999;
    double chi2_min_M1 = 40000;//200^2 > 178^2 ~ (x-top_mass)^2 //NOTE: will bound the range of M1!!!
    //pick the wanted one (with index)
    for(int i=0; i<top_fcnh_candidates.size(); ++i){ if(top_fcnh_chi2[i]<chi2_min_M1){ index_M1_min_chi2 = i; chi2_min_M1 = top_fcnh_chi2[i];} }
    TLorentzVector null;
    if(index_M1_min_chi2 == -999){ M1 = -999; return null;}//No suitable candidate
    else{
        index_q = indices[index_M1_min_chi2];
        jet_q = Jets[index_M1_min_chi2];
        M1 = top_fcnh_candidates[index_M1_min_chi2].M();
        return top_fcnh_candidates[index_M1_min_chi2];
    }
}
//}}}
//### Functions{{{
void MakePlots(TCanvas *c1, TH1D* hist, const char* outputFile){
    hist->Draw();
    hist->SetYTitle("Entries");
    hist->GetYaxis()->SetTitleOffset(1.4);
    hist->Write();
    c1->SaveAs(outputFile);
}
bool isThisDataOrNot(char* dataset){
    if((string)dataset == "DoubleEG_B") return true;
    if((string)dataset == "DoubleEG_C") return true;
    if((string)dataset == "DoubleEG_D") return true;
    if((string)dataset == "DoubleEG_E") return true;
    if((string)dataset == "DoubleEG_F") return true;
    return false;
}
bool isThisMCsignal(char* dataset){
    if((string)dataset == "ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    return false;
}
double Chi2_calculator_modified(double w_mass, double t_mass){
    TVectorD vec_mass(2);
    vec_mass(0) = w_mass - w_boson_mass;
    vec_mass(1) = t_mass - top_quark_mass;
    TMatrixD matrix(2,2);
    matrix(0,0) = 331.84; matrix(0,1) = 374.95;
    matrix(1,0) = 374.95; matrix(1,1) = 965.17;

    //--- ST ---//
    //matrix(0,0) = 305.98; matrix(0,1) = 323.97;
    //matrix(1,0) = 323.97; matrix(1,1) = 787.07;
    return matrix.Invert()*vec_mass*vec_mass;
}
//}}}
//### Class functions{{{
void myTreeClass::InitTree(const char* treeName){
    mytree = new TTree(treeName, treeName);
}
void myTreeClass::AddRootFile(TFile* input){
    mytree = (TTree*) input->Get("mytree");
    printf("[INFO] myTreeClass::AddRootFile : Finished!\n");
}
int myTreeClass::GetEntries(void){
    return mytree->GetEntries();
}
TTree* myTreeClass::GetTTree(void){ 
    return mytree;
}
void myTreeClass::SetBranchAddresses(){
    mytree -> SetBranchAddress("EvtInfo_totalEntry_before_preselection", &EvtInfo_totalEntry_before_preselection);
    mytree -> SetBranchAddress("EvtInfo_NormalizationFactor_lumi", &EvtInfo_NormalizationFactor_lumi);
    mytree -> SetBranchAddress("EvtInfo_NPu", &EvtInfo_NPu);
    mytree -> SetBranchAddress("EvtInfo_Rho", &EvtInfo_Rho);
    mytree -> SetBranchAddress("EvtInfo_NVtx", &EvtInfo_NVtx);
    mytree -> SetBranchAddress("EvtInfo_genweight", &EvtInfo_genweight);
    //------------------------
    mytree -> SetBranchAddress("DiPhoInfo_mass", &DiPhoInfo_mass);
    mytree -> SetBranchAddress("DiPhoInfo_pt", &DiPhoInfo_pt);
    mytree -> SetBranchAddress("DiPhoInfo_eta", &DiPhoInfo_eta);
    mytree -> SetBranchAddress("DiPhoInfo_phi", &DiPhoInfo_phi);
    mytree -> SetBranchAddress("DiPhoInfo_energy", &DiPhoInfo_energy);
    mytree -> SetBranchAddress("DiPhoInfo_leadPt", &DiPhoInfo_leadPt);
    mytree -> SetBranchAddress("DiPhoInfo_leadEta", &DiPhoInfo_leadEta);
    mytree -> SetBranchAddress("DiPhoInfo_leadPhi", &DiPhoInfo_leadPhi);
    mytree -> SetBranchAddress("DiPhoInfo_leadE", &DiPhoInfo_leadE);
    mytree -> SetBranchAddress("DiPhoInfo_leadhoe", &DiPhoInfo_leadhoe);
    mytree -> SetBranchAddress("DiPhoInfo_leadIDMVA", &DiPhoInfo_leadIDMVA);
    mytree -> SetBranchAddress("DiPhoInfo_subleadPt", &DiPhoInfo_subleadPt);
    mytree -> SetBranchAddress("DiPhoInfo_subleadEta", &DiPhoInfo_subleadEta);
    mytree -> SetBranchAddress("DiPhoInfo_subleadPhi", &DiPhoInfo_subleadPhi);
    mytree -> SetBranchAddress("DiPhoInfo_subleadE", &DiPhoInfo_subleadE);
    mytree -> SetBranchAddress("DiPhoInfo_subleadhoe", &DiPhoInfo_subleadhoe);
    mytree -> SetBranchAddress("DiPhoInfo_subleadIDMVA", &DiPhoInfo_subleadIDMVA);
    //------------------------
    mytree -> SetBranchAddress("ElecInfo_Size", &ElecInfo_Size);
    mytree -> SetBranchAddress("MuonInfo_Size", &MuonInfo_Size);
    mytree -> SetBranchAddress("num_leptons", &num_leptons);// # of selected objects.
    mytree -> SetBranchAddress("num_electrons", &num_electrons);// # of selected objects.
    mytree -> SetBranchAddress("num_muons", &num_muons);// # of selected objects.
    mytree -> SetBranchAddress("ElecInfo_electron_charge", &ElecInfo_electron_charge_selection);
    mytree -> SetBranchAddress("ElecInfo_electron_pt", &ElecInfo_electron_pt_selection);
    mytree -> SetBranchAddress("ElecInfo_electron_eta", &ElecInfo_electron_eta_selection);
    mytree -> SetBranchAddress("ElecInfo_electron_phi", &ElecInfo_electron_phi_selection);
    mytree -> SetBranchAddress("ElecInfo_electron_energy", &ElecInfo_electron_energy_selection);
    mytree -> SetBranchAddress("ElecInfo_electron_diphoton_deltaR", &ElecInfo_electron_diphoton_deltaR_selection);
    mytree -> SetBranchAddress("MuonInfo_muon_charge", &MuonInfo_muon_charge_selection);
    mytree -> SetBranchAddress("MuonInfo_muon_pt", &MuonInfo_muon_pt_selection);
    mytree -> SetBranchAddress("MuonInfo_muon_eta", &MuonInfo_muon_eta_selection);
    mytree -> SetBranchAddress("MuonInfo_muon_phi", &MuonInfo_muon_phi_selection);
    mytree -> SetBranchAddress("MuonInfo_muon_energy", &MuonInfo_muon_energy_selection);
    mytree -> SetBranchAddress("MuonInfo_muon_diphoton_deltaR", &MuonInfo_muon_diphoton_deltaR_selection);
    //------------------------
    //mytree -> SetBranchAddress("GenPartInfo_size", &GenPartInfo_size);
    //mytree -> SetBranchAddress("GenPartInfo_gen_Pt", &GenPartInfo_gen_Pt_selection);
    //mytree -> SetBranchAddress("GenPartInfo_gen_Eta", &GenPartInfo_gen_Eta_selection);
    //mytree -> SetBranchAddress("GenPartInfo_gen_Phi", &GenPartInfo_gen_Phi_selection);
    //mytree -> SetBranchAddress("GenPartInfo_gen_Mass", &GenPartInfo_gen_Mass_selection);
    //mytree -> SetBranchAddress("GenPartInfo_gen_PdgID", &GenPartInfo_gen_PdgID_selection);
    //mytree -> SetBranchAddress("GenPartInfo_gen_Status", &GenPartInfo_gen_Status_selection);
    //mytree -> SetBranchAddress("GenPartInfo_gen_nMo", &GenPartInfo_gen_nMo_selection);
    //mytree -> SetBranchAddress("GenPartInfo_gen_nDa", &GenPartInfo_gen_nDa_selection);
    //------------------------
    mytree -> SetBranchAddress("jets_size", &jets_size);
    mytree -> SetBranchAddress("num_jets", &num_jets);
    mytree -> SetBranchAddress("JetInfo_jet_pt", &JetInfo_jet_pt_selection);
    mytree -> SetBranchAddress("JetInfo_jet_eta", &JetInfo_jet_eta_selection);
    mytree -> SetBranchAddress("JetInfo_jet_phi", &JetInfo_jet_phi_selection);
    mytree -> SetBranchAddress("JetInfo_jet_energy", &JetInfo_jet_energy_selection);
    mytree -> SetBranchAddress("JetInfo_jet_diphoton_deltaR", &JetInfo_jet_diphoton_deltaR_selection);
    mytree -> SetBranchAddress("JetInfo_jet_pfDeepCSVJetTags_probb", &JetInfo_jet_pfDeepCSVJetTags_probb_selection);
    mytree -> SetBranchAddress("JetInfo_jet_pfDeepCSVJetTags_probbb", &JetInfo_jet_pfDeepCSVJetTags_probbb_selection);
    mytree -> SetBranchAddress("JetInfo_jet_pfDeepCSVJetTags_probc", &JetInfo_jet_pfDeepCSVJetTags_probc_selection);
    mytree -> SetBranchAddress("JetInfo_jet_pfDeepCSVJetTags_probudsg", &JetInfo_jet_pfDeepCSVJetTags_probudsg_selection);
    //------------------------
    mytree -> SetBranchAddress("num_bjets", &num_bjets);// # of selected objects.
    mytree -> SetBranchAddress("JetInfo_leading_bjet_pt", &JetInfo_leading_bjet_pt_selection);
    mytree -> SetBranchAddress("JetInfo_leading_bjet_eta", &JetInfo_leading_bjet_eta_selection);
    mytree -> SetBranchAddress("JetInfo_leading_bjet_phi", &JetInfo_leading_bjet_phi_selection);
    mytree -> SetBranchAddress("JetInfo_leading_bjet_energy", &JetInfo_leading_bjet_energy_selection);
    //------------------------
    mytree -> SetBranchAddress("MetInfo_Pt", &MetInfo_Pt);
    mytree -> SetBranchAddress("MetInfo_Phi", &MetInfo_Phi);
    mytree -> SetBranchAddress("MetInfo_Px", &MetInfo_Px);
    mytree -> SetBranchAddress("MetInfo_Py", &MetInfo_Py);
    mytree -> SetBranchAddress("MetInfo_SumET", &MetInfo_SumET);
    //------------------------
    printf("[INFO] myTreeClass::SetBranchAddresses : Finished!\n");
}
myParameters::myParameters(){
    GenPartInfo_gen_Pt_selection = new std::vector<float>;
    GenPartInfo_gen_Eta_selection = new std::vector<float>;
    GenPartInfo_gen_Phi_selection = new std::vector<float>;
    GenPartInfo_gen_Mass_selection = new std::vector<float>;
    GenPartInfo_gen_PdgID_selection = new std::vector<int>;
    GenPartInfo_gen_Status_selection = new std::vector<int>;
    GenPartInfo_gen_nMo_selection = new std::vector<int>;
    GenPartInfo_gen_nDa_selection = new std::vector<int>;
    JetInfo_jet_pt_selection = new std::vector<float>;
    JetInfo_jet_eta_selection = new std::vector<float>;
    JetInfo_jet_phi_selection = new std::vector<float>;
    JetInfo_jet_energy_selection = new std::vector<float>;
    JetInfo_jet_diphoton_deltaR_selection = new std::vector<float>;
    JetInfo_jet_pfDeepCSVJetTags_probb_selection = new std::vector<float>;
    JetInfo_jet_pfDeepCSVJetTags_probbb_selection = new std::vector<float>;
    JetInfo_jet_pfDeepCSVJetTags_probc_selection = new std::vector<float>;
    JetInfo_jet_pfDeepCSVJetTags_probudsg_selection = new std::vector<float>;
    JetInfo_leading_bjet_pt_selection = new std::vector<float>;
    JetInfo_leading_bjet_eta_selection = new std::vector<float>;
    JetInfo_leading_bjet_phi_selection = new std::vector<float>;
    JetInfo_leading_bjet_energy_selection = new std::vector<float>;
    ElecInfo_electron_charge_selection = new std::vector<int>;
    ElecInfo_electron_pt_selection = new std::vector<float>;
    ElecInfo_electron_eta_selection = new std::vector<float>;
    ElecInfo_electron_phi_selection = new std::vector<float>;
    ElecInfo_electron_energy_selection = new std::vector<float>;
    ElecInfo_electron_diphoton_deltaR_selection = new std::vector<float>;
    MuonInfo_muon_charge_selection = new std::vector<int>;
    MuonInfo_muon_pt_selection = new std::vector<float>;
    MuonInfo_muon_eta_selection = new std::vector<float>;
    MuonInfo_muon_phi_selection = new std::vector<float>;
    MuonInfo_muon_energy_selection = new std::vector<float>;
    MuonInfo_muon_diphoton_deltaR_selection = new std::vector<float>;
}
myParameters::~myParameters(){
    delete GenPartInfo_gen_Pt_selection;
    delete GenPartInfo_gen_Eta_selection;
    delete GenPartInfo_gen_Phi_selection;
    delete GenPartInfo_gen_Mass_selection;
    delete GenPartInfo_gen_PdgID_selection;
    delete GenPartInfo_gen_Status_selection;
    delete GenPartInfo_gen_nMo_selection;
    delete GenPartInfo_gen_nDa_selection;
    delete JetInfo_jet_pt_selection;
    delete JetInfo_jet_eta_selection;
    delete JetInfo_jet_phi_selection;
    delete JetInfo_jet_energy_selection;
    delete JetInfo_jet_diphoton_deltaR_selection;
    delete JetInfo_jet_pfDeepCSVJetTags_probb_selection;
    delete JetInfo_jet_pfDeepCSVJetTags_probbb_selection;
    delete JetInfo_jet_pfDeepCSVJetTags_probc_selection;
    delete JetInfo_jet_pfDeepCSVJetTags_probudsg_selection;
    delete JetInfo_leading_bjet_pt_selection;
    delete JetInfo_leading_bjet_eta_selection;
    delete JetInfo_leading_bjet_phi_selection;
    delete JetInfo_leading_bjet_energy_selection;
    delete ElecInfo_electron_charge_selection;
    delete ElecInfo_electron_pt_selection;
    delete ElecInfo_electron_eta_selection;
    delete ElecInfo_electron_phi_selection;
    delete ElecInfo_electron_energy_selection;
    delete ElecInfo_electron_diphoton_deltaR_selection;
    delete MuonInfo_muon_charge_selection;
    delete MuonInfo_muon_pt_selection;
    delete MuonInfo_muon_eta_selection;
    delete MuonInfo_muon_phi_selection;
    delete MuonInfo_muon_energy_selection;
    delete MuonInfo_muon_diphoton_deltaR_selection;
}
flashggStdTreeParameters::flashggStdTreeParameters(){
    JetInfo_Pt = new std::vector<float>;
    JetInfo_Eta = new std::vector<float>;
    JetInfo_Phi = new std::vector<float>;
    JetInfo_Mass = new std::vector<float>;
    JetInfo_Energy = new std::vector<float>;
    JetInfo_pfDeepCSVJetTags_probb = new std::vector<float>;
    JetInfo_pfDeepCSVJetTags_probbb = new std::vector<float>;
    JetInfo_pfDeepCSVJetTags_probc = new std::vector<float>;
    JetInfo_pfDeepCSVJetTags_probudsg = new std::vector<float>;
    //------------------------
    ElecInfo_Charge = new std::vector<int>;
    ElecInfo_Pt = new std::vector<float>;
    ElecInfo_Eta = new std::vector<float>;
    ElecInfo_Phi = new std::vector<float>;
    ElecInfo_Energy = new std::vector<float>;
    ElecInfo_EtaSC = new std::vector<float>;
    ElecInfo_PhiSC = new std::vector<float>;
    ElecInfo_GsfTrackDz = new std::vector<float>;
    ElecInfo_GsfTrackDxy = new std::vector<float>;
    ElecInfo_EGMCutBasedIDVeto = new std::vector<bool>;
    ElecInfo_EGMCutBasedIDLoose = new std::vector<bool>;
    ElecInfo_EGMCutBasedIDMedium = new std::vector<bool>;
    ElecInfo_EGMCutBasedIDTight = new std::vector<bool>;
    ElecInfo_fggPhoVeto = new std::vector<bool>;
    //ElecInfo_tmpPhoVeto = new std::vector<bool>;
    ElecInfo_EnergyCorrFactor = new std::vector<float>;
    ElecInfo_EnergyPostCorrErr = new std::vector<float>;
    ElecInfo_EnergyPostCorrScaleUp = new std::vector<float>;
    ElecInfo_EnergyPostCorrScaleDown = new std::vector<float>;
    ElecInfo_EnergyPostCorrSmearUp = new std::vector<float>;
    ElecInfo_EnergyPostCorrSmearDown = new std::vector<float>;
    ElecInfo_GenMatch = new std::vector<bool>;
    ElecInfo_GenPdgID = new std::vector<int>;
    ElecInfo_GenPt = new std::vector<float>;
    ElecInfo_GenEta = new std::vector<float>;
    ElecInfo_GenPhi = new std::vector<float>;
    //------------------------
    MuonInfo_Charge = new std::vector<int>;
    MuonInfo_MuonType = new std::vector<float>;
    MuonInfo_Pt = new std::vector<float>;
    MuonInfo_Eta = new std::vector<float>;
    MuonInfo_Phi = new std::vector<float>;
    MuonInfo_Energy = new std::vector<float>;
    MuonInfo_BestTrackDz = new std::vector<float>;
    MuonInfo_BestTrackDxy = new std::vector<float>;
    MuonInfo_PFIsoDeltaBetaCorrR04 = new std::vector<float>;
    MuonInfo_TrackerBasedIsoR03 = new std::vector<float>;
    MuonInfo_CutBasedIdMedium = new std::vector<bool>;
    MuonInfo_CutBasedIdTight = new std::vector<bool>;
    MuonInfo_GenMatch = new std::vector<bool>;
    MuonInfo_GenPdgID = new std::vector<int>;
    MuonInfo_GenPt = new std::vector<float>;
    MuonInfo_GenEta = new std::vector<float>;
    MuonInfo_GenPhi = new std::vector<float>;
    //------------------------
}
flashggStdTreeParameters::~flashggStdTreeParameters(){
    delete JetInfo_Pt;
    delete JetInfo_Eta;
    delete JetInfo_Phi;
    delete JetInfo_Mass;
    delete JetInfo_Energy;
    delete JetInfo_pfDeepCSVJetTags_probb;
    delete JetInfo_pfDeepCSVJetTags_probbb;
    delete JetInfo_pfDeepCSVJetTags_probc;
    delete JetInfo_pfDeepCSVJetTags_probudsg;
    //------------------------
    delete ElecInfo_Charge;
    delete ElecInfo_Pt;
    delete ElecInfo_Eta;
    delete ElecInfo_Phi;
    delete ElecInfo_Energy;
    delete ElecInfo_EtaSC;
    delete ElecInfo_PhiSC;
    delete ElecInfo_GsfTrackDz;
    delete ElecInfo_GsfTrackDxy;
    delete ElecInfo_EGMCutBasedIDVeto;
    delete ElecInfo_EGMCutBasedIDLoose;
    delete ElecInfo_EGMCutBasedIDMedium;
    delete ElecInfo_EGMCutBasedIDTight;
    delete ElecInfo_fggPhoVeto;
    //delete ElecInfo_tmpPhoVeto;
    delete ElecInfo_EnergyCorrFactor;
    delete ElecInfo_EnergyPostCorrErr;
    delete ElecInfo_EnergyPostCorrScaleUp;
    delete ElecInfo_EnergyPostCorrScaleDown;
    delete ElecInfo_EnergyPostCorrSmearUp;
    delete ElecInfo_EnergyPostCorrSmearDown;
    delete ElecInfo_GenMatch;
    delete ElecInfo_GenPdgID;
    delete ElecInfo_GenPt;
    delete ElecInfo_GenEta;
    delete ElecInfo_GenPhi;
    //------------------------
    delete MuonInfo_Charge;
    delete MuonInfo_MuonType;
    delete MuonInfo_Pt;
    delete MuonInfo_Eta;
    delete MuonInfo_Phi;
    delete MuonInfo_Energy;
    delete MuonInfo_BestTrackDz;
    delete MuonInfo_BestTrackDxy;
    delete MuonInfo_PFIsoDeltaBetaCorrR04;
    delete MuonInfo_TrackerBasedIsoR03;
    delete MuonInfo_CutBasedIdMedium;
    delete MuonInfo_CutBasedIdTight;
    delete MuonInfo_GenMatch;
    delete MuonInfo_GenPdgID;
    delete MuonInfo_GenPt;
    delete MuonInfo_GenEta;
    delete MuonInfo_GenPhi;
    //------------------------
}
//}}}
