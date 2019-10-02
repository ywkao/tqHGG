//***************************************************************************
//
// FileName    : selection.cpp
// Purpose     : Develop for top FCNH with H to two photons analysis
// Description : Applying event selection & Preparing histograms for individual dataset.
// Deetail     : PU reweighting, leptonic/hadronic channels, hists, (top reconstruction).
// Author      : Yu-Wei Kao [ykao@cern.ch]
//
//***************************************************************************
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
using namespace std;

bool bool_isHadronic;

//--- control bjet selection ---//
bool bool_bjet_is_loose  = false;
bool bool_bjet_is_medium = true;
bool bool_bjet_is_tight  = false;
bool bool_num_bjets_is_exactly_one = false;
bool bool_num_bjets_is_atleast_one = true;

void Selection(char* input_file, char* output_file, char* dataset, char* output_dir, char* channel){
    bool isData = isThisDataOrNot(dataset);
    //### I/O, histograms{{{
    if((string)channel == "hadronic") bool_isHadronic = true; else bool_isHadronic = false;
    if(bool_isHadronic) printf("[CHECK] isHadronic!\n");
    if(!bool_isHadronic) printf("[CHECK] isLeptonic!\n");
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
    //==================================================//
    //-------------------- Event Loop ------------------//
    //==================================================//
    int counter_bjet_is_bquark = 0;
    int nentries = treeReader.GetEntries();
    for(int ientry=0; ientry<nentries; ientry++){
        TTree* tmp = treeReader.GetTTree();
        tmp->GetEntry(ientry);
        if(ientry==0) printf("[INFO] N_entries = %d/%d\n", nentries, treeReader.EvtInfo_totalEntry_before_preselection);
        if((ientry+1)%10000==0 || (ientry+1)==nentries) printf("ientry = %d\r", ientry);
        //### PU & Event selection{{{
        //========= PU Reweighting =========//
        bool apply_PU_reweighting = true;
        double PU_reweighting_factor = apply_PU_reweighting ? h_pu_reweighting_factor->GetBinContent((int)treeReader.EvtInfo_NPu+1) : 1.;
        double NormalizationFactor = treeReader.EvtInfo_genweight * treeReader.EvtInfo_NormalizationFactor_lumi * PU_reweighting_factor;
        double NormalizationFactor_wopu = treeReader.EvtInfo_genweight * treeReader.EvtInfo_NormalizationFactor_lumi;
        //Reminder: EvtInfo_NormalizationFactor_lumi = 1000. * Luminosity * CrossSection * BranchingFraction / TotalGenweight;


        //========= Event Selections =========//
        bool passEvent=true;
        //--- Leptonic Channel ---//
        if(!bool_isHadronic){
            if(treeReader.num_leptons<1) passEvent = false;
            if(treeReader.num_jets<1) passEvent = false;
        //--- Hadronic Channel ---//
        } else{
            if(treeReader.num_leptons>0) passEvent = false;
            //if(treeReader.num_jets<3) passEvent = false;
            if(treeReader.num_jets<4) passEvent = false;
        }
        if(!passEvent) continue;
        //--------- check bjets ---------//
        TLorentzVector bjet, jet;
        int num_bjets = 0, index_bjet = -999;

        std::vector<TLorentzVector> bjets_tight, bjets_loose, bjets_medium;
        int num_bjets_tight = 0, num_bjets_loose = 0, num_bjets_medium = 0;
        int index_leading_bjet_tight = -999, index_leading_bjet_loose = -999, index_leading_bjet_medium = -999;

        if(!(treeReader.num_jets<1)){
            for(int i=0; i<treeReader.num_jets; ++i){
                jet.SetPtEtaPhiE(treeReader.JetInfo_jet_pt_selection->at(i),treeReader.JetInfo_jet_eta_selection->at(i),treeReader.JetInfo_jet_phi_selection->at(i),treeReader.JetInfo_jet_energy_selection->at(i));
                if(treeReader.JetInfo_jet_pfDeepCSVJetTags_probb_selection->at(i)+treeReader.JetInfo_jet_pfDeepCSVJetTags_probbb_selection->at(i) >= pfDeepCSVJetTags_loose){
                    num_bjets_loose += 1;
                    bjets_loose.push_back(jet);
                    if(num_bjets_loose == 1) index_leading_bjet_loose = i;
                }
                if(treeReader.JetInfo_jet_pfDeepCSVJetTags_probb_selection->at(i)+treeReader.JetInfo_jet_pfDeepCSVJetTags_probbb_selection->at(i) >= pfDeepCSVJetTags_medium){
                    num_bjets_medium += 1;
                    bjets_medium.push_back(jet);
                    if(num_bjets_medium == 1) index_leading_bjet_medium = i;
                }
                if(treeReader.JetInfo_jet_pfDeepCSVJetTags_probb_selection->at(i)+treeReader.JetInfo_jet_pfDeepCSVJetTags_probbb_selection->at(i) >= pfDeepCSVJetTags_tight){
                    num_bjets_tight += 1;
                    bjets_tight.push_back(jet);
                    if(num_bjets_tight == 1) index_leading_bjet_tight = i;
                }
            }//end of looping jets
        }
        if(bool_bjet_is_loose)  { index_bjet = index_leading_bjet_loose  ; if(index_bjet != -999){ bjet = bjets_loose[0]  ; num_bjets = num_bjets_loose  ;} }
        if(bool_bjet_is_medium) { index_bjet = index_leading_bjet_medium ; if(index_bjet != -999){ bjet = bjets_medium[0] ; num_bjets = num_bjets_medium ;} }
        if(bool_bjet_is_tight)  { index_bjet = index_leading_bjet_tight  ; if(index_bjet != -999){ bjet = bjets_tight[0]  ; num_bjets = num_bjets_tight  ;} }

        if(index_bjet != -999){//at least one bjet!
            h[hist_leading_bjet_pt] -> Fill(bjet.Pt(), isData ? 1. : NormalizationFactor);
            h[hist_leading_bjet_eta] -> Fill(bjet.Eta(), isData ? 1. : NormalizationFactor);
            h[hist_leading_bjet_phi] -> Fill(bjet.Phi(), isData ? 1. : NormalizationFactor);
            h[hist_leading_bjet_energy] -> Fill(bjet.E(), isData ? 1. : NormalizationFactor);
        }
        bool pass_bjets_multiplicity_selection;
        if(bool_num_bjets_is_exactly_one) pass_bjets_multiplicity_selection = num_bjets == 1;
        if(bool_num_bjets_is_atleast_one) pass_bjets_multiplicity_selection = num_bjets >= 1;
        if(!pass_bjets_multiplicity_selection) continue;

        //To be modified later...{{{
        //if(num_bjets != 1) continue;// exactly one
        //if(num_bjets == 0) continue;// at least one
        /* No gen-info yet......
        
        //--- check the rate of bjet being bqurak ---//
        bool is_b_quark = CheckBJetID(bjets_tight[0],\
                          treeReader.GenPartInfo_size, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                          treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
        if(is_b_quark) count_bjet_is_bquark_tight += 1;
        */
        /*
        if(treeReader.num_bjets == 0) continue;
        if(treeReader.num_bjets > 0){
            h[hist_leading_bjet_pt]     -> Fill(treeReader.JetInfo_leading_bjet_pt_selection->at(0), isData ? 1. : NormalizationFactor);
            h[hist_leading_bjet_eta]    -> Fill(treeReader.JetInfo_leading_bjet_eta_selection->at(0), isData ? 1. : NormalizationFactor);
            h[hist_leading_bjet_phi]    -> Fill(treeReader.JetInfo_leading_bjet_phi_selection->at(0), isData ? 1. : NormalizationFactor);
            h[hist_leading_bjet_energy] -> Fill(treeReader.JetInfo_leading_bjet_energy_selection->at(0), isData ? 1. : NormalizationFactor);
        }
        */
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
        //}}}
        //### Leptons{{{
        //--------- Leptons ---------//
        std::vector<TLorentzVector> Leptons;
        h[hist_ElecInfo_Size] -> Fill(treeReader.ElecInfo_Size, isData ? 1. : NormalizationFactor);
        h[hist_MuonInfo_Size] -> Fill(treeReader.MuonInfo_Size, isData ? 1. : NormalizationFactor);
        h[hist_num_leptons] -> Fill(treeReader.num_leptons, isData ? 1. : NormalizationFactor);// # of selected objects.
        h[hist_num_electrons] -> Fill(treeReader.num_electrons, isData ? 1. : NormalizationFactor);// # of selected objects.
        h[hist_num_muons] -> Fill(treeReader.num_muons, isData ? 1. : NormalizationFactor);// # of selected objects.
        if(treeReader.num_electrons>0){
            for(int i=0; i<treeReader.num_electrons; ++i){
                h[hist_ElecInfo_electron_pt] -> Fill(treeReader.ElecInfo_electron_pt_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_ElecInfo_electron_eta] -> Fill(treeReader.ElecInfo_electron_eta_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_ElecInfo_electron_phi] -> Fill(treeReader.ElecInfo_electron_phi_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_ElecInfo_electron_energy] -> Fill(treeReader.ElecInfo_electron_energy_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_ElecInfo_electron_diphoton_deltaR] -> Fill(treeReader.ElecInfo_electron_diphoton_deltaR_selection->at(i), isData ? 1. : NormalizationFactor);
                //------------------------
                h[hist_lepton_pt] -> Fill(treeReader.ElecInfo_electron_pt_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_eta] -> Fill(treeReader.ElecInfo_electron_eta_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_phi] -> Fill(treeReader.ElecInfo_electron_phi_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_energy] -> Fill(treeReader.ElecInfo_electron_energy_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_diphoton_deltaR] -> Fill(treeReader.ElecInfo_electron_diphoton_deltaR_selection->at(i), isData ? 1. : NormalizationFactor);
                //------------------------
                TLorentzVector electron; electron.SetPtEtaPhiE(treeReader.ElecInfo_electron_pt_selection->at(i), treeReader.ElecInfo_electron_eta_selection->at(i), treeReader.ElecInfo_electron_phi_selection->at(i), treeReader.ElecInfo_electron_energy_selection->at(i));
                Leptons.push_back(electron);
            }
        }
        if(treeReader.num_muons>0){
            for(int i=0; i<treeReader.num_muons; ++i){
                h[hist_MuonInfo_muon_pt] -> Fill(treeReader.MuonInfo_muon_pt_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_MuonInfo_muon_eta] -> Fill(treeReader.MuonInfo_muon_eta_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_MuonInfo_muon_phi] -> Fill(treeReader.MuonInfo_muon_phi_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_MuonInfo_muon_energy] -> Fill(treeReader.MuonInfo_muon_energy_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_MuonInfo_muon_diphoton_deltaR] -> Fill(treeReader.MuonInfo_muon_diphoton_deltaR_selection->at(i), isData ? 1. : NormalizationFactor);
                //------------------------
                h[hist_lepton_pt] -> Fill(treeReader.MuonInfo_muon_pt_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_eta] -> Fill(treeReader.MuonInfo_muon_eta_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_phi] -> Fill(treeReader.MuonInfo_muon_phi_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_energy] -> Fill(treeReader.MuonInfo_muon_energy_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_lepton_diphoton_deltaR] -> Fill(treeReader.MuonInfo_muon_diphoton_deltaR_selection->at(i), isData ? 1. : NormalizationFactor);
                //------------------------
                TLorentzVector muon; muon.SetPtEtaPhiE(treeReader.MuonInfo_muon_pt_selection->at(i), treeReader.MuonInfo_muon_eta_selection->at(i), treeReader.MuonInfo_muon_phi_selection->at(i), treeReader.MuonInfo_muon_energy_selection->at(i));
                Leptons.push_back(muon);
            }
        }
        //}}}
        //### Jets{{{
        //--------- Jets ---------//
        std::vector<TLorentzVector> Jets;
        std::vector<double> Jets_btag_score;
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
                double btag_score = treeReader.JetInfo_jet_pfDeepCSVJetTags_probb_selection->at(i)+treeReader.JetInfo_jet_pfDeepCSVJetTags_probbb_selection->at(i);
                Jets_btag_score.push_back(btag_score);
                if(i==0){//leading jet
                    h[hist_jet1_pt] -> Fill(treeReader.JetInfo_jet_pt_selection->at(i), isData ? 1. : NormalizationFactor);
                    h[hist_jet1_eta] -> Fill(treeReader.JetInfo_jet_eta_selection->at(i), isData ? 1. : NormalizationFactor);
                    h[hist_jet1_phi] -> Fill(treeReader.JetInfo_jet_phi_selection->at(i), isData ? 1. : NormalizationFactor);
                    h[hist_jet1_energy] -> Fill(treeReader.JetInfo_jet_energy_selection->at(i), isData ? 1. : NormalizationFactor);
                    h[hist_jet1_btag_score] -> Fill(btag_score, isData ? 1. : NormalizationFactor);
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
        //================================================//
        //-----------   top reconstruction     -----------//
        //================================================//
        std::vector<int> index_jet_chi2_modified(2);
        std::vector<TLorentzVector> jet_chi2_modified(2);
        //### chi-2 sorting{{{
        //----- chi-2 modified -----//
        double chi2_min = 999;
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
        h[hist_wjet1_pt] -> Fill(jet_chi2_modified[0].Pt(), isData ? 1. : NormalizationFactor);
        h[hist_wjet1_eta] -> Fill(jet_chi2_modified[0].Eta(), isData ? 1. : NormalizationFactor);
        h[hist_wjet1_phi] -> Fill(jet_chi2_modified[0].Phi(), isData ? 1. : NormalizationFactor);
        h[hist_wjet1_energy] -> Fill(jet_chi2_modified[0].E(), isData ? 1. : NormalizationFactor);
        h[hist_wjet1_btag_score] -> Fill(Jets_btag_score[index_jet_chi2_modified[0]], isData ? 1. : NormalizationFactor);
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
        h[hist_wjet2_diphoton_deltaR] -> Fill(jet_chi2_modified[1].DeltaR(diphoton), isData ? 1. : NormalizationFactor);
        if(treeReader.num_leptons>0){
            for(int i=0; i<treeReader.num_leptons; ++i){
                double delta_R = jet_chi2_modified[1].DeltaR(Leptons.at(i));
                h[hist_wjet2_lepton_deltaR] -> Fill(delta_R, isData ? 1. : NormalizationFactor);
            }
        }
        //}}}
        //### Reconstruct Mass (W, M2, M1){{{
        TLorentzVector w_candidate_chi2_modified = jet_chi2_modified[0] + jet_chi2_modified[1];
        TLorentzVector top_candidate_chi2_modified = w_candidate_chi2_modified + bjet;
        h[hist_w_candidate_pt]->Fill(w_candidate_chi2_modified.Pt(), isData ? 1. : NormalizationFactor);
        h[hist_w_candidate_eta]->Fill(w_candidate_chi2_modified.Eta(), isData ? 1. : NormalizationFactor);
        h[hist_w_candidate_mass]->Fill(w_candidate_chi2_modified.M(), isData ? 1. : NormalizationFactor);
        h[hist_top_tbw_pt]->Fill(top_candidate_chi2_modified.Pt(), isData ? 1. : NormalizationFactor);
        h[hist_top_tbw_eta]->Fill(top_candidate_chi2_modified.Eta(), isData ? 1. : NormalizationFactor);
        h[hist_top_tbw_mass]->Fill(top_candidate_chi2_modified.M(), isData ? 1. : NormalizationFactor);
        double M1;
        int index_q;
        TLorentzVector jet_q;
        TLorentzVector top_fcnh = GetBestM1(M1, treeReader.num_jets, index_bjet, index_jet_chi2_modified, diphoton, Jets, index_q, jet_q);
        if(M1 != -999){
            h[hist_top_tqh_pt]->Fill(top_fcnh.Pt(), isData ? 1. : NormalizationFactor);//exclude event without suitable candidate
            h[hist_top_tqh_eta]->Fill(top_fcnh.Eta(), isData ? 1. : NormalizationFactor);//exclude event without suitable candidate
            h[hist_top_tqh_mass]->Fill(top_fcnh.M(), isData ? 1. : NormalizationFactor);//exclude event without suitable candidate
            //-----
            h[hist_jetq_pt] -> Fill(jet_q.Pt(), isData ? 1. : NormalizationFactor);
            h[hist_jetq_eta] -> Fill(jet_q.Eta(), isData ? 1. : NormalizationFactor);
            h[hist_jetq_phi] -> Fill(jet_q.Phi(), isData ? 1. : NormalizationFactor);
            h[hist_jetq_energy] -> Fill(jet_q.E(), isData ? 1. : NormalizationFactor);
            h[hist_jetq_btag_score] -> Fill(Jets_btag_score[index_q], isData ? 1. : NormalizationFactor);
            h[hist_jetq_diphoton_deltaR] -> Fill(jet_q.DeltaR(diphoton), isData ? 1. : NormalizationFactor);
            if(treeReader.num_leptons>0){
                for(int i=0; i<treeReader.num_leptons; ++i){
                    double delta_R = jet_q.DeltaR(Leptons.at(i));
                    h[hist_jetq_lepton_deltaR] -> Fill(delta_R, isData ? 1. : NormalizationFactor);
                }
            }
        }
        //}}}
        //### deltaR{{{
        double deltaR;
        if(M1 != -999) deltaR = top_candidate_chi2_modified.DeltaR(top_fcnh)    ; h[hist_deltaR_top_top] -> Fill(deltaR, isData ? 1. : NormalizationFactor)       ;
        if(M1 != -999) deltaR = diphoton.DeltaR(jet_q)                          ; h[hist_deltaR_qH] -> Fill(deltaR, isData ? 1. : NormalizationFactor)            ;
        deltaR = Jets[0].DeltaR(Jets[1])                                        ; h[hist_deltaR_jet1_jet2] -> Fill(deltaR, isData ? 1. : NormalizationFactor)     ;
        deltaR = jet_chi2_modified[0].DeltaR(jet_chi2_modified[1])              ; h[hist_deltaR_wjet1_wjet2] -> Fill(deltaR, isData ? 1. : NormalizationFactor)     ;
        deltaR = leading_photon.DeltaR(subleading_photon)                       ; h[hist_deltaR_photon_photon] -> Fill(deltaR, isData ? 1. : NormalizationFactor) ;
        deltaR = bjet.DeltaR(w_candidate_chi2_modified)                         ; h[hist_deltaR_bW] -> Fill(deltaR, isData ? 1. : NormalizationFactor)            ;
        deltaR = diphoton.DeltaR(w_candidate_chi2_modified)                     ; h[hist_deltaR_HW] -> Fill(deltaR, isData ? 1. : NormalizationFactor)            ;
        deltaR = diphoton.DeltaR(top_candidate_chi2_modified)                   ; h[hist_deltaR_tH] -> Fill(deltaR, isData ? 1. : NormalizationFactor)            ;
        //}}}
    }//end of event loop

    //### yields, plots, close{{{
    //==================================================//
    //---------------------  Yields  -------------------//
    //==================================================//
    double yields = h[hist_EvtInfo_NPu] -> Integral();
    printf("[INFO] Expected yields = %f\n", yields);
    //==================================================//
    //--------------------- MakePlot -------------------//
    //==================================================//
    TCanvas *c1 = new TCanvas("c1", "c1", 700, 800);
    for(int i=0; i<totalHistNum; ++i){
        MakePlots(c1, h[i], Form("%s/%s.png", output_dir, histNames[i].c_str()));
    }
    //=================================================//
    //---------------------  Close  -------------------//
    //=================================================//
    fout -> Close();
    f_mcpu -> Close();
    fin  -> Close();
    //}}}
}

//### Main function{{{
int main(int argc, char *argv[]){
    char input_file[512]; sprintf(input_file, "%s", argv[1]); printf("[INFO] input_file  = %s\n", input_file);
    char output_file[512]; sprintf(output_file, "%s", argv[2]); printf("[INFO] output_file = %s\n", output_file);
    char dataset[512]; sprintf(dataset, "%s", argv[3]); printf("[INFO] dataset     = %s\n", dataset);
    char output_dir[512];  sprintf(output_dir, "%s", argv[4]); printf("[INFO] output_dir  = %s\n", output_dir);
    char channel[512];  sprintf(channel, "%s", argv[5]); printf("[INFO] channel     = %s\n", channel);
    Selection(input_file, output_file, dataset, output_dir, channel);
    return 1;
}
//}}}
// ### GetBestM1{{{
TLorentzVector GetBestM1(double M1, int num_jets, int index_bjet, std::vector<int> index_jet, TLorentzVector diphoton, std::vector<TLorentzVector> Jets, int &index_q, TLorentzVector &jet_q){
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
    mytree -> SetBranchAddress("ElecInfo_electron_pt", &ElecInfo_electron_pt_selection);
    mytree -> SetBranchAddress("ElecInfo_electron_eta", &ElecInfo_electron_eta_selection);
    mytree -> SetBranchAddress("ElecInfo_electron_phi", &ElecInfo_electron_phi_selection);
    mytree -> SetBranchAddress("ElecInfo_electron_energy", &ElecInfo_electron_energy_selection);
    mytree -> SetBranchAddress("ElecInfo_electron_diphoton_deltaR", &ElecInfo_electron_diphoton_deltaR_selection);
    mytree -> SetBranchAddress("MuonInfo_muon_pt", &MuonInfo_muon_pt_selection);
    mytree -> SetBranchAddress("MuonInfo_muon_eta", &MuonInfo_muon_eta_selection);
    mytree -> SetBranchAddress("MuonInfo_muon_phi", &MuonInfo_muon_phi_selection);
    mytree -> SetBranchAddress("MuonInfo_muon_energy", &MuonInfo_muon_energy_selection);
    mytree -> SetBranchAddress("MuonInfo_muon_diphoton_deltaR", &MuonInfo_muon_diphoton_deltaR_selection);
    //------------------------
    mytree -> SetBranchAddress("GenPartInfo_size", &GenPartInfo_size);
    mytree -> SetBranchAddress("GenPartInfo_gen_Pt", &GenPartInfo_gen_Pt_selection);
    mytree -> SetBranchAddress("GenPartInfo_gen_Eta", &GenPartInfo_gen_Eta_selection);
    mytree -> SetBranchAddress("GenPartInfo_gen_Phi", &GenPartInfo_gen_Phi_selection);
    mytree -> SetBranchAddress("GenPartInfo_gen_Mass", &GenPartInfo_gen_Mass_selection);
    mytree -> SetBranchAddress("GenPartInfo_gen_PdgID", &GenPartInfo_gen_PdgID_selection);
    mytree -> SetBranchAddress("GenPartInfo_gen_Status", &GenPartInfo_gen_Status_selection);
    mytree -> SetBranchAddress("GenPartInfo_gen_nMo", &GenPartInfo_gen_nMo_selection);
    mytree -> SetBranchAddress("GenPartInfo_gen_nDa", &GenPartInfo_gen_nDa_selection);
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
    //------------------------
    mytree -> SetBranchAddress("num_bjets", &num_bjets);// # of selected objects.
    mytree -> SetBranchAddress("JetInfo_leading_bjet_pt", &JetInfo_leading_bjet_pt_selection);
    mytree -> SetBranchAddress("JetInfo_leading_bjet_eta", &JetInfo_leading_bjet_eta_selection);
    mytree -> SetBranchAddress("JetInfo_leading_bjet_phi", &JetInfo_leading_bjet_phi_selection);
    mytree -> SetBranchAddress("JetInfo_leading_bjet_energy", &JetInfo_leading_bjet_energy_selection);
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
    JetInfo_leading_bjet_pt_selection = new std::vector<float>;
    JetInfo_leading_bjet_eta_selection = new std::vector<float>;
    JetInfo_leading_bjet_phi_selection = new std::vector<float>;
    JetInfo_leading_bjet_energy_selection = new std::vector<float>;
    ElecInfo_electron_pt_selection = new std::vector<float>;
    ElecInfo_electron_eta_selection = new std::vector<float>;
    ElecInfo_electron_phi_selection = new std::vector<float>;
    ElecInfo_electron_energy_selection = new std::vector<float>;
    ElecInfo_electron_diphoton_deltaR_selection = new std::vector<float>;
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
    delete JetInfo_leading_bjet_pt_selection;
    delete JetInfo_leading_bjet_eta_selection;
    delete JetInfo_leading_bjet_phi_selection;
    delete JetInfo_leading_bjet_energy_selection;
    delete ElecInfo_electron_pt_selection;
    delete ElecInfo_electron_eta_selection;
    delete ElecInfo_electron_phi_selection;
    delete ElecInfo_electron_energy_selection;
    delete ElecInfo_electron_diphoton_deltaR_selection;
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
