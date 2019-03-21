//===========================================//
//=== 27. Feb. 2019                       ===//
//=== Develop for single top tqH analysis ===//
//===========================================//
#include <stdio.h>
#include <math.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <vector>
#include <string>
#include "../include/main.h"
#include "../include/cross_section.h"
using namespace std;
//==========================//
//---  Useful Constants  ---//
//==========================//
double pfDeepCSVJetTags_tight  = 0.8001;
double pfDeepCSVJetTags_medium = 0.4941;
double pfDeepCSVJetTags_loose  = 0.1522;
double w_boson_mass = 80.385;//GeV
double top_quark_mass = 173.0 ;//GeV


int main(int argc, char *argv[]){
    char input_file[256]; sprintf(input_file, "%s", argv[1]); printf("[INFO] input_file  = %s\n", input_file);
    char output_file[256]; sprintf(output_file, "%s", argv[2]); printf("[INFO] output_file = %s\n", output_file);
    char dataset[256]; sprintf(dataset, "%s", argv[3]); printf("[INFO] dataset     = %s\n", dataset);
    bool isData = isThisDataOrNot(dataset);
    bool isMCsignal = isThisMCsignal(dataset);
    bool isMultiFile = isThisMultiFile(dataset);
    TFile *fout = new TFile(output_file, "RECREATE");
    TChain *flashggStdTree = new TChain("flashggNtuples/flashggStdTree");
    if(isMultiFile) flashggStdTree->Add(Form("%s/*.root", input_file));
    else flashggStdTree->Add(input_file);

    //==================//
    //--- histograms ---//
    //==================//
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    TH1D  *hist_num_jets = new TH1D("hist_num_jets", "hist_num_jets", 20, 0, 20);
    TH1D  *hist_num_btagged_jets = new TH1D("hist_num_btagged_jets", "hist_num_btagged_jets", 10, 0, 10);
    TH1D  *hist_num_nonbtagged_jets = new TH1D("hist_num_nonbtagged_jets", "hist_num_nonbtagged_jets", 20, 0, 20);
    //------------------------
    TH1D  *hist_bjet_pt = new TH1D("hist_bjet_pt", "hist_bjet_pt", 50, 0, 1000);
    TH1D  *hist_jet1_pt = new TH1D("hist_jet1_pt", "hist_jet1_pt", 50, 0, 1000);
    TH1D  *hist_jet2_pt = new TH1D("hist_jet2_pt", "hist_jet2_pt", 50, 0, 1000);
    TH1D  *hist_bjet_eta = new TH1D("hist_bjet_eta", "hist_bjet_eta", 40, -2.5, 2.5);
    TH1D  *hist_jet1_eta = new TH1D("hist_jet1_eta", "hist_jet1_eta", 40, -2.5, 2.5);
    TH1D  *hist_jet2_eta = new TH1D("hist_jet2_eta", "hist_jet2_eta", 40, -2.5, 2.5);
    TH1D  *hist_bjet_phi = new TH1D("hist_bjet_phi", "hist_bjet_phi", 40, -3.0, 3.0);
    TH1D  *hist_jet1_phi = new TH1D("hist_jet1_phi", "hist_jet1_phi", 40, -3.0, 3.0);
    TH1D  *hist_jet2_phi = new TH1D("hist_jet2_phi", "hist_jet2_phi", 40, -3.0, 3.0);
    //------------------------
    TH1D  *hist_dijet_eta = new TH1D("hist_dijet_eta", "hist_dijet_eta", 40, 0., 5.0);
    TH1D  *hist_dijet_phi = new TH1D("hist_dijet_phi", "hist_dijet_phi", 40, 0., 6.0);
    TH1D  *hist_dijet_angle = new TH1D("hist_dijet_angle", "hist_dijet_angle", 40, 0., 8.0);
    TH1D  *hist_inv_mass_dijet = new TH1D("hist_inv_mass_dijet", "hist_inv_mass_dijet", 50, 55, 105);
    TH1D  *hist_inv_mass_diphoton = new TH1D("hist_inv_mass_diphoton", "hist_inv_mass_diphoton", 50, 100, 150);
    TH1D  *hist_inv_mass_tbw = new TH1D("hist_inv_mass_tbw", "hist_inv_mass_tbw", 50, 0, 500);
    //------------------------
    TH1D  *hist_DiPhoInfo_leadPt = new TH1D("hist_DiPhoInfo_leadPt", "DiPhoInfo_leadPt", 50, 0., 1000);
    TH1D  *hist_DiPhoInfo_leadEta = new TH1D("hist_DiPhoInfo_leadEta", "DiPhoInfo_leadEta", 40, -2.5, 2.5);
    TH1D  *hist_DiPhoInfo_leadPhi = new TH1D("hist_DiPhoInfo_leadPhi", "DiPhoInfo_leadPhi", 40, -3., 3.);
    TH1D  *hist_DiPhoInfo_leadE = new TH1D("hist_DiPhoInfo_leadE", "DiPhoInfo_leadE", 50, 0., 1000);
    TH1D  *hist_DiPhoInfo_leadIDMVA = new TH1D("hist_DiPhoInfo_leadIDMVA", "DiPhoInfo_leadIDMVA", 50, -1., 1.);
    TH1D  *hist_DiPhoInfo_subleadPt = new TH1D("hist_DiPhoInfo_subleadPt", "DiPhoInfo_subleadPt", 50, 0., 1000);
    TH1D  *hist_DiPhoInfo_subleadEta = new TH1D("hist_DiPhoInfo_subleadEta", "DiPhoInfo_subleadEta", 40, -2.5, 2.5);
    TH1D  *hist_DiPhoInfo_subleadPhi = new TH1D("hist_DiPhoInfo_subleadPhi", "DiPhoInfo_subleadPhi", 40, -3., 3.);
    TH1D  *hist_DiPhoInfo_subleadE = new TH1D("hist_DiPhoInfo_subleadE", "DiPhoInfo_subleadE", 50, 0., 1000);
    TH1D  *hist_DiPhoInfo_subleadIDMVA = new TH1D("hist_DiPhoInfo_subleadIDMVA", "DiPhoInfo_subleadIDMVA", 50, -1., 1.);
    //------------------------
    TH1D  *hist_DiPhoInfo_leadIDMVA_ori = new TH1D("hist_DiPhoInfo_leadIDMVA_ori", "DiPhoInfo_leadIDMVA", 50, -1., 1.);
    TH1D  *hist_DiPhoInfo_subleadIDMVA_ori = new TH1D("hist_DiPhoInfo_subleadIDMVA_ori", "DiPhoInfo_subleadIDMVA", 50, -1., 1.);
    TH1D  *hist_inv_mass_diphoton_ori = new TH1D("hist_inv_mass_diphoton_ori", "hist_inv_mass_diphoton", 50, 100, 150);
    //------------------------
    hist_num_jets -> Sumw2();
    hist_num_btagged_jets -> Sumw2();
    hist_num_nonbtagged_jets -> Sumw2();
    hist_bjet_pt -> Sumw2();
    hist_jet1_pt -> Sumw2();
    hist_jet2_pt -> Sumw2();
    hist_bjet_eta -> Sumw2();
    hist_jet1_eta -> Sumw2();
    hist_jet2_eta -> Sumw2();
    hist_bjet_phi -> Sumw2();
    hist_jet1_phi -> Sumw2();
    hist_jet2_phi -> Sumw2();
    hist_dijet_eta -> Sumw2();
    hist_dijet_phi -> Sumw2();
    hist_dijet_angle -> Sumw2();
    hist_inv_mass_dijet -> Sumw2();
    hist_inv_mass_diphoton -> Sumw2();
    hist_inv_mass_tbw -> Sumw2();
    hist_DiPhoInfo_leadPt -> Sumw2();
    hist_DiPhoInfo_leadEta -> Sumw2();
    hist_DiPhoInfo_leadPhi -> Sumw2();
    hist_DiPhoInfo_leadE -> Sumw2();
    hist_DiPhoInfo_leadIDMVA -> Sumw2();
    hist_DiPhoInfo_subleadPt -> Sumw2();
    hist_DiPhoInfo_subleadEta -> Sumw2();
    hist_DiPhoInfo_subleadPhi -> Sumw2();
    hist_DiPhoInfo_subleadE -> Sumw2();
    hist_DiPhoInfo_subleadIDMVA -> Sumw2();
    hist_DiPhoInfo_leadIDMVA_ori -> Sumw2();
    hist_DiPhoInfo_subleadIDMVA_ori -> Sumw2();
    hist_inv_mass_diphoton_ori -> Sumw2();

    //==========================//
    //--- readout parameters ---//
    //==========================//
    Int_t jets_size=0;
    std::vector<float> *JetInfo_Pt=0;
    std::vector<float> *JetInfo_Eta=0;
    std::vector<float> *JetInfo_Phi=0;
    std::vector<float> *JetInfo_Mass=0;
    std::vector<float> *JetInfo_Energy=0;
    std::vector<float> *JetInfo_pfDeepCSVJetTags_probb=0;
    std::vector<float> *JetInfo_pfDeepCSVJetTags_probbb=0;
    float EvtInfo_genweight=0;
    float DiPhoInfo_mass=0;
    float DiPhoInfo_leadPt=0;
    float DiPhoInfo_leadEta=0;
    float DiPhoInfo_leadPhi=0;
    float DiPhoInfo_leadE=0;
    float DiPhoInfo_leadIDMVA=0;
    float DiPhoInfo_subleadPt=0;
    float DiPhoInfo_subleadEta=0;
    float DiPhoInfo_subleadPhi=0;
    float DiPhoInfo_subleadE=0;
    float DiPhoInfo_subleadIDMVA=0;
    //------------------------
    flashggStdTree->SetBranchAddress("EvtInfo.genweight", &EvtInfo_genweight);
    flashggStdTree->SetBranchAddress("jets_size", &jets_size);
    flashggStdTree->SetBranchAddress("JetInfo.Pt", &JetInfo_Pt);
    flashggStdTree->SetBranchAddress("JetInfo.Eta", &JetInfo_Eta);
    flashggStdTree->SetBranchAddress("JetInfo.Phi", &JetInfo_Phi);
    flashggStdTree->SetBranchAddress("JetInfo.Mass", &JetInfo_Mass);
    flashggStdTree->SetBranchAddress("JetInfo.Energy", &JetInfo_Energy);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probb", &JetInfo_pfDeepCSVJetTags_probb);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probbb", &JetInfo_pfDeepCSVJetTags_probbb);
    flashggStdTree->SetBranchAddress("DiPhoInfo.mass", &DiPhoInfo_mass);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadPt", &DiPhoInfo_leadPt);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadEta", &DiPhoInfo_leadEta);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadPhi", &DiPhoInfo_leadPhi);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadE", &DiPhoInfo_leadE);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadIDMVA", &DiPhoInfo_leadIDMVA);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadPt", &DiPhoInfo_subleadPt);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadEta", &DiPhoInfo_subleadEta);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadPhi", &DiPhoInfo_subleadPhi);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadE", &DiPhoInfo_subleadE);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadIDMVA", &DiPhoInfo_subleadIDMVA);

    std::vector<int> vec_bjet_indices;
    std::vector<TLorentzVector> vec_btagged_jets;
    std::vector<TLorentzVector> vec_nonbtagged_jets;

    //##################################################//
    //########   Event Loop [Normalization]   ##########//
    //##################################################//
    int nentries = flashggStdTree->GetEntries(); printf("[INFO] N_entries = %d\n", nentries);
    double NormalizationFactor;
    double Luminosity = 42.; //fb{-1}
    double CrossSection = GetXsec(dataset); //pb
    double BranchingFraction = GetBranchingFraction(dataset); //pb
    printf("[INFO] CrossSection = %f !\n", CrossSection);
    printf("[INFO] Equivalent lumi. = %f !\n", (double)nentries/CrossSection);
    printf("[INFO] BranchingFraction = %f !\n", BranchingFraction);
    double TotalGenweight=0;
    for(int ientry=0; ientry<nentries; ientry++){
        flashggStdTree->GetEntry(ientry);
        TotalGenweight+=EvtInfo_genweight;
    }
    NormalizationFactor = 1000. * Luminosity * CrossSection * BranchingFraction / TotalGenweight;
    printf("\n[INFO] TotalGenweight = %f!\n", TotalGenweight);
    printf("[INFO] NormalizationFactor = %f!\n", isData ? 1. : NormalizationFactor);
    //##################################################//
    //#########    Event Loop [Selection]    ###########//
    //##################################################//
    // Goal 1: t->b+W(jj), 1 bjet + 2 chi2 jets 
    // Goal 2: diphoton info
    int nevents_pass_selection = 0;
    int Nevents_CUT_num_bjets_is_not_exactly_one = 0;
    int Nevents_anti_CUT_num_bjets_is_not_exactly_one = 0;
    int Nevents_CUT_no_bjet_events = 0;
    int Nevents_anti_CUT_no_bjet_events = 0;
    int Nevents_CUT_num_nonbjets_less_than_2 = 0;
    int Nevents_anti_CUT_num_nonbjets_less_than_2 = 0;
    int Nevents_CUT_nand = 0;
    int Nevents_anti_CUT_nand = 0;
    for(int ientry=0; ientry<nentries; ientry++){
        flashggStdTree->GetEntry(ientry);//load data
        if((ientry+1)%1000==0 || (ientry+1)==nentries) printf("ientry = %d\r", ientry);
        //==================================================//
        //-------------   Reset Parameters   ---------------//
        //==================================================//
        vec_nonbtagged_jets.clear();
        //==================================================//
        //--------------   Basic Selectoin   ---------------//
        //==================================================//
        if(DiPhoInfo_mass<0) continue;
        if( !(jets_size>0) ) continue;
        if(DiPhoInfo_mass<100 || DiPhoInfo_mass>150) continue;
        if( !isMCsignal && DiPhoInfo_mass>120 && DiPhoInfo_mass<130) continue;
        //--------------------
        hist_DiPhoInfo_leadIDMVA_ori -> Fill(DiPhoInfo_leadIDMVA, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_DiPhoInfo_subleadIDMVA_ori -> Fill(DiPhoInfo_subleadIDMVA, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_inv_mass_diphoton_ori -> Fill(DiPhoInfo_mass, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        //--------------------
        if(!(DiPhoInfo_leadIDMVA>0)) continue;
        //==================================================//
        //------------   Physical Observables   ------------//
        //==================================================//
        
        //==================================================//
        //-----  Reconstruction(tbW): Select one bjet  -----//
        //==================================================//
        int bjetindex=-1, num_bjets=0;
        TLorentzVector leading_bjet;
        for(int i=0; i<jets_size; i++){
            if( fabs(JetInfo_Eta->at(i)) > 2.4 ) continue;
            if( fabs(JetInfo_Pt->at(i))  < 30  ) continue;
            if(JetInfo_pfDeepCSVJetTags_probb->at(i)+JetInfo_pfDeepCSVJetTags_probbb->at(i) >= pfDeepCSVJetTags_tight){
                if(bjetindex==-1) leading_bjet.SetPtEtaPhiE(JetInfo_Pt->at(i), JetInfo_Eta->at(i), JetInfo_Phi->at(i), JetInfo_Energy->at(i));
                bjetindex = i;
                TLorentzVector jet_lorentzvector;
                jet_lorentzvector.SetPtEtaPhiE(JetInfo_Pt->at(i), JetInfo_Eta->at(i), JetInfo_Phi->at(i), JetInfo_Energy->at(i));
                vec_btagged_jets.push_back(jet_lorentzvector);
                vec_bjet_indices.push_back(bjetindex);
                num_bjets+=1;
            }
        }
        hist_num_btagged_jets->Fill(num_bjets, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        bool CUT_no_bjet_events = (bjetindex==-1) ? true : false;
        bool CUT_num_bjets_is_not_exactly_one = (num_bjets!=1) ? true : false;
        //if(bjetindex==-1) continue;// discard no-b-jet events
        //if(num_bjets!=1)   continue;// exactly 1 b-jet criterion
        hist_bjet_pt->Fill( CUT_no_bjet_events ? -999. : leading_bjet.Pt(), isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_bjet_eta->Fill( CUT_no_bjet_events ? -999. : leading_bjet.Eta(), isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_bjet_phi->Fill( CUT_no_bjet_events ? -999. : leading_bjet.Phi(), isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        //==================================================//
        //-----  Reconstruction(tbW): Store rest jets  -----//
        //==================================================//
        bool isbjet;
        int num_nonbtagged_jets=0;
        for(int i=0; i<jets_size; i++){
            if( fabs(JetInfo_Eta->at(i)) > 2.4 ) continue;
            if( fabs(JetInfo_Pt->at(i))  < 30  ) continue;
            // Exclude bjet
            //if(i == bjetindex) continue;
            isbjet=false;
            if(!CUT_no_bjet_events){for(int j=0; j<vec_btagged_jets.size(); j++){if(i==vec_bjet_indices[j]) isbjet=true;}}
            if(isbjet) continue;
            //--------------------
            TLorentzVector jet_lorentzvector;
            jet_lorentzvector.SetPtEtaPhiE(JetInfo_Pt->at(i), JetInfo_Eta->at(i), JetInfo_Phi->at(i), JetInfo_Energy->at(i));
            vec_nonbtagged_jets.push_back(jet_lorentzvector);
            num_nonbtagged_jets+=1;
        }
        bool CUT_num_nonbjets_less_than_2 = (num_nonbtagged_jets<2) ? true : false;
        //if(num_nonbtagged_jets<2) continue;// at least 2 non-b-tagged jets for construction of w boson
        hist_num_nonbtagged_jets->Fill(num_nonbtagged_jets, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_num_jets->Fill(num_bjets+num_nonbtagged_jets, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        //==================================================//
        //---------------------  Test  ---------------------//
        //==================================================//
        if(CUT_num_bjets_is_not_exactly_one) Nevents_CUT_num_bjets_is_not_exactly_one += 1;
        else Nevents_anti_CUT_num_bjets_is_not_exactly_one += 1;
        if(CUT_no_bjet_events) Nevents_CUT_no_bjet_events += 1;
        else Nevents_anti_CUT_no_bjet_events += 1;
        if(CUT_num_nonbjets_less_than_2) Nevents_CUT_num_nonbjets_less_than_2 += 1;
        else Nevents_anti_CUT_num_nonbjets_less_than_2 += 1;
        if(!CUT_no_bjet_events && !CUT_num_nonbjets_less_than_2) Nevents_CUT_nand += 1;
        else Nevents_anti_CUT_nand += 1;


        if(bjetindex==-1) continue;
        if(num_bjets!=1)   continue;
        if(num_nonbtagged_jets<2) continue;

        int wjetindices[2]={0};//indices in non-btagged vector
        double dijet_invariant_mass, chi2, chi2_min=9999;
        TLorentzVector best_dijet_pair;
        for(int i=0; i<vec_nonbtagged_jets.size(); i++){
            for(int j=i+1; j<vec_nonbtagged_jets.size(); j++){
                TLorentzVector w_candidate_temp = vec_nonbtagged_jets[i] + vec_nonbtagged_jets[j];
                dijet_invariant_mass = w_candidate_temp.M();
                chi2 = (dijet_invariant_mass-w_boson_mass)*(dijet_invariant_mass-w_boson_mass);
                if(chi2<chi2_min){wjetindices[0]=i; wjetindices[1]=j; chi2_min=chi2;};
            }
        }
        best_dijet_pair = vec_nonbtagged_jets[wjetindices[0]] + vec_nonbtagged_jets[wjetindices[1]];
        dijet_invariant_mass = best_dijet_pair.M();

        TLorentzVector tbw_bjj = best_dijet_pair + leading_bjet;
        double inv_mass_tbw = tbw_bjj.M();
        
        hist_inv_mass_dijet->Fill(dijet_invariant_mass, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_inv_mass_tbw->Fill(inv_mass_tbw, isData ? 1. : NormalizationFactor*EvtInfo_genweight);


        /*
        //==================================================//
        //------  Reconstruction(tbW): Chi-2 sorting  ------//
        //==================================================//
        TLorentzVector best_dijet, best_trijet;//trijet=bjj
        int good_bjet_index=-1, wjetindices[2]={0};//particles indices for min chi2 
        double dijet_invariant_mass, trijet_invariant_mass, chi2, chi2_min=99999;
        if(!CUT_num_nonbjets_less_than_2){// Requirement for at least 1 w-boson!
            for(int i=0; i<vec_nonbtagged_jets.size(); i++){
                for(int j=i+1; j<vec_nonbtagged_jets.size(); j++){
                    TLorentzVector w_candidate = vec_nonbtagged_jets[i] + vec_nonbtagged_jets[j];
                    dijet_invariant_mass = w_candidate.M();
                    if(!CUT_no_bjet_events){// Requirement for at least 1 top!
                        for(int k=0; k<vec_btagged_jets.size(); k++){
                            TLorentzVector top_candidate = w_candidate + vec_btagged_jets[k];
                            trijet_invariant_mass = top_candidate.M();
                            chi2 = Chi2_calculator(dijet_invariant_mass, trijet_invariant_mass);
                            if(chi2<chi2_min){wjetindices[0]=i; wjetindices[1]=j; good_bjet_index=k; chi2_min=chi2;};
                        }// end of k loop (b-jet))
                    }
                }// end of jet2
            }// end of jet1
        }
        bool CUT_no_any_top_candidate = CUT_num_nonbjets_less_than_2 || CUT_no_bjet_events;
        //==================================================//
        //-----  Reconstruction(tbW): Store dijet info  ----//
        //==================================================//
        //If CUT is true, to prevent invalid addition of elements in vec, best_dijet is set to certain value.
        best_dijet = CUT_num_nonbjets_less_than_2 ? best_dijet : vec_nonbtagged_jets[wjetindices[0]] + vec_nonbtagged_jets[wjetindices[1]];
        dijet_invariant_mass = CUT_num_nonbjets_less_than_2 ? -999. : best_dijet.M();
        //------------------------------
        double dijet_jet1_pt = CUT_num_nonbjets_less_than_2 ? -999. : vec_nonbtagged_jets[wjetindices[0]].Pt();
        double dijet_jet1_eta = CUT_num_nonbjets_less_than_2 ? -999. : vec_nonbtagged_jets[wjetindices[0]].Eta();
        double dijet_jet1_phi = CUT_num_nonbjets_less_than_2 ? -999. : vec_nonbtagged_jets[wjetindices[0]].Phi();
        double dijet_jet2_pt = CUT_num_nonbjets_less_than_2 ? -999. : vec_nonbtagged_jets[wjetindices[1]].Pt();
        double dijet_jet2_eta = CUT_num_nonbjets_less_than_2 ? -999. : vec_nonbtagged_jets[wjetindices[1]].Eta();
        double dijet_jet2_phi = CUT_num_nonbjets_less_than_2 ? -999. : vec_nonbtagged_jets[wjetindices[1]].Phi();
        //------------------------------
        double dijet_eta_difference = CUT_num_nonbjets_less_than_2 ? -999. : fabs(dijet_jet1_eta - dijet_jet2_eta);
        double dijet_phi_difference = CUT_num_nonbjets_less_than_2 ? -999. : fabs(dijet_jet1_phi - dijet_jet2_phi);
        double dijet_angle_difference = sqrt(dijet_eta_difference*dijet_eta_difference + dijet_phi_difference*dijet_phi_difference);
        //------------------------------
        hist_jet1_pt->Fill(dijet_jet1_pt, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_jet1_eta->Fill(dijet_jet1_eta, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_jet1_phi->Fill(dijet_jet1_phi, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_jet2_pt->Fill(dijet_jet2_pt, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_jet2_eta->Fill(dijet_jet2_eta, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_jet2_phi->Fill(dijet_jet2_phi, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_dijet_eta->Fill(dijet_eta_difference, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_dijet_phi->Fill(dijet_phi_difference, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_dijet_angle->Fill(dijet_angle_difference, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_inv_mass_dijet->Fill(dijet_invariant_mass, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        //==================================================//
        //-----  Reconstruction(tbW): InvMass of top  ------//
        //==================================================//
        //TLorentzVector tbw_bjj = best_dijet + leading_bjet;
        //double inv_mass_tbw = CUT_no_any_top_candidate ? -999. : tbw_bjj.M();
        best_trijet = CUT_no_any_top_candidate ? best_trijet : best_dijet + vec_btagged_jets[good_bjet_index];
        trijet_invariant_mass = CUT_no_any_top_candidate ? -999. : best_trijet.M();
        hist_inv_mass_tbw->Fill(trijet_invariant_mass, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        //==================================================//
        //-----------   Diphoton Related Info    -----------//
        //==================================================//
        hist_DiPhoInfo_leadPt->Fill(DiPhoInfo_leadPt, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_DiPhoInfo_leadEta->Fill(DiPhoInfo_leadEta, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_DiPhoInfo_leadPhi->Fill(DiPhoInfo_leadPhi, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_DiPhoInfo_leadE->Fill(DiPhoInfo_leadE, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_DiPhoInfo_leadIDMVA->Fill(DiPhoInfo_leadIDMVA, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_DiPhoInfo_subleadPt->Fill(DiPhoInfo_subleadPt, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_DiPhoInfo_subleadEta->Fill(DiPhoInfo_subleadEta, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_DiPhoInfo_subleadPhi->Fill(DiPhoInfo_subleadPhi, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_DiPhoInfo_subleadE->Fill(DiPhoInfo_subleadE, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_DiPhoInfo_subleadIDMVA->Fill(DiPhoInfo_subleadIDMVA, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_inv_mass_diphoton->Fill(DiPhoInfo_mass, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        */
        //==================================================//
        //-------------   Event Selection    ---------------//
        //==================================================//
        nevents_pass_selection += 1;
    }// End of event loop.

    //==================================================//
    //---------------------  Test  ---------------------//
    //==================================================//
    printf("[INFO] N_nevents_pass_selection = %d\n", nevents_pass_selection);
    printf("[DEBUG] Nevents_CUT_num_bjets_is_not_exactly_one = %d (%.4f)\n", Nevents_CUT_num_bjets_is_not_exactly_one, (double)Nevents_CUT_num_bjets_is_not_exactly_one/(double)nevents_pass_selection);
    printf("[DEBUG] Nevents_anti_CUT_num_bjets_is_not_exactly_one = %d (%.4f)\n", Nevents_anti_CUT_num_bjets_is_not_exactly_one, (double)Nevents_anti_CUT_num_bjets_is_not_exactly_one/(double)nevents_pass_selection);
    printf("[DEBUG] Nevents_CUT_no_bjet_events = %d (%.4f)\n", Nevents_CUT_no_bjet_events, (double)Nevents_CUT_no_bjet_events/(double)nevents_pass_selection);
    printf("[DEBUG] Nevents_anti_CUT_no_bjet_events = %d (%.4f)\n", Nevents_anti_CUT_no_bjet_events, (double)Nevents_anti_CUT_no_bjet_events/(double)nevents_pass_selection);
    printf("[DEBUG] Nevents_CUT_num_nonbjets_less_than_2 = %d (%.4f)\n", Nevents_CUT_num_nonbjets_less_than_2, (double)Nevents_CUT_num_nonbjets_less_than_2/(double)nevents_pass_selection);
    printf("[DEBUG] Nevents_anti_CUT_num_nonbjets_less_than_2 = %d (%.4f)\n", Nevents_anti_CUT_num_nonbjets_less_than_2, (double)Nevents_anti_CUT_num_nonbjets_less_than_2/(double)nevents_pass_selection);
    printf("[DEBUG] Nevents_CUT_nand = %d (%.4f)\n", Nevents_CUT_nand, (double)Nevents_CUT_nand/(double)nevents_pass_selection);
    printf("[DEBUG] Nevents_anti_CUT_nand = %d (%.4f)\n", Nevents_anti_CUT_nand, (double)Nevents_anti_CUT_nand/(double)nevents_pass_selection);

    char output_dir[256]; sprintf(output_dir, "%s", argv[4]); printf("[INFO] output_dir = %s\n", output_dir);
    MakePlots(c1, hist_inv_mass_tbw, "M2 mass spectrum", Form("%s/hist_inv_mass_tbw.png", output_dir));
    MakePlots(c1, hist_inv_mass_dijet, "W boson mass spectrum", Form("%s/hist_inv_mass_dijet.png", output_dir));


    //##################################################//
    //##############     Make Plots !!    ##############//
    //##################################################//
    /*
    printf("[INFO] Start making plots!\n");
    char output_dir[256]; sprintf(output_dir, "%s", argv[4]); printf("[INFO] output_dir = %s\n", output_dir);
    //------------------------------
    MakePlots(c1, hist_num_jets, "Multiplicity of jets", Form("%s/hist_num_jets.png", output_dir));
    MakePlots(c1, hist_num_btagged_jets, "Multiplicity of b-tagged jets", Form("%s/hist_num_btagged_jets.png", output_dir));
    MakePlots(c1, hist_num_nonbtagged_jets, "Multiplicity of non-b-tagged jets", Form("%s/hist_num_nonbtagged_jets.png", output_dir));
    //------------------------------
    MakePlots(c1, hist_bjet_pt, "Pt of b-tagged jet", Form("%s/hist_bjet_pt.png", output_dir));
    MakePlots(c1, hist_bjet_eta, "Eta of b-tagged jet", Form("%s/hist_bjet_eta.png", output_dir));
    MakePlots(c1, hist_bjet_phi, "Phi of b-tagged jet", Form("%s/hist_bjet_phi.png", output_dir));
    MakePlots(c1, hist_jet1_pt, "Pt of hadronic jet1 candidate", Form("%s/hist_jet1_pt.png", output_dir));
    MakePlots(c1, hist_jet1_eta, "Eta of hadronic jet1 candidate", Form("%s/hist_jet1_eta.png", output_dir));
    MakePlots(c1, hist_jet1_phi, "Phi of hadronic jet1 candidate", Form("%s/hist_jet1_phi.png", output_dir));
    MakePlots(c1, hist_jet2_pt, "Pt of hadronic jet2 candidate", Form("%s/hist_jet2_pt.png", output_dir));
    MakePlots(c1, hist_jet2_eta, "Eta of hadronic jet2 candidate", Form("%s/hist_jet2_eta.png", output_dir));
    MakePlots(c1, hist_jet2_phi, "Phi of hadronic jet2 candidate", Form("%s/hist_jet2_phi.png", output_dir));
    //------------------------------
    MakePlots(c1, hist_dijet_eta, "Eta difference of hadronic dijet candidate", Form("%s/hist_dijet_eta.png", output_dir));
    MakePlots(c1, hist_dijet_phi, "Phi difference of hadronic dijet candidate", Form("%s/hist_dijet_phi.png", output_dir));
    MakePlots(c1, hist_dijet_angle, "Angle difference of hadronic dijet candidate", Form("%s/hist_dijet_angle.png", output_dir));
    //------------------------------
    MakePlots(c1, hist_DiPhoInfo_leadPt, "Diphoton_leadPt", Form("%s/hist_DiPhoInfo_leadPt.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_leadEta, "Diphoton_leadEta", Form("%s/hist_DiPhoInfo_leadEta.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_leadPhi, "Diphoton_leadPhi", Form("%s/hist_DiPhoInfo_leadPhi.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_leadE, "Diphoton_leadE", Form("%s/hist_DiPhoInfo_leadE.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_leadIDMVA, "Diphoton_leadIDMVA", Form("%s/hist_DiPhoInfo_leadIDMVA.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_subleadPt, "Diphoton_subleadPt", Form("%s/hist_DiPhoInfo_subleadPt.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_subleadEta, "Diphoton_subleadEta", Form("%s/hist_DiPhoInfo_subleadEta.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_subleadPhi, "Diphoton_subleadPhi", Form("%s/hist_DiPhoInfo_subleadPhi.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_subleadE, "Diphoton_subleadE", Form("%s/hist_DiPhoInfo_subleadE.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_subleadIDMVA, "Diphoton_subleadIDMVA", Form("%s/hist_DiPhoInfo_subleadIDMVA.png", output_dir));
    MakePlots(c1, hist_inv_mass_tbw, "M2 mass spectrum", Form("%s/hist_inv_mass_tbw.png", output_dir));
    MakePlots(c1, hist_inv_mass_dijet, "W boson mass spectrum", Form("%s/hist_inv_mass_dijet.png", output_dir));
    MakePlots(c1, hist_inv_mass_diphoton, "Diphoton mass spectrum", Form("%s/hist_inv_mass_diphoton.png", output_dir));
    //------------------------------
    MakePlots(c1, hist_DiPhoInfo_leadIDMVA_ori, "Diphoton_leadIDMVA", Form("%s/hist_DiPhoInfo_leadIDMVA_ori.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_subleadIDMVA_ori, "Diphoton_subleadIDMVA", Form("%s/hist_DiPhoInfo_subleadIDMVA_ori.png", output_dir));
    MakePlots(c1, hist_inv_mass_diphoton_ori, "Diphoton mass spectrum", Form("%s/hist_inv_mass_diphoton_ori.png", output_dir));
    //------------------------------
    */
    fout->Close();
    return 1;
}

double Chi2_calculator(double w_mass, double t_mass){
    return (w_mass-w_boson_mass)*(w_mass-w_boson_mass) + (t_mass-top_quark_mass)*(t_mass-top_quark_mass);
}
void MakePlots(TCanvas *c1, TH1D* hist, const char* title, const char* outputFile){
    hist->Draw();
    hist->SetTitle(title);
    hist->SetXTitle(title);
    hist->SetYTitle("Entries");
    hist->GetYaxis()->SetTitleOffset(1.4);
    hist->Write();
    c1->SaveAs(outputFile);
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
bool isThisDataOrNot(char* dataset){
    if((string)dataset == "DoubleEG_B") return true;
    if((string)dataset == "DoubleEG_C") return true;
    if((string)dataset == "DoubleEG_D") return true;
    if((string)dataset == "DoubleEG_E") return true;
    if((string)dataset == "DoubleEG_F") return true;
    return false;
}

bool isThisMultiFile(char* dataset){
    if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return true;
    if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return true;
    return false;
}
