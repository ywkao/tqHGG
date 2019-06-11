//===========================================//
//=== 24. April. 2019                     ===//
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
//2016, https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
double pfCombinedInclusiveSecondaryVertexV2BJetTags_tight  = 0.9535;
double pfCombinedInclusiveSecondaryVertexV2BJetTags_medium = 0.8484;
double pfCombinedInclusiveSecondaryVertexV2BJetTags_loose  = 0.5426;
//2017, https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
double pfDeepCSVJetTags_tight  = 0.8001;
double pfDeepCSVJetTags_medium = 0.4941;
double pfDeepCSVJetTags_loose  = 0.1522;
double w_boson_mass = 80.379;//GeV
double w_boson_width = 9.5;//GeV
double top_quark_mass = 173.0;//GeV
double top_quark_width = 16.3;//GeV


int main(int argc, char *argv[]){
    //============================//
    //----- Input file names -----//
    //============================//
    char input_file[512]; sprintf(input_file, "%s", argv[1]); printf("[INFO] input_file  = %s\n", input_file);
    char output_file[512]; sprintf(output_file, "%s", argv[2]); printf("[INFO] output_file = %s\n", output_file);
    char dataset[512]; sprintf(dataset, "%s", argv[3]); printf("[INFO] dataset     = %s\n", dataset);
    bool isData = isThisDataOrNot(dataset);//Determine normalization factor
    bool isMCsignal = isThisMCsignal(dataset);//Determine to show up in sig region
    TFile *fout = new TFile(output_file, "RECREATE");

    //==============================//
    //----- Read input file(s) -----//
    //==============================//
    flashggStdTreeReader treeReader;
    treeReader.InitChain("flashggNtuples/flashggStdTree");
    treeReader.AddMultiRootFile(input_file);
    treeReader.SetBranchAddresses();

    //===============================//
    //----- Prepare output file -----//
    //===============================//
    myTreeClass mytree;
    mytree.InitTree();
    mytree.MakeNewBranchAddresses();


    // Note: Normalization factors are not used during preselection.
    // It is stored and would be used during selection.
    // Question: keep total genweight before or after preselection?
    //##################################################//
    //########   Event Loop [Normalization]   ##########//
    //##################################################//
    int nentries = treeReader.GetEntries(); printf("[INFO] N_entries = %d\n", nentries);
    double NormalizationFactor;
    double Luminosity = 35.9; //fb{-1}
    double CrossSection = GetXsec(dataset); //pb
    double BranchingFraction = GetBranchingFraction(dataset); //pb
    printf("[INFO] CrossSection = %f !\n", CrossSection);
    printf("[INFO] Equivalent lumi. = %f !\n", (double)nentries/CrossSection);
    printf("[INFO] BranchingFraction = %f !\n", BranchingFraction);
    double TotalGenweight=0;
    for(int ientry=0; ientry<nentries; ientry++){
        //treeReader.GetTChain()->GetEntry(ientry);
        TChain* tmp = treeReader.GetTChain();
        tmp->GetEntry(ientry);
        TotalGenweight+=treeReader.GetGenWeight();
    }
    NormalizationFactor = 1000. * Luminosity * CrossSection * BranchingFraction / TotalGenweight;
    printf("[INFO] TotalGenweight = %f!\n", TotalGenweight);
    printf("[INFO] NormalizationFactor = %f!\n", isData ? 1. : NormalizationFactor);


    //##################################################//
    //#######    Parameters used to report     #########//
    //##################################################//
    int Nevents_pass_selection = 0;
    int Nevents_CUT_num_bjets_is_not_exactly_one = 0;
    int Nevents_anti_CUT_num_bjets_is_not_exactly_one = 0;
    int Nevents_CUT_no_bjet_events = 0;
    int Nevents_anti_CUT_no_bjet_events = 0;
    int Nevents_CUT_num_nonbjets_less_than_2 = 0;
    int Nevents_anti_CUT_num_nonbjets_less_than_2 = 0;
    int Nevents_CUT_nand = 0;
    int Nevents_anti_CUT_nand = 0;

    //##################################################//
    //#########    Event Loop [Selection]    ###########//
    //##################################################//
    // Goal 1: t->b+W(jj), 1 bjet + 2 chi2 jets 
    // Goal 2: diphoton info
    //int nentries = treeReader.flashggStdTree->GetEntries(); printf("[INFO] N_entries = %d\n", nentries);
    for(int ientry=0; ientry<nentries; ientry++){
        treeReader.flashggStdTree->GetEntry(ientry);//load data
        if((ientry+1)%1000==0 || (ientry+1)==nentries) printf("ientry = %d\r", ientry);
        //==================================================//
        //-------------   Reset Parameters   ---------------//
        //==================================================//
        mytree.Clear();
        //------------------------
        std::vector<int> vec_bjet_indices;
        std::vector<TLorentzVector> vec_btagged_jets;
        std::vector<TLorentzVector> vec_nonbtagged_jets;
        //==================================================//
        //--------------   Basic Selectoin   ---------------//
        //==================================================//
        ////if(treeReader.DiPhoInfo_mass<100) continue;
        ////if(treeReader.DiPhoInfo_mass<0) continue;
        //if(treeReader.DiPhoInfo_mass<100 || treeReader.DiPhoInfo_mass>180) continue;
        //if( !isMCsignal && treeReader.DiPhoInfo_mass>120 && treeReader.DiPhoInfo_mass<130) continue;
        //if(!(treeReader.DiPhoInfo_leadIDMVA>0)) continue;
        //if(!(treeReader.DiPhoInfo_subleadIDMVA>0)) continue;
        ////if( !(treeReader.jets_size>0) ) continue;
        //==================================================//
        //------------   Physical Observables   ------------//
        //==================================================//
        //==================================================//
        //------------   Normalization factor   ------------//
        //==================================================//
        //------------   Physical Observables   ------------//
        //==================================================//
        //==================================================//
        //------------   Normalization factor   ------------//
        //==================================================//
        mytree.EvtInfo_totalEntry_before_preselection = nentries;
        mytree.EvtInfo_NormalizationFactor_lumi = isData ? 1. : NormalizationFactor;
        
        //==================================================//
        //-----  Reconstruction(tbW): Select one bjet  -----//
        //==================================================//
        int bjetindex=-1;
        TLorentzVector leading_bjet;
        for(int i=0; i<treeReader.jets_size; i++){
            if( fabs(treeReader.JetInfo_Eta->at(i)) > 2.4 ) continue;
            if( fabs(treeReader.JetInfo_Pt->at(i))  < 30  ) continue;
            if(treeReader.JetInfo_pfCombinedInclusiveSecondaryVertexV2BJetTags->at(i) >= pfCombinedInclusiveSecondaryVertexV2BJetTags_tight){
                if(bjetindex==-1) leading_bjet.SetPtEtaPhiE(treeReader.JetInfo_Pt->at(i), treeReader.JetInfo_Eta->at(i), treeReader.JetInfo_Phi->at(i), treeReader.JetInfo_Energy->at(i));
                bjetindex = i;
                TLorentzVector jet_lorentzvector;
                jet_lorentzvector.SetPtEtaPhiE(treeReader.JetInfo_Pt->at(i), treeReader.JetInfo_Eta->at(i), treeReader.JetInfo_Phi->at(i), treeReader.JetInfo_Energy->at(i));
                vec_btagged_jets.push_back(jet_lorentzvector);
                vec_bjet_indices.push_back(bjetindex);
                mytree.num_btagged_jets+=1;
            }
        }
        bool CUT_no_bjet_events = (bjetindex==-1) ? true : false; // discard no-b-jet events
        bool CUT_num_bjets_is_not_exactly_one = (mytree.num_btagged_jets!=1) ? true : false; // exactly 1 b-jet criterion
        mytree.JetInfo_leading_bjet_pt = CUT_no_bjet_events ? -999. : leading_bjet.Pt();
        mytree.JetInfo_leading_bjet_eta = CUT_no_bjet_events ? -999. : leading_bjet.Eta();
        mytree.JetInfo_leading_bjet_phi = CUT_no_bjet_events ? -999. : leading_bjet.Phi();
        //==================================================//
        //-----  Reconstruction(tbW): Store rest jets  -----//
        //==================================================//
        bool isbjet;
        for(int i=0; i<treeReader.jets_size; i++){
            if( fabs(treeReader.JetInfo_Eta->at(i)) > 2.4 ) continue;
            if( fabs(treeReader.JetInfo_Pt->at(i))  < 30  ) continue;
            //--------------------
            // Exclude bjet
            isbjet=false;
            if(!CUT_no_bjet_events){for(int j=0; j<vec_btagged_jets.size(); j++){if(i==vec_bjet_indices[j]) isbjet=true;}}
            if(isbjet) continue;
            //--------------------
            TLorentzVector jet_lorentzvector;
            jet_lorentzvector.SetPtEtaPhiE(treeReader.JetInfo_Pt->at(i), treeReader.JetInfo_Eta->at(i), treeReader.JetInfo_Phi->at(i), treeReader.JetInfo_Energy->at(i));
            vec_nonbtagged_jets.push_back(jet_lorentzvector);
            mytree.num_nonbtagged_jets+=1;
        }
        bool CUT_num_nonbjets_less_than_2 = (mytree.num_nonbtagged_jets<2) ? true : false;
        mytree.num_jets = mytree.num_btagged_jets + mytree.num_nonbtagged_jets;
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
                            if(chi2<chi2_min){wjetindices[0]=i; wjetindices[1]=j; good_bjet_index=k; chi2_min=chi2;}
                        }// end of k loop (b-jet))
                    } else{
                            chi2 = Chi2_calculator_w_only(dijet_invariant_mass);
                            if(chi2<chi2_min){wjetindices[0]=i; wjetindices[1]=j; chi2_min=chi2;}
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
        double dijet_angle_difference = CUT_num_nonbjets_less_than_2 ? -999. : sqrt(dijet_eta_difference*dijet_eta_difference + dijet_phi_difference*dijet_phi_difference);
        //------------------------------
        double tbw_chosen_bjet_pt = CUT_no_bjet_events ? -999. : vec_btagged_jets[good_bjet_index].Pt();
        double tbw_chosen_bjet_eta = CUT_no_bjet_events ? -999. : vec_btagged_jets[good_bjet_index].Eta();
        double tbw_chosen_bjet_phi = CUT_no_bjet_events ? -999. : vec_btagged_jets[good_bjet_index].Phi();
        //------------------------------
        mytree.JetInfo_jet1_pt = dijet_jet1_pt;
        mytree.JetInfo_jet1_eta = dijet_jet1_eta;
        mytree.JetInfo_jet1_phi = dijet_jet1_phi;
        mytree.JetInfo_jet2_pt = dijet_jet2_pt;
        mytree.JetInfo_jet2_eta = dijet_jet2_eta;
        mytree.JetInfo_jet2_phi = dijet_jet2_phi;
        mytree.JetInfo_dijet_delta_eta = dijet_eta_difference;
        mytree.JetInfo_dijet_delta_phi = dijet_phi_difference;
        mytree.JetInfo_dijet_delta_angle = dijet_angle_difference;
        mytree.inv_mass_dijet = dijet_invariant_mass;
        mytree.JetInfo_chosen_bjet_is_leading_bjet = (good_bjet_index == 0) ? 1 : 0;
        mytree.JetInfo_chosen_bjet_pt = tbw_chosen_bjet_pt;
        mytree.JetInfo_chosen_bjet_eta = tbw_chosen_bjet_eta;
        mytree.JetInfo_chosen_bjet_phi = tbw_chosen_bjet_phi;
        //==================================================//
        //-----  Reconstruction(tbW): InvMass of top  ------//
        //==================================================//
        //TLorentzVector tbw_bjj = best_dijet + leading_bjet;
        best_trijet = CUT_no_any_top_candidate ? best_trijet : best_dijet + vec_btagged_jets[good_bjet_index];
        trijet_invariant_mass = CUT_no_any_top_candidate ? -999. : best_trijet.M();
        mytree.inv_mass_tbw = trijet_invariant_mass;
        //==================================================//
        //-----------   Diphoton Related Info    -----------//
        //==================================================//
        mytree.DiPhoInfo_leadPt = treeReader.DiPhoInfo_leadPt;
        mytree.DiPhoInfo_leadEta = treeReader.DiPhoInfo_leadEta;
        mytree.DiPhoInfo_leadPhi = treeReader.DiPhoInfo_leadPhi;
        mytree.DiPhoInfo_leadE = treeReader.DiPhoInfo_leadE;
        mytree.DiPhoInfo_leadIDMVA = treeReader.DiPhoInfo_leadIDMVA;
        mytree.DiPhoInfo_subleadPt = treeReader.DiPhoInfo_subleadPt;
        mytree.DiPhoInfo_subleadEta = treeReader.DiPhoInfo_subleadEta;
        mytree.DiPhoInfo_subleadPhi = treeReader.DiPhoInfo_subleadPhi;
        mytree.DiPhoInfo_subleadE = treeReader.DiPhoInfo_subleadE;
        mytree.DiPhoInfo_subleadIDMVA = treeReader.DiPhoInfo_subleadIDMVA;
        mytree.DiPhoInfo_mass = treeReader.DiPhoInfo_mass;
        mytree.inv_mass_diphoton = treeReader.DiPhoInfo_mass;
        //==================================================//
        //-----------   EventPar Related Info    -----------//
        //==================================================//
        mytree.EvtInfo_passTrigger = treeReader.EvtInfo_passTrigger;
        mytree.EvtInfo_NPu = treeReader.EvtInfo_NPu;
        mytree.EvtInfo_Rho = treeReader.EvtInfo_Rho;
        mytree.EvtInfo_NVtx = treeReader.EvtInfo_NVtx;
        mytree.EvtInfo_genweight = treeReader.EvtInfo_genweight;
        //==================================================//
        //-------------   Event Counting     ---------------//
        //==================================================//
        Nevents_pass_selection += 1;

        mytree.Fill();
    }// End of event loop.

    //==================================================//
    //---------------------  Report  -------------------//
    //==================================================//
    printf("\n[INFO] N_Nevents_pass_selection = %d/%d\n", Nevents_pass_selection, nentries);
    printf("[DEBUG] Nevents_CUT_num_bjets_is_not_exactly_one = %d (%.4f)\n", Nevents_CUT_num_bjets_is_not_exactly_one, (double)Nevents_CUT_num_bjets_is_not_exactly_one/(double)Nevents_pass_selection);
    printf("[DEBUG] Nevents_anti_CUT_num_bjets_is_not_exactly_one = %d (%.4f)\n", Nevents_anti_CUT_num_bjets_is_not_exactly_one, (double)Nevents_anti_CUT_num_bjets_is_not_exactly_one/(double)Nevents_pass_selection);
    printf("[DEBUG] Nevents_CUT_no_bjet_events = %d (%.4f)\n", Nevents_CUT_no_bjet_events, (double)Nevents_CUT_no_bjet_events/(double)Nevents_pass_selection);
    printf("[DEBUG] Nevents_anti_CUT_no_bjet_events = %d (%.4f)\n", Nevents_anti_CUT_no_bjet_events, (double)Nevents_anti_CUT_no_bjet_events/(double)Nevents_pass_selection);
    printf("[DEBUG] Nevents_CUT_num_nonbjets_less_than_2 = %d (%.4f)\n", Nevents_CUT_num_nonbjets_less_than_2, (double)Nevents_CUT_num_nonbjets_less_than_2/(double)Nevents_pass_selection);
    printf("[DEBUG] Nevents_anti_CUT_num_nonbjets_less_than_2 = %d (%.4f)\n", Nevents_anti_CUT_num_nonbjets_less_than_2, (double)Nevents_anti_CUT_num_nonbjets_less_than_2/(double)Nevents_pass_selection);
    printf("[DEBUG] Nevents_CUT_nand = %d (%.4f)\n", Nevents_CUT_nand, (double)Nevents_CUT_nand/(double)Nevents_pass_selection);
    printf("[DEBUG] Nevents_anti_CUT_nand = %d (%.4f)\n", Nevents_anti_CUT_nand, (double)Nevents_anti_CUT_nand/(double)Nevents_pass_selection);


    fout->Write();
    fout->Close();
    return 1;
}



double Chi2_calculator(double w_mass, double t_mass){
    return (w_mass-w_boson_mass)*(w_mass-w_boson_mass)/w_boson_width + (t_mass-top_quark_mass)*(t_mass-top_quark_mass)/top_quark_width;
}
double Chi2_calculator_w_only(double w_mass){
    return (w_mass-w_boson_mass)*(w_mass-w_boson_mass);
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
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016B-07Aug17_ver2-v2") return true;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016C-07Aug17-v1") return true;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016D-07Aug17-v1") return true;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016E-07Aug17-v1") return true;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016F-07Aug17-v1") return true;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016G-07Aug17-v1") return true;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016H-07Aug17-v1") return true;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016B-07Aug17_ver2-v2") return true;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016C-07Aug17-v1") return true;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016D-07Aug17-v1") return true;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016E-07Aug17-v1") return true;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016F-07Aug17-v1") return true;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016G-07Aug17-v1") return true;
    if((string)dataset == "DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016H-07Aug17-v1") return true;
    return false;
}



flashggStdTreeParameters::flashggStdTreeParameters(){
    JetInfo_Pt = new std::vector<float>;
    JetInfo_Eta = new std::vector<float>;
    JetInfo_Phi = new std::vector<float>;
    JetInfo_Mass = new std::vector<float>;
    JetInfo_Energy = new std::vector<float>;
    JetInfo_pfDeepCSVJetTags_probb = new std::vector<float>;
    JetInfo_pfDeepCSVJetTags_probbb = new std::vector<float>;
    JetInfo_pfCombinedInclusiveSecondaryVertexV2BJetTags = new std::vector<float>;
}
flashggStdTreeParameters::~flashggStdTreeParameters(){
    delete JetInfo_Pt;
    delete JetInfo_Eta;
    delete JetInfo_Phi;
    delete JetInfo_Mass;
    delete JetInfo_Energy;
    delete JetInfo_pfDeepCSVJetTags_probb;
    delete JetInfo_pfDeepCSVJetTags_probbb;
    delete JetInfo_pfCombinedInclusiveSecondaryVertexV2BJetTags;
}
flashggStdTreeReader::flashggStdTreeReader(void){
    printf("[INFO] Reading data...\n");
}
void flashggStdTreeReader::InitChain(const char* treeName){
    flashggStdTree = new TChain(treeName);
    printf("[INFO] flashggStdTreeReader::InitChain : Finished!\n");
}
void flashggStdTreeReader::AddSingleRootFile(char* input_file){
    flashggStdTree->Add(input_file);
    printf("[INFO] flashggStdTreeReader::AddSingleRootFile : Finished!\n");
}
void flashggStdTreeReader::AddMultiRootFile(char* input_file){
    flashggStdTree->Add(Form("%s/*.root", input_file));
    printf("[INFO] flashggStdTreeReader::AddMultiRootFile : Finished!\n");
}
int flashggStdTreeReader::GetEntries(void){
    printf("[INFO] flashggStdTreeReader::GetEntries : %d\n", flashggStdTree->GetEntries());
    return flashggStdTree->GetEntries();
}
double flashggStdTreeReader::GetGenWeight(void){
    return EvtInfo_genweight;
}
TChain* flashggStdTreeReader::GetTChain(void){
    return flashggStdTree;
}
void flashggStdTreeReader::SetBranchAddresses(){
    flashggStdTree->SetBranchAddress("EvtInfo.passTrigger", &EvtInfo_passTrigger);
    flashggStdTree->SetBranchAddress("EvtInfo.NPu", &EvtInfo_NPu);
    flashggStdTree->SetBranchAddress("EvtInfo.Rho", &EvtInfo_Rho);
    flashggStdTree->SetBranchAddress("EvtInfo.NVtx", &EvtInfo_NVtx);
    flashggStdTree->SetBranchAddress("EvtInfo.genweight", &EvtInfo_genweight);
    flashggStdTree->SetBranchAddress("jets_size", &jets_size);
    flashggStdTree->SetBranchAddress("JetInfo.Pt", &JetInfo_Pt);
    flashggStdTree->SetBranchAddress("JetInfo.Eta", &JetInfo_Eta);
    flashggStdTree->SetBranchAddress("JetInfo.Phi", &JetInfo_Phi);
    flashggStdTree->SetBranchAddress("JetInfo.Mass", &JetInfo_Mass);
    flashggStdTree->SetBranchAddress("JetInfo.Energy", &JetInfo_Energy);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probb", &JetInfo_pfDeepCSVJetTags_probb);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probbb", &JetInfo_pfDeepCSVJetTags_probbb);
    flashggStdTree->SetBranchAddress("JetInfo.pfCombinedInclusiveSecondaryVertexV2BJetTags", &JetInfo_pfCombinedInclusiveSecondaryVertexV2BJetTags);
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
    printf("[INFO] flashggStdTreeReader::SetBranchAddresses : Finished!\n");
}



void myTreeClass::InitTree(){
    mytree =  new TTree("mytree", "mytree");
}
void myTreeClass::MakeNewBranchAddresses(){
    mytree -> Branch("EvtInfo_totalEntry_before_preselection", &EvtInfo_totalEntry_before_preselection, "EvtInfo_totalEntry_before_preselection/I");
    mytree -> Branch("EvtInfo_NormalizationFactor_lumi", &EvtInfo_NormalizationFactor_lumi, "EvtInfo_NormalizationFactor_lumi/F");
    mytree -> Branch("EvtInfo_passTrigger", &EvtInfo_passTrigger, "EvtInfo_passTrigger/F");
    mytree -> Branch("EvtInfo_NPu", &EvtInfo_NPu, "EvtInfo_NPu/F");
    mytree -> Branch("EvtInfo_Rho", &EvtInfo_Rho, "EvtInfo_Rho/F");
    mytree -> Branch("EvtInfo_NVtx", &EvtInfo_NVtx, "EvtInfo_NVtx/I");
    mytree -> Branch("EvtInfo_genweight", &EvtInfo_genweight, "EvtInfo_genweight/F");
    mytree -> Branch("num_jets", &num_jets, "num_jets/I");
    mytree -> Branch("num_btagged_jets", &num_btagged_jets, "num_btagged_jets/I");
    mytree -> Branch("num_nonbtagged_jets", &num_nonbtagged_jets, "num_nonbtagged_jets/I");
    //------------------------
    mytree -> Branch("inv_mass_dijet", &inv_mass_dijet, "inv_mass_dijet/F");
    mytree -> Branch("inv_mass_diphoton", &inv_mass_diphoton, "inv_mass_diphoton/F");
    mytree -> Branch("inv_mass_tbw", &inv_mass_tbw, "inv_mass_tbw/F");
    //------------------------
    mytree -> Branch("JetInfo_leading_bjet_pt", &JetInfo_leading_bjet_pt, "JetInfo_leading_bjet_pt/F");
    mytree -> Branch("JetInfo_leading_bjet_eta", &JetInfo_leading_bjet_eta, "JetInfo_leading_bjet_eta/F");
    mytree -> Branch("JetInfo_leading_bjet_phi", &JetInfo_leading_bjet_phi, "JetInfo_leading_bjet_phi/F");
    mytree -> Branch("JetInfo_chosen_bjet_is_leading_bjet", &JetInfo_chosen_bjet_is_leading_bjet, "JetInfo_chosen_bjet_is_leading_bjet/I");
    mytree -> Branch("JetInfo_chosen_bjet_pt", &JetInfo_chosen_bjet_pt, "JetInfo_chosen_bjet_pt/F");
    mytree -> Branch("JetInfo_chosen_bjet_eta", &JetInfo_chosen_bjet_eta, "JetInfo_chosen_bjet_eta/F");
    mytree -> Branch("JetInfo_chosen_bjet_phi", &JetInfo_chosen_bjet_phi, "JetInfo_chosen_bjet_phi/F");
    mytree -> Branch("JetInfo_jet1_pt", &JetInfo_jet1_pt, "JetInfo_jet1_pt/F");
    mytree -> Branch("JetInfo_jet1_eta", &JetInfo_jet1_eta, "JetInfo_jet1_eta/F");
    mytree -> Branch("JetInfo_jet1_phi", &JetInfo_jet1_phi, "JetInfo_jet1_phi/F");
    mytree -> Branch("JetInfo_jet2_pt", &JetInfo_jet2_pt, "JetInfo_jet2_pt/F");
    mytree -> Branch("JetInfo_jet2_eta", &JetInfo_jet2_eta, "JetInfo_jet2_eta/F");
    mytree -> Branch("JetInfo_jet2_phi", &JetInfo_jet2_phi, "JetInfo_jet2_phi/F");
    //------------------------
    mytree -> Branch("JetInfo_dijet_delta_eta", &JetInfo_dijet_delta_eta, "JetInfo_dijet_delta_eta/F");
    mytree -> Branch("JetInfo_dijet_delta_phi", &JetInfo_dijet_delta_phi, "JetInfo_dijet_delta_phi/F");
    mytree -> Branch("JetInfo_dijet_delta_angle", &JetInfo_dijet_delta_angle, "JetInfo_dijet_delta_angle/F");
    //------------------------
    mytree -> Branch("DiPhoInfo_leadPt", &DiPhoInfo_leadPt, "DiPhoInfo_leadPt/F");
    mytree -> Branch("DiPhoInfo_leadEta", &DiPhoInfo_leadEta, "DiPhoInfo_leadEta/F");
    mytree -> Branch("DiPhoInfo_leadPhi", &DiPhoInfo_leadPhi, "DiPhoInfo_leadPhi/F");
    mytree -> Branch("DiPhoInfo_leadE", &DiPhoInfo_leadE, "DiPhoInfo_leadE/F");
    mytree -> Branch("DiPhoInfo_leadIDMVA", &DiPhoInfo_leadIDMVA, "DiPhoInfo_leadIDMVA/F");
    mytree -> Branch("DiPhoInfo_subleadPt", &DiPhoInfo_subleadPt, "DiPhoInfo_subleadPt/F");
    mytree -> Branch("DiPhoInfo_subleadEta", &DiPhoInfo_subleadEta, "DiPhoInfo_subleadEta/F");
    mytree -> Branch("DiPhoInfo_subleadPhi", &DiPhoInfo_subleadPhi, "DiPhoInfo_subleadPhi/F");
    mytree -> Branch("DiPhoInfo_subleadE", &DiPhoInfo_subleadE, "DiPhoInfo_subleadE/F");
    mytree -> Branch("DiPhoInfo_subleadIDMVA", &DiPhoInfo_subleadIDMVA, "DiPhoInfo_subleadIDMVA/F");
    mytree -> Branch("DiPhoInfo_mass", &DiPhoInfo_mass, "DiPhoInfo_mass/F");
}
void myTreeClass::Fill(){
    mytree -> Fill();
}
void myParameters::Clear(){
    EvtInfo_totalEntry_before_preselection = 0;
    EvtInfo_NormalizationFactor_lumi = 0;
    //------------------------
    num_jets = 0;
    num_btagged_jets = 0;
    num_nonbtagged_jets = 0;
    //------------------------
    inv_mass_dijet = 0;
    inv_mass_diphoton = 0;
    inv_mass_tbw = 0;
    //------------------------
    JetInfo_dijet_delta_eta = 0;
    JetInfo_dijet_delta_phi = 0;
    JetInfo_dijet_delta_angle = 0;
    //------------------------
    JetInfo_leading_bjet_pt = 0;
    JetInfo_leading_bjet_eta = 0;
    JetInfo_leading_bjet_phi = 0;
    //------------------------
    JetInfo_chosen_bjet_is_leading_bjet = -1;
    JetInfo_chosen_bjet_pt = 0;
    JetInfo_chosen_bjet_eta = 0;
    JetInfo_chosen_bjet_phi = 0;
    //------------------------
    JetInfo_jet1_pt = 0;
    JetInfo_jet1_eta = 0;
    JetInfo_jet1_phi = 0;
    //------------------------
    JetInfo_jet2_pt = 0;
    JetInfo_jet2_eta = 0;
    JetInfo_jet2_phi = 0;
}
