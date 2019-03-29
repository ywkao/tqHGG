//===========================================//
//=== 22. Mar. 2019                       ===//
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
double w_boson_width = 2.085;//GeV
double top_quark_mass = 173.0 ;//GeV
double top_quark_width = 1.41 ;//GeV


int main(int argc, char *argv[]){
    char input_file[256]; sprintf(input_file, "%s", argv[1]); printf("[INFO] input_file  = %s\n", input_file);
    char output_file[256]; sprintf(output_file, "%s", argv[2]); printf("[INFO] output_file = %s\n", output_file);
    char dataset[256]; sprintf(dataset, "%s", argv[3]); printf("[INFO] dataset     = %s\n", dataset);
    bool isData = isThisDataOrNot(dataset);
    bool isMCsignal = isThisMCsignal(dataset);
    bool isMultiFile = isThisMultiFile(dataset);
    TFile *fout = new TFile(output_file, "RECREATE");

    flashggStdTreeReader treeReader;
    treeReader.InitChain("flashggNtuples/flashggStdTree");
    if(isMultiFile) treeReader.AddSingleRootFile(input_file);
    else            treeReader.AddMultiRootFile(input_file);
    treeReader.SetBranchAddresses();

    myTreeClass mytree;
    mytree.InitTree();
    mytree.SetBranchAddresses();

    /*
    TChain *flashggStdTree = new TChain("flashggNtuples/flashggStdTree");
    if(isMultiFile) flashggStdTree->Add(Form("%s/*.root", input_file));
    else flashggStdTree->Add(input_file);
    */

    //==========================//
    //--- Readout Parameters ---//
    //==========================//
    /*
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
    */


    //==================//
    //--- MyTreeVar ----//
    //==================//
    /*
    Int_t num_jets;
    Int_t num_btagged_jets;
    Int_t num_nonbtagged_jets;
    //------------------------
    float inv_mass_dijet;
    float inv_mass_diphoton;
    float inv_mass_tbw;
    float JetInfo_dijet_delta_eta;
    float JetInfo_dijet_delta_phi;
    float JetInfo_dijet_delta_angle;
    //------------------------
    std::vector<float> JetInfo_bjet_pt;
    std::vector<float> JetInfo_bjet_eta;
    std::vector<float> JetInfo_bjet_phi;
    std::vector<float> JetInfo_jet1_pt;
    std::vector<float> JetInfo_jet1_eta;
    std::vector<float> JetInfo_jet1_phi;
    std::vector<float> JetInfo_jet2_pt;
    std::vector<float> JetInfo_jet2_eta;
    std::vector<float> JetInfo_jet2_phi;
    */

    //==================//
    //---- Branches ----//
    //==================//
    /*
    TTree *mytree = new TTree("mytree", "mytree");
    //------------------------
    mytree -> Branch("num_jets", &num_jets, "num_jets/I");
    mytree -> Branch("num_btagged_jets", &num_btagged_jets, "num_btagged_jets/I");
    mytree -> Branch("num_nonbtagged_jets", &num_nonbtagged_jets, "num_nonbtagged_jets/I");
    //------------------------
    mytree -> Branch("inv_mass_dijet", &inv_mass_dijet, "inv_mass_dijet/F");
    mytree -> Branch("inv_mass_diphoton", &inv_mass_diphoton, "inv_mass_diphoton/F");
    mytree -> Branch("inv_mass_tbw", &inv_mass_tbw, "inv_mass_tbw/F");
    //------------------------
    mytree -> Branch("JetInfo_bjet_pt", &JetInfo_bjet_pt);
    mytree -> Branch("JetInfo_bjet_eta", &JetInfo_bjet_eta);
    mytree -> Branch("JetInfo_bjet_phi", &JetInfo_bjet_phi);
    mytree -> Branch("JetInfo_jet1_pt", &JetInfo_jet1_pt);
    mytree -> Branch("JetInfo_jet1_eta", &JetInfo_jet1_eta);
    mytree -> Branch("JetInfo_jet1_phi", &JetInfo_jet1_phi);
    mytree -> Branch("JetInfo_jet2_pt", &JetInfo_jet2_pt);
    mytree -> Branch("JetInfo_jet2_eta", &JetInfo_jet2_eta);
    mytree -> Branch("JetInfo_jet2_phi", &JetInfo_jet2_phi);
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
    */

    //##################################################//
    //########   Event Loop [Normalization]   ##########//
    //##################################################//
    /*
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
    */

    //##################################################//
    //#######    Parameters used to report     #########//
    //##################################################//
    int nevents_pass_selection = 0;
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
    int nentries = treeReader.flashggStdTree->GetEntries(); printf("[INFO] N_entries = %d\n", nentries);
    for(int ientry=0; ientry<nentries; ientry++){
        treeReader.flashggStdTree->GetEntry(ientry);//load data
        if((ientry+1)%1000==0 || (ientry+1)==nentries) printf("ientry = %d\r", ientry);
        //==================================================//
        //-------------   Reset Parameters   ---------------//
        //==================================================//
        std::vector<int> vec_bjet_indices;
        std::vector<TLorentzVector> vec_btagged_jets;
        std::vector<TLorentzVector> vec_nonbtagged_jets;
        //vec_bjet_indices.clear();
        //vec_btagged_jets.clear();
        //vec_nonbtagged_jets.clear();
        //------------------------
        num_jets = 0;
        num_btagged_jets = 0;
        num_nonbtagged_jets = 0;
        inv_mass_dijet = 0;
        inv_mass_diphoton = 0;
        inv_mass_tbw = 0;
        JetInfo_dijet_delta_eta = 0;
        JetInfo_dijet_delta_phi = 0;
        JetInfo_dijet_delta_angle = 0;
        //------------------------
        JetInfo_bjet_pt.clear();
        JetInfo_bjet_eta.clear();
        JetInfo_bjet_phi.clear();
        JetInfo_jet1_pt.clear();
        JetInfo_jet1_eta.clear();
        JetInfo_jet1_phi.clear();
        JetInfo_jet2_pt.clear();
        JetInfo_jet2_eta.clear();
        JetInfo_jet2_phi.clear();
        //==================================================//
        //--------------   Basic Selectoin   ---------------//
        //==================================================//
        if(DiPhoInfo_mass<0) continue;
        if( !(jets_size>0) ) continue;
        if(DiPhoInfo_mass<100 || DiPhoInfo_mass>150) continue;
        if( !isMCsignal && DiPhoInfo_mass>120 && DiPhoInfo_mass<130) continue;
        //if(!(DiPhoInfo_leadIDMVA>0)) continue;
        //==================================================//
        //------------   Physical Observables   ------------//
        //==================================================//
        
        //==================================================//
        //-----  Reconstruction(tbW): Select one bjet  -----//
        //==================================================//
        int bjetindex=-1;
        num_btagged_jets = 0;
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
                num_btagged_jets+=1;
            }
        }
        bool CUT_no_bjet_events = (bjetindex==-1) ? true : false; // discard no-b-jet events
        bool CUT_num_bjets_is_not_exactly_one = (num_btagged_jets!=1) ? true : false; // exactly 1 b-jet criterion
        JetInfo_bjet_pt.push_back(leading_bjet.Pt());
        JetInfo_bjet_eta.push_back(leading_bjet.Eta());
        JetInfo_bjet_phi.push_back(leading_bjet.Phi());
        //==================================================//
        //-----  Reconstruction(tbW): Store rest jets  -----//
        //==================================================//
        bool isbjet;
        num_nonbtagged_jets=0;
        for(int i=0; i<jets_size; i++){
            if( fabs(JetInfo_Eta->at(i)) > 2.4 ) continue;
            if( fabs(JetInfo_Pt->at(i))  < 30  ) continue;
            //--------------------
            // Exclude bjet
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
        num_jets = num_btagged_jets + num_nonbtagged_jets;
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
        JetInfo_jet1_pt.push_back(dijet_jet1_pt);
        JetInfo_jet1_eta.push_back(dijet_jet1_eta);
        JetInfo_jet1_phi.push_back(dijet_jet1_phi);
        JetInfo_jet2_pt.push_back(dijet_jet2_pt);
        JetInfo_jet2_eta.push_back(dijet_jet2_eta);
        JetInfo_jet2_phi.push_back(dijet_jet2_phi);
        JetInfo_dijet_delta_eta = dijet_eta_difference;
        JetInfo_dijet_delta_phi = dijet_phi_difference;
        JetInfo_dijet_delta_angle = dijet_angle_difference;
        inv_mass_dijet = dijet_invariant_mass;
        //==================================================//
        //-----  Reconstruction(tbW): InvMass of top  ------//
        //==================================================//
        //TLorentzVector tbw_bjj = best_dijet + leading_bjet;
        best_trijet = CUT_no_any_top_candidate ? best_trijet : best_dijet + vec_btagged_jets[good_bjet_index];
        trijet_invariant_mass = CUT_no_any_top_candidate ? -999. : best_trijet.M();
        inv_mass_tbw = trijet_invariant_mass;
        //==================================================//
        //-----------   Diphoton Related Info    -----------//
        //==================================================//
        // Automatically stored
        //==================================================//
        //-------------   Event Counting     ---------------//
        //==================================================//
        nevents_pass_selection += 1;
        mytree->Fill();
    }// End of event loop.

    //==================================================//
    //---------------------  Report  -------------------//
    //==================================================//
    printf("[INFO] N_nevents_pass_selection = %d/%d\n", nevents_pass_selection, nentries);
    printf("[DEBUG] Nevents_CUT_num_bjets_is_not_exactly_one = %d (%.4f)\n", Nevents_CUT_num_bjets_is_not_exactly_one, (double)Nevents_CUT_num_bjets_is_not_exactly_one/(double)nevents_pass_selection);
    printf("[DEBUG] Nevents_anti_CUT_num_bjets_is_not_exactly_one = %d (%.4f)\n", Nevents_anti_CUT_num_bjets_is_not_exactly_one, (double)Nevents_anti_CUT_num_bjets_is_not_exactly_one/(double)nevents_pass_selection);
    printf("[DEBUG] Nevents_CUT_no_bjet_events = %d (%.4f)\n", Nevents_CUT_no_bjet_events, (double)Nevents_CUT_no_bjet_events/(double)nevents_pass_selection);
    printf("[DEBUG] Nevents_anti_CUT_no_bjet_events = %d (%.4f)\n", Nevents_anti_CUT_no_bjet_events, (double)Nevents_anti_CUT_no_bjet_events/(double)nevents_pass_selection);
    printf("[DEBUG] Nevents_CUT_num_nonbjets_less_than_2 = %d (%.4f)\n", Nevents_CUT_num_nonbjets_less_than_2, (double)Nevents_CUT_num_nonbjets_less_than_2/(double)nevents_pass_selection);
    printf("[DEBUG] Nevents_anti_CUT_num_nonbjets_less_than_2 = %d (%.4f)\n", Nevents_anti_CUT_num_nonbjets_less_than_2, (double)Nevents_anti_CUT_num_nonbjets_less_than_2/(double)nevents_pass_selection);
    printf("[DEBUG] Nevents_CUT_nand = %d (%.4f)\n", Nevents_CUT_nand, (double)Nevents_CUT_nand/(double)nevents_pass_selection);
    printf("[DEBUG] Nevents_anti_CUT_nand = %d (%.4f)\n", Nevents_anti_CUT_nand, (double)Nevents_anti_CUT_nand/(double)nevents_pass_selection);


    fout->Write();
    fout->Close();
    return 1;
}

double Chi2_calculator(double w_mass, double t_mass){
    return (w_mass-w_boson_mass)*(w_mass-w_boson_mass)/w_boson_width + (t_mass-top_quark_mass)*(t_mass-top_quark_mass)/top_quark_width;
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
