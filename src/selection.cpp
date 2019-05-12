#include <stdio.h>
#include <TFile.h>
#include <TH1D.h>
#include <TChain.h>
#include <TCanvas.h>
#include "../include/main.h"
#include "../include/selection.h"
using namespace std;

void Selection(){
    char input_file[256]; sprintf(input_file, "%s", "ntuples_skimmed/ntuple_ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root"); printf("[INFO] input_file  = %s\n", input_file);
    char output_file[256]; sprintf(output_file, "%s", "plots/test.root"); printf("[INFO] output_file = %s\n", output_file);
    char dataset[256]; sprintf(dataset, "%s", "ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8"); printf("[INFO] dataset     = %s\n", dataset);
    bool isData = isThisDataOrNot(dataset);

    TFile *fin  = TFile::Open(input_file);
    TFile *fout = new TFile(output_file, "RECREATE");
    myTreeClass treeReader;
    treeReader.InitTree("mytree");
    treeReader.AddRootFile(fin);
    treeReader.SetBranchAddresses();

    //==================//
    //--- histograms ---//
    //==================//
    TH1D  *hist_NPu = new TH1D("hist_NPu", "hist_NPu", 100, 0, 100);
    TH1D  *hist_NVtx = new TH1D("hist_NVtx", "hist_NVtx", 100, 0, 100);
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
    hist_NPu -> Sumw2();
    hist_NVtx -> Sumw2();
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


    int nentries = treeReader.GetEntries();
    for(int ientry=0; ientry<nentries; ientry++){
        TTree* tmp = treeReader.GetTTree();
        tmp->GetEntry(ientry);
        if(ientry==0) printf("[INFO] N_entries = %d/%d\n", nentries, treeReader.EvtInfo_totalEntry_before_preselection);
        if((ientry+1)%10000==0 || (ientry+1)==nentries) printf("ientry = %d\r", ientry);
        double NormalizationFactor = treeReader.EvtInfo_genweight * treeReader.EvtInfo_NormalizationFactor_lumi;
        //EvtInfo_NormalizationFactor_lumi = 1000. * Luminosity * CrossSection * BranchingFraction / TotalGenweight;
        //=== Selections ===//
        //=== Store Info ===//
        //for(unsigned int i=0; i<treeReader.JetInfo_bjet_pt.size(); i++){
        //    hist_bjet_pt -> Fill(treeReader.JetInfo_bjet_pt[i]);
        //    hist_bjet_eta -> Fill(treeReader.JetInfo_bjet_eta[i]);
        //    hist_bjet_phi -> Fill(treeReader.JetInfo_bjet_phi[i]);
        //}
        //hist_jet1_pt -> Fill(treeReader.jet1_pt);
        //hist_jet2_pt -> Fill(treeReader.jet2_pt);
        //hist_jet1_eta -> Fill(treeReader.jet1_eta);
        //hist_jet2_eta -> Fill(treeReader.jet2_eta);
        //hist_jet1_phi -> Fill(treeReader.jet1_phi);
        //hist_jet2_phi -> Fill(treeReader.jet2_phi);
        //hist_dijet_eta -> Fill(treeReader.dijet_eta);
        //hist_dijet_phi -> Fill(treeReader.dijet_phi);
        //hist_dijet_angle -> Fill(treeReader.dijet_angle);
        hist_NPu  -> Fill(treeReader.EvtInfo_NPu, isData ? 1. : NormalizationFactor);
        hist_NVtx -> Fill(treeReader.EvtInfo_NVtx, isData ? 1. : NormalizationFactor);
        hist_num_jets -> Fill(treeReader.num_jets, isData ? 1. : NormalizationFactor);
        hist_num_btagged_jets -> Fill(treeReader.num_btagged_jets, isData ? 1. : NormalizationFactor);
        hist_num_nonbtagged_jets -> Fill(treeReader.num_nonbtagged_jets, isData ? 1. : NormalizationFactor);
        hist_inv_mass_dijet -> Fill(treeReader.inv_mass_dijet, isData ? 1. : NormalizationFactor);
        hist_inv_mass_diphoton -> Fill(treeReader.inv_mass_diphoton, isData ? 1. : NormalizationFactor);
        hist_inv_mass_tbw -> Fill(treeReader.inv_mass_tbw, isData ? 1. : NormalizationFactor);
        hist_DiPhoInfo_leadPt -> Fill(treeReader.DiPhoInfo_leadPt, isData ? 1. : NormalizationFactor);
        hist_DiPhoInfo_leadEta -> Fill(treeReader.DiPhoInfo_leadEta, isData ? 1. : NormalizationFactor);
        hist_DiPhoInfo_leadPhi -> Fill(treeReader.DiPhoInfo_leadPhi, isData ? 1. : NormalizationFactor);
        hist_DiPhoInfo_leadE -> Fill(treeReader.DiPhoInfo_leadE, isData ? 1. : NormalizationFactor);
        hist_DiPhoInfo_leadIDMVA -> Fill(treeReader.DiPhoInfo_leadIDMVA, isData ? 1. : NormalizationFactor);
        hist_DiPhoInfo_subleadPt -> Fill(treeReader.DiPhoInfo_subleadPt, isData ? 1. : NormalizationFactor);
        hist_DiPhoInfo_subleadEta -> Fill(treeReader.DiPhoInfo_subleadEta, isData ? 1. : NormalizationFactor);
        hist_DiPhoInfo_subleadPhi -> Fill(treeReader.DiPhoInfo_subleadPhi, isData ? 1. : NormalizationFactor);
        hist_DiPhoInfo_subleadE -> Fill(treeReader.DiPhoInfo_subleadE, isData ? 1. : NormalizationFactor);
        hist_DiPhoInfo_subleadIDMVA -> Fill(treeReader.DiPhoInfo_subleadIDMVA, isData ? 1. : NormalizationFactor);
        //hist_DiPhoInfo_leadIDMVA_ori -> Fill(treeReader.DiPhoInfo_leadIDMVA_ori, isData ? 1. : NormalizationFactor);
        //hist_DiPhoInfo_subleadIDMVA_ori -> Fill(treeReader.DiPhoInfo_subleadIDMVA_ori, isData ? 1. : NormalizationFactor);
        //hist_inv_mass_diphoton_ori -> Fill(treeReader.inv_mass_diphoton_ori, isData ? 1. : NormalizationFactor);
    }

    char output_dir[256] = "plots"; printf("[INFO] output_dir = %s\n", output_dir);
    TCanvas *c1 = new TCanvas("c1", "c1", 700, 800);
    //MakePlots(c1, hist_bjet_pt, "hist_bjet_pt", Form("%s/hist_bjet_pt.png", output_dir));
    //MakePlots(c1, hist_bjet_eta, "hist_bjet_eta", Form("%s/hist_bjet_eta.png", output_dir));
    //MakePlots(c1, hist_bjet_phi, "hist_bjet_phi", Form("%s/hist_bjet_phi.png", output_dir));
    MakePlots(c1, hist_NPu, "hist_NPu", Form("%s/hist_NPu.png", output_dir));
    MakePlots(c1, hist_NVtx, "hist_NVtx", Form("%s/hist_NVtx.png", output_dir));
    MakePlots(c1, hist_num_jets, "hist_num_jets", Form("%s/hist_num_jets.png", output_dir));
    MakePlots(c1, hist_num_btagged_jets, "hist_num_btagged_jets", Form("%s/hist_num_btagged_jets.png", output_dir));
    MakePlots(c1, hist_num_nonbtagged_jets, "hist_num_nonbtagged_jets", Form("%s/hist_num_nonbtagged_jets.png", output_dir));
    MakePlots(c1, hist_inv_mass_dijet, "hist_inv_mass_dijet", Form("%s/hist_inv_mass_dijet.png", output_dir));
    MakePlots(c1, hist_inv_mass_diphoton, "hist_inv_mass_diphoton", Form("%s/hist_inv_mass_diphoton.png", output_dir));
    MakePlots(c1, hist_inv_mass_tbw, "hist_inv_mass_tbw", Form("%s/hist_inv_mass_tbw.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_leadPt, "hist_DiPhoInfo_leadPt", Form("%s/hist_DiPhoInfo_leadPt.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_leadEta, "hist_DiPhoInfo_leadEta", Form("%s/hist_DiPhoInfo_leadEta.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_leadPhi, "hist_DiPhoInfo_leadPhi", Form("%s/hist_DiPhoInfo_leadPhi.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_leadE, "hist_DiPhoInfo_leadE", Form("%s/hist_DiPhoInfo_leadE.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_leadIDMVA, "hist_DiPhoInfo_leadIDMVA", Form("%s/hist_DiPhoInfo_leadIDMVA.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_subleadPt, "hist_DiPhoInfo_subleadPt", Form("%s/hist_DiPhoInfo_subleadPt.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_subleadEta, "hist_DiPhoInfo_subleadEta", Form("%s/hist_DiPhoInfo_subleadEta.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_subleadPhi, "hist_DiPhoInfo_subleadPhi", Form("%s/hist_DiPhoInfo_subleadPhi.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_subleadE, "hist_DiPhoInfo_subleadE", Form("%s/hist_DiPhoInfo_subleadE.png", output_dir));
    MakePlots(c1, hist_DiPhoInfo_subleadIDMVA, "hist_DiPhoInfo_subleadIDMVA", Form("%s/hist_DiPhoInfo_subleadIDMVA.png", output_dir));

    fout -> Close();
}

int main(){
    Selection();
    return 1;
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
bool isThisDataOrNot(char* dataset){
    if((string)dataset == "DoubleEG_B") return true;
    if((string)dataset == "DoubleEG_C") return true;
    if((string)dataset == "DoubleEG_D") return true;
    if((string)dataset == "DoubleEG_E") return true;
    if((string)dataset == "DoubleEG_F") return true;
    return false;
}


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
    mytree -> SetBranchAddress("EvtInfo_NVtx", &EvtInfo_NVtx);
    mytree -> SetBranchAddress("EvtInfo_genweight", &EvtInfo_genweight);
    //------------------------
    mytree -> SetBranchAddress("num_jets", &num_jets);
    mytree -> SetBranchAddress("num_btagged_jets", &num_btagged_jets);
    mytree -> SetBranchAddress("num_nonbtagged_jets", &num_nonbtagged_jets);
    //------------------------
    mytree -> SetBranchAddress("inv_mass_dijet", &inv_mass_dijet);
    mytree -> SetBranchAddress("inv_mass_diphoton", &inv_mass_diphoton);
    mytree -> SetBranchAddress("inv_mass_tbw", &inv_mass_tbw);
    //------------------------
    //mytree -> SetBranchAddress("JetInfo_bjet_pt", &JetInfo_bjet_pt);
    //mytree -> SetBranchAddress("JetInfo_bjet_eta", &JetInfo_bjet_eta);
    //mytree -> SetBranchAddress("JetInfo_bjet_phi", &JetInfo_bjet_phi);
    //mytree -> SetBranchAddress("JetInfo_jet1_pt", &JetInfo_jet1_pt);
    //mytree -> SetBranchAddress("JetInfo_jet1_eta", &JetInfo_jet1_eta);
    //mytree -> SetBranchAddress("JetInfo_jet1_phi", &JetInfo_jet1_phi);
    //mytree -> SetBranchAddress("JetInfo_jet2_pt", &JetInfo_jet2_pt);
    //mytree -> SetBranchAddress("JetInfo_jet2_eta", &JetInfo_jet2_eta);
    //mytree -> SetBranchAddress("JetInfo_jet2_phi", &JetInfo_jet2_phi);
    //------------------------
    mytree -> SetBranchAddress("JetInfo_dijet_delta_eta", &JetInfo_dijet_delta_eta);
    mytree -> SetBranchAddress("JetInfo_dijet_delta_phi", &JetInfo_dijet_delta_phi);
    mytree -> SetBranchAddress("JetInfo_dijet_delta_angle", &JetInfo_dijet_delta_angle);
    //------------------------
    mytree -> SetBranchAddress("DiPhoInfo_leadPt", &DiPhoInfo_leadPt);
    mytree -> SetBranchAddress("DiPhoInfo_leadEta", &DiPhoInfo_leadEta);
    mytree -> SetBranchAddress("DiPhoInfo_leadPhi", &DiPhoInfo_leadPhi);
    mytree -> SetBranchAddress("DiPhoInfo_leadE", &DiPhoInfo_leadE);
    mytree -> SetBranchAddress("DiPhoInfo_leadIDMVA", &DiPhoInfo_leadIDMVA);
    mytree -> SetBranchAddress("DiPhoInfo_subleadPt", &DiPhoInfo_subleadPt);
    mytree -> SetBranchAddress("DiPhoInfo_subleadEta", &DiPhoInfo_subleadEta);
    mytree -> SetBranchAddress("DiPhoInfo_subleadPhi", &DiPhoInfo_subleadPhi);
    mytree -> SetBranchAddress("DiPhoInfo_subleadE", &DiPhoInfo_subleadE);
    mytree -> SetBranchAddress("DiPhoInfo_subleadIDMVA", &DiPhoInfo_subleadIDMVA);
    mytree -> SetBranchAddress("DiPhoInfo_mass", &DiPhoInfo_mass);
    printf("[INFO] myTreeClass::SetBranchAddresses : Finished!\n");
}
flashggStdTreeParameters::flashggStdTreeParameters(){
    JetInfo_Pt = new std::vector<float>;
    JetInfo_Eta = new std::vector<float>;
    JetInfo_Phi = new std::vector<float>;
    JetInfo_Mass = new std::vector<float>;
    JetInfo_Energy = new std::vector<float>;
    JetInfo_pfDeepCSVJetTags_probb = new std::vector<float>;
    JetInfo_pfDeepCSVJetTags_probbb = new std::vector<float>;
}
flashggStdTreeParameters::~flashggStdTreeParameters(){
    delete JetInfo_Pt;
    delete JetInfo_Eta;
    delete JetInfo_Phi;
    delete JetInfo_Mass;
    delete JetInfo_Energy;
    delete JetInfo_pfDeepCSVJetTags_probb;
    delete JetInfo_pfDeepCSVJetTags_probbb;
}
