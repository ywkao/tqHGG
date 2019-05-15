#include <stdio.h>
#include <TFile.h>
#include <TH1D.h>
#include <TChain.h>
#include <TCanvas.h>
#include "../include/main.h"
#include "../include/selection.h"
#include "../include/enumhist.h"
using namespace std;

void Selection(char* input_file, char* output_file, char* dataset, char* output_dir){
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
    TH1D* h[totalHistNum];
    for(int i=0; i<totalHistNum; i++){
        //printf("[%d] %s, %d, %f, %f\n", i, histNames[i].c_str(), histNbins[i], histBinLow[i], histBinHigh[i]);
        h[i] = new TH1D(histNames[i].c_str(), histNames[i].c_str(), histNbins[i], histBinLow[i], histBinHigh[i]);
        h[i] -> Sumw2();
    }

    //==================//
    //--- Event Loop ---//
    //==================//
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
        h[hist_NPu]  -> Fill(treeReader.EvtInfo_NPu, isData ? 1. : NormalizationFactor);
        h[hist_NVtx] -> Fill(treeReader.EvtInfo_NVtx, isData ? 1. : NormalizationFactor);
        h[hist_num_jets] -> Fill(treeReader.num_jets, isData ? 1. : NormalizationFactor);
        h[hist_num_btagged_jets] -> Fill(treeReader.num_btagged_jets, isData ? 1. : NormalizationFactor);
        h[hist_num_nonbtagged_jets] -> Fill(treeReader.num_nonbtagged_jets, isData ? 1. : NormalizationFactor);
        if(!(treeReader.JetInfo_leading_bjet_pt<-900.)){
            h[hist_leading_bjet_pt] -> Fill(treeReader.JetInfo_leading_bjet_pt, isData ? 1. : NormalizationFactor);
            h[hist_leading_bjet_eta] -> Fill(treeReader.JetInfo_leading_bjet_eta, isData ? 1. : NormalizationFactor);
            h[hist_leading_bjet_phi] -> Fill(treeReader.JetInfo_leading_bjet_phi, isData ? 1. : NormalizationFactor);
        }
        if(!(treeReader.JetInfo_chosen_bjet_pt<-900.)){
            h[hist_chosen_bjet_pt] -> Fill(treeReader.JetInfo_chosen_bjet_pt, isData ? 1. : NormalizationFactor);
            h[hist_chosen_bjet_eta] -> Fill(treeReader.JetInfo_chosen_bjet_eta, isData ? 1. : NormalizationFactor);
            h[hist_chosen_bjet_phi] -> Fill(treeReader.JetInfo_chosen_bjet_phi, isData ? 1. : NormalizationFactor);
        }
        if(!(treeReader.JetInfo_jet1_pt<-900.)){
            h[hist_jet1_pt] -> Fill(treeReader.JetInfo_jet1_pt, isData ? 1. : NormalizationFactor);
            h[hist_jet1_eta] -> Fill(treeReader.JetInfo_jet1_eta, isData ? 1. : NormalizationFactor);
            h[hist_jet1_phi] -> Fill(treeReader.JetInfo_jet1_phi, isData ? 1. : NormalizationFactor);
        }
        if(!(treeReader.JetInfo_jet2_pt<-900.)){
            h[hist_jet2_pt] -> Fill(treeReader.JetInfo_jet2_pt, isData ? 1. : NormalizationFactor);
            h[hist_jet2_eta] -> Fill(treeReader.JetInfo_jet2_eta, isData ? 1. : NormalizationFactor);
            h[hist_jet2_phi] -> Fill(treeReader.JetInfo_jet2_phi, isData ? 1. : NormalizationFactor);
        }
        if(!(treeReader.JetInfo_dijet_delta_eta<-900.)){
            h[hist_dijet_eta] -> Fill(treeReader.JetInfo_dijet_delta_eta, isData ? 1. : NormalizationFactor);
            h[hist_dijet_phi] -> Fill(treeReader.JetInfo_dijet_delta_phi, isData ? 1. : NormalizationFactor);
            h[hist_dijet_angle] -> Fill(treeReader.JetInfo_dijet_delta_angle, isData ? 1. : NormalizationFactor);
        }

        if(!(treeReader.inv_mass_dijet<-900.))
            h[hist_inv_mass_dijet] -> Fill(treeReader.inv_mass_dijet, isData ? 1. : NormalizationFactor);
        if(!(treeReader.inv_mass_diphoton<-900.))
            h[hist_inv_mass_diphoton] -> Fill(treeReader.inv_mass_diphoton, isData ? 1. : NormalizationFactor);
        if(!(treeReader.inv_mass_tbw<-900.))
            h[hist_inv_mass_tbw] -> Fill(treeReader.inv_mass_tbw, isData ? 1. : NormalizationFactor);

        h[hist_DiPhoInfo_leadPt] -> Fill(treeReader.DiPhoInfo_leadPt, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_leadEta] -> Fill(treeReader.DiPhoInfo_leadEta, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_leadPhi] -> Fill(treeReader.DiPhoInfo_leadPhi, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_leadE] -> Fill(treeReader.DiPhoInfo_leadE, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_leadIDMVA] -> Fill(treeReader.DiPhoInfo_leadIDMVA, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_subleadPt] -> Fill(treeReader.DiPhoInfo_subleadPt, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_subleadEta] -> Fill(treeReader.DiPhoInfo_subleadEta, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_subleadPhi] -> Fill(treeReader.DiPhoInfo_subleadPhi, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_subleadE] -> Fill(treeReader.DiPhoInfo_subleadE, isData ? 1. : NormalizationFactor);
        h[hist_DiPhoInfo_subleadIDMVA] -> Fill(treeReader.DiPhoInfo_subleadIDMVA, isData ? 1. : NormalizationFactor);
        //hist_DiPhoInfo_leadIDMVA_ori -> Fill(treeReader.DiPhoInfo_leadIDMVA_ori, isData ? 1. : NormalizationFactor);
        //hist_DiPhoInfo_subleadIDMVA_ori -> Fill(treeReader.DiPhoInfo_subleadIDMVA_ori, isData ? 1. : NormalizationFactor);
        //hist_inv_mass_diphoton_ori -> Fill(treeReader.inv_mass_diphoton_ori, isData ? 1. : NormalizationFactor);
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 700, 800);
    for(int i=0; i<totalHistNum; i++){
        MakePlots(c1, h[i], Form("%s/%s.png", output_dir, histNames[i].c_str()));
    }

    fout -> Close();
    fin  -> Close();
}

int main(int argc, char *argv[]){
    char input_file[512]; sprintf(input_file, "%s", argv[1]); printf("[INFO] input_file  = %s\n", input_file);
    char output_file[512]; sprintf(output_file, "%s", argv[2]); printf("[INFO] output_file = %s\n", output_file);
    char dataset[512]; sprintf(dataset, "%s", argv[3]); printf("[INFO] dataset     = %s\n", dataset);
    char output_dir[512];  sprintf(output_dir, "%s", argv[4]); printf("[INFO] output_dir = %s\n", output_dir);
    Selection(input_file, output_file, dataset, output_dir);
    return 1;
}


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
    mytree -> SetBranchAddress("JetInfo_leading_bjet_pt", &JetInfo_leading_bjet_pt);
    mytree -> SetBranchAddress("JetInfo_leading_bjet_eta", &JetInfo_leading_bjet_eta);
    mytree -> SetBranchAddress("JetInfo_leading_bjet_phi", &JetInfo_leading_bjet_phi);
    mytree -> SetBranchAddress("JetInfo_chosen_bjet_pt", &JetInfo_chosen_bjet_pt);
    mytree -> SetBranchAddress("JetInfo_chosen_bjet_eta", &JetInfo_chosen_bjet_eta);
    mytree -> SetBranchAddress("JetInfo_chosen_bjet_phi", &JetInfo_chosen_bjet_phi);
    mytree -> SetBranchAddress("JetInfo_jet1_pt", &JetInfo_jet1_pt);
    mytree -> SetBranchAddress("JetInfo_jet1_eta", &JetInfo_jet1_eta);
    mytree -> SetBranchAddress("JetInfo_jet1_phi", &JetInfo_jet1_phi);
    mytree -> SetBranchAddress("JetInfo_jet2_pt", &JetInfo_jet2_pt);
    mytree -> SetBranchAddress("JetInfo_jet2_eta", &JetInfo_jet2_eta);
    mytree -> SetBranchAddress("JetInfo_jet2_phi", &JetInfo_jet2_phi);
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
