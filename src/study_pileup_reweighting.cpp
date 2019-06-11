//**************************************************
//
// Study of pile up reweighting.
// Compare shape of NVtx and Rho between data & MC
// Input: ntuples from You-Ying
// Output: root file with histograms (NPu, NVtx, Rho)
//
//**************************************************
#include <stdio.h>
#include <math.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
using namespace std;

bool isThisDataOrNot(char* dataset);

int main(int argc, char *argv[]){
    //============================//
    //----- Input parameters -----//
    //============================//
    char input_file[256]; sprintf(input_file, "%s", argv[1]); printf("[INFO] input_file  = %s\n", input_file);
    char output_file[256]; sprintf(output_file, "%s", argv[2]); printf("[INFO] output_file = %s\n", output_file);
    char dataset[256]; sprintf(dataset, "%s", argv[3]); printf("[INFO] dataset     = %s\n", dataset);
    bool isData = isThisDataOrNot(dataset);//Determine normalization factor
    //==============================//
    //----- Read input file(s) -----//
    //==============================//
    TChain *treeReader;
    treeReader = new TChain("flashggNtuples/flashggStdTree");
    treeReader->Add(Form("%s/*.root", input_file));
    //--------------------------------------------------
    TFile *f_pu = TFile::Open("data/MCPileUp.root");
    TH1D  *h_pu = (TH1D*)f_pu->Get("puhist");
    //--------------------------------------------------
    Int_t EvtInfo_NVtx;
    float EvtInfo_NPu;
    float EvtInfo_Rho;
    float EvtInfo_genweight;
    treeReader->SetBranchAddress("EvtInfo.NPu", &EvtInfo_NPu);
    treeReader->SetBranchAddress("EvtInfo.Rho", &EvtInfo_Rho);
    treeReader->SetBranchAddress("EvtInfo.NVtx", &EvtInfo_NVtx);
    treeReader->SetBranchAddress("EvtInfo.genweight", &EvtInfo_genweight);
    //===============================//
    //----- Prepare output file -----//
    //===============================//
    TFile *fout = new TFile(output_file, "RECREATE");
    TH1D *hist_NPu  = new TH1D("hist_NPu", ";;", 100, 0, 100);          hist_NPu->Sumw2();
    TH1D *hist_Rho  = new TH1D("hist_Rho", ";;", 100, 0, 100);          hist_Rho->Sumw2();
    TH1D *hist_NVtx = new TH1D("hist_NVtx", ";;", 100, 0, 100);         hist_NVtx->Sumw2();
    TH1D *hist_Rho_pu  = new TH1D("hist_Rho_pu", ";;", 100, 0, 100);    hist_Rho_pu->Sumw2();
    TH1D *hist_NVtx_pu = new TH1D("hist_NVtx_pu", ";;", 100, 0, 100);   hist_NVtx_pu->Sumw2();
    TH1D *hist_Rho_pu_gen  = new TH1D("hist_Rho_pu_gen", ";;", 100, 0, 100);    hist_Rho_pu_gen->Sumw2();
    TH1D *hist_NVtx_pu_gen = new TH1D("hist_NVtx_pu_gen", ";;", 100, 0, 100);   hist_NVtx_pu_gen->Sumw2();
    //===============================//
    //--------  Event loop  ---------//
    //===============================//
    // Goal: Nvtx, Rho
    int nentries = treeReader->GetEntries(); printf("[INFO] N_entries = %d\n", nentries);
    for(int ientry=0; ientry<nentries; ientry++){
        treeReader->GetEntry(ientry);//load data
        if((ientry+1)%1000==0 || (ientry+1)==nentries) printf("ientry = %d\r", ientry);
        double PU_reweighting_factor = isData ? 1. : h_pu->GetBinContent((int)EvtInfo_NPu+1);
        double genweight_factor = isData ? 1. : EvtInfo_genweight;
        double NormalizationFactor = PU_reweighting_factor;
        double NormalizationFactor_gen = PU_reweighting_factor * genweight_factor;//consider event genweight
        double value_NPu  = EvtInfo_NPu;
        double value_Rho  = EvtInfo_Rho;
        double value_NVtx = EvtInfo_NVtx;
        hist_NPu->Fill(value_NPu, 1.);
        hist_Rho->Fill(value_Rho, 1.);
        hist_NVtx->Fill(value_NVtx, 1.);
        hist_Rho_pu->Fill(value_Rho, NormalizationFactor);
        hist_NVtx_pu->Fill(value_NVtx, NormalizationFactor);
        hist_Rho_pu_gen->Fill(value_Rho, NormalizationFactor_gen);
        hist_NVtx_pu_gen->Fill(value_NVtx, NormalizationFactor_gen);
    }// End of event loop.

    fout->Write();
    fout->Close();
    return 1;
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
