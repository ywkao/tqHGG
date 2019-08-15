//***************************************************************************
//
// FileName    : selection.cpp
// Purpose     : Develop for top FCNH with H to two photons analysis
// Description : Applying event selection & Preparing histograms for individual dataset.
// Deetail     : PU reweighting, leptonic/hadronic channels, hists, (top reconstruction).
// Author      : Yu-Wei Kao [ykao@cern.ch]
//
//***************************************************************************
#include <stdio.h>
#include <TFile.h>
#include <TH1D.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include "../include/main.h"
#include "../include/selection.h"
#include "../include/enumhist.h"
using namespace std;

bool bool_isHadronic;

void Selection(char* input_file, char* output_file, char* dataset, char* output_dir, char* channel){
    bool isData = isThisDataOrNot(dataset);
    //-----
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
    //==================================================//
    //-------------------- Event Loop ------------------//
    //==================================================//
    int nentries = treeReader.GetEntries();
    for(int ientry=0; ientry<nentries; ientry++){
        TTree* tmp = treeReader.GetTTree();
        tmp->GetEntry(ientry);
        if(ientry==0) printf("[INFO] N_entries = %d/%d\n", nentries, treeReader.EvtInfo_totalEntry_before_preselection);
        if((ientry+1)%10000==0 || (ientry+1)==nentries) printf("ientry = %d\r", ientry);
        //========= PU Reweighting =========//
        double PU_reweighting_factor = h_pu_reweighting_factor->GetBinContent((int)treeReader.EvtInfo_NPu+1);
        //double PU_reweighting_factor = 1.; //No PU
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
            if(treeReader.num_jets<3) passEvent = false;
        }
        if(!passEvent) continue;
        //--------- check bjets ---------//
        int num_bjets = 0;
        if(!(treeReader.num_jets<1)){
            for(int i=0; i<treeReader.num_jets; ++i){
                if(treeReader.JetInfo_jet_pfDeepCSVJetTags_probb_selection->at(i)+treeReader.JetInfo_jet_pfDeepCSVJetTags_probbb_selection->at(i) >= pfDeepCSVJetTags_tight){
                    num_bjets += 1;
                    if(num_bjets == 1){
                        h[hist_leading_bjet_pt] -> Fill(treeReader.JetInfo_jet_pt_selection->at(i), isData ? 1. : NormalizationFactor);
                        h[hist_leading_bjet_eta] -> Fill(treeReader.JetInfo_jet_eta_selection->at(i), isData ? 1. : NormalizationFactor);
                        h[hist_leading_bjet_phi] -> Fill(treeReader.JetInfo_jet_phi_selection->at(i), isData ? 1. : NormalizationFactor);
                        h[hist_leading_bjet_energy] -> Fill(treeReader.JetInfo_jet_energy_selection->at(i), isData ? 1. : NormalizationFactor);
                    }
                }
            }//end of looping jets
        }
        if(num_bjets == 0) continue;
        //========= Store Info =========//
        h[hist_EvtInfo_NPu] -> Fill(treeReader.EvtInfo_NPu, isData ? 1. : NormalizationFactor);
        h[hist_EvtInfo_Rho] -> Fill(treeReader.EvtInfo_Rho, isData ? 1. : NormalizationFactor);
        h[hist_EvtInfo_NVtx] -> Fill(treeReader.EvtInfo_NVtx, isData ? 1. : NormalizationFactor);
        h[hist_EvtInfo_NVtx_wopu] -> Fill(treeReader.EvtInfo_NVtx, isData ? 1. : NormalizationFactor_wopu);
        h[hist_EvtInfo_genweight] -> Fill(treeReader.EvtInfo_genweight, isData ? 1. : NormalizationFactor);
        //--------- Diphoton ---------//
        TLorentzVector Diphoton;
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
        Diphoton.SetPtEtaPhiE(treeReader.DiPhoInfo_pt, treeReader.DiPhoInfo_eta, treeReader.DiPhoInfo_phi, treeReader.DiPhoInfo_energy);
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
        //--------- Jets ---------//
        h[hist_jets_size] -> Fill(treeReader.jets_size, isData ? 1. : NormalizationFactor);
        h[hist_num_jets] -> Fill(treeReader.num_jets, isData ? 1. : NormalizationFactor);
        if(treeReader.num_jets>0){
            for(int i=0; i<treeReader.num_jets; ++i){
                h[hist_JetInfo_jet_pt] -> Fill(treeReader.JetInfo_jet_pt_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_JetInfo_jet_eta] -> Fill(treeReader.JetInfo_jet_eta_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_JetInfo_jet_phi] -> Fill(treeReader.JetInfo_jet_phi_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_JetInfo_jet_energy] -> Fill(treeReader.JetInfo_jet_energy_selection->at(i), isData ? 1. : NormalizationFactor);
                h[hist_JetInfo_jet_diphoton_deltaR] -> Fill(treeReader.JetInfo_jet_diphoton_deltaR_selection->at(i), isData ? 1. : NormalizationFactor);
                if(i==0){//leading jet
                    h[hist_jet1_pt] -> Fill(treeReader.JetInfo_jet_pt_selection->at(i), isData ? 1. : NormalizationFactor);
                    h[hist_jet1_eta] -> Fill(treeReader.JetInfo_jet_eta_selection->at(i), isData ? 1. : NormalizationFactor);
                    h[hist_jet1_phi] -> Fill(treeReader.JetInfo_jet_phi_selection->at(i), isData ? 1. : NormalizationFactor);
                    h[hist_jet1_energy] -> Fill(treeReader.JetInfo_jet_energy_selection->at(i), isData ? 1. : NormalizationFactor);
                    h[hist_jet1_diphoton_deltaR] -> Fill(treeReader.JetInfo_jet_diphoton_deltaR_selection->at(i), isData ? 1. : NormalizationFactor);
                    TLorentzVector jet; jet.SetPtEtaPhiE(treeReader.JetInfo_jet_pt_selection->at(i), treeReader.JetInfo_jet_eta_selection->at(i), treeReader.JetInfo_jet_phi_selection->at(i), treeReader.JetInfo_jet_energy_selection->at(i));
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
                    h[hist_jet2_diphoton_deltaR] -> Fill(treeReader.JetInfo_jet_diphoton_deltaR_selection->at(i), isData ? 1. : NormalizationFactor);
                    TLorentzVector jet; jet.SetPtEtaPhiE(treeReader.JetInfo_jet_pt_selection->at(i), treeReader.JetInfo_jet_eta_selection->at(i), treeReader.JetInfo_jet_phi_selection->at(i), treeReader.JetInfo_jet_energy_selection->at(i));
                    if(treeReader.num_leptons>0){
                        for(int i=0; i<treeReader.num_leptons; ++i){
                            double delta_R = jet.DeltaR(Leptons.at(i));
                            h[hist_jet2_lepton_deltaR] -> Fill(delta_R, isData ? 1. : NormalizationFactor);
                        }
                    }
                }
            }//end of looping jets
        }
        
        //------------------------
        // Reconstruction tbw
        //------------------------
        //hist_chosen_bjet_pt;
        //hist_chosen_bjet_eta;
        //hist_chosen_bjet_phi;
        ////------------------------
        //hist_inv_mass_dijet;
        //hist_inv_mass_diphoton;
        //hist_inv_mass_tbw;
    }//end of event loop
    double yields = h[hist_EvtInfo_NPu] -> Integral();
    printf("[INFO] Expected yields");
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
}

int main(int argc, char *argv[]){
    char input_file[512]; sprintf(input_file, "%s", argv[1]); printf("[INFO] input_file  = %s\n", input_file);
    char output_file[512]; sprintf(output_file, "%s", argv[2]); printf("[INFO] output_file = %s\n", output_file);
    char dataset[512]; sprintf(dataset, "%s", argv[3]); printf("[INFO] dataset     = %s\n", dataset);
    char output_dir[512];  sprintf(output_dir, "%s", argv[4]); printf("[INFO] output_dir  = %s\n", output_dir);
    char channel[512];  sprintf(channel, "%s", argv[5]); printf("[INFO] channel     = %s\n", channel);
    Selection(input_file, output_file, dataset, output_dir, channel);
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
    mytree -> SetBranchAddress("jets_size", &jets_size);
    mytree -> SetBranchAddress("num_jets", &num_jets);
    mytree -> SetBranchAddress("JetInfo_jet_pt", &JetInfo_jet_pt_selection);
    mytree -> SetBranchAddress("JetInfo_jet_eta", &JetInfo_jet_eta_selection);
    mytree -> SetBranchAddress("JetInfo_jet_phi", &JetInfo_jet_phi_selection);
    mytree -> SetBranchAddress("JetInfo_jet_energy", &JetInfo_jet_energy_selection);
    mytree -> SetBranchAddress("JetInfo_jet_diphoton_deltaR", &JetInfo_jet_diphoton_deltaR_selection);
    mytree -> SetBranchAddress("JetInfo_jet_pfDeepCSVJetTags_probb", &JetInfo_jet_pfDeepCSVJetTags_probb_selection);
    mytree -> SetBranchAddress("JetInfo_jet_pfDeepCSVJetTags_probbb", &JetInfo_jet_pfDeepCSVJetTags_probbb_selection);
    printf("[INFO] myTreeClass::SetBranchAddresses : Finished!\n");
}
myParameters::myParameters(){
    JetInfo_jet_pt_selection = new std::vector<float>;
    JetInfo_jet_eta_selection = new std::vector<float>;
    JetInfo_jet_phi_selection = new std::vector<float>;
    JetInfo_jet_energy_selection = new std::vector<float>;
    JetInfo_jet_diphoton_deltaR_selection = new std::vector<float>;
    JetInfo_jet_pfDeepCSVJetTags_probb_selection = new std::vector<float>;
    JetInfo_jet_pfDeepCSVJetTags_probbb_selection = new std::vector<float>;
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
    delete JetInfo_jet_pt_selection;
    delete JetInfo_jet_eta_selection;
    delete JetInfo_jet_phi_selection;
    delete JetInfo_jet_energy_selection;
    delete JetInfo_jet_diphoton_deltaR_selection;
    delete JetInfo_jet_pfDeepCSVJetTags_probb_selection;
    delete JetInfo_jet_pfDeepCSVJetTags_probbb_selection;
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
