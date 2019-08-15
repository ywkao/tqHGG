//***************************************************************************
//
// FileName    : preselection.cpp
// Purpose     : Develop for top FCNH with H to two photons analysis
// Description : Extracting event info & Selecting objects (diphoton, leptons, jets).
// Author      : Yu-Wei Kao [ykao@cern.ch]
//
//***************************************************************************
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
#include "../include/cross_section_v2.h"
using namespace std;

int main(int argc, char *argv[]){
    //============================//
    //----- Input file names -----//
    //============================//
    char input_file[512]; sprintf(input_file, "%s", argv[1]); printf("[INFO] input_file  = %s\n", input_file);
    char output_file[512]; sprintf(output_file, "%s", argv[2]); printf("[INFO] output_file = %s\n", output_file);
    char dataset[512]; sprintf(dataset, "%s", argv[3]); printf("[INFO] dataset     = %s\n", dataset);
    bool isData = isThisDataOrNot(dataset);
    bool isMCsignal = isThisMCsignal(dataset);
    bool isMultiFile = isThisMultiFile(dataset);
    TFile *fout = new TFile(output_file, "RECREATE");
    //==============================//
    //----- Read input file(s) -----//
    //==============================//
    flashggStdTreeReader treeReader;
    treeReader.InitChain("flashggNtuples/flashggStdTree");
    if(!isMultiFile) treeReader.AddSingleRootFile(input_file);
    else             treeReader.AddMultiRootFile(input_file);
    treeReader.SetBranchAddresses();
    //===============================//
    //----- Prepare output file -----//
    //===============================//
    myTreeClass mytree;
    mytree.InitTree();
    mytree.MakeNewBranchAddresses();
    //==================================================//
    //--------   Event Loop [Normalization]   ----------//
    //==================================================//
    // Note: Normalization factors are not used during preselection.
    // It is stored and would be used during selection.
    // Question: keep total genweight before or after preselection?
    int nentries = treeReader.GetEntries(); printf("[INFO] N_entries = %d\n", nentries);
    double NormalizationFactor;
    double Luminosity = 41.53; //fb{-1}
    double CrossSection = GetXsec(dataset); //pb
    double BranchingFraction = GetBranchingFraction(dataset); //pb
    printf("[INFO] CrossSection = %f !\n", CrossSection);
    printf("[INFO] Equivalent lumi. = %f !\n", (double)nentries/CrossSection);
    printf("[INFO] BranchingFraction = %f !\n", BranchingFraction);
    double TotalGenweight=0;
    for(int ientry=0; ientry<nentries; ientry++){
        treeReader.flashggStdTree->GetEntry(ientry);//load data
        TotalGenweight+=treeReader.GetGenWeight();
    }
    NormalizationFactor = 1000. * Luminosity * CrossSection * BranchingFraction / TotalGenweight;
    printf("[INFO] TotalGenweight = %f!\n", TotalGenweight);
    printf("[INFO] NormalizationFactor = %f!\n", isData ? 1. : NormalizationFactor);


    //##################################################//
    //#########    Event Loop [Selection]    ###########//
    //##################################################//
    // Goal: leptons, jets, diphoton; t->b+W(jj), 1 bjet + 2 chi2 jets 
    int Nevents_pass_selection = 0;
    printf("[CHECK-0] Nevents_pass_selection = %d\n", Nevents_pass_selection);
    for(int ientry=0; ientry<nentries; ientry++){
        treeReader.flashggStdTree->GetEntry(ientry);//load data
        if((ientry+1)%1000==0 || (ientry+1)==nentries) printf("ientry = %d\r", ientry);
        //==================================================//
        //-------------   Reset Parameters   ---------------//
        //==================================================//
        mytree.Clear();
        //==================================================//
        //--------------   Basic Selectoin   ---------------//
        //==================================================//
        //require MC events pass trigger.
        if(!treeReader.EvtInfo_passTrigger) continue;
        //control region
        if(treeReader.DiPhoInfo_mass<100 || treeReader.DiPhoInfo_mass>180) continue;
        if( !isMCsignal && treeReader.DiPhoInfo_mass>120 && treeReader.DiPhoInfo_mass<130) continue;
        //require the quality of photons. (PT requirement)
        bool pass_leadingPhotonPT = treeReader.DiPhoInfo_leadPt > treeReader.DiPhoInfo_mass / 2.;
        bool pass_subleadingPhotonPT = treeReader.DiPhoInfo_subleadPt > treeReader.DiPhoInfo_mass / 4.;
        bool pass_photon_criteria = pass_leadingPhotonPT && pass_subleadingPhotonPT;
        if(!pass_photon_criteria) continue;
        //Others
        //if(treeReader.DiPhoInfo_mass<0) continue;
        //if( !(treeReader.jets_size>0) ) continue;
        //if(!(DiPhoInfo_leadIDMVA>0)) continue;
        //==================================================//
        //------------   Normalization factor   ------------//
        //==================================================//
        mytree.EvtInfo_totalEntry_before_preselection = nentries;
        mytree.EvtInfo_NormalizationFactor_lumi = isData ? 1. : NormalizationFactor;
        //===============================================//
        //-----------   Store Diphoton Info   -----------//
        //===============================================//
        mytree.DiPhoInfo_leadPt = treeReader.DiPhoInfo_leadPt;
        mytree.DiPhoInfo_leadEta = treeReader.DiPhoInfo_leadEta;
        mytree.DiPhoInfo_leadPhi = treeReader.DiPhoInfo_leadPhi;
        mytree.DiPhoInfo_leadE = treeReader.DiPhoInfo_leadE;
        mytree.DiPhoInfo_leadhoe = treeReader.DiPhoInfo_leadhoe;
        mytree.DiPhoInfo_leadIDMVA = treeReader.DiPhoInfo_leadIDMVA;
        mytree.DiPhoInfo_subleadPt = treeReader.DiPhoInfo_subleadPt;
        mytree.DiPhoInfo_subleadEta = treeReader.DiPhoInfo_subleadEta;
        mytree.DiPhoInfo_subleadPhi = treeReader.DiPhoInfo_subleadPhi;
        mytree.DiPhoInfo_subleadE = treeReader.DiPhoInfo_subleadE;
        mytree.DiPhoInfo_subleadhoe = treeReader.DiPhoInfo_subleadhoe;
        mytree.DiPhoInfo_subleadIDMVA = treeReader.DiPhoInfo_subleadIDMVA;
        mytree.DiPhoInfo_mass = treeReader.DiPhoInfo_mass;
        mytree.DiPhoInfo_pt = treeReader.DiPhoInfo_pt;

        TLorentzVector leading_photon, subleading_photon, diphoton;
        leading_photon.SetPtEtaPhiE(treeReader.DiPhoInfo_leadPt, treeReader.DiPhoInfo_leadEta, treeReader.DiPhoInfo_leadPhi, treeReader.DiPhoInfo_leadE);
        subleading_photon.SetPtEtaPhiE(treeReader.DiPhoInfo_subleadPt, treeReader.DiPhoInfo_subleadEta, treeReader.DiPhoInfo_subleadPhi, treeReader.DiPhoInfo_subleadE);
        diphoton = leading_photon + subleading_photon;
        mytree.DiPhoInfo_eta = diphoton.Eta();
        mytree.DiPhoInfo_phi = diphoton.Phi();
        mytree.DiPhoInfo_energy = diphoton.Energy();
        //==============================//
        //-----  Select Electrons  -----//
        //==============================//
        mytree.ElecInfo_Size = treeReader.ElecInfo_Size;
        std::vector<TLorentzVector> Electrons;
        bool bool_AtLeastOneElectron = treeReader.ElecInfo_Size>0 ? true : false;//treeReader.ElecInfo_Size = -999 => event without diphoton candidate
        if(bool_AtLeastOneElectron){
            for(int i=0; i<treeReader.ElecInfo_Size; i++){
                if( !treeReader.ElecInfo_EGMCutBasedIDMedium->at(i) ) continue;
                if( fabs(treeReader.ElecInfo_Eta->at(i)) > 2.4 ) continue;
                if( fabs(treeReader.ElecInfo_Eta->at(i)) > 1.4442 && fabs(treeReader.ElecInfo_Eta->at(i)) < 1.566 ) continue;
                if( fabs(treeReader.ElecInfo_Pt->at(i))  < 20  ) continue;
                //--- check deltaR(electron,photon) ---//
                TLorentzVector electron; 
                electron.SetPtEtaPhiE(treeReader.ElecInfo_Pt->at(i), treeReader.ElecInfo_Eta->at(i), treeReader.ElecInfo_Phi->at(i), treeReader.ElecInfo_Energy->at(i));
                double delta_R = electron.DeltaR(leading_photon);
                if( delta_R<0.3 ) continue;
                delta_R = electron.DeltaR(subleading_photon);
                if( delta_R<0.3 ) continue;
                //--- calculate deltaR(electron,photon/diphoton) and store info ---//
                delta_R = electron.DeltaR(diphoton);          mytree.ElecInfo_electron_diphoton_deltaR.push_back(delta_R);
                delta_R = electron.DeltaR(leading_photon);    mytree.ElecInfo_electron_leadingPhoton_deltaR.push_back(delta_R);
                delta_R = electron.DeltaR(subleading_photon); mytree.ElecInfo_electron_subleadingPhoton_deltaR.push_back(delta_R);
                //--- store information ---//
                mytree.ElecInfo_electron_pt.push_back(treeReader.ElecInfo_Pt->at(i));
                mytree.ElecInfo_electron_eta.push_back(treeReader.ElecInfo_Eta->at(i));
                mytree.ElecInfo_electron_phi.push_back(treeReader.ElecInfo_Phi->at(i));
                mytree.ElecInfo_electron_energy.push_back(treeReader.ElecInfo_Phi->at(i));
                mytree.num_electrons+=1;
                Electrons.push_back(electron);
            }
        }
        else{
                mytree.num_electrons=treeReader.ElecInfo_Size;
        }
        bool bool_AtLeastOneSelectedElectron = mytree.num_electrons>0 ? true : false;//for calculation of deltaR(e,j).
        //==========================//
        //-----  Select Muons  -----//
        //==========================//
        mytree.MuonInfo_Size = treeReader.MuonInfo_Size;
        std::vector<TLorentzVector> Muons;
        bool bool_AtLeastOneMuon = treeReader.MuonInfo_Size>0 ? true : false;//treeReader.MuonInfo_Size = -999 => event without diphoton candidate
        if(bool_AtLeastOneMuon){
            for(int i=0; i<treeReader.MuonInfo_Size; i++){
                if( !treeReader.MuonInfo_CutBasedIdTight->at(i) ) continue;
                if( fabs(treeReader.MuonInfo_Eta->at(i)) > 2.4 ) continue;
                if( fabs(treeReader.MuonInfo_Eta->at(i)) > 1.4442 && fabs(treeReader.MuonInfo_Eta->at(i)) < 1.566 ) continue;
                if( fabs(treeReader.MuonInfo_Pt->at(i))  < 20  ) continue;
                if( treeReader.MuonInfo_PFIsoDeltaBetaCorrR04->at(i) > 0.25  ) continue;
                //--- check deltaR(muon,photon) ---//
                TLorentzVector muon; 
                muon.SetPtEtaPhiE(treeReader.MuonInfo_Pt->at(i), treeReader.MuonInfo_Eta->at(i), treeReader.MuonInfo_Phi->at(i), treeReader.MuonInfo_Energy->at(i));
                double delta_R = muon.DeltaR(leading_photon);
                if( delta_R<0.3 ) continue;
                delta_R = muon.DeltaR(subleading_photon);
                if( delta_R<0.3 ) continue;
                //--- calculate deltaR(muon,diphoton) and store info ---//
                delta_R = muon.DeltaR(diphoton);          mytree.MuonInfo_muon_diphoton_deltaR.push_back(delta_R);
                delta_R = muon.DeltaR(leading_photon);    mytree.MuonInfo_muon_leadingPhoton_deltaR.push_back(delta_R);
                delta_R = muon.DeltaR(subleading_photon); mytree.MuonInfo_muon_subleadingPhoton_deltaR.push_back(delta_R);
                //--- store information ---//
                mytree.MuonInfo_muon_pt.push_back(treeReader.MuonInfo_Pt->at(i));
                mytree.MuonInfo_muon_eta.push_back(treeReader.MuonInfo_Eta->at(i));
                mytree.MuonInfo_muon_phi.push_back(treeReader.MuonInfo_Phi->at(i));
                mytree.MuonInfo_muon_energy.push_back(treeReader.MuonInfo_Phi->at(i));
                
                mytree.num_muons+=1;
                Muons.push_back(muon);
            }
        }
        else{
                mytree.num_muons=treeReader.MuonInfo_Size;
        }
        bool bool_AtLeastOneSelectedMuon = mytree.num_muons>0 ? true : false;//for calculation of deltaR(mu,j).
        mytree.num_leptons = mytree.num_electrons + mytree.num_muons;
        //printf("num_leptons = %d\n", mytree.num_leptons);
        //=========================//
        //-----  Select Jets  -----//
        //=========================//
        mytree.jets_size = treeReader.jets_size;
        std::vector<TLorentzVector> Jets;
        bool bool_AtLeastOneJet = treeReader.jets_size>0 ? true : false;//treeReader.jets_size = -999 => event without diphoton candidate
        if(bool_AtLeastOneJet){
            for(int i=0; i<treeReader.jets_size; i++){
                //flashgg::Tight2017 jet ID #Already applied flashgg package.
                if( fabs(treeReader.JetInfo_Eta->at(i)) > 2.4 ) continue;
                if( fabs(treeReader.JetInfo_Pt->at(i))  < 25  ) continue;
                //--- check deltaR(jet,photon) ---//
                TLorentzVector jet; 
                jet.SetPtEtaPhiE(treeReader.JetInfo_Pt->at(i), treeReader.JetInfo_Eta->at(i), treeReader.JetInfo_Phi->at(i), treeReader.JetInfo_Energy->at(i));
                double delta_R = jet.DeltaR(leading_photon);
                if( delta_R<0.4 ) continue;
                delta_R = jet.DeltaR(subleading_photon);
                if( delta_R<0.4 ) continue;
                bool bool_passJetLeptonSeparation = true;//if no leptons selected, the jet pass the delta_R criterion automatically.
                //--- check deltaR(jet,e) ---//
                if(bool_AtLeastOneSelectedElectron){
                    for(int i=0; i<mytree.num_electrons; i++){
                        delta_R = jet.DeltaR(Electrons.at(i));
                        if( delta_R<0.4 ) bool_passJetLeptonSeparation = false;
                    }
                }
                if(!bool_passJetLeptonSeparation) continue;//if not pass, reject the jet.
                //--- check deltaR(jet,mu) ---//
                if(bool_AtLeastOneSelectedMuon){
                    for(int i=0; i<mytree.num_muons; i++){
                        delta_R = jet.DeltaR(Muons.at(i));
                        if( delta_R<0.4 ) bool_passJetLeptonSeparation = false;
                    }
                }
                if(!bool_passJetLeptonSeparation) continue;//if not pass, reject the jet.
                //--- calculate deltaR(jet,diphoton) and store info ---//
                delta_R = jet.DeltaR(diphoton);          mytree.JetInfo_jet_diphoton_deltaR.push_back(delta_R);
                delta_R = jet.DeltaR(leading_photon);    mytree.JetInfo_jet_leadingPhoton_deltaR.push_back(delta_R);
                delta_R = jet.DeltaR(subleading_photon); mytree.JetInfo_jet_subleadingPhoton_deltaR.push_back(delta_R);
                //--- store information ---//
                mytree.JetInfo_jet_pt.push_back(treeReader.JetInfo_Pt->at(i));
                mytree.JetInfo_jet_eta.push_back(treeReader.JetInfo_Eta->at(i));
                mytree.JetInfo_jet_phi.push_back(treeReader.JetInfo_Phi->at(i));
                mytree.JetInfo_jet_energy.push_back(treeReader.JetInfo_Phi->at(i));
                
                mytree.JetInfo_jet_pfDeepCSVJetTags_probb.push_back(treeReader.JetInfo_pfDeepCSVJetTags_probb->at(i));
                mytree.JetInfo_jet_pfDeepCSVJetTags_probbb.push_back(treeReader.JetInfo_pfDeepCSVJetTags_probbb->at(i));
                mytree.num_jets+=1;
                Jets.push_back(jet);
            }
        }
        else{
                mytree.num_jets=treeReader.jets_size;
        }
        //================================================//
        //-----------   Store EventPar Info    -----------//
        //================================================//
        mytree.EvtInfo_NPu = treeReader.EvtInfo_NPu;
        mytree.EvtInfo_Rho = treeReader.EvtInfo_Rho;
        mytree.EvtInfo_NVtx = treeReader.EvtInfo_NVtx;
        mytree.EvtInfo_genweight = treeReader.EvtInfo_genweight;
        //==================================================//
        //-------------   Event Counting     ---------------//
        //==================================================//
        Nevents_pass_selection += 1;
        mytree.Fill();
        if(ientry == nentries - 1) printf("[CHECK-1] Nevents_pass_selection = %d\n", Nevents_pass_selection);
    }// End of event loop.
    //==================================================//
    //---------------------  Report  -------------------//
    //==================================================//
    printf("[CHECK-2] Nevents_pass_selection = %d\n", Nevents_pass_selection);
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
    if((string)dataset == "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8") return true;
    if((string)dataset == "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8") return true;
    if((string)dataset == "WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8") return true;
    if((string)dataset == "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8") return true;
    if((string)dataset == "WW_TuneCP5_13TeV-pythia8") return true;
    if((string)dataset == "WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8") return true;
    if((string)dataset == "ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8") return true;
    if((string)dataset == "ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8") return true;
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
    flashggStdTree->SetBranchAddress("EvtInfo.NVtx", &EvtInfo_NVtx);
    flashggStdTree->SetBranchAddress("EvtInfo.Rho", &EvtInfo_Rho);
    flashggStdTree->SetBranchAddress("EvtInfo.genweight", &EvtInfo_genweight);
    //------------------------
    flashggStdTree->SetBranchAddress("DiPhoInfo.mass", &DiPhoInfo_mass);
    flashggStdTree->SetBranchAddress("DiPhoInfo.pt", &DiPhoInfo_pt);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadPt", &DiPhoInfo_leadPt);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadEta", &DiPhoInfo_leadEta);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadPhi", &DiPhoInfo_leadPhi);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadE", &DiPhoInfo_leadE);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadhoe", &DiPhoInfo_leadhoe);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadIDMVA", &DiPhoInfo_leadIDMVA);
    //------------------------
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadPt", &DiPhoInfo_subleadPt);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadEta", &DiPhoInfo_subleadEta);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadPhi", &DiPhoInfo_subleadPhi);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadE", &DiPhoInfo_subleadE);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadhoe", &DiPhoInfo_subleadhoe);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadIDMVA", &DiPhoInfo_subleadIDMVA);
    //------------------------
    flashggStdTree->SetBranchAddress("jets_size", &jets_size);
    flashggStdTree->SetBranchAddress("JetInfo.Pt", &JetInfo_Pt);
    flashggStdTree->SetBranchAddress("JetInfo.Eta", &JetInfo_Eta);
    flashggStdTree->SetBranchAddress("JetInfo.Phi", &JetInfo_Phi);
    flashggStdTree->SetBranchAddress("JetInfo.Mass", &JetInfo_Mass);
    flashggStdTree->SetBranchAddress("JetInfo.Energy", &JetInfo_Energy);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probb", &JetInfo_pfDeepCSVJetTags_probb);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probbb", &JetInfo_pfDeepCSVJetTags_probbb);
    //------------------------
    flashggStdTree->SetBranchAddress("ElecInfo.Size", &ElecInfo_Size);
    flashggStdTree->SetBranchAddress("ElecInfo.Charge", &ElecInfo_Charge);
    flashggStdTree->SetBranchAddress("ElecInfo.Pt", &ElecInfo_Pt);
    flashggStdTree->SetBranchAddress("ElecInfo.Eta", &ElecInfo_Eta);
    flashggStdTree->SetBranchAddress("ElecInfo.Phi", &ElecInfo_Phi);
    flashggStdTree->SetBranchAddress("ElecInfo.Energy", &ElecInfo_Energy);
    flashggStdTree->SetBranchAddress("ElecInfo.EtaSC", &ElecInfo_EtaSC);
    flashggStdTree->SetBranchAddress("ElecInfo.PhiSC", &ElecInfo_PhiSC);
    flashggStdTree->SetBranchAddress("ElecInfo.GsfTrackDz", &ElecInfo_GsfTrackDz);
    flashggStdTree->SetBranchAddress("ElecInfo.GsfTrackDxy", &ElecInfo_GsfTrackDxy);
    flashggStdTree->SetBranchAddress("ElecInfo.EGMCutBasedIDVeto", &ElecInfo_EGMCutBasedIDVeto);
    flashggStdTree->SetBranchAddress("ElecInfo.EGMCutBasedIDLoose", &ElecInfo_EGMCutBasedIDLoose);
    flashggStdTree->SetBranchAddress("ElecInfo.EGMCutBasedIDMedium", &ElecInfo_EGMCutBasedIDMedium);
    flashggStdTree->SetBranchAddress("ElecInfo.EGMCutBasedIDTight", &ElecInfo_EGMCutBasedIDTight);
    flashggStdTree->SetBranchAddress("ElecInfo.fggPhoVeto", &ElecInfo_fggPhoVeto);
    //flashggStdTree->SetBranchAddress("ElecInfo.tmpPhoVeto", &ElecInfo_tmpPhoVeto);
    flashggStdTree->SetBranchAddress("ElecInfo.EnergyCorrFactor", &ElecInfo_EnergyCorrFactor);
    flashggStdTree->SetBranchAddress("ElecInfo.EnergyPostCorrErr", &ElecInfo_EnergyPostCorrErr);
    flashggStdTree->SetBranchAddress("ElecInfo.EnergyPostCorrScaleUp", &ElecInfo_EnergyPostCorrScaleUp);
    flashggStdTree->SetBranchAddress("ElecInfo.EnergyPostCorrScaleDown", &ElecInfo_EnergyPostCorrScaleDown);
    flashggStdTree->SetBranchAddress("ElecInfo.EnergyPostCorrSmearUp", &ElecInfo_EnergyPostCorrSmearUp);
    flashggStdTree->SetBranchAddress("ElecInfo.EnergyPostCorrSmearDown", &ElecInfo_EnergyPostCorrSmearDown);
    flashggStdTree->SetBranchAddress("ElecInfo.GenMatch", &ElecInfo_GenMatch);
    flashggStdTree->SetBranchAddress("ElecInfo.GenPdgID", &ElecInfo_GenPdgID);
    flashggStdTree->SetBranchAddress("ElecInfo.GenPt", &ElecInfo_GenPt);
    flashggStdTree->SetBranchAddress("ElecInfo.GenEta", &ElecInfo_GenEta);
    flashggStdTree->SetBranchAddress("ElecInfo.GenPhi", &ElecInfo_GenPhi);
    //------------------------
    flashggStdTree->SetBranchAddress("MuonInfo.Size", &MuonInfo_Size);
    flashggStdTree->SetBranchAddress("MuonInfo.Charge", &MuonInfo_Charge);
    flashggStdTree->SetBranchAddress("MuonInfo.MuonType", &MuonInfo_MuonType);
    flashggStdTree->SetBranchAddress("MuonInfo.Pt", &MuonInfo_Pt);
    flashggStdTree->SetBranchAddress("MuonInfo.Eta", &MuonInfo_Eta);
    flashggStdTree->SetBranchAddress("MuonInfo.Phi", &MuonInfo_Phi);
    flashggStdTree->SetBranchAddress("MuonInfo.Energy", &MuonInfo_Energy);
    flashggStdTree->SetBranchAddress("MuonInfo.BestTrackDz", &MuonInfo_BestTrackDz);
    flashggStdTree->SetBranchAddress("MuonInfo.BestTrackDxy", &MuonInfo_BestTrackDxy);
    flashggStdTree->SetBranchAddress("MuonInfo.PFIsoDeltaBetaCorrR04", &MuonInfo_PFIsoDeltaBetaCorrR04);
    flashggStdTree->SetBranchAddress("MuonInfo.TrackerBasedIsoR03", &MuonInfo_TrackerBasedIsoR03);
    flashggStdTree->SetBranchAddress("MuonInfo.CutBasedIdMedium", &MuonInfo_CutBasedIdMedium);
    flashggStdTree->SetBranchAddress("MuonInfo.CutBasedIdTight", &MuonInfo_CutBasedIdTight);
    flashggStdTree->SetBranchAddress("MuonInfo.GenMatch", &MuonInfo_GenMatch);
    flashggStdTree->SetBranchAddress("MuonInfo.GenPdgID", &MuonInfo_GenPdgID);
    flashggStdTree->SetBranchAddress("MuonInfo.GenPt", &MuonInfo_GenPt);
    flashggStdTree->SetBranchAddress("MuonInfo.GenEta", &MuonInfo_GenEta);
    flashggStdTree->SetBranchAddress("MuonInfo.GenPhi", &MuonInfo_GenPhi);
    //------------------------
    printf("[INFO] flashggStdTreeReader::SetBranchAddresses : Finished!\n");
}


myParameters::myParameters(){
    JetInfo_jet_pt_selection = new std::vector<float>;
    JetInfo_jet_eta_selection = new std::vector<float>;
    JetInfo_jet_phi_selection = new std::vector<float>;
    JetInfo_jet_energy_selection = new std::vector<float>;
    JetInfo_jet_diphoton_deltaR_selection = new std::vector<float>;
    JetInfo_jet_leadingPhoton_deltaR_selection = new std::vector<float>;
    JetInfo_jet_subleadingPhoton_deltaR_selection = new std::vector<float>;
    JetInfo_jet_pfDeepCSVJetTags_probb_selection = new std::vector<float>;
    JetInfo_jet_pfDeepCSVJetTags_probbb_selection = new std::vector<float>;
    ElecInfo_electron_pt_selection = new std::vector<float>;
    ElecInfo_electron_eta_selection = new std::vector<float>;
    ElecInfo_electron_phi_selection = new std::vector<float>;
    ElecInfo_electron_energy_selection = new std::vector<float>;
    ElecInfo_electron_diphoton_deltaR_selection = new std::vector<float>;
    ElecInfo_electron_leadingPhoton_deltaR_selection = new std::vector<float>;
    ElecInfo_electron_subleadingPhoton_deltaR_selection = new std::vector<float>;
    MuonInfo_muon_pt_selection = new std::vector<float>;
    MuonInfo_muon_eta_selection = new std::vector<float>;
    MuonInfo_muon_phi_selection = new std::vector<float>;
    MuonInfo_muon_energy_selection = new std::vector<float>;
    MuonInfo_muon_diphoton_deltaR_selection = new std::vector<float>;
    MuonInfo_muon_leadingPhoton_deltaR_selection = new std::vector<float>;
    MuonInfo_muon_subleadingPhoton_deltaR_selection = new std::vector<float>;
}
myParameters::~myParameters(){
    delete JetInfo_jet_pt_selection;
    delete JetInfo_jet_eta_selection;
    delete JetInfo_jet_phi_selection;
    delete JetInfo_jet_energy_selection;
    delete JetInfo_jet_diphoton_deltaR_selection;
    delete JetInfo_jet_leadingPhoton_deltaR_selection;
    delete JetInfo_jet_subleadingPhoton_deltaR_selection;
    delete JetInfo_jet_pfDeepCSVJetTags_probb_selection;
    delete JetInfo_jet_pfDeepCSVJetTags_probbb_selection;
    delete ElecInfo_electron_pt_selection;
    delete ElecInfo_electron_eta_selection;
    delete ElecInfo_electron_phi_selection;
    delete ElecInfo_electron_energy_selection;
    delete ElecInfo_electron_diphoton_deltaR_selection;
    delete ElecInfo_electron_leadingPhoton_deltaR_selection;
    delete ElecInfo_electron_subleadingPhoton_deltaR_selection;
    delete MuonInfo_muon_pt_selection;
    delete MuonInfo_muon_eta_selection;
    delete MuonInfo_muon_phi_selection;
    delete MuonInfo_muon_energy_selection;
    delete MuonInfo_muon_diphoton_deltaR_selection;
    delete MuonInfo_muon_leadingPhoton_deltaR_selection;
    delete MuonInfo_muon_subleadingPhoton_deltaR_selection;
}

void myTreeClass::InitTree(){
    mytree =  new TTree("mytree", "mytree");
}
void myTreeClass::MakeNewBranchAddresses(){
    mytree -> Branch("EvtInfo_totalEntry_before_preselection", &EvtInfo_totalEntry_before_preselection, "EvtInfo_totalEntry_before_preselection/I");
    mytree -> Branch("EvtInfo_NormalizationFactor_lumi", &EvtInfo_NormalizationFactor_lumi, "EvtInfo_NormalizationFactor_lumi/F");
    mytree -> Branch("EvtInfo_NPu", &EvtInfo_NPu, "EvtInfo_NPu/I");
    //mytree -> Branch("EvtInfo_NPu", &EvtInfo_NPu, "EvtInfo_NPu/F");
    mytree -> Branch("EvtInfo_Rho", &EvtInfo_Rho, "EvtInfo_Rho/F");
    mytree -> Branch("EvtInfo_NVtx", &EvtInfo_NVtx, "EvtInfo_NVtx/I");
    mytree -> Branch("EvtInfo_genweight", &EvtInfo_genweight, "EvtInfo_genweight/F");
    //------------------------
    mytree -> Branch("DiPhoInfo_mass", &DiPhoInfo_mass, "DiPhoInfo_mass/F");
    mytree -> Branch("DiPhoInfo_pt", &DiPhoInfo_pt, "DiPhoInfo_pt/F");
    mytree -> Branch("DiPhoInfo_eta", &DiPhoInfo_eta, "DiPhoInfo_eta/F");
    mytree -> Branch("DiPhoInfo_phi", &DiPhoInfo_phi, "DiPhoInfo_phi/F");
    mytree -> Branch("DiPhoInfo_energy", &DiPhoInfo_energy, "DiPhoInfo_energy/F");
    mytree -> Branch("DiPhoInfo_leadPt", &DiPhoInfo_leadPt, "DiPhoInfo_leadPt/F");
    mytree -> Branch("DiPhoInfo_leadEta", &DiPhoInfo_leadEta, "DiPhoInfo_leadEta/F");
    mytree -> Branch("DiPhoInfo_leadPhi", &DiPhoInfo_leadPhi, "DiPhoInfo_leadPhi/F");
    mytree -> Branch("DiPhoInfo_leadE", &DiPhoInfo_leadE, "DiPhoInfo_leadE/F");
    mytree -> Branch("DiPhoInfo_leadhoe", &DiPhoInfo_leadhoe, "DiPhoInfo_leadhoe/F");
    mytree -> Branch("DiPhoInfo_leadIDMVA", &DiPhoInfo_leadIDMVA, "DiPhoInfo_leadIDMVA/F");
    mytree -> Branch("DiPhoInfo_subleadPt", &DiPhoInfo_subleadPt, "DiPhoInfo_subleadPt/F");
    mytree -> Branch("DiPhoInfo_subleadEta", &DiPhoInfo_subleadEta, "DiPhoInfo_subleadEta/F");
    mytree -> Branch("DiPhoInfo_subleadPhi", &DiPhoInfo_subleadPhi, "DiPhoInfo_subleadPhi/F");
    mytree -> Branch("DiPhoInfo_subleadE", &DiPhoInfo_subleadE, "DiPhoInfo_subleadE/F");
    mytree -> Branch("DiPhoInfo_subleadhoe", &DiPhoInfo_subleadhoe, "DiPhoInfo_subleadhoe/F");
    mytree -> Branch("DiPhoInfo_subleadIDMVA", &DiPhoInfo_subleadIDMVA, "DiPhoInfo_subleadIDMVA/F");
    //------------------------
    mytree -> Branch("ElecInfo_Size", &ElecInfo_Size, "ElecInfo_Size/I");
    mytree -> Branch("MuonInfo_Size", &MuonInfo_Size, "MuonInfo_Size/I");
    mytree -> Branch("num_leptons", &num_leptons, "num_leptons/I");// # of selected objects.
    mytree -> Branch("num_electrons", &num_electrons, "num_electrons/I");// # of selected objects.
    mytree -> Branch("num_muons", &num_muons, "num_muons/I");// # of selected objects.
    mytree -> Branch("ElecInfo_electron_pt", &ElecInfo_electron_pt);
    mytree -> Branch("ElecInfo_electron_eta", &ElecInfo_electron_eta);
    mytree -> Branch("ElecInfo_electron_phi", &ElecInfo_electron_phi);
    mytree -> Branch("ElecInfo_electron_energy", &ElecInfo_electron_energy);
    mytree -> Branch("ElecInfo_electron_diphoton_deltaR", &ElecInfo_electron_diphoton_deltaR);
    mytree -> Branch("ElecInfo_electron_leadingPhoton_deltaR", &ElecInfo_electron_leadingPhoton_deltaR);
    mytree -> Branch("ElecInfo_electron_subleadingPhoton_deltaR", &ElecInfo_electron_subleadingPhoton_deltaR);
    mytree -> Branch("MuonInfo_muon_pt", &MuonInfo_muon_pt);
    mytree -> Branch("MuonInfo_muon_eta", &MuonInfo_muon_eta);
    mytree -> Branch("MuonInfo_muon_phi", &MuonInfo_muon_phi);
    mytree -> Branch("MuonInfo_muon_energy", &MuonInfo_muon_energy);
    mytree -> Branch("MuonInfo_muon_diphoton_deltaR", &MuonInfo_muon_diphoton_deltaR);
    mytree -> Branch("MuonInfo_muon_leadingPhoton_deltaR", &MuonInfo_muon_leadingPhoton_deltaR);
    mytree -> Branch("MuonInfo_muon_subleadingPhoton_deltaR", &MuonInfo_muon_subleadingPhoton_deltaR);
    //------------------------
    mytree -> Branch("jets_size", &jets_size, "jets_size/I");
    mytree -> Branch("num_jets", &num_jets, "num_jets/I");
    mytree -> Branch("JetInfo_jet_pt", &JetInfo_jet_pt);
    mytree -> Branch("JetInfo_jet_eta", &JetInfo_jet_eta);
    mytree -> Branch("JetInfo_jet_phi", &JetInfo_jet_phi);
    mytree -> Branch("JetInfo_jet_energy", &JetInfo_jet_energy);
    mytree -> Branch("JetInfo_jet_diphoton_deltaR", &JetInfo_jet_diphoton_deltaR);
    mytree -> Branch("JetInfo_jet_leadingPhoton_deltaR", &JetInfo_jet_leadingPhoton_deltaR);
    mytree -> Branch("JetInfo_jet_subleadingPhoton_deltaR", &JetInfo_jet_subleadingPhoton_deltaR);
    mytree -> Branch("JetInfo_jet_pfDeepCSVJetTags_probb", &JetInfo_jet_pfDeepCSVJetTags_probb);
    mytree -> Branch("JetInfo_jet_pfDeepCSVJetTags_probbb", &JetInfo_jet_pfDeepCSVJetTags_probbb);
    //mytree -> Branch("num_btagged_jets", &num_btagged_jets, "num_btagged_jets/I");
    //mytree -> Branch("num_nonbtagged_jets", &num_nonbtagged_jets, "num_nonbtagged_jets/I");
    //------------------------
    //mytree -> Branch("inv_mass_dijet", &inv_mass_dijet, "inv_mass_dijet/F");
    //mytree -> Branch("inv_mass_diphoton", &inv_mass_diphoton, "inv_mass_diphoton/F");
    //mytree -> Branch("inv_mass_tbw", &inv_mass_tbw, "inv_mass_tbw/F");
}
void myTreeClass::Fill(){
    mytree -> Fill();
}
void myParameters::Clear(){
    EvtInfo_totalEntry_before_preselection = 0;
    EvtInfo_NormalizationFactor_lumi = 0;
    //------------------------
    DiPhoInfo_eta = 0;
    DiPhoInfo_phi = 0;
    DiPhoInfo_energy = 0;
    //------------------------
    num_jets = 0;
    JetInfo_jet_pt.clear();
    JetInfo_jet_eta.clear();
    JetInfo_jet_phi.clear();
    JetInfo_jet_energy.clear();
    JetInfo_jet_diphoton_deltaR.clear();
    JetInfo_jet_leadingPhoton_deltaR.clear();
    JetInfo_jet_subleadingPhoton_deltaR.clear();
    JetInfo_jet_pfDeepCSVJetTags_probb.clear();
    JetInfo_jet_pfDeepCSVJetTags_probbb.clear();
    //------------------------
    num_leptons = 0;
    num_electrons = 0;
    num_muons = 0;
    ElecInfo_electron_pt.clear();
    ElecInfo_electron_eta.clear();
    ElecInfo_electron_phi.clear();
    ElecInfo_electron_energy.clear();
    ElecInfo_electron_diphoton_deltaR.clear();
    ElecInfo_electron_leadingPhoton_deltaR.clear();
    ElecInfo_electron_subleadingPhoton_deltaR.clear();
    MuonInfo_muon_pt.clear();
    MuonInfo_muon_eta.clear();
    MuonInfo_muon_phi.clear();
    MuonInfo_muon_energy.clear();
    MuonInfo_muon_diphoton_deltaR.clear();
    MuonInfo_muon_leadingPhoton_deltaR.clear();
    MuonInfo_muon_subleadingPhoton_deltaR.clear();
    //------------------------
    //Not used in preselection stage
    //------------------------
    num_btagged_jets = 0;
    num_nonbtagged_jets = 0;
    //------------------------
    //Chi-2 sorting related
    //------------------------
    inv_mass_dijet = 0;
    inv_mass_diphoton = 0;
    inv_mass_tbw = 0;
    //------------------------
    JetInfo_dijet_delta_eta = 0;
    JetInfo_dijet_delta_phi = 0;
    JetInfo_dijet_delta_angle = 0;
    //------------------------
}
/*
        std::vector<int> vec_bjet_indices;
        std::vector<TLorentzVector> vec_btagged_jets;
        std::vector<TLorentzVector> vec_nonbtagged_jets;
        //=================================//
        //-----  Select B-tagged Jet  -----//
        //=================================//
        int bjetindex=-1;
        TLorentzVector leading_bjet;
        for(int i=0; i<treeReader.jets_size; i++){
            if( fabs(treeReader.JetInfo_Eta->at(i)) > 2.4 ) continue;
            if( fabs(treeReader.JetInfo_Pt->at(i))  < 25  ) continue;
            if(treeReader.JetInfo_pfDeepCSVJetTags_probb->at(i)+treeReader.JetInfo_pfDeepCSVJetTags_probbb->at(i) >= pfDeepCSVJetTags_tight){
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
        //==============================//
        //-----  Select Rest Jets  -----//
        //==============================//
        bool isbjet;
        for(int i=0; i<treeReader.jets_size; i++){
            if( fabs(treeReader.JetInfo_Eta->at(i)) > 2.4 ) continue;
            if( fabs(treeReader.JetInfo_Pt->at(i))  < 25  ) continue;
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
        mytree.JetInfo_jet1_pt = CUT_num_nonbjets_less_than_2 ? -999. : vec_nonbtagged_jets[0].Pt();
        mytree.JetInfo_jet1_eta = CUT_num_nonbjets_less_than_2 ? -999. : vec_nonbtagged_jets[0].Eta();
        mytree.JetInfo_jet1_phi = CUT_num_nonbjets_less_than_2 ? -999. : vec_nonbtagged_jets[0].Phi();
        mytree.JetInfo_jet2_pt = CUT_num_nonbjets_less_than_2 ? -999. : vec_nonbtagged_jets[1].Pt();
        mytree.JetInfo_jet2_eta = CUT_num_nonbjets_less_than_2 ? -999. : vec_nonbtagged_jets[1].Eta();
        mytree.JetInfo_jet2_phi = CUT_num_nonbjets_less_than_2 ? -999. : vec_nonbtagged_jets[1].Phi();
*/
