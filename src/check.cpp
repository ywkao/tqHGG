#include <stdio.h>
#include <vector>
#include <TCanvas.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <iostream>
#include "../include/cross_section.h"
using namespace std;
double pfDeepCSVJetTags_loose  = 0.1522;
double pfDeepCSVJetTags_medium = 0.4941;
double pfDeepCSVJetTags_tight  = 0.8001;

bool isThisDataOrNot(char* dataset);
bool isThisMultiFile(char* dataset);


int main(int argc, char *argv[]){
    //HelloWorld{{{
    printf("Hello World.\n");
    char input_file[512]; sprintf(input_file, "%s", argv[1]); printf("[INFO] input_file  = %s\n", input_file);
    char dataset[512]; sprintf(dataset, "%s", argv[2]);       printf("[INFO] dataset     = %s\n", dataset);
    bool isData = isThisDataOrNot(dataset);
    bool isMultiFile = isThisMultiFile(dataset);
    TChain *flashggStdTree = new TChain("flashggNtuples/flashggStdTree");
    if(!isMultiFile) flashggStdTree->Add(input_file);
    else             flashggStdTree->Add(Form("%s/*.root", input_file));

    TFile *f_mcpu = TFile::Open("data/MCPileUp.root");
    TH1D  *h_pu_reweighting_factor = (TH1D*)f_mcpu->Get("puhist");
    //}}}
    //SetBranchAddresses{{{
    bool  EvtInfo_passTrigger;
    Int_t EvtInfo_NVtx;
    Int_t EvtInfo_NPu;
    float EvtInfo_genweight;
    float DiPhoInfo_mass;
    float DiPhoInfo_pt;
    float DiPhoInfo_leadPt;
    float DiPhoInfo_leadEta;
    float DiPhoInfo_leadPhi;
    float DiPhoInfo_leadE;
    float DiPhoInfo_leadhoe;
    float DiPhoInfo_leadIDMVA;
    float DiPhoInfo_subleadPt;
    float DiPhoInfo_subleadEta;
    float DiPhoInfo_subleadPhi;
    float DiPhoInfo_subleadE;
    float DiPhoInfo_subleadhoe;
    float DiPhoInfo_subleadIDMVA;
    Int_t jets_size;
    std::vector<float> *JetInfo_Pt;
    std::vector<float> *JetInfo_Eta;
    std::vector<float> *JetInfo_Phi;
    std::vector<float> *JetInfo_Mass;
    std::vector<float> *JetInfo_Energy;
    std::vector<float> *JetInfo_pfDeepCSVJetTags_probb;
    std::vector<float> *JetInfo_pfDeepCSVJetTags_probbb;
    std::vector<float>  JetInfo_jet_pfDeepCSVJetTags_probb;//bjet purpose
    std::vector<float>  JetInfo_jet_pfDeepCSVJetTags_probbb;//bjet purpose
    Int_t           ElecInfo_Size;
    std::vector<int>     *ElecInfo_Charge;
    std::vector<float>   *ElecInfo_Pt;
    std::vector<float>   *ElecInfo_Eta;
    std::vector<float>   *ElecInfo_Phi;
    std::vector<float>   *ElecInfo_Energy;
    std::vector<bool>    *ElecInfo_EGMCutBasedIDVeto;
    std::vector<bool>    *ElecInfo_EGMCutBasedIDLoose;
    std::vector<bool>    *ElecInfo_EGMCutBasedIDMedium;
    std::vector<bool>    *ElecInfo_EGMCutBasedIDTight;
    Int_t           MuonInfo_Size;
    std::vector<int>     *MuonInfo_Charge;
    std::vector<float>   *MuonInfo_MuonType;
    std::vector<float>   *MuonInfo_Pt;
    std::vector<float>   *MuonInfo_Eta;
    std::vector<float>   *MuonInfo_Phi;
    std::vector<float>   *MuonInfo_Energy;
    std::vector<float>   *MuonInfo_PFIsoDeltaBetaCorrR04;
    std::vector<float>   *MuonInfo_TrackerBasedIsoR03;
    std::vector<bool>    *MuonInfo_CutBasedIdMedium;
    std::vector<bool>    *MuonInfo_CutBasedIdTight;
    //------------------------
    JetInfo_Pt = 0;
    JetInfo_Eta = 0;
    JetInfo_Phi = 0;
    JetInfo_Mass = 0;
    JetInfo_Energy = 0;
    JetInfo_pfDeepCSVJetTags_probb = 0;
    JetInfo_pfDeepCSVJetTags_probbb = 0;
    ElecInfo_Charge = 0;
    ElecInfo_Pt = 0;
    ElecInfo_Eta = 0;
    ElecInfo_Phi = 0;
    ElecInfo_Energy = 0;
    ElecInfo_EGMCutBasedIDVeto = 0;
    ElecInfo_EGMCutBasedIDLoose = 0;
    ElecInfo_EGMCutBasedIDMedium = 0;
    ElecInfo_EGMCutBasedIDTight = 0;
    MuonInfo_Charge = 0;
    MuonInfo_MuonType = 0;
    MuonInfo_Pt = 0;
    MuonInfo_Eta = 0;
    MuonInfo_Phi = 0;
    MuonInfo_Energy = 0;
    MuonInfo_PFIsoDeltaBetaCorrR04 = 0;
    MuonInfo_TrackerBasedIsoR03 = 0;
    MuonInfo_CutBasedIdMedium = 0;
    MuonInfo_CutBasedIdTight = 0;
    //------------------------
    flashggStdTree->SetBranchAddress("EvtInfo.passTrigger", &EvtInfo_passTrigger);
    flashggStdTree->SetBranchAddress("EvtInfo.NVtx", &EvtInfo_NVtx);
    flashggStdTree->SetBranchAddress("EvtInfo.NPu", &EvtInfo_NPu);
    flashggStdTree->SetBranchAddress("EvtInfo.genweight", &EvtInfo_genweight);
    flashggStdTree->SetBranchAddress("DiPhoInfo.mass", &DiPhoInfo_mass);
    flashggStdTree->SetBranchAddress("DiPhoInfo.pt", &DiPhoInfo_pt);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadPt", &DiPhoInfo_leadPt);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadEta", &DiPhoInfo_leadEta);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadPhi", &DiPhoInfo_leadPhi);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadE", &DiPhoInfo_leadE);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadhoe", &DiPhoInfo_leadhoe);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadIDMVA", &DiPhoInfo_leadIDMVA);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadPt", &DiPhoInfo_subleadPt);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadEta", &DiPhoInfo_subleadEta);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadPhi", &DiPhoInfo_subleadPhi);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadE", &DiPhoInfo_subleadE);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadhoe", &DiPhoInfo_subleadhoe);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadIDMVA", &DiPhoInfo_subleadIDMVA);
    flashggStdTree->SetBranchAddress("jets_size", &jets_size);
    flashggStdTree->SetBranchAddress("JetInfo.Pt", &JetInfo_Pt);
    flashggStdTree->SetBranchAddress("JetInfo.Eta", &JetInfo_Eta);
    flashggStdTree->SetBranchAddress("JetInfo.Phi", &JetInfo_Phi);
    flashggStdTree->SetBranchAddress("JetInfo.Mass", &JetInfo_Mass);
    flashggStdTree->SetBranchAddress("JetInfo.Energy", &JetInfo_Energy);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probb", &JetInfo_pfDeepCSVJetTags_probb);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probbb", &JetInfo_pfDeepCSVJetTags_probbb);
    flashggStdTree->SetBranchAddress("ElecInfo.Size", &ElecInfo_Size);
    flashggStdTree->SetBranchAddress("ElecInfo.Charge", &ElecInfo_Charge);
    flashggStdTree->SetBranchAddress("ElecInfo.Pt", &ElecInfo_Pt);
    flashggStdTree->SetBranchAddress("ElecInfo.Eta", &ElecInfo_Eta);
    flashggStdTree->SetBranchAddress("ElecInfo.Phi", &ElecInfo_Phi);
    flashggStdTree->SetBranchAddress("ElecInfo.Energy", &ElecInfo_Energy);
    flashggStdTree->SetBranchAddress("ElecInfo.EGMCutBasedIDVeto", &ElecInfo_EGMCutBasedIDVeto);
    flashggStdTree->SetBranchAddress("ElecInfo.EGMCutBasedIDLoose", &ElecInfo_EGMCutBasedIDLoose);
    flashggStdTree->SetBranchAddress("ElecInfo.EGMCutBasedIDMedium", &ElecInfo_EGMCutBasedIDMedium);
    flashggStdTree->SetBranchAddress("ElecInfo.EGMCutBasedIDTight", &ElecInfo_EGMCutBasedIDTight);
    flashggStdTree->SetBranchAddress("MuonInfo.Size", &MuonInfo_Size);
    flashggStdTree->SetBranchAddress("MuonInfo.Charge", &MuonInfo_Charge);
    flashggStdTree->SetBranchAddress("MuonInfo.MuonType", &MuonInfo_MuonType);
    flashggStdTree->SetBranchAddress("MuonInfo.Pt", &MuonInfo_Pt);
    flashggStdTree->SetBranchAddress("MuonInfo.Eta", &MuonInfo_Eta);
    flashggStdTree->SetBranchAddress("MuonInfo.Phi", &MuonInfo_Phi);
    flashggStdTree->SetBranchAddress("MuonInfo.Energy", &MuonInfo_Energy);
    flashggStdTree->SetBranchAddress("MuonInfo.PFIsoDeltaBetaCorrR04", &MuonInfo_PFIsoDeltaBetaCorrR04);
    flashggStdTree->SetBranchAddress("MuonInfo.TrackerBasedIsoR03", &MuonInfo_TrackerBasedIsoR03);
    flashggStdTree->SetBranchAddress("MuonInfo.CutBasedIdMedium", &MuonInfo_CutBasedIdMedium);
    flashggStdTree->SetBranchAddress("MuonInfo.CutBasedIdTight", &MuonInfo_CutBasedIdTight);
    //}}}End of SetBranchAddresses
    //Set up parameters{{{
    Int_t num_jets;
    Int_t num_leptons;

    int nentries = flashggStdTree->GetEntries(); printf("[INFO] nentries = %d\n", nentries);

    double TotalGenweight=0;
    double total=0;
    for(int ientry=0; ientry<nentries; ientry++){
        flashggStdTree->GetEntry(ientry);//load data
        TotalGenweight += EvtInfo_genweight;
    }
    
    double Luminosity = 41.53; //fb{-1}
    double CrossSection = GetXsec(dataset); //pb
    double BranchingFraction = 1.; //pb
    double NormalizationFactor = 1000. * Luminosity * CrossSection * BranchingFraction / TotalGenweight;
    printf("[INFO] Luminosity = %f\n", Luminosity);
    printf("[INFO] CrossSection = %f\n", CrossSection);
    printf("[INFO] TotalGenweight = %.7f\n", TotalGenweight);
    printf("[INFO] NormalizationFactor = %.7f\n", NormalizationFactor);
    //}}}
    //Counters {{{
    double counter_nocut = 0;
    double counter_cut0  = 0;
    double counter_cut1  = 0;
    double counter_cut2  = 0;
    double counter_cut3  = 0;
    double counter_cut4  = 0;
    double counter_cut5  = 0;

    double counter_cut_hadronic  = 0;

    double counter_nocut_norm = 0;
    double counter_cut0_norm  = 0;
    double counter_cut1_norm  = 0;
    double counter_cut2_norm  = 0;
    double counter_cut3_norm  = 0;
    double counter_cut4_norm  = 0;
    double counter_cut5_norm  = 0;

    double counter_nocut_successive = 0;
    double counter_cut0_successive  = 0;
    double counter_cut1_successive  = 0;
    double counter_cut2_successive  = 0;
    double counter_cut3_successive  = 0;
    double counter_cut4_successive  = 0;
    double counter_cut5_successive  = 0;

    double counter_nocut_norm_successive = 0;
    double counter_cut0_norm_successive  = 0;
    double counter_cut1_norm_successive  = 0;
    double counter_cut2_norm_successive  = 0;
    double counter_cut3_norm_successive  = 0;
    double counter_cut4_norm_successive  = 0;
    double counter_cut5_norm_successive  = 0;
    //}}}

    TH1D *hist = new TH1D("hist", "; Number of b-tagged jets; Entries", 10, 0, 10);
    TH1D *hist_probb = new TH1D("hist_probb", "; JetInfo_pfDeepCSVJetTags_probb; Entries", 10, 0, 1.);
    TH1D *hist_probbb = new TH1D("hist_probbb", "; JetInfo_pfDeepCSVJetTags_probbb; Entries", 10, 0, 1.);
    TH1D *hist_jet_probb  = new TH1D("hist_jet_probb",  "; JetInfo_jet_pfDeepCSVJetTags_probb; Entries", 10, 0, 1.);
    TH1D *hist_jet_probbb = new TH1D("hist_jet_probbb", "; JetInfo_jet_pfDeepCSVJetTags_probbb; Entries", 10, 0, 1.);


    for(int ientry=0; ientry<nentries; ientry++){
        //printf("ientry = %d\n", ientry);
        flashggStdTree->GetEntry(ientry);//load data

        //========= PU Reweighting =========//
        bool apply_PU_reweighting = true;
        double PU_reweighting_factor = apply_PU_reweighting ? h_pu_reweighting_factor->GetBinContent((int)EvtInfo_NPu+1) : 1.;
        double effectiveYields = EvtInfo_genweight * NormalizationFactor * PU_reweighting_factor; //### the final weighted value will be different....
        
        counter_nocut += 1;
        counter_nocut_norm += effectiveYields;
        counter_nocut_successive += 1;
        counter_nocut_norm_successive += effectiveYields;
        JetInfo_jet_pfDeepCSVJetTags_probb.clear();//bjet purpose
        JetInfo_jet_pfDeepCSVJetTags_probbb.clear();//bjet purpose
        //if(ientry<1000) { total += EvtInfo_genweight*NormalizationFactor; printf("ie = %d, yields = %11.7f\n", ientry+1, EvtInfo_genweight * NormalizationFactor); }
        //Store photon info{{{
        //==============================//
        //-----  Store PhotonInfo  -----//
        //==============================//
        //bool pass_leadingPhotonPT = DiPhoInfo_leadPt > DiPhoInfo_mass / 2.;
        //bool pass_subleadingPhotonPT = DiPhoInfo_subleadPt > DiPhoInfo_mass / 4.;
        bool pass_leadingPhotonPT = DiPhoInfo_leadPt > 35.;
        bool pass_subleadingPhotonPT = DiPhoInfo_subleadPt > 25.;
        bool pass_photon_criteria_pt = pass_leadingPhotonPT && pass_subleadingPhotonPT;
        bool pass_leadingPhotonEta =  (DiPhoInfo_leadEta < 1.4442) || (DiPhoInfo_leadEta > 1.566 && DiPhoInfo_leadEta < 2.5);
        bool pass_subleadingPhotonEta = (DiPhoInfo_subleadEta < 1.4442) || (DiPhoInfo_subleadEta > 1.566 && DiPhoInfo_subleadEta < 2.5);
        bool pass_photon_criteria_eta = pass_leadingPhotonEta && pass_subleadingPhotonEta;
        bool pass_photon_criteria = pass_photon_criteria_pt && pass_photon_criteria_eta;

        TLorentzVector leading_photon, subleading_photon, diphoton;
        leading_photon.SetPtEtaPhiE(DiPhoInfo_leadPt, DiPhoInfo_leadEta, DiPhoInfo_leadPhi, DiPhoInfo_leadE);
        subleading_photon.SetPtEtaPhiE(DiPhoInfo_subleadPt, DiPhoInfo_subleadEta, DiPhoInfo_subleadPhi, DiPhoInfo_subleadE);
        diphoton = leading_photon + subleading_photon;
        //}}}
        //Select leptons{{{
        //==============================//
        //-----  Select Electrons  -----//
        //==============================//
        int num_electrons = 0;
        std::vector<TLorentzVector> Leptons;
        std::vector<TLorentzVector> Electrons;
        bool bool_AtLeastOneElectron = ElecInfo_Size>0 ? true : false;//ElecInfo_Size = -999 => event without diphoton candidate
        if(bool_AtLeastOneElectron){
            for(int i=0; i<ElecInfo_Size; i++){
                if( !ElecInfo_EGMCutBasedIDMedium->at(i) ) continue;
                if( fabs(ElecInfo_Eta->at(i)) > 2.4 ) continue;
                if( fabs(ElecInfo_Eta->at(i)) > 1.4442 && fabs(ElecInfo_Eta->at(i)) < 1.566 ) continue;
                if( fabs(ElecInfo_Pt->at(i))  < 20  ) continue;
                //--- check deltaR(electron,photon) ---//
                TLorentzVector electron; 
                electron.SetPtEtaPhiE(ElecInfo_Pt->at(i), ElecInfo_Eta->at(i), ElecInfo_Phi->at(i), ElecInfo_Energy->at(i));
                double delta_R = electron.DeltaR(leading_photon);
                if( delta_R<0.3 ) continue;
                delta_R = electron.DeltaR(subleading_photon);
                if( delta_R<0.3 ) continue;
                //--- store information ---//
                num_electrons+=1;
                Electrons.push_back(electron);
                Leptons.push_back(electron);

                /*
                if(ientry==10642){
                    printf("(%d) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f, Energy = %6.2f, Mass = %6.2f (electron)\n", i, electron.Pt(), electron.Eta(), electron.Phi(), electron.Energy(), electron.M());
                    printf("(%d) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f, Energy = %6.2f, Mass = %6.2f (leadingPhoton)\n", i, leading_photon.Pt(), leading_photon.Eta(), leading_photon.Phi(), leading_photon.Energy(), leading_photon.M());
                    printf("(%d) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f, Energy = %6.2f, Mass = %6.2f (subleadingPhoton)\n", i, subleading_photon.Pt(), subleading_photon.Eta(), subleading_photon.Phi(), subleading_photon.Energy(), subleading_photon.M());
                }
                */
            }
        }
        bool bool_AtLeastOneSelectedElectron = num_electrons>0 ? true : false;//for calculation of deltaR(e,j).
        //==========================//
        //-----  Select Muons  -----//
        //==========================//
        int num_muons = 0;
        std::vector<TLorentzVector> Muons;
        bool bool_AtLeastOneMuon = MuonInfo_Size>0 ? true : false;//MuonInfo_Size = -999 => event without diphoton candidate
        if(bool_AtLeastOneMuon){
            for(int i=0; i<MuonInfo_Size; i++){
                if( !MuonInfo_CutBasedIdTight->at(i) ) continue;
                if( fabs(MuonInfo_Eta->at(i)) > 2.4 ) continue;
                if( fabs(MuonInfo_Pt->at(i))  < 20  ) continue;
                if( MuonInfo_PFIsoDeltaBetaCorrR04->at(i)  > 0.25  ) continue;
                //--- check deltaR(muon,photon) ---//
                TLorentzVector muon; 
                muon.SetPtEtaPhiE(MuonInfo_Pt->at(i), MuonInfo_Eta->at(i), MuonInfo_Phi->at(i), MuonInfo_Energy->at(i));
                double delta_R = muon.DeltaR(leading_photon);
                if( delta_R<0.3 ) continue;
                delta_R = muon.DeltaR(subleading_photon);
                if( delta_R<0.3 ) continue;
                //--- store information ---//
                num_muons+=1;
                Muons.push_back(muon);
                Leptons.push_back(muon);
                /*
                if(ientry==10642){
                    printf("(%d) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f, Energy = %6.2f, Mass = %6.2f (muon)\n", i, muon.Pt(), muon.Eta(), muon.Phi(), muon.Energy(), muon.M());
                }
                */
            }
        }
        bool bool_AtLeastOneSelectedMuon = num_muons>0 ? true : false;//for calculation of deltaR(mu,j).
        num_leptons = num_electrons + num_muons;
        //}}}
        //Select jets{{{
        //=========================//
        //-----  Select Jets  -----//
        //=========================//
        num_jets = 0;
        std::vector<TLorentzVector> Jets;
        bool bool_AtLeastOneJet = jets_size>0 ? true : false;//jets_size = -999 => event without diphoton candidate
        //--------------------------------------------------
        bool bool_cut02;
        if(bool_AtLeastOneJet){
            for(int i=0; i<jets_size; i++){
                //flashgg::Tight2017 jet ID #Already applied flashgg package.
                if( fabs(JetInfo_Eta->at(i)) > 2.4 ) continue;
                if( fabs(JetInfo_Pt->at(i))  < 25  ) continue;
                num_jets+=1;
            }
        }
        if(num_jets>0) bool_cut02 = true;
        else bool_cut02 = false;
        //--------------------------------------------------
        num_jets = 0;
        if(bool_AtLeastOneJet){
            for(int i=0; i<jets_size; i++){
                //flashgg::Tight2017 jet ID #Already applied flashgg package.
                if( fabs(JetInfo_Eta->at(i)) > 2.4 ) continue;
                if( fabs(JetInfo_Pt->at(i))  < 25  ) continue;
                //--- check deltaR(jet,photon) ---//
                TLorentzVector jet; 
                jet.SetPtEtaPhiE(JetInfo_Pt->at(i), JetInfo_Eta->at(i), JetInfo_Phi->at(i), JetInfo_Energy->at(i));
                double delta_R = jet.DeltaR(leading_photon);
                if( delta_R<0.4 ) continue;
                delta_R = jet.DeltaR(subleading_photon);
                if( delta_R<0.4 ) continue;
                //--- check deltaR(jet,lep) ---//
                bool bool_passJetLeptonSeparation = true;//if no leptons selected, the jet pass the delta_R criterion automatically.
                if(num_leptons>0){
                    for(int j=0; j<num_leptons; j++){
                        TLorentzVector lepton = Leptons.at(j);
                        delta_R = jet.DeltaR(lepton);
                        if( delta_R<0.4 ) bool_passJetLeptonSeparation = false;
                    }
                }
                if(!bool_passJetLeptonSeparation) continue;//if not pass, reject the jet.
                /*
                //--- check deltaR(jet,e) ---//
                if(bool_AtLeastOneSelectedElectron){
                    for(int i=0; i<num_electrons; i++){
                        delta_R = jet.DeltaR(Electrons.at(i));
                        if( delta_R<0.4 ) bool_passJetLeptonSeparation = false;
                    }
                }
                if(!bool_passJetLeptonSeparation) continue;//if not pass, reject the jet.
                //--- check deltaR(jet,mu) ---//
                if(bool_AtLeastOneSelectedMuon){
                    for(int i=0; i<num_muons; i++){
                        delta_R = jet.DeltaR(Muons.at(i));
                        if( delta_R<0.4 ) bool_passJetLeptonSeparation = false;
                    }
                }
                if(!bool_passJetLeptonSeparation) continue;//if not pass, reject the jet.
                */
                //--- store information ---//
                num_jets+=1;
                Jets.push_back(jet);
                JetInfo_jet_pfDeepCSVJetTags_probb.push_back(JetInfo_pfDeepCSVJetTags_probb->at(i));//bjet purpose
                JetInfo_jet_pfDeepCSVJetTags_probbb.push_back(JetInfo_pfDeepCSVJetTags_probbb->at(i));//bjet purpose
                //printf("[check] JetInfo_pfDeepCSVJetTags_probb = %f\n", JetInfo_pfDeepCSVJetTags_probb->at(i));
                //printf("[check] JetInfo_jet_pfDeepCSVJetTags_probb = %f\n", JetInfo_jet_pfDeepCSVJetTags_probb.at(num_jets-1));
                //printf("[check] JetInfo_pfDeepCSVJetTags_probbb = %f\n", JetInfo_pfDeepCSVJetTags_probbb->at(i));
                //printf("[check] JetInfo_jet_pfDeepCSVJetTags_probbb = %f\n\n", JetInfo_jet_pfDeepCSVJetTags_probbb.at(num_jets-1));
                hist_probb -> Fill(JetInfo_pfDeepCSVJetTags_probb->at(i));
                hist_probbb -> Fill(JetInfo_pfDeepCSVJetTags_probbb->at(i));
                hist_jet_probb  -> Fill(JetInfo_jet_pfDeepCSVJetTags_probb.at(num_jets-1));
                hist_jet_probbb -> Fill(JetInfo_jet_pfDeepCSVJetTags_probbb.at(num_jets-1));
                //printf("%2d JetInfo_pfDeepCSVJetTags_probb = %f\n", i, JetInfo_pfDeepCSVJetTags_probb->at(i));
                //printf("%2d JetInfo_pfDeepCSVJetTags_probbb = %f\n", i, JetInfo_pfDeepCSVJetTags_probbb->at(i));
            }
        }
        //}}}
        //Select bjet{{{
        //----- b-jet -----//
        int num_bjets = 0;
        for(int i=0; i<num_jets; ++i){
            if(JetInfo_jet_pfDeepCSVJetTags_probb.at(i)+JetInfo_jet_pfDeepCSVJetTags_probbb.at(i) >= pfDeepCSVJetTags_medium){
                //printf("JetInfo_jet_pfDeepCSVJetTags_probb = %f\n", JetInfo_jet_pfDeepCSVJetTags_probb.at(i));
                //printf("JetInfo_jet_pfDeepCSVJetTags_probbb = %f\n", JetInfo_jet_pfDeepCSVJetTags_probbb.at(i));
                num_bjets += 1;
            }
        }//end of looping jets
        if(!(num_jets==0) && pass_photon_criteria && EvtInfo_passTrigger){
            //printf("[check] %7d num_bjets = %d, num_jets = %d\n", ientry, num_bjets, num_jets);
            hist->Fill(num_bjets);
        }
        //}}}
        //Event selection (counters){{{
        //==========================//
        //----- EventSelection -----//
        //==========================//
        //Individual cuts{{{
        if(pass_photon_criteria && EvtInfo_passTrigger) counter_cut0 += 1;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_leptons>0 ) counter_cut1 += 1;
        if(pass_photon_criteria && EvtInfo_passTrigger && bool_cut02 ) counter_cut2 += 1;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_jets>0 ) counter_cut3 += 1;
        if(pass_photon_criteria && EvtInfo_passTrigger && ((DiPhoInfo_mass > 100 && DiPhoInfo_mass < 120) || (DiPhoInfo_mass > 130 && DiPhoInfo_mass < 180)) ) counter_cut4 += 1;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_bjets>0 ) counter_cut5 += 1;
        //---
        if(pass_photon_criteria && EvtInfo_passTrigger) counter_cut0_norm += effectiveYields;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_leptons>0 ) counter_cut1_norm += effectiveYields;
        if(pass_photon_criteria && EvtInfo_passTrigger && bool_cut02 ) counter_cut2_norm += effectiveYields;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_jets>0 ) counter_cut3_norm += effectiveYields;
        if(pass_photon_criteria && EvtInfo_passTrigger && ((DiPhoInfo_mass > 100 && DiPhoInfo_mass < 120) || (DiPhoInfo_mass > 130 && DiPhoInfo_mass < 180)) ) counter_cut4_norm += effectiveYields;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_bjets>0 ) counter_cut5_norm += effectiveYields;
        //}}}
        //Successive cuts{{{
        if(pass_photon_criteria && EvtInfo_passTrigger) counter_cut0_successive += 1;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_leptons>0 ) counter_cut1_successive += 1;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_leptons>0 && bool_cut02 ) counter_cut2_successive += 1;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_leptons>0 && bool_cut02 && num_jets>0 ) counter_cut3_successive += 1;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_leptons>0 && bool_cut02 && num_jets>0 && ((DiPhoInfo_mass > 100 && DiPhoInfo_mass < 120) || (DiPhoInfo_mass > 130 && DiPhoInfo_mass < 180)) ) counter_cut4_successive += 1;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_leptons>0 && bool_cut02 && num_jets>0 && ((DiPhoInfo_mass > 100 && DiPhoInfo_mass < 120) || (DiPhoInfo_mass > 130 && DiPhoInfo_mass < 180)) && num_bjets>0) counter_cut5_successive += 1;
        //---
        if(pass_photon_criteria && EvtInfo_passTrigger) counter_cut0_norm_successive += effectiveYields;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_leptons>0 ) counter_cut1_norm_successive += effectiveYields;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_leptons>0 && bool_cut02 ) counter_cut2_norm_successive += effectiveYields;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_leptons>0 && bool_cut02 && num_jets>0 ) counter_cut3_norm_successive += effectiveYields;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_leptons>0 && bool_cut02 && num_jets>0 && ((DiPhoInfo_mass > 100 && DiPhoInfo_mass < 120) || (DiPhoInfo_mass > 130 && DiPhoInfo_mass < 180)) ) counter_cut4_norm_successive += effectiveYields;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_leptons>0 && bool_cut02 && num_jets>0 && ((DiPhoInfo_mass > 100 && DiPhoInfo_mass < 120) || (DiPhoInfo_mass > 130 && DiPhoInfo_mass < 180)) && num_bjets>0) counter_cut5_norm_successive += effectiveYields;
        //}}}
        // MoreInfo (hadronic){{{
        if(pass_photon_criteria && EvtInfo_passTrigger && num_leptons<1 && bool_cut02 && num_jets>3 && ((DiPhoInfo_mass > 100 && DiPhoInfo_mass < 120) || (DiPhoInfo_mass > 130 && DiPhoInfo_mass < 180)) && num_bjets>0) counter_cut_hadronic += 1;
        
        // }}}
        //}}}
    }//end of event loop.

    printf("\nIndividual selection\n");
    printf("--------------------\n");
    printf("[INFO] Hadronic: Entries = %.0f\n", counter_cut_hadronic);

    //report individual{{{
    printf("\nIndividual selection\n");
    printf("--------------------\n");
    printf("[INFO] Individual NoCut: Entries = %.0f\n", counter_nocut);
    printf("[INFO] Individual Cut00: Entries = %.0f\n", counter_cut0);
    printf("[INFO] Individual Cut01: Entries = %.0f\n", counter_cut1);
    printf("[INFO] Individual Cut02: Entries = %.0f\n", counter_cut2);
    printf("[INFO] Individual Cut03: Entries = %.0f\n", counter_cut3);
    printf("[INFO] Individual Cut04: Entries = %.0f\n", counter_cut4);
    printf("[INFO] Individual Cut05: Entries = %.0f\n", counter_cut5);
    printf("--------------------\n");
    printf("[INFO] Individual NoCut: Yields = %7.3f\n", counter_nocut_norm);
    printf("[INFO] Individual Cut00: Yields = %7.3f\n", counter_cut0_norm);
    printf("[INFO] Individual Cut01: Yields = %7.3f\n", counter_cut1_norm);
    printf("[INFO] Individual Cut02: Yields = %7.3f\n", counter_cut2_norm);
    printf("[INFO] Individual Cut03: Yields = %7.3f\n", counter_cut3_norm);
    printf("[INFO] Individual Cut04: Yields = %7.3f\n", counter_cut4_norm);
    printf("[INFO] Individual Cut05: Yields = %7.3f\n", counter_cut5_norm);
    //}}}
    //report successive{{{
    printf("\nSuccessive selection\n");
    printf("--------------------\n");
    printf("[INFO] Successive NoCut: Entries = %.0f\n", counter_nocut_successive);
    printf("[INFO] Successive Cut00: Entries = %.0f\n", counter_cut0_successive);
    printf("[INFO] Successive Cut01: Entries = %.0f\n", counter_cut1_successive);
    printf("[INFO] Successive Cut02: Entries = %.0f\n", counter_cut2_successive);
    printf("[INFO] Successive Cut03: Entries = %.0f\n", counter_cut3_successive);
    printf("[INFO] Successive Cut04: Entries = %.0f\n", counter_cut4_successive);
    printf("[INFO] Successive Cut05: Entries = %.0f\n", counter_cut5_successive);
    printf("--------------------\n");
    printf("[INFO] Successive NoCut: Yields = %7.3f\n", counter_nocut_norm_successive);
    printf("[INFO] Successive Cut00: Yields = %7.3f\n", counter_cut0_norm_successive);
    printf("[INFO] Successive Cut01: Yields = %7.3f\n", counter_cut1_norm_successive);
    printf("[INFO] Successive Cut02: Yields = %7.3f\n", counter_cut2_norm_successive);
    printf("[INFO] Successive Cut03: Yields = %7.3f\n", counter_cut3_norm_successive);
    printf("[INFO] Successive Cut04: Yields = %7.3f\n", counter_cut4_norm_successive);
    printf("[INFO] Successive Cut05: Yields = %7.3f\n", counter_cut5_norm_successive);
    //}}}

    /*
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    hist->Draw();
    c1->SaveAs("check_num_bjets.png");
    hist_probb->SetLineColor(kBlue);
    hist_jet_probb->SetLineColor(kRed);
    hist_probb->Draw();
    hist_jet_probb->Draw("same");
    c1->SaveAs("check_probb.png");

    hist_probbb->SetLineColor(kBlue);
    hist_jet_probbb->SetLineColor(kRed);
    hist_probbb->Draw();
    hist_jet_probbb->Draw("same");
    c1->SaveAs("check_probbb.png");
    */
    return 0;
}

// functions{{{
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
//}}}
