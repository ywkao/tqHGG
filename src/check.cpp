#include <stdio.h>
#include <vector>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>
using namespace std;

int main(int argc, char *argv[]){

    printf("Hello World.\n");
    char input_file[512]; sprintf(input_file, "%s", argv[1]); printf("[INFO] input_file  = %s\n", input_file);
    TChain *flashggStdTree = new TChain("flashggNtuples/flashggStdTree");
    flashggStdTree->Add(input_file);
    //flashggStdTree->Add("/wk_cms2/youying/public/2017_94X_3_1_X_and_3_2_0/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root");
    //flashggStdTree->Add("/wk_cms2/youying/public/2017_94X_3_1_X_and_3_2_0/DoubleEG_B.root");
    //flashggStdTree->Add("/wk_cms2/youying/public/2017_94X_3_1_X_and_3_2_0/DoubleEG_C.root");
    //flashggStdTree->Add("/wk_cms2/youying/public/2017_94X_3_1_X_and_3_2_0/DoubleEG_D.root");
    //flashggStdTree->Add("/wk_cms2/youying/public/2017_94X_3_1_X_and_3_2_0/DoubleEG_E.root");
    //flashggStdTree->Add("/wk_cms2/youying/public/2017_94X_3_1_X_and_3_2_0/DoubleEG_F.root");
    //SetBranchAddresses{{{
    bool  EvtInfo_passTrigger;
    Int_t EvtInfo_NVtx;
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

    Int_t num_jets;
    Int_t num_leptons;
    /*
    TFile *fout   = new TFile("check.root", "RECREATE");
    TTree *mytree = new TTree("mytree", "mytree");
    
    //Branch{{{
    Int_t num_jets;
    Int_t num_leptons;
    std::vector<float> jet_pt;
    std::vector<float> jet_eta;
    std::vector<float> jet_phi;
    std::vector<float> jet_energy;
    std::vector<float> lepton_pt;
    std::vector<float> lepton_eta;
    std::vector<float> lepton_phi;
    std::vector<float> lepton_energy;
    //------------------------
    mytree -> Branch("num_jets", &num_jets, "num_jets/I");
    mytree -> Branch("jet_pt", &jet_pt);
    mytree -> Branch("jet_eta", &jet_eta);
    mytree -> Branch("jet_phi", &jet_phi);
    mytree -> Branch("jet_energy", &jet_energy);
    mytree -> Branch("num_leptons", &num_leptons, "num_leptons/I");
    mytree -> Branch("lepton_pt", &lepton_pt);
    mytree -> Branch("lepton_eta", &lepton_eta);
    mytree -> Branch("lepton_phi", &lepton_phi);
    mytree -> Branch("lepton_energy", &lepton_energy);
    */
    //}}}


    int nentries = flashggStdTree->GetEntries(); printf("[INFO] nentries = %d\n", nentries);

    double TotalGenweight=0;
    for(int ientry=0; ientry<nentries; ientry++){
        flashggStdTree->GetEntry(ientry);//load data
        //TotalGenweight += EvtInfo_genweight;
        TotalGenweight += 1.;
    }
    
    double Luminosity = 41.5; //fb{-1}
    double CrossSection = 0.0016; //pb
    double BranchingFraction = 1.; //pb
    double NormalizationFactor = 1000. * Luminosity * CrossSection * BranchingFraction / TotalGenweight;
    printf("NormalizationFactor = %.5f\n", NormalizationFactor);
    printf("TotalGenweight = %.5f\n", TotalGenweight);

    float counter_nocut = 0;
    float counter_cut0 = 0;
    float counter_cut1 = 0;
    float counter_cut2 = 0;
    float counter_cut3 = 0;
    float counter_cut4 = 0;


    for(int ientry=0; ientry<nentries; ientry++){
        //counter_nocut+=EvtInfo_genweight * NormalizationFactor;
        counter_nocut+=1;
        flashggStdTree->GetEntry(ientry);//load data
        //Reset{{{
        /*
        jet_pt.clear();
        jet_eta.clear();
        jet_phi.clear();
        jet_energy.clear();
        lepton_pt.clear();
        lepton_eta.clear();
        lepton_phi.clear();
        lepton_energy.clear();
        */
        //}}}
        //Store photon info{{{
        //==============================//
        //-----  Store PhotonInfo  -----//
        //==============================//
        bool pass_leadingPhotonPT = DiPhoInfo_leadPt > DiPhoInfo_mass / 2.;
        bool pass_subleadingPhotonPT = DiPhoInfo_subleadPt > DiPhoInfo_mass / 4.;
        bool pass_photon_criteria = pass_leadingPhotonPT && pass_subleadingPhotonPT;
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
                if(ientry==9740){
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
                if(ientry==9740){
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
        int num_jets = 0;
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
                if(ientry==9740){
                    printf("(%d) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f, Energy = %6.2f, Mass = %6.2f (jet)\n", i, jet.Pt(), jet.Eta(), jet.Phi(), jet.Energy(), jet.M());
                }
                if(num_leptons>0){
                    for(int j=0; j<num_leptons; j++){
                        TLorentzVector lepton = Leptons.at(j);
                        delta_R = jet.DeltaR(lepton);
                        if( delta_R<0.4 ) bool_passJetLeptonSeparation = false;
                        if(ientry==9740){
                            printf("(%d) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f, Energy = %6.2f, Mass = %6.2f, deltaR = %6.2f (lepton)\n", j, lepton.Pt(), lepton.Eta(), lepton.Phi(), lepton.Energy(), lepton.M(), delta_R);
                        }
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
            }
        }
        //}}}
        //Event selection (counters){{{
        //==========================//
        //----- EventSelection -----//
        //==========================//
        //Individual cuts
        if(pass_photon_criteria && EvtInfo_passTrigger) counter_cut0 += 1;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_leptons>0 ) counter_cut1 += 1;
        if(pass_photon_criteria && EvtInfo_passTrigger && bool_cut02 ) counter_cut2 += 1;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_jets>0 ) counter_cut3 += 1;
        //if(pass_photon_criteria && EvtInfo_passTrigger && num_jets>0 ) printf("Entry passing cut3: %d\n", ientry+1);
        if(pass_photon_criteria && EvtInfo_passTrigger && ((DiPhoInfo_mass > 100 && DiPhoInfo_mass < 120) || (DiPhoInfo_mass > 130 && DiPhoInfo_mass < 180)) ) counter_cut4 += 1;
        /*
        if(pass_photon_criteria && EvtInfo_passTrigger) counter_cut0 += EvtInfo_genweight * NormalizationFactor;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_leptons>0 ) counter_cut1 += EvtInfo_genweight * NormalizationFactor;
        if(pass_photon_criteria && EvtInfo_passTrigger && bool_cut02 ) counter_cut2 += EvtInfo_genweight * NormalizationFactor;
        if(pass_photon_criteria && EvtInfo_passTrigger && num_jets>0 ) counter_cut3 += EvtInfo_genweight * NormalizationFactor;
        if(pass_photon_criteria && EvtInfo_passTrigger && ((DiPhoInfo_mass > 100 && DiPhoInfo_mass < 120) || (DiPhoInfo_mass > 130 && DiPhoInfo_mass < 180)) ) counter_cut4 += EvtInfo_genweight * NormalizationFactor;
        */
        //}}}

        //if(pass_photon_criteria && EvtInfo_passTrigger && num_jets>0 ) mytree->Fill();//store cut3 info
    }//end of event loop.

    //mytree->Draw(">>eventList_ywk");
    //fout->Write();
    //fout->Close();

    printf("[INFO] NoCut: Entries = %7.3f\n", counter_nocut);
    printf("[INFO] Cut00: Entries = %7.3f\n", counter_cut0);
    printf("[INFO] Cut01: Entries = %7.3f\n", counter_cut1);
    printf("[INFO] Cut02: Entries = %7.3f\n", counter_cut2);
    printf("[INFO] Cut03: Entries = %7.3f\n", counter_cut3);
    printf("[INFO] Cut04: Entries = %7.3f\n", counter_cut4);
    return 0;
}
