//***************************************************************************
//
// FileName    : preselection.cpp
// Purpose     : Develop for top FCNH with H to two photons analysis
// Description : Extracting event info & Selecting objects (diphoton, leptons, jets).
// Author      : Yu-Wei Kao [ykao@cern.ch]
//
//***************************************************************************
// includes{{{
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
#include "../include/cross_section_2016.h"
#include "../include/cross_section_2017.h"
#include "../include/cross_section_2018.h"
#include "../include/preselection_criteria.h"
using namespace std;
//}}}

int main(int argc, char *argv[]){
    // I/O{{{
    //============================//
    //----- Input file names -----//
    //============================//
    char input_file[512] ; sprintf(input_file , "%s", argv[1]); printf("[INFO] input_file  = %s\n", input_file) ;
    char output_file[512]; sprintf(output_file, "%s", argv[2]); printf("[INFO] output_file = %s\n", output_file);
    char dataset[512]    ; sprintf(dataset    , "%s", argv[3]); printf("[INFO] dataset     = %s\n", dataset)    ;
    char tag[512]        ; sprintf(tag        , "%s", argv[4]); printf("[INFO] tag         = %s\n", tag)        ;
    char year[512]       ; sprintf(year       , "%s", argv[5]); printf("[INFO] year        = %s\n", year)       ;
    bool isData = isThisDataOrNot(dataset);
    bool isMCsignal = isThisMCsignal(dataset);
    bool isMultiFile = isDirectory(tag);
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
    fout->cd();
    myTreeClass mytree;
    mytree.InitTree();
    mytree.MakeNewBranchAddresses();
    //}}}
    // luminosity, cross section, total genweight{{{
    //==================================//
    //--------   Basic Info   ----------//
    //==================================//
    // Note: Normalization factors will be used in selection stage instead of preselection.
    int nentries = treeReader.GetEntries(); printf("[INFO] N_entries = %d\n", nentries);
    float Luminosity = GetLuminosity(year); //fb{-1}
    float CrossSection = GetXsec(year, dataset); //pb
    float BranchingFraction = 1.;
    printf("[INFO] Luminosity = %.2f\n", Luminosity);
    printf("[INFO] CrossSection = %f\n", CrossSection);
    printf("[INFO] Equivalent lumi. = %.2f\n", (float)nentries/CrossSection);
    printf("[INFO] BranchingFraction = %.2f\n", BranchingFraction);

    float TotalGenweight = GetTotalGenweight(year, dataset);
    // 2017old get TotalGenweight{{{
    if((string)year == "2017old"){
        float _TotalGenweight =0.;
        int nentries = treeReader.flashggStdTree->GetEntries();
        for(int ientry=0; ientry<nentries; ientry++){
            treeReader.flashggStdTree->GetEntry(ientry);//load data
            _TotalGenweight+=treeReader.GetGenWeight();
        }
        TotalGenweight = _TotalGenweight;
    }
    //}}}
    printf("[INFO] TotalGenweight = %.2f\n", TotalGenweight);

    float NormalizationFactor = 1000. * Luminosity * CrossSection * BranchingFraction / TotalGenweight;
    printf("[INFO] NormalizationFactor = %.2f\n", isData ? 1. : NormalizationFactor);
    //}}}

    //##################################################//
    //#########    Event Loop [Selection]    ###########//
    //##################################################//
    int Nevents_pass_selection = 0;
    for(int ientry=0; ientry<nentries; ientry++){
        treeReader.flashggStdTree->GetEntry(ientry); // load data
        mytree.Clear(); // reset parameters
        // # Basic selection{{{
        //==================================================//
        //--------------   Basic Selectoin   ---------------//
        //==================================================//
        //bool pass_leadingPhotonPT = treeReader.DiPhoInfo_leadPt > treeReader.DiPhoInfo_mass / 2.;
        //bool pass_subleadingPhotonPT = treeReader.DiPhoInfo_subleadPt > treeReader.DiPhoInfo_mass / 4.;
        bool pass_leadingPhotonPT = treeReader.DiPhoInfo_leadPt > CRITERION_LEADING_PHOTON_PT;
        bool pass_subleadingPhotonPT = treeReader.DiPhoInfo_subleadPt > CRITERION_SUBLEADING_PHOTON_PT;
        bool pass_photon_criteria_pt = pass_leadingPhotonPT && pass_subleadingPhotonPT;
        //--------------------------------------------------
        bool pass_leadingPhotonEta =  (treeReader.DiPhoInfo_leadEta < CRITERION_EE_EB_GAP_LOWER) || \
                                      (treeReader.DiPhoInfo_leadEta > CRITERION_EE_EB_GAP_UPPER && \
                                       treeReader.DiPhoInfo_leadEta < CRITERION_PHOTON_ETA);
        bool pass_subleadingPhotonEta = (treeReader.DiPhoInfo_subleadEta < CRITERION_EE_EB_GAP_LOWER) || \
                                        (treeReader.DiPhoInfo_subleadEta > CRITERION_EE_EB_GAP_UPPER && \
                                         treeReader.DiPhoInfo_subleadEta < CRITERION_PHOTON_ETA);
        bool pass_photon_criteria_eta = pass_leadingPhotonEta && pass_subleadingPhotonEta;
        //--------------------------------------------------
        //bool pass_interested_region = treeReader.DiPhoInfo_mass > 100 && treeReader.DiPhoInfo_mass < 180;
        bool pass_interested_region = treeReader.DiPhoInfo_mass > 100;
        bool pass_signal_region = treeReader.DiPhoInfo_mass>120 && treeReader.DiPhoInfo_mass<130;

        //require MC events pass trigger.
        if(!treeReader.EvtInfo_passTrigger) continue;
        //require the quality of photons.
        if(!pass_photon_criteria_pt) continue;
        if(!pass_photon_criteria_eta) continue;
        //control region
        if(!pass_interested_region) continue;
        if(isData && pass_signal_region) continue;
        //Others: photon id mva -> to be cut at src/selection.cpp
        //end of basic selection}}}
        // # Store event info{{{
        //=============================================//
        //-----------   Store Event Info    -----------//
        //=============================================//
        mytree.EvtInfo_totalEntry_before_preselection = nentries;
        mytree.EvtInfo_NormalizationFactor_lumi = isData ? 1. : NormalizationFactor;
        mytree.EvtInfo_NPu = treeReader.EvtInfo_NPu;
        mytree.EvtInfo_Rho = treeReader.EvtInfo_Rho;
        mytree.EvtInfo_NVtx = treeReader.EvtInfo_NVtx;
        mytree.EvtInfo_genweight = treeReader.EvtInfo_genweight;
        //end of store event info}}}
        // # Store photon info{{{
        //===============================================//
        //-----------   Store Diphoton Info   -----------//
        //===============================================//
        mytree.DiPhoInfo_leadPt              = treeReader.DiPhoInfo_leadPt;
        mytree.DiPhoInfo_leadEta             = treeReader.DiPhoInfo_leadEta;
        mytree.DiPhoInfo_leadPhi             = treeReader.DiPhoInfo_leadPhi;
        mytree.DiPhoInfo_leadE               = treeReader.DiPhoInfo_leadE;
        mytree.DiPhoInfo_leadhoe             = treeReader.DiPhoInfo_leadhoe;
        mytree.DiPhoInfo_leadIDMVA           = treeReader.DiPhoInfo_leadIDMVA;
        mytree.DiPhoInfo_leadhasPixelSeed    = treeReader.DiPhoInfo_leadhasPixelSeed;
        mytree.DiPhoInfo_subleadPt           = treeReader.DiPhoInfo_subleadPt;
        mytree.DiPhoInfo_subleadEta          = treeReader.DiPhoInfo_subleadEta;
        mytree.DiPhoInfo_subleadPhi          = treeReader.DiPhoInfo_subleadPhi;
        mytree.DiPhoInfo_subleadE            = treeReader.DiPhoInfo_subleadE;
        mytree.DiPhoInfo_subleadhoe          = treeReader.DiPhoInfo_subleadhoe;
        mytree.DiPhoInfo_subleadIDMVA        = treeReader.DiPhoInfo_subleadIDMVA;
        mytree.DiPhoInfo_subleadhasPixelSeed = treeReader.DiPhoInfo_subleadhasPixelSeed;
        mytree.DiPhoInfo_mass                = treeReader.DiPhoInfo_mass;
        mytree.DiPhoInfo_pt                  = treeReader.DiPhoInfo_pt;

        TLorentzVector leading_photon, subleading_photon, diphoton;
        leading_photon.SetPtEtaPhiE(treeReader.DiPhoInfo_leadPt, treeReader.DiPhoInfo_leadEta, treeReader.DiPhoInfo_leadPhi, treeReader.DiPhoInfo_leadE);
        subleading_photon.SetPtEtaPhiE(treeReader.DiPhoInfo_subleadPt, treeReader.DiPhoInfo_subleadEta, treeReader.DiPhoInfo_subleadPhi, treeReader.DiPhoInfo_subleadE);
        diphoton = leading_photon + subleading_photon;
        mytree.DiPhoInfo_eta = diphoton.Eta();
        mytree.DiPhoInfo_phi = diphoton.Phi();
        mytree.DiPhoInfo_energy = diphoton.Energy();
        //end of diphoton info}}}
        // # Select electrons{{{
        //==============================//
        //-----  Select Electrons  -----//
        //==============================//
        mytree.ElecInfo_Size = treeReader.ElecInfo_Size;
        std::vector<TLorentzVector> Electrons;
        bool bool_AtLeastOneElectron = treeReader.ElecInfo_Size>0 ? true : false;//treeReader.ElecInfo_Size = -999 => event without diphoton candidate
        if(bool_AtLeastOneElectron){
            for(int i=0; i<treeReader.ElecInfo_Size; i++){
                if( !treeReader.ElecInfo_EGMCutBasedIDMedium->at(i) ) continue;
                if( fabs(treeReader.ElecInfo_Eta->at(i)) > CRITERION_ELECTRON_ETA ) continue;
                if( fabs(treeReader.ElecInfo_Eta->at(i)) > CRITERION_EE_EB_GAP_LOWER && fabs(treeReader.ElecInfo_Eta->at(i)) < CRITERION_EE_EB_GAP_UPPER ) continue;
                if( fabs(treeReader.ElecInfo_Pt->at(i))  < CRITERION_ELECTRON_PT  ) continue;
                //--- check deltaR(electron,photon) ---//
                TLorentzVector electron; 
                electron.SetPtEtaPhiE(treeReader.ElecInfo_Pt->at(i), treeReader.ElecInfo_Eta->at(i), treeReader.ElecInfo_Phi->at(i), treeReader.ElecInfo_Energy->at(i));
                float delta_R = electron.DeltaR(leading_photon);
                if( delta_R < CRITERION_ELECTRON_ISO ) continue;
                delta_R = electron.DeltaR(subleading_photon);
                if( delta_R < CRITERION_ELECTRON_ISO ) continue;
                //--- calculate deltaR(electron,photon/diphoton) and store info ---//
                delta_R = electron.DeltaR(diphoton);          mytree.ElecInfo_electron_diphoton_deltaR.push_back(delta_R);
                delta_R = electron.DeltaR(leading_photon);    mytree.ElecInfo_electron_leadingPhoton_deltaR.push_back(delta_R);
                delta_R = electron.DeltaR(subleading_photon); mytree.ElecInfo_electron_subleadingPhoton_deltaR.push_back(delta_R);
                //--- store information ---//
                mytree.ElecInfo_electron_charge.push_back(treeReader.ElecInfo_Charge->at(i));
                mytree.ElecInfo_electron_pt.push_back(treeReader.ElecInfo_Pt->at(i));
                mytree.ElecInfo_electron_eta.push_back(treeReader.ElecInfo_Eta->at(i));
                mytree.ElecInfo_electron_phi.push_back(treeReader.ElecInfo_Phi->at(i));
                mytree.ElecInfo_electron_energy.push_back(treeReader.ElecInfo_Energy->at(i));
                mytree.num_electrons+=1;
                Electrons.push_back(electron);
            }
        }
        else{
                mytree.num_electrons=treeReader.ElecInfo_Size;
        }
        bool bool_AtLeastOneSelectedElectron = mytree.num_electrons>0 ? true : false;//for calculation of deltaR(e,j).
        //end of select electrons}}}
        // # Select muons{{{
        //==========================//
        //-----  Select Muons  -----//
        //==========================//
        mytree.MuonInfo_Size = treeReader.MuonInfo_Size;
        std::vector<TLorentzVector> Muons;
        bool bool_AtLeastOneMuon = treeReader.MuonInfo_Size>0 ? true : false;//treeReader.MuonInfo_Size = -999 => event without diphoton candidate
        if(bool_AtLeastOneMuon){
            for(int i=0; i<treeReader.MuonInfo_Size; i++){
                if( !treeReader.MuonInfo_CutBasedIdTight->at(i) ) continue;
                if( fabs(treeReader.MuonInfo_Eta->at(i)) > CRITERION_MUON_ETA ) continue;
                if( fabs(treeReader.MuonInfo_Pt->at(i))  < CRITERION_MUON_PT  ) continue;
                if( treeReader.MuonInfo_PFIsoDeltaBetaCorrR04->at(i) > CRITERION_MUON_RELATIVE_PFISO_R04  ) continue;
                //--- check deltaR(muon,photon) ---//
                TLorentzVector muon; 
                muon.SetPtEtaPhiE(treeReader.MuonInfo_Pt->at(i), treeReader.MuonInfo_Eta->at(i), treeReader.MuonInfo_Phi->at(i), treeReader.MuonInfo_Energy->at(i));
                float delta_R = muon.DeltaR(leading_photon);
                if( delta_R < CRITERION_MUON_ISO ) continue;
                delta_R = muon.DeltaR(subleading_photon);
                if( delta_R < CRITERION_MUON_ISO ) continue;
                //--- calculate deltaR(muon,diphoton) and store info ---//
                delta_R = muon.DeltaR(diphoton);          mytree.MuonInfo_muon_diphoton_deltaR.push_back(delta_R);
                delta_R = muon.DeltaR(leading_photon);    mytree.MuonInfo_muon_leadingPhoton_deltaR.push_back(delta_R);
                delta_R = muon.DeltaR(subleading_photon); mytree.MuonInfo_muon_subleadingPhoton_deltaR.push_back(delta_R);
                //--- store information ---//
                mytree.MuonInfo_muon_charge.push_back(treeReader.MuonInfo_Charge->at(i));
                mytree.MuonInfo_muon_pt.push_back(treeReader.MuonInfo_Pt->at(i));
                mytree.MuonInfo_muon_eta.push_back(treeReader.MuonInfo_Eta->at(i));
                mytree.MuonInfo_muon_phi.push_back(treeReader.MuonInfo_Phi->at(i));
                mytree.MuonInfo_muon_energy.push_back(treeReader.MuonInfo_Energy->at(i));
                
                mytree.num_muons+=1;
                Muons.push_back(muon);
            }
        }
        else{
                mytree.num_muons=treeReader.MuonInfo_Size;
        }
        bool bool_AtLeastOneSelectedMuon = mytree.num_muons>0 ? true : false;//for calculation of deltaR(mu,j).
        mytree.num_leptons = mytree.num_electrons + mytree.num_muons;
        //end of select muons}}}
        // # Select jets{{{
        //=========================//
        //-----  Select Jets  -----//
        //=========================//
        mytree.jets_size = treeReader.jets_size;
        std::vector<TLorentzVector> Jets;
        bool bool_AtLeastOneJet = treeReader.jets_size>0 ? true : false;//treeReader.jets_size = -999 => event without diphoton candidate
        if(bool_AtLeastOneJet){
            for(int i=0; i<treeReader.jets_size; i++){
                //flashgg::Tight2017 jet ID #Already applied flashgg package.
                if( fabs(treeReader.JetInfo_Eta->at(i)) > CRITERION_JET_ETA ) continue;
                if( fabs(treeReader.JetInfo_Pt->at(i))  < CRITERION_JET_PT  ) continue;
                //--- check deltaR(jet,photon) ---//
                TLorentzVector jet; 
                jet.SetPtEtaPhiE(treeReader.JetInfo_Pt->at(i), treeReader.JetInfo_Eta->at(i), treeReader.JetInfo_Phi->at(i), treeReader.JetInfo_Energy->at(i));
                float delta_R = jet.DeltaR(leading_photon);
                if( delta_R < CRITERION_JET_ISO ) continue;
                delta_R = jet.DeltaR(subleading_photon);
                if( delta_R < CRITERION_JET_ISO ) continue;
                bool bool_passJetLeptonSeparation = true;//if no leptons selected, the jet pass the delta_R criterion automatically.
                //--- check deltaR(jet,e) ---//
                if(bool_AtLeastOneSelectedElectron){
                    for(int i=0; i<mytree.num_electrons; i++){
                        delta_R = jet.DeltaR(Electrons.at(i));
                        if( delta_R < CRITERION_JET_ISO ) bool_passJetLeptonSeparation = false;
                    }
                }
                if(!bool_passJetLeptonSeparation) continue;//if not pass, reject the jet.
                //--- check deltaR(jet,mu) ---//
                if(bool_AtLeastOneSelectedMuon){
                    for(int i=0; i<mytree.num_muons; i++){
                        delta_R = jet.DeltaR(Muons.at(i));
                        if( delta_R < CRITERION_JET_ISO ) bool_passJetLeptonSeparation = false;
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
                mytree.JetInfo_jet_energy.push_back(treeReader.JetInfo_Energy->at(i));
                
                mytree.JetInfo_jet_pfDeepCSVJetTags_probb.push_back(treeReader.JetInfo_pfDeepCSVJetTags_probb->at(i));
                mytree.JetInfo_jet_pfDeepCSVJetTags_probbb.push_back(treeReader.JetInfo_pfDeepCSVJetTags_probbb->at(i));
                mytree.JetInfo_jet_pfDeepCSVJetTags_probc.push_back(treeReader.JetInfo_pfDeepCSVJetTags_probc->at(i));
                mytree.JetInfo_jet_pfDeepCSVJetTags_probudsg.push_back(treeReader.JetInfo_pfDeepCSVJetTags_probudsg->at(i));
                mytree.num_jets+=1;
                Jets.push_back(jet);
            }
        }
        else{
                mytree.num_jets=treeReader.jets_size;
        }
        //end of select jets}}}
        // # Store met info{{{
        //============================//
        //-----  Store MET Info  -----//
        //============================//
        mytree.MetInfo_Pt = treeReader.MetInfo_Pt;
        mytree.MetInfo_Phi = treeReader.MetInfo_Phi;
        mytree.MetInfo_Px = treeReader.MetInfo_Px;
        mytree.MetInfo_Py = treeReader.MetInfo_Py;
        mytree.MetInfo_SumET = treeReader.MetInfo_SumET;
        //end of store met}}}
        // # Store GenInfo (skipped){{{ 
        //===========================//
        //-----  Store GenInfo  -----//
        //===========================//
        ////### [warnning] will take too much space!!! more trials are needed!!!
        //mytree.GenPartInfo_size = treeReader.GenPartInfo_size;
        //for(int i=0; i<treeReader.GenPartInfo_size; i++){
        //        mytree.GenPartInfo_gen_Pt.push_back(treeReader.GenPartInfo_Pt->at(i));
        //        mytree.GenPartInfo_gen_Eta.push_back(treeReader.GenPartInfo_Eta->at(i));
        //        mytree.GenPartInfo_gen_Phi.push_back(treeReader.GenPartInfo_Phi->at(i));
        //        mytree.GenPartInfo_gen_Phi.push_back(treeReader.GenPartInfo_Phi->at(i));
        //        mytree.GenPartInfo_gen_Mass.push_back(treeReader.GenPartInfo_Mass->at(i));
        //        mytree.GenPartInfo_gen_PdgID.push_back(treeReader.GenPartInfo_PdgID->at(i));
        //        mytree.GenPartInfo_gen_Status.push_back(treeReader.GenPartInfo_Status->at(i));
        //        mytree.GenPartInfo_gen_nMo.push_back(treeReader.GenPartInfo_nMo->at(i));
        //        mytree.GenPartInfo_gen_nDa.push_back(treeReader.GenPartInfo_nDa->at(i));
        //}
        //end of store gen info}}}
        // # Event counting{{{
        //==================================================//
        //-------------   Event Counting     ---------------//
        //==================================================//
        Nevents_pass_selection += 1;
        mytree.Fill();
        //}}}
    }// End of event loop.
    // report & close{{{
    //==================================================//
    //---------------------  Report  -------------------//
    //==================================================//
    printf("[INFO] Nevents_pass_selection = %d\n", Nevents_pass_selection);
    fout->Write();
    fout->Close();
    //}}}
    return 0;
}//end

// functions{{{
float GetTotalGenweight(char* year, char* dataset){
    if((string)year == "2016") return GetTotalGenweight_2016(dataset);
    if((string)year == "2017") return GetTotalGenweight_2017(dataset);
    if((string)year == "2018") return GetTotalGenweight_2018(dataset);
    return -1.;
}

float GetLuminosity(char* year){
    if((string)year == "2016") return GetLuminosity_2016();
    if((string)year == "2017") return GetLuminosity_2017();
    if((string)year == "2018") return GetLuminosity_2018();
    if((string)year == "2017old") return GetLuminosity_2017old();
    return -999.;
}

float GetXsec(char* year, char* dataset){
    if((string)year == "2016") return GetXsec_2016(dataset);
    if((string)year == "2017") return GetXsec_2017(dataset);
    if((string)year == "2018") return GetXsec_2018(dataset);
    if((string)year == "2017old") return GetXsec_2017old(dataset);
    return 1.;
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
    if((string)dataset == "DoubleEG_G") return true;
    if((string)dataset == "DoubleEG_H") return true;
    if((string)dataset == "DoubleEG") return true;
    if((string)dataset == "EGamma") return true;
    return false;
}
bool isDirectory(char* tag){
    if((string)tag == "directory") return true;
    else return false;
}
//}}}
// class functions{{{
flashggStdTreeParameters::flashggStdTreeParameters(){
    GenPartInfo_Pt = new std::vector<float>;
    GenPartInfo_Eta = new std::vector<float>;
    GenPartInfo_Phi = new std::vector<float>;
    GenPartInfo_Mass = new std::vector<float>;
    GenPartInfo_PdgID = new std::vector<int>;
    GenPartInfo_Status = new std::vector<int>;
    GenPartInfo_nMo = new std::vector<int>;
    GenPartInfo_nDa = new std::vector<int>;
    //------------------------
    JetInfo_Pt = new std::vector<float>;
    JetInfo_Eta = new std::vector<float>;
    JetInfo_Phi = new std::vector<float>;
    JetInfo_Mass = new std::vector<float>;
    JetInfo_Energy = new std::vector<float>;
    JetInfo_pfDeepCSVJetTags_probb = new std::vector<float>;
    JetInfo_pfDeepCSVJetTags_probbb = new std::vector<float>;
    JetInfo_pfDeepCSVJetTags_probc = new std::vector<float>;
    JetInfo_pfDeepCSVJetTags_probudsg = new std::vector<float>;
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
    delete GenPartInfo_Pt;
    delete GenPartInfo_Eta;
    delete GenPartInfo_Phi;
    delete GenPartInfo_Mass;
    delete GenPartInfo_PdgID;
    delete GenPartInfo_Status;
    delete GenPartInfo_nMo;
    delete GenPartInfo_nDa;
    //------------------------
    delete JetInfo_Pt;
    delete JetInfo_Eta;
    delete JetInfo_Phi;
    delete JetInfo_Mass;
    delete JetInfo_Energy;
    delete JetInfo_pfDeepCSVJetTags_probb;
    delete JetInfo_pfDeepCSVJetTags_probbb;
    delete JetInfo_pfDeepCSVJetTags_probc;
    delete JetInfo_pfDeepCSVJetTags_probudsg;
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
    //printf("[INFO] flashggStdTreeReader::GetEntries : %d\n", flashggStdTree->GetEntries());
    return flashggStdTree->GetEntries();
}
float flashggStdTreeReader::GetGenWeight(void){
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
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadhasPixelSeed", &DiPhoInfo_leadhasPixelSeed);
    //------------------------
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadPt", &DiPhoInfo_subleadPt);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadEta", &DiPhoInfo_subleadEta);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadPhi", &DiPhoInfo_subleadPhi);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadE", &DiPhoInfo_subleadE);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadhoe", &DiPhoInfo_subleadhoe);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadIDMVA", &DiPhoInfo_subleadIDMVA);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadhasPixelSeed", &DiPhoInfo_subleadhasPixelSeed);
    //------------------------
    flashggStdTree->SetBranchAddress("GenPartInfo.size", &GenPartInfo_size);
    flashggStdTree->SetBranchAddress("GenPartInfo.Pt", &GenPartInfo_Pt);
    flashggStdTree->SetBranchAddress("GenPartInfo.Eta", &GenPartInfo_Eta);
    flashggStdTree->SetBranchAddress("GenPartInfo.Phi", &GenPartInfo_Phi);
    flashggStdTree->SetBranchAddress("GenPartInfo.Mass", &GenPartInfo_Mass);
    flashggStdTree->SetBranchAddress("GenPartInfo.PdgID", &GenPartInfo_PdgID);
    flashggStdTree->SetBranchAddress("GenPartInfo.Status", &GenPartInfo_Status);
    flashggStdTree->SetBranchAddress("GenPartInfo.nMo", &GenPartInfo_nMo);
    flashggStdTree->SetBranchAddress("GenPartInfo.nDa", &GenPartInfo_nDa);
    //------------------------
    flashggStdTree->SetBranchAddress("jets_size", &jets_size);
    flashggStdTree->SetBranchAddress("JetInfo.Pt", &JetInfo_Pt);
    flashggStdTree->SetBranchAddress("JetInfo.Eta", &JetInfo_Eta);
    flashggStdTree->SetBranchAddress("JetInfo.Phi", &JetInfo_Phi);
    flashggStdTree->SetBranchAddress("JetInfo.Mass", &JetInfo_Mass);
    flashggStdTree->SetBranchAddress("JetInfo.Energy", &JetInfo_Energy);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probb", &JetInfo_pfDeepCSVJetTags_probb);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probbb", &JetInfo_pfDeepCSVJetTags_probbb);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probc", &JetInfo_pfDeepCSVJetTags_probc);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probudsg", &JetInfo_pfDeepCSVJetTags_probudsg);
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
    flashggStdTree->SetBranchAddress("MetInfo.Pt", &MetInfo_Pt);
    flashggStdTree->SetBranchAddress("MetInfo.Phi", &MetInfo_Phi);
    flashggStdTree->SetBranchAddress("MetInfo.Px", &MetInfo_Px);
    flashggStdTree->SetBranchAddress("MetInfo.Py", &MetInfo_Py);
    flashggStdTree->SetBranchAddress("MetInfo.SumET", &MetInfo_SumET);
    //------------------------
    printf("[INFO] flashggStdTreeReader::SetBranchAddresses : Finished!\n");
}


myParameters::myParameters(){
    GenPartInfo_gen_Pt_selection = new std::vector<float>;
    GenPartInfo_gen_Eta_selection = new std::vector<float>;
    GenPartInfo_gen_Phi_selection = new std::vector<float>;
    GenPartInfo_gen_Mass_selection = new std::vector<float>;
    GenPartInfo_gen_PdgID_selection = new std::vector<int>;
    GenPartInfo_gen_Status_selection = new std::vector<int>;
    GenPartInfo_gen_nMo_selection = new std::vector<int>;
    GenPartInfo_gen_nDa_selection = new std::vector<int>;
    //------------------------
    JetInfo_jet_pt_selection = new std::vector<float>;
    JetInfo_jet_eta_selection = new std::vector<float>;
    JetInfo_jet_phi_selection = new std::vector<float>;
    JetInfo_jet_energy_selection = new std::vector<float>;
    JetInfo_jet_diphoton_deltaR_selection = new std::vector<float>;
    JetInfo_jet_leadingPhoton_deltaR_selection = new std::vector<float>;
    JetInfo_jet_subleadingPhoton_deltaR_selection = new std::vector<float>;
    JetInfo_jet_pfDeepCSVJetTags_probb_selection = new std::vector<float>;
    JetInfo_jet_pfDeepCSVJetTags_probbb_selection = new std::vector<float>;
    JetInfo_jet_pfDeepCSVJetTags_probc_selection = new std::vector<float>;
    JetInfo_jet_pfDeepCSVJetTags_probudsg_selection = new std::vector<float>;
    ElecInfo_electron_charge_selection = new std::vector<int>;
    ElecInfo_electron_pt_selection = new std::vector<float>;
    ElecInfo_electron_eta_selection = new std::vector<float>;
    ElecInfo_electron_phi_selection = new std::vector<float>;
    ElecInfo_electron_energy_selection = new std::vector<float>;
    ElecInfo_electron_diphoton_deltaR_selection = new std::vector<float>;
    ElecInfo_electron_leadingPhoton_deltaR_selection = new std::vector<float>;
    ElecInfo_electron_subleadingPhoton_deltaR_selection = new std::vector<float>;
    MuonInfo_muon_charge_selection = new std::vector<int>;
    MuonInfo_muon_pt_selection = new std::vector<float>;
    MuonInfo_muon_eta_selection = new std::vector<float>;
    MuonInfo_muon_phi_selection = new std::vector<float>;
    MuonInfo_muon_energy_selection = new std::vector<float>;
    MuonInfo_muon_diphoton_deltaR_selection = new std::vector<float>;
    MuonInfo_muon_leadingPhoton_deltaR_selection = new std::vector<float>;
    MuonInfo_muon_subleadingPhoton_deltaR_selection = new std::vector<float>;
}
myParameters::~myParameters(){
    delete GenPartInfo_gen_Pt_selection;
    delete GenPartInfo_gen_Eta_selection;
    delete GenPartInfo_gen_Phi_selection;
    delete GenPartInfo_gen_Mass_selection;
    delete GenPartInfo_gen_PdgID_selection;
    delete GenPartInfo_gen_Status_selection;
    delete GenPartInfo_gen_nMo_selection;
    delete GenPartInfo_gen_nDa_selection;
    //------------------------
    delete JetInfo_jet_pt_selection;
    delete JetInfo_jet_eta_selection;
    delete JetInfo_jet_phi_selection;
    delete JetInfo_jet_energy_selection;
    delete JetInfo_jet_diphoton_deltaR_selection;
    delete JetInfo_jet_leadingPhoton_deltaR_selection;
    delete JetInfo_jet_subleadingPhoton_deltaR_selection;
    delete JetInfo_jet_pfDeepCSVJetTags_probb_selection;
    delete JetInfo_jet_pfDeepCSVJetTags_probbb_selection;
    delete JetInfo_jet_pfDeepCSVJetTags_probc_selection;
    delete JetInfo_jet_pfDeepCSVJetTags_probudsg_selection;
    delete ElecInfo_electron_charge_selection;
    delete ElecInfo_electron_pt_selection;
    delete ElecInfo_electron_eta_selection;
    delete ElecInfo_electron_phi_selection;
    delete ElecInfo_electron_energy_selection;
    delete ElecInfo_electron_diphoton_deltaR_selection;
    delete ElecInfo_electron_leadingPhoton_deltaR_selection;
    delete ElecInfo_electron_subleadingPhoton_deltaR_selection;
    delete MuonInfo_muon_charge_selection;
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
    mytree -> Branch("EvtInfo_Rho", &EvtInfo_Rho, "EvtInfo_Rho/F");
    mytree -> Branch("EvtInfo_NVtx", &EvtInfo_NVtx, "EvtInfo_NVtx/I");
    mytree -> Branch("EvtInfo_genweight", &EvtInfo_genweight, "EvtInfo_genweight/F");
    //------------------------
    //mytree -> Branch("GenPartInfo_size", &GenPartInfo_size, "GenPartInfo_size/I");
    //mytree -> Branch("GenPartInfo_gen_Pt", &GenPartInfo_gen_Pt);
    //mytree -> Branch("GenPartInfo_gen_Eta", &GenPartInfo_gen_Eta);
    //mytree -> Branch("GenPartInfo_gen_Phi", &GenPartInfo_gen_Phi);
    //mytree -> Branch("GenPartInfo_gen_Mass", &GenPartInfo_gen_Mass);
    //mytree -> Branch("GenPartInfo_gen_PdgID", &GenPartInfo_gen_PdgID);
    //mytree -> Branch("GenPartInfo_gen_Status", &GenPartInfo_gen_Status);
    //mytree -> Branch("GenPartInfo_gen_nMo", &GenPartInfo_gen_nMo);
    //mytree -> Branch("GenPartInfo_gen_nDa", &GenPartInfo_gen_nDa);
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
    mytree -> Branch("DiPhoInfo_leadhasPixelSeed", &DiPhoInfo_leadhasPixelSeed, "DiPhoInfo_leadhasPixelSeed/O");
    mytree -> Branch("DiPhoInfo_subleadPt", &DiPhoInfo_subleadPt, "DiPhoInfo_subleadPt/F");
    mytree -> Branch("DiPhoInfo_subleadEta", &DiPhoInfo_subleadEta, "DiPhoInfo_subleadEta/F");
    mytree -> Branch("DiPhoInfo_subleadPhi", &DiPhoInfo_subleadPhi, "DiPhoInfo_subleadPhi/F");
    mytree -> Branch("DiPhoInfo_subleadE", &DiPhoInfo_subleadE, "DiPhoInfo_subleadE/F");
    mytree -> Branch("DiPhoInfo_subleadhoe", &DiPhoInfo_subleadhoe, "DiPhoInfo_subleadhoe/F");
    mytree -> Branch("DiPhoInfo_subleadIDMVA", &DiPhoInfo_subleadIDMVA, "DiPhoInfo_subleadIDMVA/F");
    mytree -> Branch("DiPhoInfo_subleadhasPixelSeed", &DiPhoInfo_subleadhasPixelSeed, "DiPhoInfo_subleadhasPixelSeed/O");
    //------------------------
    mytree -> Branch("ElecInfo_Size", &ElecInfo_Size, "ElecInfo_Size/I");
    mytree -> Branch("MuonInfo_Size", &MuonInfo_Size, "MuonInfo_Size/I");
    mytree -> Branch("num_leptons", &num_leptons, "num_leptons/I");// # of selected objects.
    mytree -> Branch("num_electrons", &num_electrons, "num_electrons/I");// # of selected objects.
    mytree -> Branch("num_muons", &num_muons, "num_muons/I");// # of selected objects.
    mytree -> Branch("ElecInfo_electron_charge", &ElecInfo_electron_charge);
    mytree -> Branch("ElecInfo_electron_pt", &ElecInfo_electron_pt);
    mytree -> Branch("ElecInfo_electron_eta", &ElecInfo_electron_eta);
    mytree -> Branch("ElecInfo_electron_phi", &ElecInfo_electron_phi);
    mytree -> Branch("ElecInfo_electron_energy", &ElecInfo_electron_energy);
    mytree -> Branch("ElecInfo_electron_diphoton_deltaR", &ElecInfo_electron_diphoton_deltaR);
    mytree -> Branch("ElecInfo_electron_leadingPhoton_deltaR", &ElecInfo_electron_leadingPhoton_deltaR);
    mytree -> Branch("ElecInfo_electron_subleadingPhoton_deltaR", &ElecInfo_electron_subleadingPhoton_deltaR);
    mytree -> Branch("MuonInfo_muon_charge", &MuonInfo_muon_charge);
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
    mytree -> Branch("JetInfo_jet_pfDeepCSVJetTags_probc", &JetInfo_jet_pfDeepCSVJetTags_probc);
    mytree -> Branch("JetInfo_jet_pfDeepCSVJetTags_probudsg", &JetInfo_jet_pfDeepCSVJetTags_probudsg);
    mytree -> Branch("num_bjets", &num_bjets, "num_bjets/I");
    //------------------------
    mytree -> Branch("MetInfo_Pt", &MetInfo_Pt, "MetInfo_Pt/F");
    mytree -> Branch("MetInfo_Phi", &MetInfo_Phi, "MetInfo_Phi/F");
    mytree -> Branch("MetInfo_Px", &MetInfo_Px, "MetInfo_Px/F");
    mytree -> Branch("MetInfo_Py", &MetInfo_Py, "MetInfo_Py/F");
    mytree -> Branch("MetInfo_SumET", &MetInfo_SumET, "MetInfo_SumET/F");
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
    JetInfo_jet_pfDeepCSVJetTags_probc.clear();
    JetInfo_jet_pfDeepCSVJetTags_probudsg.clear();
    num_bjets = 0;// # of selected objects.
    //------------------------
    num_leptons = 0;
    num_electrons = 0;
    num_muons = 0;
    ElecInfo_electron_charge.clear();
    ElecInfo_electron_pt.clear();
    ElecInfo_electron_eta.clear();
    ElecInfo_electron_phi.clear();
    ElecInfo_electron_energy.clear();
    ElecInfo_electron_diphoton_deltaR.clear();
    ElecInfo_electron_leadingPhoton_deltaR.clear();
    ElecInfo_electron_subleadingPhoton_deltaR.clear();
    MuonInfo_muon_charge.clear();
    MuonInfo_muon_pt.clear();
    MuonInfo_muon_eta.clear();
    MuonInfo_muon_phi.clear();
    MuonInfo_muon_energy.clear();
    MuonInfo_muon_diphoton_deltaR.clear();
    MuonInfo_muon_leadingPhoton_deltaR.clear();
    MuonInfo_muon_subleadingPhoton_deltaR.clear();
}
//}}}
