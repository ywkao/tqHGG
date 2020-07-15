// vim: set fdm=marker:{{{
//***************************************************************************
//
// FileName    : generalChiSquareStudy_leptonic.cpp
// Purpose     : Develop top reconstruction methods in leptonic channel (test code) 
// Description : Neutrino Pz is evaluated by quadratic method and TopKinFit method
//             : the performance of these 2 methods is under investigation
// Author      : Yu-Wei Kao [ykao@cern.ch]
//
//***************************************************************************
#include "../include/generalChiSquareStudy.C"
#include "../include/generalChiSquareStudy_leptonic.C"
#include "../include/cross_section.h"
#include "../../TopKinFit/kinfit.h"
#include "../../TopKinFit/TopLep.h"
#include "../include/hist_components.cpp"
#include "../include/plotHelper.C"
using namespace std;
bool printSelectedInfo = false;
bool bool_bjet_is_loose  = true;
bool bool_bjet_is_medium = false;
bool bool_bjet_is_tight  = false;
bool bool_num_bjets_is_exactly_one = false;
bool bool_num_bjets_is_atleast_one = !bool_num_bjets_is_exactly_one;
//}}}

int main(int argc, char *argv[]){
    initHists();
    // I/O and event info{{{
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

    char dir[512] = "/wk_cms2/youying/public/tH_FCNC/Era2017_RR-31Mar2018_v2";
    //sprintf(input_file, "%s", Form("%s/ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root", dir)     ); treeReader.AddSingleRootFile(input_file) ;
    sprintf(input_file, "%s", Form("%s/ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root", dir)     ); treeReader.AddSingleRootFile(input_file) ;
    //sprintf(input_file, "%s", Form("%s/TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root", dir) ); treeReader.AddSingleRootFile(input_file) ;
    //sprintf(input_file, "%s", Form("%s/TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8.root", dir) ); treeReader.AddSingleRootFile(input_file) ;
    //sprintf(input_file, "%s", Form("%s/TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root", dir) ); treeReader.AddSingleRootFile(input_file) ;
    //sprintf(input_file, "%s", Form("%s/TT_FCNC-T2HJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8.root", dir)  ); treeReader.AddSingleRootFile(input_file) ;

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
    double CrossSection = GetXsec_2017old(dataset); //pb
    double BranchingFraction = GetBranchingFraction_2017old(dataset); //pb
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
    printf("[INFO] NormalizationFactor = %f!\n", NormalizationFactor);
    //}}}
    // # topKinFit method (skipped){{{
    /*
    KINFIT::kfit *kf = new KINFIT::kfit();
    kf->Init(TOPLEP); // Initialize tool for ttbar with FCNC top decay to Higgs(->bb)+u/c hypothesis
    kf->SetNToy(1); // Set number of toys for minimization
    // Define PDFs{{{
    std::string pdfFileName = "../TopKinFit/test/GenAnalysis/TopLep/pdf.root";
    //kf->SetPDF("TopWMass",pdfFileName.c_str(),"TopLepWM_Fit");
    kf->SetPDF("TopWMass",pdfFileName.c_str(),"TopWM_Fit");
    kf->SetPDF("TopMass",pdfFileName.c_str(),"TopLepRecM_Fit");
    //kf->SetPDF("HiggsMass",pdfFileName.c_str(),"HiggsRecM_Fit");
    //kf->SetPDF("TopHadMass",pdfFileName.c_str(),"TopHadRecM_Fit");
    kf->SetPDF("MetPx",pdfFileName.c_str(),"dMetPx_Gaus");
    kf->SetPDF("MetPy",pdfFileName.c_str(),"dMetPy_Gaus");
    kf->SetPDF("BJetPx",pdfFileName.c_str(),"dBJetPx_Fit");
    kf->SetPDF("BJetPy",pdfFileName.c_str(),"dBJetPy_Fit");
    kf->SetPDF("BJetPz",pdfFileName.c_str(),"dBJetPz_Fit");
    kf->SetPDF("BJetE",pdfFileName.c_str(),"dBJetE_Fit");
    kf->SetPDF("NonBJetPx",pdfFileName.c_str(),"dNonBJetPx_Fit");
    kf->SetPDF("NonBJetPy",pdfFileName.c_str(),"dNonBJetPy_Fit");
    kf->SetPDF("NonBJetPz",pdfFileName.c_str(),"dNonBJetPz_Fit");
    kf->SetPDF("NonBJetE",pdfFileName.c_str(),"dNonBJetE_Fit");
    kf->SetPDF("ElecPx",pdfFileName.c_str(),"dElecPx_Fit");
    kf->SetPDF("ElecPy",pdfFileName.c_str(),"dElecPy_Fit");
    kf->SetPDF("ElecPz",pdfFileName.c_str(),"dElecPz_Fit");
    kf->SetPDF("ElecE",pdfFileName.c_str(),"dElecE_Fit");
    kf->SetPDF("MuonPx",pdfFileName.c_str(),"dMuonPx_Fit");
    kf->SetPDF("MuonPy",pdfFileName.c_str(),"dMuonPy_Fit");
    kf->SetPDF("MuonPz",pdfFileName.c_str(),"dMuonPz_Fit");
    kf->SetPDF("MuonE",pdfFileName.c_str(),"dMuonE_Fit");
    //}}}
    */
    //}}}
    
    //===============================================//
    //----------------   Event Loop  ----------------//
    //===============================================//
    for(int ientry=0; ientry<nentries; ientry++){
    //for(int ientry=0; ientry<100; ientry++){
        treeReader.flashggStdTree->GetEntry(ientry);//load data
        //# Pre-process{{{
        // reset, selection, Normalization factor{{{ 
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
        //if( !isMCsignal && treeReader.DiPhoInfo_mass>120 && treeReader.DiPhoInfo_mass<130) continue;
        //require the quality of photons. (PT requirement)
        //bool pass_leadingPhotonPT = treeReader.DiPhoInfo_leadPt > treeReader.DiPhoInfo_mass / 2.;
        //bool pass_subleadingPhotonPT = treeReader.DiPhoInfo_subleadPt > treeReader.DiPhoInfo_mass / 4.;
        bool pass_leadingPhotonPT = treeReader.DiPhoInfo_leadPt > 35.;
        bool pass_subleadingPhotonPT = treeReader.DiPhoInfo_subleadPt > 25.;
        bool pass_photon_criteria_pt = pass_leadingPhotonPT && pass_subleadingPhotonPT;

        bool pass_leadingPhotonEta =  (treeReader.DiPhoInfo_leadEta < 1.4442) || (treeReader.DiPhoInfo_leadEta > 1.566 && treeReader.DiPhoInfo_leadEta < 2.4);
        bool pass_subleadingPhotonEta = (treeReader.DiPhoInfo_subleadEta < 1.4442) || (treeReader.DiPhoInfo_subleadEta > 1.566 && treeReader.DiPhoInfo_subleadEta < 2.4);
        bool pass_photon_criteria_eta = pass_leadingPhotonEta && pass_subleadingPhotonEta;

        bool pass_photon_criteria = pass_photon_criteria_pt && pass_photon_criteria_eta;
        if(!pass_photon_criteria) continue;

        bool pass_photon_IDMVA = treeReader.DiPhoInfo_leadIDMVA>-0.7 && treeReader.DiPhoInfo_subleadIDMVA>-0.7;
        if(!pass_photon_IDMVA) continue;

        //Others
        //if(treeReader.DiPhoInfo_mass<0) continue;
        //if( !(treeReader.jets_size>0) ) continue;
        //if(!(DiPhoInfo_leadIDMVA>0)) continue;
        //==================================================//
        //------------   Normalization factor   ------------//
       //==================================================//
        mytree.EvtInfo_totalEntry_before_preselection = nentries;
        mytree.EvtInfo_NormalizationFactor_lumi = isData ? 1. : NormalizationFactor;
        // }}} 
        // diphoton info{{{
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
        // }}}
        // electron selection {{{
        std::vector<TLorentzVector> Leptons;
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
                if( fabs(treeReader.ElecInfo_Pt->at(i))  < 10.  ) continue;
                //--- check deltaR(electron,photon) ---//
                TLorentzVector electron; 
                electron.SetPtEtaPhiE(treeReader.ElecInfo_Pt->at(i), treeReader.ElecInfo_Eta->at(i), treeReader.ElecInfo_Phi->at(i), treeReader.ElecInfo_Energy->at(i));
                double delta_R = electron.DeltaR(leading_photon);
                if( delta_R<0.2 ) continue;
                delta_R = electron.DeltaR(subleading_photon);
                if( delta_R<0.2 ) continue;
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
                Leptons.push_back(electron);
            }
        }
        else{
                mytree.num_electrons=treeReader.ElecInfo_Size;
        }
        bool bool_AtLeastOneSelectedElectron = mytree.num_electrons>0 ? true : false;//for calculation of deltaR(e,j).
        // }}}
        // muon selection {{{
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
                if( fabs(treeReader.MuonInfo_Pt->at(i))  < 5.  ) continue;
                if( treeReader.MuonInfo_PFIsoDeltaBetaCorrR04->at(i)  > 0.25  ) continue;
                //--- check deltaR(muon,photon) ---//
                TLorentzVector muon; 
                muon.SetPtEtaPhiE(treeReader.MuonInfo_Pt->at(i), treeReader.MuonInfo_Eta->at(i), treeReader.MuonInfo_Phi->at(i), treeReader.MuonInfo_Energy->at(i));
                double delta_R = muon.DeltaR(leading_photon);
                if( delta_R<0.2 ) continue;
                delta_R = muon.DeltaR(subleading_photon);
                if( delta_R<0.2 ) continue;
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
                Leptons.push_back(muon);
            }
        }
        else{
                mytree.num_muons=treeReader.MuonInfo_Size;
        }
        bool bool_AtLeastOneSelectedMuon = mytree.num_muons>0 ? true : false;//for calculation of deltaR(mu,j).
        mytree.num_leptons = mytree.num_electrons + mytree.num_muons;
        //printf("num_leptons = %d\n", mytree.num_leptons);
        // }}}
        // jet selection {{{
        //=========================//
        //-----  Select Jets  -----//
        //=========================//
        mytree.jets_size = treeReader.jets_size;
        std::vector<TLorentzVector> Jets;
        std::vector<int> Jets_GenFalavor;
        bool bool_AtLeastOneJet = treeReader.jets_size>0 ? true : false;//treeReader.jets_size = -999 => event without diphoton candidate
        if(bool_AtLeastOneJet){
            for(int i=0; i<treeReader.jets_size; i++){
                //flashgg::Tight2017 jet ID #Already applied flashgg package.
                if( fabs(treeReader.JetInfo_Eta->at(i)) > 2.4 ) continue;
                if( fabs(treeReader.JetInfo_Pt->at(i))  < 25.  ) continue;
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
                mytree.JetInfo_jet_energy.push_back(treeReader.JetInfo_Energy->at(i));
                
                mytree.JetInfo_jet_pfDeepCSVJetTags_probb.push_back(treeReader.JetInfo_pfDeepCSVJetTags_probb->at(i));
                mytree.JetInfo_jet_pfDeepCSVJetTags_probbb.push_back(treeReader.JetInfo_pfDeepCSVJetTags_probbb->at(i));
                mytree.num_jets+=1;
                Jets.push_back(jet);
                Jets_GenFalavor.push_back(treeReader.JetInfo_GenFlavor->at(i));
            }
        }
        else{
                mytree.num_jets=treeReader.jets_size;
        }
        // }}}
        // Check b quark info, store diphoton and jets{{{
        int num_bquark=0, num_light_quark=0;
        for(int i=0; i<treeReader.GenPartInfo_size; i++){
            if(treeReader.GenPartInfo_Status->at(i)==23 && abs(treeReader.GenPartInfo_PdgID->at(i))==5) num_bquark+=1;
            if(treeReader.GenPartInfo_Status->at(i)==23 && abs(treeReader.GenPartInfo_PdgID->at(i))<5)  num_light_quark+=1;
        }
        if(num_bquark==1) check_num_bquak_is_one += 1;
        if(num_bquark==2) check_num_bquak_is_two += 1;
        if(num_bquark==3) check_num_bquak_is_three += 1;

        hist_num_gen_bquark->Fill(num_bquark);
        hist_num_gen_light_quark->Fill(num_light_quark);
        hist_num_gen_quarks->Fill(num_bquark + num_light_quark);
        hist_num_selected_jets->Fill(mytree.num_jets);
        hist_mass_diphoton->Fill(diphoton.M());
        //}}}
        // Leptonic event selection{{{
        if(mytree.num_leptons<1) continue;
        if(mytree.num_jets<1) continue;
        //}}}
        hist_num_leptons -> Fill(mytree.num_leptons);
        //}}}
        ////check geninfo(skipped){{{
        //int _count_wboson = 0;
        //int _count_lepton = 0;
        //for(int i=0; i<treeReader.GenPartInfo_size; i++){
        //    int pdgID = treeReader.GenPartInfo_PdgID->at(i);
        //    bool isPromptFinalState = treeReader.GenPartInfo_isPromptFinalState->at(i);
        //    bool isNeutrino = (abs(pdgID) == 12 || abs(pdgID) == 14 || abs(pdgID) == 16) && isPromptFinalState;
        //    bool isChargedLepton = (abs(pdgID) == 11 || abs(pdgID) == 13 || abs(pdgID) == 15) && isPromptFinalState;
        //    bool isWboson = (abs(pdgID) == 24);
        //    bool isTop = (abs(pdgID) == 6);
        //    
        //    if(isWboson) _count_wboson += 1;
        //    if(isChargedLepton) _count_lepton += 1;
        //}

        //bool wboson_is_more_than_one = _count_wboson > 1;
        //if(wboson_is_more_than_one){
        //    counter_events_with_two_wbosons += 1;
        //}

        //bool lepton_is_more_than_one = _count_lepton > 1;
        //if(lepton_is_more_than_one){
        //    counter_events_with_morethan_genleptons += 1;
        //    for(int i=0; i<treeReader.GenPartInfo_size; i++){
        //        printf("(%2d) ", i);
        //        printf("Status = %3d, ", treeReader.GenPartInfo_Status->at(i));
        //        printf("PdgID = %3d, ", treeReader.GenPartInfo_PdgID->at(i));
        //        printf("Pt = %6.2f, ", treeReader.GenPartInfo_Pt->at(i));
        //        printf("Eta = %9.2f, ", treeReader.GenPartInfo_Eta->at(i));
        //        printf("Phi = %6.2f, ", treeReader.GenPartInfo_Phi->at(i));
        //        printf("Mass = %6.2f, ", treeReader.GenPartInfo_Mass->at(i));
        //        printf("isHardProcess = %3d, ", treeReader.GenPartInfo_isHardProcess->at(i) ? 1 : 0);
        //        printf("isPromptFinalState = %3d, ", treeReader.GenPartInfo_isPromptFinalState->at(i) ? 1 : 0);
        //        printf("MomPdgID = %5d, ", treeReader.GenPartInfo_MomPdgID->at(i));
        //        printf("MomStatus = %3d\n", treeReader.GenPartInfo_MomStatus->at(i));
        //    }
        //    printf("\n");
        //}
        ////}}}


        // pick up indices of jets/leptons from hard process according to MC truth{{{
        if(printSelectedInfo){
            printf("ientry = %d/%d, ", ientry+1, nentries+1);
            printf("num_jets = %d, ", mytree.num_jets);
            printf("num_leptons = %d\n", mytree.num_leptons);
        }

        int counter_is_bquarkFromSMtop = 0, jetIndex_is_bquarkFromSMtop = -999;
        int counter_is_quarkFromFCNtop = 0, jetIndex_is_quarkFromFCNtop = -999;
        std::vector<int> jetIndex_is_bquarkFromSMtop_moreThanOne(2, -999);
        std::vector<int> jetIndex_is_quarkFromFCNtop_moreThanOne(2, -999);
        for(int i=0; i<mytree.num_jets; ++i){
            //get matched ID
            int index_gen, Matched_PdgID, Matched_MomPdgID;
            if(printSelectedInfo) kinematics_info("jet", Jets[i], i);
            obtain_gen_matched_ID(printSelectedInfo, Jets[i], index_gen, Matched_PdgID, Matched_MomPdgID, treeReader.GenPartInfo_size,\
                                  treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                                  treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
            //counting, index
            if( abs(Matched_PdgID)==5 && abs(Matched_MomPdgID)==6 ){counter_is_bquarkFromSMtop += 1; jetIndex_is_bquarkFromSMtop = i;}
            if( abs(Matched_PdgID)< 5 && abs(Matched_MomPdgID)==6 ){counter_is_quarkFromFCNtop += 1; jetIndex_is_quarkFromFCNtop = i;}

            // fix bug when bquark is matched twice
            if( abs(Matched_PdgID)==5 && abs(Matched_MomPdgID)==6 ){
                if(counter_is_bquarkFromSMtop==1) jetIndex_is_bquarkFromSMtop_moreThanOne[0] = i;
                if(counter_is_bquarkFromSMtop==2) jetIndex_is_bquarkFromSMtop_moreThanOne[1] = i;
            }

            // fix bug when q-tqh is matched twice
            if( abs(Matched_PdgID)!=5 && abs(Matched_MomPdgID)==6 ){
                if(counter_is_quarkFromFCNtop==1) jetIndex_is_quarkFromFCNtop_moreThanOne[0] = i;
                if(counter_is_quarkFromFCNtop==2) jetIndex_is_quarkFromFCNtop_moreThanOne[1] = i;
            }
        }

        int counter_momPdgID_is_wboson = 0, lepIndex_is_leptonFromSMtop = -999;
        for(int i=0; i<mytree.num_leptons; ++i){
            //get matched ID
            int index_gen, Matched_PdgID, Matched_MomPdgID;
            if(printSelectedInfo) kinematics_info("lep", Leptons[i], i);
            obtain_gen_matched_ID(printSelectedInfo, Leptons[i], index_gen, Matched_PdgID, Matched_MomPdgID, treeReader.GenPartInfo_size,\
                                  treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                                  treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
            //counting, index
            bool isChargedLepton = (abs(Matched_PdgID)==11 || abs(Matched_PdgID)==13 || abs(Matched_PdgID)==15);
            if( isChargedLepton && abs(Matched_MomPdgID)==24 ){counter_momPdgID_is_wboson += 1; lepIndex_is_leptonFromSMtop = i;}
        }

        // fix bug when bquark is matched twice{{{
        double deltaR_jet_bquark[2] = {0};
        if( counter_is_bquarkFromSMtop == 2 ){ 
            //pick the jet with smaller deltaR (jet, bquark)
            int Matched_PdgID, Matched_MomPdgID;
            deltaR_jet_bquark[0] = obtain_deltaR(printSelectedInfo, Jets[jetIndex_is_bquarkFromSMtop_moreThanOne[0]],\
                                  Matched_PdgID, Matched_MomPdgID, treeReader.GenPartInfo_size,\
                                  treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                                  treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
            deltaR_jet_bquark[1] = obtain_deltaR(printSelectedInfo, Jets[jetIndex_is_bquarkFromSMtop_moreThanOne[1]],\
                                  Matched_PdgID, Matched_MomPdgID, treeReader.GenPartInfo_size,\
                                  treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                                  treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
            jetIndex_is_bquarkFromSMtop = (deltaR_jet_bquark[0] < deltaR_jet_bquark[1]) ? jetIndex_is_bquarkFromSMtop_moreThanOne[0] : jetIndex_is_bquarkFromSMtop_moreThanOne[1];
            //printf("[check] chosen correct index of bjet: %d\n", jetIndex_is_bquarkFromSMtop);
            //jetIndex_is_bquarkFromSMtop = -999;
        }//}}}
        // fix bug when q-tqh is matched twice{{{
        double deltaR_jet_tqh_quark[2] = {0};
        if( counter_is_quarkFromFCNtop == 2 ){ 
            //pick the jet with smaller deltaR (jet, bquark)
            int Matched_PdgID, Matched_MomPdgID;
            deltaR_jet_tqh_quark[0] = obtain_deltaR(printSelectedInfo, Jets[jetIndex_is_quarkFromFCNtop_moreThanOne[0]],\
                                  Matched_PdgID, Matched_MomPdgID, treeReader.GenPartInfo_size,\
                                  treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                                  treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
            deltaR_jet_tqh_quark[1] = obtain_deltaR(printSelectedInfo, Jets[jetIndex_is_quarkFromFCNtop_moreThanOne[1]],\
                                  Matched_PdgID, Matched_MomPdgID, treeReader.GenPartInfo_size,\
                                  treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                                  treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
            jetIndex_is_quarkFromFCNtop = (deltaR_jet_tqh_quark[0] < deltaR_jet_tqh_quark[1]) ? jetIndex_is_quarkFromFCNtop_moreThanOne[0] : jetIndex_is_quarkFromFCNtop_moreThanOne[1];
        }//}}}
        
        bool tbwCanBeReconstructed = (counter_momPdgID_is_wboson >= 1) && (counter_is_bquarkFromSMtop >= 1);
        if(tbwCanBeReconstructed) counter_selectedJets_tbwCanBeReconstructed += 1;
        bool tqhCanBeReconstructed = (counter_is_quarkFromFCNtop >= 1);
        if(tqhCanBeReconstructed) counter_selectedJets_tqhCanBeReconstructed += 1;
        bool sigCanBeReconstructed = tbwCanBeReconstructed && tqhCanBeReconstructed;
        if(sigCanBeReconstructed) counter_selectedJets_sigCanBeReconstructed += 1;
        //}}}
//        //check geninfo(skipped){{{
//        printf("\n");
//        printf("[GenCheck] GenPartInfo_size = %d\n", treeReader.GenPartInfo_size);
//        for(int i=0; i<treeReader.GenPartInfo_size; i++){
//            int pdgID = treeReader.GenPartInfo_PdgID->at(i);
//            bool isPromptFinalState = treeReader.GenPartInfo_isPromptFinalState->at(i);
//            bool isNeutrino = (abs(pdgID) == 12 || abs(pdgID) == 14 || abs(pdgID) == 16) && isPromptFinalState;
//            bool isChargedLepton = (abs(pdgID) == 11 || abs(pdgID) == 13 || abs(pdgID) == 15) && isPromptFinalState;
//            bool isWboson = (abs(pdgID) == 24);
//            bool isTop = (abs(pdgID) == 6);
//
//            printf("(%d) ", i);
//            printf("Status = %3d, ", treeReader.GenPartInfo_Status->at(i));
//            printf("PdgID = %3d, ", treeReader.GenPartInfo_PdgID->at(i));
//            printf("Pt = %6.2f, ", treeReader.GenPartInfo_Pt->at(i));
//            printf("Eta = %9.2f, ", treeReader.GenPartInfo_Eta->at(i));
//            printf("Phi = %6.2f, ", treeReader.GenPartInfo_Phi->at(i));
//            printf("Mass = %6.2f, ", treeReader.GenPartInfo_Mass->at(i));
//            printf("isHardProcess = %3d, ", treeReader.GenPartInfo_isHardProcess->at(i) ? 1 : 0);
//            printf("isPromptFinalState = %3d, ", treeReader.GenPartInfo_isPromptFinalState->at(i) ? 1 : 0);
//            printf("MomPdgID = %5d, ", treeReader.GenPartInfo_MomPdgID->at(i));
//            printf("MomStatus = %3d\n", treeReader.GenPartInfo_MomStatus->at(i));
//        }
//        //}}}
//        //check index of b = -999{{{
//        if(jetIndex_is_bquarkFromSMtop == -999 && mytree.num_jets > 1)
//        {
//            TLorentzVector gen_bquark;
//            printf("[check] event = %d, num_jets = %d\n", ientry, mytree.num_jets);
//            //check geninfo(skipped){{{
//            for(int i=0; i<treeReader.GenPartInfo_size; i++){
//                int pdgID = treeReader.GenPartInfo_PdgID->at(i);
//                int momPdgID = treeReader.GenPartInfo_MomPdgID->at(i);
//                bool isPromptFinalState = treeReader.GenPartInfo_isPromptFinalState->at(i);
//                bool isNeutrino = (abs(pdgID) == 12 || abs(pdgID) == 14 || abs(pdgID) == 16) && isPromptFinalState;
//                bool isChargedLepton = (abs(pdgID) == 11 || abs(pdgID) == 13 || abs(pdgID) == 15) && isPromptFinalState;
//                bool isWboson = (abs(pdgID) == 24);
//                bool isTop = (abs(pdgID) == 6);
//                bool isBquarkFromTop = (abs(pdgID) == 5) && (abs(momPdgID) == 6);
//                bool isLquarkFromTop = (abs(pdgID) < 5) && (abs(momPdgID) == 6);
//
//                if(isBquarkFromTop || isLquarkFromTop)
//                {
//                    gen_bquark.SetPtEtaPhiM(treeReader.GenPartInfo_Pt->at(i), treeReader.GenPartInfo_Eta->at(i), treeReader.GenPartInfo_Phi->at(i), treeReader.GenPartInfo_Mass->at(i));
//                    printf("(%2d) ", i);
//                    printf("Status = %3d, ", treeReader.GenPartInfo_Status->at(i));
//                    printf("PdgID = %3d, ", treeReader.GenPartInfo_PdgID->at(i));
//                    printf("Pt = %6.2f, ", treeReader.GenPartInfo_Pt->at(i));
//                    printf("Eta = %9.2f, ", treeReader.GenPartInfo_Eta->at(i));
//                    printf("Phi = %6.2f, ", treeReader.GenPartInfo_Phi->at(i));
//                    printf("Mass = %6.2f, ", treeReader.GenPartInfo_Mass->at(i));
//                    printf("isHardProcess = %3d, ", treeReader.GenPartInfo_isHardProcess->at(i) ? 1 : 0);
//                    printf("isPromptFinalState = %3d, ", treeReader.GenPartInfo_isPromptFinalState->at(i) ? 1 : 0);
//                    printf("MomPdgID = %5d, ", treeReader.GenPartInfo_MomPdgID->at(i));
//                    printf("MomStatus = %3d\n", treeReader.GenPartInfo_MomStatus->at(i));
//                }
//            }
//            //}}}
//            for(int i=0; i<mytree.num_jets; ++i)
//                kinematics_info(i, Jets[i], gen_bquark);
//            printf("\n");
//        }
//        //}}}

        // b-jet{{{
        TLorentzVector bjet;
        int num_bjets = 0, index_bjet = -999;
        int num_bjets_tight = 0, num_bjets_loose = 0, num_bjets_medium = 0;
        int index_leading_bjet_tight = -999, index_leading_bjet_loose = -999, index_leading_bjet_medium = -999;
        std::vector<int> indices_tight, indices_loose, indices_medium;
        std::vector<TLorentzVector> bjets_tight, bjets_loose, bjets_medium;
        std::vector<double> btag_scores;
        //categorize bjet according to deepCSV score{{{
        for(int i=0; i<mytree.num_jets; ++i){
            double btag_score = mytree.JetInfo_jet_pfDeepCSVJetTags_probb.at(i)+mytree.JetInfo_jet_pfDeepCSVJetTags_probbb.at(i);
            btag_scores.push_back(btag_score);

            if(btag_score >= pfDeepCSVJetTags_loose){
                num_bjets_loose += 1;
                indices_loose.push_back(i);
                bjets_loose.push_back(Jets[i]);
                if(num_bjets_loose == 1) index_leading_bjet_loose = i;
            }
            if(btag_score >= pfDeepCSVJetTags_medium){
                num_bjets_medium += 1;
                indices_medium.push_back(i);
                bjets_medium.push_back(Jets[i]);
                if(num_bjets_medium == 1) index_leading_bjet_medium = i;
            }
            if(btag_score >= pfDeepCSVJetTags_tight){
                num_bjets_tight += 1;
                indices_tight.push_back(i);
                bjets_tight.push_back(Jets[i]);
                if(num_bjets_tight == 1) index_leading_bjet_tight = i;
            }
        }//end of looping jets
        //}}}
        hist_num_bjets_loose->Fill(num_bjets_loose);
        hist_num_bjets_medium->Fill(num_bjets_medium);
        hist_num_bjets_tight->Fill(num_bjets_tight);
        //counting correct rate in when num_bjet==1 {{{
        if(num_bjets_loose==1){
            bool is_b_quark = CheckBJetID(bjets_loose[0], treeReader.GenPartInfo_size,\
                              treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                              treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
            if(is_b_quark){
                count_bjet_is_bquark_loose += 1;
                hist_num_bjets_loose_matched -> Fill(num_bjets_loose);
            }
        }
        if(num_bjets_medium==1){
            bool is_b_quark = CheckBJetID(bjets_medium[0], treeReader.GenPartInfo_size,\
                              treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                              treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
            if(is_b_quark){
                count_bjet_is_bquark_medium += 1;
                hist_num_bjets_medium_matched -> Fill(num_bjets_medium);
            }
        }
        if(num_bjets_tight==1){
            bool is_b_quark = CheckBJetID(bjets_tight[0], treeReader.GenPartInfo_size,\
                              treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                              treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
            if(is_b_quark){
                count_bjet_is_bquark_tight += 1;
                hist_num_bjets_tight_matched -> Fill(num_bjets_tight);
            }
        }
        //}}}
        //determine leading bjet according to chosen WP{{{
        if(bool_bjet_is_loose){
            index_bjet = index_leading_bjet_loose;
            if(index_bjet != -999){
                bjet = bjets_loose[0];
                num_bjets = num_bjets_loose;
            }
        }
        if(bool_bjet_is_medium){
            index_bjet = index_leading_bjet_medium;
            if(index_bjet != -999){
                bjet = bjets_medium[0];
                num_bjets = num_bjets_medium;
            }
        }
        if(bool_bjet_is_tight){
            index_bjet = index_leading_bjet_tight;
            if(index_bjet != -999){
                bjet = bjets_tight[0];
                num_bjets = num_bjets_tight;
            }
        }
        //}}}
        bool pass_bjets_multiplicity_selection;
        if(bool_num_bjets_is_exactly_one) pass_bjets_multiplicity_selection = num_bjets == 1;
        if(bool_num_bjets_is_atleast_one) pass_bjets_multiplicity_selection = num_bjets >= 1;
        if(!pass_bjets_multiplicity_selection) continue;
        // }}}
        //b-jet rates-WP study with gen-info{{{
        //printf("[check] loose   = "); for(int i=0; i<indices_loose.size(); ++i) printf("%d ", indices_loose[i]); printf("\n");
        //printf("[check] medium  = "); for(int i=0; i<indices_medium.size(); ++i) printf("%d ", indices_medium[i]); printf("\n");
        //printf("[check] tight   = "); for(int i=0; i<indices_tight.size(); ++i) printf("%d ", indices_tight[i]); printf("\n");
        //printf("[check] gen matched idx (tbw, tqh) = (%d, %d)\n", jetIndex_is_bquarkFromSMtop, jetIndex_is_quarkFromFCNtop);
        for(std::size_t i=0; i!=indices_loose.size(); ++i){
            if(indices_loose[0] == jetIndex_is_bquarkFromSMtop){ hist_rate_leading_bjet_is_bquark->Fill(0); }
            if(indices_loose[i] == jetIndex_is_bquarkFromSMtop){ hist_rate_bquark_is_included->Fill(0); break; }
        }
        for(std::size_t i=0; i!=indices_medium.size(); ++i){
            if(indices_medium[0] == jetIndex_is_bquarkFromSMtop){ hist_rate_leading_bjet_is_bquark->Fill(1); }
            if(indices_medium[i] == jetIndex_is_bquarkFromSMtop){ hist_rate_bquark_is_included->Fill(1); break; }
        }
        for(std::size_t i=0; i!=indices_tight.size(); ++i){
            if(indices_tight[0] == jetIndex_is_bquarkFromSMtop){ hist_rate_leading_bjet_is_bquark->Fill(2); }
            if(indices_tight[i] == jetIndex_is_bquarkFromSMtop){ hist_rate_bquark_is_included->Fill(2); break; }
        }
        //}}}
        
        int counter_num_genLepton = 0;
        //geninfo: extract 4-momentum{{{
        int mom_pdgID_bquark; // used to identified t(bw) and t(qH) in a signal event
        TLorentzVector chargedLepton, neutrino, wboson, bquark, topquark, antitopquark, topquark_tbw, topquark_tqh;
        for(int i=0; i<treeReader.GenPartInfo_size; i++){
            int pdgID = treeReader.GenPartInfo_PdgID->at(i);
            int mom_pdgID = treeReader.GenPartInfo_MomPdgID->at(i);
            bool isPromptFinalState = treeReader.GenPartInfo_isPromptFinalState->at(i);
            bool isNeutrino = (abs(pdgID) == 12 || abs(pdgID) == 14 || abs(pdgID) == 16) && isPromptFinalState;
            bool isChargedLepton = (abs(pdgID) == 11 || abs(pdgID) == 13 || abs(pdgID) == 15) && isPromptFinalState;
            bool isWboson = (abs(pdgID) == 24);
            bool isbquark = (abs(pdgID) == 5 && abs(mom_pdgID) == 6);
            bool isTopQuark = ( pdgID == 6 );
            bool isAntiTopQuark = ( pdgID == 6 );
            if(isNeutrino){
                neutrino.SetPtEtaPhiM(treeReader.GenPartInfo_Pt->at(i), treeReader.GenPartInfo_Eta->at(i), treeReader.GenPartInfo_Phi->at(i), treeReader.GenPartInfo_Mass->at(i));
                //kinematics_info("v", neutrino);
            }
            if(isChargedLepton){
                chargedLepton.SetPtEtaPhiM(treeReader.GenPartInfo_Pt->at(i), treeReader.GenPartInfo_Eta->at(i), treeReader.GenPartInfo_Phi->at(i), treeReader.GenPartInfo_Mass->at(i));
                counter_num_genLepton += 1;
                //kinematics_info("l", chargedLepton);
            }
            if(isbquark){
                mom_pdgID_bquark = mom_pdgID;
                bquark.SetPtEtaPhiM(treeReader.GenPartInfo_Pt->at(i), treeReader.GenPartInfo_Eta->at(i), treeReader.GenPartInfo_Phi->at(i), treeReader.GenPartInfo_Mass->at(i));
                //kinematics_info("b", bquark);
            }
            if(isWboson){
                wboson.SetPtEtaPhiM(treeReader.GenPartInfo_Pt->at(i), treeReader.GenPartInfo_Eta->at(i), treeReader.GenPartInfo_Phi->at(i), treeReader.GenPartInfo_Mass->at(i));
                //kinematics_info("w", wboson);
            }
            if(isTopQuark){
                topquark.SetPtEtaPhiM(treeReader.GenPartInfo_Pt->at(i), treeReader.GenPartInfo_Eta->at(i), treeReader.GenPartInfo_Phi->at(i), treeReader.GenPartInfo_Mass->at(i));
                //kinematics_info("t", topquark, i);
            }
            if(isAntiTopQuark){
                antitopquark.SetPtEtaPhiM(treeReader.GenPartInfo_Pt->at(i), treeReader.GenPartInfo_Eta->at(i), treeReader.GenPartInfo_Phi->at(i), treeReader.GenPartInfo_Mass->at(i));
                //kinematics_info("t", topquark, i);
            }
        }
        if(mom_pdgID_bquark ==  6){ topquark_tbw = topquark    ; topquark_tqh = antitopquark; }
        if(mom_pdgID_bquark == -6){ topquark_tbw = antitopquark; topquark_tqh = topquark    ; }
        //}}}
        hist_num_genLeptons -> Fill(counter_num_genLepton);

        // init reco leading lepton and MET{{{
        TLorentzVector lepton = Leptons[0]; // leading lepton
        float met_pt    = treeReader.MetInfo_Pt;
        float met_phi   = treeReader.MetInfo_Phi;
        float met_px    = treeReader.MetInfo_Px;
        float met_py    = treeReader.MetInfo_Py;
        float met_sumET = treeReader.MetInfo_SumET;
        vector<double> met_info = { met_pt, met_px, met_py };
        //double neutrino_pz = evaluate_neutrino_pz(lepton, met_info); // quadratic method
        double met_pz_solution_2 = evaluate_neutrino_pz(lepton, met_info); // quadratic method
        //}}}
        //--- quality of the leading lepton ---//
        double deltaR = chargedLepton.DeltaR(lepton);
        bool is_the_leading_leptong_good = deltaR < 0.01;
        if(is_the_leading_leptong_good) counter_the_leading_lepton_isGood += 1;

        TLorentzVector reco_neutrino = derive_reco_neutrino(lepton, met_info); // quadratic method
        TLorentzVector reco_wboson   = derive_reco_wboson(lepton, reco_neutrino);

        // get indices{{{
        // Question: How should a bjet be determined?
        // Method 1: leading bjet
        int index_bjet_method_1 = index_bjet;
        // Method 2: highest b-tagged score
        // --- test purporse ---//
        //int index_bjet_method_2 = std::max_element(btag_scores.begin(), btag_scores.end()) - btag_scores.begin();
        int index_bjet_method_2 = jetIndex_is_bquarkFromSMtop;

        // Method 3: top info

        vector<double> chi2_masses;
        //printf("[check] event = %-4d: Jets.size() = %d, indices_loose.size() = %d; ", Jets.size(), ientry, indices_loose.size());
        for(size_t i=0; i!=indices_loose.size(); ++i)
        {
            TLorentzVector bjet = Jets[indices_loose[i]];
            TLorentzVector reco_tbw = derive_reco_tbw(reco_wboson, bjet);
            double reco_mass = reco_tbw.M();
            double chi2_mass = (reco_mass - top_quark_mass)*(reco_mass - top_quark_mass);
            chi2_masses.push_back(chi2_mass);
            //printf("[check] (%d) chi2 = %f\n", i, chi2_mass);
        }

        int smallest_chi2_index = std::min_element(chi2_masses.begin(),chi2_masses.end()) - chi2_masses.begin();
        double smallest_chi2 = *std::min_element(chi2_masses.begin(),chi2_masses.end());
        //printf("[check] chosen: (%d) chi2 = %f\n", smallest_chi2_index, smallest_chi2);


        int index_bjet_method_3 = indices_loose[smallest_chi2_index];
        //printf("[result] true, method-1, method-2, method-3 = ");
        //printf("%d, ", jetIndex_is_bquarkFromSMtop);
        //printf("%d, ", index_bjet_method_1);
        //printf("%d, ", index_bjet_method_2);
        //printf("%d\n", index_bjet_method_3);

        if(isMatched(index_bjet_method_1, jetIndex_is_bquarkFromSMtop)) counter_index_bjet_method_1 += 1;
        if(isMatched(index_bjet_method_2, jetIndex_is_bquarkFromSMtop)) counter_index_bjet_method_2 += 1;
        if(isMatched(index_bjet_method_3, jetIndex_is_bquarkFromSMtop)) counter_index_bjet_method_3 += 1;

        if(jetIndex_is_bquarkFromSMtop==-999) counter_no_gen_matched_bjet += 1;


        int index_q = get_q_index_min_chi2(Jets, index_bjet_method_2, diphoton);
        if(isMatched(index_q, jetIndex_is_quarkFromFCNtop)) counter_index_qjet_matched += 1;
        //}}}

        bool bool_bjet_is_matched = isMatched(index_bjet_method_2, jetIndex_is_bquarkFromSMtop);
        bool bool_qjet_is_matched = isMatched(index_q, jetIndex_is_quarkFromFCNtop);
        if(bool_bjet_is_matched && is_the_leading_leptong_good) counter_tbw_matched += 1;
        if(bool_bjet_is_matched && bool_qjet_is_matched && is_the_leading_leptong_good) counter_overall_matched += 1;

        /*
        double deltaR;
        //--- M1{{{
        double M1;
        int index_q;
        std::vector<int> index_jet_chi2_modified(2, -999);
        TLorentzVector jet_q;
        TLorentzVector top_fcnh = GetBestM1(M1, mytree.num_jets, index_bjet, index_jet_chi2_modified, diphoton, Jets, index_q, jet_q);
        bool tqhIsCorrectlyMatched = check_tqhIsCorrectlyMatched(index_q, jetIndex_is_quarkFromFCNtop);
        bool bjetIsCorrectOne = index_bjet == jetIndex_is_bquarkFromSMtop;

        deltaR = chargedLepton.DeltaR(lepton);
        if(deltaR<0.2) counter_lepton_is_correctly_chosen += 1;
        if(bjetIsCorrectOne) counter_bjet_is_correct += 1;
        if(bjetIsCorrectOne && deltaR<0.2) counter_tbw_is_correct += 1;
        if(tqhIsCorrectlyMatched) counter_tqhIsCorrectlyMatched += 1;
        if(tqhIsCorrectlyMatched && bjetIsCorrectOne && deltaR<0.2) counter_signal_is_correct += 1;

        if(M1>0){
            hf_tqh_quadratic.Fill_hist(top_fcnh, topquark_tqh);
        }

        if(jetIndex_is_quarkFromFCNtop != -999){
            TLorentzVector tqh_matched = diphoton + Jets[jetIndex_is_quarkFromFCNtop];
            hf_tqh_matched.Fill_hist(tqh_matched, topquark_tqh);
        }
        //}}}
        */

        // original(skipped){{{
        /*
        // geninfo{{{
        TLorentzVector chargedLepton, neutrino, wboson, bquark, topquark_bw, topquark;
        for(int i=0; i<treeReader.GenPartInfo_size; i++){
            int pdgID = treeReader.GenPartInfo_PdgID->at(i);
            int mom_pdgID = treeReader.GenPartInfo_MomPdgID->at(i);
            bool isPromptFinalState = treeReader.GenPartInfo_isPromptFinalState->at(i);
            bool isNeutrino = (abs(pdgID) == 12 || abs(pdgID) == 14 || abs(pdgID) == 16) && isPromptFinalState;
            bool isChargedLepton = (abs(pdgID) == 11 || abs(pdgID) == 13 || abs(pdgID) == 15) && isPromptFinalState;
            bool isWboson = (abs(pdgID) == 24);
            bool isbquark = (abs(pdgID) == 5 && abs(mom_pdgID) == 6);
            bool isTopQuark = ( abs(pdgID) == 6 && i==3 ); //pick up sm top
            if(isNeutrino){
                neutrino.SetPtEtaPhiM(treeReader.GenPartInfo_Pt->at(i), treeReader.GenPartInfo_Eta->at(i), treeReader.GenPartInfo_Phi->at(i), treeReader.GenPartInfo_Mass->at(i));
                //kinematics_info("v", neutrino);
            }
            if(isChargedLepton){
                chargedLepton.SetPtEtaPhiM(treeReader.GenPartInfo_Pt->at(i), treeReader.GenPartInfo_Eta->at(i), treeReader.GenPartInfo_Phi->at(i), treeReader.GenPartInfo_Mass->at(i));
                //kinematics_info("l", chargedLepton);
            }
            if(isbquark){
                bquark.SetPtEtaPhiM(treeReader.GenPartInfo_Pt->at(i), treeReader.GenPartInfo_Eta->at(i), treeReader.GenPartInfo_Phi->at(i), treeReader.GenPartInfo_Mass->at(i));
                //kinematics_info("b", bquark);
            }
            if(isWboson){
                wboson.SetPtEtaPhiM(treeReader.GenPartInfo_Pt->at(i), treeReader.GenPartInfo_Eta->at(i), treeReader.GenPartInfo_Phi->at(i), treeReader.GenPartInfo_Mass->at(i));
                //kinematics_info("w", wboson);
                topquark_bw = wboson + bquark;
                //kinematics_info("t", topquark_bw, i);
                hist_mass_gen_wboson_leptonic->Fill(wboson.M());
                hist_mass_gen_topquark_leptonic->Fill(topquark_bw.M());
            }
            if(isTopQuark){
                topquark.SetPtEtaPhiM(treeReader.GenPartInfo_Pt->at(i), treeReader.GenPartInfo_Eta->at(i), treeReader.GenPartInfo_Phi->at(i), treeReader.GenPartInfo_Mass->at(i));
                //kinematics_info("t", topquark, i);
            }
        }
        //TLorentzVector wboson_check = neutrino + chargedLepton;
        //kinematics_info("W", wboson_check);
        hist_gen_neutrino_pz->Fill(neutrino.Pz());
        //}}}
        // init MET{{{
        double deltaR;
        float met_pt = treeReader.MetInfo_Pt;
        float met_phi = treeReader.MetInfo_Phi;
        float met_px = treeReader.MetInfo_Px;
        float met_py = treeReader.MetInfo_Py;
        float met_sumET = treeReader.MetInfo_SumET;
        //}}}
        //Quadratic Method (W, M2, M1){{{
        //--- solve met_pz{{{
        TLorentzVector lepton = Leptons[0]; // leading lepton
        float lepton_px = lepton.Px();
        float lepton_py = lepton.Py();
        float lepton_pz = lepton.Pz();
        float lepton_energy = lepton.E();
        //float coefficient_factor = ( w_boson_mass*w_boson_mass + 2*lepton_px*met_px + 2*lepton_py*met_py ) / (2.*lepton_energy);
        float coefficient_factor = ( 80.375*80.375 + 2.*lepton_px*met_px + 2.*lepton_py*met_py ) / (2.*lepton_energy);
        //float coefficient_A = 1. - (lepton_pz/lepton_energy)*(lepton_pz/lepton_energy);
        float coefficient_A = 1. - (lepton_pz*lepton_pz)/(lepton_energy*lepton_energy);
        float coefficient_B = 2.*coefficient_factor*lepton_pz/lepton_energy;
        float coefficient_C = met_pt*met_pt - coefficient_factor*coefficient_factor;
        float coefficient_D = coefficient_B*coefficient_B - 4.*coefficient_A*coefficient_C;
        
        //printf("[CHECK-coeff] coefficient_A = %6.2f, ", coefficient_A);
        //printf("lepton_pz = %6.2f, ", lepton_pz);
        //printf("lepton_energy = %6.2f \n", lepton_energy);

        hist_MetInfo_coeff_A -> Fill(coefficient_A);
        hist_MetInfo_coeff_B -> Fill(coefficient_B);
        hist_MetInfo_coeff_C -> Fill(coefficient_C);
        hist_MetInfo_coeff_B2A -> Fill(-coefficient_B / (2*coefficient_A));

        float met_pz_solution_1 = 0.0;
        float met_pz_solution_2 = 0.0;

        counter_coeff_D += 1;
        if(coefficient_D < 0){
            counter_coeff_D_isNegative += 1;
            //printf("[check] coefficient_D = %f\n", coefficient_D);
            met_pz_solution_1 = coefficient_B / (2.*coefficient_A);
            met_pz_solution_2 = coefficient_B / (2.*coefficient_A);
            //met_pz_solution_1 = sqrt( coefficient_C / coefficient_A);
            //met_pz_solution_2 = sqrt( coefficient_C / coefficient_A);
            hist_MetInfo_coeff_D -> Fill(-sqrt(-coefficient_D));//keep tracking negative value
            hist_MetInfo_coeff_D2A -> Fill(-sqrt(-coefficient_D) / (2*coefficient_A));
            if(neutrino.Pz() >  20.)      counter_coeff_D_isNegative_gen_largerThan20 += 1;
            else if(neutrino.Pz() > -20.) counter_coeff_D_isNegative_gen_smallerThan20 += 1;
            else                          counter_coeff_D_isNegative_gen_between20 += 1;
        } else{
            met_pz_solution_1 = (coefficient_B + TMath::Sqrt(coefficient_D))/(2.*coefficient_A);
            met_pz_solution_2 = (coefficient_B - TMath::Sqrt(coefficient_D))/(2.*coefficient_A);
            hist_MetInfo_coeff_D -> Fill(sqrt(coefficient_D));
            hist_MetInfo_coeff_D2A -> Fill(sqrt(coefficient_D) / (2*coefficient_A));
        }
        ////ordering
        float larger_pz  = (abs(met_pz_solution_1) > abs(met_pz_solution_2) ) ? met_pz_solution_1 : met_pz_solution_2;
        float smaller_pz = (abs(met_pz_solution_1) < abs(met_pz_solution_2) ) ? met_pz_solution_1 : met_pz_solution_2;
        met_pz_solution_1 = larger_pz;
        met_pz_solution_2 = smaller_pz;
        
        hist_MetInfo_Pz_solution_1 -> Fill(met_pz_solution_1);
        hist_MetInfo_Pz_solution_2 -> Fill(met_pz_solution_2);
        if(coefficient_D>=0) hist_MetInfo_Pz_solution_1_positiveD -> Fill(met_pz_solution_1);
        if(coefficient_D>=0) hist_MetInfo_Pz_solution_2_positiveD -> Fill(met_pz_solution_2);
        if(coefficient_D<0) hist_MetInfo_Pz_solution_1_negativeD -> Fill(met_pz_solution_1);
        if(coefficient_D<0) hist_MetInfo_Pz_solution_2_negativeD -> Fill(met_pz_solution_2);
        //}}}
        //--- W, M2{{{
        TLorentzVector L_b_lep = bjet;
        TLorentzVector L_met_lep[2];
        TLorentzVector L_w_lep[2];
        TLorentzVector L_bw_lep[2];

        float met_energy_solution_1 = TMath::Sqrt(met_pt*met_pt + met_pz_solution_1*met_pz_solution_1);
        float met_energy_solution_2 = TMath::Sqrt(met_pt*met_pt + met_pz_solution_2*met_pz_solution_2);
        L_met_lep[0].SetPxPyPzE( met_px, met_py, met_pz_solution_1, met_energy_solution_1 );
        L_met_lep[1].SetPxPyPzE( met_px, met_py, met_pz_solution_2, met_energy_solution_2 );

        L_w_lep[0].SetPxPyPzE(  (lepton_px + met_px), (lepton_py + met_py), (lepton_pz + met_pz_solution_1), (lepton_energy + met_energy_solution_1)  );
        L_w_lep[1].SetPxPyPzE(  (lepton_px + met_px), (lepton_py + met_py), (lepton_pz + met_pz_solution_2), (lepton_energy + met_energy_solution_2)  );

        L_bw_lep[0] = L_b_lep + L_w_lep[0];
        L_bw_lep[1] = L_b_lep + L_w_lep[1];

        hf_chargedLepton.Fill_hist(chargedLepton, lepton);
        hf_neutrino_sol1.Fill_hist(L_met_lep[0], neutrino);
        hf_neutrino_sol2.Fill_hist(L_met_lep[1], neutrino);
        if(coefficient_D>=0) hf_neutrino_sol2_positive.Fill_hist(L_met_lep[1], neutrino);
        else                 hf_neutrino_sol2_negative.Fill_hist(L_met_lep[1], neutrino);
        hf_wboson_quadratic.Fill_hist(L_w_lep[1], wboson);
        hf_top_quadratic.Fill_hist(L_bw_lep[1], topquark);

        deltaR = diphoton.DeltaR(L_bw_lep[1])  ; hist_deltaR_reco_top_higgs_leptonic->Fill(deltaR)  ;
        bool canFoundSolution_quadratic = coefficient_D>=0;


        //}}}
        //--- solve met_pz (gen study){{{
        TLorentzVector gen_lepton = chargedLepton; // leading gen_lepton
        TLorentzVector gen_neutrino = neutrino; // neutrino
        float gen_lepton_px = gen_lepton.Px();
        float gen_lepton_py = gen_lepton.Py();
        float gen_lepton_pz = gen_lepton.Pz();
        float gen_lepton_energy = gen_lepton.E();

        float gen_met_pt = gen_neutrino.Pt();
        float gen_met_phi = gen_neutrino.Phi();
        float gen_met_px = gen_neutrino.Px();
        float gen_met_py = gen_neutrino.Py();
        float gen_coefficient_factor = ( wboson.M()*wboson.M() + 2.*gen_lepton_px*gen_met_px + 2.*gen_lepton_py*gen_met_py ) / (2.*gen_lepton_energy);
        float gen_coefficient_A = 1. - (gen_lepton_pz*gen_lepton_pz)/(gen_lepton_energy*gen_lepton_energy);
        float gen_coefficient_B = 2.*gen_coefficient_factor*gen_lepton_pz/gen_lepton_energy;
        float gen_coefficient_C = gen_met_pt*gen_met_pt - gen_coefficient_factor*gen_coefficient_factor;
        float gen_coefficient_D = gen_coefficient_B*gen_coefficient_B - 4.*gen_coefficient_A*gen_coefficient_C;
        counter_coeff_D_gen += 1;
        if(gen_coefficient_D < 0){
            counter_coeff_D_isNegative_gen += 1;
            hist_MetInfo_coeff_D_gen -> Fill(-sqrt(-gen_coefficient_D));//keep tracking negative value
        } else{
            hist_MetInfo_coeff_D_gen -> Fill(sqrt(gen_coefficient_D));
        }
        //}}}
        //--- M1{{{
        double M1;
        int index_q;
        std::vector<int> index_jet_chi2_modified(2, -999);
        TLorentzVector jet_q;
        TLorentzVector top_fcnh = GetBestM1(M1, mytree.num_jets, index_bjet, index_jet_chi2_modified, diphoton, Jets, index_q, jet_q);
        bool tqhIsCorrectlyMatched = check_tqhIsCorrectlyMatched(index_q, jetIndex_is_quarkFromFCNtop);
        bool bjetIsCorrectOne = index_bjet == jetIndex_is_bquarkFromSMtop;
        if(bjetIsCorrectOne && tqhIsCorrectlyMatched) counter_tqhIsCorrectlyMatched += 1;
        if(bjetIsCorrectOne && deltaR<0.4) counter_lepton_is_correctly_chosen += 1;
        //}}}
        //}}}
        //TopKinFit Method{{{
        // set up input parameters{{{
        std::vector<float> BJetPt;
        std::vector<float> BJetEta;
        std::vector<float> BJetPhi;
        std::vector<float> BJetE;
        std::vector<float> NonBJetFilteredPt;
        std::vector<float> NonBJetFilteredEta;
        std::vector<float> NonBJetFilteredPhi;
        std::vector<float> NonBJetFilteredE;
        std::vector<float> ElectronPt;
        std::vector<float> ElectronEta;
        std::vector<float> ElectronPhi;
        std::vector<float> ElectronE;
        std::vector<float> MuonPt;
        std::vector<float> MuonEta;
        std::vector<float> MuonPhi;
        std::vector<float> MuonE;
        float MetRecPx = treeReader.MetInfo_Px;
        float MetRecPy = treeReader.MetInfo_Py;

        //exactly one tight bjet
        BJetPt.push_back(bjet.Pt());
        BJetEta.push_back(bjet.Eta());
        BJetPhi.push_back(bjet.Phi());
        BJetE.push_back(bjet.E());
        for(std::size_t i=0; i<Jets.size(); ++i){
            if(i==index_bjet) continue;
            NonBJetFilteredPt.push_back(Jets[i].Pt());
            NonBJetFilteredEta.push_back(Jets[i].Eta());
            NonBJetFilteredPhi.push_back(Jets[i].Phi());
            NonBJetFilteredE.push_back(Jets[i].E());
        }
        for(std::size_t i=0; i<Electrons.size(); ++i){
            ElectronPt.push_back(Electrons[i].Pt());
            ElectronEta.push_back(Electrons[i].Eta());
            ElectronPhi.push_back(Electrons[i].Phi());
            ElectronE.push_back(Electrons[i].E());
        }
        for(std::size_t i=0; i<Muons.size(); ++i){
            MuonPt.push_back(Muons[i].Pt());
            MuonEta.push_back(Muons[i].Eta());
            MuonPhi.push_back(Muons[i].Phi());
            MuonE.push_back(Muons[i].E());
        }

        // Pass reconstructed objects to the tool (std::vector<float> for all objects except float for Met)
        kf->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
        kf->SetNonBJet(NonBJetFilteredPt,NonBJetFilteredEta,NonBJetFilteredPhi,NonBJetFilteredE);
        kf->SetElectron(ElectronPt,ElectronEta,ElectronPhi,ElectronE);
        kf->SetMuon(MuonPt,MuonEta,MuonPhi,MuonE);
        kf->SetMet(MetRecPx,MetRecPy);
        //}}}
        //run{{{
        kf->Run(); // Run the tool
        int NPerm = kf->GetNPerm(); // Get number of permutations
        std::vector<float> NuPz, disc;

        bool isMoreThanOnePermutation = false;
        if(NPerm!=1)
        {
            counter_nperm += 1;
            isMoreThanOnePermutation = true;
            printf("[WARNNING] NPerm = %d\n", NPerm);
        }

        for(int ip=0;ip<NPerm;ip++) // Loop over permutations - already sorted in likelihood value from min to max
        {
            //printf("[INFO-kinfit] disc = %f\n", disc);
            disc.push_back( kf->GetDisc(ip) ); // Get minimized likelihood value
            // Get reconstructed neutrino
            float NuPx = kf->GetNuPx(ip,0);
            float NuPy = kf->GetNuPy(ip,0);
            NuPz.push_back( kf->GetNuPz(ip,0) );
        }
        //}}}

        if(disc[0] > 100000.) counter_irregular_disc += 1;
        hist_disc_topKinFit->Fill(disc[0]);

        //reconstruct W, top
        TLorentzVector L_met_topKinFit;
        TLorentzVector L_w_topKinFit;
        TLorentzVector L_bw_topKinFit;
        float met_pz_topKinFit = NuPz[0];
        float met_energy_topKinFit = TMath::Sqrt(met_pt*met_pt + met_pz_topKinFit*met_pz_topKinFit);
        L_met_topKinFit.SetPxPyPzE( met_px, met_py, met_pz_topKinFit, met_energy_topKinFit );
        L_w_topKinFit.SetPxPyPzE( (lepton_px + met_px), (lepton_py + met_py), (lepton_pz + met_pz_topKinFit), (lepton_energy + met_energy_topKinFit) );
        L_bw_topKinFit = bjet + L_w_topKinFit;

        hf_neutrino_topKinFit.Fill_hist(L_met_topKinFit, neutrino);
        hf_wboson_topKinFit.Fill_hist(L_w_topKinFit, wboson);
        hf_top_topKinFit.Fill_hist(L_bw_topKinFit, topquark);

        bool canFoundSolution_topKinFit = !(disc[0] > 100000.);
        bool is_reg_and_positive = canFoundSolution_topKinFit && met_pz_topKinFit >= 0;

        //}}}
        //# event counting{{{
        if(isMoreThanOnePermutation && canFoundSolution_topKinFit) counter_yy += 1;
        if(!isMoreThanOnePermutation && canFoundSolution_topKinFit) counter_ny += 1;
        if(!isMoreThanOnePermutation && !canFoundSolution_topKinFit) counter_nn += 1;
        if(isMoreThanOnePermutation && !canFoundSolution_topKinFit) counter_yn += 1;
        //if(canFoundSolution_quadratic && canFoundSolution_topKinFit) counter_yy += 1;
        //if(!canFoundSolution_quadratic && canFoundSolution_topKinFit) counter_ny += 1;
        //if(!canFoundSolution_quadratic && !canFoundSolution_topKinFit) counter_nn += 1;
        //if(canFoundSolution_quadratic && !canFoundSolution_topKinFit) counter_yn += 1;
        //}}}
        */
        //}}}
        Nevents_pass_selection += 1;
    }// End of event loop.
    //==================================================//
    //---------------------  Report  -------------------//
    //==================================================//
    //PrintCountsAndRatio("counter_events_with_two_wbosons", counter_events_with_two_wbosons, Nevents_pass_selection);
    //PrintCountsAndRatio("counter_events_with_morethan_genleptons", counter_events_with_morethan_genleptons, Nevents_pass_selection);

    PrintCountsAndRatio("Nevents_pass_selection", Nevents_pass_selection, nentries);

    PrintCountsAndRatio("counter_selectedJets_tbwCanBeReconstructed", counter_selectedJets_tbwCanBeReconstructed, Nevents_pass_selection);
    PrintCountsAndRatio("counter_selectedJets_tqhCanBeReconstructed", counter_selectedJets_tqhCanBeReconstructed, Nevents_pass_selection);
    PrintCountsAndRatio("counter_selectedJets_sigCanBeReconstructed", counter_selectedJets_sigCanBeReconstructed, Nevents_pass_selection);

    PrintCountsAndRatio("counter_index_bjet_method_1"       , counter_index_bjet_method_1       , Nevents_pass_selection);
    PrintCountsAndRatio("counter_index_bjet_method_2"       , counter_index_bjet_method_2       , Nevents_pass_selection);
    PrintCountsAndRatio("counter_index_bjet_method_3"       , counter_index_bjet_method_3       , Nevents_pass_selection);
    PrintCountsAndRatio("counter_no_gen_matched_bjet"       , counter_no_gen_matched_bjet       , Nevents_pass_selection);
    PrintCountsAndRatio("counter_the_leading_lepton_isGood" , counter_the_leading_lepton_isGood , Nevents_pass_selection);
    PrintCountsAndRatio("counter_index_qjet_matched"        , counter_index_qjet_matched        , Nevents_pass_selection);
    PrintCountsAndRatio("counter_tbw_matched"               , counter_tbw_matched           , Nevents_pass_selection);
    PrintCountsAndRatio("counter_tqh_matched"               , counter_index_qjet_matched        , Nevents_pass_selection);
    PrintCountsAndRatio("counter_overall_matched"           , counter_overall_matched           , Nevents_pass_selection);

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    ph_makePlot(c1 , hist_num_gen_bquark      , Form("%s/hist_num_gen_bquark.png"      , output_histDir));
    ph_makePlot(c1 , hist_num_gen_light_quark , Form("%s/hist_num_gen_light_quark.png" , output_histDir));
    ph_makePlot(c1 , hist_num_gen_quarks      , Form("%s/hist_num_gen_quarks.png"      , output_histDir));
    ph_makePlot(c1 , hist_num_leptons         , Form("%s/hist_num_leptons.png"         , output_histDir));
    ph_makePlot(c1 , hist_num_genLeptons      , Form("%s/hist_num_genLeptons.png"      , output_histDir));
    ph_makePlot(c1 , hist_num_selected_jets   , Form("%s/hist_num_selected_jets.png"   , output_histDir));
    ph_makePlot(c1 , hist_mass_diphoton       , Form("%s/hist_mass_diphoton.png"       , output_histDir));
    //hist_rate_bquark_is_included_atleastOne.png{{{
    //relative counting, not rate yet
    for(int i=0; i<3; ++i){
        //loop over 3 different b-tagged wp
        //the Nevents_pass_selection depends on wp because the multiplicity selection on b-jet!
        double content;
        content = hist_rate_bquark_is_included -> GetBinContent(i+1) / (double)Nevents_pass_selection;
        hist_rate_bquark_is_included -> SetBinContent(i+1, content);
        content = hist_rate_leading_bjet_is_bquark -> GetBinContent(i+1) / (double)Nevents_pass_selection;
        hist_rate_leading_bjet_is_bquark -> SetBinContent(i+1, content);
    }
    ph_bjet_study(c1, hist_rate_bquark_is_included, hist_rate_leading_bjet_is_bquark, Form("%s/hist_rate_bquark_is_included_atleastOne.png", output_histDir));
    //}}}
    //hist_num_bjets_bquark.png{{{
    string          xtitle  = "Number of particles";
    vector<TH1D*>   hists   = {hist_num_gen_bquark, hist_num_bjets_tight, hist_num_bjets_medium, hist_num_bjets_loose};
    vector<string>  styles  = { "l", "l", "l", "l" };
    vector<int>     colors  = {kRed, kBlue, kGreen+2, kGray+2};
    vector<int>     orders  = {0, 1, 2, 3};
    vector<string>  legends = {"b quark (gen-level)", "tight b-tagged jets", "medium b-tagged jets", "loose b-tagged jets"};
    ph_makePlots(c1, hists, styles, colors, orders, legends, xtitle, Form("%s/hist_num_bjets_bquark.png", output_histDir));
    //}}}


    // original(skipped){{{
    /*
    // report{{{
    PrintCountsAndRatio("counter_lepton_is_correctly_chosen", counter_lepton_is_correctly_chosen, Nevents_pass_selection);
    PrintCountsAndRatio("counter_bjet_is_correct", counter_bjet_is_correct, Nevents_pass_selection);
    PrintCountsAndRatio("counter_tbw_is_correct", counter_tbw_is_correct, Nevents_pass_selection);
    PrintCountsAndRatio("counter_tqhIsCorrectlyMatched", counter_tqhIsCorrectlyMatched, Nevents_pass_selection);
    PrintCountsAndRatio("counter_signal_is_correct", counter_signal_is_correct, Nevents_pass_selection);

    //TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    //hf_chargedLepton.Draw_all_hist(c1);
    //hf_neutrino_sol1.Draw_all_hist(c1);
    //hf_neutrino_sol2.Draw_all_hist(c1);
    //hf_neutrino_sol2_positive.Draw_all_hist(c1);
    //hf_neutrino_sol2_negative.Draw_all_hist(c1);
    //hf_wboson_quadratic.Draw_all_hist(c1);
    //hf_top_quadratic.Draw_all_hist(c1);
    //hf_tqh_quadratic.Draw_all_hist(c1);
    //hf_tqh_matched.Draw_all_hist(c1);
    //}}}
    // printf{{{
    PrintCountsAndRatio("counter_nperm", counter_nperm, Nevents_pass_selection);
    PrintCountsAndRatio("counter_yy", counter_yy, Nevents_pass_selection);
    PrintCountsAndRatio("counter_yn", counter_yn, Nevents_pass_selection);
    PrintCountsAndRatio("counter_ny", counter_ny, Nevents_pass_selection);
    PrintCountsAndRatio("counter_nn", counter_nn, Nevents_pass_selection);
    PrintCountsAndRatio("Nevents_tqhIsCorrectlyMatched", counter_tqhIsCorrectlyMatched, Nevents_pass_selection);
    PrintCountsAndRatio("Nevents_leptonIsCorrectlyMatched", counter_lepton_is_correctly_chosen, Nevents_pass_selection);
    PrintCountsAndRatio("number of irregular disc value", counter_irregular_disc, counter_coeff_D);
    PrintCountsAndRatio("number of counter_coeff_D_isNegative", counter_coeff_D_isNegative, counter_coeff_D);
    PrintCountsAndRatio("number of counter_coeff_D_isNegative_gen_largerThan20", counter_coeff_D_isNegative_gen_largerThan20, counter_coeff_D);
    PrintCountsAndRatio("number of counter_coeff_D_isNegative_gen_smallerThan20", counter_coeff_D_isNegative_gen_smallerThan20, counter_coeff_D);
    PrintCountsAndRatio("number of counter_coeff_D_isNegative_gen_between20", counter_coeff_D_isNegative_gen_between20, counter_coeff_D);
    PrintCountsAndRatio("number of counter_coeff_D_isNegative_gen", counter_coeff_D_isNegative_gen, counter_coeff_D_gen);
    //printf("\n//--- report of quadratic ---//\n");
    //printf("\n//--- report of topKinFit ---//\n");
    //}}}
    // hist_factory{{{
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    hf_chargedLepton.Draw_all_hist(c1);
    hf_neutrino_sol1.Draw_all_hist(c1);
    hf_neutrino_sol2.Draw_all_hist(c1);
    hf_neutrino_sol2_positive.Draw_all_hist(c1);
    hf_neutrino_sol2_negative.Draw_all_hist(c1);
    hf_neutrino_topKinFit.Draw_all_hist(c1);
    hf_wboson_quadratic.Draw_all_hist(c1);
    hf_wboson_topKinFit.Draw_all_hist(c1);
    hf_top_quadratic.Draw_all_hist(c1);
    hf_top_topKinFit.Draw_all_hist(c1);

    MakeComparisonPlots(output_histDir, c1, hf_neutrino_sol2, hf_neutrino_topKinFit);
    MakeComparisonPlots(output_histDir, c1, hf_wboson_quadratic, hf_wboson_topKinFit);
    MakeComparisonPlots(output_histDir, c1, hf_top_quadratic, hf_top_topKinFit);

    printf("\n//--- Correlation Factor of quadratic vs. topKinFit ---//\n");
    MakeComparison_CorrelationFactors(hf_neutrino_sol2, hf_neutrino_topKinFit);
    MakeComparison_CorrelationFactors(hf_neutrino_sol2_positive, hf_neutrino_topKinFit);
    MakeComparison_CorrelationFactors(hf_neutrino_sol2_negative, hf_neutrino_topKinFit);
    MakeComparison_CorrelationFactors(hf_wboson_quadratic, hf_wboson_topKinFit);
    MakeComparison_CorrelationFactors(hf_top_quadratic, hf_top_topKinFit);
    //}}}
    //NOTE: need to modify SaveAs ntuples_skimmed!!!
    //1-D plots with legend{{{
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.10);
    TLegend *legend = new TLegend(0.60,0.65,0.85,0.85);
    //# coeff_D at gen-level{{{
    double max;
    double scale = 1.2;
    double max_01 = hist_MetInfo_coeff_D -> GetMaximum();
    double max_02 = hist_MetInfo_coeff_D_gen  -> GetMaximum();
    max = (max_01 > max_02) ? max_01 : max_02;

    hist_MetInfo_coeff_D -> SetStats(0)    ;
    hist_MetInfo_coeff_D -> SetMaximum(max*scale);
    hist_MetInfo_coeff_D -> SetLineWidth(2)    ;
    hist_MetInfo_coeff_D -> Draw("hist")    ;
    hist_MetInfo_coeff_D_gen -> Draw("hist;same")    ;
    hist_MetInfo_coeff_D_gen -> SetLineWidth(0)    ;
    hist_MetInfo_coeff_D_gen -> SetFillStyle(3001)    ;
    hist_MetInfo_coeff_D_gen -> SetFillColor(kRed)    ;
    hist_MetInfo_coeff_D -> Draw("hist;same")    ;
    legend->Clear();
    legend->AddEntry(hist_MetInfo_coeff_D,  "Reco.", "l");
    legend->AddEntry(hist_MetInfo_coeff_D_gen,  "Gen-level", "f");
    legend->SetLineColor(0);
    legend->Draw("same");
    c1->SaveAs("ntuples_skimmed/hist_MetInfo_coeff_D_gen.png")    ;
    //}}}
    //# highlight D>0 and D<0 in the distribution of Pz solution{{{
    TLegend *legend_sol = new TLegend(0.15,0.65,0.45,0.85);
    hist_MetInfo_Pz_solution_1 -> Draw("hist")                  ;
    hist_MetInfo_Pz_solution_1 -> SetLineWidth(2)                  ;
    hist_MetInfo_Pz_solution_1_positiveD -> Draw("hist;same")                  ;
    hist_MetInfo_Pz_solution_1_positiveD -> SetLineColor(kRed)                  ;
    hist_MetInfo_Pz_solution_1_positiveD -> SetLineWidth(2)                  ;
    hist_MetInfo_Pz_solution_1_negativeD -> Draw("hist;same")                  ;
    hist_MetInfo_Pz_solution_1_negativeD -> SetLineColor(kGreen)                  ;
    hist_MetInfo_Pz_solution_1_negativeD -> SetLineWidth(2)                  ;
    legend_sol->Clear();
    legend_sol->AddEntry(hist_MetInfo_Pz_solution_1, "All", "l");
    legend_sol->AddEntry(hist_MetInfo_Pz_solution_1_positiveD, "D > 0", "l");
    legend_sol->AddEntry(hist_MetInfo_Pz_solution_1_negativeD, "D < 0", "l");
    legend_sol->SetLineColor(0);
    legend_sol->Draw("same");
    c1->SaveAs("ntuples_skimmed/hist_MetInfo_Pz_solution_1.png")                  ;
    //---
    hist_MetInfo_Pz_solution_2 -> Draw("hist")                  ;
    hist_MetInfo_Pz_solution_2 -> SetLineWidth(2)                  ;
    hist_MetInfo_Pz_solution_2_positiveD -> Draw("hist;same")                  ;
    hist_MetInfo_Pz_solution_2_positiveD -> SetLineColor(kRed)                  ;
    hist_MetInfo_Pz_solution_2_positiveD -> SetLineWidth(2)                  ;
    hist_MetInfo_Pz_solution_2_negativeD -> Draw("hist;same")                  ;
    hist_MetInfo_Pz_solution_2_negativeD -> SetLineColor(kGreen)                  ;
    hist_MetInfo_Pz_solution_2_negativeD -> SetLineWidth(2)                  ;
    legend_sol->Clear();
    legend_sol->AddEntry(hist_MetInfo_Pz_solution_2, "All", "l");
    legend_sol->AddEntry(hist_MetInfo_Pz_solution_2_positiveD, "D > 0", "l");
    legend_sol->AddEntry(hist_MetInfo_Pz_solution_2_negativeD, "D < 0", "l");
    legend_sol->SetLineColor(0);
    legend_sol->Draw("same");
    c1->SaveAs("ntuples_skimmed/hist_MetInfo_Pz_solution_2.png")                  ;
    //}}}
    //--------------------
    TLegend *legend_ratio = new TLegend(0.65,0.53,0.85,0.73);
    //}}}
    // 1-D plots{{{
    const char* dir = "ntuples_skimmed/chi2_study_leptonic_1D_plots";
    hist_mass_gen_wboson_leptonic -> Draw("hist")               ; c1->SaveAs( Form("%s/hist_mass_gen_wboson_leptonic.png", dir) )               ;
    hist_mass_gen_topquark_leptonic -> Draw("hist")             ; c1->SaveAs( Form("%s/hist_mass_gen_topquark_leptonic.png", dir) )             ;
    hist_deltaR_reco_top_higgs_leptonic -> Draw("hist")         ; c1->SaveAs( Form("%s/hist_deltaR_reco_top_higgs_leptonic.png", dir) )         ;
    //---
    hist_gen_neutrino_pz -> Draw("hist")                        ; c1->SaveAs( Form("%s/hist_gen_neutrino_pz.png", dir) )                        ;
    hist_MetInfo_Pz_solution_1 -> Draw("hist")                  ; c1->SaveAs( Form("%s/hist_MetInfo_Pz_solution_1.png", dir) )                  ;
    hist_MetInfo_Pz_solution_2 -> Draw("hist")                  ; c1->SaveAs( Form("%s/hist_MetInfo_Pz_solution_2.png", dir) )                  ;
    hist_MetInfo_coeff_A -> Draw("hist")                        ; c1->SaveAs( Form("%s/hist_MetInfo_coeff_A.png", dir) )                        ;
    hist_MetInfo_coeff_B -> Draw("hist")                        ; c1->SaveAs( Form("%s/hist_MetInfo_coeff_B.png", dir) )                        ;
    hist_MetInfo_coeff_C -> Draw("hist")                        ; c1->SaveAs( Form("%s/hist_MetInfo_coeff_C.png", dir) )                        ;
    hist_MetInfo_coeff_D -> Draw("hist")                        ; c1->SaveAs( Form("%s/hist_MetInfo_coeff_D.png", dir) )                        ;
    hist_MetInfo_coeff_D_gen -> Draw("hist")                    ; c1->SaveAs( Form("%s/hist_MetInfo_coeff_D_gen.png", dir) )                    ;
    hist_MetInfo_coeff_B2A -> Draw("hist")                      ; c1->SaveAs( Form("%s/hist_MetInfo_coeff_B2A.png", dir) )                      ;
    hist_MetInfo_coeff_D2A -> Draw("hist")                      ; c1->SaveAs( Form("%s/hist_MetInfo_coeff_D2A.png", dir) )                      ;
    //---
    hist_disc_topKinFit->GetXaxis()->SetRange(1, hist_disc_topKinFit->GetNbinsX() + 1);
    hist_disc_topKinFit -> Draw("hist"); c1 -> SaveAs( Form("%s/hist_disc_topKinFit.png", dir) )  ;
    //}}}
    */
        //}}}
    // close{{{
    printf("End!!\n");
    fout->Write();
    fout->Close();
    return 0;
    //}}}
}
