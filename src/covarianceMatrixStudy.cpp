//***************************************************************************
//
// FileName    : covarianceMatrixStudy.cc
// Purpose     : Study covariance matrix of Mjj and Mbjj.
// Description : 
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
#include "../include/covarianceMatrixStudy.h"
#include "../include/cross_section.h"
using namespace std;
double w_boson_mass = 80.379;//GeV
double top_quark_mass = 173.0;//GeV

int main(int argc, char *argv[]){
    // ### I/O {{{
    //============================//
    //----- Input file names -----//
    //============================//
    char input_file[512]; sprintf(input_file, "%s", argv[1]); printf("[INFO] input_file  = %s\n", input_file);
    char output_file[512]; sprintf(output_file, "%s", argv[2]); printf("[INFO] output_file = %s\n", output_file);
    char dataset[512]; sprintf(dataset, "%s", argv[3]); printf("[INFO] dataset     = %s\n", dataset);
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
    TH1D *hist_gen_w_all = new TH1D("hist_gen_w_all", "hist_gen_w_all; Mass [GeV]; Entries", 160, 0, 160);
    TH1D *hist_gen_w_sel = new TH1D("hist_gen_w_sel", "hist_gen_w_sel; Mass [GeV]; Entries", 160, 0, 160);
    //}}}
    //==================================================//
    //=========    Event Loop [Selection]    ===========//
    //==================================================//
    double covarianceMatrix[2][2] = {0};
    double mean[2] = {0};
    int Nevents_pass_selection = 0;
    int nentries = treeReader.GetEntries(); printf("[INFO] N_entries = %d\n", nentries);
    for(int ientry=0; ientry<nentries; ientry++){
        treeReader.flashggStdTree->GetEntry(ientry);//load data
        //if((ientry+1)%1000==0 || (ientry+1)==nentries) printf("ientry = %d\r", ientry);
        // ### Basic Selections{{{
        //==================================================//
        //-------------   Reset Parameters   ---------------//
        //==================================================//
        mytree.Clear();
        //==================================================//
        //--------------   Basic Selectoin   ---------------//
        //==================================================//
        //require MC events pass trigger.
        if(!treeReader.EvtInfo_passTrigger) continue;// require MC events pass trigger.
        //control region
        if(treeReader.DiPhoInfo_mass<100 || treeReader.DiPhoInfo_mass>180) continue;
        //if(treeReader.DiPhoInfo_mass>120 && treeReader.DiPhoInfo_mass<130) continue;
        //require the quality of photons. (PT requirement)
        bool pass_leadingPhotonPT = treeReader.DiPhoInfo_leadPt > treeReader.DiPhoInfo_mass / 2.;
        bool pass_subleadingPhotonPT = treeReader.DiPhoInfo_subleadPt > treeReader.DiPhoInfo_mass / 4.;
        bool pass_photon_criteria_pt = pass_leadingPhotonPT && pass_subleadingPhotonPT;

        bool pass_leadingPhotonEta =  (treeReader.DiPhoInfo_leadEta < 1.4442) || (treeReader.DiPhoInfo_leadEta > 1.566 && treeReader.DiPhoInfo_leadEta < 2.4);
        bool pass_subleadingPhotonEta = (treeReader.DiPhoInfo_subleadEta < 1.4442) || (treeReader.DiPhoInfo_subleadEta > 1.566 && treeReader.DiPhoInfo_subleadEta < 2.4);
        bool pass_photon_criteria_eta = pass_leadingPhotonEta && pass_subleadingPhotonEta;

        bool pass_photon_criteria = pass_photon_criteria_pt && pass_photon_criteria_eta;
        if(!pass_photon_criteria) continue;
        //====== Hadronic Channel Criteria ======//
        if(treeReader.jets_size<3) continue;//quick skimmed
        // ### photon && lepton criteria{{{
        //---  Store Diphoton Inf{{{
        TLorentzVector leading_photon, subleading_photon, diphoton;
        leading_photon.SetPtEtaPhiE(treeReader.DiPhoInfo_leadPt, treeReader.DiPhoInfo_leadEta, treeReader.DiPhoInfo_leadPhi, treeReader.DiPhoInfo_leadE);
        subleading_photon.SetPtEtaPhiE(treeReader.DiPhoInfo_subleadPt, treeReader.DiPhoInfo_subleadEta, treeReader.DiPhoInfo_subleadPhi, treeReader.DiPhoInfo_subleadE);
        diphoton = leading_photon + subleading_photon;
        //}}}
        //---  Select Electron{{{
        int num_electrons = 0;
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
                num_electrons+=1;
                Electrons.push_back(electron);
            }
        }
        else{
                num_electrons=treeReader.ElecInfo_Size;
        }
        bool bool_AtLeastOneSelectedElectron = mytree.num_electrons>0 ? true : false;//for calculation of deltaR(e,j).
        //}}}
        //---  Select Muon{{{
        int num_muons = 0, num_leptons = 0;
        std::vector<TLorentzVector> Muons;
        bool bool_AtLeastOneMuon = treeReader.MuonInfo_Size>0 ? true : false;//treeReader.MuonInfo_Size = -999 => event without diphoton candidate
        if(bool_AtLeastOneMuon){
            for(int i=0; i<treeReader.MuonInfo_Size; i++){
                if( !treeReader.MuonInfo_CutBasedIdTight->at(i) ) continue;
                if( fabs(treeReader.MuonInfo_Eta->at(i)) > 2.4 ) continue;
                if( fabs(treeReader.MuonInfo_Pt->at(i))  < 20  ) continue;
                if( treeReader.MuonInfo_PFIsoDeltaBetaCorrR04->at(i) > 0.25  ) continue;
                //--- check deltaR(muon,photon) ---//
                TLorentzVector muon; 
                muon.SetPtEtaPhiE(treeReader.MuonInfo_Pt->at(i), treeReader.MuonInfo_Eta->at(i), treeReader.MuonInfo_Phi->at(i), treeReader.MuonInfo_Energy->at(i));
                double delta_R = muon.DeltaR(leading_photon);
                if( delta_R<0.3 ) continue;
                delta_R = muon.DeltaR(subleading_photon);
                if( delta_R<0.3 ) continue;
                num_muons+=1;
                Muons.push_back(muon);
            }
        }
        else{
                num_muons=treeReader.MuonInfo_Size;
        }
        num_leptons = num_electrons + num_muons;
        bool bool_AtLeastOneSelectedMuon = mytree.num_muons>0 ? true : false;//for calculation of deltaR(mu,j).
        //}}}
        //}}}
        if(num_leptons > 0) continue;
        //}}}

        //### Start GenMatching for Jets{{{
        //========================================//
        //-----  Start GenMatching for Jets  -----//
        //========================================//
        //### GenMatching: find the gen particle (MC truth) for each jet (reconstructed). 
        //### pdgID: (1, 2, 3, 4, 5, 6) = (d, u, s, c, b, t)
        //### This is the simplest version. Identify the corresponding gen particle by selecting smallest deltaR(gen, jet).
        //### One can try to print out the info of pt, eta, phi, energy, and deltaR of jet and corresponding gen particle to see if they are matched.
        std::vector<int> index_GenParticles, GenParticles_PdgID;
        std::vector<TLorentzVector> Jets, GenParticles;
        //--- GenMatching for each reco jet ---//
        //printf("ientry = %d\n", ientry);
        for(int i=0; i<treeReader.jets_size; i++){
            //--- basic jet selections{{{
            if( fabs(treeReader.JetInfo_Eta->at(i)) > 2.4 ) continue;
            if( fabs(treeReader.JetInfo_Pt->at(i))  < 25  ) continue;
            TLorentzVector jet; jet.SetPtEtaPhiE(treeReader.JetInfo_Pt->at(i), treeReader.JetInfo_Eta->at(i), treeReader.JetInfo_Phi->at(i), treeReader.JetInfo_Energy->at(i));
            double delta_R = jet.DeltaR(leading_photon);
            if( delta_R<0.4 ) continue;
            delta_R = jet.DeltaR(subleading_photon);
            if( delta_R<0.4 ) continue;
            //--------------------------------------------------//
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
            //--------------------------------------------------//
            //}}}
            TLorentzVector truelove;//we are looking for the right genParticle to match the jet.
            int index = -1, truelove_PdgID = -999; double delta_R_min = 999.;
            for(int j=0; j<treeReader.GenPartInfo_size; j++){
                if( abs(treeReader.GenPartInfo_Eta->at(j)) > 10000 ) continue;//remove incoming particle
                if( abs(treeReader.GenPartInfo_PdgID->at(j)) == 6 ) continue;//exclude top quark
                bool isAvailable = checkAvailability(j, index_GenParticles); if(!isAvailable) continue;//the truelove of previous jets shall not become another one's truelove.
                //--------------------
                TLorentzVector genParticle;
                genParticle.SetPtEtaPhiM(treeReader.GenPartInfo_Pt->at(j), treeReader.GenPartInfo_Eta->at(j), treeReader.GenPartInfo_Phi->at(j), treeReader.GenPartInfo_Mass->at(j));
                double delta_R = jet.DeltaR(genParticle);
                //select quark & require min(deltaR)
                if( abs(treeReader.GenPartInfo_PdgID->at(j)) < 7 && delta_R < 0.4 && delta_R < delta_R_min){
                    index = j;//record the matched genParticle
                    delta_R_min = delta_R;
                    truelove = genParticle;
                    truelove_PdgID = treeReader.GenPartInfo_PdgID->at(j);
                }
            }//end of gen loop
            index_GenParticles.push_back(index);
            //store particle info
            Jets.push_back(jet);
            GenParticles.push_back(truelove);
            GenParticles_PdgID.push_back(truelove_PdgID);
            //printf("(%d) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f\n", i, jet.Pt(), jet.Eta(), jet.Phi());
            //printf("(%d) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f (PdgID:%2d, min_deltaR = %5.3f)\n", i, truelove.Pt(), truelove.Eta(), truelove.Phi(), truelove_PdgID, delta_R_min);
        }//end of jet loop
        //}}}
        //### Selection on gen-matched quarks{{{
        //--- require at least 3 quarks ---//
        int  count_quarks = 0;
        for(std::size_t i=0; i<GenParticles_PdgID.size(); ++i){
            if( abs(GenParticles_PdgID.at(i)) < 7){ count_quarks += 1; }
        }
        //if(count_quarks<3) continue; // st FCNH
        if(count_quarks<4) continue; // tt FCNH


        //--- reject event without bottom quark ---//
        bool has_bottom_quark;
        for(std::size_t i=0; i<GenParticles_PdgID.size(); ++i){
            if( abs(GenParticles_PdgID[i]) == 5 ){ has_bottom_quark = true; break; }
            else has_bottom_quark = false;
        }
        if(!has_bottom_quark) continue;
        //}}}
        //### GenMother: W boson info{{{
        //========================================//
        //------  GenMother: W boson info   ------//
        //========================================//
        //### Purpose: determine w-jets according to correct combination of gen particles.
        //### 1. calculate the mass of all combination of gen particles.
        //### 2. pickup the best one, and use corresponding jets to reconstruct w boson.
        int best_indices[2]; std::fill(std::begin(best_indices), std::end(best_indices), -999);
        double chi_min = 999;
        for(std::size_t i=0; i<GenParticles_PdgID.size(); ++i){
            if( abs(GenParticles_PdgID[i]) == 5 ) continue;
            for(std::size_t j=0; j<GenParticles_PdgID.size(); ++j){
                if(j<i) continue;
                if( abs(GenParticles_PdgID[j]) == 5 ) continue;
                TLorentzVector gen_w = GenParticles[i] + GenParticles[j];
                hist_gen_w_all->Fill(gen_w.M());
                double chi = (gen_w.M() - w_boson_mass) * (gen_w.M() - w_boson_mass);
                if(chi < chi_min){
                    chi_min = chi;
                    best_indices[0] = i;
                    best_indices[1] = j;
                }
            }
        }
        if(best_indices[0] == -999 || best_indices[1] == -999) continue;
        TLorentzVector gen_w_sel = GenParticles[best_indices[0]] + GenParticles[best_indices[1]];

        //--- check gen w info ---//
        double delta_R_gen_w = 999;
        //printf("\n[CHECK] ientry = %d\n", ientry);
        for(int j=0; j<treeReader.GenPartInfo_size; j++){
            if( abs(treeReader.GenPartInfo_Eta->at(j)) > 10000 ) continue;//remove incoming particle
            if( abs(treeReader.GenPartInfo_PdgID->at(j)) == 6 ) continue;//exclude top quark
            if( abs(treeReader.GenPartInfo_PdgID->at(j)) == 24 ){
                TLorentzVector genParticle;
                genParticle.SetPtEtaPhiM(treeReader.GenPartInfo_Pt->at(j), treeReader.GenPartInfo_Eta->at(j), treeReader.GenPartInfo_Phi->at(j), treeReader.GenPartInfo_Mass->at(j));
                delta_R_gen_w = gen_w_sel.DeltaR(genParticle);
                //printf("(%d) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f\n", j, gen_w_sel.Pt(), gen_w_sel.Eta(), gen_w_sel.Phi());
                //printf("(%d) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f (PdgID:%2d, deltaR = %5.3f)\n", j, genParticle.Pt(), genParticle.Eta(), genParticle.Phi(), treeReader.GenPartInfo_PdgID->at(j), delta_R);
            }
        }//end of gen loop

        //--- remove event with bad combination ---//
        if(delta_R_gen_w > 0.005) continue;
        if(gen_w_sel.M() < 20) continue; 
        hist_gen_w_sel->Fill(gen_w_sel.M());
        //}}}

        //--- hadronic top reconstructon ---//
        TLorentzVector bjet, bquark;
        std::vector<TLorentzVector> jets, gens;
        for(std::size_t i=0; i<GenParticles_PdgID.size(); ++i){
            if( abs(GenParticles_PdgID[i]) == 5 ){ bjet = Jets[i]; bquark = GenParticles[i]; }
            //else if( abs(GenParticles_PdgID.at(i)) < 5){ jets.push_back(Jets[i]); gens.push_back(GenParticles[i]); }
            //else continue;
        }
        jets.push_back(Jets[best_indices[0]]);
        jets.push_back(Jets[best_indices[1]]);
        gens.push_back(GenParticles[best_indices[0]]);
        gens.push_back(GenParticles[best_indices[1]]);


        TLorentzVector w_candidate_gen, top_candidate_gen;
        w_candidate_gen = gens[0] + gens[1];
        top_candidate_gen = bquark + w_candidate_gen;
        mytree.GenMatching_Mass_wboson_gen = w_candidate_gen.M();
        mytree.GenMatching_Mass_top_gen    = top_candidate_gen.M();

        TLorentzVector w_candidate_reco, top_candidate_reco;
        w_candidate_reco = jets[0] + jets[1];
        top_candidate_reco = bjet + w_candidate_reco;
        mytree.GenMatching_Mass_wboson_reco = w_candidate_reco.M();
        mytree.GenMatching_Mass_top_reco    = top_candidate_reco.M();

        //--- remove unphysical value ---//
        //### comment: caused by best_index of -999
        //if(w_candidate_reco.M()>999) continue;
        //if(top_candidate_reco.M()>999) continue;

        //==================================================//
        //-------------   Covariant Matrix   ---------------//
        //==================================================//
        //if(ientry>2390 && ientry<2398){//debug purpose
        //    printf("\n[CHECK] ientry = %d\n", ientry);
        //    printf("[CHECK] |%6.2f %6.2f|\n", covarianceMatrix[0][0], covarianceMatrix[0][1]);
        //    printf("[CHECK] |%6.2f %6.2f|\n", covarianceMatrix[1][0], covarianceMatrix[1][1]);
        //    printf("[CHECK] mean = %6.2f, sigma_prime_jj  = %6.2f\n", mean[0], sqrt(covarianceMatrix[0][0]));
        //    printf("[CHECK] mean = %6.2f, sigma_prime_bjj = %6.2f\n", mean[1], sqrt(covarianceMatrix[1][1]));
        //    printf("[CHECK] best_index = %d\n", best_indices[0]);
        //    printf("[CHECK] best_index = %d\n", best_indices[1]);
        //    printf("[CHECK] jet0: pt = %6.2f, eta = %6.2f, phi = %6.2f, e = %6.2f\n", jets[0].Pt(), jets[0].Eta(), jets[0].Phi(), jets[0].Energy());
        //    printf("[CHECK] jet1: pt = %6.2f, eta = %6.2f, phi = %6.2f, e = %6.2f\n", jets[1].Pt(), jets[1].Eta(), jets[1].Phi(), jets[1].Energy());
        //    printf("[CHECK] bjet: pt = %6.2f, eta = %6.2f, phi = %6.2f, e = %6.2f\n", bjet.Pt(), bjet.Eta(), bjet.Phi(), bjet.Energy());
        //    printf("[CHECK] w reco = %6.2f\n", w_candidate_reco.M());
        //    printf("[CHECK] t reco = %6.2f\n", top_candidate_reco.M());
        //}
        mean[0] += w_candidate_reco.M();
        mean[1] += top_candidate_reco.M();
        covarianceMatrix[0][0] += (w_candidate_reco.M() - w_boson_mass) * (w_candidate_reco.M() - w_boson_mass);
        covarianceMatrix[0][1] += (w_candidate_reco.M() - w_boson_mass) * (top_candidate_reco.M() - top_quark_mass);
        covarianceMatrix[1][0] += (top_candidate_reco.M() - top_quark_mass) * (w_candidate_reco.M() - w_boson_mass);
        covarianceMatrix[1][1] += (top_candidate_reco.M() - top_quark_mass) * (top_candidate_reco.M() - top_quark_mass);
        //==================================================//
        //-------------   Event Counting     ---------------//
        //==================================================//
        Nevents_pass_selection += 1;
        mytree.Fill();
    }// End of event loop.

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    MakePlots(c1, hist_gen_w_all, "no effect", "cov_hist_gen_w_all.png");
    MakePlots(c1, hist_gen_w_sel, "no effect", "cov_hist_gen_w_sel.png");


    mean[0] /= Nevents_pass_selection;
    mean[1] /= Nevents_pass_selection;
    covarianceMatrix[0][0] /= Nevents_pass_selection;
    covarianceMatrix[0][1] /= Nevents_pass_selection;
    covarianceMatrix[1][0] /= Nevents_pass_selection;
    covarianceMatrix[1][1] /= Nevents_pass_selection;

    printf("matrix(0,0) = %6.2f; matrix(0,1) = %6.2f;\n", covarianceMatrix[0][0], 0.);
    printf("matrix(1,0) = %6.2f; matrix(1,1) = %6.2f;\n", 0., covarianceMatrix[1][1]);
    printf("---------------\n");
    printf("matrix(0,0) = %6.2f; matrix(0,1) = %6.2f;\n", covarianceMatrix[0][0], covarianceMatrix[0][1]);
    printf("matrix(1,0) = %6.2f; matrix(1,1) = %6.2f;\n", covarianceMatrix[1][0], covarianceMatrix[1][1]);
    //printf("|%6.2f %6.2f|\n", covarianceMatrix[0][0], covarianceMatrix[0][1]);
    //printf("|%6.2f %6.2f|\n", covarianceMatrix[1][0], covarianceMatrix[1][1]);
    printf("mean = %6.2f, sigma_prime_jj  = %6.2f\n", mean[0], sqrt(covarianceMatrix[0][0]));
    printf("mean = %6.2f, sigma_prime_bjj = %6.2f\n", mean[1], sqrt(covarianceMatrix[1][1]));
    //==================================================//
    //---------------------  Report  -------------------//
    //==================================================//
    printf("[INFO] Nevents_pass_selection = %d\n", Nevents_pass_selection);
    fout->Write();
    fout->Close();
    return 1;
}


// functions{{{
bool checkAvailability(int index, std::vector<int> ID_IsChosen){
    bool result = true;//if pass the following for loop, the genParticle of the index is available.
    for(std::size_t i=0; i<ID_IsChosen.size(); ++i){
        if(index == ID_IsChosen[i]){ result = false; break; } 
    }
    return result;
}

void MakePlots(TCanvas *c1, TH1D* hist, const char* title, const char* outputFile){
    hist->Draw();
    //hist->SetTitle(title);
    //hist->SetXTitle(title);
    //hist->SetYTitle("Entries");
    hist->GetYaxis()->SetTitleOffset(1.4);
    hist->Write();
    c1->SaveAs(outputFile);
}
bool isThisMultiFile(char* dataset){
    if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return true;
    if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return true;
    return false;
}
//}}}
//class functions{{{
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
    ElecInfo_tmpPhoVeto = new std::vector<bool>;
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
    delete ElecInfo_tmpPhoVeto;
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
    flashggStdTree->SetBranchAddress("ElecInfo.tmpPhoVeto", &ElecInfo_tmpPhoVeto);
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
    //------------------------
    mytree -> Branch("GenMatching_Mass_top_gen", &GenMatching_Mass_top_gen, "GenMatching_Mass_top_gen/D");
    mytree -> Branch("GenMatching_Mass_wboson_gen", &GenMatching_Mass_wboson_gen, "GenMatching_Mass_wboson_gen/D");
    mytree -> Branch("GenMatching_Mass_top_reco", &GenMatching_Mass_top_reco, "GenMatching_Mass_top_reco/D");
    mytree -> Branch("GenMatching_Mass_wboson_reco", &GenMatching_Mass_wboson_reco, "GenMatching_Mass_wboson_reco/D");
    //------------------------
}
void myTreeClass::Fill(){
    mytree -> Fill();
}
void myParameters::Clear(){
    //------------------------
    GenMatching_Mass_top_gen = 0;
    GenMatching_Mass_wboson_gen = 0;
    GenMatching_Mass_top_reco = 0;
    GenMatching_Mass_wboson_reco = 0;
    //------------------------
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
//}}}
