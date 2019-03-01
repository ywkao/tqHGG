//=== 10. Jan. 2019 ===//
//=== ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X?fbclid=IwAR0QAekcVDaD2SL6lI7xFXkHxQigtkpuPiUHiP14t_i9pKvwfZ__v92MYiE
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
using namespace std;

int main(int argc, char *argv[]){
    //const char *input_file = "/wk_cms2/youying/public/2017_94X_3_1_X_and_3_2_0/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root";
    char input_file[256]; sprintf(input_file, "%s", argv[1]); printf("[INFO] input_file  = %s\n", input_file);
    char output_file[256]; sprintf(output_file, "%s", argv[2]); printf("[INFO] output_file = %s\n", output_file);
    char dataset[256]; sprintf(dataset, "%s", argv[3]); printf("[INFO] dataset     = %s\n", dataset);
    bool isData = isThisDataOrNot(dataset);
    bool isMultiFile = isThisMultiFile(dataset);
    TFile *fout = new TFile(output_file, "RECREATE");
    TChain *flashggStdTree = new TChain("flashggNtuples/flashggStdTree");
    //if(isMultiFile) flashggStdTree->Add("/home/ykao/legacy/CMSSW_9_4_10/src/t2cH/test/*.root"); // For test only.
    if(isMultiFile) flashggStdTree->Add(Form("%s/*.root", input_file));
    else flashggStdTree->Add(input_file);
    //=== input old ver. ===//
    //TFile *fin  = TFile::Open(input_file);
    //TTree *flashggStdTree = (TTree*)fin->Get("flashggNtuples/flashggStdTree");
    //=== test ===//
    //TFile *fin  = TFile::Open("/wk_cms2/youying/public/forYuWei/tthTest.root");
    //TFile *fout = new TFile("plots/hist_tth.root", "RECREATE");
    //TTree *flashggStdTree = (TTree*)fin->Get("flashggNtuple/flashggStdTree");
    //TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

    //==================//
    //--- histograms ---//
    //==================//
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    TH1D  *hist_num_jets = new TH1D("hist_num_jets", "hist_num_jets", 20, 0, 20);
    TH1D  *hist_num_btagged_jets = new TH1D("hist_num_btagged_jets", "hist_num_btagged_jets", 10, 0, 10);
    TH1D  *hist_num_nonbtagged_jets = new TH1D("hist_num_nonbtagged_jets", "hist_num_nonbtagged_jets", 20, 0, 20);
    TH1D  *hist_wboson_mass_spectrum = new TH1D("hist_wboson_mass_spectrum", "hist_wboson_mass_spectrum", 50, 55, 105);
    TH1D  *hist_diphoton_mass_spectrum = new TH1D("hist_diphoton_mass_spectrum", "hist_diphoton_mass_spectrum", 50, 100, 150);
    TH1D  *hist_bjet_pt = new TH1D("hist_bjet_pt", "hist_bjet_pt", 50, 0, 1000);
    TH1D  *hist_jet1_pt = new TH1D("hist_jet1_pt", "hist_jet1_pt", 50, 0, 1000);
    TH1D  *hist_jet2_pt = new TH1D("hist_jet2_pt", "hist_jet2_pt", 50, 0, 1000);
    TH1D  *hist_cjet_pt = new TH1D("hist_cjet_pt", "hist_cjet_pt", 50, 0, 1000);
    TH1D  *hist_bjet_eta = new TH1D("hist_bjet_eta", "hist_bjet_eta", 40, -2.5, 2.5);
    TH1D  *hist_jet1_eta = new TH1D("hist_jet1_eta", "hist_jet1_eta", 40, -2.5, 2.5);
    TH1D  *hist_jet2_eta = new TH1D("hist_jet2_eta", "hist_jet2_eta", 40, -2.5, 2.5);
    TH1D  *hist_cjet_eta = new TH1D("hist_cjet_eta", "hist_cjet_eta", 40, -2.5, 2.5);
    TH1D  *hist_bjet_phi = new TH1D("hist_bjet_phi", "hist_bjet_phi", 40, -3.0, 3.0);
    TH1D  *hist_jet1_phi = new TH1D("hist_jet1_phi", "hist_jet1_phi", 40, -3.0, 3.0);
    TH1D  *hist_jet2_phi = new TH1D("hist_jet2_phi", "hist_jet2_phi", 40, -3.0, 3.0);
    TH1D  *hist_cjet_phi = new TH1D("hist_cjet_phi", "hist_cjet_phi", 40, -3.0, 3.0);
    TH1D  *hist_inv_mass_tch = new TH1D("hist_inv_mass_tch", "hist_inv_mass_tch", 50, 0, 500);
    TH1D  *hist_inv_mass_tbw = new TH1D("hist_inv_mass_tbw", "hist_inv_mass_tbw", 50, 0, 500);
    //TH1D  *hist = new TH1D("hist", "hist", 50, 0, 1000);
    //TH1D  *hist = new TH1D("hist", "hist", 50, -10, 10);
    //------------------------
    hist_num_jets -> Sumw2();
    hist_num_btagged_jets -> Sumw2();
    hist_num_nonbtagged_jets -> Sumw2();
    hist_wboson_mass_spectrum -> Sumw2();
    hist_diphoton_mass_spectrum -> Sumw2();
    hist_bjet_pt -> Sumw2();
    hist_jet1_pt -> Sumw2();
    hist_jet2_pt -> Sumw2();
    hist_cjet_pt -> Sumw2();
    hist_bjet_eta -> Sumw2();
    hist_jet1_eta -> Sumw2();
    hist_jet2_eta -> Sumw2();
    hist_cjet_eta -> Sumw2();
    hist_bjet_phi -> Sumw2();
    hist_jet1_phi -> Sumw2();
    hist_jet2_phi -> Sumw2();
    hist_cjet_phi -> Sumw2();
    hist_inv_mass_tch -> Sumw2();
    hist_inv_mass_tbw -> Sumw2();

    //==========================//
    //--- readout parameters ---//
    //==========================//
    Int_t jets_size=0;
    std::vector<float> *JetInfo_Pt=0;
    std::vector<float> *JetInfo_Eta=0;
    std::vector<float> *JetInfo_Phi=0;
    std::vector<float> *JetInfo_Mass=0;
    std::vector<float> *JetInfo_Energy=0;
    std::vector<float> *JetInfo_pfDeepCSVJetTags_probb=0;
    std::vector<float> *JetInfo_pfDeepCSVJetTags_probbb=0;
    float EvtInfo_genweight=0;
    float DiPhoInfo_mass=0;
    float DiPhoInfo_leadPt=0;
    float DiPhoInfo_leadEta=0;
    float DiPhoInfo_leadPhi=0;
    float DiPhoInfo_leadE=0;
    //------------------------
    flashggStdTree->SetBranchAddress("EvtInfo.genweight", &EvtInfo_genweight);
    flashggStdTree->SetBranchAddress("jets_size", &jets_size);
    flashggStdTree->SetBranchAddress("JetInfo.Pt", &JetInfo_Pt);
    flashggStdTree->SetBranchAddress("JetInfo.Eta", &JetInfo_Eta);
    flashggStdTree->SetBranchAddress("JetInfo.Phi", &JetInfo_Phi);
    flashggStdTree->SetBranchAddress("JetInfo.Mass", &JetInfo_Mass);
    flashggStdTree->SetBranchAddress("JetInfo.Energy", &JetInfo_Energy);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probb", &JetInfo_pfDeepCSVJetTags_probb);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probbb", &JetInfo_pfDeepCSVJetTags_probbb);
    flashggStdTree->SetBranchAddress("DiPhoInfo.mass", &DiPhoInfo_mass);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadPt", &DiPhoInfo_leadPt);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadEta", &DiPhoInfo_leadEta);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadPhi", &DiPhoInfo_leadPhi);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadE", &DiPhoInfo_leadE);

    //==========================//
    //---  Useful Constants  ---//
    //==========================//
    double pfDeepCSVJetTags_tight  = 0.8001;
    double pfDeepCSVJetTags_medium = 0.4941;
    double pfDeepCSVJetTags_loose  = 0.1522;
    double w_boson_mass = 80.385;//GeV

    std::vector<TLorentzVector> vec_nonbtagged_jets;

    //##################################################//
    //########   Event Loop [Normalization]   ##########//
    //##################################################//
    int nentries = flashggStdTree->GetEntries(); printf("[INFO] N_entries = %d\n", nentries);
    double NormalizationFactor;
    double Luminosity = 42.; //fb{-1}
    double CrossSection = GetXsec(dataset); //pb
    double BranchingFraction = GetBranchingFraction(dataset); //pb
    printf("[INFO] CrossSection = %f !\n", CrossSection);
    printf("[INFO] Equivalent lumi. = %f !\n", (double)nentries/CrossSection);
    printf("[INFO] BranchingFraction = %f !\n", BranchingFraction);
    double TotalGenweight=0;
    for(int ientry=0; ientry<nentries; ientry++){
        flashggStdTree->GetEntry(ientry);
        TotalGenweight+=EvtInfo_genweight;
    }
    NormalizationFactor = 1000. * Luminosity * CrossSection * BranchingFraction / TotalGenweight;
    printf("\n[INFO] TotalGenweight = %f!\n", TotalGenweight);
    printf("[INFO] NormalizationFactor = %f!\n", isData ? 1. : NormalizationFactor);
    //##################################################//
    //#########    Event Loop [Selection]    ###########//
    //##################################################//
    int nevents_pass_selection = 0;
    for(int ientry=0; ientry<nentries; ientry++){
        flashggStdTree->GetEntry(ientry);//load data
        if((ientry+1)%10000==0 || (ientry+1)==nentries) printf("ientry = %d\r", ientry);
        //==================================================//
        //-------------   Reset Parameters   ---------------//
        //==================================================//
        vec_nonbtagged_jets.clear();
        //==================================================//
        //--------------   Basic Selectoin   ---------------//
        //==================================================//
        if(DiPhoInfo_mass<0) continue;
        if( !(jets_size>0) ) continue;
        if(DiPhoInfo_mass<100 || DiPhoInfo_mass>150) continue;
        if(DiPhoInfo_mass>120 && DiPhoInfo_mass<130) continue;
        //==================================================//
        //------------   Physical Observables   ------------//
        //==================================================//
        //==================================================//
        //---------   Diphoton Jet Mass Spectrum   ---------//
        //==================================================//
        //=== t-cH, H-gg; t-bW, W-jj ===//
        TLorentzVector leading_diphoton;
        leading_diphoton.SetPtEtaPhiE(DiPhoInfo_leadPt, DiPhoInfo_leadEta, DiPhoInfo_leadPhi, DiPhoInfo_leadE);
        //#######################################################//
        //### 1 bjet + 2 chi2 jets + the leading jet of rest ones
        //#######################################################//
        //=== Select one bjet ===//
        int bjetindex=-1, num_bjet=0;
        TLorentzVector leading_bjet;
        for(int i=0; i<jets_size; i++){
            if( fabs(JetInfo_Eta->at(i)) > 2.4 ) continue;
            if( fabs(JetInfo_Pt->at(i))  < 30  ) continue;
            if(JetInfo_pfDeepCSVJetTags_probb->at(i)+JetInfo_pfDeepCSVJetTags_probbb->at(i) >= pfDeepCSVJetTags_tight){
                bjetindex = i;
                leading_bjet.SetPtEtaPhiE(JetInfo_Pt->at(i), JetInfo_Eta->at(i), JetInfo_Phi->at(i), JetInfo_Energy->at(i));
                num_bjet+=1;
            }
            //if(num_bjet==2) break;
        }
        hist_num_btagged_jets->Fill(num_bjet, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        if(bjetindex==-1) continue;
        if(num_bjet!=1)   continue;
        hist_bjet_pt->Fill(leading_bjet.Pt(), isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_bjet_eta->Fill(leading_bjet.Eta(), isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_bjet_phi->Fill(leading_bjet.Phi(), isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        //=== Store rest jets ===//
        int num_nonbtagged_jets=0;
        for(int i=0; i<jets_size; i++){
            if( fabs(JetInfo_Eta->at(i)) > 2.4 ) continue;
            if( fabs(JetInfo_Pt->at(i))  < 30  ) continue;
            if(i == bjetindex) continue;
            TLorentzVector jet_lorentzvector;
            jet_lorentzvector.SetPtEtaPhiE(JetInfo_Pt->at(i), JetInfo_Eta->at(i), JetInfo_Phi->at(i), JetInfo_Energy->at(i));
            vec_nonbtagged_jets.push_back(jet_lorentzvector);
            num_nonbtagged_jets+=1;
        }
        if(num_nonbtagged_jets<3) continue;
        hist_num_nonbtagged_jets->Fill(num_nonbtagged_jets, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_num_jets->Fill(num_bjet+num_nonbtagged_jets, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        //=== Chi-2 sorting ===//
        int wjetindices[2]={0};//indices in non-btagged vector
        double dijet_invariant_mass, chi2, chi2_min=9999;
        TLorentzVector best_dijet_pair;
        for(int i=0; i<vec_nonbtagged_jets.size(); i++){
            for(int j=i+1; j<vec_nonbtagged_jets.size(); j++){
                TLorentzVector w_candidate_temp = vec_nonbtagged_jets[i] + vec_nonbtagged_jets[j];
                dijet_invariant_mass = w_candidate_temp.M();
                chi2 = (dijet_invariant_mass-w_boson_mass)*(dijet_invariant_mass-w_boson_mass);
                if(chi2<chi2_min){wjetindices[0]=i; wjetindices[1]=j; chi2_min=chi2;};
            }
        }
        best_dijet_pair = vec_nonbtagged_jets[wjetindices[0]] + vec_nonbtagged_jets[wjetindices[1]];
        dijet_invariant_mass = best_dijet_pair.M();
        hist_wboson_mass_spectrum->Fill(dijet_invariant_mass, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_jet1_pt->Fill(vec_nonbtagged_jets[wjetindices[0]].Pt(), isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_jet1_eta->Fill(vec_nonbtagged_jets[wjetindices[0]].Eta(), isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_jet1_phi->Fill(vec_nonbtagged_jets[wjetindices[0]].Phi(), isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_jet2_pt->Fill(vec_nonbtagged_jets[wjetindices[1]].Pt(), isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_jet2_eta->Fill(vec_nonbtagged_jets[wjetindices[1]].Eta(), isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_jet2_phi->Fill(vec_nonbtagged_jets[wjetindices[1]].Phi(), isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        //=== Diphoton inv mass ===//
        hist_diphoton_mass_spectrum->Fill(DiPhoInfo_mass, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        //=== Determine cjet ===//
        int cjetindex=-1;//indices in non-btagged vector
        TLorentzVector leading_cjet;
        for(int i=0; i<vec_nonbtagged_jets.size(); i++){
            if(i==wjetindices[0] || i==wjetindices[1]) continue;
            cjetindex=i;
            leading_cjet.SetPtEtaPhiE(JetInfo_Pt->at(i), JetInfo_Eta->at(i), JetInfo_Phi->at(i), JetInfo_Energy->at(i));
            if(cjetindex!=-1) break;
        }
        hist_cjet_pt->Fill(leading_cjet.Pt(), isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_cjet_eta->Fill(leading_cjet.Eta(), isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_cjet_phi->Fill(leading_cjet.Phi(), isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        //=== InvMass of tops ===//
        TLorentzVector tch_jgg = leading_diphoton + leading_cjet;
        TLorentzVector tbw_bjj = best_dijet_pair + leading_bjet;
        double inv_mass_tch = tch_jgg.M();
        double inv_mass_tbw = tbw_bjj.M();
        hist_inv_mass_tch->Fill(inv_mass_tch, isData ? 1. : NormalizationFactor*EvtInfo_genweight);
        hist_inv_mass_tbw->Fill(inv_mass_tbw, isData ? 1. : NormalizationFactor*EvtInfo_genweight);

        //==================================================//
        //-------------   Event Selection    ---------------//
        //==================================================//
        nevents_pass_selection += 1;
    }// End of event loop.
    printf("[INFO] N_nevents_pass_selection = %d\n", nevents_pass_selection);


    //##################################################//
    //##############     Make Plots !!    ##############//
    //##################################################//
    printf("[INFO] Start making plots!\n");
    char output_dir[256]; sprintf(output_dir, "%s", argv[4]); printf("[INFO] output_dir = %s\n", output_dir);
    //------------------------------
    hist_num_jets->Draw();
    hist_num_jets->SetTitle("Multiplicity of jets");
    hist_num_jets->SetXTitle("# of jets");
    hist_num_jets->SetYTitle("Entries");
    hist_num_jets->GetYaxis()->SetTitleOffset(1.4);
    hist_num_jets->Write();
    c1->SaveAs(Form("%s/hist_num_jets.png", output_dir));
    //------------------------------
    hist_num_btagged_jets->Draw();
    hist_num_btagged_jets->SetTitle("Multiplicity of b-tagged jets");
    hist_num_btagged_jets->SetXTitle("# of b-tagged jets");
    hist_num_btagged_jets->SetYTitle("Entries");
    hist_num_btagged_jets->GetYaxis()->SetTitleOffset(1.4);
    hist_num_btagged_jets->Write();
    c1->SaveAs(Form("%s/hist_num_btagged_jets.png", output_dir));
    //------------------------------
    hist_num_nonbtagged_jets->Draw();
    hist_num_nonbtagged_jets->SetTitle("Multiplicity of non-b-tagged jets");
    hist_num_nonbtagged_jets->SetXTitle("# of non-b-tagged jets");
    hist_num_nonbtagged_jets->SetYTitle("Entries");
    hist_num_nonbtagged_jets->GetYaxis()->SetTitleOffset(1.4);
    hist_num_nonbtagged_jets->Write();
    c1->SaveAs(Form("%s/hist_num_nonbtagged_jets.png", output_dir));
    //------------------------------
    hist_bjet_pt->Draw();
    hist_bjet_pt->SetTitle("Pt of b-tagged jet");
    hist_bjet_pt->SetXTitle("Pt of b-tagged jet [GeV]");
    hist_bjet_pt->SetYTitle("Entries / 20 [GeV]");
    hist_bjet_pt->GetYaxis()->SetTitleOffset(1.4);
    hist_bjet_pt->Write();
    c1->SaveAs(Form("%s/hist_bjet_pt.png", output_dir));
    //------------------------------
    hist_bjet_eta->Draw();
    hist_bjet_eta->SetTitle("Eta of b-tagged jet");
    hist_bjet_eta->SetXTitle("Eta of b-tagged jet");
    hist_bjet_eta->SetYTitle("Entries / 0.5 unit");
    hist_bjet_eta->GetYaxis()->SetTitleOffset(1.4);
    hist_bjet_eta->Write();
    c1->SaveAs(Form("%s/hist_bjet_eta.png", output_dir));
    //------------------------------
    hist_bjet_phi->Draw();
    hist_bjet_phi->SetTitle("Phi of b-tagged jet");
    hist_bjet_phi->SetXTitle("Phi of b-tagged jet");
    hist_bjet_phi->SetYTitle("Entries / 0.5 [rad]");
    hist_bjet_phi->GetYaxis()->SetTitleOffset(1.4);
    hist_bjet_phi->Write();
    c1->SaveAs(Form("%s/hist_bjet_phi.png", output_dir));
    //------------------------------
    hist_jet1_pt->Draw();
    hist_jet1_pt->SetTitle("Pt of hadronic jet1 candidate");
    hist_jet1_pt->SetXTitle("Pt of hadronic jet1 candidate [GeV]");
    hist_jet1_pt->SetYTitle("Entries / 20 [GeV]");
    hist_jet1_pt->GetYaxis()->SetTitleOffset(1.4);
    hist_jet1_pt->Write();
    c1->SaveAs(Form("%s/hist_jet1_pt.png", output_dir));
    //------------------------------
    hist_jet1_eta->Draw();
    hist_jet1_eta->SetTitle("Eta of hadronic jet1 candidate");
    hist_jet1_eta->SetXTitle("Eta of hadronic jet1 candidate");
    hist_jet1_eta->SetYTitle("Entries / 0.5 unit");
    hist_jet1_eta->GetYaxis()->SetTitleOffset(1.4);
    hist_jet1_eta->Write();
    c1->SaveAs(Form("%s/hist_jet1_eta.png", output_dir));
    //------------------------------
    hist_jet1_phi->Draw();
    hist_jet1_phi->SetTitle("Phi of hadronic jet1 candidate");
    hist_jet1_phi->SetXTitle("Phi of hadronic jet1 candidate");
    hist_jet1_phi->SetYTitle("Entries / 0.5 [rad]");
    hist_jet1_phi->GetYaxis()->SetTitleOffset(1.4);
    hist_jet1_phi->Write();
    c1->SaveAs(Form("%s/hist_jet1_phi.png", output_dir));
    //------------------------------
    hist_jet2_pt->Draw();
    hist_jet2_pt->SetTitle("Pt of hadronic jet2 candidate");
    hist_jet2_pt->SetXTitle("Pt of hadronic jet2 candidate [GeV]");
    hist_jet2_pt->SetYTitle("Entries / 20 [GeV]");
    hist_jet2_pt->GetYaxis()->SetTitleOffset(1.4);
    hist_jet2_pt->Write();
    c1->SaveAs(Form("%s/hist_jet2_pt.png", output_dir));
    //------------------------------
    hist_jet2_eta->Draw();
    hist_jet2_eta->SetTitle("Eta of hadronic jet2 candidate");
    hist_jet2_eta->SetXTitle("Eta of hadronic jet2 candidate");
    hist_jet2_eta->SetYTitle("Entries / 0.5 unit");
    hist_jet2_eta->GetYaxis()->SetTitleOffset(1.4);
    hist_jet2_eta->Write();
    c1->SaveAs(Form("%s/hist_jet2_eta.png", output_dir));
    //------------------------------
    hist_jet2_phi->Draw();
    hist_jet2_phi->SetTitle("Phi of hadronic jet2 candidate");
    hist_jet2_phi->SetXTitle("Phi of hadronic jet2 candidate");
    hist_jet2_phi->SetYTitle("Entries / 0.5 [rad]");
    hist_jet2_phi->GetYaxis()->SetTitleOffset(1.4);
    hist_jet2_phi->Write();
    c1->SaveAs(Form("%s/hist_jet2_phi.png", output_dir));
    //------------------------------
    hist_cjet_pt->Draw();
    hist_cjet_pt->SetTitle("Pt of cjet candidates");
    hist_cjet_pt->SetXTitle("Pt of cjet candidates [GeV]");
    hist_cjet_pt->SetYTitle("Entries / 20 [GeV]");
    hist_cjet_pt->GetYaxis()->SetTitleOffset(1.4);
    hist_cjet_pt->Write();
    c1->SaveAs(Form("%s/hist_cjet_pt.png", output_dir));
    //------------------------------
    hist_cjet_eta->Draw();
    hist_cjet_eta->SetTitle("Eta of cjet candidates");
    hist_cjet_eta->SetXTitle("Eta of cjet candidates");
    hist_cjet_eta->SetYTitle("Entries / 0.5 unit");
    hist_cjet_eta->GetYaxis()->SetTitleOffset(1.4);
    hist_cjet_eta->Write();
    c1->SaveAs(Form("%s/hist_cjet_eta.png", output_dir));
    //------------------------------
    hist_cjet_phi->Draw();
    hist_cjet_phi->SetTitle("Phi of cjet candidates");
    hist_cjet_phi->SetXTitle("Phi of cjet candidates");
    hist_cjet_phi->SetYTitle("Entries / 0.5 [rad]");
    hist_cjet_phi->GetYaxis()->SetTitleOffset(1.4);
    hist_cjet_phi->Write();
    c1->SaveAs(Form("%s/hist_cjet_phi.png", output_dir));
    //------------------------------
    hist_wboson_mass_spectrum->Draw();
    hist_wboson_mass_spectrum->SetTitle("W boson mass spectrum");
    hist_wboson_mass_spectrum->SetXTitle("dijet invariant mass [GeV / c^{2}]");
    hist_wboson_mass_spectrum->SetYTitle("Entries / 1 GeV");
    hist_wboson_mass_spectrum->GetYaxis()->SetTitleOffset(1.4);
    hist_wboson_mass_spectrum->Write();
    c1->SaveAs(Form("%s/hist_wboson_mass_spectrum.png", output_dir));
    //------------------------------
    hist_diphoton_mass_spectrum->Draw();
    hist_diphoton_mass_spectrum->SetTitle("Diphoton mass spectrum");
    hist_diphoton_mass_spectrum->SetXTitle("Diphoton invariant mass [GeV / c^{2}]");
    hist_diphoton_mass_spectrum->SetYTitle("Entries / 1 GeV");
    hist_diphoton_mass_spectrum->GetYaxis()->SetTitleOffset(1.4);
    hist_diphoton_mass_spectrum->Write();
    c1->SaveAs(Form("%s/hist_diphoton_mass_spectrum.png", output_dir));
    //------------------------------
    hist_inv_mass_tch->Draw();
    hist_inv_mass_tch->SetTitle("M1 mass spectrum");
    hist_inv_mass_tch->SetXTitle("diphoton + cjet candidate invariant_mass [GeV/c^{2}]");
    hist_inv_mass_tch->SetYTitle("Entries / 10 GeV");
    hist_inv_mass_tch->GetYaxis()->SetTitleOffset(1.4);
    hist_inv_mass_tch->Write();
    c1->SaveAs(Form("%s/hist_inv_mass_tch.png", output_dir));
    //------------------------------
    hist_inv_mass_tbw->Draw();
    hist_inv_mass_tbw->SetTitle("M2 mass spectrum");
    hist_inv_mass_tbw->SetXTitle("dijet + b-tagged jet invariant_mass [GeV / c^{2}]");
    hist_inv_mass_tbw->SetYTitle("Entries / 10 GeV");
    hist_inv_mass_tbw->GetYaxis()->SetTitleOffset(1.4);
    hist_inv_mass_tbw->Write();
    c1->SaveAs(Form("%s/hist_inv_mass_tbw.png", output_dir));

    fout->Close();
    return 1;
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
    return false;
}
