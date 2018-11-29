#include <stdio.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <vector>
//ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X?fbclid=IwAR0QAekcVDaD2SL6lI7xFXkHxQigtkpuPiUHiP14t_i9pKvwfZ__v92MYiE

int main(){
    TFile *fin = TFile::Open("/wk_cms2/youying/public/forYuWei/tthTest.root");
    TTree *flashggStdTree = (TTree*)fin->Get("flashggNtuple/flashggStdTree");
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    TH1D  *hist_num_bjets = new TH1D("hist_num_bjets", "hist_num_bjets", 10, 0, 10);
    TH1D  *hist_wboson_mass_spectrum = new TH1D("hist_wboson_mass_spectrum", "hist_wboson_mass_spectrum", 50, 55, 105);
    TH1D  *hist_diphoton_mass = new TH1D("hist_diphoton_mass", "hist_diphoton_mass", 50, 100, 150);
    // Declaration of leaf types
    Int_t jets_size=0;
    std::vector<float> *JetInfo_Pt=0;
    std::vector<float> *JetInfo_Eta=0;
    std::vector<float> *JetInfo_Phi=0;
    std::vector<float> *JetInfo_Mass=0;
    std::vector<float> *JetInfo_Energy=0;
    std::vector<float> *JetInfo_pfDeepCSVJetTags_probb=0;
    std::vector<float> *JetInfo_pfDeepCSVJetTags_probbb=0;
    float DiPhoInfo_mass=0;
    // Set branch addresses
    flashggStdTree->SetBranchAddress("jets_size", &jets_size);
    flashggStdTree->SetBranchAddress("JetInfo.Pt", &JetInfo_Pt);
    flashggStdTree->SetBranchAddress("JetInfo.Eta", &JetInfo_Eta);
    flashggStdTree->SetBranchAddress("JetInfo.Phi", &JetInfo_Phi);
    flashggStdTree->SetBranchAddress("JetInfo.Mass", &JetInfo_Mass);
    flashggStdTree->SetBranchAddress("JetInfo.Energy", &JetInfo_Energy);
    flashggStdTree->SetBranchAddress("JetInfo.Energy", &JetInfo_Energy);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probb", &JetInfo_pfDeepCSVJetTags_probb);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probbb", &JetInfo_pfDeepCSVJetTags_probbb);
    flashggStdTree->SetBranchAddress("DiPhoInfo.mass", &DiPhoInfo_mass);

    double pfDeepCSVJetTags_tight  = 0.8001;
    double pfDeepCSVJetTags_medium = 0.4941;
    double pfDeepCSVJetTags_loose  = 0.1522;
    double w_boson_mass = 80.385;//GeV

    std::vector<TLorentzVector> vec_btagged_jets;
    std::vector<TLorentzVector> vec_nonbtagged_jets;

    //===== Loop Events =====//
    int nentries = flashggStdTree->GetEntries();
    printf("N = %d\n", nentries);
    for(int ientry=0; ientry<nentries; ientry++){
        flashggStdTree->GetEntry(ientry);//load data
        if((ientry+1)%100==0 || (ientry+1)==nentries) printf("ientry = %d\r", ientry);

        //##################################################//
        //##########     W-boson Mass Spectrum    ##########//
        //##################################################//
        if( !(jets_size>0) ) continue;
        //=== Select 2 b-tagged jets ===//
        int id_bjet1=0, id_bjet2=0, counter=0;
        vec_btagged_jets.clear();
        for(int i=0; i<jets_size; i++){
            TLorentzVector jet_lorentzvector;
            if(JetInfo_pfDeepCSVJetTags_probb->at(i)+JetInfo_pfDeepCSVJetTags_probbb->at(i) >= pfDeepCSVJetTags_tight){
                jet_lorentzvector.SetPtEtaPhiE(JetInfo_Pt->at(i), JetInfo_Eta->at(i), JetInfo_Phi->at(i), JetInfo_Energy->at(i));
                vec_btagged_jets.push_back(jet_lorentzvector);
                counter+=1;
            }
            if(counter==1) {id_bjet1=i;}
            if(counter==2) {id_bjet2=i; break;}
        }
        if(counter<2) continue;
        //=== Store rest jets ===//
        counter=0;
        vec_nonbtagged_jets.clear();
        for(int i=0; i<jets_size; i++){
            TLorentzVector jet_lorentzvector;
            if(i==id_bjet1 || i==id_bjet2) continue;
            jet_lorentzvector.SetPtEtaPhiE(JetInfo_Pt->at(i), JetInfo_Eta->at(i), JetInfo_Phi->at(i), JetInfo_Energy->at(i));
            vec_nonbtagged_jets.push_back(jet_lorentzvector);
            counter+=1;
        }
        if(counter<2) continue;
        hist_num_bjets->Fill(counter);
        //=== Chi-2 sorting ===//
        int id_jet1=0, id_jet2=0;
        double dijet_invariant_mass, chi2, chi2_min=999;
        for(int i=0; i<vec_nonbtagged_jets.size(); i++){
            for(int j=i+1; j<vec_nonbtagged_jets.size(); j++){
                TLorentzVector w_candidate_temp = vec_nonbtagged_jets[i] + vec_nonbtagged_jets[j];
                dijet_invariant_mass = w_candidate_temp.M();
                chi2 = (dijet_invariant_mass-w_boson_mass)*(dijet_invariant_mass-w_boson_mass);
                if(chi2<chi2_min){id_jet1=i; id_jet2=j; chi2_min=chi2;};
            }
        }
        TLorentzVector best_dijet_pair = vec_nonbtagged_jets[id_jet1] + vec_nonbtagged_jets[id_jet2];
        dijet_invariant_mass = best_dijet_pair.M();
        hist_wboson_mass_spectrum->Fill(dijet_invariant_mass);

        //##################################################//
        //##########    DiPhoton Mass Spectrum    ##########//
        //##################################################//
        if(DiPhoInfo_mass>0) hist_diphoton_mass->Fill(DiPhoInfo_mass);
    }// End of event loop.


    //##################################################//
    //##############     Make Plots !!    ##############//
    //##################################################//
    //hist_num_bjets->Draw();
    hist_wboson_mass_spectrum->Draw();
    hist_wboson_mass_spectrum->SetTitle("W boson mass spectrum");
    hist_wboson_mass_spectrum->SetXTitle("dijet_invariant_mass [GeV / c^{2}]");
    hist_wboson_mass_spectrum->SetYTitle("Entries / 1 GeV");
    hist_wboson_mass_spectrum->GetYaxis()->SetTitleOffset(1.4);
    c1->SaveAs("hist_wboson_mass_spectrum.png");

    hist_diphoton_mass->Draw();
    hist_diphoton_mass->SetTitle("Diphoton mass spectrum");
    hist_diphoton_mass->SetXTitle("DiPhoInfo_mass [GeV / c^{2}]");
    hist_diphoton_mass->SetYTitle("Entries / 1 GeV");
    hist_diphoton_mass->GetYaxis()->SetTitleOffset(1.4);
    c1->SaveAs("hist_diphoton_mass.png");
    return 1;
}
