#include <stdio.h>
#include <math.h>
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

    //==================//
    //--- histograms ---//
    //==================//
    TH1D  *hist_num_nonbtagged_jets = new TH1D("hist_num_nonbtagged_jets", "hist_num_nonbtagged_jets", 20, 0, 20);
    TH1D  *hist_wboson_mass_spectrum = new TH1D("hist_wboson_mass_spectrum", "hist_wboson_mass_spectrum", 50, 55, 105);
    TH1D  *hist_cjet_pt = new TH1D("hist_cjet_pt", "hist_cjet_pt", 50, 0, 1000);
    TH1D  *hist_inv_mass_tch = new TH1D("hist_inv_mass_tch", "hist_inv_mass_tch", 50, 0, 500);
    TH1D  *hist_inv_mass_tbw = new TH1D("hist_inv_mass_tbw", "hist_inv_mass_tbw", 50, 0, 500);
    TH1D  *hist = new TH1D("hist", "hist", 50, 0, 1000);
    //TH1D  *hist = new TH1D("hist", "hist", 50, -10, 10);

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
    float DiPhoInfo_mass=0;
    float DiPhoInfo_leadPt=0;
    float DiPhoInfo_leadEta=0;
    float DiPhoInfo_leadPhi=0;
    float DiPhoInfo_leadE=0;
    //------------------------
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
    //################    Event Loop     ###############//
    //##################################################//
    int nentries = flashggStdTree->GetEntries(); printf("N = %d\n", nentries);
    for(int ientry=0; ientry<nentries; ientry++){
        flashggStdTree->GetEntry(ientry);//load data
        if((ientry+1)%100==0 || (ientry+1)==nentries) printf("ientry = %d\r", ientry);
        //==================================================//
        //-------------   Reset Parameters   ---------------//
        //==================================================//
        vec_nonbtagged_jets.clear();
        //==================================================//
        //--------------   Basic Selectoin   ---------------//
        //==================================================//
        if(DiPhoInfo_mass<0) continue;
        if( !(jets_size>0) ) continue;
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
            }
            num_bjet+=1;
            if(num_bjet==2) break;
        }
        if(bjetindex==-1) continue;
        //if(num_bjet!=1)   continue;
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
        hist_num_nonbtagged_jets->Fill(num_nonbtagged_jets);
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
        hist_wboson_mass_spectrum->Fill(dijet_invariant_mass);
        //=== Determine cjet ===//
        int cjetindex=-1;//indices in non-btagged vector
        TLorentzVector leading_cjet;
        for(int i=0; i<vec_nonbtagged_jets.size(); i++){
            if(i==wjetindices[0] || i==wjetindices[1]) continue;
            cjetindex=i;
            leading_cjet.SetPtEtaPhiE(JetInfo_Pt->at(i), JetInfo_Eta->at(i), JetInfo_Phi->at(i), JetInfo_Energy->at(i));
            if(cjetindex!=-1) break;
        }
        hist_cjet_pt->Fill(leading_cjet.Pt());
        //=== InvMass of tops ===//
        TLorentzVector tch_jgg = leading_diphoton + leading_cjet;
        TLorentzVector tbw_bjj = best_dijet_pair + leading_bjet;
        double inv_mass_tch = tch_jgg.M();
        double inv_mass_tbw = tbw_bjj.M();
        hist_inv_mass_tch->Fill(inv_mass_tch);
        hist_inv_mass_tbw->Fill(inv_mass_tbw);

    }// End of event loop.


    //##################################################//
    //##############     Make Plots !!    ##############//
    //##################################################//
    printf("Start making plots!\n");
    //hist->Draw();
    //c1->SaveAs("plots/hist.png");
    hist_num_nonbtagged_jets->Draw();
    c1->SaveAs("plots/hist_num_nonbtagged_jets.png");
    hist_cjet_pt->Draw();
    c1->SaveAs("plots/hist_cjet_pt.png");
    hist_wboson_mass_spectrum->Draw();
    hist_wboson_mass_spectrum->SetTitle("W boson mass spectrum");
    hist_wboson_mass_spectrum->SetXTitle("dijet_invariant_mass [GeV / c^{2}]");
    hist_wboson_mass_spectrum->SetYTitle("Entries / 1 GeV");
    hist_wboson_mass_spectrum->GetYaxis()->SetTitleOffset(1.4);
    c1->SaveAs("plots/hist_wboson_mass_spectrum.png");
    hist_inv_mass_tch->Draw();
    c1->SaveAs("plots/hist_inv_mass_tch.png");
    hist_inv_mass_tbw->Draw();
    c1->SaveAs("plots/hist_inv_mass_tbw.png");

    return 1;
}
