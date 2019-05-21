#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>

const int N_dataset = 27;
double SumContents(TH1D* &h);
void Normalization(TH1D* &h, double Ntotal);

void puComparison2017(){
    TH1D *h_winter = new TH1D("h_winter", "", 100, 0, 100);
    TH1D *h_UltraLegacy = new TH1D("h_UltraLegacy", ";Number of pileup; probValue", 100, 0, 100);

    TFile *fin = TFile::Open("plots/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/hist_GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8.root");
    TFile *fin_pu = TFile::Open("plots_pu/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/hist_GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8.root");

    TH1D *h_mc_Npu = (TH1D*)fin->Get("hist_NPu");
    TH1D *h_mc_Npu_afterPUreweighting = (TH1D*)fin_pu->Get("hist_NPu");

    double Ntotal[2]={0};
    Ntotal[0] = SumContents(h_mc_Npu);
    Ntotal[1] = SumContents(h_mc_Npu_afterPUreweighting);
    printf("Ntotal = %f\n", Ntotal[0]);
    Normalization(h_mc_Npu, Ntotal[0]);
    Normalization(h_mc_Npu_afterPUreweighting, Ntotal[1]);

    double Nbins[99] = {
               0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
							   15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 
							   27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 
							   39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 
							   51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 
							   63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 
							   75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 
							   87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98
    };

    //https://raw.githubusercontent.com/cms-sw/cmssw/master/SimGeneral/MixingModule/python/mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi.py
    double pu_winter[99] = {
              3.39597497605e-05,
              6.63688402133e-06,
              1.39533611284e-05,
              3.64963078209e-05,
              6.00872171664e-05,
              9.33932578027e-05,
              0.000120591524486,
              0.000128694546198,
              0.000361697233219,
              0.000361796847553,
              0.000702474896113,
              0.00133766053707,
              0.00237817050805,
              0.00389825605651,
              0.00594546732588,
              0.00856825906255,
              0.0116627396044,
              0.0148793350787,
              0.0179897368379,
              0.0208723871946,
              0.0232564170641,
              0.0249826433945,
              0.0262245860346,
              0.0272704617569,
              0.0283301107549,
              0.0294006137386,
              0.0303026836965,
              0.0309692426278,
              0.0308818046328,
              0.0310566806228,
              0.0309692426278,
              0.0310566806228,
              0.0310566806228,
              0.0310566806228,
              0.0307696426944,
              0.0300103336052,
              0.0288355370103,
              0.0273233309106,
              0.0264343533951,
              0.0255453758796,
              0.0235877272306,
              0.0215627588047,
              0.0195825559393,
              0.0177296309658,
              0.0160560731931,
              0.0146022004183,
              0.0134080690078,
              0.0129586991411,
              0.0125093292745,
              0.0124360740539,
              0.0123547104433,
              0.0123953922486,
              0.0124360740539,
              0.0124360740539,
              0.0123547104433,
              0.0124360740539,
              0.0123387597772,
              0.0122414455005,
              0.011705203844,
              0.0108187105305,
              0.00963985508986,
              0.00827210065136,
              0.00683770076341,
              0.00545237697118,
              0.00420456901556,
              0.00367513566191,
              0.00314570230825,
              0.0022917978982,
              0.00163221454973,
              0.00114065309494,
              0.000784838366118,
              0.000533204105387,
              0.000358474034915,
              0.000238881117601,
              0.0001984254989,
              0.000157969880198,
              0.00010375646169,
              6.77366175538e-05,
              4.39850477645e-05,
              2.84298066026e-05,
              1.83041729561e-05,
              1.17473542058e-05,
              7.51982735129e-06,
              6.16160108867e-06,
              4.80337482605e-06,
              3.06235473369e-06,
              1.94863396999e-06,
              1.23726800704e-06,
              7.83538083774e-07,
              4.94602064224e-07,
              3.10989480331e-07,
              1.94628487765e-07,
              1.57888581037e-07,
              1.2114867431e-07,
              7.49518929908e-08,
              4.6060444984e-08,
              2.81008884326e-08,
              1.70121486128e-08,
              1.02159894812e-08
    };

    //https://raw.githubusercontent.com/cms-sw/cmssw/master/SimGeneral/MixingModule/python/mix_2017_25ns_UltraLegacy_PoissonOOTPU_cfi.py
    double pu_ultralegacy[99] = {
        1.1840841518e-05, 3.46661037703e-05, 8.98772521472e-05, 7.47400487733e-05, 0.000123005176624,
        0.000156501700614, 0.000154660478659, 0.000177496185603, 0.000324149805611, 0.000737524009713,
        0.00140432980253, 0.00244424508696, 0.00380027898037, 0.00541093042612, 0.00768803501793,
        0.010828224552, 0.0146608623707, 0.01887739113, 0.0228418813823, 0.0264817796874,
        0.0294637401336, 0.0317960986171, 0.0336645950831, 0.0352638818387, 0.036869429333,
        0.0382797316998, 0.039386705577, 0.0398389681346, 0.039646211131, 0.0388392805703,
        0.0374195678161, 0.0355377892706, 0.0333383902828, 0.0308286549265, 0.0282914440969,
        0.0257860718304, 0.02341635055, 0.0213126338243, 0.0195035612803, 0.0181079838989,
        0.0171991315458, 0.0166377598339, 0.0166445341361, 0.0171943735369, 0.0181980997278,
        0.0191339792146, 0.0198518804356, 0.0199714909193, 0.0194616474094, 0.0178626975229,
        0.0153296785464, 0.0126789254325, 0.0100766041988, 0.00773867100481, 0.00592386091874,
        0.00434706240169, 0.00310217013427, 0.00213213401899, 0.0013996000761, 0.000879148859271,
        0.000540866009427, 0.000326115560156, 0.000193965828516, 0.000114607606623, 6.74262828734e-05,
        3.97805301078e-05, 2.19948704638e-05, 9.72007976207e-06, 4.26179259146e-06, 2.80015581327e-06,
        1.14675436465e-06, 2.52452411995e-07, 9.08394910044e-08, 1.14291987912e-08, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0
    };


    for(int i=0; i<99; ++i){
        h_winter -> SetBinContent(Nbins[i]+1, pu_winter[i]);
        h_UltraLegacy -> SetBinContent(Nbins[i]+1, pu_ultralegacy[i]);
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

    h_UltraLegacy->GetYaxis()->SetTitleOffset(1.4);
    h_UltraLegacy->SetStats(0);
    h_UltraLegacy->SetLineWidth(2);
    h_UltraLegacy->SetLineColor(kRed);
    h_UltraLegacy->Draw();

    h_mc_Npu_afterPUreweighting->SetStats(0);
    h_mc_Npu_afterPUreweighting->SetLineWidth(2);
    h_mc_Npu_afterPUreweighting->SetLineColor(kGray);
    h_mc_Npu_afterPUreweighting->Draw("h,same");

    h_winter->SetStats(0);
    h_winter->SetLineWidth(2);
    h_winter->SetLineColor(kBlue);
    h_winter->Draw("same");

    h_mc_Npu->SetStats(0);
    h_mc_Npu->SetLineWidth(2);
    h_mc_Npu->SetLineColor(kBlack);
    h_mc_Npu->Draw("h,same");

    TLegend *legend = new TLegend(0.50,0.55,0.85,0.85);
    legend->AddEntry(h_winter, "Winter MC", "l");
    legend->AddEntry(h_UltraLegacy, "UtraLegacy MC", "l");
    legend->AddEntry(h_mc_Npu, "ggH(#gamma#gamma) before pu", "l");
    legend->AddEntry(h_mc_Npu_afterPUreweighting, "ggH(#gamma#gamma) after pu", "l");
    legend->SetLineColor(0);
    legend->Draw("same");

    c1->SaveAs("plots/puComparison.png");
    fin_pu->Close();
    fin->Close();



    //=== Check all 2017 dataset ===//
    TFile *f[N_dataset] = {
        //#--------------- TT Signal ---------------#
        TFile::Open("plots/TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8/hist_TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8.root"),
        TFile::Open("plots/TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8/hist_TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8.root"),
        TFile::Open("plots/TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8/hist_TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8.root"),
        TFile::Open("plots/TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8/hist_TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8.root"),
        TFile::Open("plots/TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8/hist_TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root"),
        TFile::Open("plots/TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8/hist_TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root"),
        TFile::Open("plots/TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8/hist_TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root"),
        TFile::Open("plots/TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8/hist_TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root"),
        //#--------------- Resonant bkg ---------------#
        TFile::Open("plots/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/hist_GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8.root"),
        TFile::Open("plots/VBFHToGG_M125_13TeV_amcatnlo_pythia8/hist_VBFHToGG_M125_13TeV_amcatnlo_pythia8.root"),
        TFile::Open("plots/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/hist_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root"),
        TFile::Open("plots/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/hist_ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root"),
        //#--------------- non-Resonant bkg ---------------#
        TFile::Open("plots/DiPhotonJetsBox_M40_80-Sherpa/hist_DiPhotonJetsBox_M40_80-Sherpa.root"),
        TFile::Open("plots/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa/hist_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa.root"),
        TFile::Open("plots/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/hist_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root"),
        TFile::Open("plots/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/hist_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root"),
        TFile::Open("plots/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8/hist_GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8.root"),
        TFile::Open("plots/QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/hist_QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root"),
        TFile::Open("plots/QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8/hist_QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8.root"),
        TFile::Open("plots/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/hist_QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root"),
        TFile::Open("plots/TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8/hist_TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8.root"),
        TFile::Open("plots/TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8/hist_TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8.root"),
        TFile::Open("plots/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/hist_TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root"),
        //#--------------- ST Signal ---------------#
        TFile::Open("plots/ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8/hist_ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root"),
        TFile::Open("plots/ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8/hist_ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8.root"),
        TFile::Open("plots/ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8/hist_ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root"),
        TFile::Open("plots/ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8/hist_ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root")
    };

    Int_t color[N_dataset] = {
        kMagenta, kMagenta-1, kMagenta-2, kMagenta-3, kMagenta-4, kMagenta-5,  kMagenta-6, kMagenta-7,
        kGray, kGray-1, kGray-2, kGray-3,
        kGreen, kGreen-1, kGreen-2, kGreen-3, kGreen-4, kGreen-5,  kGreen-6, kGreen-7, kGreen-8,  kGreen-9, kGreen-10,
        kOrange, kOrange-1, kOrange-2, kOrange-3
    };

    Int_t color_nonres[N_dataset] = {
        kMagenta, kMagenta-1, kMagenta-2, kMagenta-3, kMagenta-4, kMagenta-5,  kMagenta-6, kMagenta-7,//not used
        kGray, kGray-1, kGray-2, kGray-3,//not used
        kViolet, kViolet-1, kGreen, kGreen+1, kGreen+2, kOrange+1, kOrange+2, kOrange+3, kGray+1,  kGray+2, kGray+3,
        kOrange, kOrange-1, kOrange-2, kOrange-3//not used
    };


    //--------------------------------------------------
    printf("\n[INFO] Start to make plots (all)!\n");
    double total[N_dataset];
    TH1D* h[N_dataset];
    for(int i=0; i<N_dataset; ++i){
        printf("%d/%d\r", i+1, N_dataset);
        h[i] = (TH1D*)f[i]->Get("hist_NPu");
        total[i] = SumContents(h[i]);
        Normalization(h[i], total[i]);
        h[i]->SetLineColor(color[i]);
    }

    h_UltraLegacy->SetMaximum(0.1);
    h_UltraLegacy->Draw();
    h_winter->Draw("same");
    for(int i=0; i<N_dataset; ++i) h[i]->Draw("h,same");
    legend->Clear();
    legend->AddEntry(h[0],  "TT signal", "l");
    legend->AddEntry(h[23], "ST signal", "l");
    legend->AddEntry(h[8],  "Res. bkg", "l");
    legend->AddEntry(h[12], "Non-Res. bkg", "l");
    legend->Draw("same");

    TLegend *leg = new TLegend(0.15,0.70,0.50,0.85);
    leg->SetLineColor(0);
    leg->AddEntry(h_winter, "Winter MC", "l");
    leg->AddEntry(h_UltraLegacy, "UtraLegacy MC", "l");
    leg->Draw("same");

    h_UltraLegacy->Draw("same");
    h_winter->Draw("same");
    c1->SaveAs("plots/puComparison_all.png");

    
    //--------------------------------------------------
    printf("\n[INFO] Start to make plots (non-res)!\n");
    for(int i=12; i<23; ++i) h[i]->SetLineColor(color_nonres[i]);
    h_UltraLegacy->SetMaximum(0.1);
    h_UltraLegacy->Draw();
    h_winter->Draw("same");
    for(int i=12; i<23; ++i) h[i]->Draw("h,same");
    legend->Clear();
    legend->AddEntry(h[12], "DiphotonJetBox", "l");
    legend->AddEntry(h[14], "GJet", "l");
    legend->AddEntry(h[17], "QCD", "l");
    legend->AddEntry(h[20], "Others non-res", "l");
    legend->Draw("same");

    leg->Draw("same");

    h_UltraLegacy->Draw("same");
    h_winter->Draw("same");
    c1->SaveAs("plots/puComparison_nonres.png");


    for(int i=0; i<N_dataset; ++i) f[i]->Close();
}


double SumContents(TH1D* &h){
    double Ntotal = 0.;
    for(int i=0; i<99; ++i){
        Ntotal += h->GetBinContent(i+1);
    }
    return Ntotal;
}


void Normalization(TH1D* &h, double Ntotal){
    for(int i=0; i<99; ++i){
        double content = h->GetBinContent(i+1);
        double error = h->GetBinError(i+1);
        h->SetBinContent(i+1, content/Ntotal);
        h->SetBinError(i+1, error/Ntotal);
    }
}
