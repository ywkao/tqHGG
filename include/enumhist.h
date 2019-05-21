#ifndef __ENUMHIST_H__
#define __ENUMHIST_H__
#include <string>
using namespace std;

enum histList {
    hist_NPu,
    hist_NVtx,
    hist_NVtx_wopu,
    hist_num_jets,
    hist_num_btagged_jets,
    hist_num_nonbtagged_jets,
    hist_leading_bjet_pt,
    hist_leading_bjet_eta,
    hist_leading_bjet_phi,
    hist_chosen_bjet_pt,
    hist_chosen_bjet_eta,
    hist_chosen_bjet_phi,
    hist_jet1_pt,
    hist_jet1_eta,
    hist_jet1_phi,
    hist_jet2_pt,
    hist_jet2_eta,
    hist_jet2_phi,
    hist_dijet_eta,
    hist_dijet_phi,
    hist_dijet_angle,
    hist_inv_mass_dijet,
    hist_inv_mass_diphoton,
    hist_inv_mass_tbw,
    hist_DiPhoInfo_leadPt,
    hist_DiPhoInfo_leadEta,
    hist_DiPhoInfo_leadPhi,
    hist_DiPhoInfo_leadE,
    hist_DiPhoInfo_leadIDMVA,
    hist_DiPhoInfo_subleadPt,
    hist_DiPhoInfo_subleadEta,
    hist_DiPhoInfo_subleadPhi,
    hist_DiPhoInfo_subleadE,
    hist_DiPhoInfo_subleadIDMVA,
    hist_DiPhoInfo_leadIDMVA_ori,
    hist_DiPhoInfo_subleadIDMVA_ori,
    hist_inv_mass_diphoton_ori,
    totalHistNum
};

std::string histNames[totalHistNum]{
    "hist_NPu",
    "hist_NVtx",
    "hist_NVtx_wopu",
    "hist_num_jets",
    "hist_num_btagged_jets",
    "hist_num_nonbtagged_jets",
    "hist_leading_bjet_pt",
    "hist_leading_bjet_eta",
    "hist_leading_bjet_phi",
    "hist_chosen_bjet_pt",
    "hist_chosen_bjet_eta",
    "hist_chosen_bjet_phi",
    "hist_jet1_pt",
    "hist_jet1_eta",
    "hist_jet1_phi",
    "hist_jet2_pt",
    "hist_jet2_eta",
    "hist_jet2_phi",
    "hist_dijet_eta",
    "hist_dijet_phi",
    "hist_dijet_angle",
    "hist_inv_mass_dijet",
    "hist_inv_mass_diphoton",
    "hist_inv_mass_tbw",
    "hist_DiPhoInfo_leadPt",
    "hist_DiPhoInfo_leadEta",
    "hist_DiPhoInfo_leadPhi",
    "hist_DiPhoInfo_leadE",
    "hist_DiPhoInfo_leadIDMVA",
    "hist_DiPhoInfo_subleadPt",
    "hist_DiPhoInfo_subleadEta",
    "hist_DiPhoInfo_subleadPhi",
    "hist_DiPhoInfo_subleadE",
    "hist_DiPhoInfo_subleadIDMVA",
    "hist_DiPhoInfo_leadIDMVA_ori",
    "hist_DiPhoInfo_subleadIDMVA_ori",
    "hist_inv_mass_diphoton_ori"
};
int histNbins[totalHistNum]{
     100,//"hist_NPu"
     100,//"hist_NVtx"
     100,//"hist_NVtx_wopu"
     20,//"hist_num_jets"
     10,//"hist_num_btagged_jets"
     20,//"hist_num_nonbtagged_jets"
     50,//"hist_leading_bjet_pt"
     40,//"hist_leading_bjet_eta"
     40,//"hist_leading_bjet_phi"
     50,//"hist_chosen_bjet_pt"
     40,//"hist_chosen_bjet_eta"
     40,//"hist_chosen_bjet_phi"
     50,//"hist_jet1_pt"
     40,//"hist_jet1_eta"
     40,//"hist_jet1_phi"
     50,//"hist_jet2_pt"
     40,//"hist_jet2_eta"
     40,//"hist_jet2_phi"
     40,//"hist_dijet_eta"
     40,//"hist_dijet_phi"
     40,//"hist_dijet_angle"
     50,//"hist_inv_mass_dijet"
     50,//"hist_inv_mass_diphoton"
     50,//"hist_inv_mass_tbw"
     50,//"hist_DiPhoInfo_leadPt"
     40,//"hist_DiPhoInfo_leadEta"
     40,//"hist_DiPhoInfo_leadPhi"
     50,//"hist_DiPhoInfo_leadE"
     50,//"hist_DiPhoInfo_leadIDMVA"
     50,//"hist_DiPhoInfo_subleadPt"
     40,//"hist_DiPhoInfo_subleadEta"
     40,//"hist_DiPhoInfo_subleadPhi"
     50,//"hist_DiPhoInfo_subleadE"
     50,//"hist_DiPhoInfo_subleadIDMVA"
     50,//"hist_DiPhoInfo_leadIDMVA_ori"
     50,//"hist_DiPhoInfo_subleadIDMVA_ori"
     50 //"hist_inv_mass_diphoton_ori"
};
double histBinLow[totalHistNum]{
     0,//"hist_NPu"
     0,//"hist_NVtx"
     0,//"hist_NVtx_wopu"
     0,//"hist_num_jets"
     0,//"hist_num_btagged_jets"
     0,//"hist_num_nonbtagged_jets"
     0,//"hist_leading_bjet_pt"
     -2.5,//"hist_leading_bjet_eta"
     -3.0,//"hist_leading_bjet_phi"
     0,//"hist_chosen_bjet_pt"
     -2.5,//"hist_chosen_bjet_eta"
     -3.0,//"hist_chosen_bjet_phi"
     0,//"hist_jet1_pt"
     -2.5,//"hist_jet1_eta"
     -3.0,//"hist_jet1_phi"
     0,//"hist_jet2_pt"
     -2.5,//"hist_jet2_eta"
     -3.0,//"hist_jet2_phi"
     0.,//"hist_dijet_eta"
     0.,//"hist_dijet_phi"
     0.,//"hist_dijet_angle"
     55,//"hist_inv_mass_dijet"
     100,//"hist_inv_mass_diphoton"
     0,//"hist_inv_mass_tbw"
     0.,//"hist_DiPhoInfo_leadPt"
     -2.5,//"hist_DiPhoInfo_leadEta"
     -3.,//"hist_DiPhoInfo_leadPhi"
     0.,//"hist_DiPhoInfo_leadE"
     -1.,//"hist_DiPhoInfo_leadIDMVA"
     0.,//"hist_DiPhoInfo_subleadPt"
     -2.5,//"hist_DiPhoInfo_subleadEta"
     -3.,//"hist_DiPhoInfo_subleadPhi"
     0.,//"hist_DiPhoInfo_subleadE"
     -1.,//"hist_DiPhoInfo_subleadIDMVA"
     -1.,//"hist_DiPhoInfo_leadIDMVA_ori"
     -1.,//"hist_DiPhoInfo_subleadIDMVA_ori"
     100 //"hist_inv_mass_diphoton_ori"
};
double histBinHigh[totalHistNum]{
     100,//"hist_NPu"
     100,//"hist_NVtx"
     100,//"hist_NVtx_wopu"
     20,//"hist_num_jets"
     10,//"hist_num_btagged_jets"
     20,//"hist_num_nonbtagged_jets"
     1000,//"hist_leading_bjet_pt"
     2.5,//"hist_leading_bjet_eta"
     3.0,//"hist_leading_bjet_phi"
     1000,//"hist_chosen_bjet_pt"
     2.5,//"hist_chosen_bjet_eta"
     3.0,//"hist_chosen_bjet_phi"
     1000,//"hist_jet1_pt"
     2.5,//"hist_jet1_eta"
     3.0,//"hist_jet1_phi"
     1000,//"hist_jet2_pt"
     2.5,//"hist_jet2_eta"
     3.0,//"hist_jet2_phi"
     5.0,//"hist_dijet_eta"
     6.0,//"hist_dijet_phi"
     8.0,//"hist_dijet_angle"
     105,//"hist_inv_mass_dijet"
     150,//"hist_inv_mass_diphoton"
     500,//"hist_inv_mass_tbw"
     1000,//"hist_DiPhoInfo_leadPt"
     2.5,//"hist_DiPhoInfo_leadEta"
     3.,//"hist_DiPhoInfo_leadPhi"
     1000,//"hist_DiPhoInfo_leadE"
     1.,//"hist_DiPhoInfo_leadIDMVA"
     1000,//"hist_DiPhoInfo_subleadPt"
     2.5,//"hist_DiPhoInfo_subleadEta"
     3.,//"hist_DiPhoInfo_subleadPhi"
     1000,//"hist_DiPhoInfo_subleadE"
     1.,//"hist_DiPhoInfo_subleadIDMVA"
     1.,//"hist_DiPhoInfo_leadIDMVA_ori"
     1.,//"hist_DiPhoInfo_subleadIDMVA_ori"
     150 //"hist_inv_mass_diphoton_ori"
};
#endif
/*
    TH1D  *hist_NPu = new TH1D("hist_NPu", "hist_NPu", 100, 0, 100);
    TH1D  *hist_NVtx = new TH1D("hist_NVtx", "hist_NVtx", 100, 0, 100);
    TH1D  *hist_NVtx_wopu = new TH1D("hist_NVtx_wopu", "hist_NVtx_wopu", 100, 0, 100);
    TH1D  *hist_num_jets = new TH1D("hist_num_jets", "hist_num_jets", 20, 0, 20);
    TH1D  *hist_num_btagged_jets = new TH1D("hist_num_btagged_jets", "hist_num_btagged_jets", 10, 0, 10);
    TH1D  *hist_num_nonbtagged_jets = new TH1D("hist_num_nonbtagged_jets", "hist_num_nonbtagged_jets", 20, 0, 20);
    //------------------------
    TH1D  *hist_leading_bjet_pt = new TH1D("hist_leading_bjet_pt", "hist_leading_bjet_pt", 50, 0, 1000);
    TH1D  *hist_leading_bjet_eta = new TH1D("hist_leading_bjet_eta", "hist_leading_bjet_eta", 40, -2.5, 2.5);
    TH1D  *hist_leading_bjet_phi = new TH1D("hist_leading_bjet_phi", "hist_leading_bjet_phi", 40, -3.0, 3.0);
    TH1D  *hist_chosen_bjet_pt = new TH1D("hist_chosen_bjet_pt", "hist_chosen_bjet_pt", 50, 0, 1000);
    TH1D  *hist_chosen_bjet_eta = new TH1D("hist_chosen_bjet_eta", "hist_chosen_bjet_eta", 40, -2.5, 2.5);
    TH1D  *hist_chosen_bjet_phi = new TH1D("hist_chosen_bjet_phi", "hist_chosen_bjet_phi", 40, -3.0, 3.0);
    TH1D  *hist_jet1_pt = new TH1D("hist_jet1_pt", "hist_jet1_pt", 50, 0, 1000);
    TH1D  *hist_jet1_eta = new TH1D("hist_jet1_eta", "hist_jet1_eta", 40, -2.5, 2.5);
    TH1D  *hist_jet1_phi = new TH1D("hist_jet1_phi", "hist_jet1_phi", 40, -3.0, 3.0);
    TH1D  *hist_jet2_pt = new TH1D("hist_jet2_pt", "hist_jet2_pt", 50, 0, 1000);
    TH1D  *hist_jet2_eta = new TH1D("hist_jet2_eta", "hist_jet2_eta", 40, -2.5, 2.5);
    TH1D  *hist_jet2_phi = new TH1D("hist_jet2_phi", "hist_jet2_phi", 40, -3.0, 3.0);
    //------------------------
    TH1D  *hist_dijet_eta = new TH1D("hist_dijet_eta", "hist_dijet_eta", 40, 0., 5.0);
    TH1D  *hist_dijet_phi = new TH1D("hist_dijet_phi", "hist_dijet_phi", 40, 0., 6.0);
    TH1D  *hist_dijet_angle = new TH1D("hist_dijet_angle", "hist_dijet_angle", 40, 0., 8.0);
    TH1D  *hist_inv_mass_dijet = new TH1D("hist_inv_mass_dijet", "hist_inv_mass_dijet", 50, 55, 105);
    TH1D  *hist_inv_mass_diphoton = new TH1D("hist_inv_mass_diphoton", "hist_inv_mass_diphoton", 50, 100, 150);
    TH1D  *hist_inv_mass_tbw = new TH1D("hist_inv_mass_tbw", "hist_inv_mass_tbw", 50, 0, 500);
    //------------------------
    TH1D  *hist_DiPhoInfo_leadPt = new TH1D("hist_DiPhoInfo_leadPt", "DiPhoInfo_leadPt", 50, 0., 1000);
    TH1D  *hist_DiPhoInfo_leadEta = new TH1D("hist_DiPhoInfo_leadEta", "DiPhoInfo_leadEta", 40, -2.5, 2.5);
    TH1D  *hist_DiPhoInfo_leadPhi = new TH1D("hist_DiPhoInfo_leadPhi", "DiPhoInfo_leadPhi", 40, -3., 3.);
    TH1D  *hist_DiPhoInfo_leadE = new TH1D("hist_DiPhoInfo_leadE", "DiPhoInfo_leadE", 50, 0., 1000);
    TH1D  *hist_DiPhoInfo_leadIDMVA = new TH1D("hist_DiPhoInfo_leadIDMVA", "DiPhoInfo_leadIDMVA", 50, -1., 1.);
    TH1D  *hist_DiPhoInfo_subleadPt = new TH1D("hist_DiPhoInfo_subleadPt", "DiPhoInfo_subleadPt", 50, 0., 1000);
    TH1D  *hist_DiPhoInfo_subleadEta = new TH1D("hist_DiPhoInfo_subleadEta", "DiPhoInfo_subleadEta", 40, -2.5, 2.5);
    TH1D  *hist_DiPhoInfo_subleadPhi = new TH1D("hist_DiPhoInfo_subleadPhi", "DiPhoInfo_subleadPhi", 40, -3., 3.);
    TH1D  *hist_DiPhoInfo_subleadE = new TH1D("hist_DiPhoInfo_subleadE", "DiPhoInfo_subleadE", 50, 0., 1000);
    TH1D  *hist_DiPhoInfo_subleadIDMVA = new TH1D("hist_DiPhoInfo_subleadIDMVA", "DiPhoInfo_subleadIDMVA", 50, -1., 1.);
    //------------------------
    TH1D  *hist_DiPhoInfo_leadIDMVA_ori = new TH1D("hist_DiPhoInfo_leadIDMVA_ori", "DiPhoInfo_leadIDMVA", 50, -1., 1.);
    TH1D  *hist_DiPhoInfo_subleadIDMVA_ori = new TH1D("hist_DiPhoInfo_subleadIDMVA_ori", "DiPhoInfo_subleadIDMVA", 50, -1., 1.);
    TH1D  *hist_inv_mass_diphoton_ori = new TH1D("hist_inv_mass_diphoton_ori", "hist_inv_mass_diphoton", 50, 100, 150);
*/
