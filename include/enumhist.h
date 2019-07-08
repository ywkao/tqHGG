#ifndef __ENUMHIST_H__
#define __ENUMHIST_H__
#include <string>
using namespace std;

enum histList {
    hist_EvtInfo_NPu,
    hist_EvtInfo_Rho,
    hist_EvtInfo_NVtx,
    hist_EvtInfo_NVtx_wopu,
    hist_EvtInfo_genweight,
    //------------------------
    hist_DiPhoInfo_mass,
    hist_DiPhoInfo_pt,
    hist_DiPhoInfo_eta,
    hist_DiPhoInfo_phi,
    hist_DiPhoInfo_energy,
    hist_DiPhoInfo_leadPt,
    hist_DiPhoInfo_leadEta,
    hist_DiPhoInfo_leadPhi,
    hist_DiPhoInfo_leadE,
    hist_DiPhoInfo_leadhoe,
    hist_DiPhoInfo_leadIDMVA,
    hist_DiPhoInfo_subleadPt,
    hist_DiPhoInfo_subleadEta,
    hist_DiPhoInfo_subleadPhi,
    hist_DiPhoInfo_subleadE,
    hist_DiPhoInfo_subleadhoe,
    hist_DiPhoInfo_subleadIDMVA,
    //------------------------
    hist_ElecInfo_Size,
    hist_MuonInfo_Size,
    hist_num_leptons,// # of selected objects.
    hist_num_electrons,// # of selected objects.
    hist_num_muons,// # of selected objects.
    hist_ElecInfo_electron_pt,
    hist_ElecInfo_electron_eta,
    hist_ElecInfo_electron_phi,
    hist_ElecInfo_electron_energy,
    hist_ElecInfo_electron_diphoton_deltaR,
    hist_MuonInfo_muon_pt,
    hist_MuonInfo_muon_eta,
    hist_MuonInfo_muon_phi,
    hist_MuonInfo_muon_energy,
    hist_MuonInfo_muon_diphoton_deltaR,
    //------------------------
    hist_jets_size,
    hist_num_jets,
    hist_JetInfo_jet_pt,
    hist_JetInfo_jet_eta,
    hist_JetInfo_jet_phi,
    hist_JetInfo_jet_energy,
    hist_JetInfo_jet_diphoton_deltaR,
    //------------------------
    //Define in selection stage
    //------------------------
    hist_lepton_pt,
    hist_lepton_eta,
    hist_lepton_phi,
    hist_lepton_energy,
    hist_lepton_diphoton_deltaR,
    //------------------------
    hist_jet1_pt,
    hist_jet1_eta,
    hist_jet1_phi,
    hist_jet1_energy,
    hist_jet1_diphoton_deltaR,
    hist_jet2_pt,
    hist_jet2_eta,
    hist_jet2_phi,
    hist_jet2_energy,
    hist_jet2_diphoton_deltaR,
    //------------------------
    hist_jet1_lepton_deltaR,
    hist_jet2_lepton_deltaR,
    //------------------------
    hist_leading_bjet_pt,
    hist_leading_bjet_eta,
    hist_leading_bjet_phi,
    hist_leading_bjet_energy,
    hist_chosen_bjet_pt,
    hist_chosen_bjet_eta,
    hist_chosen_bjet_phi,
    hist_chosen_bjet_energy,
    //------------------------
    hist_inv_mass_dijet,
    hist_inv_mass_diphoton,
    hist_inv_mass_tbw,
    totalHistNum
};
std::string histNames[totalHistNum]{
    "hist_EvtInfo_NPu",
    "hist_EvtInfo_Rho",
    "hist_EvtInfo_NVtx",
    "hist_EvtInfo_NVtx_wopu",
    "hist_EvtInfo_genweight",
    //------------------------
    "hist_DiPhoInfo_mass",
    "hist_DiPhoInfo_pt",
    "hist_DiPhoInfo_eta",
    "hist_DiPhoInfo_phi",
    "hist_DiPhoInfo_energy",
    "hist_DiPhoInfo_leadPt",
    "hist_DiPhoInfo_leadEta",
    "hist_DiPhoInfo_leadPhi",
    "hist_DiPhoInfo_leadE",
    "hist_DiPhoInfo_leadhoe",
    "hist_DiPhoInfo_leadIDMVA",
    "hist_DiPhoInfo_subleadPt",
    "hist_DiPhoInfo_subleadEta",
    "hist_DiPhoInfo_subleadPhi",
    "hist_DiPhoInfo_subleadE",
    "hist_DiPhoInfo_subleadhoe",
    "hist_DiPhoInfo_subleadIDMVA",
    //------------------------
    "hist_ElecInfo_Size",
    "hist_MuonInfo_Size",
    "hist_num_leptons",// # of selected objects.
    "hist_num_electrons",// # of selected objects.
    "hist_num_muons",// # of selected objects.
    "hist_ElecInfo_electron_pt",
    "hist_ElecInfo_electron_eta",
    "hist_ElecInfo_electron_phi",
    "hist_ElecInfo_electron_energy",
    "hist_ElecInfo_electron_diphoton_deltaR",
    "hist_MuonInfo_muon_pt",
    "hist_MuonInfo_muon_eta",
    "hist_MuonInfo_muon_phi",
    "hist_MuonInfo_muon_energy",
    "hist_MuonInfo_muon_diphoton_deltaR",
    //------------------------
    "hist_jets_size",
    "hist_num_jets",
    "hist_JetInfo_jet_pt",
    "hist_JetInfo_jet_eta",
    "hist_JetInfo_jet_phi",
    "hist_JetInfo_jet_energy",
    "hist_JetInfo_jet_diphoton_deltaR",
    //------------------------
    //Define in selection stage
    //------------------------
    "hist_lepton_pt",
    "hist_lepton_eta",
    "hist_lepton_phi",
    "hist_lepton_energy",
    "hist_lepton_diphoton_deltaR",
    //------------------------
    "hist_jet1_pt",
    "hist_jet1_eta",
    "hist_jet1_phi",
    "hist_jet1_energy",
    "hist_jet1_diphoton_deltaR",
    "hist_jet2_pt",
    "hist_jet2_eta",
    "hist_jet2_phi",
    "hist_jet2_energy",
    "hist_jet2_diphoton_deltaR",
    //------------------------
    "hist_jet1_lepton_deltaR",
    "hist_jet2_lepton_deltaR",
    //------------------------
    "hist_leading_bjet_pt",
    "hist_leading_bjet_eta",
    "hist_leading_bjet_phi",
    "hist_leading_bjet_energy",
    "hist_chosen_bjet_pt",
    "hist_chosen_bjet_eta",
    "hist_chosen_bjet_phi",
    "hist_chosen_bjet_energy",
    //------------------------
    "hist_inv_mass_dijet",
    "hist_inv_mass_diphoton",
    "hist_inv_mass_tbw"
};
int histNbins[totalHistNum]{
    100,//hist_EvtInfo_NPu
    50,//hist_EvtInfo_Rho
    100,//hist_EvtInfo_NVtx
    100,//hist_EvtInfo_NVtx_wopu
    100,//hist_EvtInfo_genweight
    //------------------------
    40,//hist_DiPhoInfo_mass
    40,//hist_DiPhoInfo_pt
    40,//hist_DiPhoInfo_eta
    40,//hist_DiPhoInfo_phi
    40,//hist_DiPhoInfo_energy
    40,//hist_DiPhoInfo_leadPt
    40,//hist_DiPhoInfo_leadEta
    40,//hist_DiPhoInfo_leadPhi
    40,//hist_DiPhoInfo_leadE
    50,//hist_DiPhoInfo_leadhoe
    50,//hist_DiPhoInfo_leadIDMVA
    40,//hist_DiPhoInfo_subleadPt
    40,//hist_DiPhoInfo_subleadEta
    40,//hist_DiPhoInfo_subleadPhi
    40,//hist_DiPhoInfo_subleadE
    50,//hist_DiPhoInfo_subleadhoe
    50,//hist_DiPhoInfo_subleadIDMVA
    //------------------------
    10,//hist_ElecInfo_Size
    10,//hist_MuonInfo_Size
    10,//hist_num_leptons// # of selected objects.
    10,//hist_num_electrons// # of selected objects.
    10,//hist_num_muons// # of selected objects.
    40,//hist_ElecInfo_electron_pt
    40,//hist_ElecInfo_electron_eta
    40,//hist_ElecInfo_electron_phi
    40,//hist_ElecInfo_electron_energy
    40,//hist_ElecInfo_electron_diphoton_deltaR
    40,//hist_MuonInfo_muon_pt
    40,//hist_MuonInfo_muon_eta
    40,//hist_MuonInfo_muon_phi
    40,//hist_MuonInfo_muon_energy
    40,//hist_MuonInfo_muon_diphoton_deltaR
    //------------------------
    10,//hist_jets_size
    10,//hist_num_jets
    40,//hist_JetInfo_jet_pt
    40,//hist_JetInfo_jet_eta
    40,//hist_JetInfo_jet_phi
    40,//hist_JetInfo_jet_energy
    40,//hist_JetInfo_jet_diphoton_deltaR
    //------------------------
    //Define in selection stage
    //------------------------
    40,//hist_lepton_pt
    40,//hist_lepton_eta
    40,//hist_lepton_phi
    40,//hist_lepton_energy
    40,//hist_lepton_diphoton_deltaR
    //------------------------
    40,//hist_jet1_pt
    40,//hist_jet1_eta
    40,//hist_jet1_phi
    40,//hist_jet1_energy
    40,//hist_jet1_diphoton_deltaR
    40,//hist_jet2_pt
    40,//hist_jet2_eta
    40,//hist_jet2_phi
    40,//hist_jet2_energy
    40,//hist_jet2_diphoton_deltaR
    //------------------------
    40,//hist_jet1_lepton_deltaR
    40,//hist_jet2_lepton_deltaR
    //------------------------
    40,//hist_leading_bjet_pt
    40,//hist_leading_bjet_eta
    40,//hist_leading_bjet_phi
    40,//hist_leading_bjet_energy
    40,//hist_chosen_bjet_pt
    40,//hist_chosen_bjet_eta
    40,//hist_chosen_bjet_phi
    40,//hist_chosen_bjet_energy
    //------------------------
    50,//hist_inv_mass_dijet
    50,//hist_inv_mass_diphoton
    50//hist_inv_mass_tbw
};
double histBinLow[totalHistNum]{
    0,//hist_EvtInfo_NPu
    0,//hist_EvtInfo_Rho
    0,//hist_EvtInfo_NVtx
    0,//hist_EvtInfo_NVtx_wopu
    0,//hist_EvtInfo_genweight
    //------------------------
    100,//hist_DiPhoInfo_mass
    0,//hist_DiPhoInfo_pt
    -2.5,//hist_DiPhoInfo_eta
    -3.0,//hist_DiPhoInfo_phi
    100,//hist_DiPhoInfo_energy
    0,//hist_DiPhoInfo_leadPt
    -2.5,//hist_DiPhoInfo_leadEta
    -3.0,//hist_DiPhoInfo_leadPhi
    0,//hist_DiPhoInfo_leadE
    0,//hist_DiPhoInfo_leadhoe
    -1.,//hist_DiPhoInfo_leadIDMVA
    0,//hist_DiPhoInfo_subleadPt
    -2.5,//hist_DiPhoInfo_subleadEta
    -3.0,//hist_DiPhoInfo_subleadPhi
    0,//hist_DiPhoInfo_subleadE
    0,//hist_DiPhoInfo_subleadhoe
    -1.,//hist_DiPhoInfo_subleadIDMVA
    //------------------------
    0,//hist_ElecInfo_Size
    0,//hist_MuonInfo_Size
    0,//hist_num_leptons// # of selected objects.
    0,//hist_num_electrons// # of selected objects.
    0,//hist_num_muons// # of selected objects.
    0,//hist_ElecInfo_electron_pt
    -2.5,//hist_ElecInfo_electron_eta
    -3.0,//hist_ElecInfo_electron_phi
    0,//hist_ElecInfo_electron_energy
    0,//hist_ElecInfo_electron_diphoton_deltaR
    0,//hist_MuonInfo_muon_pt
    -2.5,//hist_MuonInfo_muon_eta
    -3.0,//hist_MuonInfo_muon_phi
    0,//hist_MuonInfo_muon_energy
    0,//hist_MuonInfo_muon_diphoton_deltaR
    //------------------------
    0,//hist_jets_size
    0,//hist_num_jets
    0,//hist_JetInfo_jet_pt
    -2.5,//hist_JetInfo_jet_eta
    -3.0,//hist_JetInfo_jet_phi
    0,//hist_JetInfo_jet_energy
    0,//hist_JetInfo_jet_diphoton_deltaR
    //------------------------
    //Define in selection stage
    //------------------------
    0,//hist_lepton_pt
    -2.5,//hist_lepton_eta
    -3.0,//hist_lepton_phi
    0,//hist_lepton_energy
    0,//hist_lepton_diphoton_deltaR
    //------------------------
    0,//hist_jet1_pt
    -2.5,//hist_jet1_eta
    -3.0,//hist_jet1_phi
    0,//hist_jet1_energy
    0,//hist_jet1_diphoton_deltaR
    0,//hist_jet2_pt
    -2.5,//hist_jet2_eta
    -3.0,//hist_jet2_phi
    0,//hist_jet2_energy
    0,//hist_jet2_diphoton_deltaR
    //------------------------
    0,//hist_jet1_lepton_deltaR
    0,//hist_jet2_lepton_deltaR
    //------------------------
    0,//hist_leading_bjet_pt
    -2.5,//hist_leading_bjet_eta
    -3.0,//hist_leading_bjet_phi
    0,//hist_leading_bjet_energy
    0,//hist_chosen_bjet_pt
    -2.5,//hist_chosen_bjet_eta
    -3.0,//hist_chosen_bjet_phi
    0,//hist_chosen_bjet_energy
    //------------------------
    55,//hist_inv_mass_dijet
    100,//hist_inv_mass_diphoton
    0//hist_inv_mass_tbw
};
double histBinHigh[totalHistNum]{
    100,//hist_EvtInfo_NPu
    200,//hist_EvtInfo_Rho
    100,//hist_EvtInfo_NVtx
    100,//hist_EvtInfo_NVtx_wopu
    100,//hist_EvtInfo_genweight
    //------------------------
    180,//hist_DiPhoInfo_mass
    200,//hist_DiPhoInfo_pt
    2.5,//hist_DiPhoInfo_eta
    3.0,//hist_DiPhoInfo_phi
    300,//hist_DiPhoInfo_energy
    200,//hist_DiPhoInfo_leadPt
    2.5,//hist_DiPhoInfo_leadEta
    3.0,//hist_DiPhoInfo_leadPhi
    200,//hist_DiPhoInfo_leadE
    0.10,//hist_DiPhoInfo_leadhoe
    1.,//hist_DiPhoInfo_leadIDMVA
    200,//hist_DiPhoInfo_subleadPt
    2.5,//hist_DiPhoInfo_subleadEta
    3.0,//hist_DiPhoInfo_subleadPhi
    200,//hist_DiPhoInfo_subleadE
    0.10,//hist_DiPhoInfo_subleadhoe
    1.,//hist_DiPhoInfo_subleadIDMVA
    //------------------------
    10,//hist_ElecInfo_Size
    10,//hist_MuonInfo_Size
    10,//hist_num_leptons// # of selected objects.
    10,//hist_num_electrons// # of selected objects.
    10,//hist_num_muons// # of selected objects.
    200,//hist_ElecInfo_electron_pt
    2.5,//hist_ElecInfo_electron_eta
    3.0,//hist_ElecInfo_electron_phi
    200,//hist_ElecInfo_electron_energy
    6,//hist_ElecInfo_electron_diphoton_deltaR
    200,//hist_MuonInfo_muon_pt
    2.5,//hist_MuonInfo_muon_eta
    3.0,//hist_MuonInfo_muon_phi
    200,//hist_MuonInfo_muon_energy
    6,//hist_MuonInfo_muon_diphoton_deltaR
    //------------------------
    10,//hist_jets_size
    10,//hist_num_jets
    200,//hist_JetInfo_jet_pt
    2.5,//hist_JetInfo_jet_eta
    3.0,//hist_JetInfo_jet_phi
    200,//hist_JetInfo_jet_energy
    6,//hist_JetInfo_jet_diphoton_deltaR
    //------------------------
    //Define in selection stage
    //------------------------
    200,//hist_lepton_pt
    2.5,//hist_lepton_eta
    3.0,//hist_lepton_phi
    200,//hist_lepton_energy
    6,//hist_lepton_diphoton_deltaR
    //------------------------
    200,//hist_jet1_pt
    2.5,//hist_jet1_eta
    3.0,//hist_jet1_phi
    200,//hist_jet1_energy
    6,//hist_jet1_diphoton_deltaR
    200,//hist_jet2_pt
    2.5,//hist_jet2_eta
    3.0,//hist_jet2_phi
    200,//hist_jet2_energy
    6,//hist_jet2_diphoton_deltaR
    //------------------------
    6,//hist_jet1_lepton_deltaR
    6,//hist_jet2_lepton_deltaR
    //------------------------
    200,//hist_leading_bjet_pt
    2.5,//hist_leading_bjet_eta
    3.0,//hist_leading_bjet_phi
    200,//hist_leading_bjet_energy
    200,//hist_chosen_bjet_pt
    2.5,//hist_chosen_bjet_eta
    3.0,//hist_chosen_bjet_phi
    200,//hist_chosen_bjet_energy
    //------------------------
    105,//hist_inv_mass_dijet
    150,//hist_inv_mass_diphoton
    500//hist_inv_mass_tbw
};
#endif
/*
    TH1D *hist_EvtInfo_NPu = new TH1D("hist_EvtInfo_NPu", "", 100, 0, 100);
    TH1D *hist_EvtInfo_Rho = new TH1D("hist_EvtInfo_Rho", "", 50, 0, 1000);
    TH1D *hist_EvtInfo_NVtx = new TH1D("hist_EvtInfo_NVtx", "", 100, 0, 100);
    TH1D *hist_EvtInfo_NVtx_wopu = new TH1D("hist_EvtInfo_NVtx_wopu", "", 100, 0, 100);
    TH1D *hist_EvtInfo_genweight = new TH1D("hist_EvtInfo_genweight", "", 100, 0, 100);
    //------------------------
    TH1D *hist_DiPhoInfo_mass = new TH1D("hist_DiPhoInfo_mass", "", 40, 100, 180);
    TH1D *hist_DiPhoInfo_pt = new TH1D("hist_DiPhoInfo_pt", "", 30, 0, 600);
    TH1D *hist_DiPhoInfo_eta = new TH1D("hist_DiPhoInfo_eta", "", 40, -2.5, 2.5);
    TH1D *hist_DiPhoInfo_phi = new TH1D("hist_DiPhoInfo_phi", "", 40, -3.0, 3.0);
    TH1D *hist_DiPhoInfo_energy = new TH1D("hist_DiPhoInfo_energy", "", 40, 0, 200);
    TH1D *hist_DiPhoInfo_leadPt = new TH1D("hist_DiPhoInfo_leadPt", "", 30, 0, 600);
    TH1D *hist_DiPhoInfo_leadEta = new TH1D("hist_DiPhoInfo_leadEta", "", 40, -2.5, 2.5);
    TH1D *hist_DiPhoInfo_leadPhi = new TH1D("hist_DiPhoInfo_leadPhi", "", 40, -3.0, 3.0);
    TH1D *hist_DiPhoInfo_leadE = new TH1D("hist_DiPhoInfo_leadE", "", 40, 0, 200);
    TH1D *hist_DiPhoInfo_leadhoe = new TH1D("hist_DiPhoInfo_leadhoe", "", 50, 0, 0.25);
    TH1D *hist_DiPhoInfo_leadIDMVA = new TH1D("hist_DiPhoInfo_leadIDMVA", "", 50, -1., 1.);
    TH1D *hist_DiPhoInfo_subleadPt = new TH1D("hist_DiPhoInfo_subleadPt", "", 30, 0, 600);
    TH1D *hist_DiPhoInfo_subleadEta = new TH1D("hist_DiPhoInfo_subleadEta", "", 40, -2.5, 2.5);
    TH1D *hist_DiPhoInfo_subleadPhi = new TH1D("hist_DiPhoInfo_subleadPhi", "", 40, -3.0, 3.0);
    TH1D *hist_DiPhoInfo_subleadE = new TH1D("hist_DiPhoInfo_subleadE", "", 40, 0, 200);
    TH1D *hist_DiPhoInfo_subleadhoe = new TH1D("hist_DiPhoInfo_subleadhoe", "", 50, 0, 0.25);
    TH1D *hist_DiPhoInfo_subleadIDMVA = new TH1D("hist_DiPhoInfo_subleadIDMVA", "", 50, -1., 1.);
    //------------------------
    TH1D *hist_ElecInfo_Size = new TH1D("hist_ElecInfo_Size", "", 10, 0, 10);
    TH1D *hist_MuonInfo_Size = new TH1D("hist_MuonInfo_Size", "", 10, 0, 10);
    TH1D *hist_num_leptons = new TH1D("hist_num_leptons", "", 10, 0, 10);// # of selected objects.
    TH1D *hist_num_electrons = new TH1D("hist_num_electrons", "", 10, 0, 10);// # of selected objects.
    TH1D *hist_num_muons = new TH1D("hist_num_muons", "", 10, 0, 10);// # of selected objects.
    TH1D *hist_ElecInfo_electron_pt = new TH1D("hist_ElecInfo_electron_pt", "", 50, 0, 1000);
    TH1D *hist_ElecInfo_electron_eta = new TH1D("hist_ElecInfo_electron_eta", "", 40, -2.5, 2.5);
    TH1D *hist_ElecInfo_electron_phi = new TH1D("hist_ElecInfo_electron_phi", "", 40, -3.0, 3.0);
    TH1D *hist_ElecInfo_electron_energy = new TH1D("hist_ElecInfo_electron_energy", "", 40, 0, 1000);
    TH1D *hist_ElecInfo_electron_diphoton_deltaR = new TH1D("hist_ElecInfo_electron_diphoton_deltaR", "", 60, 0, 60);
    TH1D *hist_MuonInfo_muon_pt = new TH1D("hist_MuonInfo_muon_pt", "", 50, 0, 1000);
    TH1D *hist_MuonInfo_muon_eta = new TH1D("hist_MuonInfo_muon_eta", "", 40, -2.5, 2.5);
    TH1D *hist_MuonInfo_muon_phi = new TH1D("hist_MuonInfo_muon_phi", "", 40, -3.0, 3.0);
    TH1D *hist_MuonInfo_muon_energy = new TH1D("hist_MuonInfo_muon_energy", "", 40, 0, 1000);
    TH1D *hist_MuonInfo_muon_diphoton_deltaR = new TH1D("hist_MuonInfo_muon_diphoton_deltaR", "", 60, 0, 60);
    //------------------------
    TH1D *hist_jets_size = new TH1D("hist_jets_size", "", 10, 0, 10);
    TH1D *hist_num_jets = new TH1D("hist_num_jets", "", 10, 0, 10);
    TH1D *hist_JetInfo_jet_pt = new TH1D("hist_JetInfo_jet_pt", "", 50, 0, 1000);
    TH1D *hist_JetInfo_jet_eta = new TH1D("hist_JetInfo_jet_eta", "", 40, -2.5, 2.5);
    TH1D *hist_JetInfo_jet_phi = new TH1D("hist_JetInfo_jet_phi", "", 40, -3.0, 3.0);
    TH1D *hist_JetInfo_jet_energy = new TH1D("hist_JetInfo_jet_energy", "", 40, 0, 1000);
    TH1D *hist_JetInfo_jet_diphoton_deltaR = new TH1D("hist_JetInfo_jet_diphoton_deltaR", "", 60, 0, 60);
    //------------------------
    //Define in selection stage
    //------------------------
    TH1D *hist_lepton_pt = new TH1D("hist_lepton_pt", "", 50, 0, 1000);
    TH1D *hist_lepton_eta = new TH1D("hist_lepton_eta", "", 40, -2.5, 2.5);
    TH1D *hist_lepton_phi = new TH1D("hist_lepton_phi", "", 40, -3.0, 3.0);
    TH1D *hist_lepton_energy = new TH1D("hist_lepton_energy", "", 40, 0, 1000);
    TH1D *hist_lepton_diphoton_deltaR = new TH1D("hist_lepton_diphoton_deltaR", "", 60, 0, 60);
    //------------------------
    TH1D *hist_jet1_diphoton_deltaR = new TH1D("hist_jet1_diphoton_deltaR", "", 60, 0, 60);
    TH1D *hist_jet2_diphoton_deltaR = new TH1D("hist_jet2_diphoton_deltaR", "", 60, 0, 60);
    TH1D *hist_jet1_lepton_deltaR = new TH1D("hist_jet1_lepton_deltaR", "", 60, 0, 60);
    TH1D *hist_jet2_lepton_deltaR = new TH1D("hist_jet2_lepton_deltaR", "", 60, 0, 60);
    //------------------------
    TH1D *hist_jet1_pt = new TH1D("hist_jet1_pt", "", 50, 0, 1000);
    TH1D *hist_jet1_eta = new TH1D("hist_jet1_eta", "", 40, -2.5, 2.5);
    TH1D *hist_jet1_phi = new TH1D("hist_jet1_phi", "", 40, -3.0, 3.0);
    TH1D *hist_jet1_energy = new TH1D("hist_jet1_energy", "", 40, 0, 1000);
    TH1D *hist_jet2_pt = new TH1D("hist_jet2_pt", "", 50, 0, 1000);
    TH1D *hist_jet2_eta = new TH1D("hist_jet2_eta", "", 40, -2.5, 2.5);
    TH1D *hist_jet2_phi = new TH1D("hist_jet2_phi", "", 40, -3.0, 3.0);
    TH1D *hist_jet2_energy = new TH1D("hist_jet2_energy", "", 40, 0, 1000);
    //------------------------
    TH1D  *hist_leading_bjet_pt = new TH1D("hist_leading_bjet_pt", "hist_leading_bjet_pt", 50, 0, 1000);
    TH1D  *hist_leading_bjet_eta = new TH1D("hist_leading_bjet_eta", "hist_leading_bjet_eta", 40, -2.5, 2.5);
    TH1D  *hist_leading_bjet_phi = new TH1D("hist_leading_bjet_phi", "hist_leading_bjet_phi", 40, -3.0, 3.0);
    TH1D  *hist_leading_bjet_energy = new TH1D("hist_leading_bjet_energy", "hist_leading_bjet_energy", 40, 0, 1000);
    TH1D  *hist_chosen_bjet_pt = new TH1D("hist_chosen_bjet_pt", "hist_chosen_bjet_pt", 50, 0, 1000);
    TH1D  *hist_chosen_bjet_eta = new TH1D("hist_chosen_bjet_eta", "hist_chosen_bjet_eta", 40, -2.5, 2.5);
    TH1D  *hist_chosen_bjet_phi = new TH1D("hist_chosen_bjet_phi", "hist_chosen_bjet_phi", 40, -3.0, 3.0);
    TH1D  *hist_chosen_bjet_energy = new TH1D("hist_chosen_bjet_energy", "hist_chosen_bjet_energy", 40, 0, 1000);
    //------------------------
    TH1D  *hist_inv_mass_dijet = new TH1D("hist_inv_mass_dijet", "hist_inv_mass_dijet", 50, 55, 105);
    TH1D  *hist_inv_mass_diphoton = new TH1D("hist_inv_mass_diphoton", "hist_inv_mass_diphoton", 50, 100, 150);
    TH1D  *hist_inv_mass_tbw = new TH1D("hist_inv_mass_tbw", "hist_inv_mass_tbw", 50, 0, 500);
*/
