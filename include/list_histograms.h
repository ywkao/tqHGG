//***************************************************************************
//
// FileName    : histogramList.h
// Purpose     : Simplify source code & Maintain histogram in an easier way.
// Usage       : Use ./script/createEnumHist to convert this file into ./include/enumhist.h
// Author      : Yu-Wei Kao [ykao@cern.ch]
//
//***************************************************************************
TH1D *hist_EvtInfo_NPu = new TH1D("hist_EvtInfo_NPu", "", 100, 0, 100);
TH1D *hist_EvtInfo_Rho = new TH1D("hist_EvtInfo_Rho", "", 40, 0, 80);
TH1D *hist_EvtInfo_Rho_wopu = new TH1D("hist_EvtInfo_Rho_wopu", "", 40, 0, 80);
TH1D *hist_EvtInfo_NVtx = new TH1D("hist_EvtInfo_NVtx", "", 50, 0, 100);
TH1D *hist_EvtInfo_NVtx_wopu = new TH1D("hist_EvtInfo_NVtx_wopu", "", 50, 0, 100);
TH1D *hist_EvtInfo_genweight = new TH1D("hist_EvtInfo_genweight", "", 100, 0, 100);
//------------------------
TH1D *hist_DiPhoInfo_mass = new TH1D("hist_DiPhoInfo_mass", "", 40, 100, 180);
TH1D *hist_DiPhoInfo_pt = new TH1D("hist_DiPhoInfo_pt", "", 40, 0, 200);
TH1D *hist_DiPhoInfo_eta = new TH1D("hist_DiPhoInfo_eta", "", 40, -2.5, 2.5);
TH1D *hist_DiPhoInfo_phi = new TH1D("hist_DiPhoInfo_phi", "", 40, -3.0, 3.0);
TH1D *hist_DiPhoInfo_energy = new TH1D("hist_DiPhoInfo_energy", "", 40, 0, 500);
TH1D *hist_DiPhoInfo_leadPt = new TH1D("hist_DiPhoInfo_leadPt", "", 40, 0, 200);
TH1D *hist_DiPhoInfo_leadEta = new TH1D("hist_DiPhoInfo_leadEta", "", 40, -2.5, 2.5);
TH1D *hist_DiPhoInfo_leadPhi = new TH1D("hist_DiPhoInfo_leadPhi", "", 40, -3.0, 3.0);
TH1D *hist_DiPhoInfo_leadE = new TH1D("hist_DiPhoInfo_leadE", "", 40, 0, 200);
TH1D *hist_DiPhoInfo_leadhoe = new TH1D("hist_DiPhoInfo_leadhoe", "", 50, 0, 0.10);
TH1D *hist_DiPhoInfo_leadIDMVA = new TH1D("hist_DiPhoInfo_leadIDMVA", "", 50, -1., 1.);
TH1D *hist_DiPhoInfo_subleadPt = new TH1D("hist_DiPhoInfo_subleadPt", "", 40, 0, 200);
TH1D *hist_DiPhoInfo_subleadEta = new TH1D("hist_DiPhoInfo_subleadEta", "", 40, -2.5, 2.5);
TH1D *hist_DiPhoInfo_subleadPhi = new TH1D("hist_DiPhoInfo_subleadPhi", "", 40, -3.0, 3.0);
TH1D *hist_DiPhoInfo_subleadE = new TH1D("hist_DiPhoInfo_subleadE", "", 40, 0, 200);
TH1D *hist_DiPhoInfo_subleadhoe = new TH1D("hist_DiPhoInfo_subleadhoe", "", 50, 0, 0.10);
TH1D *hist_DiPhoInfo_subleadIDMVA = new TH1D("hist_DiPhoInfo_subleadIDMVA", "", 50, -1., 1.);
//------------------------
TH1D *hist_ElecInfo_Size = new TH1D("hist_ElecInfo_Size", "", 10, 0, 10);
TH1D *hist_MuonInfo_Size = new TH1D("hist_MuonInfo_Size", "", 10, 0, 10);
TH1D *hist_num_leptons = new TH1D("hist_num_leptons", "", 10, 0, 10);
TH1D *hist_num_electrons = new TH1D("hist_num_electrons", "", 10, 0, 10);
TH1D *hist_num_muons = new TH1D("hist_num_muons", "", 10, 0, 10);
TH1D *hist_ElecInfo_electron_charge = new TH1D("hist_ElecInfo_electron_charge", "", 2, -2, 2);
TH1D *hist_ElecInfo_electron_pt = new TH1D("hist_ElecInfo_electron_pt", "", 40, 0, 200);
TH1D *hist_ElecInfo_electron_eta = new TH1D("hist_ElecInfo_electron_eta", "", 40, -2.5, 2.5);
TH1D *hist_ElecInfo_electron_phi = new TH1D("hist_ElecInfo_electron_phi", "", 40, -3.0, 3.0);
TH1D *hist_ElecInfo_electron_energy = new TH1D("hist_ElecInfo_electron_energy", "", 40, 0, 200);
TH1D *hist_ElecInfo_electron_diphoton_deltaR = new TH1D("hist_ElecInfo_electron_diphoton_deltaR", "", 40, 0, 6);
TH1D *hist_MuonInfo_muon_charge = new TH1D("hist_MuonInfo_muon_charge", "", 2, -2, 2);
TH1D *hist_MuonInfo_muon_pt = new TH1D("hist_MuonInfo_muon_pt", "", 40, 0, 200);
TH1D *hist_MuonInfo_muon_eta = new TH1D("hist_MuonInfo_muon_eta", "", 40, -2.5, 2.5);
TH1D *hist_MuonInfo_muon_phi = new TH1D("hist_MuonInfo_muon_phi", "", 40, -3.0, 3.0);
TH1D *hist_MuonInfo_muon_energy = new TH1D("hist_MuonInfo_muon_energy", "", 40, 0, 200);
TH1D *hist_MuonInfo_muon_diphoton_deltaR = new TH1D("hist_MuonInfo_muon_diphoton_deltaR", "", 40, 0, 6);
//------------------------
TH1D *hist_jets_size = new TH1D("hist_jets_size", "", 10, 0, 10);
TH1D *hist_num_jets = new TH1D("hist_num_jets", "", 15, 0, 15);
TH1D *hist_JetInfo_jet_pt = new TH1D("hist_JetInfo_jet_pt", "", 40, 0, 200);
TH1D *hist_JetInfo_jet_eta = new TH1D("hist_JetInfo_jet_eta", "", 40, -2.5, 2.5);
TH1D *hist_JetInfo_jet_phi = new TH1D("hist_JetInfo_jet_phi", "", 40, -3.0, 3.0);
TH1D *hist_JetInfo_jet_energy = new TH1D("hist_JetInfo_jet_energy", "", 40, 0, 200);
TH1D *hist_JetInfo_jet_diphoton_deltaR = new TH1D("hist_JetInfo_jet_diphoton_deltaR", "", 40, 0, 6);
//------------------------
TH1D *hist_MetInfo_Pt = new TH1D("hist_MetInfo_Pt", "", 40, 0, 200);
TH1D *hist_MetInfo_Phi = new TH1D("hist_MetInfo_Phi", "", 40, -3.0, 3.0);
TH1D *hist_MetInfo_Px = new TH1D("hist_MetInfo_Px", "", 40, -200, 200);
TH1D *hist_MetInfo_Py = new TH1D("hist_MetInfo_Py", "", 40, -200, 200);
TH1D *hist_MetInfo_SumET = new TH1D("hist_MetInfo_SumET", "", 40, 0, 200);
TH1D *hist_MetInfo_Pz_solution_1 = new TH1D("hist_MetInfo_Pz_solution_1", "", 40, -200, 200);
TH1D *hist_MetInfo_Pz_solution_2 = new TH1D("hist_MetInfo_Pz_solution_2", "", 40, -200, 200);
//TH1D *hist_MetInfo_coeff_D = new TH1D("hist_MetInfo_coeff_D", "", 40, -1.e+6, 1e+6);
TH1D *hist_MetInfo_coeff_D = new TH1D("hist_MetInfo_coeff_D", "", 40, -1000, 1000);
TH1D *hist_MetInfo_coeff_A = new TH1D("hist_MetInfo_coeff_A", "", 40, -1000., 1000.);
TH1D *hist_MetInfo_coeff_B2A = new TH1D("hist_MetInfo_coeff_B2A", "", 40, -600, 600);
TH1D *hist_MetInfo_coeff_D2A = new TH1D("hist_MetInfo_coeff_D2A", "", 40, -600, 600);
//------------------------
//Define in selection stage
//------------------------
TH1D *hist_lepton_charge = new TH1D("hist_lepton_charge", "", 2, -2, 2);
TH1D *hist_lepton_pt = new TH1D("hist_lepton_pt", "", 40, 0, 200);
TH1D *hist_lepton_eta = new TH1D("hist_lepton_eta", "", 40, -2.5, 2.5);
TH1D *hist_lepton_phi = new TH1D("hist_lepton_phi", "", 40, -3.0, 3.0);
TH1D *hist_lepton_energy = new TH1D("hist_lepton_energy", "", 40, 0, 200);
TH1D *hist_lepton_diphoton_deltaR = new TH1D("hist_lepton_diphoton_deltaR", "", 40, 0, 6);
//------------------------
//jet1: leading jet
//jet2: subleading jet
//wjet1(2): jet constitute w_candidate // had process
//jetq: jet chosen to be the one in tqh coupling //sig tt samples
//------------------------
TH1D *hist_jet1_pt = new TH1D("hist_jet1_pt", "", 40, 0, 200);
TH1D *hist_jet1_eta = new TH1D("hist_jet1_eta", "", 40, -2.5, 2.5);
TH1D *hist_jet1_phi = new TH1D("hist_jet1_phi", "", 40, -3.0, 3.0);
TH1D *hist_jet1_energy = new TH1D("hist_jet1_energy", "", 40, 0, 200);
TH1D *hist_jet1_btag_score = new TH1D("hist_jet1_btag_score", "", 10, 0, 1);
TH1D *hist_jet1_CvsL_score = new TH1D("hist_jet1_CvsL_score", "", 10, 0, 1);
TH1D *hist_jet1_CvsB_score = new TH1D("hist_jet1_CvsB_score", "", 10, 0, 1);
TH1D *hist_jet1_diphoton_deltaR = new TH1D("hist_jet1_diphoton_deltaR", "", 40, 0, 6);
TH1D *hist_jet1_lepton_deltaR = new TH1D("hist_jet1_lepton_deltaR", "", 40, 0, 6);
//---
TH1D *hist_jet2_pt = new TH1D("hist_jet2_pt", "", 40, 0, 200);
TH1D *hist_jet2_eta = new TH1D("hist_jet2_eta", "", 40, -2.5, 2.5);
TH1D *hist_jet2_phi = new TH1D("hist_jet2_phi", "", 40, -3.0, 3.0);
TH1D *hist_jet2_energy = new TH1D("hist_jet2_energy", "", 40, 0, 200);
TH1D *hist_jet2_btag_score = new TH1D("hist_jet2_btag_score", "", 10, 0, 1);
TH1D *hist_jet2_CvsL_score = new TH1D("hist_jet2_CvsL_score", "", 10, 0, 1);
TH1D *hist_jet2_CvsB_score = new TH1D("hist_jet2_CvsB_score", "", 10, 0, 1);
TH1D *hist_jet2_diphoton_deltaR = new TH1D("hist_jet2_diphoton_deltaR", "", 40, 0, 6);
TH1D *hist_jet2_lepton_deltaR = new TH1D("hist_jet2_lepton_deltaR", "", 40, 0, 6);
//---
TH1D *hist_wjet1_pt = new TH1D("hist_wjet1_pt", "", 40, 0, 200);
TH1D *hist_wjet1_eta = new TH1D("hist_wjet1_eta", "", 40, -2.5, 2.5);
TH1D *hist_wjet1_phi = new TH1D("hist_wjet1_phi", "", 40, -3.0, 3.0);
TH1D *hist_wjet1_energy = new TH1D("hist_wjet1_energy", "", 40, 0, 200);
TH1D *hist_wjet1_btag_score = new TH1D("hist_wjet1_btag_score", "", 10, 0, 1);
TH1D *hist_wjet1_CvsL_score = new TH1D("hist_wjet1_CvsL_score", "", 10, 0, 1);
TH1D *hist_wjet1_CvsB_score = new TH1D("hist_wjet1_CvsB_score", "", 10, 0, 1);
TH1D *hist_wjet1_diphoton_deltaR = new TH1D("hist_wjet1_diphoton_deltaR", "", 40, 0, 6);
TH1D *hist_wjet1_lepton_deltaR = new TH1D("hist_wjet1_lepton_deltaR", "", 40, 0, 6);
//---
TH1D *hist_wjet2_pt = new TH1D("hist_wjet2_pt", "", 40, 0, 200);
TH1D *hist_wjet2_eta = new TH1D("hist_wjet2_eta", "", 40, -2.5, 2.5);
TH1D *hist_wjet2_phi = new TH1D("hist_wjet2_phi", "", 40, -3.0, 3.0);
TH1D *hist_wjet2_energy = new TH1D("hist_wjet2_energy", "", 40, 0, 200);
TH1D *hist_wjet2_btag_score = new TH1D("hist_wjet2_btag_score", "", 10, 0, 1);
TH1D *hist_wjet2_CvsL_score = new TH1D("hist_wjet2_CvsL_score", "", 10, 0, 1);
TH1D *hist_wjet2_CvsB_score = new TH1D("hist_wjet2_CvsB_score", "", 10, 0, 1);
TH1D *hist_wjet2_diphoton_deltaR = new TH1D("hist_wjet2_diphoton_deltaR", "", 40, 0, 6);
TH1D *hist_wjet2_lepton_deltaR = new TH1D("hist_wjet2_lepton_deltaR", "", 40, 0, 6);
//---
TH1D *hist_jetq_pt = new TH1D("hist_jetq_pt", "", 40, 0, 200);
TH1D *hist_jetq_eta = new TH1D("hist_jetq_eta", "", 40, -2.5, 2.5);
TH1D *hist_jetq_phi = new TH1D("hist_jetq_phi", "", 40, -3.0, 3.0);
TH1D *hist_jetq_energy = new TH1D("hist_jetq_energy", "", 40, 0, 200);
TH1D *hist_jetq_btag_score = new TH1D("hist_jetq_btag_score", "", 10, 0, 1);
TH1D *hist_jetq_CvsL_score = new TH1D("hist_jetq_CvsL_score", "", 10, 0, 1);
TH1D *hist_jetq_CvsB_score = new TH1D("hist_jetq_CvsB_score", "", 10, 0, 1);
TH1D *hist_jetq_diphoton_deltaR = new TH1D("hist_jetq_diphoton_deltaR", "", 40, 0, 6);
TH1D *hist_jetq_lepton_deltaR = new TH1D("hist_jetq_lepton_deltaR", "", 40, 0, 6);
//------------------------
TH1D  *hist_leading_bjet_pt = new TH1D("hist_leading_bjet_pt", "hist_leading_bjet_pt", 40, 0, 200);
TH1D  *hist_leading_bjet_eta = new TH1D("hist_leading_bjet_eta", "hist_leading_bjet_eta", 40, -2.5, 2.5);
TH1D  *hist_leading_bjet_phi = new TH1D("hist_leading_bjet_phi", "hist_leading_bjet_phi", 40, -3.0, 3.0);
TH1D  *hist_leading_bjet_energy = new TH1D("hist_leading_bjet_energy", "hist_leading_bjet_energy", 40, 0, 200);
TH1D  *hist_leading_bjet_btag_score = new TH1D("hist_leading_bjet_btag_score", "", 10, 0, 1);
TH1D  *hist_leading_bjet_CvsL_score = new TH1D("hist_leading_bjet_CvsL_score", "", 10, 0, 1);
TH1D  *hist_leading_bjet_CvsB_score = new TH1D("hist_leading_bjet_CvsB_score", "", 10, 0, 1);
//------------------------
TH1D *hist_deltaR_top_top = new TH1D("hist_deltaR_top_top", "", 40, 0, 6);
TH1D *hist_deltaR_qH = new TH1D("hist_deltaR_qH", "", 40, 0, 6);
TH1D *hist_deltaR_photon_photon = new TH1D("hist_deltaR_photon_photon", "", 40, 0, 6);
//---
TH1D *hist_deltaR_bW = new TH1D("hist_deltaR_bW", "", 40, 0, 6);
TH1D *hist_deltaR_HW = new TH1D("hist_deltaR_HW", "", 40, 0, 6);
TH1D *hist_deltaR_tH = new TH1D("hist_deltaR_tH", "", 40, 0, 6);
//---
TH1D *hist_deltaR_lep_met = new TH1D("hist_deltaR_lep_met", "", 40, 0, 6);
TH1D *hist_deltaR_jet1_jet2 = new TH1D("hist_deltaR_jet1_jet2", "", 40, 0, 6);
TH1D *hist_deltaR_wjet1_wjet2 = new TH1D("hist_deltaR_wjet1_wjet2", "", 40, 0, 6);
//------------------------
TH1D *hist_top_tqh_pt = new TH1D("hist_top_tqh_pt", "", 40, 0, 400);
TH1D *hist_top_tqh_eta = new TH1D("hist_top_tqh_eta", "", 40, -2.5, 2.5);
TH1D *hist_top_tqh_mass = new TH1D("hist_top_tqh_mass", "", 40, 0, 400);
//---
TH1D *hist_hadronic_w_candidate_pt = new TH1D("hist_hadronic_w_candidate_pt", "", 40, 0, 200);
TH1D *hist_hadronic_w_candidate_eta = new TH1D("hist_hadronic_w_candidate_eta", "", 40, -2.5, 2.5);
TH1D *hist_hadronic_w_candidate_mass = new TH1D("hist_hadronic_w_candidate_mass", "", 50, 0, 150);
TH1D *hist_hadronic_top_tbw_pt = new TH1D("hist_hadronic_top_tbw_pt", "", 40, 0, 400);
TH1D *hist_hadronic_top_tbw_eta = new TH1D("hist_hadronic_top_tbw_eta", "", 40, -2.5, 2.5);
TH1D *hist_hadronic_top_tbw_mass = new TH1D("hist_hadronic_top_tbw_mass", "", 40, 0, 400);
//---
TH1D *hist_leptonic_w_candidate_solution1_pt = new TH1D("hist_leptonic_w_candidate_solution1_pt", "", 40, 0, 200);
TH1D *hist_leptonic_w_candidate_solution1_eta = new TH1D("hist_leptonic_w_candidate_solution1_eta", "", 40, -2.5, 2.5);
TH1D *hist_leptonic_w_candidate_solution1_mass = new TH1D("hist_leptonic_w_candidate_solution1_mass", "", 50, 0, 150);
TH1D *hist_leptonic_top_tbw_solution1_pt = new TH1D("hist_leptonic_top_tbw_solution1_pt", "", 40, 0, 400);
TH1D *hist_leptonic_top_tbw_solution1_eta = new TH1D("hist_leptonic_top_tbw_solution1_eta", "", 40, -2.5, 2.5);
TH1D *hist_leptonic_top_tbw_solution1_mass = new TH1D("hist_leptonic_top_tbw_solution1_mass", "", 40, 0, 400);
//---
TH1D *hist_leptonic_w_candidate_solution2_pt = new TH1D("hist_leptonic_w_candidate_solution2_pt", "", 40, 0, 200);
TH1D *hist_leptonic_w_candidate_solution2_eta = new TH1D("hist_leptonic_w_candidate_solution2_eta", "", 40, -2.5, 2.5);
TH1D *hist_leptonic_w_candidate_solution2_mass = new TH1D("hist_leptonic_w_candidate_solution2_mass", "", 50, 0, 150);
TH1D *hist_leptonic_top_tbw_solution2_pt = new TH1D("hist_leptonic_top_tbw_solution2_pt", "", 40, 0, 400);
TH1D *hist_leptonic_top_tbw_solution2_eta = new TH1D("hist_leptonic_top_tbw_solution2_eta", "", 40, -2.5, 2.5);
TH1D *hist_leptonic_top_tbw_solution2_mass = new TH1D("hist_leptonic_top_tbw_solution2_mass", "", 40, 0, 400);
//---
TH1D *hist_leptonic_w_candidate_topKinFit_pt = new TH1D("hist_leptonic_w_candidate_topKinFit_pt", "", 40, 0, 200);
TH1D *hist_leptonic_w_candidate_topKinFit_eta = new TH1D("hist_leptonic_w_candidate_topKinFit_eta", "", 40, -2.5, 2.5);
TH1D *hist_leptonic_w_candidate_topKinFit_mass = new TH1D("hist_leptonic_w_candidate_topKinFit_mass", "", 50, 0, 150);
TH1D *hist_leptonic_top_tbw_topKinFit_pt = new TH1D("hist_leptonic_top_tbw_topKinFit_pt", "", 40, 0, 400);
TH1D *hist_leptonic_top_tbw_topKinFit_eta = new TH1D("hist_leptonic_top_tbw_topKinFit_eta", "", 40, -2.5, 2.5);
TH1D *hist_leptonic_top_tbw_topKinFit_mass = new TH1D("hist_leptonic_top_tbw_topKinFit_mass", "", 40, 0, 400);
//------------------------
//TH1D  *hist_mass_w_candidate = new TH1D("hist_mass_w_candidate", "hist_mass_w_candidate", 50, 0, 150);
//TH1D  *hist_mass_top_candidate = new TH1D("hist_mass_top_candidate", "hist_mass_top_candidate", 40, 0, 400);
//TH1D  *hist_mass_top_fcnh_candidate = new TH1D("hist_mass_top_fcnh_candidate", "hist_mass_top_fcnh_candidate", 40, 0, 400);
