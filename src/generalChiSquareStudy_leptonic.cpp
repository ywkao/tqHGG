//***************************************************************************
//
// FileName    : generalChiSquareStudy.cpp
// Purpose     : Develop for top FCNH with H to two photons analysis
// Description : Extracting event info & Selecting objects (diphoton, leptons, jets).
// Author      : Yu-Wei Kao [ykao@cern.ch]
//
//***************************************************************************
//### includes{{{
#include <algorithm> //sts::fill_n(array, N_elements, -999)
#include <stdio.h>
#include <math.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TVectorD.h>
#include <TStyle.h>
#include <TMatrixD.h>
#include <TTree.h>
#include <vector>
#include <string>
#include "../include/generalChiSquareStudy.h"
#include "../include/cross_section.h"
#include "../../TopKinFit/kinfit.h"
#include "../../TopKinFit/TopLep.h"
using namespace std;
//}}}

int main(int argc, char *argv[]){
    //### I/O and event info{{{
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
    if(!isMultiFile) treeReader.AddSingleRootFile(input_file);
    else             treeReader.AddMultiRootFile(input_file);
    //treeReader.flashggStdTree->Add("/wk_cms2/youying/public/tH_FCNC/Era2017_RR-31Mar2018_v2/TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8.root");
    //treeReader.flashggStdTree->Add("/wk_cms2/youying/public/tH_FCNC/Era2017_RR-31Mar2018_v2/TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8.root");
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
    double CrossSection = GetXsec(dataset); //pb
    double BranchingFraction = GetBranchingFraction(dataset); //pb
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
    printf("[INFO] NormalizationFactor = %f!\n");
    //}}}
    //### histograms{{{
    //==================================================//
    //--------   Histograms for mass spectra  ----------//
    //==================================================//
    TH1D *hist_num_bjets_tight = new TH1D("hist_num_bjets_tight", ";Number of tight b-tagged jets;Entries", 10, 0, 10);
    TH1D *hist_num_bjets_medium = new TH1D("hist_num_bjets_medium", ";Number of medium b-tagged jets;Entries", 10, 0, 10);
    TH1D *hist_num_bjets_loose = new TH1D("hist_num_bjets_loose", ";Number of loose b-tagged jets;Entries", 10, 0, 10);
    TH1D *hist_num_bjets_tight_matched = new TH1D("hist_num_bjets_tight_matched", ";Number of tight b-tagged jets;Entries", 10, 0, 10);
    TH1D *hist_num_bjets_medium_matched = new TH1D("hist_num_bjets_medium_matched", ";Number of medium b-tagged jets;Entries", 10, 0, 10);
    TH1D *hist_num_bjets_loose_matched = new TH1D("hist_num_bjets_loose_matched", ";Number of loose b-tagged jets;Entries", 10, 0, 10);
    TH1D *hist_num_gen_bquark = new TH1D("hist_num_gen_bquark", ";Number of b quark (gen-level);Entries", 10, 0, 10);
    TH1D *hist_num_gen_light_quark = new TH1D("hist_num_gen_light_quark", ";Number of light quarks (gen-level);Entries", 10, 0, 10);
    TH1D *hist_num_selected_jets = new TH1D("hist_num_selected_jets", ";Number of selected jets;Entries", 10, 0, 10);
    TH1D *hist_deltaR_W_genW_simple = new TH1D("hist_deltaR_W_genW_simple", ";deltaR;Entries", 20, 0, 0.01);
    TH1D *hist_deltaR_W_genW_modified = new TH1D("hist_deltaR_W_genW_modified", ";deltaR;Entries", 20, 0, 0.01);
    TH1D *hist_deltaR_W_genW_yfyj = new TH1D("hist_deltaR_W_genW_yfyj", ";deltaR;Entries", 20, 0, 0.01);

    TH1D *hist_mass_diphoton = new TH1D("hist_mass_diphoton", ";Mass [GeV/c^2];Entries", 50, 0, 250);
    TH1D *hist_mass_top_fcnh_simple = new TH1D("hist_mass_top_fcnh_simple", ";Mass [GeV/c^2];Entries", 70, 0, 350);
    TH1D *hist_mass_top_fcnh_modified = new TH1D("hist_mass_top_fcnh_modified", ";Mass [GeV/c^2];Entries", 70, 0, 350);

    TH1D *hist_mass_w_candidate_yfyj = new TH1D("hist_mass_w_candidate_yfyj", ";Mass [GeV/c^2];Entries", 60, 0, 180);
    TH1D *hist_mass_t1_candidate_yfyj = new TH1D("hist_mass_t1_candidate_yfyj", ";Mass [GeV/c^2];Entries", 70, 0, 350);
    TH1D *hist_mass_t2_candidate_yfyj = new TH1D("hist_mass_t2_candidate_yfyj", ";Mass [GeV/c^2];Entries", 70, 0, 350);
    TH1D *hist_mass_gen_w_candidate_yfyj = new TH1D("hist_mass_gen_w_candidate_yfyj", ";Mass [GeV/c^2];Entries", 50, 0, 150);
    TH1D *hist_mass_gen_t2_candidate_yfyj = new TH1D("hist_mass_gen_t2_candidate_yfyj", ";Mass [GeV/c^2];Entries", 70, 0, 350);
    TH1D *hist_mass_conditioned_gen_w_candidate_yfyj    = new TH1D("hist_mass_conditioned_gen_w_candidate_yfyj", ";Mass [GeV/c^2];Entries", 50, 0, 150);
    TH1D *hist_mass_conditioned_gen_t2_candidate_yfyj   = new TH1D("hist_mass_conditioned_gen_t2_candidate_yfyj", ";Mass [GeV/c^2];Entries", 70, 0, 350);

    TH1D *hist_mass_w_candidate_chi2_simple = new TH1D("hist_mass_w_candidate_chi2_simple", ";Mass [GeV/c^2];Entries", 60, 0, 180);
    TH1D *hist_mass_top_candidate_chi2_simple = new TH1D("hist_mass_top_candidate_chi2_simple", ";Mass [GeV/c^2];Entries", 70, 0, 350);
    TH1D *hist_mass_w_candidate_chi2_modified = new TH1D("hist_mass_w_candidate_chi2_modified", ";Mass [GeV/c^2];Entries", 60, 0, 180);
    TH1D *hist_mass_top_candidate_chi2_modified = new TH1D("hist_mass_top_candidate_chi2_modified", ";Mass [GeV/c^2];Entries", 70, 0, 350);
    TH1D *hist_mass_gen_w_candidate_chi2_simple     = new TH1D("hist_mass_gen_w_candidate_chi2_simple", ";Mass [GeV/c^2];Entries", 50, 0, 150);
    TH1D *hist_mass_gen_top_candidate_chi2_simple   = new TH1D("hist_mass_gen_top_candidate_chi2_simple", ";Mass [GeV/c^2];Entries", 70, 0, 350);
    TH1D *hist_mass_gen_w_candidate_chi2_modified   = new TH1D("hist_mass_gen_w_candidate_chi2_modified", ";Mass [GeV/c^2];Entries", 50, 0, 150);
    TH1D *hist_mass_gen_top_candidate_chi2_modified = new TH1D("hist_mass_gen_top_candidate_chi2_modified", ";Mass [GeV/c^2];Entries", 70, 0, 350);

    TH1D *hist_mass_conditioned_gen_w_candidate_chi2_simple     = new TH1D("hist_mass_conditioned_gen_w_candidate_chi2_simple", ";Mass [GeV/c^2];Entries", 50, 0, 150);
    TH1D *hist_mass_conditioned_gen_top_candidate_chi2_simple   = new TH1D("hist_mass_conditioned_gen_top_candidate_chi2_simple", ";Mass [GeV/c^2];Entries", 70, 0, 350);
    TH1D *hist_mass_conditioned_gen_w_candidate_chi2_modified   = new TH1D("hist_mass_conditioned_gen_w_candidate_chi2_modified", ";Mass [GeV/c^2];Entries", 50, 0, 150);
    TH1D *hist_mass_conditioned_gen_top_candidate_chi2_modified = new TH1D("hist_mass_conditioned_gen_top_candidate_chi2_modified", ";Mass [GeV/c^2];Entries", 70, 0, 350);
    //}}}
    //### hist leptonic{{{
    //# deltaR{{{
    TH1D *hist_mass_gen_wboson_leptonic = new TH1D("hist_mass_gen_wboson_leptonic", ";Mass [GeV/c^{2}];Entries", 32, 0, 160);
    TH1D *hist_mass_gen_topquark_leptonic = new TH1D("hist_mass_gen_topquark_leptonic", ";Mass [GeV/c^{2}];Entries", 34, 0, 340);
    TH1D *hist_deltaR_reco_top_higgs_leptonic = new TH1D("hist_deltaR_reco_top_higgs_leptonic", "", 40, 0, 6);
    TH1D *hist_deltaR_gen_reco_chargedLepton = new TH1D("hist_deltaR_gen_reco_chargedLepton", "", 40, 0, 0.1);
    TH1D *hist_deltaR_gen_reco_neutrino_sol0 = new TH1D("hist_deltaR_gen_reco_neutrino_sol0", "", 40, 0, 6);
    TH1D *hist_deltaR_gen_reco_wboson_sol0 = new TH1D("hist_deltaR_gen_reco_wboson_sol0", "", 40, 0, 6);
    TH1D *hist_deltaR_gen_reco_neutrino_sol1 = new TH1D("hist_deltaR_gen_reco_neutrino_sol1", "", 40, 0, 6);
    TH1D *hist_deltaR_gen_reco_wboson_sol1 = new TH1D("hist_deltaR_gen_reco_wboson_sol1", "", 40, 0, 6);
    //}}}
    //# deltaPT{{{
    TH1D *hist_deltaPT_gen_reco_chargedLepton = new TH1D("hist_deltaPT_gen_reco_chargedLepton", "", 40, -1, 1);
    TH1D *hist_deltaPT_gen_reco_neutrino_sol0 = new TH1D("hist_deltaPT_gen_reco_neutrino_sol0", "", 40, -1, 1);
    TH1D *hist_deltaPT_gen_reco_wboson_sol0 = new TH1D("hist_deltaPT_gen_reco_wboson_sol0", "", 40, -1, 1);
    TH1D *hist_deltaPT_gen_reco_neutrino_sol1 = new TH1D("hist_deltaPT_gen_reco_neutrino_sol1", "", 40, -1, 1);
    TH1D *hist_deltaPT_gen_reco_wboson_sol1 = new TH1D("hist_deltaPT_gen_reco_wboson_sol1", "", 40, -1, 1);
    TH1D *hist_deltaPT_reco_neutrino_sol0_doubleCheck = new TH1D("hist_deltaPT_gen_reco_neutrino_sol0_doubleCheck", "", 40, -1, 1);
    TH1D *hist_deltaPT_reco_neutrino_sol1_doubleCheck = new TH1D("hist_deltaPT_gen_reco_neutrino_sol1_doubleCheck", "", 40, -1, 1);
    //}}}
    //---
    TH1D *hist_deltaPz_gen_reco_chargedLepton = new TH1D("hist_deltaPz_gen_reco_chargedLepton", "", 40, -100, 100);
    TH1D *hist_deltaPz_gen_reco_neutrino_sol0 = new TH1D("hist_deltaPz_gen_reco_neutrino_sol0", "", 40, -200, 200);
    TH1D *hist_deltaPz_gen_reco_neutrino_sol1 = new TH1D("hist_deltaPz_gen_reco_neutrino_sol1", "", 40, -200, 200);
    TH1D *hist_deltaPz_gen_reco_wboson_sol0 = new TH1D("hist_deltaPz_gen_reco_wboson_sol0", "", 40, -100, 100);
    TH1D *hist_deltaPz_gen_reco_wboson_sol1 = new TH1D("hist_deltaPz_gen_reco_wboson_sol1", "", 40, -100, 100);
    //---
    TH1D *hist_deltaPzRatio_gen_reco_neutrino_topKinFit = new TH1D("hist_deltaPzRatio_gen_reco_neutrino_topKinFit", ";[Pz(reco) - Pz(gen)] / Pz(gen); Entries", 40, -2, 2);
    TH1D *hist_deltaPzRatio_gen_reco_neutrino_topKinFit_regDisc = new TH1D("hist_deltaPzRatio_gen_reco_neutrino_topKinFit_regDisc", "", 40, -2, 2);
    TH1D *hist_deltaPzRatio_gen_reco_neutrino_topKinFit_irrDisc = new TH1D("hist_deltaPzRatio_gen_reco_neutrino_topKinFit_irrDisc", "", 40, -2, 2);
    TH2D *hist_reco_Pz_gen_Pz_topKinFit = new TH2D("hist_reco_Pz_gen_Pz_topKinFit", ";Gen-level Pz [GeV/c];Reco Pz [GeV/c]", 40, -200, 200, 40, -200, 200);
    TH2D *hist_reco_Pz_gen_Pz_topKinFit_regDisc = new TH2D("hist_reco_Pz_gen_Pz_topKinFit_regDisc", ";Gen-level Pz [GeV/c];Reco Pz [GeV/c]", 40, -200, 200, 40, -200, 200);
    TH2D *hist_reco_Pz_gen_Pz_topKinFit_irrDisc = new TH2D("hist_reco_Pz_gen_Pz_topKinFit_irrDisc", ";Gen-level Pz [GeV/c];Reco Pz [GeV/c]", 40, -200, 200, 40, -200, 200);
    //---
    //TH1D *hist_deltaPz_gen_reco_chargedLepton = new TH1D("hist_deltaPz_gen_reco_chargedLepton", "", 40, -1, 1);
    TH1D *hist_deltaPzRatio_gen_reco_neutrino_sol0 = new TH1D("hist_deltaPzRatio_gen_reco_neutrino_sol0", ";[Pz(reco) - Pz(gen)] / Pz(gen); Entries", 40, -2, 2);
    TH1D *hist_deltaPzRatio_gen_reco_neutrino_sol0_positiveD = new TH1D("hist_deltaPzRatio_gen_reco_neutrino_sol0_positiveD", "", 40, -2, 2);
    TH1D *hist_deltaPzRatio_gen_reco_neutrino_sol0_negativeD = new TH1D("hist_deltaPzRatio_gen_reco_neutrino_sol0_negativeD", "", 40, -2, 2);

    TH1D *hist_deltaPzRatio_gen_reco_neutrino_sol1 = new TH1D("hist_deltaPzRatio_gen_reco_neutrino_sol1", ";[Pz(reco) - Pz(gen)] / Pz(gen); Entries", 40, -2, 2);
    TH1D *hist_deltaPzRatio_gen_reco_neutrino_sol1_positiveD = new TH1D("hist_deltaPzRatio_gen_reco_neutrino_sol1_positiveD", "", 40, -2, 2);
    TH1D *hist_deltaPzRatio_gen_reco_neutrino_sol1_negativeD = new TH1D("hist_deltaPzRatio_gen_reco_neutrino_sol1_negativeD", "", 40, -2, 2);

    TH2D *hist_reco_Pz_gen_Pz = new TH2D("hist_reco_Pz_gen_Pz", ";Gen-level Pz [GeV/c];Reco Pz [GeV/c]", 40, -200, 200, 40, -200, 200);
    TH2D *hist_reco_Pz_gen_Pz_negativeD = new TH2D("hist_reco_Pz_gen_Pz_negativeD", ";Gen-level Pz [GeV/c];Reco Pz [GeV/c]", 40, -200, 200, 40, -200, 200);
    TH2D *hist_reco_Pz_gen_Pz_positiveD = new TH2D("hist_reco_Pz_gen_Pz_positiveD", ";Gen-level Pz [GeV/c];Reco Pz [GeV/c]", 40, -200, 200, 40, -200, 200);

    TH1D *hist_deltaPzRatio_gen_reco_wboson_sol0 = new TH1D("hist_deltaPzRatio_gen_reco_wboson_sol0", "", 40, -1, 1);
    TH1D *hist_deltaPzRatio_gen_reco_wboson_sol1 = new TH1D("hist_deltaPzRatio_gen_reco_wboson_sol1", "", 40, -1, 1);
    //---
    //# hist of neutrino pz{{{
    TH1D *hist_MetInfo_Pz_solution_1 = new TH1D("hist_MetInfo_Pz_solution_1", "", 40, -200, 200);
    TH1D *hist_MetInfo_Pz_solution_2 = new TH1D("hist_MetInfo_Pz_solution_2", "", 40, -200, 200);
    TH1D *hist_MetInfo_Pz_solution_1_negativeD = new TH1D("hist_MetInfo_Pz_solution_1_negativeD", "", 40, -200, 200);
    TH1D *hist_MetInfo_Pz_solution_2_negativeD = new TH1D("hist_MetInfo_Pz_solution_2_negativeD", "", 40, -200, 200);
    TH1D *hist_MetInfo_Pz_solution_1_positiveD = new TH1D("hist_MetInfo_Pz_solution_1_positiveD", "", 40, -200, 200);
    TH1D *hist_MetInfo_Pz_solution_2_positiveD = new TH1D("hist_MetInfo_Pz_solution_2_positiveD", "", 40, -200, 200);
    TH1D *hist_MetInfo_coeff_A = new TH1D("hist_MetInfo_coeff_A", "", 40, 0., 1.5);
    TH1D *hist_MetInfo_coeff_B = new TH1D("hist_MetInfo_coeff_B", "", 40, -1000, 1000);
    TH1D *hist_MetInfo_coeff_C = new TH1D("hist_MetInfo_coeff_C", "", 40, -1000, 1000);
    TH1D *hist_MetInfo_coeff_D = new TH1D("hist_MetInfo_coeff_D", "", 40, -1000, 1000);
    TH1D *hist_MetInfo_coeff_D_gen = new TH1D("hist_MetInfo_coeff_D_gen", "", 40, -1000, 1000);
    TH1D *hist_MetInfo_coeff_B2A = new TH1D("hist_MetInfo_coeff_B2A", "", 40, -600, 600);
    TH1D *hist_MetInfo_coeff_D2A = new TH1D("hist_MetInfo_coeff_D2A", "", 40, -600, 600);
    //---
    TH1D *hist_gen_neutrino_pz = new TH1D("hist_gen_neutrino_pz", "", 40, -200, 200);
    //---
    //}}}
    //# hist of top, w{{{
    TH1D *hist_leptonic_w_candidate_solution1_pt = new TH1D("hist_leptonic_w_candidate_solution1_pt", "", 40, 0, 200);
    TH1D *hist_leptonic_w_candidate_solution1_eta = new TH1D("hist_leptonic_w_candidate_solution1_eta", "", 40, -2.5, 2.5);
    TH1D *hist_leptonic_w_candidate_solution1_mass = new TH1D("hist_leptonic_w_candidate_solution1_mass", ";Mass [GeV/c^{2}];Entries", 50, 0, 150);
    TH1D *hist_leptonic_w_candidate_solution1_mass_negativeD = new TH1D("hist_leptonic_w_candidate_solution1_mass_negativeD", ";Mass [GeV/c^{2}];Entries", 50, 0, 150);
    TH1D *hist_leptonic_top_tbw_solution1_pt = new TH1D("hist_leptonic_top_tbw_solution1_pt", "", 40, 0, 400);
    TH1D *hist_leptonic_top_tbw_solution1_eta = new TH1D("hist_leptonic_top_tbw_solution1_eta", "", 40, -2.5, 2.5);
    TH1D *hist_leptonic_top_tbw_solution1_mass = new TH1D("hist_leptonic_top_tbw_solution1_mass", ";Mass [GeV/c^{2}];Entries", 70, 0, 350);
    TH1D *hist_leptonic_top_tbw_solution1_mass_negativeD = new TH1D("hist_leptonic_top_tbw_solution1_mass_negativeD", ";Mass [GeV/c^{2}];Entries", 70, 0, 350);
    //---
    TH1D *hist_leptonic_w_candidate_solution2_pt = new TH1D("hist_leptonic_w_candidate_solution2_pt", "", 40, 0, 200);
    TH1D *hist_leptonic_w_candidate_solution2_eta = new TH1D("hist_leptonic_w_candidate_solution2_eta", "", 40, -2.5, 2.5);
    TH1D *hist_leptonic_w_candidate_solution2_mass = new TH1D("hist_leptonic_w_candidate_solution2_mass", ";Mass [GeV/c^{2}];Entries", 50, 0, 150);
    TH1D *hist_leptonic_w_candidate_solution2_mass_negativeD = new TH1D("hist_leptonic_w_candidate_solution2_mass_negativeD", ";Mass [GeV/c^{2}];Entries", 50, 0, 150);
    TH1D *hist_leptonic_top_tbw_solution2_pt = new TH1D("hist_leptonic_top_tbw_solution2_pt", "", 40, 0, 400);
    TH1D *hist_leptonic_top_tbw_solution2_eta = new TH1D("hist_leptonic_top_tbw_solution2_eta", "", 40, -2.5, 2.5);
    TH1D *hist_leptonic_top_tbw_solution2_mass = new TH1D("hist_leptonic_top_tbw_solution2_mass", ";Mass [GeV/c^{2}];Entries", 70, 0, 350);
    TH1D *hist_leptonic_top_tbw_solution2_mass_negativeD = new TH1D("hist_leptonic_top_tbw_solution2_mass_negativeD", ";Mass [GeV/c^{2}];Entries", 70, 0, 350);
    //}}}

    //}}}
    //### Counters{{{
    int counter_irregular_disc = 0;
    int counter_coeff_D = 0;
    int counter_coeff_D_gen = 0;
    int counter_coeff_D_isNegative = 0;
    int counter_coeff_D_isNegative_gen = 0;
    int Nevents_pass_selection = 0;
    int check_num_bquak_is_one = 0;
    int check_num_bquak_is_two = 0;
    int check_num_bquak_is_three = 0;
    int count_bjet_is_bquark_loose = 0;
    int count_bjet_is_bquark_medium = 0;
    int count_bjet_is_bquark_tight = 0;
    int count_bjet_is_bquark = 0;
    int check_M1_20 = 0;
    int check_M20_simple = 0;
    int check_M20_modified = 0;
    int check_gen_exclude_id_999 = 0;
    int check_gen_exclude_the_same = 0;
    double accuracy_chi2_improved = 0;
    double accuracy_chi2_simple = 0, accuracy_chi2_modified = 0, accuracy_yfyj = 0;
    double accuracy_tbw_chi2_simple = 0, accuracy_tbw_chi2_modified = 0, accuracy_tbw_chi2_improved = 0, accuracy_tbw_yfyj = 0;
    double accuracy_tqh_chi2_simple = 0, accuracy_tqh_chi2_modified = 0, accuracy_tqh_chi2_improved = 0, accuracy_tqh_yfyj = 0;
    //}}}
    
    //##################################################//
    //#########    Event Loop [Selection]    ###########//
    //##################################################//
    // Goal: leptons, jets, diphoton; t->b+W(jj), 1 bjet + 2 chi2 jets 
    // # topKinFit method{{{
    KINFIT::kfit *kf = new KINFIT::kfit();
    kf->Init(TOPLEP); // Initialize tool for ttbar with FCNC top decay to Higgs(->bb)+u/c hypothesis
    kf->SetNToy(100); // Set number of toys for minimization
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
    //}}}
    for(int ientry=0; ientry<nentries; ientry++){
        treeReader.flashggStdTree->GetEntry(ientry);//load data
        if((ientry+1)%1000==0) printf("ientry = %d/%d\r", ientry+1, nentries+1);
        if((ientry+1)==nentries) printf("ientry = %d/%d\n", ientry+1, nentries+1);
        //# Pre-process{{{
        //### reset, selection, Normalization factor{{{ 
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
        bool pass_leadingPhotonPT = treeReader.DiPhoInfo_leadPt > treeReader.DiPhoInfo_mass / 2.;
        bool pass_subleadingPhotonPT = treeReader.DiPhoInfo_subleadPt > treeReader.DiPhoInfo_mass / 4.;
        bool pass_photon_criteria_pt = pass_leadingPhotonPT && pass_subleadingPhotonPT;

        bool pass_leadingPhotonEta =  (treeReader.DiPhoInfo_leadEta < 1.4442) || (treeReader.DiPhoInfo_leadEta > 1.566 && treeReader.DiPhoInfo_leadEta < 2.4);
        bool pass_subleadingPhotonEta = (treeReader.DiPhoInfo_subleadEta < 1.4442) || (treeReader.DiPhoInfo_subleadEta > 1.566 && treeReader.DiPhoInfo_subleadEta < 2.4);
        bool pass_photon_criteria_eta = pass_leadingPhotonEta && pass_subleadingPhotonEta;

        bool pass_photon_criteria = pass_photon_criteria_pt && pass_photon_criteria_eta;
        if(!pass_photon_criteria) continue;
        //Others
        //if(treeReader.DiPhoInfo_mass<0) continue;
        //if( !(treeReader.jets_size>0) ) continue;
        //if(!(DiPhoInfo_leadIDMVA>0)) continue;
        //==================================================//
        //------------   Normalization factor   ------------//
       //==================================================//
        mytree.EvtInfo_totalEntry_before_preselection = nentries;
        mytree.EvtInfo_NormalizationFactor_lumi = isData ? 1. : NormalizationFactor;
        //### }}} 
        //### diphoton info{{{
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
        //### }}}
        //### electron selection {{{
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
                if( fabs(treeReader.ElecInfo_Pt->at(i))  < 20  ) continue;
                //--- check deltaR(electron,photon) ---//
                TLorentzVector electron; 
                electron.SetPtEtaPhiE(treeReader.ElecInfo_Pt->at(i), treeReader.ElecInfo_Eta->at(i), treeReader.ElecInfo_Phi->at(i), treeReader.ElecInfo_Energy->at(i));
                double delta_R = electron.DeltaR(leading_photon);
                if( delta_R<0.3 ) continue;
                delta_R = electron.DeltaR(subleading_photon);
                if( delta_R<0.3 ) continue;
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
        //### }}}
        //### muon selection {{{
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
                if( fabs(treeReader.MuonInfo_Pt->at(i))  < 20  ) continue;
                if( treeReader.MuonInfo_PFIsoDeltaBetaCorrR04->at(i)  > 0.25  ) continue;
                //--- check deltaR(muon,photon) ---//
                TLorentzVector muon; 
                muon.SetPtEtaPhiE(treeReader.MuonInfo_Pt->at(i), treeReader.MuonInfo_Eta->at(i), treeReader.MuonInfo_Phi->at(i), treeReader.MuonInfo_Energy->at(i));
                double delta_R = muon.DeltaR(leading_photon);
                if( delta_R<0.3 ) continue;
                delta_R = muon.DeltaR(subleading_photon);
                if( delta_R<0.3 ) continue;
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
        //### }}}
        //### jet selection {{{
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
                if( fabs(treeReader.JetInfo_Pt->at(i))  < 25  ) continue;
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
                mytree.JetInfo_jet_energy.push_back(treeReader.JetInfo_Phi->at(i));
                
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
        //### }}}
        //### Check b quark info, store diphoton and jets{{{
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
        hist_num_selected_jets->Fill(mytree.num_jets);
        hist_mass_diphoton->Fill(diphoton.M());
        //}}}
        //### Leptonic event selection{{{
        if(mytree.num_leptons<1) continue;
        if(mytree.num_jets<1) continue;
        //###}}}
        //}}}

        //================================================//
        //-----------   top reconstruction     -----------//
        //================================================//
        TLorentzVector bjet;
        TLorentzVector genParticle_bjet;
        //### b-jet{{{
        //----- b-jet -----//
        int num_bjets_loose = 0, num_bjets_medium = 0;
        int num_bjets_tight = 0, index_bjet = -999;
        //TLorentzVector bjet_tight, bjet_loose, bjet_medium;
        std::vector<TLorentzVector> bjets_tight, bjets_loose, bjets_medium;
        for(int i=0; i<mytree.num_jets; ++i){
            if(mytree.JetInfo_jet_pfDeepCSVJetTags_probb.at(i)+mytree.JetInfo_jet_pfDeepCSVJetTags_probbb.at(i) >= pfDeepCSVJetTags_loose){
                num_bjets_loose += 1;
                bjets_loose.push_back(Jets[i]);
            }
            if(mytree.JetInfo_jet_pfDeepCSVJetTags_probb.at(i)+mytree.JetInfo_jet_pfDeepCSVJetTags_probbb.at(i) >= pfDeepCSVJetTags_medium){
                num_bjets_medium += 1;
                bjets_medium.push_back(Jets[i]);
            }
            if(mytree.JetInfo_jet_pfDeepCSVJetTags_probb.at(i)+mytree.JetInfo_jet_pfDeepCSVJetTags_probbb.at(i) >= pfDeepCSVJetTags_tight){
                num_bjets_tight += 1;
                bjets_tight.push_back(Jets[i]);
                index_bjet = i;
            }
        }//end of looping jets
        hist_num_bjets_loose->Fill(num_bjets_loose);
        hist_num_bjets_medium->Fill(num_bjets_medium);
        hist_num_bjets_tight->Fill(num_bjets_tight);

        //--- counting correct rate in when num_bjet==1 ---//
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

        //if(num_bjets_loose != 1) continue;//exactly 1 loose bjet candidate
        //if(num_bjets_loose < 1) continue;//at least 1 loose bjet candidate
        //bjet = bjets_loose[0];

        if(num_bjets_tight != 1) continue;//exactly 1 tight bjet candidate
        //if(num_bjets_tight < 1) continue;//at least 1 tight bjet candidate
        bjet = bjets_tight[0];
        
        //if(num_bjets_medium != 1) continue;//exactly 1 medium bjet candidate
        //if(num_bjets_medium < 1) continue;//at least 1 medium bjet candidate
        //bjet = bjets_medium[0];
        //### }}}
        //# LEPTONIC Channel{{{
        // geninfo{{{
        TLorentzVector chargedLepton, neutrino, wboson, bquark, topquark;
        for(int i=0; i<treeReader.GenPartInfo_size; i++){
            int pdgID = treeReader.GenPartInfo_PdgID->at(i);
            int mom_pdgID = treeReader.GenPartInfo_MomPdgID->at(i);
            bool isPromptFinalState = treeReader.GenPartInfo_isPromptFinalState->at(i);
            bool isNeutrino = (abs(pdgID) == 12 || abs(pdgID) == 14 || abs(pdgID) == 16) && isPromptFinalState;
            bool isChargedLepton = (abs(pdgID) == 11 || abs(pdgID) == 13 || abs(pdgID) == 15) && isPromptFinalState;
            bool isWboson = (abs(pdgID) == 24);
            bool isbquark = (abs(pdgID) == 5 && abs(mom_pdgID) == 6);
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
                topquark = wboson + bquark;
                //kinematics_info("t", topquark);
                hist_mass_gen_wboson_leptonic->Fill(wboson.M());
                hist_mass_gen_topquark_leptonic->Fill(topquark.M());
            }
        }
        //TLorentzVector wboson_check = neutrino + chargedLepton;
        //kinematics_info("W", wboson_check);
        //}}}
        //Reconstruct Mass (W, M2, M1){{{
        float met_pt = treeReader.MetInfo_Pt;
        float met_phi = treeReader.MetInfo_Phi;
        float met_px = treeReader.MetInfo_Px;
        float met_py = treeReader.MetInfo_Py;
        float met_sumET = treeReader.MetInfo_SumET;

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
        
        hist_reco_Pz_gen_Pz -> Fill(neutrino.Pz(), met_pz_solution_2);
        if(coefficient_D<0) hist_reco_Pz_gen_Pz_negativeD -> Fill(neutrino.Pz(), met_pz_solution_2);
        if(coefficient_D>0) hist_reco_Pz_gen_Pz_positiveD -> Fill(neutrino.Pz(), met_pz_solution_2);
        hist_MetInfo_Pz_solution_1 -> Fill(met_pz_solution_1);
        hist_MetInfo_Pz_solution_2 -> Fill(met_pz_solution_2);
        if(coefficient_D>0) hist_MetInfo_Pz_solution_1_positiveD -> Fill(met_pz_solution_1);
        if(coefficient_D>0) hist_MetInfo_Pz_solution_2_positiveD -> Fill(met_pz_solution_2);
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

        hist_leptonic_w_candidate_solution1_pt->Fill(L_w_lep[0].Pt());
        hist_leptonic_w_candidate_solution1_eta->Fill(L_w_lep[0].Eta());
        hist_leptonic_w_candidate_solution1_mass->Fill(L_w_lep[0].M());
        if(coefficient_D<0) hist_leptonic_w_candidate_solution1_mass_negativeD->Fill(L_w_lep[0].M());
        hist_leptonic_top_tbw_solution1_pt->Fill(L_bw_lep[0].Pt());
        hist_leptonic_top_tbw_solution1_eta->Fill(L_bw_lep[0].Eta());
        hist_leptonic_top_tbw_solution1_mass->Fill(L_bw_lep[0].M());
        if(coefficient_D<0) hist_leptonic_top_tbw_solution1_mass_negativeD->Fill(L_bw_lep[0].M());

        hist_leptonic_w_candidate_solution2_pt->Fill(L_w_lep[1].Pt());
        hist_leptonic_w_candidate_solution2_eta->Fill(L_w_lep[1].Eta());
        hist_leptonic_w_candidate_solution2_mass->Fill(L_w_lep[1].M());
        if(coefficient_D<0)hist_leptonic_w_candidate_solution2_mass_negativeD->Fill(L_w_lep[1].M());
        hist_leptonic_top_tbw_solution2_pt->Fill(L_bw_lep[1].Pt());
        hist_leptonic_top_tbw_solution2_eta->Fill(L_bw_lep[1].Eta());
        hist_leptonic_top_tbw_solution2_mass->Fill(L_bw_lep[1].M());
        if(coefficient_D<0) hist_leptonic_top_tbw_solution2_mass_negativeD->Fill(L_bw_lep[1].M());

        //## comments: further MC truth study is needed
        //if(L_w_lep[1].M()>90){
        //    printf("[check - wboson] mass = %f\n", L_w_lep[1].M());
        //    printf("[check - wboson] coefficient_D = %f\n", coefficient_D);
        //    printf("[check - wboson] lepton px = %6.2f, ", lepton.Px());
        //    printf("py = %6.2f, ", lepton.Py());
        //    printf("pz = %6.2f, ", lepton.Pz());
        //    printf("energy = %6.2f\n", lepton.E());
        //    printf("[check - wboson] metinf px = %6.2f, ", met_px);
        //    printf("py = %6.2f, ", met_py);
        //    printf("pz = %6.2f, ", met_pz_solution_2);
        //    printf("energy = %6.2f\n", met_energy_solution_2);
        //    printf("[check - wboson] wboson px = %6.2f, ", L_w_lep[1].Px());
        //    printf("py = %6.2f, ", L_w_lep[1].Py());
        //    printf("pz = %6.2f, ", L_w_lep[1].Pz());
        //    printf("energy = %6.2f\n", L_w_lep[1].E());
        //}
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
        //--- deltaR/Pz/PT(gen, reco){{{
        double deltaR;
        deltaR = diphoton.DeltaR(L_bw_lep[1])                             ; hist_deltaR_reco_top_higgs_leptonic->Fill(deltaR)   ;
        deltaR = chargedLepton.DeltaR(lepton)                             ; hist_deltaR_gen_reco_chargedLepton->Fill(deltaR)   ;
        deltaR = neutrino.DeltaR(L_met_lep[0])                            ; hist_deltaR_gen_reco_neutrino_sol0->Fill(deltaR)   ;
        deltaR = neutrino.DeltaR(L_met_lep[1])                            ; hist_deltaR_gen_reco_neutrino_sol1->Fill(deltaR)   ;
        deltaR = wboson.DeltaR(L_w_lep[0])                                ; hist_deltaR_gen_reco_wboson_sol0->Fill(deltaR)     ;
        deltaR = wboson.DeltaR(L_w_lep[1])                                ; hist_deltaR_gen_reco_wboson_sol1->Fill(deltaR)     ;
        double deltaPz;
        deltaPz = (lepton.Pz() - chargedLepton.Pz()) ; hist_deltaPz_gen_reco_chargedLepton->Fill(deltaPz) ;
        deltaPz = (L_met_lep[0].Pz() - neutrino.Pz()) ; hist_deltaPz_gen_reco_neutrino_sol0->Fill(deltaPz) ;
        deltaPz = (L_met_lep[1].Pz() - neutrino.Pz()) ; hist_deltaPz_gen_reco_neutrino_sol1->Fill(deltaPz) ;
        deltaPz = (L_w_lep[0].Pz() - wboson.Pz()) ; hist_deltaPz_gen_reco_wboson_sol0->Fill(deltaPz)   ;
        deltaPz = (L_w_lep[1].Pz() - wboson.Pz()) ; hist_deltaPz_gen_reco_wboson_sol1->Fill(deltaPz)   ;
        //deltaPz = (lepton.Pz() - chargedLepton.Pz()) / chargedLepton.Pz() ; hist_deltaPz_gen_reco_chargedLepton->Fill(deltaPz) ;
        deltaPz = (L_met_lep[0].Pz() - neutrino.Pz()) / neutrino.Pz()     ; hist_deltaPzRatio_gen_reco_neutrino_sol0->Fill(deltaPz) ;
        deltaPz = (L_met_lep[1].Pz() - neutrino.Pz()) / neutrino.Pz()     ; hist_deltaPzRatio_gen_reco_neutrino_sol1->Fill(deltaPz) ;
        if(coefficient_D>0){
            deltaPz = (L_met_lep[0].Pz() - neutrino.Pz()) / neutrino.Pz()     ; hist_deltaPzRatio_gen_reco_neutrino_sol0_positiveD->Fill(deltaPz) ;
            deltaPz = (L_met_lep[1].Pz() - neutrino.Pz()) / neutrino.Pz()     ; hist_deltaPzRatio_gen_reco_neutrino_sol1_positiveD->Fill(deltaPz) ;
        }
        if(coefficient_D<0){
            deltaPz = (L_met_lep[0].Pz() - neutrino.Pz()) / neutrino.Pz()     ; hist_deltaPzRatio_gen_reco_neutrino_sol0_negativeD->Fill(deltaPz) ;
            deltaPz = (L_met_lep[1].Pz() - neutrino.Pz()) / neutrino.Pz()     ; hist_deltaPzRatio_gen_reco_neutrino_sol1_negativeD->Fill(deltaPz) ;
        }

        deltaPz = (L_w_lep[0].Pz() - wboson.Pz()) / wboson.Pz()           ; hist_deltaPzRatio_gen_reco_wboson_sol0->Fill(deltaPz)   ;
        deltaPz = (L_w_lep[1].Pz() - wboson.Pz()) / wboson.Pz()           ; hist_deltaPzRatio_gen_reco_wboson_sol1->Fill(deltaPz)   ;
        double deltaPT;
        deltaPT = (lepton.Pt() - chargedLepton.Pt()) / chargedLepton.Pt() ; hist_deltaPT_gen_reco_chargedLepton->Fill(deltaPT) ;
        deltaPT = (L_met_lep[0].Pt() - neutrino.Pt()) / neutrino.Pt()     ; hist_deltaPT_gen_reco_neutrino_sol0->Fill(deltaPT) ;
        deltaPT = (L_met_lep[1].Pt() - neutrino.Pt()) / neutrino.Pt()     ; hist_deltaPT_gen_reco_neutrino_sol1->Fill(deltaPT) ;
        deltaPT = (L_w_lep[0].Pt() - wboson.Pt()) / wboson.Pt()           ; hist_deltaPT_gen_reco_wboson_sol0->Fill(deltaPT)   ;
        deltaPT = (L_w_lep[1].Pt() - wboson.Pt()) / wboson.Pt()           ; hist_deltaPT_gen_reco_wboson_sol1->Fill(deltaPT)   ;
        deltaPT = (met_pt - L_met_lep[0].Pt()) / met_pt                   ; hist_deltaPT_reco_neutrino_sol0_doubleCheck->Fill(deltaPT) ;
        deltaPT = (met_pt - L_met_lep[1].Pt()) / met_pt                   ; hist_deltaPT_reco_neutrino_sol1_doubleCheck->Fill(deltaPT) ;
        //}}}
//        //--- M1{{{
//        double M1;
//        int index_q;
//        TLorentzVector jet_q;
//        TLorentzVector top_fcnh = GetBestM1(M1, treeReader.num_jets, index_bjet, index_jet_chi2_modified, diphoton, Jets, index_q, jet_q);
//        if(M1 != -999){
//            //-----
//            h[hist_top_tqh_pt]->Fill(top_fcnh.Pt());//exclude event without suitable candidate
//            h[hist_top_tqh_eta]->Fill(top_fcnh.Eta());//exclude event without suitable candidate
//            h[hist_top_tqh_mass]->Fill(top_fcnh.M());//exclude event without suitable candidate
//            //-----
//            h[hist_jetq_pt] -> Fill(jet_q.Pt());
//            h[hist_jetq_eta] -> Fill(jet_q.Eta());
//            h[hist_jetq_phi] -> Fill(jet_q.Phi());
//            h[hist_jetq_energy] -> Fill(jet_q.E());
//            h[hist_jetq_btag_score] -> Fill(Jets_btag_score[index_q]);
//            h[hist_jetq_diphoton_deltaR] -> Fill(jet_q.DeltaR(diphoton));
//            if(treeReader.num_leptons>0){
//                for(int i=0; i<treeReader.num_leptons; ++i){
//                    double delta_R = jet_q.DeltaR(Leptons.at(i));
//                    h[hist_jetq_lepton_deltaR] -> Fill(delta_R);
//                }
//            }
//        }
//
//        //}}}
//        //deltaR{{{
//        //deltaR 
//        double deltaR;
//        if(M1 != -999){ deltaR = diphoton.DeltaR(jet_q)                       ; h[hist_deltaR_qH] -> Fill(deltaR)            ; }
//        deltaR = leading_photon.DeltaR(subleading_photon)                     ; h[hist_deltaR_photon_photon] -> Fill(deltaR) ;
//        //if(M1 != -999){ deltaR = top_candidate_chi2_modified.DeltaR(top_fcnh) ; h[hist_deltaR_top_top] -> Fill(deltaR)       ; }
//        //deltaR = bjet.DeltaR(w_candidate_chi2_modified)                       ; h[hist_deltaR_bW] -> Fill(deltaR)            ;
//        //deltaR = diphoton.DeltaR(w_candidate_chi2_modified)                   ; h[hist_deltaR_HW] -> Fill(deltaR)            ;
//        //deltaR = diphoton.DeltaR(top_candidate_chi2_modified)                 ; h[hist_deltaR_tH] -> Fill(deltaR)            ;
//
//        //deltaR = Jets[0].DeltaR(Jets[1])                                      ; h[hist_deltaR_jet1_jet2] -> Fill(deltaR)     ;
//        //}}}
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
        kf->Run(); // Run the tool
        int NPerm = kf->GetNPerm(); // Get number of permutations
        std::vector<float> NuPz, disc;

        if(NPerm!=1) printf("[WARNNING] NPerm = %d\n", NPerm);

        for(int ip=0;ip<NPerm;ip++) // Loop over permutations - already sorted in likelihood value from min to max
        {
            //printf("[INFO-kinfit] disc = %f\n", disc);
            disc.push_back( kf->GetDisc(ip) ); // Get minimized likelihood value
                       
            // Get reconstructed neutrino
            float NuPx = kf->GetNuPx(ip,0);
            float NuPy = kf->GetNuPy(ip,0);
            NuPz.push_back( kf->GetNuPz(ip,0) );
                                                                           
            //// Get index of object in input collection
            //int idxTopLepWElecFit = kf->GetIndex(ELECTRON_TOPLEP,ip);
            //int idxTopLepWMuonFit = kf->GetIndex(MUON_TOPLEP,ip);
            //int idxTopLepBJetFit = kf->GetIndex(BJETLEP_TOPLEP,ip);
            //int idxTopHadNonBJetFit = kf->GetIndex(NONBJETHAD_TOPLEP,ip);
            //int idxHiggsBJet1Fit = kf->GetIndex(BJET1_TOPLEP,ip);
            //int idxHiggsBJet2Fit = kf->GetIndex(BJET2_TOPLEP,ip);
                                                       
            //// Build up b jet from leptonic top quark decay
            //float TopLepBJetFitPt = BJetPt[idxTopLepBJetFit];
            //float TopLepBJetFitEta = BJetEta[idxTopLepBJetFit];
            //float TopLepBJetFitPhi = BJetPhi[idxTopLepBJetFit];
            //float TopLepBJetFitE = BJetE[idxTopLepBJetFit];
        }

        if(disc[0] > 100000.) counter_irregular_disc += 1;

        //2D plots
        hist_reco_Pz_gen_Pz_topKinFit -> Fill(neutrino.Pz(), NuPz[0]);
        if(disc[0] > 100000.) hist_reco_Pz_gen_Pz_topKinFit_irrDisc -> Fill(neutrino.Pz(), NuPz[0]);
        else                  hist_reco_Pz_gen_Pz_topKinFit_regDisc -> Fill(neutrino.Pz(), NuPz[0]);

        //deltaPz/Pz
        deltaPz = (NuPz[0]-neutrino.Pz()) / neutrino.Pz();
        hist_deltaPzRatio_gen_reco_neutrino_topKinFit->Fill(deltaPz);
        if(disc[0] > 100000.) hist_deltaPzRatio_gen_reco_neutrino_topKinFit_irrDisc -> Fill(deltaPz);
        else                  hist_deltaPzRatio_gen_reco_neutrino_topKinFit_regDisc -> Fill(deltaPz);

        //}}}

        hist_gen_neutrino_pz->Fill(neutrino.Pz());
        //}}}
    }// End of event loop.
    //==================================================//
    //---------------------  Report  -------------------//
    //==================================================//
    //# LEPTONIC reco study{{{
    printf("[INFO] correlation factor all = %f\n", hist_reco_Pz_gen_Pz->GetCorrelationFactor());
    printf("[INFO] correlation factor all (topKinFit) = %f\n", hist_reco_Pz_gen_Pz_topKinFit->GetCorrelationFactor());
    printf("-----\n");
    printf("[INFO] correlation factor positive D = %f\n", hist_reco_Pz_gen_Pz_positiveD->GetCorrelationFactor());
    printf("[INFO] correlation factor negative D = %f\n", hist_reco_Pz_gen_Pz_negativeD->GetCorrelationFactor());
    printf("[INFO] correlation factor regular Disc   = %f\n", hist_reco_Pz_gen_Pz_topKinFit_regDisc->GetCorrelationFactor());
    printf("[INFO] correlation factor irregular Disc = %f\n", hist_reco_Pz_gen_Pz_topKinFit_irrDisc->GetCorrelationFactor());
    printf("-----\n");
    printf("[CHECK] number of irregular disc value = %d / %d\n", counter_irregular_disc, counter_coeff_D);
    printf("[CHECK] number of counter_coeff_D_isNegative = %d / %d\n", counter_coeff_D_isNegative, counter_coeff_D);
    printf("[CHECK] number of counter_coeff_D_isNegative_gen = %d / %d\n", counter_coeff_D_isNegative_gen, counter_coeff_D_gen);

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    // 1-D plots{{{
    const char* dir = "ntuples_skimmed/chi2_study_leptonic_1D_plots";
    hist_mass_gen_wboson_leptonic -> Draw("hist")               ; c1->SaveAs( Form("%s/hist_mass_gen_wboson_leptonic.png", dir) )               ;
    hist_mass_gen_topquark_leptonic -> Draw("hist")             ; c1->SaveAs( Form("%s/hist_mass_gen_topquark_leptonic.png", dir) )             ;
    hist_deltaR_reco_top_higgs_leptonic -> Draw("hist")         ; c1->SaveAs( Form("%s/hist_deltaR_reco_top_higgs_leptonic.png", dir) )         ;
    hist_deltaR_gen_reco_chargedLepton -> Draw("hist")          ; c1->SaveAs( Form("%s/hist_deltaR_gen_reco_chargedLepton.png", dir) )          ;
    hist_deltaR_gen_reco_neutrino_sol0 -> Draw("hist")          ; c1->SaveAs( Form("%s/hist_deltaR_gen_reco_neutrino_sol0.png", dir) )          ;
    hist_deltaR_gen_reco_wboson_sol0 -> Draw("hist")            ; c1->SaveAs( Form("%s/hist_deltaR_gen_reco_wboson_sol0.png", dir) )            ;
    hist_deltaR_gen_reco_neutrino_sol1 -> Draw("hist")          ; c1->SaveAs( Form("%s/hist_deltaR_gen_reco_neutrino_sol1.png", dir) )          ;
    hist_deltaR_gen_reco_wboson_sol1 -> Draw("hist")            ; c1->SaveAs( Form("%s/hist_deltaR_gen_reco_wboson_sol1.png", dir) )            ;
    //---
    hist_deltaPz_gen_reco_chargedLepton -> Draw("hist")         ; c1->SaveAs( Form("%s/hist_deltaPz_gen_reco_chargedLepton.png", dir) )         ;
    hist_deltaPz_gen_reco_neutrino_sol0 -> Draw("hist")         ; c1->SaveAs( Form("%s/hist_deltaPz_gen_reco_neutrino_sol0.png", dir) )         ;
    hist_deltaPz_gen_reco_wboson_sol0 -> Draw("hist")           ; c1->SaveAs( Form("%s/hist_deltaPz_gen_reco_wboson_sol0.png", dir) )           ;
    hist_deltaPz_gen_reco_neutrino_sol1 -> Draw("hist")         ; c1->SaveAs( Form("%s/hist_deltaPz_gen_reco_neutrino_sol1.png", dir) )         ;
    hist_deltaPz_gen_reco_wboson_sol1 -> Draw("hist")           ; c1->SaveAs( Form("%s/hist_deltaPz_gen_reco_wboson_sol1.png", dir) )           ;
    //---
    hist_deltaPzRatio_gen_reco_neutrino_sol0 -> Draw("hist")    ; c1->SaveAs( Form("%s/hist_deltaPzRatio_gen_reco_neutrino_sol0.png", dir) )    ;
    hist_deltaPzRatio_gen_reco_wboson_sol0 -> Draw("hist")      ; c1->SaveAs( Form("%s/hist_deltaPzRatio_gen_reco_wboson_sol0.png", dir) )      ;
    hist_deltaPzRatio_gen_reco_neutrino_sol1 -> Draw("hist")    ; c1->SaveAs( Form("%s/hist_deltaPzRatio_gen_reco_neutrino_sol1.png", dir) )    ;
    hist_deltaPzRatio_gen_reco_wboson_sol1 -> Draw("hist")      ; c1->SaveAs( Form("%s/hist_deltaPzRatio_gen_reco_wboson_sol1.png", dir) )      ;
    //---
    hist_deltaPT_gen_reco_chargedLepton -> Draw("hist")         ; c1->SaveAs( Form("%s/hist_deltaPT_gen_reco_chargedLepton.png", dir) )         ;
    hist_deltaPT_gen_reco_neutrino_sol0 -> Draw("hist")         ; c1->SaveAs( Form("%s/hist_deltaPT_gen_reco_neutrino_sol0.png", dir) )         ;
    hist_deltaPT_gen_reco_wboson_sol0 -> Draw("hist")           ; c1->SaveAs( Form("%s/hist_deltaPT_gen_reco_wboson_sol0.png", dir) )           ;
    hist_deltaPT_gen_reco_neutrino_sol1 -> Draw("hist")         ; c1->SaveAs( Form("%s/hist_deltaPT_gen_reco_neutrino_sol1.png", dir) )         ;
    hist_deltaPT_gen_reco_wboson_sol1 -> Draw("hist")           ; c1->SaveAs( Form("%s/hist_deltaPT_gen_reco_wboson_sol1.png", dir) )           ;
    //---
    hist_deltaPT_reco_neutrino_sol0_doubleCheck -> Draw("hist") ; c1->SaveAs( Form("%s/hist_deltaPT_reco_neutrino_sol0_doubleCheck.png", dir) ) ;
    hist_deltaPT_reco_neutrino_sol1_doubleCheck -> Draw("hist") ; c1->SaveAs( Form("%s/hist_deltaPT_reco_neutrino_sol1_doubleCheck.png", dir) ) ;
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
    hist_leptonic_w_candidate_solution1_pt -> Draw("hist")      ; c1->SaveAs( Form("%s/hist_leptonic_w_candidate_solution1_pt.png", dir) )      ;
    hist_leptonic_w_candidate_solution1_eta -> Draw("hist")     ; c1->SaveAs( Form("%s/hist_leptonic_w_candidate_solution1_eta.png", dir) )     ;
    hist_leptonic_w_candidate_solution1_mass -> Draw("hist")    ; c1->SaveAs( Form("%s/hist_leptonic_w_candidate_solution1_mass.png", dir) )    ;
    hist_leptonic_top_tbw_solution1_pt -> Draw("hist")          ; c1->SaveAs( Form("%s/hist_leptonic_top_tbw_solution1_pt.png", dir) )          ;
    hist_leptonic_top_tbw_solution1_eta -> Draw("hist")         ; c1->SaveAs( Form("%s/hist_leptonic_top_tbw_solution1_eta.png", dir) )         ;
    hist_leptonic_top_tbw_solution1_mass -> Draw("hist")        ; c1->SaveAs( Form("%s/hist_leptonic_top_tbw_solution1_mass.png", dir) )        ;
    //---
    hist_leptonic_w_candidate_solution2_pt -> Draw("hist")      ; c1->SaveAs( Form("%s/hist_leptonic_w_candidate_solution2_pt.png", dir) )      ;
    hist_leptonic_w_candidate_solution2_eta -> Draw("hist")     ; c1->SaveAs( Form("%s/hist_leptonic_w_candidate_solution2_eta.png", dir) )     ;
    hist_leptonic_w_candidate_solution2_mass -> Draw("hist")    ; c1->SaveAs( Form("%s/hist_leptonic_w_candidate_solution2_mass.png", dir) )    ;
    hist_leptonic_top_tbw_solution2_pt -> Draw("hist")          ; c1->SaveAs( Form("%s/hist_leptonic_top_tbw_solution2_pt.png", dir) )          ;
    hist_leptonic_top_tbw_solution2_eta -> Draw("hist")         ; c1->SaveAs( Form("%s/hist_leptonic_top_tbw_solution2_eta.png", dir) )         ;
    hist_leptonic_top_tbw_solution2_mass -> Draw("hist")        ; c1->SaveAs( Form("%s/hist_leptonic_top_tbw_solution2_mass.png", dir) )        ;
    //}}}
    //1-D plots with legend{{{
    TLegend *legend = new TLegend(0.60,0.65,0.85,0.85);
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

    //# highlight the portion of events with D < 0{{{
    //void makeComparisonPlots_2(TH1D* hist, TH1D* hist_a, TLegend *legend, TCanvas *c1, ){} // # to be finished later
    hist_leptonic_w_candidate_solution1_mass -> SetStats(0)    ;
    hist_leptonic_w_candidate_solution1_mass -> SetLineWidth(2)    ;
    hist_leptonic_w_candidate_solution1_mass -> Draw("hist")    ;
    hist_leptonic_w_candidate_solution1_mass_negativeD -> Draw("hist;same")    ;
    hist_leptonic_w_candidate_solution1_mass_negativeD -> SetLineWidth(0)    ;
    hist_leptonic_w_candidate_solution1_mass_negativeD -> SetFillStyle(3001)    ;
    hist_leptonic_w_candidate_solution1_mass_negativeD -> SetFillColor(kRed)    ;
    legend->Clear();
    legend->AddEntry(hist_leptonic_w_candidate_solution1_mass,  "Every reco. W boson", "l");
    legend->AddEntry(hist_leptonic_w_candidate_solution1_mass_negativeD,  "D < 0", "f");
    legend->SetLineColor(0);
    legend->Draw("same");
    c1->SaveAs("ntuples_skimmed/hist_leptonic_w_candidate_solution1_mass_negativeD.png")    ;

    hist_leptonic_w_candidate_solution2_mass -> SetStats(0)    ;
    hist_leptonic_w_candidate_solution2_mass -> SetLineWidth(2)    ;
    hist_leptonic_w_candidate_solution2_mass -> Draw("hist")    ;
    hist_leptonic_w_candidate_solution2_mass_negativeD -> Draw("hist;same")    ;
    hist_leptonic_w_candidate_solution2_mass_negativeD -> SetLineWidth(0)    ;
    hist_leptonic_w_candidate_solution2_mass_negativeD -> SetFillStyle(3001)    ;
    hist_leptonic_w_candidate_solution2_mass_negativeD -> SetFillColor(kRed)    ;
    legend->Clear();
    legend->AddEntry(hist_leptonic_w_candidate_solution2_mass,  "Every reco. W boson", "l");
    legend->AddEntry(hist_leptonic_w_candidate_solution2_mass_negativeD,  "D < 0", "f");
    legend->SetLineColor(0);
    legend->Draw("same");
    c1->SaveAs("ntuples_skimmed/hist_leptonic_w_candidate_solution2_mass_negativeD.png")    ;

    //---
    hist_leptonic_top_tbw_solution1_mass -> SetStats(0)    ;
    hist_leptonic_top_tbw_solution1_mass -> SetLineWidth(2)    ;
    hist_leptonic_top_tbw_solution1_mass -> Draw("hist")    ;
    hist_leptonic_top_tbw_solution1_mass_negativeD -> Draw("hist;same")    ;
    hist_leptonic_top_tbw_solution1_mass_negativeD -> SetLineWidth(0)    ;
    hist_leptonic_top_tbw_solution1_mass_negativeD -> SetFillStyle(3001)    ;
    hist_leptonic_top_tbw_solution1_mass_negativeD -> SetFillColor(kRed)    ;
    legend->Clear();
    legend->AddEntry(hist_leptonic_top_tbw_solution1_mass,  "Every reco. top quark", "l");
    legend->AddEntry(hist_leptonic_top_tbw_solution1_mass_negativeD,  "D < 0", "f");
    legend->SetLineColor(0);
    legend->Draw("same");
    c1->SaveAs("ntuples_skimmed/hist_leptonic_top_tbw_solution1_mass_negativeD.png")    ;

    hist_leptonic_top_tbw_solution2_mass -> SetStats(0)    ;
    hist_leptonic_top_tbw_solution2_mass -> SetLineWidth(2)    ;
    hist_leptonic_top_tbw_solution2_mass -> Draw("hist")    ;
    hist_leptonic_top_tbw_solution2_mass_negativeD -> Draw("hist;same")    ;
    hist_leptonic_top_tbw_solution2_mass_negativeD -> SetLineWidth(0)    ;
    hist_leptonic_top_tbw_solution2_mass_negativeD -> SetFillStyle(3001)    ;
    hist_leptonic_top_tbw_solution2_mass_negativeD -> SetFillColor(kRed)    ;
    legend->Clear();
    legend->AddEntry(hist_leptonic_top_tbw_solution2_mass,  "Every reco. top quark", "l");
    legend->AddEntry(hist_leptonic_top_tbw_solution2_mass_negativeD,  "D < 0", "f");
    legend->SetLineColor(0);
    legend->Draw("same");
    c1->SaveAs("ntuples_skimmed/hist_leptonic_top_tbw_solution2_mass_negativeD.png")    ;
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
    //---
    TLegend *legend_ratio = new TLegend(0.65,0.53,0.85,0.73);
    hist_deltaPzRatio_gen_reco_neutrino_topKinFit -> Draw("hist")                  ;
    hist_deltaPzRatio_gen_reco_neutrino_topKinFit -> SetLineWidth(2)                  ;
    hist_deltaPzRatio_gen_reco_neutrino_topKinFit_regDisc -> Draw("hist;same")                  ;
    hist_deltaPzRatio_gen_reco_neutrino_topKinFit_regDisc -> SetLineColor(kRed)                  ;
    hist_deltaPzRatio_gen_reco_neutrino_topKinFit_regDisc -> SetLineWidth(2)                  ;
    hist_deltaPzRatio_gen_reco_neutrino_topKinFit_regDisc -> SetLineStyle(2)                  ;
    hist_deltaPzRatio_gen_reco_neutrino_topKinFit_irrDisc -> Draw("hist;same")                  ;
    hist_deltaPzRatio_gen_reco_neutrino_topKinFit_irrDisc -> SetLineColor(kGreen)                  ;
    hist_deltaPzRatio_gen_reco_neutrino_topKinFit_irrDisc -> SetLineWidth(2)                  ;
    legend_ratio->Clear();
    legend_ratio->AddEntry(hist_deltaPzRatio_gen_reco_neutrino_topKinFit, "All", "l");
    legend_ratio->AddEntry(hist_deltaPzRatio_gen_reco_neutrino_topKinFit_regDisc, "regular Disc.", "l");
    legend_ratio->AddEntry(hist_deltaPzRatio_gen_reco_neutrino_topKinFit_irrDisc, "irregular Disc.", "l");
    legend_ratio->SetLineColor(0);
    legend_ratio->Draw("same");
    c1->SaveAs("ntuples_skimmed/hist_deltaPzRatio_gen_reco_neutrino_topKinFit.png")                  ;
    //---
    hist_deltaPzRatio_gen_reco_neutrino_sol0 -> Draw("hist")                  ;
    hist_deltaPzRatio_gen_reco_neutrino_sol0 -> SetLineWidth(2)                  ;
    hist_deltaPzRatio_gen_reco_neutrino_sol0_positiveD -> Draw("hist;same")                  ;
    hist_deltaPzRatio_gen_reco_neutrino_sol0_positiveD -> SetLineColor(kRed)                  ;
    hist_deltaPzRatio_gen_reco_neutrino_sol0_positiveD -> SetLineWidth(2)                  ;
    hist_deltaPzRatio_gen_reco_neutrino_sol0_positiveD -> SetLineStyle(2)                  ;
    hist_deltaPzRatio_gen_reco_neutrino_sol0_negativeD -> Draw("hist;same")                  ;
    hist_deltaPzRatio_gen_reco_neutrino_sol0_negativeD -> SetLineColor(kGreen)                  ;
    hist_deltaPzRatio_gen_reco_neutrino_sol0_negativeD -> SetLineWidth(2)                  ;
    legend_ratio->Clear();
    legend_ratio->AddEntry(hist_deltaPzRatio_gen_reco_neutrino_sol0, "All", "l");
    legend_ratio->AddEntry(hist_deltaPzRatio_gen_reco_neutrino_sol0_positiveD, "D > 0", "l");
    legend_ratio->AddEntry(hist_deltaPzRatio_gen_reco_neutrino_sol0_negativeD, "D < 0", "l");
    legend_ratio->SetLineColor(0);
    legend_ratio->Draw("same");
    c1->SaveAs("ntuples_skimmed/hist_deltaPzRatio_gen_reco_neutrino_sol0.png")                  ;
    //---
    hist_deltaPzRatio_gen_reco_neutrino_sol1 -> Draw("hist")                  ;
    hist_deltaPzRatio_gen_reco_neutrino_sol1 -> SetLineWidth(2)                  ;
    hist_deltaPzRatio_gen_reco_neutrino_sol1_positiveD -> Draw("hist;same")                  ;
    hist_deltaPzRatio_gen_reco_neutrino_sol1_positiveD -> SetLineColor(kRed)                  ;
    hist_deltaPzRatio_gen_reco_neutrino_sol1_positiveD -> SetLineWidth(2)                  ;
    hist_deltaPzRatio_gen_reco_neutrino_sol1_positiveD -> SetLineStyle(2)                  ;
    hist_deltaPzRatio_gen_reco_neutrino_sol1_negativeD -> Draw("hist;same")                  ;
    hist_deltaPzRatio_gen_reco_neutrino_sol1_negativeD -> SetLineColor(kGreen)                  ;
    hist_deltaPzRatio_gen_reco_neutrino_sol1_negativeD -> SetLineWidth(2)                  ;
    legend_ratio->Clear();
    legend_ratio->AddEntry(hist_deltaPzRatio_gen_reco_neutrino_sol1, "All", "l");
    legend_ratio->AddEntry(hist_deltaPzRatio_gen_reco_neutrino_sol1_positiveD, "D > 0", "l");
    legend_ratio->AddEntry(hist_deltaPzRatio_gen_reco_neutrino_sol1_negativeD, "D < 0", "l");
    legend_ratio->SetLineColor(0);
    legend_ratio->Draw("same");
    c1->SaveAs("ntuples_skimmed/hist_deltaPzRatio_gen_reco_neutrino_sol1.png")                  ;
    //}}}
    //2-D plots{{{
    hist_reco_Pz_gen_Pz_topKinFit -> SetStats(0);
    hist_reco_Pz_gen_Pz_topKinFit_regDisc -> SetStats(0);
    hist_reco_Pz_gen_Pz_topKinFit_irrDisc -> SetStats(0);
    hist_reco_Pz_gen_Pz_topKinFit -> Draw("COLZ"); c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_topKinFit.png")           ;
    hist_reco_Pz_gen_Pz_topKinFit_regDisc -> Draw("COLZ"); c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_topKinFit_regDisc.png")           ;
    hist_reco_Pz_gen_Pz_topKinFit_irrDisc -> Draw("COLZ"); c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_topKinFit_irrDisc.png")           ;
    //---
    hist_reco_Pz_gen_Pz -> SetStats(0);
    hist_reco_Pz_gen_Pz_negativeD -> SetStats(0);
    hist_reco_Pz_gen_Pz_positiveD -> SetStats(0);
    //hist_reco_Pz_gen_Pz -> SetMaximum(16);
    //hist_reco_Pz_gen_Pz_negativeD -> SetMaximum(16);
    //hist_reco_Pz_gen_Pz_positiveD -> SetMaximum(16);
    //---
    hist_reco_Pz_gen_Pz -> Draw("COLZ")           ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz.png")           ;
    hist_reco_Pz_gen_Pz_negativeD -> Draw("COLZ") ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_negativeD.png") ;
    hist_reco_Pz_gen_Pz_positiveD -> Draw("COLZ") ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_positiveD.png") ;
    gStyle->SetPalette(kRainBow);
    hist_reco_Pz_gen_Pz -> Draw("COLZ")           ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_kRainBow.png")           ;
    hist_reco_Pz_gen_Pz_negativeD -> Draw("COLZ") ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_negativeD_kRainBow.png") ;
    hist_reco_Pz_gen_Pz_positiveD -> Draw("COLZ") ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_positiveD_kRainBow.png") ;
    gStyle->SetPalette(kVisibleSpectrum);
    hist_reco_Pz_gen_Pz -> Draw("COLZ")           ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_kVisibleSpectrum.png")           ;
    hist_reco_Pz_gen_Pz_negativeD -> Draw("COLZ") ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_negativeD_kVisibleSpectrum.png") ;
    hist_reco_Pz_gen_Pz_positiveD -> Draw("COLZ") ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_positiveD_kVisibleSpectrum.png") ;
    gStyle->SetPalette(kInvertedDarkBodyRadiator);
    hist_reco_Pz_gen_Pz -> Draw("COLZ")           ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_kInvertedDarkBodyRadiator.png")           ;
    hist_reco_Pz_gen_Pz_negativeD -> Draw("COLZ") ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_negativeD_kInvertedDarkBodyRadiator.png") ;
    hist_reco_Pz_gen_Pz_positiveD -> Draw("COLZ") ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_positiveD_kInvertedDarkBodyRadiator.png") ;
    gStyle->SetPalette(kGistEarth);
    hist_reco_Pz_gen_Pz -> Draw("COLZ")           ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_kGistEarth.png")           ;
    hist_reco_Pz_gen_Pz_negativeD -> Draw("COLZ") ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_negativeD_kGistEarth.png") ;
    hist_reco_Pz_gen_Pz_positiveD -> Draw("COLZ") ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_positiveD_kGistEarth.png") ;
    gStyle->SetPalette(kCherry)                   ;
    hist_reco_Pz_gen_Pz -> Draw("COLZ")           ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_kCherry.png")              ;
    hist_reco_Pz_gen_Pz_negativeD -> Draw("COLZ") ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_negativeD_kCherry.png") ;
    hist_reco_Pz_gen_Pz_positiveD -> Draw("COLZ") ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_positiveD_kCherry.png") ;
    gStyle->SetPalette(kViridis)                  ;
    hist_reco_Pz_gen_Pz -> Draw("COLZ")           ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_kViridis.png")             ;
    hist_reco_Pz_gen_Pz_negativeD -> Draw("COLZ") ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_negativeD_kViridis.png") ;
    hist_reco_Pz_gen_Pz_positiveD -> Draw("COLZ") ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_positiveD_kViridis.png") ;
    //gStyle->SetPalette(kOcean);
    //TColor::InvertPalette();
    //TColor::InvertPalette();
    //hist_reco_Pz_gen_Pz -> Draw("COLZ")                  ; c1->SaveAs("ntuples_skimmed/hist_reco_Pz_gen_Pz_kOcean_inverted.png")                  ;
    //}}}
    //}}}
    //---
    printf("End!!");
    fout->Write();
    fout->Close();
    return 1;
}

// functions{{{
double GetM1M2_ratio(double M1, double M2){
    double ratio = abs(M1/M2-1) + abs(M2/M1-1);
    return ratio;
}
// ### GetBestM1{{{
TLorentzVector GetBestM1(double &M1, int num_jets, int index_bjet, std::vector<int> index_jet, TLorentzVector diphoton, std::vector<TLorentzVector> Jets, int &index_q, TLorentzVector &jet_q){
    //record all the possible combinations
    std::vector<int> indices;
    std::vector<double> top_fcnh_chi2;
    std::vector<TLorentzVector> top_fcnh_candidates;
    for(int i=0; i<num_jets; ++i){
        if(i==index_bjet || i==index_jet[0] || i==index_jet[1]) continue;//skip the jets for bjj
        TLorentzVector top_fcnh_tmp = diphoton + Jets[i];
        double chi2 = (top_fcnh_tmp.M() - top_quark_mass) * (top_fcnh_tmp.M() - top_quark_mass);
        indices.push_back(i);
        top_fcnh_chi2.push_back(chi2);
        top_fcnh_candidates.push_back(top_fcnh_tmp);
    }
    //choose the candidate with Mass closest to top mass
    int index_M1_min_chi2 = -999;
    double chi2_min_M1 = 40000;//200^2 > 178^2 ~ (x-top_mass)^2 //NOTE: will bound the range of M1!!!
    //pick the wanted one (with index)
    for(int i=0; i<top_fcnh_candidates.size(); ++i){ if(top_fcnh_chi2[i]<chi2_min_M1){ index_M1_min_chi2 = i; chi2_min_M1 = top_fcnh_chi2[i];} }
    TLorentzVector null;
    if(index_M1_min_chi2 == -999){ M1 = -999; return null;}//No suitable candidate
    else{
        index_q = indices[index_M1_min_chi2];
        jet_q = Jets[index_M1_min_chi2];
        M1 = top_fcnh_candidates[index_M1_min_chi2].M();
        return top_fcnh_candidates[index_M1_min_chi2];
    }
}
//}}}
// ### Alternative Mother info{{{
bool isMatched_with_Gen_tbw(std::vector<int> *GenPartInfo_MomPdgID, int index_bjet, int index_jet1, int index_jet2){
    bool isCorrect = \
    index_bjet > 0 &&\
    index_jet1 > 0 &&\
    index_jet2 > 0 &&\
    abs(GenPartInfo_MomPdgID->at(index_bjet)) == 6 &&\
    abs(GenPartInfo_MomPdgID->at(index_jet1)) == 24 &&\
    abs(GenPartInfo_MomPdgID->at(index_jet2)) == 24 &&\
    (index_jet1 != index_jet2);

    return isCorrect;
}
//old code{{{
//bool isMatched_with_Gen_W_Boson(TLorentzVector gen_w_sel, TH1D *&hist, Int_t GenPartInfo_size, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_PdgID){
//    double delta_R_gen_w = 999;
//    for(int j=0; j<GenPartInfo_size; j++){
//        if( abs(GenPartInfo_PdgID->at(j)) == 24 ){
//            TLorentzVector genParticle;
//            genParticle.SetPtEtaPhiM(GenPartInfo_Pt->at(j), GenPartInfo_Eta->at(j), GenPartInfo_Phi->at(j), GenPartInfo_Mass->at(j));
//            delta_R_gen_w = gen_w_sel.DeltaR(genParticle);
//        }
//    }//end of gen loop
//    //--- remove event with bad combination ---//
//    //if(gen_w_sel.M() < 20) continue; 
//    hist->Fill(delta_R_gen_w);
//    if(delta_R_gen_w > 0.005 || gen_w_sel.M() < 20) return false;
//    else return true;
//}
//}}}
//}}}
//### genMatching{{{
bool is_this_tqh_quark(TLorentzVector jet, Int_t GenPartInfo_size, std::vector<int> *GenPartInfo_MomPdgID, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_Status, std::vector<int> *GenPartInfo_PdgID){
    TLorentzVector truelove;//we are looking for the right genParticle to match the jet.
    int index = -999, truelove_PdgID = -999; double delta_R_min = 999.;
    for(int i=0; i<GenPartInfo_size; i++){
        if( abs(GenPartInfo_Status->at(i)) != 23 ) continue;//remove incoming/intermediate particles
        if( abs(GenPartInfo_PdgID->at(i)) > 6 ) continue;//exclude top quark & other particles
        //--------------------
        TLorentzVector genParticle;
        genParticle.SetPtEtaPhiM(GenPartInfo_Pt->at(i), GenPartInfo_Eta->at(i), GenPartInfo_Phi->at(i), GenPartInfo_Mass->at(i));
        double delta_R = jet.DeltaR(genParticle);
        //select quark & require min(deltaR)
        if( delta_R < 0.4 && delta_R < delta_R_min){
            index = i;//record the matched genParticle
            delta_R_min = delta_R;
            truelove = genParticle;
            truelove_PdgID = GenPartInfo_PdgID->at(i);
        }
    }//end of gen loop

    /*
    if(index>0){
        printf("[check-tqh] \n");
            printf("Status = %3d, ", GenPartInfo_Status->at(index));
            printf("PdgID = %3d, ", GenPartInfo_PdgID->at(index));
            printf("Pt = %6.2f, ", GenPartInfo_Pt->at(index));
            printf("Eta = %9.2f, ", GenPartInfo_Eta->at(index));
            printf("Phi = %6.2f, ", GenPartInfo_Phi->at(index));
            printf("Mass = %6.2f, ", GenPartInfo_Mass->at(index));
            printf("MomPdgID = %5d\n", GenPartInfo_MomPdgID->at(index));
    } else{
        printf("[check-tqh] index<0\n");
    }
    */

    bool is_tqh_quark = ( index > 0 && abs(GenPartInfo_MomPdgID->at(index)) == 6 );
    return is_tqh_quark;

}
bool CheckBJetID(TLorentzVector jet, Int_t GenPartInfo_size, std::vector<int> *GenPartInfo_MomPdgID, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_Status, std::vector<int> *GenPartInfo_PdgID){
    TLorentzVector truelove;//we are looking for the right genParticle to match the jet.
    int index = -999, truelove_PdgID = -999; double delta_R_min = 999.;
    for(int i=0; i<GenPartInfo_size; i++){
        if( abs(GenPartInfo_Status->at(i)) != 23 ) continue;//remove incoming/intermediate particles
        if( abs(GenPartInfo_PdgID->at(i)) == 6 ) continue;//exclude top quark
        //--------------------
        TLorentzVector genParticle;
        genParticle.SetPtEtaPhiM(GenPartInfo_Pt->at(i), GenPartInfo_Eta->at(i), GenPartInfo_Phi->at(i), GenPartInfo_Mass->at(i));
        double delta_R = jet.DeltaR(genParticle);
        //select quark & require min(deltaR)
        if( abs(GenPartInfo_PdgID->at(i)) < 7 && delta_R < 0.4 && delta_R < delta_R_min){
        //if( abs(GenPartInfo_PdgID->at(i)) < 7 && delta_R < delta_R_min){
            index = i;//record the matched genParticle
            delta_R_min = delta_R;
            truelove = genParticle;
            truelove_PdgID = GenPartInfo_PdgID->at(i);
        }
    }//end of gen loop
    // bjet is bquark coming from top
    bool bjet_is_bquark = ( abs(truelove_PdgID) == 5 && abs(GenPartInfo_MomPdgID->at(index)) == 6);
    return bjet_is_bquark;
}
TLorentzVector GetGenParticle(TLorentzVector jet, Int_t GenPartInfo_size, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_Status, std::vector<int> *GenPartInfo_PdgID, std::vector<int> &index_GenParticles, int &genParticle_PdgID){
    TLorentzVector truelove;//we are looking for the right genParticle to match the jet.
    int index = -1, truelove_PdgID = -999; double delta_R_min = 999.;
    for(int i=0; i<GenPartInfo_size; i++){
        if( abs(GenPartInfo_Status->at(i)) != 23 ) continue;//remove incoming/intermediate particles
        if( abs(GenPartInfo_PdgID->at(i)) == 6 ) continue;//exclude top quark
        //bool isAvailable = checkAvailability(i, index_GenParticles); if(!isAvailable) continue;//the truelove of previous jets shall not become another one's truelove.
        //--------------------
        TLorentzVector genParticle;
        genParticle.SetPtEtaPhiM(GenPartInfo_Pt->at(i), GenPartInfo_Eta->at(i), GenPartInfo_Phi->at(i), GenPartInfo_Mass->at(i));
        double delta_R = jet.DeltaR(genParticle);
        //select quark & require min(deltaR)
        if( abs(GenPartInfo_PdgID->at(i)) < 7 && delta_R < 0.4 && delta_R < delta_R_min){
        //if( abs(GenPartInfo_PdgID->at(i)) < 7 && delta_R < delta_R_min){
            index = i;//record the matched genParticle
            delta_R_min = delta_R;
            truelove = genParticle;
            truelove_PdgID = GenPartInfo_PdgID->at(i);
        }
    }//end of gen loop
    index_GenParticles.push_back(index);
    genParticle_PdgID = truelove_PdgID;
    //store particle info
    return truelove;
}

bool checkAvailability(int index, std::vector<int> ID_IsChosen){
    bool result = true;//if pass the following for loop, the genParticle of the index is available.
    for(std::size_t i=0; i<ID_IsChosen.size(); ++i){
        if(index == ID_IsChosen[i]){ result = false; break; } 
    }
    return result;
}
//}}}
// ##### print gen info{{{
void print_matched_gen_info(TLorentzVector jet, Int_t GenPartInfo_size, std::vector<int> *GenPartInfo_MomPdgID, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_Status, std::vector<int> *GenPartInfo_PdgID){
    TLorentzVector truelove;//we are looking for the right genParticle to match the jet.
    int index = -999, truelove_PdgID = -999; double delta_R_min = 999.;
    for(int i=0; i<GenPartInfo_size; i++){
        if( abs(GenPartInfo_Status->at(i)) != 23 ) continue;//remove incoming/intermediate particles
        if( abs(GenPartInfo_PdgID->at(i)) > 6 ) continue;//exclude top quark & other particles
        //--------------------
        TLorentzVector genParticle;
        genParticle.SetPtEtaPhiM(GenPartInfo_Pt->at(i), GenPartInfo_Eta->at(i), GenPartInfo_Phi->at(i), GenPartInfo_Mass->at(i));
        double delta_R = jet.DeltaR(genParticle);
        //select quark & require min(deltaR)
        if( delta_R < 0.4 && delta_R < delta_R_min){
            index = i;//record the matched genParticle
            delta_R_min = delta_R;
            truelove = genParticle;
            truelove_PdgID = GenPartInfo_PdgID->at(i);
        }
    }//end of gen loop

    if(index>0){
        printf("[check-tqh] \n");
            printf("Status = %3d, ", GenPartInfo_Status->at(index));
            printf("PdgID = %3d, ", GenPartInfo_PdgID->at(index));
            printf("Pt = %6.2f, ", GenPartInfo_Pt->at(index));
            printf("Eta = %9.2f, ", GenPartInfo_Eta->at(index));
            printf("Phi = %6.2f, ", GenPartInfo_Phi->at(index));
            printf("Mass = %6.2f, ", GenPartInfo_Mass->at(index));
            printf("MomPdgID = %5d\n", GenPartInfo_MomPdgID->at(index));
    } else{
        printf("[check-tqh] Not found.\n");
    }
}
//}}}
//### report{{{
void kinematics_info(const char* Title, TLorentzVector Particle){
        printf("(%s) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f, Energy = %6.2f, Mass = %6.2f\n", Title, Particle.Pt(), Particle.Eta(), Particle.Phi(), Particle.Energy(), Particle.M());
}
void kinematics_report(const char* recoTitle, TLorentzVector recoParticle, int id_recoParticle, const char* genTitle, TLorentzVector genParticle, int genParticle_PdgID){
        double delta_R = genParticle.DeltaR(recoParticle);
        printf("(%s) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f, Energy = %6.2f, Mass = %6.2f, id = %d\n", recoTitle, recoParticle.Pt(), recoParticle.Eta(), recoParticle.Phi(), recoParticle.Energy(), recoParticle.M(), id_recoParticle);
        printf("(%s) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f, Energy = %6.2f, Mass = %6.2f, id = %d, delta_R = %6.2f\n", genTitle, genParticle.Pt(), genParticle.Eta(), genParticle.Phi(), genParticle.Energy(), genParticle.M(), genParticle_PdgID, delta_R);
}
void hist_bin_fraction(TH1D *hist, const char* title, int entries_matched_in_bin2){
    int total_entries = hist->GetEntries();
    printf("[INFO-hist] %s (%d)\n", title, total_entries);
    printf("[INFO-hist] fractions = ");
    for(int i=0; i<3; ++i){
        if(i<2) printf("%6.4f, ", (double) hist->GetBinContent(i+1) / (double)total_entries);
        else{
            int accumulation = 0;
            for(int j=i; j<10; ++j) accumulation += (double) hist->GetBinContent(j+1);
            printf("%6.4f; ", i, (double) accumulation / (double)total_entries);
            printf("%6.4f, ", i, (double) entries_matched_in_bin2 / (double)total_entries);
            printf("%6.4f\n", i, (double) entries_matched_in_bin2 / (double)hist->GetBinContent(2));
        }
    }
}
void hist_report(TH1D *hist, const char* chi2_type){
    double mean = hist->GetMean();
    double sigma = hist->GetMeanError();
    double sigma2 = hist->GetRMS();
    int bin_maximum = hist->GetMaximumBin();
    double maxbin_lower_edge = hist->GetBinLowEdge(bin_maximum);
    double bin_width = hist->GetBinWidth(bin_maximum);
    double maxbin_upper_edge = maxbin_lower_edge + bin_width;

    //printf("[INFO] chi2 %s: mean = %6.2f, sigma = %6.2f, %6.2f\n", chi2_type, mean, sigma, sigma2);
    printf("[INFO] chi2 %s: mean = %6.2f, sigma = %6.2f, mode lies in [%6.2f, %6.2f]\n", chi2_type, mean, sigma2, maxbin_lower_edge, maxbin_upper_edge);
}
//}}}
//### chi2{{{
double Chi2_calculator_simple(double w_mass, double t_mass){
    TVectorD vec_mass(2);
    vec_mass(0) = w_mass - w_boson_mass;
    vec_mass(1) = t_mass - top_quark_mass;
    TMatrixD matrix(2,2);
    matrix(0,0) = 345.84; matrix(0,1) =   0.00;
    matrix(1,0) =   0.00; matrix(1,1) = 961.79;
    //--- ST ---//
    //matrix(0,0) = 305.98; matrix(0,1) = 0;
    //matrix(1,0) = 0;      matrix(1,1) = 787.07;
    return matrix.Invert()*vec_mass*vec_mass;
}
double Chi2_calculator_modified(double w_mass, double t_mass){
    TVectorD vec_mass(2);
    vec_mass(0) = w_mass - w_boson_mass;
    vec_mass(1) = t_mass - top_quark_mass;
    TMatrixD matrix(2,2);
    matrix(0,0) = 345.84; matrix(0,1) = 378.21;
    matrix(1,0) = 378.21; matrix(1,1) = 961.79;
    //--- ST ---//
    //matrix(0,0) = 305.98; matrix(0,1) = 323.97;
    //matrix(1,0) = 323.97; matrix(1,1) = 787.07;
    return matrix.Invert()*vec_mass*vec_mass;
}
double Chi2_calculator_improved(double w_mass, double t_mass, double fcnc_top_mass){
    TVectorD vec_mass(3);
    vec_mass(0) = w_mass - w_boson_mass;
    vec_mass(1) = t_mass - top_quark_mass;
    vec_mass(2) = fcnc_top_mass - top_quark_mass;
    TMatrixD matrix(3,3);
    matrix(0,0) = 345.84; matrix(0,1) = 378.21; matrix(0,2) =   3.26;
    matrix(1,0) = 378.21; matrix(1,1) = 961.79; matrix(1,2) =   3.78;
    matrix(2,0) =   3.26; matrix(2,1) =   3.78; matrix(2,2) = 272.53;
    //--- ST ---//
    //matrix(0,0) = 305.98; matrix(0,1) = 323.97;
    //matrix(1,0) = 323.97; matrix(1,1) = 787.07;
    return matrix.Invert()*vec_mass*vec_mass;
}
//###}}}
//### mkplots{{{
void MakeTwoPlots(TCanvas *c1, TH1D* hist_gen, TH1D* hist_reco, TLegend *legend, const char* name){
    hist_reco->SetStats(0);
    hist_reco->SetLineWidth(2);
    hist_reco->SetLineColor(kBlue);
    hist_reco->Draw();
    hist_gen->SetLineWidth(2);
    hist_gen->SetLineColor(kRed);
    hist_gen->Draw("same;hist");

    legend->Clear();
    legend->AddEntry(hist_gen, "gen-level", "l");
    legend->AddEntry(hist_reco, "reco", "l");
    legend->SetLineColor(0);
    legend->Draw("same");

    c1->SaveAs(Form("ntuples_skimmed/%s", name));
}
void MakePlots(TCanvas *c1, TH1D* hist, const char* title, const char* outputFile){
    hist->Draw();
    hist->SetTitle(title);
    hist->SetXTitle(title);
    hist->SetYTitle("Entries");
    hist->GetYaxis()->SetTitleOffset(1.4);
    hist->Write();
    c1->SaveAs(outputFile);
}
void MakeFinalPlots(TCanvas *c1, TH1D* hist_simple, TH1D* hist_modified, TH1D*hist_yfyj, TLegend *legend, const char* name){
    double max;
    double scale = 1.2;
    double max_simple = hist_simple -> GetMaximum();
    double max_modified = hist_modified -> GetMaximum();
    double max_yfyj = hist_yfyj -> GetMaximum();
    max = (max_modified > max_simple) ? max_modified : max_simple;
    max = (max_yfyj > max_modified) ? max_yfyj : max_modified;
    double overflow_simple = hist_simple->GetBinContent(hist_simple->GetNbinsX() + 1);
    double overflow_modified = hist_modified->GetBinContent(hist_modified->GetNbinsX() + 1);
    double overflow_yfyj = hist_yfyj->GetBinContent(hist_yfyj->GetNbinsX() + 1);
    //max = (max > overflow_simple) ? max : overflow_simple;
    //max = (max > overflow_modified) ? max : overflow_modified;
    //max = (max > overflow_yfyj) ? max : overflow_yfyj;

    hist_simple->SetStats(0);
    hist_simple->SetMaximum(max*scale);
    hist_simple->SetLineWidth(2);
    hist_simple->SetLineColor(kBlue);
    hist_simple->GetXaxis()->SetRange(1, hist_simple->GetNbinsX() + 1);
    hist_simple->Draw();
    hist_modified->SetLineWidth(2);
    hist_modified->SetLineColor(kRed);
    hist_modified->GetXaxis()->SetRange(1, hist_modified->GetNbinsX() + 1);
    hist_modified->Draw("same");
    hist_yfyj->SetLineWidth(2);
    hist_yfyj->SetLineColor(kGreen+4);
    hist_yfyj->GetXaxis()->SetRange(1, hist_yfyj->GetNbinsX() + 1);
    hist_yfyj->Draw("same");

    legend->Clear();
    legend->SetTextSize(0.03);
    legend->AddEntry(hist_yfyj, "leading-jets method", "l");
    legend->AddEntry(hist_simple, "simple #chi^{2} method", "l");
    legend->AddEntry(hist_modified, "modified #chi^{2} method", "l");
    legend->SetLineColor(0);
    legend->Draw("same");

    c1->SaveAs(Form("ntuples_skimmed/%s", name));
}
//###}}}
//### bool functions{{{
bool isThisMCsignal(char* dataset){
    if((string)dataset == "ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-T2HJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8") return true;
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
    return false;
}
bool isThisMultiFile(char* dataset){
    if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return true;
    if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return true;
    return false;
}
//###}}}
//### flashggStdTree {{{
flashggStdTreeParameters::flashggStdTreeParameters(){
    GenPartInfo_Pt = new std::vector<float>;
    GenPartInfo_Eta = new std::vector<float>;
    GenPartInfo_Phi = new std::vector<float>;
    GenPartInfo_Mass = new std::vector<float>;
    GenPartInfo_PdgID = new std::vector<int>;
    GenPartInfo_Status = new std::vector<int>;
    GenPartInfo_nMo = new std::vector<int>;
    GenPartInfo_nDa = new std::vector<int>;
    GenPartInfo_isHardProcess = new std::vector<bool>;
    GenPartInfo_fromHardProcessFinalState = new std::vector<bool>;
    GenPartInfo_isPromptFinalState = new std::vector<bool>;
    GenPartInfo_isDirectPromptTauDecayProductFinalState = new std::vector<bool>;
    GenPartInfo_MomPdgID = new std::vector<int>;
    GenPartInfo_MomStatus = new std::vector<int>;
    GenPartInfo_MomPt = new std::vector<float>;
    GenPartInfo_MomEta = new std::vector<float>;
    GenPartInfo_MomPhi = new std::vector<float>;
    GenPartInfo_MomMass = new std::vector<float>;
    //------------------------
    JetInfo_Pt = new std::vector<float>;
    JetInfo_Eta = new std::vector<float>;
    JetInfo_Phi = new std::vector<float>;
    JetInfo_Mass = new std::vector<float>;
    JetInfo_Energy = new std::vector<float>;
    JetInfo_GenFlavor = new std::vector<int>;
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
    delete GenPartInfo_isHardProcess;
    delete GenPartInfo_fromHardProcessFinalState;
    delete GenPartInfo_isPromptFinalState;
    delete GenPartInfo_isDirectPromptTauDecayProductFinalState;
    delete GenPartInfo_MomPdgID;
    delete GenPartInfo_MomStatus;
    delete GenPartInfo_MomPt;
    delete GenPartInfo_MomEta;
    delete GenPartInfo_MomPhi;
    delete GenPartInfo_MomMass;
    //------------------------
    delete JetInfo_Pt;
    delete JetInfo_Eta;
    delete JetInfo_Phi;
    delete JetInfo_Mass;
    delete JetInfo_Energy;
    delete JetInfo_GenFlavor;
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
    flashggStdTree->SetBranchAddress("GenPartInfo.isHardProcess", &GenPartInfo_isHardProcess);
    flashggStdTree->SetBranchAddress("GenPartInfo.fromHardProcessFinalState", &GenPartInfo_fromHardProcessFinalState);
    flashggStdTree->SetBranchAddress("GenPartInfo.isPromptFinalState", &GenPartInfo_isPromptFinalState);
    flashggStdTree->SetBranchAddress("GenPartInfo.isDirectPromptTauDecayProductFinalState", &GenPartInfo_isDirectPromptTauDecayProductFinalState);
    flashggStdTree->SetBranchAddress("GenPartInfo.MomPdgID", &GenPartInfo_MomPdgID);
    flashggStdTree->SetBranchAddress("GenPartInfo.MomStatus", &GenPartInfo_MomStatus);
    flashggStdTree->SetBranchAddress("GenPartInfo.MomPt", &GenPartInfo_MomPt);
    flashggStdTree->SetBranchAddress("GenPartInfo.MomEta", &GenPartInfo_MomEta);
    flashggStdTree->SetBranchAddress("GenPartInfo.MomPhi", &GenPartInfo_MomPhi);
    flashggStdTree->SetBranchAddress("GenPartInfo.MomMass", &GenPartInfo_MomMass);
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
    flashggStdTree->SetBranchAddress("JetInfo.GenFlavor", &JetInfo_GenFlavor);
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
///### }}}
//### myTreeClass {{{
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
    mytree -> Branch("Mass_w_candidate_chi2_simple", &Mass_w_candidate_chi2_simple, "Mass_w_candidate_chi2_simple/F");
    mytree -> Branch("Mass_top_candidate_chi2_simple", &Mass_top_candidate_chi2_simple, "Mass_top_candidate_chi2_simple/F");
    mytree -> Branch("Mass_w_candidate_chi2_modified", &Mass_w_candidate_chi2_modified, "Mass_w_candidate_chi2_modified/F");
    mytree -> Branch("Mass_top_candidate_chi2_modified", &Mass_top_candidate_chi2_modified, "Mass_top_candidate_chi2_modified/F");
    mytree -> Branch("Mass_gen_w_candidate_chi2_simple", &Mass_gen_w_candidate_chi2_simple, "Mass_gen_w_candidate_chi2_simple/F");
    mytree -> Branch("Mass_gen_top_candidate_chi2_simple", &Mass_gen_top_candidate_chi2_simple, "Mass_gen_top_candidate_chi2_simple/F");
    mytree -> Branch("Mass_gen_w_candidate_chi2_modified", &Mass_gen_w_candidate_chi2_modified, "Mass_gen_w_candidate_chi2_modified/F");
    mytree -> Branch("Mass_gen_top_candidate_chi2_modified", &Mass_gen_top_candidate_chi2_modified, "Mass_gen_top_candidate_chi2_modified/F");
    //------------------------
    mytree -> Branch("EvtInfo_totalEntry_before_preselection", &EvtInfo_totalEntry_before_preselection, "EvtInfo_totalEntry_before_preselection/I");
    mytree -> Branch("EvtInfo_NormalizationFactor_lumi", &EvtInfo_NormalizationFactor_lumi, "EvtInfo_NormalizationFactor_lumi/F");
    mytree -> Branch("EvtInfo_NPu", &EvtInfo_NPu, "EvtInfo_NPu/I");
    mytree -> Branch("EvtInfo_Rho", &EvtInfo_Rho, "EvtInfo_Rho/F");
    mytree -> Branch("EvtInfo_NVtx", &EvtInfo_NVtx, "EvtInfo_NVtx/I");
    mytree -> Branch("EvtInfo_genweight", &EvtInfo_genweight, "EvtInfo_genweight/F");
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
    mytree -> Branch("DiPhoInfo_subleadPt", &DiPhoInfo_subleadPt, "DiPhoInfo_subleadPt/F");
    mytree -> Branch("DiPhoInfo_subleadEta", &DiPhoInfo_subleadEta, "DiPhoInfo_subleadEta/F");
    mytree -> Branch("DiPhoInfo_subleadPhi", &DiPhoInfo_subleadPhi, "DiPhoInfo_subleadPhi/F");
    mytree -> Branch("DiPhoInfo_subleadE", &DiPhoInfo_subleadE, "DiPhoInfo_subleadE/F");
    mytree -> Branch("DiPhoInfo_subleadhoe", &DiPhoInfo_subleadhoe, "DiPhoInfo_subleadhoe/F");
    mytree -> Branch("DiPhoInfo_subleadIDMVA", &DiPhoInfo_subleadIDMVA, "DiPhoInfo_subleadIDMVA/F");
    //------------------------
    mytree -> Branch("ElecInfo_Size", &ElecInfo_Size, "ElecInfo_Size/I");
    mytree -> Branch("MuonInfo_Size", &MuonInfo_Size, "MuonInfo_Size/I");
    mytree -> Branch("num_leptons", &num_leptons, "num_leptons/I");// # of selected objects.
    mytree -> Branch("num_electrons", &num_electrons, "num_electrons/I");// # of selected objects.
    mytree -> Branch("num_muons", &num_muons, "num_muons/I");// # of selected objects.
    mytree -> Branch("ElecInfo_electron_pt", &ElecInfo_electron_pt);
    mytree -> Branch("ElecInfo_electron_eta", &ElecInfo_electron_eta);
    mytree -> Branch("ElecInfo_electron_phi", &ElecInfo_electron_phi);
    mytree -> Branch("ElecInfo_electron_energy", &ElecInfo_electron_energy);
    mytree -> Branch("ElecInfo_electron_diphoton_deltaR", &ElecInfo_electron_diphoton_deltaR);
    mytree -> Branch("ElecInfo_electron_leadingPhoton_deltaR", &ElecInfo_electron_leadingPhoton_deltaR);
    mytree -> Branch("ElecInfo_electron_subleadingPhoton_deltaR", &ElecInfo_electron_subleadingPhoton_deltaR);
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
    //mytree -> Branch("num_btagged_jets", &num_btagged_jets, "num_btagged_jets/I");
    //mytree -> Branch("num_nonbtagged_jets", &num_nonbtagged_jets, "num_nonbtagged_jets/I");
    //------------------------
    //mytree -> Branch("inv_mass_dijet", &inv_mass_dijet, "inv_mass_dijet/F");
    //mytree -> Branch("inv_mass_diphoton", &inv_mass_diphoton, "inv_mass_diphoton/F");
    //mytree -> Branch("inv_mass_tbw", &inv_mass_tbw, "inv_mass_tbw/F");
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
    Mass_w_candidate_chi2_simple = 0;
    Mass_top_candidate_chi2_simple = 0;
    Mass_w_candidate_chi2_modified = 0;
    Mass_top_candidate_chi2_modified = 0;
    Mass_gen_w_candidate_chi2_simple = 0;
    Mass_gen_top_candidate_chi2_simple = 0;
    Mass_gen_w_candidate_chi2_modified = 0;
    Mass_gen_top_candidate_chi2_modified = 0;

    inv_mass_dijet = 0;
    inv_mass_diphoton = 0;
    inv_mass_tbw = 0;
    //------------------------
    JetInfo_dijet_delta_eta = 0;
    JetInfo_dijet_delta_phi = 0;
    JetInfo_dijet_delta_angle = 0;
    //------------------------
}
//### }}}
//}}}