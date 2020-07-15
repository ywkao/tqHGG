// vim: set fdm=marker: {{{
//***************************************************************************
//
// FileName    : generalChiSquareStudy_hadronic.cpp
// Purpose     : Develop top reconstruction methods in hadronic channel (test code) 
// Description : several methods for top reconstruction are tested
//             : (1) leading jets method (only high pt jets are selected for reconstruction)
//             : (2) simple chi-2 method (without consider covariance term)
//             : (3) modified chi-2 method (2X2 matrix)
//             : (4) improved chi-2 method (3X3 matrix)
// Author      : Yu-Wei Kao [ykao@cern.ch]
//
//***************************************************************************
#include "../include/generalChiSquareStudy.C"
#include "../include/cross_section.h"
#include "../include/hist_factory.h"
#include "../include/print_tool.C"
#include "../include/plotHelper.C"
//--- control bjet selection ---//
bool printSelectedJetsInfo = false;
//bool bool_bjet_is_loose  = true;
//bool bool_bjet_is_medium = false;
//bool bool_bjet_is_tight  = false;
//bool bool_num_bjets_is_exactly_one = false;
//bool bool_num_bjets_is_atleast_one = !bool_num_bjets_is_exactly_one;
//}}}
int CONSTRAINT_NUM_JETS = 3; //TT=4, ST=3
int main(int argc, char *argv[]){
    // I/O and event info
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
    //sprintf(input_file, "%s", Form("%s/ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root", dir))     ; treeReader.AddSingleRootFile(input_file) ;
    //sprintf(input_file, "%s", Form("%s/TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root", dir))   ; treeReader.AddSingleRootFile(input_file) ;
    //sprintf(input_file, "%s", Form("%s/TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8.root", dir))   ; treeReader.AddSingleRootFile(input_file) ;
    sprintf(input_file, "%s", Form("%s/ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root", dir))     ; treeReader.AddSingleRootFile(input_file) ;
    //sprintf(input_file, "%s", Form("%s/TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root", dir)) ; treeReader.AddSingleRootFile(input_file) ;
    //sprintf(input_file, "%s", Form("%s/TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8.root", dir)) ; treeReader.AddSingleRootFile(input_file) ;
    
    // Opening{{{
    // init{{{
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
    printf("[INFO] NormalizationFactor = %f!\n");
    //}}}
    // histograms{{{
    //==================================================//
    //--------   Histograms for bjet rate-WP  ----------//
    //==================================================//
    const char* wp[3] = {"Loose", "Medium", "Tight"};
    TH1D *hist_rate_bquark_is_included           = new TH1D("hist_rate_bquark_is_included", ";;rate", 3, 0, 3);
    TH1D *hist_rate_leading_bjet_is_bquark       = new TH1D("hist_rate_leading_bjet_is_bquark", ";;rate", 3, 0, 3);
    TH1D *hist_rate_permuting_bjet_is_bquark_sim = new TH1D("hist_rate_permuting_bjet_is_bquark_sim", ";;rate", 3, 0, 3);
    TH1D *hist_rate_permuting_bjet_is_bquark_mod = new TH1D("hist_rate_permuting_bjet_is_bquark_mod", ";;rate", 3, 0, 3);
    TH1D *hist_rate_permuting_bjet_is_bquark_imp = new TH1D("hist_rate_permuting_bjet_is_bquark_imp", ";;rate", 3, 0, 3);
    TH1D *hist_entries_bjet_wp_selected          = new TH1D("hist_entries_bjet_wp_selected", ";;Entries", 3, 0, 3);
    TH1D *hist_entries_bjet_wp_matched           = new TH1D("hist_entries_bjet_wp_matched", ";;Entries", 3, 0, 3);
    for(int i=0; i<3; ++i) hist_rate_bquark_is_included           -> GetXaxis() -> SetBinLabel(i+1, wp[i]);
    for(int i=0; i<3; ++i) hist_rate_leading_bjet_is_bquark       -> GetXaxis() -> SetBinLabel(i+1, wp[i]);
    for(int i=0; i<3; ++i) hist_rate_permuting_bjet_is_bquark_sim -> GetXaxis() -> SetBinLabel(i+1, wp[i]);
    for(int i=0; i<3; ++i) hist_rate_permuting_bjet_is_bquark_mod -> GetXaxis() -> SetBinLabel(i+1, wp[i]);
    for(int i=0; i<3; ++i) hist_rate_permuting_bjet_is_bquark_imp -> GetXaxis() -> SetBinLabel(i+1, wp[i]);
    for(int i=0; i<3; ++i) hist_entries_bjet_wp_selected          -> GetXaxis() -> SetBinLabel(i+1, wp[i]);
    for(int i=0; i<3; ++i) hist_entries_bjet_wp_matched           -> GetXaxis() -> SetBinLabel(i+1, wp[i]);
    //==================================================//
    //--------   Histograms for mass spectra  ----------//
    //==================================================//
    TH1D *hist_num_bjets_tight = new TH1D("hist_num_bjets_tight", ";Number of tight b-tagged jets;Entries", 10, 0, 10);
    TH1D *hist_num_bjets_medium = new TH1D("hist_num_bjets_medium", ";Number of medium b-tagged jets;Entries", 10, 0, 10);
    TH1D *hist_num_bjets_loose = new TH1D("hist_num_bjets_loose", ";Number of loose b-tagged jets;Entries", 10, 0, 10);
    TH1D *hist_num_bjets_tight_matched = new TH1D("hist_num_bjets_tight_matched", ";Number of tight b-tagged jets;Entries", 10, 0, 10);
    TH1D *hist_num_bjets_medium_matched = new TH1D("hist_num_bjets_medium_matched", ";Number of medium b-tagged jets;Entries", 10, 0, 10);
    TH1D *hist_num_bjets_loose_matched = new TH1D("hist_num_bjets_loose_matched", ";Number of loose b-tagged jets;Entries", 10, 0, 10);
    TH1D *hist_num_gen_bquark = new TH1D("hist_num_gen_bquark", ";Number of b quarks (gen-level);Entries", 10, 0, 10);
    TH1D *hist_num_gen_light_quark = new TH1D("hist_num_gen_light_quark", ";Number of light quarks (gen-level);Entries", 10, 0, 10);
    TH1D *hist_num_gen_quarks = new TH1D("hist_num_gen_quarks", ";Number of quarks (gen-level);Entries", 10, 0, 10);
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
    //hist_factory{{{
    char output_histDir[128] = "result_top_reco_study/hist_factory_hadronic";
    hist_factory hf_wboson (output_histDir, "wboson", "modified", 200);
    hist_factory hf_tbw (output_histDir, "tbw", "modified", 350);
    hist_factory hf_tqh (output_histDir, "tqh", "modified", 350);
    hist_factory hf_tqh_matched (output_histDir, "tqh", "matched", 350);
    //}}}
    // counters{{{
    int counter_irregular_disc = 0;
    int counter_coeff_D = 0;
    int counter_coeff_D_gen = 0;
    int counter_coeff_D_isNegative = 0;
    int counter_coeff_D_isNegative_gen = 0;
    //int Nevents_pass_selection = 0;
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
    int counter_selectedJets_tbwCanBeReconstructed = 0;
    int counter_selectedJets_tqhCanBeReconstructed = 0;
    int counter_selectedJets_sigCanBeReconstructed = 0;

    // for wp = loose, medium, tight
    int Nevents_pass_selection[3][2] = {};
    int counter_bjet_is_matched_improved[3][2] = {};
    int counter_bjet_is_matched_modified[3][2] = {};
    int counter_bjet_is_matched_simple[3][2] = {};
    int counter_tbwIsCorrectlyMatched_yfyj[3][2] = {};
    int counter_tbwIsCorrectlyMatched_simple[3][2] = {};
    int counter_tbwIsCorrectlyMatched_modified[3][2] = {};
    int counter_tbwIsCorrectlyMatched_improved[3][2] = {};
    int counter_tqhIsCorrectlyMatched_yfyj[3][2] = {};
    int counter_tqhIsCorrectlyMatched_simple[3][2] = {};
    int counter_tqhIsCorrectlyMatched_modified[3][2] = {};
    int counter_tqhIsCorrectlyMatched_improved[3][2] = {};
    int counter_sigIsCorrectlyReco_yfyj[3][2] = {};
    int counter_sigIsCorrectlyReco_simple[3][2] = {};
    int counter_sigIsCorrectlyReco_modified[3][2] = {};
    int counter_sigIsCorrectlyReco_improved[3][2] = {};
    printf("[curiosity] Nevents_pass_selection ");
    for(int i=0; i<2; ++i){
        for(int j=0; j<3; ++j){
            printf("%d ", Nevents_pass_selection[j][i]);
        }
    }
    printf("\n");

    int counter_mismatched = 0;
    int counter_intersection = 0;
    int counter_M_notC = 0;
    int counter_notM_C = 0;
    int counter_notM_notC = 0;
    int counter_bjj_can_be_matched = 0;
    //}}}
        int counter_index_bjet_method_2 = 0;
        int anti_counter_index_bjet_method_2 = 0;
        int counter_genmatched_index_bjet_unphysical = 0;
        int anti_counter_genmatched_index_bjet_unphysical = 0;
        int counter_not_matched_but_is_in_the_list = 0;
        int counter_difference = 0;
    //}}}
    //===============================================//
    //----------------   Event Loop  ----------------//
    //===============================================//
    //for(int ientry=0; ientry<nentries; ientry++){
    for(int ientry=0; ientry<1000; ientry++){
        //preselection, b-jet, geninfo{{{
        treeReader.flashggStdTree->GetEntry(ientry);//load data
        // Pre-process{{{
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
        // hadronic event selection{{{
        if(mytree.num_leptons>0) continue;
        if(mytree.num_jets<CONSTRAINT_NUM_JETS) continue;
        //}}}
        // store EventInfo{{{
        mytree.EvtInfo_NPu = treeReader.EvtInfo_NPu;
        mytree.EvtInfo_Rho = treeReader.EvtInfo_Rho;
        mytree.EvtInfo_NVtx = treeReader.EvtInfo_NVtx;
        mytree.EvtInfo_genweight = treeReader.EvtInfo_genweight;
        //}}}
        //### check b quark info, store diphoton and jets{{{
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
        //}}}
        // b-jet{{{
        TLorentzVector bjet;
        TLorentzVector genParticle_bjet;
        int num_bjets = 0, index_bjet = -999;
        int num_bjets_tight = 0, num_bjets_loose = 0, num_bjets_medium = 0;
        int index_leading_bjet_tight = -999, index_leading_bjet_loose = -999, index_leading_bjet_medium = -999;
        //std::vector<int> indices_bjet, indices_tight, indices_loose, indices_medium;
        //std::vector<TLorentzVector> bjets, bjets_tight, bjets_loose, bjets_medium;
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
        ////determine leading bjet according to chosen WP (skipped){{{
        //if(bool_bjet_is_loose){
        //    indices_bjet = indices_loose;
        //    bjets = bjets_loose;
        //    index_bjet = index_leading_bjet_loose;
        //    if(index_bjet != -999){
        //        bjet = bjets_loose[0];
        //        num_bjets = num_bjets_loose;
        //    }
        //}
        //if(bool_bjet_is_medium){
        //    indices_bjet = indices_medium;
        //    bjets = bjets_medium;
        //    index_bjet = index_leading_bjet_medium;
        //    if(index_bjet != -999){
        //        bjet = bjets_medium[0];
        //        num_bjets = num_bjets_medium;
        //    }
        //}
        //if(bool_bjet_is_tight){
        //    indices_bjet = indices_tight;
        //    bjets = bjets_tight;
        //    index_bjet = index_leading_bjet_tight;
        //    if(index_bjet != -999){
        //        bjet = bjets_tight[0];
        //        num_bjets = num_bjets_tight;
        //    }
        //}
        ////}}}
        
        //bool pass_bjets_multiplicity_selection;
        //if(bool_num_bjets_is_exactly_one) pass_bjets_multiplicity_selection = num_bjets == 1;
        //if(bool_num_bjets_is_atleast_one) pass_bjets_multiplicity_selection = num_bjets >= 1;
        //if(!pass_bjets_multiplicity_selection) continue;
        
        //### }}}
        //check geninfo(skipped){{{
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
        //}}}
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
        // pick up indices of jets from hard process according to MC truth{{{
        if(printSelectedJetsInfo){
            printf("\nientry = %d/%d\n", ientry+1, nentries+1);
            printf("[CHECK] num_jets = %d\n", mytree.num_jets);
        }
        int counter_momPdgID_is_wboson = 0;
        int counter_is_bquarkFromSMtop = 0;
        int counter_is_quarkFromFCNtop = 0;

        int jetIndex_is_bquarkFromSMtop = -999;
        int jetIndex_is_quarkFromFCNtop = -999;
        std::vector<int> jetIndex_momPdgID_is_wboson(2, -999);
        
        std::vector<int> jetIndex_is_bquarkFromSMtop_moreThanOne(2, -999);
        std::vector<int> jetIndex_is_quarkFromFCNtop_moreThanOne(2, -999);
        std::vector<int> genIndex_is_quarkFromWboson(2, -999);
        for(int i=0; i<mytree.num_jets; ++i){
            //retrieve index of genParticle, matched ID, and matched mom ID
            int index_gen, Matched_PdgID, Matched_MomPdgID;
            if(printSelectedJetsInfo) kinematics_info("jet", Jets[i], i);
            obtain_gen_matched_ID(printSelectedJetsInfo, Jets[i], index_gen, Matched_PdgID, Matched_MomPdgID, treeReader.GenPartInfo_size,\
                                  treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                                  treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
            //counting
            if( abs(Matched_MomPdgID)==24 ) counter_momPdgID_is_wboson += 1;
            if( abs(Matched_PdgID)==5 && abs(Matched_MomPdgID)==6 ) counter_is_bquarkFromSMtop += 1;
            if( abs(Matched_PdgID)!=5 && abs(Matched_MomPdgID)==6 ) counter_is_quarkFromFCNtop += 1;

            //index
            if( abs(Matched_MomPdgID)==24 ){
                if(counter_momPdgID_is_wboson==1) jetIndex_momPdgID_is_wboson[0] = i;
                if(counter_momPdgID_is_wboson==2) jetIndex_momPdgID_is_wboson[1] = i;
            }
            if( abs(Matched_PdgID)==5 && abs(Matched_MomPdgID)==6 ) jetIndex_is_bquarkFromSMtop = i;
            if( abs(Matched_PdgID)!=5 && abs(Matched_MomPdgID)==6 ) jetIndex_is_quarkFromFCNtop = i;

            // fix bug when w quark is matched twice
            if( abs(Matched_MomPdgID)==24 ){
                if(counter_momPdgID_is_wboson==1) genIndex_is_quarkFromWboson[0] = index_gen;
                if(counter_momPdgID_is_wboson==2) genIndex_is_quarkFromWboson[1] = index_gen;
            }

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
        //--------------------------------------------------
        // fix bug when w quark is matched twice{{{
        bool w_quark_is_matched_twice = genIndex_is_quarkFromWboson[0] == genIndex_is_quarkFromWboson[1];
        //}}}
        // fix bug when bquark is matched twice{{{
        double deltaR_jet_bquark[2] = {0};
        if( counter_is_bquarkFromSMtop == 2 ){ 
            //pick the jet with smaller deltaR (jet, bquark)
            int Matched_PdgID, Matched_MomPdgID;
            deltaR_jet_bquark[0] = obtain_deltaR(printSelectedJetsInfo, Jets[jetIndex_is_bquarkFromSMtop_moreThanOne[0]],\
                                  Matched_PdgID, Matched_MomPdgID, treeReader.GenPartInfo_size,\
                                  treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                                  treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
            deltaR_jet_bquark[1] = obtain_deltaR(printSelectedJetsInfo, Jets[jetIndex_is_bquarkFromSMtop_moreThanOne[1]],\
                                  Matched_PdgID, Matched_MomPdgID, treeReader.GenPartInfo_size,\
                                  treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                                  treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
            jetIndex_is_bquarkFromSMtop = (deltaR_jet_bquark[0] < deltaR_jet_bquark[1]) ? jetIndex_is_bquarkFromSMtop_moreThanOne[0] : jetIndex_is_bquarkFromSMtop_moreThanOne[1];
            //printf("[check] chosen correct index of bjet: %d\n", jetIndex_is_bquarkFromSMtop);
            //jetIndex_is_bquarkFromSMtop = -999;
        }
        //}}}
        // fix bug when q-tqh is matched twice{{{
        double deltaR_jet_tqh_quark[2] = {0};
        if( counter_is_quarkFromFCNtop == 2 ){ 
            //pick the jet with smaller deltaR (jet, bquark)
            int Matched_PdgID, Matched_MomPdgID;
            deltaR_jet_tqh_quark[0] = obtain_deltaR(printSelectedJetsInfo, Jets[jetIndex_is_quarkFromFCNtop_moreThanOne[0]],\
                                  Matched_PdgID, Matched_MomPdgID, treeReader.GenPartInfo_size,\
                                  treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                                  treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
            deltaR_jet_tqh_quark[1] = obtain_deltaR(printSelectedJetsInfo, Jets[jetIndex_is_quarkFromFCNtop_moreThanOne[1]],\
                                  Matched_PdgID, Matched_MomPdgID, treeReader.GenPartInfo_size,\
                                  treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                                  treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
            jetIndex_is_quarkFromFCNtop = (deltaR_jet_tqh_quark[0] < deltaR_jet_tqh_quark[1]) ? jetIndex_is_quarkFromFCNtop_moreThanOne[0] : jetIndex_is_quarkFromFCNtop_moreThanOne[1];
        }
        //}}}

        bool exists_unphysical_value = ( jetIndex_is_bquarkFromSMtop < 0 ) || ( jetIndex_momPdgID_is_wboson[0] < 0 ) || ( jetIndex_momPdgID_is_wboson[1] < 0 );
        if(!exists_unphysical_value) counter_bjj_can_be_matched += 1;
        bool tbwCanBeReconstructed = (counter_momPdgID_is_wboson >= 2) && (counter_is_bquarkFromSMtop >= 1);
        if(tbwCanBeReconstructed) counter_selectedJets_tbwCanBeReconstructed += 1;
        bool tqhCanBeReconstructed = (counter_is_quarkFromFCNtop >= 1);
        if(tqhCanBeReconstructed) counter_selectedJets_tqhCanBeReconstructed += 1;
        bool sigCanBeReconstructed = tbwCanBeReconstructed && tqhCanBeReconstructed && !w_quark_is_matched_twice;
        if(sigCanBeReconstructed) counter_selectedJets_sigCanBeReconstructed += 1;
        //}}}
        //}}}
        //b-jet rates-WP study with gen-info{{{
        //printf("[check] indices = "); for(int i=0; i<indices_bjet.size(); ++i) printf("%d ", indices_bjet[i]); printf("\n");
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
        
        // gen-matched indices
        std::vector<int> index_jet_gen_matched = { jetIndex_is_bquarkFromSMtop, jetIndex_momPdgID_is_wboson[0], jetIndex_momPdgID_is_wboson[1], jetIndex_is_quarkFromFCNtop };

        // Multiplicity & WP
        std::vector<bool> bool_multiplicity_condition = {false, true};
        for(std::size_t ibool=0; ibool!=bool_multiplicity_condition.size(); ++ibool){
            bool bool_num_bjets_is_exactly_one = bool_multiplicity_condition[ibool];
            bool bool_num_bjets_is_atleast_one = !bool_num_bjets_is_exactly_one;
            std::vector<std::vector<int>> bjet_indices_wp = {indices_loose, indices_medium, indices_tight};
            for(std::size_t iwp=0; iwp!=bjet_indices_wp.size(); ++iwp){
                if(ibool!=0) continue; // consider loose WP only

                std::vector<int> indices_bjet = bjet_indices_wp[iwp];
                // bjet multiplicity selection{{{
                bool pass_bjets_multiplicity_selection;
                if(bool_num_bjets_is_exactly_one) pass_bjets_multiplicity_selection = indices_bjet.size() == 1;
                if(bool_num_bjets_is_atleast_one) pass_bjets_multiplicity_selection = indices_bjet.size() >= 1;
                if(!pass_bjets_multiplicity_selection) continue;
                //}}}
                hist_entries_bjet_wp_selected -> Fill(iwp);

                //--- test: mis-id bjet leads to the Hut and Hct difference ---//
                int index_bjet_method_2 = std::max_element(btag_scores.begin(), btag_scores.end()) - btag_scores.begin();
                //int index_bjet_method_2 = jetIndex_is_bquarkFromSMtop; //test purpose
                // to be further studied{{{
                if(ibool==0 && iwp==0){
                    if(index_bjet_method_2 == jetIndex_is_bquarkFromSMtop) counter_index_bjet_method_2 += 1;
                    else
                    {
                        anti_counter_index_bjet_method_2 += 1;
                        if(jetIndex_is_bquarkFromSMtop == -999) counter_genmatched_index_bjet_unphysical += 1;
                        else
                        {
                            anti_counter_genmatched_index_bjet_unphysical += 1;
                            //printf("[check-bjet] jetIndex_is_bquarkFromSMtop = %d\n", jetIndex_is_bquarkFromSMtop);
                            bool isIncluded = false;
                            vector<double> vec_delta_R;
                            for(std::size_t i=0; i!=indices_bjet.size(); ++i){
                                vec_delta_R.push_back(bquark.DeltaR(Jets[indices_bjet[i]]));
                                if(indices_bjet[i] == jetIndex_is_bquarkFromSMtop) isIncluded = true;
                            }

                            if(isIncluded)
                            {
                                counter_not_matched_but_is_in_the_list += 1;
                                
                                //for(std::size_t i=0; i!=indices_bjet.size(); ++i){
                                //    if(i==0) printf("[check-bjet] indices of loose b-jets = %d ", indices_bjet[i]);
                                //    else     printf(" %d", indices_bjet[i]);
                                //}
                                //printf("\n");

                                //for(std::size_t i=0; i!=indices_bjet.size(); ++i){
                                //    if(i==0) printf("[check-bjet] deltaR of loose b-jets = %7.4f ", vec_delta_R[i]);
                                //    else     printf(" %7.4f", vec_delta_R[i]);
                                //}
                                //printf("\n");

                                //for(std::size_t i=0; i!=indices_bjet.size(); ++i){
                                //    if(i==0) printf("[check-bjet] score of loose b-jets = %7.4f ", btag_scores[indices_bjet[i]]);
                                //    else     printf(" %7.4f", btag_scores[i]);
                                //}
                                //printf("\n");

                            }

                        }
                    }
                }
                //}}}
                if(index_bjet_method_2 == jetIndex_is_bquarkFromSMtop)
                    hist_entries_bjet_wp_matched  -> Fill(iwp);


                // four methods for top reconstruction(ST){{{
                //### chi-2 Study
                indices_bjet = {index_bjet_method_2}; // b-jet with the highest b-tag scores
                std::vector<int> index_jet_chi2_simple   = get_bjj_indices_min_chi2(Jets, indices_bjet, false); index_jet_chi2_simple.push_back(-1);
                std::vector<int> index_jet_chi2_modified = get_bjj_indices_min_chi2(Jets, indices_bjet, true);  index_jet_chi2_modified.push_back(-1);
                std::vector<int> index_jet_chi2_improved = {-1, -1, -1, -1}; // b j j q

                bool double_check = (index_jet_chi2_improved[0] != jetIndex_is_bquarkFromSMtop) && (index_bjet_method_2 == jetIndex_is_bquarkFromSMtop);
                if(double_check) counter_difference += 1;
                //}}}
                //// four methods for top reconstruction(TT){{{
                ////### chi-2 Study
                //indices_bjet = {index_bjet_method_2}; // b-jet with the highest b-tag scores
                //std::vector<int> index_jet_chi2_simple = get_bjjq_indices_min_chi2(Jets, indices_bjet, diphoton, false);
                //std::vector<int> index_jet_chi2_modified = get_bjjq_indices_min_chi2(Jets, indices_bjet, diphoton, true);
                //std::vector<int> index_jet_chi2_improved = get_bjjq_indices_min_chi2_3x3(Jets, indices_bjet, diphoton);

                //bool double_check = (index_jet_chi2_improved[0] != jetIndex_is_bquarkFromSMtop) && (index_bjet_method_2 == jetIndex_is_bquarkFromSMtop);
                //if(double_check) counter_difference += 1;
                ////}}}
                //when the higest b-tag score not work, how about other methods? (skipped){{{
                //if(index_bjet_method_2 != jetIndex_is_bquarkFromSMtop && jetIndex_is_bquarkFromSMtop != -999){
                //printf("[check-bjet] index_bjet (gen, sim, mod, imp, max_score) = (%d, ", index_jet_gen_matched[0]);
                //printf("%d, ", index_jet_chi2_simple[0]);
                //printf("%d, ", index_jet_chi2_modified[0]);
                //printf("%d, ", index_jet_chi2_improved[0]);
                //printf("%d)", index_bjet_method_2);
                //printf("\n");
                //}
                //}}}
                //printf(debug){{{
                //printf("[check-ywk] reco (gen-matching): (%d, %d, %d, %d)\n", jetIndex_is_bquarkFromSMtop, jetIndex_momPdgID_is_wboson[0], jetIndex_momPdgID_is_wboson[1], jetIndex_is_quarkFromFCNtop);
                //printf("[check-ywk] chi-2 sim          : (%d, %d, %d, %d)\n", index_jet_chi2_simple[0], index_jet_chi2_simple[1], index_jet_chi2_simple[2], index_jet_chi2_simple[3]);
                //printf("[check-ywk] chi-2 mod          : (%d, %d, %d, %d)\n", index_jet_chi2_modified[0], index_jet_chi2_modified[1], index_jet_chi2_modified[2], index_jet_chi2_modified[3]);
                //printf("[check-ywk] chi-2 imp          : (%d, %d, %d, %d)\n\n", index_jet_chi2_improved[0], index_jet_chi2_improved[1], index_jet_chi2_improved[2], index_jet_chi2_improved[3]);
                //}}}
                // Evaluation: counting for efficiency{{{
                std::vector<bool> matching_result;
    
                matching_result = get_matching_result(index_jet_chi2_simple, index_jet_gen_matched); //{tbwIsCorrectlyMatched, tqhIsCorrectlyMatched, sigIsCorrectlyReco};
                if(matching_result[0]) counter_tbwIsCorrectlyMatched_simple[iwp][ibool] += 1;
                if(matching_result[1]) counter_tqhIsCorrectlyMatched_simple[iwp][ibool] += 1;
                if(matching_result[2]) counter_sigIsCorrectlyReco_simple[iwp][ibool]    += 1;
                if(matching_result[3]) counter_bjet_is_matched_simple[iwp][ibool]       += 1;
    
                matching_result = get_matching_result(index_jet_chi2_modified, index_jet_gen_matched); //{tbwIsCorrectlyMatched, tqhIsCorrectlyMatched, sigIsCorrectlyReco};
                if(matching_result[0]) counter_tbwIsCorrectlyMatched_modified[iwp][ibool] += 1;
                if(matching_result[1]) counter_tqhIsCorrectlyMatched_modified[iwp][ibool] += 1;
                if(matching_result[2]) counter_sigIsCorrectlyReco_modified[iwp][ibool]    += 1;
                if(matching_result[3]) counter_bjet_is_matched_modified[iwp][ibool]       += 1;
    
                matching_result = get_matching_result(index_jet_chi2_improved, index_jet_gen_matched); //{tbwIsCorrectlyMatched, tqhIsCorrectlyMatched, sigIsCorrectlyReco};
                if(matching_result[0]) counter_tbwIsCorrectlyMatched_improved[iwp][ibool] += 1;
                if(matching_result[1]) counter_tqhIsCorrectlyMatched_improved[iwp][ibool] += 1;
                if(matching_result[2]) counter_sigIsCorrectlyReco_improved[iwp][ibool]    += 1;
                if(matching_result[3]) counter_bjet_is_matched_improved[iwp][ibool]       += 1;
                //}}}
                // Evaluation: hist for bjet matching efficiency{{{
                if( index_jet_chi2_simple[0]   == jetIndex_is_bquarkFromSMtop && ibool == 0 ){ hist_rate_permuting_bjet_is_bquark_sim -> Fill(iwp); }
                if( index_jet_chi2_modified[0] == jetIndex_is_bquarkFromSMtop && ibool == 0 ){ hist_rate_permuting_bjet_is_bquark_mod -> Fill(iwp); }
                if( index_jet_chi2_improved[0] == jetIndex_is_bquarkFromSMtop && ibool == 0 ){ hist_rate_permuting_bjet_is_bquark_imp -> Fill(iwp); }
                //}}}
                Nevents_pass_selection[iwp][ibool] += 1;
            } // end of wp for loop
        } // end of mulplicity condition for loop
        //}}}
        
        /* skip tmp{{{
        // Define variables / interface{{{
        bool isMatched, is_tqh_quark=false;//boolean used to check matching
        bool isMatched_tbw_simple, is_tqh_quark_simple=false;//boolean used to check matching
        bool isMatched_tbw_modified, is_tqh_quark_modified=false;//boolean used to check matching
        bool isMatched_tbw_improved, is_tqh_quark_improved=false;//boolean used to check matching
        //-------------------------//
        int index_q_simple = -999;

        std::vector<int> id_jet_chi2_simple(2);
        id_jet_chi2_simple[0] = Jets_GenFalavor[index_jet_chi2_simple[1]];
        id_jet_chi2_simple[1] = Jets_GenFalavor[index_jet_chi2_simple[2]];

        std::vector<TLorentzVector> jet_chi2_simple(2);
        jet_chi2_simple[0] = Jets[index_jet_chi2_simple[1]];
        jet_chi2_simple[1] = Jets[index_jet_chi2_simple[2]];

        std::vector<TLorentzVector> genParticle_jet_chi2_simple(2);
        //-------------------------//
        int index_q_modified = -999;

        std::vector<int> id_jet_chi2_modified(2);
        id_jet_chi2_modified[0] = Jets_GenFalavor[index_jet_chi2_modified[1]];
        id_jet_chi2_modified[1] = Jets_GenFalavor[index_jet_chi2_modified[2]];

        std::vector<TLorentzVector> jet_chi2_modified(2);
        jet_chi2_modified[0] = Jets[index_jet_chi2_modified[1]];
        jet_chi2_modified[1] = Jets[index_jet_chi2_modified[2]];

        std::vector<TLorentzVector> genParticle_jet_chi2_modified(2);
        //-------------------------//
        //}}}
        // GenMatching{{{
        std::vector<int> index_GenParticles;
        int id_genParticle_bjet=-999;
        int id_genParticle_jet_chi2_simple[2]; std::fill_n(id_genParticle_jet_chi2_simple, 2, -999);
        int id_genParticle_jet_chi2_modified[2]; std::fill_n(id_genParticle_jet_chi2_modified, 2, -999);
        double delta_R = 0;
        genParticle_bjet = GetGenParticle(bjet, treeReader.GenPartInfo_size, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass, treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID, index_GenParticles, id_genParticle_bjet);

        if( abs(id_genParticle_bjet)==5 && abs(treeReader.GenPartInfo_MomPdgID->at(index_GenParticles[0]))==6 ) count_bjet_is_bquark+=1;
        genParticle_jet_chi2_simple[0] = GetGenParticle(jet_chi2_simple[0], treeReader.GenPartInfo_size, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass, treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID, index_GenParticles, id_genParticle_jet_chi2_simple[0]);
        genParticle_jet_chi2_simple[1] = GetGenParticle(jet_chi2_simple[1], treeReader.GenPartInfo_size, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass, treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID, index_GenParticles, id_genParticle_jet_chi2_simple[1]);
        genParticle_jet_chi2_modified[0] = GetGenParticle(jet_chi2_modified[0], treeReader.GenPartInfo_size, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass, treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID, index_GenParticles, id_genParticle_jet_chi2_modified[0]);
        genParticle_jet_chi2_modified[1] = GetGenParticle(jet_chi2_modified[1], treeReader.GenPartInfo_size, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass, treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID, index_GenParticles, id_genParticle_jet_chi2_modified[1]);
        // }}}
        // Reconstruct Mass (M2){{{
        //----- reco -----//
        TLorentzVector w_candidate_chi2_simple = jet_chi2_simple[0] + jet_chi2_simple[1];
        TLorentzVector top_candidate_chi2_simple = w_candidate_chi2_simple + bjet;
        mytree.Mass_w_candidate_chi2_simple = w_candidate_chi2_simple.M();
        mytree.Mass_top_candidate_chi2_simple = top_candidate_chi2_simple.M();
        hist_mass_w_candidate_chi2_simple->Fill(w_candidate_chi2_simple.M());
        hist_mass_top_candidate_chi2_simple->Fill(top_candidate_chi2_simple.M());
        //--------------------
        TLorentzVector w_candidate_chi2_modified = jet_chi2_modified[0] + jet_chi2_modified[1];
        TLorentzVector top_candidate_chi2_modified = w_candidate_chi2_modified + bjet;
        mytree.Mass_w_candidate_chi2_modified = w_candidate_chi2_modified.M();
        mytree.Mass_top_candidate_chi2_modified = top_candidate_chi2_modified.M();
        hist_mass_w_candidate_chi2_modified->Fill(w_candidate_chi2_modified.M());
        hist_mass_top_candidate_chi2_modified->Fill(top_candidate_chi2_modified.M());
        //----- gen -----//
        //reject candidates according to the MC truth
        //1. reject candidate with id_genParticle = -999 (soft jet, not from hard process)
        //2. reject candidate with id_genParticle_jet == 5 (w_boson cannot decay to b jet)
        //3. reject candidate with 2 jets having same genParticle
        TLorentzVector gen_w_candidate_chi2_simple = genParticle_jet_chi2_simple[0] + genParticle_jet_chi2_simple[1];
        TLorentzVector gen_top_candidate_chi2_simple = gen_w_candidate_chi2_simple + genParticle_bjet;
        hist_mass_gen_w_candidate_chi2_simple->Fill(gen_w_candidate_chi2_simple.M());
        hist_mass_gen_top_candidate_chi2_simple->Fill(gen_top_candidate_chi2_simple.M());
        //----- check -----//
        if(gen_w_candidate_chi2_simple.M() < 20){
            check_M20_simple += 1;
            if( !(id_genParticle_jet_chi2_simple[0]==-999 || id_genParticle_jet_chi2_simple[1]==-999) ){
                check_gen_exclude_id_999 += 1;
                if( !(genParticle_jet_chi2_simple[0].Pt() == genParticle_jet_chi2_simple[1].Pt()) ){
                    check_gen_exclude_the_same += 1;
                    //kinematics_report("jet[0]", jet_chi2_simple[0], id_jet_chi2_simple[0], "gen[0]", genParticle_jet_chi2_simple[0], id_genParticle_jet_chi2_simple[0]);
                    //kinematics_report("jet[1]", jet_chi2_simple[1], id_jet_chi2_simple[1], "gen[1]", genParticle_jet_chi2_simple[1], id_genParticle_jet_chi2_simple[1]);
                    //printf("\n");
                }
            }
        }

        isMatched_tbw_simple = isMatched_with_Gen_tbw(treeReader.GenPartInfo_PdgID, treeReader.GenPartInfo_MomPdgID, index_GenParticles[0], index_GenParticles[1], index_GenParticles[2] );//bjet, jet12(simple)
        if( isMatched_tbw_simple ){
            mytree.Mass_gen_w_candidate_chi2_simple = gen_w_candidate_chi2_simple.M();
            mytree.Mass_gen_top_candidate_chi2_simple = gen_top_candidate_chi2_simple.M();
            hist_mass_conditioned_gen_w_candidate_chi2_simple->Fill(gen_w_candidate_chi2_simple.M());
            hist_mass_conditioned_gen_top_candidate_chi2_simple->Fill(gen_top_candidate_chi2_simple.M());
            //old test code{{{
            if(gen_w_candidate_chi2_simple.M() < 20){// debug purpose
                //printf("\nCheck:\n");
                //kinematics_report("bjet", bjet, id_bjet, "genP", genParticle_bjet, id_genParticle_bjet);
                //kinematics_report("jet_chi2_simple[0]", jet_chi2_simple[0], id_jet_chi2_simple[0], "MatchedGenParticle", genParticle_jet_chi2_simple[0], id_genParticle_jet_chi2_simple[0]);
                //kinematics_report("jet_chi2_simple[1]", jet_chi2_simple[1], id_jet_chi2_simple[1], "MatchedGenParticle", genParticle_jet_chi2_simple[1], id_genParticle_jet_chi2_simple[1]);
                //printf("[INFO] M_gen_jj  = %6.2f\n", gen_w_candidate_chi2_simple.M());
                //printf("[INFO] M_gen_bjj = %6.2f\n", gen_top_candidate_chi2_simple.M());
            }
            //}}}
        }
        //--------------------
        TLorentzVector gen_w_candidate_chi2_modified = genParticle_jet_chi2_modified[0] + genParticle_jet_chi2_modified[1];
        TLorentzVector gen_top_candidate_chi2_modified = gen_w_candidate_chi2_modified + genParticle_bjet;
        hist_mass_gen_w_candidate_chi2_modified->Fill(gen_w_candidate_chi2_modified.M());
        hist_mass_gen_top_candidate_chi2_modified->Fill(gen_top_candidate_chi2_modified.M());

        isMatched_tbw_modified = isMatched_with_Gen_tbw(treeReader.GenPartInfo_PdgID, treeReader.GenPartInfo_MomPdgID, index_GenParticles[0], index_GenParticles[3], index_GenParticles[4] );//bjet, jet12(modified)
        if( isMatched_tbw_modified ){
            mytree.Mass_gen_w_candidate_chi2_modified = gen_w_candidate_chi2_modified.M();
            mytree.Mass_gen_top_candidate_chi2_modified = gen_top_candidate_chi2_modified.M();
            hist_mass_conditioned_gen_w_candidate_chi2_modified->Fill(gen_w_candidate_chi2_modified.M());
            hist_mass_conditioned_gen_top_candidate_chi2_modified->Fill(gen_top_candidate_chi2_modified.M());
        }
        //--------------------
        // }}}
        // Reconstruct Mass (M1){{{
        double M1_simple = -999;
        double M1_modified = -999;
        TLorentzVector top_fcnh_simple, top_fcnh_modified, jet_q_simple, jet_q_modified;

        //--- simple chi2 ---//
        top_fcnh_simple = GetBestM1(M1_simple, mytree.num_jets, index_bjet, index_jet_chi2_simple, diphoton, Jets, index_q_simple, jet_q_simple);
        if(M1_simple != -999 && M1_simple<20) check_M1_20+=1;
        hist_mass_top_fcnh_simple->Fill(M1_simple);
        is_tqh_quark_simple = is_this_tqh_quark(jet_q_simple, treeReader.GenPartInfo_size,\
                          treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                          treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
        
        //--- modified chi2 ---//
        top_fcnh_modified = GetBestM1(M1_modified, mytree.num_jets, index_bjet, index_jet_chi2_modified, diphoton, Jets, index_q_modified, jet_q_modified);
        hist_mass_top_fcnh_modified->Fill(M1_modified);
        is_tqh_quark_modified = is_this_tqh_quark(jet_q_modified, treeReader.GenPartInfo_size,\
                            treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                            treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
        //}}}
        //### Yueh-Feng's method{{{
        // combinations from leading jets{{{
        std::vector<int> index_leading_jets(3);// store the correct indices in meta-indices
        std::vector<TLorentzVector> leading_jets(3);
        // use meta-indices to store the leading-jets info
        if(index_bjet == 0){
            leading_jets[0] = Jets[1]; index_leading_jets[0] = 1;
            leading_jets[1] = Jets[2]; index_leading_jets[1] = 2;
            leading_jets[2] = Jets[3]; index_leading_jets[2] = 3;
        } else if(index_bjet == 1){
            leading_jets[0] = Jets[0]; index_leading_jets[0] = 0;
            leading_jets[1] = Jets[2]; index_leading_jets[1] = 2;
            leading_jets[2] = Jets[3]; index_leading_jets[2] = 3;
        } else if(index_bjet == 2){
            leading_jets[0] = Jets[0]; index_leading_jets[0] = 0;
            leading_jets[1] = Jets[1]; index_leading_jets[1] = 1;
            leading_jets[2] = Jets[3]; index_leading_jets[2] = 3;
        } else{
            leading_jets[0] = Jets[0]; index_leading_jets[0] = 0;
            leading_jets[1] = Jets[1]; index_leading_jets[1] = 1;
            leading_jets[2] = Jets[2]; index_leading_jets[2] = 2;
        }
        // the meta-index of the jet from tqh coupling is tracked
        std::vector<TLorentzVector> ggj(3);
        std::vector<TLorentzVector> wjj(3);
        std::vector<TLorentzVector> bjj(3);
        for(int i=0; i<3; ++i) ggj[i] = diphoton + leading_jets[i];
        wjj[0] = leading_jets[1] + leading_jets[2];
        wjj[1] = leading_jets[0] + leading_jets[2];
        wjj[2] = leading_jets[0] + leading_jets[1];
        bjj[0] = wjj[0] + bjet;
        bjj[1] = wjj[1] + bjet;
        bjj[2] = wjj[2] + bjet;
        //}}}
        // invariant mass of best candidate{{{
        std::vector<double> invm_w(3);
        std::vector<double> invm_ggj(3);
        std::vector<double> invm_bjj(3);
        std::vector<double> t1t2_angle(3);
        std::vector<double> M1M2_criteria(3);

        for(int i=0; i<3; ++i){
            invm_w[i] = wjj[i].M();
            invm_ggj[i] = ggj[i].M();
            invm_bjj[i] = bjj[i].M();
            t1t2_angle[i] = ggj[i].Angle(bjj[i].Vect());
            M1M2_criteria[i] = GetM1M2_ratio(invm_ggj[i], invm_bjj[i]);
        }

        int minimum_index = -999;
        double minimum_mass_ratio = 9999;
        for(int i=0; i<3; ++i){
            if(M1M2_criteria[i]<minimum_mass_ratio){
                minimum_mass_ratio = M1M2_criteria[i];
                minimum_index = i;
            }
        }

        hist_mass_w_candidate_yfyj->Fill(invm_w[minimum_index]);
        hist_mass_t1_candidate_yfyj->Fill(invm_ggj[minimum_index]);
        hist_mass_t2_candidate_yfyj->Fill(invm_bjj[minimum_index]);
        //}}}
        // decode the meta-indices{{{
        int index_q_yfyj = -999;
        std::vector<int> index_jet_yfyj(2, -999);
        // decode the meta-indices to correct indices of jets
        if     (minimum_index == 0) {
            index_q_yfyj = index_leading_jets[0];
            index_jet_yfyj[0] = index_leading_jets[1];
            index_jet_yfyj[1] = index_leading_jets[2];
        }
        else if(minimum_index == 1) {
            index_q_yfyj = index_leading_jets[1];
            index_jet_yfyj[0] = index_leading_jets[0];
            index_jet_yfyj[1] = index_leading_jets[2];
        }
        else                        {
            index_q_yfyj = index_leading_jets[2];
            index_jet_yfyj[0] = index_leading_jets[0];
            index_jet_yfyj[1] = index_leading_jets[1];
        }
        //}}}
        //Estimate accuracy{{{
        //reject candidates according to the MC truth
        //1. reject candidate with id_genParticle = -999 (soft jet, not from hard process)
        //2. reject candidate with id_genParticle_jet == 5 (w_boson cannot decay to b jet)
        //3. reject candidate with 2 jets having same genParticle
        int id_gen_wjets_yfyj[2];
        TLorentzVector leading_jets_genParticles[2];
        leading_jets_genParticles[0] = GetGenParticle(Jets[index_jet_yfyj[0]], treeReader.GenPartInfo_size, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass, treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID, index_GenParticles, id_gen_wjets_yfyj[0]);
        leading_jets_genParticles[1] = GetGenParticle(Jets[index_jet_yfyj[1]], treeReader.GenPartInfo_size, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass, treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID, index_GenParticles, id_gen_wjets_yfyj[1]);

        TLorentzVector gen_w_candidate_yfyj = leading_jets_genParticles[0] + leading_jets_genParticles[1];
        TLorentzVector gen_top_candidate_yfyj = gen_w_candidate_yfyj + genParticle_bjet;
        hist_mass_gen_w_candidate_yfyj->Fill(gen_w_candidate_yfyj.M());
        hist_mass_gen_t2_candidate_yfyj->Fill(gen_top_candidate_yfyj.M());

        isMatched = isMatched_with_Gen_tbw(treeReader.GenPartInfo_PdgID, treeReader.GenPartInfo_MomPdgID, index_GenParticles[0], index_GenParticles[5], index_GenParticles[6] );//bjet, jet12(yfyj)
        if( isMatched ){
            hist_mass_conditioned_gen_w_candidate_yfyj->Fill(gen_w_candidate_yfyj.M());
            hist_mass_conditioned_gen_t2_candidate_yfyj->Fill(gen_top_candidate_yfyj.M());
        }

        is_tqh_quark = is_this_tqh_quark(leading_jets[minimum_index], treeReader.GenPartInfo_size,\
                            treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                            treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);

        //}}}
        //}}}
        //### Improved chi-2 method{{{
        // 1. pick up 3 jets
        // 2. q-jj three combinations 
        std::vector<int> index_wjets(2);
        int index_tqh_qjet;
        TLorentzVector wjets[2], tqh_qjet;
        TLorentzVector w_boson, sm_top, fcnc_top; 
        chi2_min = 99999;
        for(int i=0; i<mytree.num_jets; ++i){
            if(i==index_bjet) continue;//bypass bjet
            for(int j=i+1; j<mytree.num_jets; ++j){
                if(j==index_bjet) continue;//bypass bjet
                for(int k=j+1; k<mytree.num_jets; ++k){
                    if(k==index_bjet) continue;//bypass bjet
                    //--- combinations{{{
                    TLorentzVector jets_chosen[3];
                    jets_chosen[0] = Jets[i];
                    jets_chosen[1] = Jets[j];
                    jets_chosen[2] = Jets[k];

                    TLorentzVector tqh_q_chosen[3];
                    tqh_q_chosen[0] = jets_chosen[0];
                    tqh_q_chosen[1] = jets_chosen[1];
                    tqh_q_chosen[2] = jets_chosen[2];

                    TLorentzVector w_candidate[3];
                    w_candidate[0] = jets_chosen[1] + jets_chosen[2];
                    w_candidate[1] = jets_chosen[0] + jets_chosen[2];
                    w_candidate[2] = jets_chosen[0] + jets_chosen[1];

                    TLorentzVector top_candidate[3];
                    TLorentzVector fcnc_top_candidate[3];
                    //}}}
                    //--- calculation{{{
                    std::vector<double> chi2;
                    double w_mass[3], t_mass[3], fcnc_top_mass[3];
                    for(int x = 0; x<3; ++x){
                        w_mass[x] =  w_candidate[x].M();

                        top_candidate[x] = w_candidate[x] + bjet;
                        t_mass[x] =  w_candidate[x].M();

                        fcnc_top_candidate[x] = tqh_q_chosen[x] + diphoton;
                        fcnc_top_mass[x] = fcnc_top_candidate[x].M();

                        chi2.push_back( Chi2_calculator_improved(w_mass[x], t_mass[x], fcnc_top_mass[x]) );
                    }
                    //}}}
                    //--- sorting{{{
                    int smallest_chi2_index = std::min_element(chi2.begin(),chi2.end()) - chi2.begin();
                    double smallest_chi2 = *std::min_element(chi2.begin(),chi2.end());

                    if(smallest_chi2 < chi2_min){
                        if(smallest_chi2_index == 0){
                            tqh_qjet = tqh_q_chosen[0];
                            wjets[0] = jets_chosen[1];
                            wjets[1] = jets_chosen[2];
                            index_tqh_qjet = i;
                            index_wjets[0] = j;
                            index_wjets[1] = k;
                        } else if(smallest_chi2_index == 1){
                            tqh_qjet = tqh_q_chosen[1];
                            wjets[0] = jets_chosen[0];
                            wjets[1] = jets_chosen[2];
                            index_tqh_qjet = j;
                            index_wjets[0] = i;
                            index_wjets[1] = k;
                        } else{
                            tqh_qjet = tqh_q_chosen[2];
                            wjets[0] = jets_chosen[0];
                            wjets[1] = jets_chosen[1];
                            index_tqh_qjet = k;
                            index_wjets[0] = i;
                            index_wjets[1] = j;
                        }

                        chi2_min = smallest_chi2;
                    }
                    //}}}
                }
            }
        }//end of looping jets
        w_boson = wjets[0] + wjets[1];
        sm_top = w_boson + bjet;
        fcnc_top = tqh_qjet + diphoton;


        int id_gen_qjet;
        int id_gen_wjets[2];
        //index_GenParticles.clear();
        TLorentzVector gen_wjets[2], gen_tqh_qjet;

        gen_wjets[0] = GetGenParticle(wjets[0], treeReader.GenPartInfo_size, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass, treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID, index_GenParticles, id_gen_wjets[0]);
        gen_wjets[1] = GetGenParticle(wjets[1], treeReader.GenPartInfo_size, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass, treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID, index_GenParticles, id_gen_wjets[1]);
        gen_tqh_qjet = GetGenParticle(tqh_qjet, treeReader.GenPartInfo_size, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass, treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID, index_GenParticles, id_gen_qjet);
        //---
        isMatched_tbw_improved = isMatched_with_Gen_tbw(treeReader.GenPartInfo_PdgID, treeReader.GenPartInfo_MomPdgID, index_GenParticles[0], index_GenParticles[7], index_GenParticles[8] );
        is_tqh_quark_improved = is_this_tqh_quark(tqh_qjet, treeReader.GenPartInfo_size,\
                            treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
                            treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
        bool isFinalMatched = isMatched_tbw_improved && is_tqh_quark_improved;

        //if(isFinalMatched){
        //    printf("[Check-mod] jets index = %d, %d\n", index_GenParticles[3], index_GenParticles[4]);
        //    printf("[Check-imp] jets index = %d, %d\n\n", index_GenParticles[7], index_GenParticles[8]);
        //}
        //printf("[Check-mod] jets index = %d, %d\n", index_jet_chi2_modified[0], index_jet_chi2_modified[1]);
        //printf("[Check-imp] jets index = %d, %d, %d\n", index_wjets[0], index_wjets[1], index_tqh_qjet);

        //}}}
        
        //fill hist_factories{{{
        hf_wboson.Fill_hist( w_candidate_chi2_modified, wboson );
        hf_tbw   .Fill_hist( top_candidate_chi2_modified, topquark_tbw );
        hf_tqh   .Fill_hist( top_fcnh_modified, topquark_tqh );

        if(jetIndex_is_quarkFromFCNtop != -999){
            TLorentzVector tqh_matched = diphoton + Jets[jetIndex_is_quarkFromFCNtop];
            hf_tqh_matched.Fill_hist(tqh_matched, topquark_tqh);
        }
        //}}}
        // Evaluation{{{
        //### print indices (skipped){{{
        //if(sigCanBeReconstructed){
        //    int _index_0 = jetIndex_momPdgID_is_wboson[0];
        //    int _index_1 = jetIndex_momPdgID_is_wboson[1];
        //    printf("[INFO] indices (q, b, j, j) correct  = %2d, %2d, %2d, %2d\n", jetIndex_is_quarkFromFCNtop, jetIndex_is_bquarkFromSMtop, _index_0, _index_1);
        //    printf("[INFO] indices (q, b, j, j) yfyj     = %2d, %2d, %2d, %2d\n", index_q_yfyj, index_bjet, index_jet_yfyj[0], index_jet_yfyj[1]);
        //    printf("[INFO] indices (q, b, j, j) simple   = %2d, %2d, %2d, %2d\n", index_q_simple, index_bjet, index_jet_chi2_simple[0], index_jet_chi2_simple[1]);
        //    printf("[INFO] indices (q, b, j, j) modified = %2d, %2d, %2d, %2d\n", index_q_modified, index_bjet, index_jet_chi2_modified[0], index_jet_chi2_modified[1]);
        //    printf("[INFO] indices (q, b, j, j) improved = %2d, %2d, %2d, %2d\n", index_tqh_qjet, index_bjet, index_wjets[0], index_wjets[1]);
        //}
        //}}}
        //### check if indices match
        bool tbwIsCorrectlyMatched;
        if(tbwCanBeReconstructed){
            //bool are_From_W_boson = (  ((index_jet_chi2_modified[0] == jetIndex_momPdgID_is_wboson[0]) && (index_jet_chi2_modified[1] == jetIndex_momPdgID_is_wboson[1]))\
            //                        || ((index_jet_chi2_modified[0] == jetIndex_momPdgID_is_wboson[1]) && (index_jet_chi2_modified[1] == jetIndex_momPdgID_is_wboson[0])) );
            //bool tbwIsCorrectlyMatched = (index_bjet == jetIndex_is_bquarkFromSMtop) && are_From_W_boson;
            tbwIsCorrectlyMatched = check_tbwIsCorrectlyMatched(index_jet_yfyj, jetIndex_momPdgID_is_wboson, index_bjet, jetIndex_is_bquarkFromSMtop);
            if(tbwIsCorrectlyMatched) counter_tbwIsCorrectlyMatched_yfyj += 1;
            tbwIsCorrectlyMatched = check_tbwIsCorrectlyMatched(index_jet_chi2_simple, jetIndex_momPdgID_is_wboson, index_bjet, jetIndex_is_bquarkFromSMtop);
            if(tbwIsCorrectlyMatched) counter_tbwIsCorrectlyMatched_simple += 1;
            tbwIsCorrectlyMatched = check_tbwIsCorrectlyMatched(index_jet_chi2_modified, jetIndex_momPdgID_is_wboson, index_bjet, jetIndex_is_bquarkFromSMtop);
            if(tbwIsCorrectlyMatched) counter_tbwIsCorrectlyMatched_modified += 1;
            tbwIsCorrectlyMatched = check_tbwIsCorrectlyMatched(index_wjets, jetIndex_momPdgID_is_wboson, index_bjet, jetIndex_is_bquarkFromSMtop);
            if(tbwIsCorrectlyMatched) counter_tbwIsCorrectlyMatched_improved += 1;
        }

        bool tqhIsCorrectlyMatched;
        if(tqhCanBeReconstructed){
            //bool tqhIsCorrectlyMatched = (index_qjet == jetIndex_is_quarkFromFCNtop);
            tqhIsCorrectlyMatched = check_tqhIsCorrectlyMatched(index_q_yfyj, jetIndex_is_quarkFromFCNtop);
            if(tqhIsCorrectlyMatched) counter_tqhIsCorrectlyMatched_yfyj += 1;
            tqhIsCorrectlyMatched = check_tqhIsCorrectlyMatched(index_q_simple, jetIndex_is_quarkFromFCNtop);
            if(tqhIsCorrectlyMatched) counter_tqhIsCorrectlyMatched_simple += 1;
            tqhIsCorrectlyMatched = check_tqhIsCorrectlyMatched(index_q_modified, jetIndex_is_quarkFromFCNtop);
            if(tqhIsCorrectlyMatched) counter_tqhIsCorrectlyMatched_modified += 1;
            tqhIsCorrectlyMatched = check_tqhIsCorrectlyMatched(index_tqh_qjet, jetIndex_is_quarkFromFCNtop);
            if(tqhIsCorrectlyMatched) counter_tqhIsCorrectlyMatched_improved += 1;
        }

        bool sigIsCorrectlyReco;
        if(tbwCanBeReconstructed && tqhCanBeReconstructed){
            tbwIsCorrectlyMatched = check_tbwIsCorrectlyMatched(index_jet_yfyj, jetIndex_momPdgID_is_wboson, index_bjet, jetIndex_is_bquarkFromSMtop);
            tqhIsCorrectlyMatched = check_tqhIsCorrectlyMatched(index_q_yfyj, jetIndex_is_quarkFromFCNtop);
            sigIsCorrectlyReco = tbwIsCorrectlyMatched && tqhIsCorrectlyMatched;
            if(sigIsCorrectlyReco) counter_sigIsCorrectlyReco_yfyj += 1;

            tbwIsCorrectlyMatched = check_tbwIsCorrectlyMatched(index_jet_chi2_simple, jetIndex_momPdgID_is_wboson, index_bjet, jetIndex_is_bquarkFromSMtop);
            tqhIsCorrectlyMatched = check_tqhIsCorrectlyMatched(index_q_simple, jetIndex_is_quarkFromFCNtop);
            sigIsCorrectlyReco = tbwIsCorrectlyMatched && tqhIsCorrectlyMatched;
            if(sigIsCorrectlyReco) counter_sigIsCorrectlyReco_simple += 1;

            tbwIsCorrectlyMatched = check_tbwIsCorrectlyMatched(index_jet_chi2_modified, jetIndex_momPdgID_is_wboson, index_bjet, jetIndex_is_bquarkFromSMtop);
            tqhIsCorrectlyMatched = check_tqhIsCorrectlyMatched(index_q_modified, jetIndex_is_quarkFromFCNtop);
            sigIsCorrectlyReco = tbwIsCorrectlyMatched && tqhIsCorrectlyMatched;
            if(sigIsCorrectlyReco) counter_sigIsCorrectlyReco_modified += 1;

            tbwIsCorrectlyMatched = check_tbwIsCorrectlyMatched(index_wjets, jetIndex_momPdgID_is_wboson, index_bjet, jetIndex_is_bquarkFromSMtop);
            tqhIsCorrectlyMatched = check_tqhIsCorrectlyMatched(index_tqh_qjet, jetIndex_is_quarkFromFCNtop);
            sigIsCorrectlyReco = tbwIsCorrectlyMatched && tqhIsCorrectlyMatched;
            if(sigIsCorrectlyReco) counter_sigIsCorrectlyReco_improved += 1;
        }
        //}}}
        //### Details of selected quark from tqh (skipped){{{
        //printf("--------------------------------------------------\n");
        //printf("[debug-tqh] index (true, sim, mod) = (%d, %d, %d) \n", jetIndex_is_quarkFromFCNtop, index_q_simple, index_q_modified);
        //printf("--------------------------------------------------\n");
        //printf("[check] tqh kinematics\n");
        //printf("[sim] ");
        //kinematics_info("", jet_q_simple, index_q_simple);
        //printf("[mod] ");
        //kinematics_info("", jet_q_modified, index_q_modified);
        //printf("[yjy] ");
        //kinematics_info("", leading_jets[minimum_index], index_q_yfyj);
        //printf("[imp] ");
        //kinematics_info("", tqh_qjet, index_tqh_qjet);
        //printf("--------------------------------------------------\n");
        //}}}
        //### Check GenInfo (skipped){{{
        ////if(num_bquark>1){
        //printf("\n");
        //printf("[GenCheck] GenPartInfo_size = %d\n", treeReader.GenPartInfo_size);
        //for(int i=0; i<treeReader.GenPartInfo_size; i++){
        //    int pdgID = treeReader.GenPartInfo_PdgID->at(i);
        //    bool isPromptFinalState = treeReader.GenPartInfo_isPromptFinalState->at(i);
        //    bool isNeutrino = (abs(pdgID) == 12 || abs(pdgID) == 14 || abs(pdgID) == 16) && isPromptFinalState;
        //    bool isChargedLepton = (abs(pdgID) == 11 || abs(pdgID) == 13 || abs(pdgID) == 15) && isPromptFinalState;
        //    bool isWboson = (abs(pdgID) == 24);
        //    bool isTop = (abs(pdgID) == 6);

        //    //if(isNeutrino){
        //    //if(isChargedLepton){
        //    //if(isWboson){
        //    if(isTop){
        //    printf("Status = %3d, ", treeReader.GenPartInfo_Status->at(i));
        //    printf("PdgID = %3d, ", treeReader.GenPartInfo_PdgID->at(i));
        //    printf("Pt = %6.2f, ", treeReader.GenPartInfo_Pt->at(i));
        //    printf("Eta = %9.2f, ", treeReader.GenPartInfo_Eta->at(i));
        //    printf("Phi = %6.2f, ", treeReader.GenPartInfo_Phi->at(i));
        //    printf("Mass = %6.2f, ", treeReader.GenPartInfo_Mass->at(i));
        //    printf("isHardProcess = %3d, ", treeReader.GenPartInfo_isHardProcess->at(i) ? 1 : 0);
        //    printf("isPromptFinalState = %3d, ", treeReader.GenPartInfo_isPromptFinalState->at(i) ? 1 : 0);
        //    printf("MomPdgID = %5d, ", treeReader.GenPartInfo_MomPdgID->at(i));
        //    printf("MomStatus = %3d\n", treeReader.GenPartInfo_MomStatus->at(i));
        //    }
        //}
        ////}
        //}}}
        //### Case study with Venn diagram{{{
        //tbwIsCorrectlyMatched = check_tbwIsCorrectlyMatched(index_jet_yfyj, jetIndex_momPdgID_is_wboson, index_bjet, jetIndex_is_bquarkFromSMtop);
        //tqhIsCorrectlyMatched = check_tqhIsCorrectlyMatched(index_q_yfyj, jetIndex_is_quarkFromFCNtop);
        tqhIsCorrectlyMatched = check_tqhIsCorrectlyMatched(index_q_simple, jetIndex_is_quarkFromFCNtop);
        //------------------------------
        if( is_tqh_quark_simple &&  tqhIsCorrectlyMatched) counter_intersection += 1;
        if( is_tqh_quark_simple && !tqhIsCorrectlyMatched) counter_M_notC += 1;
        if(!is_tqh_quark_simple &&  tqhIsCorrectlyMatched) counter_notM_C += 1;
        if(!is_tqh_quark_simple && !tqhIsCorrectlyMatched) counter_notM_notC += 1;

        //### TO BE INVESTIGATE!!! ###//
//        if(is_tqh_quark_simple && !tqhIsCorrectlyMatched){
//            //if( w_quark_is_matched_twice ) counter_mismatched += 1;
//            counter_mismatched += 1;
//            printf("\nientry = %d/%d\n", ientry+1, nentries+1);
//            printf("[CHECK-debug] counter_is_bquarkFromSMtop = %d\n", counter_is_bquarkFromSMtop);
//            printf("[CHECK-debug] num_jets = %d\n", mytree.num_jets);
//            for(int i=0; i<mytree.num_jets; ++i){
//                int index_gen, Matched_PdgID, Matched_MomPdgID;
//                kinematics_info("jet", Jets[i], i);
//                obtain_gen_matched_ID(true, Jets[i], index_gen, Matched_PdgID, Matched_MomPdgID, treeReader.GenPartInfo_size,\
//                                       treeReader.GenPartInfo_MomPdgID, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass,\
//                                       treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID);
//            }
//            int _index_0 = jetIndex_momPdgID_is_wboson[0];
//            int _index_1 = jetIndex_momPdgID_is_wboson[1];
//            printf("[INFO] indices (q, b, j, j) correct  = %2d, %2d, %2d, %2d\n", jetIndex_is_quarkFromFCNtop, jetIndex_is_bquarkFromSMtop, _index_0, _index_1);
//            printf("[INFO] indices (q, b, j, j) yfyj     = %2d, %2d, %2d, %2d\n", index_q_yfyj, index_bjet, index_jet_yfyj[0], index_jet_yfyj[1]);
//            printf("[INFO] indices (q, b, j, j) simple   = %2d, %2d, %2d, %2d\n", index_q_simple, index_bjet, index_jet_chi2_simple[0], index_jet_chi2_simple[1]);
//            printf("[INFO] indices (q, b, j, j) modified = %2d, %2d, %2d, %2d\n", index_q_modified, index_bjet, index_jet_chi2_modified[0], index_jet_chi2_modified[1]);
//            printf("[INFO] indices (q, b, j, j) improved = %2d, %2d, %2d, %2d\n", index_tqh_qjet, index_bjet, index_wjets[0], index_wjets[1]);
//            printf("[check] chosen correct index of bjet: %d (counter = %d)\n", jetIndex_is_bquarkFromSMtop, counter_is_bquarkFromSMtop);
//            printf("[check]");
//            printf("(%d) deltaR = %5.3f, " , jetIndex_is_bquarkFromSMtop_moreThanOne[0], deltaR_jet_bquark[0]);
//            printf("(%d) deltaR = %5.3f\n", jetIndex_is_bquarkFromSMtop_moreThanOne[1], deltaR_jet_bquark[1]);
//        }
        //}}}
        //Accuracy counting{{{
        if( isMatched_tbw_simple ) accuracy_tbw_chi2_simple += 1;
        if( is_tqh_quark_simple ) accuracy_tqh_chi2_simple += 1;
        if(isMatched_tbw_simple && is_tqh_quark_simple) accuracy_chi2_simple += 1.;

        if( isMatched_tbw_modified ) accuracy_tbw_chi2_modified += 1;
        if( is_tqh_quark_modified ) accuracy_tqh_chi2_modified += 1;
        if(isMatched_tbw_modified && is_tqh_quark_modified) accuracy_chi2_modified += 1.;
        //if(M1_simple!=-999 && isMatched_tbw_simple && is_tqh_quark_simple) accuracy_chi2_simple += 1.;
        //if(M1_modified!=-999 && isMatched_tbw_modified && is_tqh_quark_modified) accuracy_chi2_modified += 1.;

        if( isMatched )  accuracy_tbw_yfyj += 1;
        if(is_tqh_quark) accuracy_tqh_yfyj += 1;
        if(isMatched && is_tqh_quark) accuracy_yfyj += 1;

        if(isMatched_tbw_improved) accuracy_tbw_chi2_improved += 1;
        if(is_tqh_quark_improved) accuracy_tqh_chi2_improved += 1;
        if(isFinalMatched) accuracy_chi2_improved += 1;
        //}}}
        */
        /*
        // two methods for top reconstruction(ST){{{
        //### chi-2 Study{{{
        // Define variables{{{
        bool isMatched, is_tqh_quark=false;//boolean used to check matching
        bool isMatched_tbw_simple, is_tqh_quark_simple=false;//boolean used to check matching
        bool isMatched_tbw_modified, is_tqh_quark_modified=false;//boolean used to check matching
        bool isMatched_tbw_improved, is_tqh_quark_improved=false;//boolean used to check matching
        //-------------------------//
        int index_q_simple = -999;
        std::vector<int> id_jet_chi2_simple(2, -999);
        std::vector<int> index_jet_chi2_simple(2, -999);
        std::vector<TLorentzVector> jet_chi2_simple(2);
        std::vector<TLorentzVector> genParticle_jet_chi2_simple(2);
        //-------------------------//
        int index_q_modified = -999;
        std::vector<int> id_jet_chi2_modified(2, -999);
        std::vector<int> index_jet_chi2_modified(2, -999);
        std::vector<TLorentzVector> jet_chi2_modified(2);
        std::vector<TLorentzVector> genParticle_jet_chi2_modified(2);
        //-------------------------//
        double chi2_min = 99999;
        //}}}
        // simple & modified chi-2 methods{{{
        // chi-2 simple {{{
        chi2_min = 99999;
        for(int i=0; i<mytree.num_jets; ++i){
            if(i==index_bjet) continue;//bypass bjet
            for(int j=i+1; j<mytree.num_jets; ++j){
                if(j==index_bjet) continue;//bypass bjet
                TLorentzVector w_candidate = Jets[i] + Jets[j];
                double w_mass = w_candidate.M();
                TLorentzVector top_candidate = w_candidate + bjet;
                double t_mass = top_candidate.M();
                double chi2 = Chi2_calculator_simple(w_mass, t_mass);
                if(chi2 < chi2_min){
                    index_jet_chi2_simple[0] = i;
                    index_jet_chi2_simple[1] = j;
                    jet_chi2_simple[0] = Jets[i];
                    jet_chi2_simple[1] = Jets[j];
                    id_jet_chi2_simple[0] = Jets_GenFalavor[i];
                    id_jet_chi2_simple[1] = Jets_GenFalavor[j];
                    chi2_min = chi2;
                }
            }
        }//end of looping jets
        //}}}
        // chi-2 modified {{{
        chi2_min = 99999;
        for(int i=0; i<mytree.num_jets; ++i){
            if(i==index_bjet) continue;//bypass bjet
            for(int j=i+1; j<mytree.num_jets; ++j){
                if(j==index_bjet) continue;//bypass bjet
                TLorentzVector w_candidate = Jets[i] + Jets[j];
                double w_mass = w_candidate.M();
                TLorentzVector top_candidate = w_candidate + bjet;
                double t_mass = top_candidate.M();
                double chi2 = Chi2_calculator_modified(w_mass, t_mass);
                if(chi2 < chi2_min){
                    index_jet_chi2_modified[0] = i;
                    index_jet_chi2_modified[1] = j;
                    jet_chi2_modified[0] = Jets[i];
                    jet_chi2_modified[1] = Jets[j];
                    id_jet_chi2_modified[0] = Jets_GenFalavor[i];
                    id_jet_chi2_modified[1] = Jets_GenFalavor[j];
                    chi2_min = chi2;
                }
            }
        }//end of looping jets
        //}}}
        // }}}
        // GenMatching{{{
        //----- MC truth -----//
        //========================================//
        //-----  Start GenMatching for Jets  -----//
        //========================================//
        // GenMatching: find the gen particle (MC truth) for each jet (reconstructed). 
        // pdgID: (1, 2, 3, 4, 5, 6) = (d, u, s, c, b, t)
        // This is the simplest version. Identify the corresponding gen particle by selecting smallest deltaR(gen, jet).
        // One can try to print out the info of pt, eta, phi, energy, and deltaR of jet and corresponding gen particle to see if they are matched.
        std::vector<int> index_GenParticles;
        int id_genParticle_bjet=-999;
        int id_genParticle_jet_chi2_simple[2]; std::fill_n(id_genParticle_jet_chi2_simple, 2, -999);
        int id_genParticle_jet_chi2_modified[2]; std::fill_n(id_genParticle_jet_chi2_modified, 2, -999);
        double delta_R = 0;
        //printf("ientry = %d\n", ientry);
        genParticle_bjet = GetGenParticle(bjet, treeReader.GenPartInfo_size, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass, treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID, index_GenParticles, id_genParticle_bjet);

        //int momPdgID = abs(treeReader.GenPartInfo_MomPdgID->at(index_GenParticles[0])) ;
        //if( abs(id_genParticle_bjet)==5 ) printf("[debug] MomPdgID = %d\n", abs(treeReader.GenPartInfo_MomPdgID->at(index_GenParticles[0]))  );
        //if( abs(id_genParticle_bjet)==5 && abs(treeReader.GenPartInfo_MomPdgID->at(index_GenParticles[0]))==6 )
        //    printf("[debug] checked!\n");

        if( abs(id_genParticle_bjet)==5 && abs(treeReader.GenPartInfo_MomPdgID->at(index_GenParticles[0]))==6 ) count_bjet_is_bquark+=1;
        //if( abs(id_genParticle_bjet)==5 ) count_bjet_is_bquark+=1;
        //---
        genParticle_jet_chi2_simple[0] = GetGenParticle(jet_chi2_simple[0], treeReader.GenPartInfo_size, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass, treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID, index_GenParticles, id_genParticle_jet_chi2_simple[0]);
        //---
        genParticle_jet_chi2_simple[1] = GetGenParticle(jet_chi2_simple[1], treeReader.GenPartInfo_size, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass, treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID, index_GenParticles, id_genParticle_jet_chi2_simple[1]);
        //---
        genParticle_jet_chi2_modified[0] = GetGenParticle(jet_chi2_modified[0], treeReader.GenPartInfo_size, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass, treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID, index_GenParticles, id_genParticle_jet_chi2_modified[0]);
        //---
        genParticle_jet_chi2_modified[1] = GetGenParticle(jet_chi2_modified[1], treeReader.GenPartInfo_size, treeReader.GenPartInfo_Pt, treeReader.GenPartInfo_Eta, treeReader.GenPartInfo_Phi, treeReader.GenPartInfo_Mass, treeReader.GenPartInfo_Status, treeReader.GenPartInfo_PdgID, index_GenParticles, id_genParticle_jet_chi2_modified[1]);
        // }}}
        // Reconstruct Mass (M2){{{
        //----- reco -----//
        TLorentzVector w_candidate_chi2_simple = jet_chi2_simple[0] + jet_chi2_simple[1];
        TLorentzVector top_candidate_chi2_simple = w_candidate_chi2_simple + bjet;
        mytree.Mass_w_candidate_chi2_simple = w_candidate_chi2_simple.M();
        mytree.Mass_top_candidate_chi2_simple = top_candidate_chi2_simple.M();
        hist_mass_w_candidate_chi2_simple->Fill(w_candidate_chi2_simple.M());
        hist_mass_top_candidate_chi2_simple->Fill(top_candidate_chi2_simple.M());
        //--------------------
        TLorentzVector w_candidate_chi2_modified = jet_chi2_modified[0] + jet_chi2_modified[1];
        TLorentzVector top_candidate_chi2_modified = w_candidate_chi2_modified + bjet;
        mytree.Mass_w_candidate_chi2_modified = w_candidate_chi2_modified.M();
        mytree.Mass_top_candidate_chi2_modified = top_candidate_chi2_modified.M();
        hist_mass_w_candidate_chi2_modified->Fill(w_candidate_chi2_modified.M());
        hist_mass_top_candidate_chi2_modified->Fill(top_candidate_chi2_modified.M());
        //----- gen -----//
        //reject candidates according to the MC truth
        //1. reject candidate with id_genParticle = -999 (soft jet, not from hard process)
        //2. reject candidate with id_genParticle_jet == 5 (w_boson cannot decay to b jet)
        //3. reject candidate with 2 jets having same genParticle
        TLorentzVector gen_w_candidate_chi2_simple = genParticle_jet_chi2_simple[0] + genParticle_jet_chi2_simple[1];
        TLorentzVector gen_top_candidate_chi2_simple = gen_w_candidate_chi2_simple + genParticle_bjet;
        hist_mass_gen_w_candidate_chi2_simple->Fill(gen_w_candidate_chi2_simple.M());
        hist_mass_gen_top_candidate_chi2_simple->Fill(gen_top_candidate_chi2_simple.M());
        //----- check -----//
        if(gen_w_candidate_chi2_simple.M() < 20){
            check_M20_simple += 1;
            if( !(id_genParticle_jet_chi2_simple[0]==-999 || id_genParticle_jet_chi2_simple[1]==-999) ){
                check_gen_exclude_id_999 += 1;
                if( !(genParticle_jet_chi2_simple[0].Pt() == genParticle_jet_chi2_simple[1].Pt()) ){
                    check_gen_exclude_the_same += 1;
                    //kinematics_report("jet[0]", jet_chi2_simple[0], id_jet_chi2_simple[0], "gen[0]", genParticle_jet_chi2_simple[0], id_genParticle_jet_chi2_simple[0]);
                    //kinematics_report("jet[1]", jet_chi2_simple[1], id_jet_chi2_simple[1], "gen[1]", genParticle_jet_chi2_simple[1], id_genParticle_jet_chi2_simple[1]);
                    //printf("\n");
                }
            }
        }

        isMatched_tbw_simple = isMatched_with_Gen_tbw(treeReader.GenPartInfo_PdgID, treeReader.GenPartInfo_MomPdgID, index_GenParticles[0], index_GenParticles[1], index_GenParticles[2] );//bjet, jet12(simple)
        if( isMatched_tbw_simple ){
            mytree.Mass_gen_w_candidate_chi2_simple = gen_w_candidate_chi2_simple.M();
            mytree.Mass_gen_top_candidate_chi2_simple = gen_top_candidate_chi2_simple.M();
            hist_mass_conditioned_gen_w_candidate_chi2_simple->Fill(gen_w_candidate_chi2_simple.M());
            hist_mass_conditioned_gen_top_candidate_chi2_simple->Fill(gen_top_candidate_chi2_simple.M());
            //old test code{{{
            if(gen_w_candidate_chi2_simple.M() < 20){// debug purpose
                //printf("\nCheck:\n");
                //kinematics_report("bjet", bjet, id_bjet, "genP", genParticle_bjet, id_genParticle_bjet);
                //kinematics_report("jet_chi2_simple[0]", jet_chi2_simple[0], id_jet_chi2_simple[0], "MatchedGenParticle", genParticle_jet_chi2_simple[0], id_genParticle_jet_chi2_simple[0]);
                //kinematics_report("jet_chi2_simple[1]", jet_chi2_simple[1], id_jet_chi2_simple[1], "MatchedGenParticle", genParticle_jet_chi2_simple[1], id_genParticle_jet_chi2_simple[1]);
                //printf("[INFO] M_gen_jj  = %6.2f\n", gen_w_candidate_chi2_simple.M());
                //printf("[INFO] M_gen_bjj = %6.2f\n", gen_top_candidate_chi2_simple.M());
            }
            //}}}
        }
        //--------------------
        TLorentzVector gen_w_candidate_chi2_modified = genParticle_jet_chi2_modified[0] + genParticle_jet_chi2_modified[1];
        TLorentzVector gen_top_candidate_chi2_modified = gen_w_candidate_chi2_modified + genParticle_bjet;
        hist_mass_gen_w_candidate_chi2_modified->Fill(gen_w_candidate_chi2_modified.M());
        hist_mass_gen_top_candidate_chi2_modified->Fill(gen_top_candidate_chi2_modified.M());

        isMatched_tbw_modified = isMatched_with_Gen_tbw(treeReader.GenPartInfo_PdgID, treeReader.GenPartInfo_MomPdgID, index_GenParticles[0], index_GenParticles[3], index_GenParticles[4] );//bjet, jet12(modified)
        if( isMatched_tbw_modified ){
            mytree.Mass_gen_w_candidate_chi2_modified = gen_w_candidate_chi2_modified.M();
            mytree.Mass_gen_top_candidate_chi2_modified = gen_top_candidate_chi2_modified.M();
            hist_mass_conditioned_gen_w_candidate_chi2_modified->Fill(gen_w_candidate_chi2_modified.M());
            hist_mass_conditioned_gen_top_candidate_chi2_modified->Fill(gen_top_candidate_chi2_modified.M());
        }
        //--------------------
        // }}}
        //}}}
        //}}}
        // Evaluation{{{
        bool tbwIsCorrectlyMatched;
        if(tbwCanBeReconstructed){
            tbwIsCorrectlyMatched = check_tbwIsCorrectlyMatched(index_jet_chi2_simple, jetIndex_momPdgID_is_wboson, index_bjet, jetIndex_is_bquarkFromSMtop);
            if(tbwIsCorrectlyMatched) counter_tbwIsCorrectlyMatched_simple += 1;
            tbwIsCorrectlyMatched = check_tbwIsCorrectlyMatched(index_jet_chi2_modified, jetIndex_momPdgID_is_wboson, index_bjet, jetIndex_is_bquarkFromSMtop);
            if(tbwIsCorrectlyMatched) counter_tbwIsCorrectlyMatched_modified += 1;
        }
        //}}}
        */
        //}}}
        mytree.Fill();
    }// End of event loop.
    //==================================================//
    //---------------------  Report  -------------------//
    //==================================================//
    hist_entries_bjet_wp_selected -> Draw("b");
    hist_entries_bjet_wp_matched  -> Draw("b same");
    PrintCountsAndRatio("[re-visit] selection eff. loose",  hist_entries_bjet_wp_selected->GetBinContent(1), nentries);
    PrintCountsAndRatio("[re-visit] selection eff. medium", hist_entries_bjet_wp_selected->GetBinContent(2), nentries);
    PrintCountsAndRatio("[re-visit] selection eff. tight",  hist_entries_bjet_wp_selected->GetBinContent(3), nentries);
    PrintCountsAndRatio("[re-visit] matching  eff. loose",   hist_entries_bjet_wp_matched->GetBinContent(1), hist_entries_bjet_wp_selected->GetBinContent(1));
    PrintCountsAndRatio("[re-visit] matching  eff. medium",  hist_entries_bjet_wp_matched->GetBinContent(2), hist_entries_bjet_wp_selected->GetBinContent(2));
    PrintCountsAndRatio("[re-visit] matching  eff. tight",   hist_entries_bjet_wp_matched->GetBinContent(3), hist_entries_bjet_wp_selected->GetBinContent(3));

    //report{{{
    PrintCountsAndRatio("counter_index_bjet_method_2", counter_index_bjet_method_2, Nevents_pass_selection[0][0]);
    PrintCountsAndRatio("anti_counter_index_bjet_method_2", anti_counter_index_bjet_method_2, Nevents_pass_selection[0][0]);
    PrintCountsAndRatio("counter_genmatched_index_bjet_unphysical", counter_genmatched_index_bjet_unphysical, Nevents_pass_selection[0][0]);
    PrintCountsAndRatio("anti_counter_genmatched_index_bjet_unphysical", anti_counter_genmatched_index_bjet_unphysical, Nevents_pass_selection[0][0]);
    PrintCountsAndRatio("counter_not_matched_but_is_in_the_list", counter_not_matched_but_is_in_the_list, Nevents_pass_selection[0][0]);
    PrintCountsAndRatio("counter_difference", counter_difference, Nevents_pass_selection[0][0]);
    

    TString string_bool[2] = {"atLeastOne", "exactlyOne"};
    TString string_wp[3] = {"loose", "medium", "tight"};

    for(int ibool=0; ibool<2; ++ibool){
    for(int iwp=0; iwp<3; ++iwp){
        if(ibool==0 && iwp==0){
        report_rate( "Nevents_pass_selection" + string_bool[ibool] + "_" + string_wp[iwp], Nevents_pass_selection[iwp][ibool], nentries);
        report_rate( "Nevents_tbwCanBeReconstructed", counter_selectedJets_tbwCanBeReconstructed, Nevents_pass_selection[iwp][ibool]);
        report_rate( "Nevents_bjj_can_be_matched", counter_bjj_can_be_matched, Nevents_pass_selection[iwp][ibool]);

        printf("------------------------------\n");
        report_rate( "Nevents_sigIsCorrectlyReco_simple_" + string_bool[ibool] + "_" + string_wp[iwp]  , counter_sigIsCorrectlyReco_simple[iwp][ibool]      , Nevents_pass_selection[iwp][ibool]);
        report_rate( "Nevents_sigIsCorrectlyReco_modified_" + string_bool[ibool] + "_" + string_wp[iwp]  , counter_sigIsCorrectlyReco_modified[iwp][ibool]    , Nevents_pass_selection[iwp][ibool]);
        report_rate( "Nevents_sigIsCorrectlyReco_improved_" + string_bool[ibool] + "_" + string_wp[iwp]  , counter_sigIsCorrectlyReco_improved[iwp][ibool]    , Nevents_pass_selection[iwp][ibool]);

        report_rate( "Nevents_tqhIsCorrectlyMatched_simple_" + string_bool[ibool] + "_" + string_wp[iwp]  , counter_tqhIsCorrectlyMatched_simple[iwp][ibool]   , Nevents_pass_selection[iwp][ibool]);
        report_rate( "Nevents_tqhIsCorrectlyMatched_modified_" + string_bool[ibool] + "_" + string_wp[iwp]  , counter_tqhIsCorrectlyMatched_modified[iwp][ibool] , Nevents_pass_selection[iwp][ibool]);
        report_rate( "Nevents_tqhIsCorrectlyMatched_improved_" + string_bool[ibool] + "_" + string_wp[iwp]  , counter_tqhIsCorrectlyMatched_improved[iwp][ibool] , Nevents_pass_selection[iwp][ibool]);

        report_rate( "Nevents_tbwIsCorrectlyMatched_simple_" + string_bool[ibool] + "_" + string_wp[iwp]  , counter_tbwIsCorrectlyMatched_simple[iwp][ibool]   , Nevents_pass_selection[iwp][ibool]);
        report_rate( "Nevents_tbwIsCorrectlyMatched_modified_" + string_bool[ibool] + "_" + string_wp[iwp]  , counter_tbwIsCorrectlyMatched_modified[iwp][ibool] , Nevents_pass_selection[iwp][ibool]);
        report_rate( "Nevents_tbwIsCorrectlyMatched_improved_" + string_bool[ibool] + "_" + string_wp[iwp]  , counter_tbwIsCorrectlyMatched_improved[iwp][ibool] , Nevents_pass_selection[iwp][ibool]);

        report_rate( "counter_bjet_is_matched_simple_" + string_bool[ibool] + "_" + string_wp[iwp]  , counter_bjet_is_matched_simple[iwp][ibool] , Nevents_pass_selection[iwp][ibool]);
        report_rate( "counter_bjet_is_matched_modified_" + string_bool[ibool] + "_" + string_wp[iwp]  , counter_bjet_is_matched_modified[iwp][ibool] , Nevents_pass_selection[iwp][ibool]);
        report_rate( "counter_bjet_is_matched_improved_" + string_bool[ibool] + "_" + string_wp[iwp]  , counter_bjet_is_matched_improved[iwp][ibool] , Nevents_pass_selection[iwp][ibool]);
        printf("------------------------------\n\n");
        }
    }
    }
    //}}}

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    //TLegend *legend = new TLegend(0.50,0.20,0.88,0.45);
    TLegend *legend = new TLegend(0.60, 0.70, 0.85, 0.85);
    legend->SetTextFont(43);
    legend->SetLineColor(0);

    printf("Hello World\n");
    //hist_entries_bjet_wp.png{{{
    gPad->SetGrid();
    gPad->SetTicks(0,1);

    ph_setBarPlot(hist_entries_bjet_wp_selected, "f", 4, 0.1);
    ph_setBarPlot(hist_entries_bjet_wp_matched, "f", 2, 0.6);
    vector<TH1D*> vec_hists = {hist_entries_bjet_wp_selected, hist_entries_bjet_wp_matched};
    ph_setMaximumScope(vec_hists);
    ph_setMaximumScope(vec_hists, nentries);


    hist_entries_bjet_wp_selected -> Draw("b");
    hist_entries_bjet_wp_matched  -> Draw("b same");

    //ph_setLegendPosition(legend, 0.50, 0.60, 0.88, 0.88);
    legend->AddEntry(hist_entries_bjet_wp_selected, "selected", "f");
    legend->AddEntry(hist_entries_bjet_wp_matched,  "bjet is matched", "f");
    legend->Draw("same");

    c1->SaveAs(Form("%s/hist_entries_bjet_wp.png", output_histDir));
    //}}}

    /*
    //bjet rate-WP plots{{{
    // Normalization{{{
    for(int i=0; i<3; ++i){
        //loop over 3 different b-tagged wp
        //the Nevents_pass_selection depends on wp because the multiplicity selection on b-jet!
        double content;

        content = hist_rate_bquark_is_included -> GetBinContent(i+1) / (double)Nevents_pass_selection[i][0];
        hist_rate_bquark_is_included -> SetBinContent(i+1, content);

        content = hist_rate_leading_bjet_is_bquark -> GetBinContent(i+1) / (double)Nevents_pass_selection[i][0];
        hist_rate_leading_bjet_is_bquark -> SetBinContent(i+1, content);

        content = hist_rate_permuting_bjet_is_bquark_sim -> GetBinContent(i+1) / (double)Nevents_pass_selection[i][0];
        hist_rate_permuting_bjet_is_bquark_sim -> SetBinContent(i+1, content);

        content = hist_rate_permuting_bjet_is_bquark_mod -> GetBinContent(i+1) / (double)Nevents_pass_selection[i][0];
        hist_rate_permuting_bjet_is_bquark_mod -> SetBinContent(i+1, content);

        content = hist_rate_permuting_bjet_is_bquark_imp -> GetBinContent(i+1) / (double)Nevents_pass_selection[i][0];
        hist_rate_permuting_bjet_is_bquark_imp -> SetBinContent(i+1, content);
    }
    //}}}   

    gPad->SetTicks(0,1);
    hist_rate_bquark_is_included->SetStats(0);
    hist_rate_bquark_is_included->SetMaximum(1.);
    hist_rate_bquark_is_included->SetMinimum(0.);
    hist_rate_bquark_is_included->SetLineWidth(2);
    hist_rate_bquark_is_included->SetLineColor(kBlue);
    hist_rate_bquark_is_included->GetXaxis()->SetLabelSize(25);
    hist_rate_bquark_is_included->GetXaxis()->SetLabelFont(43);
    hist_rate_bquark_is_included->GetYaxis()->SetTitleSize(25);
    hist_rate_bquark_is_included->GetYaxis()->SetTitleFont(43);
    hist_rate_bquark_is_included->Draw();

    hist_rate_leading_bjet_is_bquark->SetLineWidth(2);
    hist_rate_leading_bjet_is_bquark->SetLineColor(2);
    hist_rate_leading_bjet_is_bquark->Draw("same");

    hist_rate_permuting_bjet_is_bquark_sim->SetLineWidth(2);
    hist_rate_permuting_bjet_is_bquark_sim->SetLineColor(kGreen);
    hist_rate_permuting_bjet_is_bquark_sim->Draw("same");

    hist_rate_permuting_bjet_is_bquark_mod->SetLineWidth(2);
    hist_rate_permuting_bjet_is_bquark_mod->SetLineColor(kGreen+2);
    hist_rate_permuting_bjet_is_bquark_mod->Draw("same");

    hist_rate_permuting_bjet_is_bquark_imp->SetLineWidth(2);
    hist_rate_permuting_bjet_is_bquark_imp->SetLineColor(kGreen+4);
    hist_rate_permuting_bjet_is_bquark_imp->Draw("same");

    legend->Clear();
    legend->AddEntry(hist_rate_bquark_is_included,      "the bjet is included", "l");
    legend->AddEntry(hist_rate_leading_bjet_is_bquark,  "the bjet is the leading one", "l");
    legend->AddEntry(hist_rate_permuting_bjet_is_bquark_sim,  "simple chi-2 method", "l");
    legend->AddEntry(hist_rate_permuting_bjet_is_bquark_mod,  "modified chi-2 method", "l");
    legend->AddEntry(hist_rate_permuting_bjet_is_bquark_imp,  "improved chi-2 method", "l");
    legend->Draw("same");

    c1->SaveAs(Form("%s/hist_rate_bquark_is_included_atleastOne.png", output_histDir))      ;
    //}}}
    */


    /* skip tmp{{{
    //--- report(TT) ---//
    // top correct rate{{{
    printf("--------------------------------------------------\n");
    report_rate( "Nevents_pass_selection", Nevents_pass_selection, nentries);
    printf("--------------------------------------------------\n");
    report_rate( "Nevents_tbwCanBeReconstructed", counter_selectedJets_tbwCanBeReconstructed, Nevents_pass_selection);
    report_rate( "Nevents_tbwIsCorrectlyMatched_yfyj", counter_tbwIsCorrectlyMatched_yfyj, Nevents_pass_selection);
    report_rate( "Nevents_tbwIsCorrectlyMatched_simple", counter_tbwIsCorrectlyMatched_simple, Nevents_pass_selection);
    report_rate( "Nevents_tbwIsCorrectlyMatched_modified", counter_tbwIsCorrectlyMatched_modified, Nevents_pass_selection);
    report_rate( "Nevents_tbwIsCorrectlyMatched_improved", counter_tbwIsCorrectlyMatched_improved, Nevents_pass_selection);
    printf("--------------------------------------------------\n");
    report_rate( "Nevents_tqhCanBeReconstructed", counter_selectedJets_tqhCanBeReconstructed, Nevents_pass_selection);
    report_rate( "Nevents_tqhIsCorrectlyMatched_yfyj", counter_tqhIsCorrectlyMatched_yfyj, Nevents_pass_selection);
    report_rate( "Nevents_tqhIsCorrectlyMatched_simple", counter_tqhIsCorrectlyMatched_simple, Nevents_pass_selection);
    report_rate( "Nevents_tqhIsCorrectlyMatched_modified", counter_tqhIsCorrectlyMatched_modified, Nevents_pass_selection);
    report_rate( "Nevents_tqhIsCorrectlyMatched_improved", counter_tqhIsCorrectlyMatched_improved, Nevents_pass_selection);
    printf("--------------------------------------------------\n");
    report_rate( "Nevents_sigCanBeReconstructed", counter_selectedJets_sigCanBeReconstructed, Nevents_pass_selection);
    report_rate( "Nevents_sigIsCorrectlyReco_yfyj", counter_sigIsCorrectlyReco_yfyj, Nevents_pass_selection);
    report_rate( "Nevents_sigIsCorrectlyReco_simple", counter_sigIsCorrectlyReco_simple, Nevents_pass_selection);
    report_rate( "Nevents_sigIsCorrectlyReco_modified", counter_sigIsCorrectlyReco_modified, Nevents_pass_selection);
    report_rate( "Nevents_sigIsCorrectlyReco_improved", counter_sigIsCorrectlyReco_improved, Nevents_pass_selection);
    printf("--------------------------------------------------\n");
    report_rate( "counter_mismatched", counter_mismatched, Nevents_pass_selection);
    report_rate( "counter_intersection", counter_intersection, Nevents_pass_selection);
    report_rate( "counter_M_notC", counter_M_notC, Nevents_pass_selection);
    report_rate( "counter_notM_C", counter_notM_C, Nevents_pass_selection);
    report_rate( "counter_notM_notC", counter_notM_notC, Nevents_pass_selection);
    printf("--------------------------------------------------\n");
    //}}}
    // bjet correct rate {{{
    printf("[INFO] bin1 bin2 bin3; correcat rate in bin2; accuracy with exactly one particle cut\n");
    hist_bin_fraction(hist_num_gen_bquark, "gen_bquark", hist_num_gen_bquark->GetBinContent(2));
    hist_bin_fraction(hist_num_bjets_loose, "bjets_loose", count_bjet_is_bquark_loose);
    hist_bin_fraction(hist_num_bjets_medium, "bjets_medium", count_bjet_is_bquark_medium);
    hist_bin_fraction(hist_num_bjets_tight, "bjets_tight", count_bjet_is_bquark_tight);
    double percentage_bjet_is_bquark = (double) count_bjet_is_bquark / (double) Nevents_pass_selection;
    printf("--------------------------------------------------\n");
    printf("[INFO] percentage_bjet_is_bquark = %f\n", percentage_bjet_is_bquark);
    //}}}
    // print accuracy{{{
    // calculate accuracy{{{
    accuracy_chi2_simple /= (double) Nevents_pass_selection;
    accuracy_chi2_modified /= (double) Nevents_pass_selection;
    accuracy_chi2_improved /= (double) Nevents_pass_selection;
    accuracy_yfyj /= (double) Nevents_pass_selection;
    accuracy_tbw_yfyj /= (double) Nevents_pass_selection;
    accuracy_tbw_chi2_simple /= (double) Nevents_pass_selection;
    accuracy_tbw_chi2_modified /= (double) Nevents_pass_selection;
    accuracy_tbw_chi2_improved /= (double) Nevents_pass_selection;
    accuracy_tqh_yfyj /= (double) Nevents_pass_selection;
    accuracy_tqh_chi2_simple /= (double) Nevents_pass_selection;
    accuracy_tqh_chi2_modified /= (double) Nevents_pass_selection;
    accuracy_tqh_chi2_improved /= (double) Nevents_pass_selection;
    //}}}
    printf("--------------------------------------------------\n");
    printf("[INFO] accuracy_tbw_yfyj = %6.4f\n", accuracy_tbw_yfyj);
    printf("[INFO] accuracy_tbw_chi2_simple = %6.4f\n", accuracy_tbw_chi2_simple);
    printf("[INFO] accuracy_tbw_chi2_modified = %6.4f\n", accuracy_tbw_chi2_modified);
    printf("[INFO] accuracy_tbw_chi2_improved = %6.4f\n", accuracy_tbw_chi2_improved);
    printf("--------------------------------------------------\n");
    printf("[INFO] accuracy_tqh_yfyj = %6.4f\n", accuracy_tqh_yfyj);
    printf("[INFO] accuracy_tqh_chi2_simple = %6.4f\n", accuracy_tqh_chi2_simple);
    printf("[INFO] accuracy_tqh_chi2_modified = %6.4f\n", accuracy_tqh_chi2_modified);
    printf("[INFO] accuracy_tqh_chi2_improved = %6.4f\n", accuracy_tqh_chi2_improved);
    printf("--------------------------------------------------\n");
    printf("[INFO] accuracy_yfyj = %6.4f\n", accuracy_yfyj);
    printf("[INFO] accuracy_chi2_simple = %6.4f\n", accuracy_chi2_simple);
    printf("[INFO] accuracy_chi2_modified = %6.4f\n", accuracy_chi2_modified);
    printf("[INFO] accuracy_chi2_improved = %6.4f\n", accuracy_chi2_improved);
    printf("--------------------------------------------------\n");
    //}}}
    // hadronic chi-2 Study (skipped){{{
//    //### performance: width, accuracy{{{
//    printf("[CHECK-2] Nevents_pass_selection = %d (%6.2f)\n", Nevents_pass_selection, 100 * (double)Nevents_pass_selection / (double)Nevents_pass_selection);
//    int entries_bquark = hist_num_gen_bquark->GetEntries();
//    printf("[CHECK] num_bquark = 1 (%6.3f%%)\n", 100*(double)check_num_bquak_is_one/(double)entries_bquark);
//    printf("[CHECK] num_bquark = 2 (%6.3f%%)\n", 100*(double)check_num_bquak_is_two/(double)entries_bquark);
//    printf("[CHECK] num_bquark = 3 (%6.3f%%)\n", 100*(double)check_num_bquak_is_three/(double)entries_bquark);
//    printf("[CHECK] gen_reco_W M < 20:    %6d (%5.2f %%) (%6.2f %%)\n", check_M20_simple, 100 * check_M20_simple/(double)Nevents_pass_selection, 100 * (double)check_M20_simple/(double)check_M20_simple);
//    printf("[CHECK] gen id != -999:       %6d (%5.2f %%) (%5.2f %%)\n", check_gen_exclude_id_999, 100 * check_gen_exclude_id_999/(double)Nevents_pass_selection, 100 * (double)check_gen_exclude_id_999/(double)check_M20_simple);
//    printf("[CHECK] exclude the same gen: %6d (%5.2f %%) (%5.2f %%)\n", check_gen_exclude_the_same, 100 * check_gen_exclude_the_same/(double)Nevents_pass_selection, 100 * (double)check_gen_exclude_the_same/(double)check_M20_simple);
//    printf("[CHECK] M1 < 20: %d\n", check_M1_20);
//    //}}}
//    //# mean, error of hists{{{
//    hist_report(hist_mass_w_candidate_yfyj, "w yfyj");
//    hist_report(hist_mass_w_candidate_chi2_simple, "w simple");
//    hist_report(hist_mass_w_candidate_chi2_modified, "w modified");
//    hist_report(hist_mass_t2_candidate_yfyj, "t2 yfyj");
//    hist_report(hist_mass_top_candidate_chi2_simple, "top simple");
//    hist_report(hist_mass_top_candidate_chi2_modified, "top modified");
//    printf("//--------------------//\n");
//    hist_report(hist_mass_gen_w_candidate_yfyj, "gen w yfyj");
//    hist_report(hist_mass_gen_w_candidate_chi2_simple, "gen w simple");
//    hist_report(hist_mass_gen_w_candidate_chi2_modified, "gen w modified");
//    hist_report(hist_mass_gen_t2_candidate_yfyj, "gen t2 yfyj");
//    hist_report(hist_mass_gen_top_candidate_chi2_simple, "gen top simple");
//    hist_report(hist_mass_gen_top_candidate_chi2_modified, "gen top modified");
//    printf("//--------------------//\n");
//    hist_report(hist_mass_conditioned_gen_w_candidate_yfyj, "conditioned_gen w yfyj");
//    hist_report(hist_mass_conditioned_gen_w_candidate_chi2_simple, "conditioned_gen w simple");
//    hist_report(hist_mass_conditioned_gen_w_candidate_chi2_modified, "conditioned_gen w modified");
//    hist_report(hist_mass_conditioned_gen_t2_candidate_yfyj, "conditioned_gen t2 yfyj");
//    hist_report(hist_mass_conditioned_gen_top_candidate_chi2_simple, "conditioned_gen top simple");
//    hist_report(hist_mass_conditioned_gen_top_candidate_chi2_modified, "conditioned_gen top modified");
//    printf("//--------------------//\n");
//    //}}}
    //}}}
    //### Make plots{{{
    //TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    //--------------------------------------------------
    hf_wboson.Draw_all_hist(c1);
    hf_tbw   .Draw_all_hist(c1);
    hf_tqh   .Draw_all_hist(c1);
    hf_tqh_matched   .Draw_all_hist(c1);
    //--------------------------------------------------
    c1->SetLogy(1);
    hist_num_gen_quarks -> Draw()      ; c1->SaveAs(Form("%s/hist_num_gen_quarks.png", output_histDir))      ;
    hist_num_gen_light_quark -> Draw() ; c1->SaveAs(Form("%s/hist_num_gen_light_quark.png", output_histDir)) ;
    hist_num_gen_bquark -> Draw()      ; c1->SaveAs(Form("%s/hist_num_gen_bquark.png", output_histDir))      ;
    hist_num_selected_jets -> Draw()   ; c1->SaveAs(Form("%s/hist_num_selected_jets.png", output_histDir))   ;
    hist_mass_diphoton -> Draw()       ; c1->SaveAs(Form("%s/hist_mass_diphoton.png", output_histDir))       ;
    // hist_num_bjets_bquark.png{{{
    c1->SetLogy(0);
    TLegend *legend_bjet = new TLegend(0.50,0.55,0.85,0.85);
    hist_num_gen_bquark           -> SetStats(0);
    hist_num_gen_bquark           -> GetXaxis() -> SetTitle("Number of particles");
    hist_num_gen_bquark           -> SetLineWidth(2);
    hist_num_bjets_loose          -> SetLineWidth(2);
    hist_num_bjets_medium         -> SetLineWidth(2);
    hist_num_bjets_tight          -> SetLineWidth(2);
    hist_num_gen_bquark           -> SetLineColor(kRed);
    hist_num_bjets_tight          -> SetLineColor(kBlue);
    hist_num_bjets_medium         -> SetLineColor(kGreen+2);
    hist_num_bjets_loose          -> SetLineColor(kGray+2);
    hist_num_bjets_tight_matched  -> SetLineWidth(0);
    hist_num_bjets_medium_matched -> SetLineWidth(0);
    hist_num_bjets_loose_matched  -> SetLineWidth(0);
    hist_num_bjets_tight_matched  -> SetFillColor(kBlue);
    hist_num_bjets_medium_matched -> SetFillColor(kGreen+2);
    hist_num_bjets_loose_matched  -> SetFillColor(kGray+2);
    hist_num_bjets_tight_matched  -> SetFillStyle(3001);
    hist_num_bjets_medium_matched -> SetFillStyle(3001);
    hist_num_bjets_loose_matched  -> SetFillStyle(3001);
    hist_num_gen_bquark           -> Draw();
    hist_num_bjets_medium_matched -> Draw("same");
    hist_num_bjets_tight_matched  -> Draw("same");
    hist_num_bjets_loose_matched  -> Draw("same");
    hist_num_bjets_loose          -> Draw("same");
    hist_num_bjets_medium         -> Draw("same");
    hist_num_bjets_tight          -> Draw("same");
    hist_num_gen_bquark           -> Draw("same");
    legend_bjet->Clear();
    legend_bjet->SetTextSize(0.03);
    legend_bjet->AddEntry(hist_num_gen_bquark,   "b quark (gen-level)", "l");
    legend_bjet->AddEntry(hist_num_bjets_tight,  "tight b-tagged jets", "l");
    legend_bjet->AddEntry(hist_num_bjets_medium, "medium b-tagged jets", "l");
    legend_bjet->AddEntry(hist_num_bjets_loose,  "loose b-tagged jets", "l");
    legend_bjet->AddEntry(hist_num_bjets_tight_matched,  "tight b-tagged jets (matched)", "f");
    legend_bjet->AddEntry(hist_num_bjets_medium_matched, "medium b-tagged jets (matched)", "f");
    legend_bjet->AddEntry(hist_num_bjets_loose_matched,  "loose b-tagged jets (matched)", "f");
    legend_bjet->SetLineColor(0);
    legend_bjet->Draw("same");
    c1->SaveAs(Form("%s/hist_num_bjets_bquark.png", output_dir));
    //}}}
    //MakeFinalPlots(skipped){{{
    //TLegend *legend = new TLegend(0.60,0.55,0.85,0.85);
    //MakeFinalPlots(c1, hist_mass_w_candidate_chi2_simple, hist_mass_w_candidate_chi2_modified, hist_mass_w_candidate_yfyj, legend, "chi2_study_w_candidate.png");
    //MakeFinalPlots(c1, hist_mass_gen_w_candidate_chi2_simple, hist_mass_gen_w_candidate_chi2_modified, hist_mass_gen_w_candidate_yfyj, legend, "chi2_study_gen_w_candidate.png");
    //MakeFinalPlots(c1, hist_mass_conditioned_gen_w_candidate_chi2_simple, hist_mass_conditioned_gen_w_candidate_chi2_modified, hist_mass_conditioned_gen_w_candidate_yfyj, legend, "chi2_study_conditioned_gen_w_candidate.png");

    //MakeFinalPlots(c1, hist_mass_top_candidate_chi2_simple, hist_mass_top_candidate_chi2_modified, hist_mass_t2_candidate_yfyj, legend, "chi2_study_top_candidate.png");
    //MakeFinalPlots(c1, hist_mass_gen_top_candidate_chi2_simple, hist_mass_gen_top_candidate_chi2_modified, hist_mass_gen_t2_candidate_yfyj, legend, "chi2_study_gen_top_candidate.png");
    //MakeFinalPlots(c1, hist_mass_conditioned_gen_top_candidate_chi2_simple, hist_mass_conditioned_gen_top_candidate_chi2_modified, hist_mass_conditioned_gen_t2_candidate_yfyj, legend, "chi2_study_conditioned_gen_top_candidate.png");

    //MakeFinalPlots(c1, hist_mass_top_fcnh_simple, hist_mass_top_fcnh_modified, hist_mass_t1_candidate_yfyj, legend, "chi2_study_top_fcnh.png");
    //MakeFinalPlots(c1, hist_deltaR_W_genW_simple, hist_deltaR_W_genW_modified, hist_deltaR_W_genW_yfyj, legend, "chi2_study_deltaR_w_genw.png");
    //}}}
    //###}}}

    //--- report(ST) ---//
    // top correct rate{{{
    printf("--------------------------------------------------\n");
    report_rate( "Nevents_pass_selection", Nevents_pass_selection, nentries);
    printf("--------------------------------------------------\n");
    report_rate( "Nevents_tbwCanBeReconstructed", counter_selectedJets_tbwCanBeReconstructed, Nevents_pass_selection);
    report_rate( "Nevents_tbwIsCorrectlyMatched_yfyj", counter_tbwIsCorrectlyMatched_yfyj, Nevents_pass_selection);
    report_rate( "Nevents_tbwIsCorrectlyMatched_simple", counter_tbwIsCorrectlyMatched_simple, Nevents_pass_selection);
    report_rate( "Nevents_tbwIsCorrectlyMatched_modified", counter_tbwIsCorrectlyMatched_modified, Nevents_pass_selection);
    report_rate( "Nevents_tbwIsCorrectlyMatched_improved", counter_tbwIsCorrectlyMatched_improved, Nevents_pass_selection);
    printf("--------------------------------------------------\n");
    //}}}
    // bjet correct rate {{{
    printf("[INFO] bin1 bin2 bin3; correcat rate in bin2; accuracy with exactly one particle cut\n");
    hist_bin_fraction(hist_num_gen_bquark, "gen_bquark", hist_num_gen_bquark->GetBinContent(2));
    hist_bin_fraction(hist_num_bjets_loose, "bjets_loose", count_bjet_is_bquark_loose);
    hist_bin_fraction(hist_num_bjets_medium, "bjets_medium", count_bjet_is_bquark_medium);
    hist_bin_fraction(hist_num_bjets_tight, "bjets_tight", count_bjet_is_bquark_tight);
    double percentage_bjet_is_bquark = (double) count_bjet_is_bquark / (double) Nevents_pass_selection;
    printf("--------------------------------------------------\n");
    printf("[INFO] percentage_bjet_is_bquark = %f\n", percentage_bjet_is_bquark);
    //}}}
    */
        //}}}
    // the end{{{
    fout->Write();
    fout->Close();
    return 0;
    //}}}
}
