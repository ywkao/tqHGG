#include "generalChiSquareStudy_leptonic.h"
#include "generalChiSquareStudy.h"
#include "hist_factory.h"
// histograms{{{
//==================================================//
//--------   Histograms for bjet rate-WP  ----------//
//==================================================//
const char* wp[3] = {"Loose", "Medium", "Tight"};
TH1D *hist_rate_bquark_is_included = new TH1D("hist_rate_bquark_is_included", ";;rate", 3, 0, 3);
TH1D *hist_rate_leading_bjet_is_bquark = new TH1D("hist_rate_leading_bjet_is_bquark", ";;rate", 3, 0, 3);
void initHists()
{
    for(int i=0; i<3; ++i) hist_rate_bquark_is_included -> GetXaxis() -> SetBinLabel(i+1, wp[i]);
    for(int i=0; i<3; ++i) hist_rate_leading_bjet_is_bquark -> GetXaxis() -> SetBinLabel(i+1, wp[i]);
}
//==================================================//
//--------   Histograms for mass spectra  ----------//
//==================================================//
TH1D *hist_num_leptons = new TH1D("hist_num_leptons", ";Number of leptons;Entries", 10, 0, 10);
TH1D *hist_num_genLeptons = new TH1D("hist_num_genLeptons", ";Number of genLeptons;Entries", 10, 0, 10);
TH1D *hist_num_bjets_tight = new TH1D("hist_num_bjets_tight", ";Number of tight b-tagged jets;Entries", 10, 0, 10);
TH1D *hist_num_bjets_medium = new TH1D("hist_num_bjets_medium", ";Number of medium b-tagged jets;Entries", 10, 0, 10);
TH1D *hist_num_bjets_loose = new TH1D("hist_num_bjets_loose", ";Number of loose b-tagged jets;Entries", 10, 0, 10);
TH1D *hist_num_bjets_tight_matched = new TH1D("hist_num_bjets_tight_matched", ";Number of tight b-tagged jets;Entries", 10, 0, 10);
TH1D *hist_num_bjets_medium_matched = new TH1D("hist_num_bjets_medium_matched", ";Number of medium b-tagged jets;Entries", 10, 0, 10);
TH1D *hist_num_bjets_loose_matched = new TH1D("hist_num_bjets_loose_matched", ";Number of loose b-tagged jets;Entries", 10, 0, 10);
TH1D *hist_num_gen_bquark = new TH1D("hist_num_gen_bquark", ";Number of b quark (gen-level);Entries", 10, 0, 10);
TH1D *hist_num_gen_light_quark = new TH1D("hist_num_gen_light_quark", ";Number of light quarks (gen-level);Entries", 10, 0, 10);
TH1D *hist_num_gen_quarks = new TH1D("hist_num_gen_quarks", ";Number of quarks (gen-level);Entries", 10, 0, 10);
TH1D *hist_num_selected_jets = new TH1D("hist_num_selected_jets", ";Number of selected jets;Entries", 10, 0, 10);
TH1D *hist_mass_diphoton = new TH1D("hist_mass_diphoton", ";Mass [GeV/c^2];Entries", 50, 0, 250);
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
//# other hist{{{
TH1D *hist_disc_topKinFit = new TH1D("hist_disc_topKinFit", "; discriminant; Entries", 40, 0, 1000);
TH1D *hist_mass_gen_wboson_leptonic = new TH1D("hist_mass_gen_wboson_leptonic", ";Mass [GeV/c^{2}];Entries", 32, 0, 160);
TH1D *hist_mass_gen_topquark_leptonic = new TH1D("hist_mass_gen_topquark_leptonic", ";Mass [GeV/c^{2}];Entries", 34, 0, 340);
TH1D *hist_deltaR_reco_top_higgs_leptonic = new TH1D("hist_deltaR_reco_top_higgs_leptonic", "", 40, 0, 6);
//}}}
//}}}
// hist leptonic{{{
char output_histDir[128] = "result_top_reco_study/hist_factory_leptonic";
hist_factory hf_chargedLepton (output_histDir, "chargedLepton", "", 100);
hist_factory hf_neutrino_sol1 (output_histDir, "neutrino_sol1", "quadratic", 100);
hist_factory hf_neutrino_sol2 (output_histDir, "neutrino_sol2", "quadratic", 100);
hist_factory hf_neutrino_sol2_positive (output_histDir, "neutrino_sol2_positive", "quadratic", 100);
hist_factory hf_neutrino_sol2_negative (output_histDir, "neutrino_sol2_negative", "quadratic", 100);
hist_factory hf_neutrino_topKinFit (output_histDir, "neutrino", "topKinFit", 100);
//---
hist_factory hf_wboson_quadratic (output_histDir, "wboson", "quadratic", 200);
hist_factory hf_wboson_topKinFit (output_histDir, "wboson", "topKinFit", 200);
hist_factory hf_top_quadratic (output_histDir, "top", "quadratic", 350);
hist_factory hf_tqh_quadratic (output_histDir, "tqh", "quadratic", 350);
hist_factory hf_tqh_matched (output_histDir, "tqh", "matched", 350);
hist_factory hf_top_topKinFit (output_histDir, "top", "topKinFit", 350);
//}}}
// Counters{{{
int counter_irregular_disc = 0;
int counter_coeff_D = 0;
int counter_coeff_D_gen = 0;
int counter_coeff_D_isNegative = 0;
int counter_coeff_D_isNegative_gen = 0;
int counter_coeff_D_isNegative_gen_largerThan20 = 0;
int counter_coeff_D_isNegative_gen_smallerThan20 = 0;
int counter_coeff_D_isNegative_gen_between20 = 0;
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

//int counter_tqhIsCorrectlyMatched = 0;
int counter_selectedJets_tbwCanBeReconstructed = 0;
int counter_selectedJets_tqhCanBeReconstructed = 0;
int counter_selectedJets_sigCanBeReconstructed = 0;
//int counter_lepton_is_correctly_chosen = 0;

int counter_yy = 0;
int counter_yn = 0;
int counter_ny = 0;
int counter_nn = 0;
int counter_nperm = 0;

int counter_lepton_is_correctly_chosen = 0;
int counter_bjet_is_correct = 0;
int counter_tbw_is_correct = 0;
int counter_tqhIsCorrectlyMatched = 0;
int counter_signal_is_correct = 0;
//}}}
int counter_index_bjet_method_1 = 0;
int counter_index_bjet_method_2 = 0;
int counter_index_bjet_method_3 = 0;
int counter_no_gen_matched_bjet = 0;
int counter_index_qjet_matched  = 0;
int counter_tbw_matched         = 0;
int counter_overall_matched     = 0;
int counter_the_leading_lepton_isGood = 0;


int counter_events_with_two_wbosons = 0;
int counter_events_with_morethan_genleptons = 0;

//double evaluate_neutrino_pz(TLorentzVector lepton, vector<double> met_info){{{
double evaluate_neutrino_pz(TLorentzVector lepton, vector<double> met_info)
{
    float met_pt = met_info[0];
    float met_px = met_info[1];
    float met_py = met_info[2];

    float lepton_px = lepton.Px();
    float lepton_py = lepton.Py();
    float lepton_pz = lepton.Pz();
    float lepton_energy = lepton.E();
    float coefficient_factor = ( w_boson_mass*w_boson_mass + 2*lepton_px*met_px + 2*lepton_py*met_py ) / (2.*lepton_energy);
    float coefficient_A = 1. - (lepton_pz*lepton_pz)/(lepton_energy*lepton_energy);
    float coefficient_B = 2.*coefficient_factor*lepton_pz/lepton_energy;
    float coefficient_C = met_pt*met_pt - coefficient_factor*coefficient_factor;
    float coefficient_D = coefficient_B*coefficient_B - 4.*coefficient_A*coefficient_C;
    
    float met_pz_solution_1 = 0.0;
    float met_pz_solution_2 = 0.0;

    if(coefficient_D < 0){
        counter_coeff_D_isNegative += 1;
        met_pz_solution_1 = coefficient_B / (2.*coefficient_A);
        met_pz_solution_2 = coefficient_B / (2.*coefficient_A);
    } else{
        met_pz_solution_1 = (coefficient_B + TMath::Sqrt(coefficient_D))/(2.*coefficient_A);
        met_pz_solution_2 = (coefficient_B - TMath::Sqrt(coefficient_D))/(2.*coefficient_A);
    }
    //ordering
    float larger_pz  = (abs(met_pz_solution_1) > abs(met_pz_solution_2) ) ? met_pz_solution_1 : met_pz_solution_2;
    float smaller_pz = (abs(met_pz_solution_1) < abs(met_pz_solution_2) ) ? met_pz_solution_1 : met_pz_solution_2;
    met_pz_solution_1 = larger_pz;
    met_pz_solution_2 = smaller_pz;

    //hists{{{
    hist_MetInfo_coeff_A -> Fill(coefficient_A);
    hist_MetInfo_coeff_B -> Fill(coefficient_B);
    hist_MetInfo_coeff_C -> Fill(coefficient_C);
    hist_MetInfo_coeff_B2A -> Fill(-coefficient_B / (2*coefficient_A));
    
    if(coefficient_D < 0){
        hist_MetInfo_coeff_D -> Fill(-sqrt(-coefficient_D));//keep tracking negative value
        hist_MetInfo_coeff_D2A -> Fill(-sqrt(-coefficient_D) / (2*coefficient_A));
    } else{
        hist_MetInfo_coeff_D -> Fill(sqrt(coefficient_D));
        hist_MetInfo_coeff_D2A -> Fill(sqrt(coefficient_D) / (2*coefficient_A));
    }

    hist_MetInfo_Pz_solution_1 -> Fill(met_pz_solution_1);
    hist_MetInfo_Pz_solution_2 -> Fill(met_pz_solution_2);
    if(coefficient_D>=0) hist_MetInfo_Pz_solution_1_positiveD -> Fill(met_pz_solution_1);
    if(coefficient_D>=0) hist_MetInfo_Pz_solution_2_positiveD -> Fill(met_pz_solution_2);
    if(coefficient_D<0) hist_MetInfo_Pz_solution_1_negativeD -> Fill(met_pz_solution_1);
    if(coefficient_D<0) hist_MetInfo_Pz_solution_2_negativeD -> Fill(met_pz_solution_2);
    //}}}

    return met_pz_solution_2;
}
//}}}
//TLorentzVector derive_reco_neutrino(TLorentzVector lepton, vector<double> met_info){{{
TLorentzVector derive_reco_neutrino(TLorentzVector lepton, vector<double> met_info)
{
    double neutrino_pz = evaluate_neutrino_pz(lepton, met_info);
    double met_pt = met_info[0];
    double met_px = met_info[1];
    double met_py = met_info[2];
    double neutrino_energy = TMath::Sqrt(met_pt*met_pt + neutrino_pz*neutrino_pz);

    TLorentzVector reco_neutrino;
    reco_neutrino.SetPxPyPzE( met_px, met_py, neutrino_pz, neutrino_energy );
    return reco_neutrino;
}
//}}}
//TLorentzVector derive_reco_wboson(TLorentzVector lepton, TLorentzVector reco_neutrino){{{
TLorentzVector derive_reco_wboson(TLorentzVector lepton, TLorentzVector reco_neutrino)
{
    TLorentzVector reco_wboson = reco_neutrino + lepton;
    return reco_wboson;
}
//}}}
//TLorentzVector derive_reco_tbw(TLorentzVector reco_wboson, TLorentzVector bjet){{{
TLorentzVector derive_reco_tbw(TLorentzVector reco_wboson, TLorentzVector bjet)
{
    TLorentzVector reco_tbw = reco_wboson + bjet;
    return reco_tbw;
}
//}}}
//int get_q_index_min_chi2(std::vector<TLorentzVector> Jets, int index_bjet, TLorentzVector diphoton){{{
int get_q_index_min_chi2(std::vector<TLorentzVector> Jets, int index_bjet, TLorentzVector diphoton)
{
    std::vector<int> indices;
    std::vector<double> top_fcnh_chi2;
    for(std::size_t i=0; i!=Jets.size(); ++i){
        if(i==index_bjet) continue; //skip the selected jets for bjet
        TLorentzVector top_fcnh_tmp = diphoton + Jets[i];
        double chi2 = (top_fcnh_tmp.M() - top_quark_mass) * (top_fcnh_tmp.M() - top_quark_mass);
        indices.push_back(i);
        top_fcnh_chi2.push_back(chi2);
        //printf("[check-ywk] q = ");
        //printf("%d, " , i);
        //printf("chi2 = %7.3f\n", chi2);
    }

    int min_index =  std::min_element(top_fcnh_chi2.begin(), top_fcnh_chi2.end()) - top_fcnh_chi2.begin();
    double min    = *std::min_element(top_fcnh_chi2.begin(), top_fcnh_chi2.end());
    //printf("[check-ywk] Jets.size() = %d, min: q = ", (int)Jets.size());
    //printf("%d, " , indices[min_index]);
    //printf("chi2 = %7.3f\n", top_fcnh_chi2[min_index]);

    int result = Jets.size() > 1 ? indices[min_index] : -1;
    return result;
}
//}}}
bool isMatched(int a, int b)
{
    if( b!= -999 ) return a==b;
    else           return false;
}

