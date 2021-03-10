#ifndef __GENERALCHISQUARESTUDY_H__
#define __GENERALCHISQUARESTUDY_H__
//includes{{{
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
#include <TString.h>
#include <TMatrixD.h>
#include <TTree.h>
#include <vector>
#include <string>
//}}}
using namespace std;

//==========================//
//---  Useful Constants  ---//
//==========================//
double pfDeepCSVJetTags_tight  = 0.8001;
double pfDeepCSVJetTags_medium = 0.4941;
double pfDeepCSVJetTags_loose  = 0.1522;
double w_boson_mass = 80.379;//GeV
double w_boson_width = 9.5;//GeV
double top_quark_mass = 173.0;//GeV
double top_quark_width = 16.3;//GeV

char output_dir[512] = "result_top_reco_study";
//==========================//
//---  Class & Function  ---//
//==========================//
void chi2_study(char *input_file, char *output_file, char* dataset);

vector<int> get_bjjq_indices_min_chi2(std::vector<TLorentzVector> Jets, std::vector<int> indices_bjet, TLorentzVector diphoton, bool is_chi2_modified);
vector<int> get_bjj_indices_min_chi2(std::vector<TLorentzVector> Jets, std::vector<int> indices_bjet, bool is_chi2_modified);
vector<int> get_indices_chi2(std::vector<TLorentzVector> Jets, int index_bjet, double &chi2_min, bool is_chi2_modified);
int get_q_index_min_chi2(std::vector<TLorentzVector> Jets, std::vector<int> indices_bjj, TLorentzVector diphoton);

std::vector<int> get_bjjq_indices_min_chi2_3x3(std::vector<TLorentzVector> Jets, std::vector<int> indices_bjet, TLorentzVector diphoton);
std::vector<int> get_bjj_indices_min_chi2_3x3(std::vector<TLorentzVector> Jets, int index_bjet, TLorentzVector diphoton, double &chi2_min);

std::vector<bool> get_matching_result(std::vector<int> reco, std::vector<int> gen);
bool check_tbwIsCorrectlyMatched(std::vector<int> index_jet_chi2, std::vector<int> index_gen_matched);

bool check_bjet_matching(int index_qjet, int index_gen_matched);
bool check_tqhIsCorrectlyMatched(int index_qjet, int jetIndex_is_quarkFromFCNtop);
bool check_tbwIsCorrectlyMatched(std::vector<int> index_jet_chi2, std::vector<int> jetIndex_momPdgID_is_wboson, int index_bjet, int jetIndex_is_bquarkFromSMtop);

// others{{{
void report_rate(const char* name, int counter, int N);
double GetM1M2_ratio(double M1, double M2);
//double GetBestM1(int num_jets, int index_bjet, std::vector<int> index_jet, TLorentzVector diphoton, std::vector<TLorentzVector> Jets);
TLorentzVector GetBestM1(double &M1, int num_jets, int index_bjet, std::vector<int> index_jet, TLorentzVector diphoton, std::vector<TLorentzVector> Jets, int &index_q, TLorentzVector &jet_q);
bool is_this_tqh_quark(TLorentzVector jet, Int_t GenPartInfo_size, std::vector<int> *GenPartInfo_MomPdgID, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_Status, std::vector<int> *GenPartInfo_PdgID);
bool isMatched_with_Gen_tbw(std::vector<int> *GenPartInfo_PdgID, std::vector<int> *GenPartInfo_MomPdgID, int index_bjet, int index_jet1, int index_jet2);
bool isMatched_with_Gen_W_Boson(TLorentzVector gen_w_sel, TH1D *&hist, Int_t GenPartInfo_size, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_PdgID);

bool CheckBJetID(TLorentzVector jet, Int_t GenPartInfo_size, std::vector<int> *GenPartInfo_MomPdgID, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_Status, std::vector<int> *GenPartInfo_PdgID);
//bool CheckBJetID(TLorentzVector jet, Int_t GenPartInfo_size, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_Status, std::vector<int> *GenPartInfo_PdgID);
TLorentzVector GetGenParticle(TLorentzVector jet, Int_t GenPartInfo_size, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_Status, std::vector<int> *GenPartInfo_PdgID, std::vector<int> &index_GenParticles, int &genParticle_PdgID);

bool checkAvailability(int index, std::vector<int> ID_IsChosen);
//void kinematics_info(const char* Title, TLorentzVector Particle);
void kinematics_info(int index, TLorentzVector Particle);
void kinematics_info(const char* Title, TLorentzVector Particle, int index);
void kinematics_report(const char* recoTitle, TLorentzVector recoParticle, int id_recoParticle, const char* genTitle, TLorentzVector genParticle, int genParticle_PdgID);
void hist_bin_fraction(TH1D *hist, const char* title, int entries_in_bin2);
void hist_report(TH1D *hist, const char* chi2_type);
double Chi2_calculator_simple(double w_mass, double t_mass);
double Chi2_calculator_modified(double w_mass, double t_mass);
double Chi2_calculator_improved(double w_mass, double t_mass, double fcnc_top_mass);
void Set2DPlot(TH2D *&h);
void MakePlots(TCanvas *c1, TH1D* hist, const char* title, const char* outputFile);
void MakeTwoPlots(TCanvas *c1, TH1D* hist_gen, TH1D* hist_reco, TLegend *legend, const char* name);
void MakeFinalPlots(TCanvas *c1, TH1D* hist_simple, TH1D* hist_modified, TH1D*hist_yfyj, TLegend *legend, const char* name);
bool isThisDataOrNot(char* dataset);
bool isThisMultiFile(char* dataset);
bool isThisMCsignal(char* dataset);

double obtain_deltaR(bool bool_print, TLorentzVector jet, int &Matched_PdgID, int &Matched_MomPdgID, Int_t GenPartInfo_size, std::vector<int> *GenPartInfo_MomPdgID, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_Status, std::vector<int> *GenPartInfo_PdgID);
void obtain_gen_matched_ID(bool bool_print, TLorentzVector jet, int &index_gen, int &Matched_PdgID, int &Matched_MomPdgID, Int_t GenPartInfo_size, std::vector<int> *GenPartInfo_MomPdgID, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_Status, std::vector<int> *GenPartInfo_PdgID);

//--- leptonic ---//
void MakePlots_coeffD(TCanvas* c1, TLegend* legend_ratio, const char* histName, const char* label1, const char* label2, TH1D* hist, TH1D* hist_positiveD, TH1D* hist_negativeD);
void PrintCountsAndRatio(const char* title, int a, int b);
void plots_comparison_two_methods(const char* dirName, const char* plotName, TCanvas *c1, TLegend *legend, TH1D* hist_quadratic, TH1D* hist_topKinFit);
//--- end of leptonic ---//

class flashggStdTreeParameters{
public:
    flashggStdTreeParameters();
    ~flashggStdTreeParameters();
    //------------------------
    Int_t           GenPartInfo_size;
    std::vector<float>   *GenPartInfo_Pt;
    std::vector<float>   *GenPartInfo_Eta;
    std::vector<float>   *GenPartInfo_Phi;
    std::vector<float>   *GenPartInfo_Mass;
    std::vector<int>     *GenPartInfo_PdgID;
    std::vector<int>     *GenPartInfo_Status;
    std::vector<int>     *GenPartInfo_nMo;
    std::vector<int>     *GenPartInfo_nDa;
    std::vector<bool>    *GenPartInfo_isHardProcess;
    std::vector<bool>    *GenPartInfo_fromHardProcessFinalState;
    std::vector<bool>    *GenPartInfo_isPromptFinalState;
    std::vector<bool>    *GenPartInfo_isDirectPromptTauDecayProductFinalState;
    std::vector<int>     *GenPartInfo_MomPdgID;
    std::vector<int>     *GenPartInfo_MomStatus;
    std::vector<float>   *GenPartInfo_MomPt;
    std::vector<float>   *GenPartInfo_MomEta;
    std::vector<float>   *GenPartInfo_MomPhi;
    std::vector<float>   *GenPartInfo_MomMass;
    //------------------------
    bool  EvtInfo_passTrigger;
    Int_t EvtInfo_NPu;
    Int_t EvtInfo_NVtx;
    float EvtInfo_Rho;
    float EvtInfo_genweight;
    //------------------------
    float DiPhoInfo_mass;
    float DiPhoInfo_pt;
    float DiPhoInfo_leadPt;
    float DiPhoInfo_leadEta;
    float DiPhoInfo_leadPhi;
    float DiPhoInfo_leadE;
    float DiPhoInfo_leadhoe;
    float DiPhoInfo_leadIDMVA;
    //------------------------
    float DiPhoInfo_subleadPt;
    float DiPhoInfo_subleadEta;
    float DiPhoInfo_subleadPhi;
    float DiPhoInfo_subleadE;
    float DiPhoInfo_subleadhoe;
    float DiPhoInfo_subleadIDMVA;
    //------------------------
    Int_t jets_size;
    std::vector<float> *JetInfo_Pt;
    std::vector<float> *JetInfo_Eta;
    std::vector<float> *JetInfo_Phi;
    std::vector<float> *JetInfo_Mass;
    std::vector<float> *JetInfo_Energy;
    std::vector<int>   *JetInfo_GenFlavor;
    std::vector<float> *JetInfo_pfDeepCSVJetTags_probb;
    std::vector<float> *JetInfo_pfDeepCSVJetTags_probbb;
    //------------------------
    Int_t           ElecInfo_Size;
    std::vector<int>     *ElecInfo_Charge;
    std::vector<float>   *ElecInfo_Pt;
    std::vector<float>   *ElecInfo_Eta;
    std::vector<float>   *ElecInfo_Phi;
    std::vector<float>   *ElecInfo_Energy;
    std::vector<float>   *ElecInfo_EtaSC;
    std::vector<float>   *ElecInfo_PhiSC;
    std::vector<float>   *ElecInfo_GsfTrackDz;
    std::vector<float>   *ElecInfo_GsfTrackDxy;
    std::vector<bool>    *ElecInfo_EGMCutBasedIDVeto;
    std::vector<bool>    *ElecInfo_EGMCutBasedIDLoose;
    std::vector<bool>    *ElecInfo_EGMCutBasedIDMedium;
    std::vector<bool>    *ElecInfo_EGMCutBasedIDTight;
    std::vector<bool>    *ElecInfo_fggPhoVeto;
    //std::vector<bool>    *ElecInfo_tmpPhoVeto;
    std::vector<float>   *ElecInfo_EnergyCorrFactor;
    std::vector<float>   *ElecInfo_EnergyPostCorrErr;
    std::vector<float>   *ElecInfo_EnergyPostCorrScaleUp;
    std::vector<float>   *ElecInfo_EnergyPostCorrScaleDown;
    std::vector<float>   *ElecInfo_EnergyPostCorrSmearUp;
    std::vector<float>   *ElecInfo_EnergyPostCorrSmearDown;
    std::vector<bool>    *ElecInfo_GenMatch;
    std::vector<int>     *ElecInfo_GenPdgID;
    std::vector<float>   *ElecInfo_GenPt;
    std::vector<float>   *ElecInfo_GenEta;
    std::vector<float>   *ElecInfo_GenPhi;
    //------------------------
    Int_t           MuonInfo_Size;
    std::vector<int>     *MuonInfo_Charge;
    std::vector<float>   *MuonInfo_MuonType;
    std::vector<float>   *MuonInfo_Pt;
    std::vector<float>   *MuonInfo_Eta;
    std::vector<float>   *MuonInfo_Phi;
    std::vector<float>   *MuonInfo_Energy;
    std::vector<float>   *MuonInfo_BestTrackDz;
    std::vector<float>   *MuonInfo_BestTrackDxy;
    std::vector<float>   *MuonInfo_PFIsoDeltaBetaCorrR04;
    std::vector<float>   *MuonInfo_TrackerBasedIsoR03;
    std::vector<bool>    *MuonInfo_CutBasedIdMedium;
    std::vector<bool>    *MuonInfo_CutBasedIdTight;
    std::vector<bool>    *MuonInfo_GenMatch;
    std::vector<int>     *MuonInfo_GenPdgID;
    std::vector<float>   *MuonInfo_GenPt;
    std::vector<float>   *MuonInfo_GenEta;
    std::vector<float>   *MuonInfo_GenPhi;
    //------------------------
    float         MetInfo_Pt;
    float         MetInfo_Phi;
    float         MetInfo_Px;
    float         MetInfo_Py;
    float         MetInfo_SumET;
};
class myParameters{
public:
    myParameters();
    ~myParameters();
    Int_t EvtInfo_totalEntry_before_preselection;
    float EvtInfo_NormalizationFactor_lumi;
    //------------------------
    float DiPhoInfo_eta;
    float DiPhoInfo_phi;
    float DiPhoInfo_energy;
    //------------------------
    Int_t num_jets;// # of selected objects.
    std::vector<float> JetInfo_jet_pt;
    std::vector<float> JetInfo_jet_eta;
    std::vector<float> JetInfo_jet_phi;
    std::vector<float> JetInfo_jet_energy;
    std::vector<float> JetInfo_jet_diphoton_deltaR;
    std::vector<float> JetInfo_jet_leadingPhoton_deltaR;
    std::vector<float> JetInfo_jet_subleadingPhoton_deltaR;
    std::vector<float> JetInfo_jet_pfDeepCSVJetTags_probb;
    std::vector<float> JetInfo_jet_pfDeepCSVJetTags_probbb;
    //------------------------
    Int_t num_leptons;// # of selected objects.
    Int_t num_electrons;// # of selected objects.
    Int_t num_muons;// # of selected objects.
    std::vector<float> ElecInfo_electron_pt;
    std::vector<float> ElecInfo_electron_eta;
    std::vector<float> ElecInfo_electron_phi;
    std::vector<float> ElecInfo_electron_energy;
    std::vector<float> ElecInfo_electron_diphoton_deltaR;
    std::vector<float> ElecInfo_electron_leadingPhoton_deltaR;
    std::vector<float> ElecInfo_electron_subleadingPhoton_deltaR;
    std::vector<float> MuonInfo_muon_pt;
    std::vector<float> MuonInfo_muon_eta;
    std::vector<float> MuonInfo_muon_phi;
    std::vector<float> MuonInfo_muon_energy;
    std::vector<float> MuonInfo_muon_diphoton_deltaR;
    std::vector<float> MuonInfo_muon_leadingPhoton_deltaR;
    std::vector<float> MuonInfo_muon_subleadingPhoton_deltaR;
    //------------------------
    //Not used in preselection stage
    //------------------------
    Int_t num_btagged_jets;// # of selected objects.
    Int_t num_nonbtagged_jets;// # of selected objects.
    //------------------------
    //Chi-2 sorting related
    //------------------------
    float Mass_w_candidate_chi2_simple;
    float Mass_top_candidate_chi2_simple;
    float Mass_w_candidate_chi2_modified;
    float Mass_top_candidate_chi2_modified;
    float Mass_gen_w_candidate_chi2_simple;
    float Mass_gen_top_candidate_chi2_simple;
    float Mass_gen_w_candidate_chi2_modified;
    float Mass_gen_top_candidate_chi2_modified;
    float inv_mass_dijet;
    float inv_mass_diphoton;
    float inv_mass_tbw;
    //------------------------
    float JetInfo_dijet_delta_eta;
    float JetInfo_dijet_delta_phi;
    float JetInfo_dijet_delta_angle;
    //------------------------
    //containers used in selection stage
    //------------------------
    std::vector<float> *JetInfo_jet_pt_selection;
    std::vector<float> *JetInfo_jet_eta_selection;
    std::vector<float> *JetInfo_jet_phi_selection;
    std::vector<float> *JetInfo_jet_energy_selection;
    std::vector<float> *JetInfo_jet_diphoton_deltaR_selection;
    std::vector<float> *JetInfo_jet_leadingPhoton_deltaR_selection;
    std::vector<float> *JetInfo_jet_subleadingPhoton_deltaR_selection;
    std::vector<float> *JetInfo_jet_pfDeepCSVJetTags_probb_selection;
    std::vector<float> *JetInfo_jet_pfDeepCSVJetTags_probbb_selection;
    std::vector<float> *ElecInfo_electron_pt_selection;
    std::vector<float> *ElecInfo_electron_eta_selection;
    std::vector<float> *ElecInfo_electron_phi_selection;
    std::vector<float> *ElecInfo_electron_energy_selection;
    std::vector<float> *ElecInfo_electron_diphoton_deltaR_selection;
    std::vector<float> *ElecInfo_electron_leadingPhoton_deltaR_selection;
    std::vector<float> *ElecInfo_electron_subleadingPhoton_deltaR_selection;
    std::vector<float> *MuonInfo_muon_pt_selection;
    std::vector<float> *MuonInfo_muon_eta_selection;
    std::vector<float> *MuonInfo_muon_phi_selection;
    std::vector<float> *MuonInfo_muon_energy_selection;
    std::vector<float> *MuonInfo_muon_diphoton_deltaR_selection;
    std::vector<float> *MuonInfo_muon_leadingPhoton_deltaR_selection;
    std::vector<float> *MuonInfo_muon_subleadingPhoton_deltaR_selection;
    void Clear();
};
class flashggStdTreeReader: public flashggStdTreeParameters{
public:
    flashggStdTreeReader();
    TChain *flashggStdTree;
    void InitChain(const char* treeName);
    void AddSingleRootFile(char* input_file);
    void AddMultiRootFile(char* input_file);
    void SetBranchAddresses(void);
    int GetEntries(void);
    double GetGenWeight(void);
    TChain *GetTChain(void);
};
class myTreeClass: public myParameters, public flashggStdTreeParameters{
public:
    TTree *mytree;
    //For creating new one
    void InitTree(void);
    void MakeNewBranchAddresses(void);
    void Fill(void);
    //For reading out
    void InitTree(const char* treeName);
    void AddRootFile(TFile* input);
    void SetBranchAddresses(void);
    int GetEntries(void);
    TTree *GetTTree(void);
};
//}}}
#endif
