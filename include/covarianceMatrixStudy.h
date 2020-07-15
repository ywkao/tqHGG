#ifndef __COVARIANCEMATRIXSTUDY_H__
#define __COVARIANCEMATRIXSTUDY_H__
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TTree.h>
#include <vector>
#include <string>

using namespace std;

//==========================//
//---  Useful Constants  ---//
//==========================//
double pfDeepCSVJetTags_tight  = 0.8001;
double pfDeepCSVJetTags_medium = 0.4941;
double pfDeepCSVJetTags_loose  = 0.1522;

//==========================//
//---  Class & Function  ---//
//==========================//
int covarianceMatrixStudy(std::vector<char*> input_file, std::string tag);
int get_constraint_on_num_jets(std::string tag);
bool get_particle_non_void_check_result(string tag, std::vector<int> indices);
bool checkAvailability(int index, std::vector<int> ID_IsChosen);
double Chi2_calculator(double w_mass, double t_mass);
double Chi2_calculator_w_only(double w_mass);
void MakePlots(TCanvas *c1, TH1D* hist, const char* title, const char* outputFile);
bool isThisMultiFile(char* dataset);

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
};
class myParameters{
public:
    myParameters();
    ~myParameters();
    //------------------------
    double GenMatching_Mass_top_gen;
    double GenMatching_Mass_wboson_gen;
    double GenMatching_Mass_top_reco;
    double GenMatching_Mass_wboson_reco;
    //------------------------
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
    void AddSingleRootFile(std::string input_file);
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

#endif
