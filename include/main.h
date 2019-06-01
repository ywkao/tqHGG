#ifndef __MAIN_H__
#define __MAIN_H__

double Chi2_calculator(double w_mass, double t_mass);
double Chi2_calculator_w_only(double w_mass);
void MakePlots(TCanvas *c1, TH1D* hist, const char* title, const char* outputFile);
bool isThisDataOrNot(char* dataset);
bool isThisMCsignal(char* dataset);

class flashggStdTreeParameters{
public:
    flashggStdTreeParameters();
    ~flashggStdTreeParameters();
    float EvtInfo_NPu;
    Int_t EvtInfo_NVtx;
    float EvtInfo_Rho;
    Int_t jets_size;
    std::vector<float> *JetInfo_Pt;
    std::vector<float> *JetInfo_Eta;
    std::vector<float> *JetInfo_Phi;
    std::vector<float> *JetInfo_Mass;
    std::vector<float> *JetInfo_Energy;
    std::vector<float> *JetInfo_pfDeepCSVJetTags_probb;
    std::vector<float> *JetInfo_pfDeepCSVJetTags_probbb;
    std::vector<float> *JetInfo_pfCombinedInclusiveSecondaryVertexV2BJetTags;
    float EvtInfo_genweight;
    float DiPhoInfo_mass;
    float DiPhoInfo_leadPt;
    float DiPhoInfo_leadEta;
    float DiPhoInfo_leadPhi;
    float DiPhoInfo_leadE;
    float DiPhoInfo_leadIDMVA;
    float DiPhoInfo_subleadPt;
    float DiPhoInfo_subleadEta;
    float DiPhoInfo_subleadPhi;
    float DiPhoInfo_subleadE;
    float DiPhoInfo_subleadIDMVA;
};
class myParameters{
public:
    Int_t EvtInfo_totalEntry_before_preselection;
    float EvtInfo_NormalizationFactor_lumi;
    //------------------------
    Int_t num_jets;
    Int_t num_btagged_jets;
    Int_t num_nonbtagged_jets;
    //------------------------
    float inv_mass_dijet;
    float inv_mass_diphoton;
    float inv_mass_tbw;
    //------------------------
    float JetInfo_dijet_delta_eta;
    float JetInfo_dijet_delta_phi;
    float JetInfo_dijet_delta_angle;
    //------------------------
    float JetInfo_leading_bjet_pt;
    float JetInfo_leading_bjet_eta;
    float JetInfo_leading_bjet_phi;
    //------------------------
    Int_t JetInfo_chosen_bjet_is_leading_bjet;
    float JetInfo_chosen_bjet_pt;
    float JetInfo_chosen_bjet_eta;
    float JetInfo_chosen_bjet_phi;
    //------------------------
    float JetInfo_jet1_pt;
    float JetInfo_jet1_eta;
    float JetInfo_jet1_phi;
    //------------------------
    float JetInfo_jet2_pt;
    float JetInfo_jet2_eta;
    float JetInfo_jet2_phi;

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

#endif
