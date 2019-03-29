#ifndef __MAIN_H__
#define __MAIN_H__

double Chi2_calculator(double w_mass, double t_mass);
void MakePlots(TCanvas *c1, TH1D* hist, const char* title, const char* outputFile);
bool isThisDataOrNot(char* dataset);
bool isThisMultiFile(char* dataset);
bool isThisMCsignal(char* dataset);

class flashggStdTreeParameters{
    Int_t jets_size;
    std::vector<float> *JetInfo_Pt;
    std::vector<float> *JetInfo_Eta;
    std::vector<float> *JetInfo_Phi;
    std::vector<float> *JetInfo_Mass;
    std::vector<float> *JetInfo_Energy;
    std::vector<float> *JetInfo_pfDeepCSVJetTags_probb;
    std::vector<float> *JetInfo_pfDeepCSVJetTags_probbb;
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
    Int_t num_jets;
    Int_t num_btagged_jets;
    Int_t num_nonbtagged_jets;
    //------------------------
    float inv_mass_dijet;
    float inv_mass_diphoton;
    float inv_mass_tbw;
    float JetInfo_dijet_delta_eta;
    float JetInfo_dijet_delta_phi;
    float JetInfo_dijet_delta_angle;
    //------------------------
    std::vector<float> JetInfo_bjet_pt;
    std::vector<float> JetInfo_bjet_eta;
    std::vector<float> JetInfo_bjet_phi;
    std::vector<float> JetInfo_jet1_pt;
    std::vector<float> JetInfo_jet1_eta;
    std::vector<float> JetInfo_jet1_phi;
    std::vector<float> JetInfo_jet2_pt;
    std::vector<float> JetInfo_jet2_eta;
    std::vector<float> JetInfo_jet2_phi;

    void Clear();
};
class flashggStdTreeReader: public flashggStdTreeParameters{
public:
    flashggStdTreeReader();
    TChain *flashggStdTree;
    void InitChain(const char* treeName);
    void AddSingleRootFile(char* input_file);
    void AddMultiRootFile(char* input_file);
    void SetBranchAddresses();
};
class myTreeClass: public myParameters, public flashggStdTreeParameters{
    TTree *mytree;
public:
    void InitTree(const char* treeName);
    void SetBranchAddresses();
    void Fill();
};


void flashggStdTreeReader::InitChain(const char* treeName){
    printf("From flashggStdTreeReader: %s\n", treeName);
    flashggStdTree = new TChain(treeName);
}
void flashggStdTreeReader::AddSingleRootFile(char* input_file){
    flashggStdTree->Add(input_file);
}
void flashggStdTreeReader::AddMultiRootFile(char* input_file){
    flashggStdTree->Add(Form("%s/*.root", input_file));
}


void flashggStdTreeReader::SetBranchAddresses(){
    flashggStdTree->SetBranchAddress("EvtInfo.genweight", &EvtInfo_genweight);
    flashggStdTree->SetBranchAddress("jets_size", &jets_size);
    flashggStdTree->SetBranchAddress("JetInfo.Pt", &JetInfo_Pt);
    flashggStdTree->SetBranchAddress("JetInfo.Eta", &JetInfo_Eta);
    flashggStdTree->SetBranchAddress("JetInfo.Phi", &JetInfo_Phi);
    flashggStdTree->SetBranchAddress("JetInfo.Mass", &JetInfo_Mass);
    flashggStdTree->SetBranchAddress("JetInfo.Energy", &JetInfo_Energy);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probb", &JetInfo_pfDeepCSVJetTags_probb);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probbb", &JetInfo_pfDeepCSVJetTags_probbb);
    flashggStdTree->SetBranchAddress("DiPhoInfo.mass", &DiPhoInfo_mass);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadPt", &DiPhoInfo_leadPt);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadEta", &DiPhoInfo_leadEta);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadPhi", &DiPhoInfo_leadPhi);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadE", &DiPhoInfo_leadE);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadIDMVA", &DiPhoInfo_leadIDMVA);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadPt", &DiPhoInfo_subleadPt);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadEta", &DiPhoInfo_subleadEta);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadPhi", &DiPhoInfo_subleadPhi);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadE", &DiPhoInfo_subleadE);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadIDMVA", &DiPhoInfo_subleadIDMVA);
}
void myTreeClass::InitTree(const char* treeName){
    mytree =  new TTree("mytree", "mytree");
}
void myTreeClass::SetBranchAddresses(){
    mytree -> Branch("num_jets", &num_jets, "num_jets/I");
    mytree -> Branch("num_btagged_jets", &num_btagged_jets, "num_btagged_jets/I");
    mytree -> Branch("num_nonbtagged_jets", &num_nonbtagged_jets, "num_nonbtagged_jets/I");
    //------------------------
    mytree -> Branch("inv_mass_dijet", &inv_mass_dijet, "inv_mass_dijet/F");
    mytree -> Branch("inv_mass_diphoton", &inv_mass_diphoton, "inv_mass_diphoton/F");
    mytree -> Branch("inv_mass_tbw", &inv_mass_tbw, "inv_mass_tbw/F");
    //------------------------
    mytree -> Branch("JetInfo_bjet_pt", &JetInfo_bjet_pt);
    mytree -> Branch("JetInfo_bjet_eta", &JetInfo_bjet_eta);
    mytree -> Branch("JetInfo_bjet_phi", &JetInfo_bjet_phi);
    mytree -> Branch("JetInfo_jet1_pt", &JetInfo_jet1_pt);
    mytree -> Branch("JetInfo_jet1_eta", &JetInfo_jet1_eta);
    mytree -> Branch("JetInfo_jet1_phi", &JetInfo_jet1_phi);
    mytree -> Branch("JetInfo_jet2_pt", &JetInfo_jet2_pt);
    mytree -> Branch("JetInfo_jet2_eta", &JetInfo_jet2_eta);
    mytree -> Branch("JetInfo_jet2_phi", &JetInfo_jet2_phi);
    //------------------------
    mytree -> Branch("JetInfo_dijet_delta_eta", &JetInfo_dijet_delta_eta, "JetInfo_dijet_delta_eta/F");
    mytree -> Branch("JetInfo_dijet_delta_phi", &JetInfo_dijet_delta_phi, "JetInfo_dijet_delta_phi/F");
    mytree -> Branch("JetInfo_dijet_delta_angle", &JetInfo_dijet_delta_angle, "JetInfo_dijet_delta_angle/F");
    //------------------------
    mytree -> Branch("DiPhoInfo_leadPt", &DiPhoInfo_leadPt, "DiPhoInfo_leadPt/F");
    mytree -> Branch("DiPhoInfo_leadEta", &DiPhoInfo_leadEta, "DiPhoInfo_leadEta/F");
    mytree -> Branch("DiPhoInfo_leadPhi", &DiPhoInfo_leadPhi, "DiPhoInfo_leadPhi/F");
    mytree -> Branch("DiPhoInfo_leadE", &DiPhoInfo_leadE, "DiPhoInfo_leadE/F");
    mytree -> Branch("DiPhoInfo_leadIDMVA", &DiPhoInfo_leadIDMVA, "DiPhoInfo_leadIDMVA/F");
    mytree -> Branch("DiPhoInfo_subleadPt", &DiPhoInfo_subleadPt, "DiPhoInfo_subleadPt/F");
    mytree -> Branch("DiPhoInfo_subleadEta", &DiPhoInfo_subleadEta, "DiPhoInfo_subleadEta/F");
    mytree -> Branch("DiPhoInfo_subleadPhi", &DiPhoInfo_subleadPhi, "DiPhoInfo_subleadPhi/F");
    mytree -> Branch("DiPhoInfo_subleadE", &DiPhoInfo_subleadE, "DiPhoInfo_subleadE/F");
    mytree -> Branch("DiPhoInfo_subleadIDMVA", &DiPhoInfo_subleadIDMVA, "DiPhoInfo_subleadIDMVA/F");
    mytree -> Branch("DiPhoInfo_mass", &DiPhoInfo_mass, "DiPhoInfo_mass/F");
}
void myTreeClass::Fill(){
    mytree -> Fill();
}

void myParameters::Clear(){
    num_jets = 0;
    num_btagged_jets = 0;
    num_nonbtagged_jets = 0;
    //------------------------
    inv_mass_dijet = 0;
    inv_mass_diphoton = 0;
    inv_mass_tbw = 0;
    JetInfo_dijet_delta_eta = 0;
    JetInfo_dijet_delta_phi = 0;
    JetInfo_dijet_delta_angle = 0;
    //------------------------
    JetInfo_bjet_pt.clear();
    JetInfo_bjet_eta.clear();
    JetInfo_bjet_phi.clear();
    JetInfo_jet1_pt.clear();
    JetInfo_jet1_eta.clear();
    JetInfo_jet1_phi.clear();
    JetInfo_jet2_pt.clear();
    JetInfo_jet2_eta.clear();
    JetInfo_jet2_phi.clear();
}
#endif
