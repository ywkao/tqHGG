//Opening{{{
/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides a simple example on how to use the trained classifiers
/// within an analysis module
/// - Project    : TMVA - a Root-integrated toolkit for multivariate data analysis
/// - Package    : TMVA
/// - Exectuable: TMVAClassificationApplication
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;
//}}}

void TMVAClassificationApplication( TString myMethodList = "" )
{
    // # Begining: set up MVA methods{{{
    // This loads the library
    TMVA::Tools::Instance();
    // Default MVA methods to be trained + tested
    std::map<std::string,int> Use;
    //MVA method
    // Cut optimisation
    Use["Cuts"]            = 0;
    Use["CutsD"]           = 0;
    Use["CutsPCA"]         = 0;
    Use["CutsGA"]          = 0;
    Use["CutsSA"]          = 0;
    //
    // 1-dimensional likelihood ("naive Bayes estimator")
    Use["Likelihood"]      = 0;
    Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
    Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
    Use["LikelihoodKDE"]   = 0;
    Use["LikelihoodMIX"]   = 0;
    //
    // Mutidimensional likelihood and Nearest-Neighbour methods
    Use["PDERS"]           = 0;
    Use["PDERSD"]          = 0;
    Use["PDERSPCA"]        = 0;
    Use["PDEFoam"]         = 0;
    Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
    Use["KNN"]             = 0; // k-nearest neighbour method
    //
    // Linear Discriminant Analysis
    Use["LD"]              = 0; // Linear Discriminant identical to Fisher
    Use["Fisher"]          = 0;
    Use["FisherG"]         = 0;
    Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
    Use["HMatrix"]         = 0;
    //
    // Function Discriminant analysis
    Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
    Use["FDA_SA"]          = 0;
    Use["FDA_MC"]          = 0;
    Use["FDA_MT"]          = 0;
    Use["FDA_GAMT"]        = 0;
    Use["FDA_MCMT"]        = 0;
    //
    // Neural Networks (all are feed-forward Multilayer Perceptrons)
    Use["MLP"]             = 0; // Recommended ANN
    Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
    Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
    Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
    Use["TMlpANN"]         = 0; // ROOT's own ANN
    Use["DNN_CPU"] = 0;         // CUDA-accelerated DNN training.
    Use["DNN_GPU"] = 0;         // Multi-core accelerated DNN.
    //
    // Support Vector Machine
    Use["SVM"]             = 0;
    //
    // Boosted Decision Trees
    Use["BDT"]             = 1; // uses Adaptive Boost
    Use["BDTG"]            = 0; // uses Gradient Boost
    Use["BDTB"]            = 0; // uses Bagging
    Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
    Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
    //
    // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
    Use["RuleFit"]         = 0;
    // ---------------------------------------------------------------
    Use["Plugin"]          = 0;
    Use["Category"]        = 0;
    Use["SVM_Gauss"]       = 0;
    Use["SVM_Poly"]        = 0;
    Use["SVM_Lin"]         = 0;
    std::cout << std::endl;
    std::cout << "==> Start TMVAClassificationApplication" << std::endl;
    // Select methods (don't look at this code - not of interest)
    if (myMethodList != "") {
       for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

       std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
       for (UInt_t i=0; i<mlist.size(); i++) {
          std::string regMethod(mlist[i]);

          if (Use.find(regMethod) == Use.end()) {
             std::cout << "Method \"" << regMethod
                       << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
             for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
                std::cout << it->first << " ";
             }
             std::cout << std::endl;
             return;
          }
          Use[regMethod] = 1;
       }
    }
    //}}}

    const char Directory[64]  = "./plots_leptonic_latest/161718/mva";
    const char signal_dir[64] = "./plots_leptonic_latest/2017old/mva";
    const char dataset[256] = "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8";
    const char tag[64]        = "testAll_24_both_hut";
    //TString I/O{{{
    TString dir               = Form("%s/dataset_%s/weights/", Directory, tag); // booking of MVA method
    TString fname = Form("%s/tree_%s.root", signal_dir, dataset);
    TString output_fname      = Form("%s/output_score_%s.root", Directory, dataset); // output file
    //}}}
    // # Create the Reader object{{{
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
    //variables{{{
    Float_t tree_event_weight;
    Float_t tree_diphoton_mass;
    Float_t tree_tH_deltaR;
    Float_t tree_top_eta;
    Float_t tree_top_pt;
    Float_t tree_top_mass;
    Float_t tree_neutrino_pz;
    Float_t tree_met;
    Float_t tree_jet4_CvsB;
    Float_t tree_jet4_CvsL;
    Float_t tree_jet4_btag;
    Float_t tree_jet4_eta;
    Float_t tree_jet4_pt;
    Float_t tree_jet3_CvsB;
    Float_t tree_jet3_CvsL;
    Float_t tree_jet3_btag;
    Float_t tree_jet3_eta;
    Float_t tree_jet3_pt;
    Float_t tree_jet2_CvsB;
    Float_t tree_jet2_CvsL;
    Float_t tree_jet2_btag;
    Float_t tree_jet2_eta;
    Float_t tree_jet2_pt;
    Float_t tree_jet1_CvsB;
    Float_t tree_jet1_CvsL;
    Float_t tree_jet1_btag;
    Float_t tree_jet1_eta;
    Float_t tree_jet1_pt;
    Float_t tree_num_jets;
    Float_t tree_lepton_diphoton_deltaTheta;
    Float_t tree_lepton_subleadingPhoton_deltaR;
    Float_t tree_lepton_leadingPhoton_deltaR;
    Float_t tree_lepton_charge; // int -> float
    Float_t tree_lepton_eta;
    Float_t tree_lepton_pt;
    Float_t tree_diphoton_cos_deltaPhi;
    Float_t tree_diphoton_eta;
    Float_t tree_diphoton_pt;
    Float_t tree_min_IDMVA;
    Float_t tree_max_IDMVA;
    Float_t tree_subleadingPhoton_hasPixelSeed; // bool -> float
    Float_t tree_subleadingPhoton_IDMVA;
    Float_t tree_subleadingPhoton_eta;
    Float_t tree_subleadingPhoton_pt_overM;
    Float_t tree_leadingPhoton_hasPixelSeed; // bool -> float
    Float_t tree_leadingPhoton_IDMVA;
    Float_t tree_leadingPhoton_eta;
    Float_t tree_leadingPhoton_pt_overM;
    //}}}
    reader->AddVariable( "tree_leadingPhoton_pt_overM", &tree_leadingPhoton_pt_overM);
    reader->AddVariable( "tree_leadingPhoton_eta", &tree_leadingPhoton_eta);
    reader->AddVariable( "tree_leadingPhoton_IDMVA", &tree_leadingPhoton_IDMVA);
    reader->AddVariable( "tree_leadingPhoton_hasPixelSeed", &tree_leadingPhoton_hasPixelSeed); // bool -> float
    reader->AddVariable( "tree_subleadingPhoton_pt_overM", &tree_subleadingPhoton_pt_overM);
    reader->AddVariable( "tree_subleadingPhoton_eta", &tree_subleadingPhoton_eta);
    reader->AddVariable( "tree_subleadingPhoton_IDMVA", &tree_subleadingPhoton_IDMVA);
    reader->AddVariable( "tree_subleadingPhoton_hasPixelSeed", &tree_subleadingPhoton_hasPixelSeed); // bool -> float
    reader->AddVariable( "tree_max_IDMVA", &tree_max_IDMVA);
    reader->AddVariable( "tree_min_IDMVA", &tree_min_IDMVA);
    reader->AddVariable( "tree_diphoton_pt", &tree_diphoton_pt);
    reader->AddVariable( "tree_diphoton_eta", &tree_diphoton_eta);
    reader->AddVariable( "tree_diphoton_cos_deltaPhi", &tree_diphoton_cos_deltaPhi);
    reader->AddVariable( "tree_lepton_pt", &tree_lepton_pt);
    reader->AddVariable( "tree_lepton_eta", &tree_lepton_eta);
    reader->AddVariable( "tree_lepton_charge", &tree_lepton_charge); // int -> float
    reader->AddVariable( "tree_lepton_leadingPhoton_deltaR", &tree_lepton_leadingPhoton_deltaR);
    reader->AddVariable( "tree_lepton_subleadingPhoton_deltaR", &tree_lepton_subleadingPhoton_deltaR);
    reader->AddVariable( "tree_lepton_diphoton_deltaTheta", &tree_lepton_diphoton_deltaTheta);
    reader->AddVariable( "tree_num_jets", &tree_num_jets);
    reader->AddVariable( "tree_jet1_pt", &tree_jet1_pt);
    reader->AddVariable( "tree_jet1_eta", &tree_jet1_eta);
    reader->AddVariable( "tree_jet1_btag", &tree_jet1_btag);
    reader->AddVariable( "tree_jet1_CvsL", &tree_jet1_CvsL);
    reader->AddVariable( "tree_jet1_CvsB", &tree_jet1_CvsB);
    reader->AddVariable( "tree_jet2_pt", &tree_jet2_pt);
    reader->AddVariable( "tree_jet2_eta", &tree_jet2_eta);
    reader->AddVariable( "tree_jet2_btag", &tree_jet2_btag);
    reader->AddVariable( "tree_jet2_CvsL", &tree_jet2_CvsL);
    reader->AddVariable( "tree_jet2_CvsB", &tree_jet2_CvsB);
    //reader->AddVariable( "tree_jet3_pt", &tree_jet3_pt);
    //reader->AddVariable( "tree_jet3_eta", &tree_jet3_eta);
    //reader->AddVariable( "tree_jet3_btag", &tree_jet3_btag);
    //reader->AddVariable( "tree_jet3_CvsL", &tree_jet3_CvsL);
    //reader->AddVariable( "tree_jet3_CvsB", &tree_jet3_CvsB);
    //reader->AddVariable( "tree_jet4_pt", &tree_jet4_pt);
    //reader->AddVariable( "tree_jet4_eta", &tree_jet4_eta);
    //reader->AddVariable( "tree_jet4_btag", &tree_jet4_btag);
    //reader->AddVariable( "tree_jet4_CvsL", &tree_jet4_CvsL);
    //reader->AddVariable( "tree_jet4_CvsB", &tree_jet4_CvsB);
    reader->AddVariable( "tree_met", &tree_met);
    reader->AddVariable( "tree_neutrino_pz", &tree_neutrino_pz);
    reader->AddVariable( "tree_top_mass", &tree_top_mass);
    reader->AddVariable( "tree_top_pt", &tree_top_pt);
    reader->AddVariable( "tree_top_eta", &tree_top_eta);
    reader->AddVariable( "tree_tH_deltaR", &tree_tH_deltaR);
    //}}}
    //Book the MVA methods{{{
    TString prefix = "TMVAClassification";
    // Book method(s)
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
       if (it->second) {
          TString methodName = TString(it->first) + TString(" method");
          TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
          reader->BookMVA( methodName, weightfile );
       }
    }
    //}}}
    // Book output histograms{{{
    UInt_t nbin = 100;
    TH1F *histLk(0);
    TH1F *histLkD(0);
    TH1F *histLkPCA(0);
    TH1F *histLkKDE(0);
    TH1F *histLkMIX(0);
    TH1F *histPD(0);
    TH1F *histPDD(0);
    TH1F *histPDPCA(0);
    TH1F *histPDEFoam(0);
    TH1F *histPDEFoamErr(0);
    TH1F *histPDEFoamSig(0);
    TH1F *histKNN(0);
    TH1F *histHm(0);
    TH1F *histFi(0);
    TH1F *histFiG(0);
    TH1F *histFiB(0);
    TH1F *histLD(0);
    TH1F *histNn(0);
    TH1F *histNnbfgs(0);
    TH1F *histNnbnn(0);
    TH1F *histNnC(0);
    TH1F *histNnT(0);
    TH1F *histBdt(0);
    TH1F *histBdt_signalRegion(0);
    TH1F *histBdt_signalRegion_final(0);
    TH1F *histBdtG(0);
    TH1F *histBdtB(0);
    TH1F *histBdtD(0);
    TH1F *histBdtF(0);
    TH1F *histRf(0);
    TH1F *histSVMG(0);
    TH1F *histSVMP(0);
    TH1F *histSVML(0);
    TH1F *histFDAMT(0);
    TH1F *histFDAGA(0);
    TH1F *histCat(0);
    TH1F *histPBdt(0);
    TH1F *histDnnGpu(0);
    TH1F *histDnnCpu(0);

    if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
    if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
    if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
    if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
    if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
    if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
    if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
    if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
    if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
    if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
    if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
    if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
    if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
    if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
    if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
    if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
    if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
    if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
    if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
    if (Use["DNN_GPU"]) histDnnGpu = new TH1F("MVA_DNN_GPU", "MVA_DNN_GPU", nbin, -0.1, 1.1);
    if (Use["DNN_CPU"]) histDnnCpu = new TH1F("MVA_DNN_CPU", "MVA_DNN_CPU", nbin, -0.1, 1.1);
    if (Use["BDT"])           histBdt                    = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
    if (Use["BDT"])           histBdt_signalRegion       = new TH1F( "MVA_BDT_signalRegion",           "MVA_BDT_signalRegion",           nbin, -0.8, 0.8 );
    if (Use["BDT"])           histBdt_signalRegion_final = new TH1F( "MVA_BDT_signalRegion_final",           "MVA_BDT_signalRegion_final",           nbin, -0.8, 0.8 );
    if (Use["BDTG"])          histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );
    if (Use["BDTB"])          histBdtB    = new TH1F( "MVA_BDTB",          "MVA_BDTB",          nbin, -1.0, 1.0 );
    if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
    if (Use["BDTF"])          histBdtF    = new TH1F( "MVA_BDTF",          "MVA_BDTF",          nbin, -1.0, 1.0 );
    if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
    if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
    if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
    if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
    if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
    if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
    if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
    if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

    // PDEFoam also returns per-event error, fill in histogram, and also fill significance
    if (Use["PDEFoam"]) {
       histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
       histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
       histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
    }

    // Book example histogram for probability (the other methods are done similarly)
    TH1F *probHistFi(0), *rarityHistFi(0);
    if (Use["Fisher"]) {
       probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
       rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
    }
    //}}}
    //Prepare input tree {{{
    TFile *input(0);
    if (!gSystem->AccessPathName( fname )) {
       input = TFile::Open( fname ); // check if file in local directory exists
    }
    else {
       std::cout << "--- Warning: file does not exists!" << std::endl;
    }
    std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;
    //}}}
    // # Prepare the event tree{{{
    std::cout << "--- Select signal sample" << std::endl;
    TTree* theTree = (TTree*)input->Get("myAnalysisTree");
    theTree->SetBranchAddress( "tree_event_weight", &tree_event_weight );
    theTree->SetBranchAddress( "tree_diphoton_mass", &tree_diphoton_mass );
         //float NormalizationFactor = treeReader.EvtInfo_genweight * treeReader.EvtInfo_NormalizationFactor_lumi * PU_reweighting_factor;
         //tree_event_weight = NormalizationFactor;
         //Reminder: EvtInfo_NormalizationFactor_lumi = 1000. * Luminosity * CrossSection * BranchingFraction / TotalGenweight;
    theTree->SetBranchAddress( "tree_leadingPhoton_pt_overM", &tree_leadingPhoton_pt_overM );
    theTree->SetBranchAddress( "tree_leadingPhoton_eta", &tree_leadingPhoton_eta );
    theTree->SetBranchAddress( "tree_leadingPhoton_IDMVA", &tree_leadingPhoton_IDMVA );
    theTree->SetBranchAddress( "tree_leadingPhoton_hasPixelSeed", &tree_leadingPhoton_hasPixelSeed ); // bool -> float
    theTree->SetBranchAddress( "tree_subleadingPhoton_pt_overM", &tree_subleadingPhoton_pt_overM );
    theTree->SetBranchAddress( "tree_subleadingPhoton_eta", &tree_subleadingPhoton_eta );
    theTree->SetBranchAddress( "tree_subleadingPhoton_IDMVA", &tree_subleadingPhoton_IDMVA );
    theTree->SetBranchAddress( "tree_subleadingPhoton_hasPixelSeed", &tree_subleadingPhoton_hasPixelSeed ); // bool -> float
    theTree->SetBranchAddress( "tree_max_IDMVA", &tree_max_IDMVA );
    theTree->SetBranchAddress( "tree_min_IDMVA", &tree_min_IDMVA );
    theTree->SetBranchAddress( "tree_diphoton_pt", &tree_diphoton_pt );
    theTree->SetBranchAddress( "tree_diphoton_eta", &tree_diphoton_eta );
    theTree->SetBranchAddress( "tree_diphoton_cos_deltaPhi", &tree_diphoton_cos_deltaPhi );
    theTree->SetBranchAddress( "tree_lepton_pt", &tree_lepton_pt );
    theTree->SetBranchAddress( "tree_lepton_eta", &tree_lepton_eta );
    theTree->SetBranchAddress( "tree_lepton_charge", &tree_lepton_charge ); // int -> float
    theTree->SetBranchAddress( "tree_lepton_leadingPhoton_deltaR", &tree_lepton_leadingPhoton_deltaR );
    theTree->SetBranchAddress( "tree_lepton_subleadingPhoton_deltaR", &tree_lepton_subleadingPhoton_deltaR );
    theTree->SetBranchAddress( "tree_lepton_diphoton_deltaTheta", &tree_lepton_diphoton_deltaTheta );
    theTree->SetBranchAddress( "tree_num_jets", &tree_num_jets );
    theTree->SetBranchAddress( "tree_jet1_pt", &tree_jet1_pt );
    theTree->SetBranchAddress( "tree_jet1_eta", &tree_jet1_eta );
    theTree->SetBranchAddress( "tree_jet1_btag", &tree_jet1_btag );
    theTree->SetBranchAddress( "tree_jet1_CvsL", &tree_jet1_CvsL );
    theTree->SetBranchAddress( "tree_jet1_CvsB", &tree_jet1_CvsB );
    theTree->SetBranchAddress( "tree_jet2_pt", &tree_jet2_pt );
    theTree->SetBranchAddress( "tree_jet2_eta", &tree_jet2_eta );
    theTree->SetBranchAddress( "tree_jet2_btag", &tree_jet2_btag );
    theTree->SetBranchAddress( "tree_jet2_CvsL", &tree_jet2_CvsL );
    theTree->SetBranchAddress( "tree_jet2_CvsB", &tree_jet2_CvsB );
    //theTree->SetBranchAddress( "tree_jet3_pt", &tree_jet3_pt );
    //theTree->SetBranchAddress( "tree_jet3_eta", &tree_jet3_eta );
    //theTree->SetBranchAddress( "tree_jet3_btag", &tree_jet3_btag );
    //theTree->SetBranchAddress( "tree_jet3_CvsL", &tree_jet3_CvsL );
    //theTree->SetBranchAddress( "tree_jet3_CvsB", &tree_jet3_CvsB );
    //theTree->SetBranchAddress( "tree_jet4_pt", &tree_jet4_pt );
    //theTree->SetBranchAddress( "tree_jet4_eta", &tree_jet4_eta );
    //theTree->SetBranchAddress( "tree_jet4_btag", &tree_jet4_btag );
    //theTree->SetBranchAddress( "tree_jet4_CvsL", &tree_jet4_CvsL );
    //theTree->SetBranchAddress( "tree_jet4_CvsB", &tree_jet4_CvsB );
    theTree->SetBranchAddress( "tree_met", &tree_met );
    theTree->SetBranchAddress( "tree_neutrino_pz", &tree_neutrino_pz );
    theTree->SetBranchAddress( "tree_top_mass", &tree_top_mass );
    theTree->SetBranchAddress( "tree_top_pt", &tree_top_pt );
    theTree->SetBranchAddress( "tree_top_eta", &tree_top_eta );
    theTree->SetBranchAddress( "tree_tH_deltaR", &tree_tH_deltaR );
    //}}}
    int   n_final_entry = 0;
    float n_final_yield = 0;
    float FINAL_SELECTION = 0.184; //25_st_hut, determined by ZA
    // Event loop{{{
    // Efficiency calculator for cut method
    Int_t    nSelCutsGA = 0;
    Double_t effS       = 0.7;
    std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
    TStopwatch sw;
    sw.Start();
    for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
       if (ievt%2000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
       theTree->GetEntry(ievt);
       bool pass_signal_region = (tree_diphoton_mass > 120 && tree_diphoton_mass < 130);

       double mva_score = reader->EvaluateMVA( "BDT method" );
       histBdt->Fill( mva_score , tree_event_weight);
       if(pass_signal_region) histBdt_signalRegion->Fill( mva_score , tree_event_weight);
       if(pass_signal_region && mva_score > FINAL_SELECTION ){
            histBdt_signalRegion_final->Fill( mva_score , tree_event_weight);
            n_final_entry += 1;
            n_final_yield += tree_event_weight;
       }
       // Return the MVA outputs and fill into histograms (skipped){{{
       /*
       if (Use["CutsGA"]) {
          // Cuts is a special case: give the desired signal efficienciy
          Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
          if (passed) nSelCutsGA++;
       }
       if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) , tree_event_weight);
       if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) , tree_event_weight);
       if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) , tree_event_weight);
       if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) , tree_event_weight);
       if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) , tree_event_weight);
       if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) , tree_event_weight);
       if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) , tree_event_weight);
       if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) , tree_event_weight);
       if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) , tree_event_weight);
       if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) , tree_event_weight);
       if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) , tree_event_weight);
       if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) , tree_event_weight);
       if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) , tree_event_weight);
       if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) , tree_event_weight);
       if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) , tree_event_weight);
       if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) , tree_event_weight);
       if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) , tree_event_weight);
       if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) , tree_event_weight);
       if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) , tree_event_weight);
       if (Use["DNN_GPU"]) histDnnGpu->Fill(reader->EvaluateMVA("DNN_GPU method"), tree_event_weight);
       if (Use["DNN_CPU"]) histDnnCpu->Fill(reader->EvaluateMVA("DNN_CPU method"), tree_event_weight);
       if (Use["BDT"          ])   histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) , tree_event_weight);
       if (Use["BDTG"         ])   histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"          ) , tree_event_weight);
       if (Use["BDTB"         ])   histBdtB   ->Fill( reader->EvaluateMVA( "BDTB method"          ) , tree_event_weight);
       if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) , tree_event_weight);
       if (Use["BDTF"         ])   histBdtF   ->Fill( reader->EvaluateMVA( "BDTF method"          ) , tree_event_weight);
       if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) , tree_event_weight);
       if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) , tree_event_weight);
       if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) , tree_event_weight);
       if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) , tree_event_weight);
       if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) , tree_event_weight);
       if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) , tree_event_weight);
       if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) , tree_event_weight);
       if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) , tree_event_weight);

       // Retrieve also per-event error
       if (Use["PDEFoam"]) {
          Double_t val = reader->EvaluateMVA( "PDEFoam method" );
          Double_t err = reader->GetMVAError();
          histPDEFoam   ->Fill( val , tree_event_weight);
          histPDEFoamErr->Fill( err , tree_event_weight);
          if (err>1.e-50) histPDEFoamSig->Fill( val/err , tree_event_weight);
       }

       // Retrieve probability instead of MVA output
       if (Use["Fisher"])   {
          probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) , tree_event_weight);
          rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) , tree_event_weight);
       }
       */
       //}}}
    }
    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();
    //}}}
    // Evaluate yields{{{
    int nbins = histBdt_signalRegion_final->GetNbinsX();
    printf("[INFO] nbins = %d\n", nbins);
    float s = 0, unc_s = 0;
    for(int i=0; i<nbins; ++i){
        s += histBdt_signalRegion_final->GetBinContent(i+1);
        unc_s += pow(histBdt_signalRegion_final->GetBinError(i+1), 2);
    }
    unc_s = sqrt(unc_s);
    //}}}
    printf("[INFO] n_final_entry = %d\n", n_final_entry);
    printf("[INFO] n_final_yield = %f\n", n_final_yield);
    printf("[INFO] yields (hist) = %.2f #pm %.2f \n", s, unc_s);
    // Get efficiency for cuts classifier{{{
    if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                 << " (for a required signal efficiency of " << effS << ")" << std::endl;

    if (Use["CutsGA"]) {

       // test: retrieve cuts for particular signal efficiency
       // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer
       TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

       if (mcuts) {
          std::vector<Double_t> cutsMin;
          std::vector<Double_t> cutsMax;
          mcuts->GetCuts( 0.7, cutsMin, cutsMax );
          std::cout << "--- -------------------------------------------------------------" << std::endl;
          std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
          for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
             std::cout << "... Cut: "
                       << cutsMin[ivar]
                       << " < \""
                       << mcuts->GetInputVar(ivar)
                       << "\" <= "
                       << cutsMax[ivar] << std::endl;
          }
          std::cout << "--- -------------------------------------------------------------" << std::endl;
       }
    }
    //}}}
    // Write histograms{{{
    TFile *target  = new TFile( output_fname, "RECREATE" );
    histBdt                    -> Write();
    histBdt_signalRegion       -> Write();
    histBdt_signalRegion_final -> Write();
    //}}}
    // Write histograms (skipped){{{
    /*
    if (Use["Likelihood"   ])   histLk     ->Write();
    if (Use["LikelihoodD"  ])   histLkD    ->Write();
    if (Use["LikelihoodPCA"])   histLkPCA  ->Write();
    if (Use["LikelihoodKDE"])   histLkKDE  ->Write();
    if (Use["LikelihoodMIX"])   histLkMIX  ->Write();
    if (Use["PDERS"        ])   histPD     ->Write();
    if (Use["PDERSD"       ])   histPDD    ->Write();
    if (Use["PDERSPCA"     ])   histPDPCA  ->Write();
    if (Use["KNN"          ])   histKNN    ->Write();
    if (Use["HMatrix"      ])   histHm     ->Write();
    if (Use["Fisher"       ])   histFi     ->Write();
    if (Use["FisherG"      ])   histFiG    ->Write();
    if (Use["BoostedFisher"])   histFiB    ->Write();
    if (Use["LD"           ])   histLD     ->Write();
    if (Use["MLP"          ])   histNn     ->Write();
    if (Use["MLPBFGS"      ])   histNnbfgs ->Write();
    if (Use["MLPBNN"       ])   histNnbnn  ->Write();
    if (Use["CFMlpANN"     ])   histNnC    ->Write();
    if (Use["TMlpANN"      ])   histNnT    ->Write();
    if (Use["DNN_GPU"]) histDnnGpu->Write();
    if (Use["DNN_CPU"]) histDnnCpu->Write();
    if (Use["BDT"          ])   histBdt    ->Write();
    if (Use["BDTG"         ])   histBdtG   ->Write();
    if (Use["BDTB"         ])   histBdtB   ->Write();
    if (Use["BDTD"         ])   histBdtD   ->Write();
    if (Use["BDTF"         ])   histBdtF   ->Write();
    if (Use["RuleFit"      ])   histRf     ->Write();
    if (Use["SVM_Gauss"    ])   histSVMG   ->Write();
    if (Use["SVM_Poly"     ])   histSVMP   ->Write();
    if (Use["SVM_Lin"      ])   histSVML   ->Write();
    if (Use["FDA_MT"       ])   histFDAMT  ->Write();
    if (Use["FDA_GA"       ])   histFDAGA  ->Write();
    if (Use["Category"     ])   histCat    ->Write();
    if (Use["Plugin"       ])   histPBdt   ->Write();

    // Write also error and significance histos
    if (Use["PDEFoam"]) { histPDEFoam->Write(); histPDEFoamErr->Write(); histPDEFoamSig->Write(); }

    // Write also probability hists
    if (Use["Fisher"]) { if (probHistFi != 0) probHistFi->Write(); if (rarityHistFi != 0) rarityHistFi->Write(); }
    */
    //}}}
    // Close{{{
    target->Close();
    std::cout << Form("--- Output directory : \"%s\"", Directory) << std::endl;
    std::cout << Form("--- Created root file: \"output_score_%s.root\" containing the MVA output histograms", dataset) << std::endl;
    delete reader;
    std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
    //}}}
}

int main( int argc, char** argv )
{
    TString methodList;
    for (int i=1; i<argc; i++) {
       TString regMethod(argv[i]);
       if(regMethod=="-b" || regMethod=="--batch") continue;
       if (!methodList.IsNull()) methodList += TString(",");
       methodList += regMethod;
    }
    TMVAClassificationApplication(methodList);
    return 0;
}
