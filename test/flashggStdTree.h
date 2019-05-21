//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 28 14:12:54 2018 by ROOT version 5.34/32
// from TTree flashggStdTree/
// found on file: ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root
//////////////////////////////////////////////////////////

#ifndef flashggStdTree_h
#define flashggStdTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class flashggStdTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           EvtInfo_NPu;
   Int_t           EvtInfo_NVtx;
   Bool_t          EvtInfo_passTrigger;
   Float_t         EvtInfo_genweight;
   Float_t         EvtInfo_Rho;
   Float_t         EvtInfo_PVz;
   Float_t         EvtInfo_BSsigmaz;
   Bool_t          EvtInfo_Flag_HBHENoiseFilter;
   Bool_t          EvtInfo_Flag_HBHENoiseIsoFilter;
   Bool_t          EvtInfo_Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          EvtInfo_Flag_goodVertices;
   Bool_t          EvtInfo_Flag_globalSuperTightHalo2016Filter;
   Bool_t          EvtInfo_Flag_BadPFMuonFilter;
   Bool_t          EvtInfo_Flag_BadChargedCandidateFilter;
   Bool_t          EvtInfo_Flag_ecalBadCalibFilter;
   Bool_t          EvtInfo_Flag_eeBadScFilter;
   Float_t         DiPhoInfo_mass;
   Float_t         DiPhoInfo_pt;
   Float_t         DiPhoInfo_leadPt;
   Float_t         DiPhoInfo_leadEta;
   Float_t         DiPhoInfo_leadPhi;
   Float_t         DiPhoInfo_leadE;
   Float_t         DiPhoInfo_leadsigEOverE;
   Float_t         DiPhoInfo_leadR9;
   Float_t         DiPhoInfo_leadsieie;
   Float_t         DiPhoInfo_leadhoe;
   Float_t         DiPhoInfo_leadIDMVA;
   Bool_t          DiPhoInfo_leadhasPixelSeed;
   Bool_t          DiPhoInfo_leadGenMatch;
   Int_t           DiPhoInfo_leadGenMatchType;
   Float_t         DiPhoInfo_subleadPt;
   Float_t         DiPhoInfo_subleadEta;
   Float_t         DiPhoInfo_subleadPhi;
   Float_t         DiPhoInfo_subleadE;
   Float_t         DiPhoInfo_subleadsigEOverE;
   Float_t         DiPhoInfo_subleadR9;
   Float_t         DiPhoInfo_subleadsieie;
   Float_t         DiPhoInfo_subleadhoe;
   Float_t         DiPhoInfo_subleadIDMVA;
   Bool_t          DiPhoInfo_subleadhasPixelSeed;
   Bool_t          DiPhoInfo_subleadGenMatch;
   Int_t           DiPhoInfo_subleadGenMatchType;
   Float_t         DiPhoInfo_SelectedVz;
   Float_t         DiPhoInfo_GenVz;
   Float_t         DiPhoInfo_centralWeight;
   Float_t         DiPhoInfo_MvaLinearSystUp;
   Float_t         DiPhoInfo_MvaLinearSystDown;
   Float_t         DiPhoInfo_LooseMvaSFUp;
   Float_t         DiPhoInfo_LooseMvaSFDown;
   Float_t         DiPhoInfo_PreselSFUp;
   Float_t         DiPhoInfo_PreselSFDown;
   Float_t         DiPhoInfo_electronVetoSFUp;
   Float_t         DiPhoInfo_electronVetoSFDown;
   Float_t         DiPhoInfo_TriggerWeightUp;
   Float_t         DiPhoInfo_TriggerWeightDown;
   Float_t         DiPhoInfo_FracRVWeightUp;
   Float_t         DiPhoInfo_FracRVWeightDown;
   Float_t         DiPhoInfo_FracRVNvtxWeightUp;
   Float_t         DiPhoInfo_FracRVNvtxWeightDown;
   Int_t           ElecInfo_Size;
   vector<int>     *ElecInfo_Charge;
   vector<float>   *ElecInfo_Pt;
   vector<float>   *ElecInfo_Eta;
   vector<float>   *ElecInfo_Phi;
   vector<float>   *ElecInfo_Energy;
   vector<float>   *ElecInfo_EtaSC;
   vector<float>   *ElecInfo_PhiSC;
   vector<float>   *ElecInfo_GsfTrackDz;
   vector<float>   *ElecInfo_GsfTrackDxy;
   vector<bool>    *ElecInfo_EGMCutBasedIDVeto;
   vector<bool>    *ElecInfo_EGMCutBasedIDLoose;
   vector<bool>    *ElecInfo_EGMCutBasedIDMedium;
   vector<bool>    *ElecInfo_EGMCutBasedIDTight;
   vector<bool>    *ElecInfo_fggPhoVeto;
   vector<bool>    *ElecInfo_tmpPhoVeto;
   vector<float>   *ElecInfo_EnergyCorrFactor;
   vector<float>   *ElecInfo_EnergyPostCorrErr;
   vector<float>   *ElecInfo_EnergyPostCorrScaleUp;
   vector<float>   *ElecInfo_EnergyPostCorrScaleDown;
   vector<float>   *ElecInfo_EnergyPostCorrSmearUp;
   vector<float>   *ElecInfo_EnergyPostCorrSmearDown;
   vector<bool>    *ElecInfo_GenMatch;
   vector<int>     *ElecInfo_GenPdgID;
   vector<float>   *ElecInfo_GenPt;
   vector<float>   *ElecInfo_GenEta;
   vector<float>   *ElecInfo_GenPhi;
   Int_t           MuonInfo_Size;
   vector<int>     *MuonInfo_Charge;
   vector<float>   *MuonInfo_MuonType;
   vector<float>   *MuonInfo_Pt;
   vector<float>   *MuonInfo_Eta;
   vector<float>   *MuonInfo_Phi;
   vector<float>   *MuonInfo_Energy;
   vector<float>   *MuonInfo_BestTrackDz;
   vector<float>   *MuonInfo_BestTrackDxy;
   vector<float>   *MuonInfo_PFIsoDeltaBetaCorrR04;
   vector<float>   *MuonInfo_TrackerBasedIsoR03;
   vector<bool>    *MuonInfo_CutBasedIdMedium;
   vector<bool>    *MuonInfo_CutBasedIdTight;
   vector<bool>    *MuonInfo_GenMatch;
   vector<int>     *MuonInfo_GenPdgID;
   vector<float>   *MuonInfo_GenPt;
   vector<float>   *MuonInfo_GenEta;
   vector<float>   *MuonInfo_GenPhi;
   Int_t           jets_size;
   vector<float>   *JetInfo_Pt;
   vector<float>   *JetInfo_Eta;
   vector<float>   *JetInfo_Phi;
   vector<float>   *JetInfo_Mass;
   vector<float>   *JetInfo_Energy;
   vector<float>   *JetInfo_PtRaw;
   vector<float>   *JetInfo_QGL;
   vector<float>   *JetInfo_RMS;
   vector<float>   *JetInfo_puJetIdMVA;
   vector<bool>    *JetInfo_GenJetMatch;
   vector<float>   *JetInfo_pfCombinedInclusiveSecondaryVertexV2BJetTags;
   vector<float>   *JetInfo_pfCombinedMVAV2BJetTags;
   vector<float>   *JetInfo_pfDeepCSVJetTags_probb;
   vector<float>   *JetInfo_pfDeepCSVJetTags_probbb;
   vector<float>   *JetInfo_pfDeepCSVJetTags_probc;
   vector<float>   *JetInfo_pfDeepCSVJetTags_probudsg;
   vector<float>   *JetInfo_JECScale;
   vector<float>   *JetInfo_JERScale;
   vector<float>   *JetInfo_JECUnc;
   vector<float>   *JetInfo_JERUp;
   vector<float>   *JetInfo_JERDown;
   vector<bool>    *JetInfo_GenPartonMatch;
   vector<float>   *JetInfo_GenPt;
   vector<float>   *JetInfo_GenEta;
   vector<float>   *JetInfo_GenPhi;
   vector<int>     *JetInfo_GenPdgID;
   vector<int>     *JetInfo_GenFlavor;
   vector<int>     *JetInfo_GenHadronFlavor;
   Float_t         MetInfo_Pt;
   Float_t         MetInfo_Phi;
   Float_t         MetInfo_Px;
   Float_t         MetInfo_Py;
   Float_t         MetInfo_SumET;
   Float_t         MetInfo_CorrPtShiftJetEnUp;
   Float_t         MetInfo_CorrPtShiftJetEnDown;
   Float_t         MetInfo_CorrPtShiftJetResUp;
   Float_t         MetInfo_CorrPtShiftJetResDown;
   Float_t         MetInfo_CorrPtShiftUncEnUp;
   Float_t         MetInfo_CorrPtShiftUncEnDown;
   Float_t         MetInfo_CorrPtShiftPhoEnUp;
   Float_t         MetInfo_CorrPtShiftPhoEnDown;
   Float_t         MetInfo_CorrPhiShiftJetEnUp;
   Float_t         MetInfo_CorrPhiShiftJetEnDown;
   Float_t         MetInfo_CorrPhiShiftJetResUp;
   Float_t         MetInfo_CorrPhiShiftJetResDown;
   Float_t         MetInfo_CorrPhiShiftUncEnUp;
   Float_t         MetInfo_CorrPhiShiftUncEnDown;
   Float_t         MetInfo_CorrPhiShiftPhoEnUp;
   Float_t         MetInfo_CorrPhiShiftPhoEnDown;
   Int_t           GenPartInfo_size;
   vector<float>   *GenPartInfo_Pt;
   vector<float>   *GenPartInfo_Eta;
   vector<float>   *GenPartInfo_Phi;
   vector<float>   *GenPartInfo_Mass;
   vector<int>     *GenPartInfo_PdgID;
   vector<int>     *GenPartInfo_Status;
   vector<int>     *GenPartInfo_nMo;
   vector<int>     *GenPartInfo_nDa;
   Int_t           HTXSstage0cat;
   Int_t           HTXSstage1cat;
   Int_t           HTXSnjets;
   Float_t         HTXSpTH;
   Float_t         HTXSpTV;

   // List of branches
   TBranch        *b_EvtInfo_NPu;   //!
   TBranch        *b_EvtInfo_NVtx;   //!
   TBranch        *b_EvtInfo_passTrigger;   //!
   TBranch        *b_EvtInfo_genweight;   //!
   TBranch        *b_EvtInfo_Rho;   //!
   TBranch        *b_EvtInfo_PVz;   //!
   TBranch        *b_EvtInfo_BSsigmaz;   //!
   TBranch        *b_EvtInfo_Flag_HBHENoiseFilter;   //!
   TBranch        *b_EvtInfo_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_EvtInfo_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_EvtInfo_Flag_goodVertices;   //!
   TBranch        *b_EvtInfo_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_EvtInfo_Flag_BadPFMuonFilter;   //!
   TBranch        *b_EvtInfo_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_EvtInfo_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_EvtInfo_Flag_eeBadScFilter;   //!
   TBranch        *b_DiPhoInfo_mass;   //!
   TBranch        *b_DiPhoInfo_pt;   //!
   TBranch        *b_DiPhoInfo_leadPt;   //!
   TBranch        *b_DiPhoInfo_leadEta;   //!
   TBranch        *b_DiPhoInfo_leadPhi;   //!
   TBranch        *b_DiPhoInfo_leadE;   //!
   TBranch        *b_DiPhoInfo_leadsigEOverE;   //!
   TBranch        *b_DiPhoInfo_leadR9;   //!
   TBranch        *b_DiPhoInfo_leadsieie;   //!
   TBranch        *b_DiPhoInfo_leadhoe;   //!
   TBranch        *b_DiPhoInfo_leadIDMVA;   //!
   TBranch        *b_DiPhoInfo_leadhasPixelSeed;   //!
   TBranch        *b_DiPhoInfo_leadGenMatch;   //!
   TBranch        *b_DiPhoInfo_leadGenMatchType;   //!
   TBranch        *b_DiPhoInfo_subleadPt;   //!
   TBranch        *b_DiPhoInfo_subleadEta;   //!
   TBranch        *b_DiPhoInfo_subleadPhi;   //!
   TBranch        *b_DiPhoInfo_subleadE;   //!
   TBranch        *b_DiPhoInfo_subleadsigEOverE;   //!
   TBranch        *b_DiPhoInfo_subleadR9;   //!
   TBranch        *b_DiPhoInfo_subleadsieie;   //!
   TBranch        *b_DiPhoInfo_subleadhoe;   //!
   TBranch        *b_DiPhoInfo_subleadIDMVA;   //!
   TBranch        *b_DiPhoInfo_subleadhasPixelSeed;   //!
   TBranch        *b_DiPhoInfo_subleadGenMatch;   //!
   TBranch        *b_DiPhoInfo_subleadGenMatchType;   //!
   TBranch        *b_DiPhoInfo_SelectedVz;   //!
   TBranch        *b_DiPhoInfo_GenVz;   //!
   TBranch        *b_DiPhoInfo_centralWeight;   //!
   TBranch        *b_DiPhoInfo_MvaLinearSystUp;   //!
   TBranch        *b_DiPhoInfo_MvaLinearSystDown;   //!
   TBranch        *b_DiPhoInfo_LooseMvaSFUp;   //!
   TBranch        *b_DiPhoInfo_LooseMvaSFDown;   //!
   TBranch        *b_DiPhoInfo_PreselSFUp;   //!
   TBranch        *b_DiPhoInfo_PreselSFDown;   //!
   TBranch        *b_DiPhoInfo_electronVetoSFUp;   //!
   TBranch        *b_DiPhoInfo_electronVetoSFDown;   //!
   TBranch        *b_DiPhoInfo_TriggerWeightUp;   //!
   TBranch        *b_DiPhoInfo_TriggerWeightDown;   //!
   TBranch        *b_DiPhoInfo_FracRVWeightUp;   //!
   TBranch        *b_DiPhoInfo_FracRVWeightDown;   //!
   TBranch        *b_DiPhoInfo_FracRVNvtxWeightUp;   //!
   TBranch        *b_DiPhoInfo_FracRVNvtxWeightDown;   //!
   TBranch        *b_ElecInfo_Size;   //!
   TBranch        *b_ElecInfo_Charge;   //!
   TBranch        *b_ElecInfo_Pt;   //!
   TBranch        *b_ElecInfo_Eta;   //!
   TBranch        *b_ElecInfo_Phi;   //!
   TBranch        *b_ElecInfo_Energy;   //!
   TBranch        *b_ElecInfo_EtaSC;   //!
   TBranch        *b_ElecInfo_PhiSC;   //!
   TBranch        *b_ElecInfo_GsfTrackDz;   //!
   TBranch        *b_ElecInfo_GsfTrackDxy;   //!
   TBranch        *b_ElecInfo_EGMCutBasedIDVeto;   //!
   TBranch        *b_ElecInfo_EGMCutBasedIDLoose;   //!
   TBranch        *b_ElecInfo_EGMCutBasedIDMedium;   //!
   TBranch        *b_ElecInfo_EGMCutBasedIDTight;   //!
   TBranch        *b_ElecInfo_fggPhoVeto;   //!
   TBranch        *b_ElecInfo_tmpPhoVeto;   //!
   TBranch        *b_ElecInfo_EnergyCorrFactor;   //!
   TBranch        *b_ElecInfo_EnergyPostCorrErr;   //!
   TBranch        *b_ElecInfo_EnergyPostCorrScaleUp;   //!
   TBranch        *b_ElecInfo_EnergyPostCorrScaleDown;   //!
   TBranch        *b_ElecInfo_EnergyPostCorrSmearUp;   //!
   TBranch        *b_ElecInfo_EnergyPostCorrSmearDown;   //!
   TBranch        *b_ElecInfo_GenMatch;   //!
   TBranch        *b_ElecInfo_GenPdgID;   //!
   TBranch        *b_ElecInfo_GenPt;   //!
   TBranch        *b_ElecInfo_GenEta;   //!
   TBranch        *b_ElecInfo_GenPhi;   //!
   TBranch        *b_MuonInfo_Size;   //!
   TBranch        *b_MuonInfo_Charge;   //!
   TBranch        *b_MuonInfo_MuonType;   //!
   TBranch        *b_MuonInfo_Pt;   //!
   TBranch        *b_MuonInfo_Eta;   //!
   TBranch        *b_MuonInfo_Phi;   //!
   TBranch        *b_MuonInfo_Energy;   //!
   TBranch        *b_MuonInfo_BestTrackDz;   //!
   TBranch        *b_MuonInfo_BestTrackDxy;   //!
   TBranch        *b_MuonInfo_PFIsoDeltaBetaCorrR04;   //!
   TBranch        *b_MuonInfo_TrackerBasedIsoR03;   //!
   TBranch        *b_MuonInfo_CutBasedIdMedium;   //!
   TBranch        *b_MuonInfo_CutBasedIdTight;   //!
   TBranch        *b_MuonInfo_GenMatch;   //!
   TBranch        *b_MuonInfo_GenPdgID;   //!
   TBranch        *b_MuonInfo_GenPt;   //!
   TBranch        *b_MuonInfo_GenEta;   //!
   TBranch        *b_MuonInfo_GenPhi;   //!
   TBranch        *b_jets_size;   //!
   TBranch        *b_JetInfo_Pt;   //!
   TBranch        *b_JetInfo_Eta;   //!
   TBranch        *b_JetInfo_Phi;   //!
   TBranch        *b_JetInfo_Mass;   //!
   TBranch        *b_JetInfo_Energy;   //!
   TBranch        *b_JetInfo_PtRaw;   //!
   TBranch        *b_JetInfo_QGL;   //!
   TBranch        *b_JetInfo_RMS;   //!
   TBranch        *b_JetInfo_puJetIdMVA;   //!
   TBranch        *b_JetInfo_GenJetMatch;   //!
   TBranch        *b_JetInfo_pfCombinedInclusiveSecondaryVertexV2BJetTags;   //!
   TBranch        *b_JetInfo_pfCombinedMVAV2BJetTags;   //!
   TBranch        *b_JetInfo_pfDeepCSVJetTags_probb;   //!
   TBranch        *b_JetInfo_pfDeepCSVJetTags_probbb;   //!
   TBranch        *b_JetInfo_pfDeepCSVJetTags_probc;   //!
   TBranch        *b_JetInfo_pfDeepCSVJetTags_probudsg;   //!
   TBranch        *b_JetInfo_JECScale;   //!
   TBranch        *b_JetInfo_JERScale;   //!
   TBranch        *b_JetInfo_JECUnc;   //!
   TBranch        *b_JetInfo_JERUp;   //!
   TBranch        *b_JetInfo_JERDown;   //!
   TBranch        *b_JetInfo_GenPartonMatch;   //!
   TBranch        *b_JetInfo_GenPt;   //!
   TBranch        *b_JetInfo_GenEta;   //!
   TBranch        *b_JetInfo_GenPhi;   //!
   TBranch        *b_JetInfo_GenPdgID;   //!
   TBranch        *b_JetInfo_GenFlavor;   //!
   TBranch        *b_JetInfo_GenHadronFlavor;   //!
   TBranch        *b_MetInfo_Pt;   //!
   TBranch        *b_MetInfo_Phi;   //!
   TBranch        *b_MetInfo_Px;   //!
   TBranch        *b_MetInfo_Py;   //!
   TBranch        *b_MetInfo_SumET;   //!
   TBranch        *b_MetInfo_CorrPtShiftJetEnUp;   //!
   TBranch        *b_MetInfo_CorrPtShiftJetEnDown;   //!
   TBranch        *b_MetInfo_CorrPtShiftJetResUp;   //!
   TBranch        *b_MetInfo_CorrPtShiftJetResDown;   //!
   TBranch        *b_MetInfo_CorrPtShiftUncEnUp;   //!
   TBranch        *b_MetInfo_CorrPtShiftUncEnDown;   //!
   TBranch        *b_MetInfo_CorrPtShiftPhoEnUp;   //!
   TBranch        *b_MetInfo_CorrPtShiftPhoEnDown;   //!
   TBranch        *b_MetInfo_CorrPhiShiftJetEnUp;   //!
   TBranch        *b_MetInfo_CorrPhiShiftJetEnDown;   //!
   TBranch        *b_MetInfo_CorrPhiShiftJetResUp;   //!
   TBranch        *b_MetInfo_CorrPhiShiftJetResDown;   //!
   TBranch        *b_MetInfo_CorrPhiShiftUncEnUp;   //!
   TBranch        *b_MetInfo_CorrPhiShiftUncEnDown;   //!
   TBranch        *b_MetInfo_CorrPhiShiftPhoEnUp;   //!
   TBranch        *b_MetInfo_CorrPhiShiftPhoEnDown;   //!
   TBranch        *b_GenPartInfo_size;   //!
   TBranch        *b_GenPartInfo_Pt;   //!
   TBranch        *b_GenPartInfo_Eta;   //!
   TBranch        *b_GenPartInfo_Phi;   //!
   TBranch        *b_GenPartInfo_Mass;   //!
   TBranch        *b_GenPartInfo_PdgID;   //!
   TBranch        *b_GenPartInfo_Status;   //!
   TBranch        *b_GenPartInfo_nMo;   //!
   TBranch        *b_GenPartInfo_nDa;   //!
   TBranch        *b_HTXSstage0cat;   //!
   TBranch        *b_HTXSstage1cat;   //!
   TBranch        *b_HTXSnjets;   //!
   TBranch        *b_HTXSpTH;   //!
   TBranch        *b_HTXSpTV;   //!

   flashggStdTree(TTree *tree=0);
   virtual ~flashggStdTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef flashggStdTree_cxx
flashggStdTree::flashggStdTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root:/flashggNtuples");
      dir->GetObject("flashggStdTree",tree);

   }
   Init(tree);
}

flashggStdTree::~flashggStdTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t flashggStdTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t flashggStdTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void flashggStdTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ElecInfo_Charge = 0;
   ElecInfo_Pt = 0;
   ElecInfo_Eta = 0;
   ElecInfo_Phi = 0;
   ElecInfo_Energy = 0;
   ElecInfo_EtaSC = 0;
   ElecInfo_PhiSC = 0;
   ElecInfo_GsfTrackDz = 0;
   ElecInfo_GsfTrackDxy = 0;
   ElecInfo_EGMCutBasedIDVeto = 0;
   ElecInfo_EGMCutBasedIDLoose = 0;
   ElecInfo_EGMCutBasedIDMedium = 0;
   ElecInfo_EGMCutBasedIDTight = 0;
   ElecInfo_fggPhoVeto = 0;
   ElecInfo_tmpPhoVeto = 0;
   ElecInfo_EnergyCorrFactor = 0;
   ElecInfo_EnergyPostCorrErr = 0;
   ElecInfo_EnergyPostCorrScaleUp = 0;
   ElecInfo_EnergyPostCorrScaleDown = 0;
   ElecInfo_EnergyPostCorrSmearUp = 0;
   ElecInfo_EnergyPostCorrSmearDown = 0;
   ElecInfo_GenMatch = 0;
   ElecInfo_GenPdgID = 0;
   ElecInfo_GenPt = 0;
   ElecInfo_GenEta = 0;
   ElecInfo_GenPhi = 0;
   MuonInfo_Charge = 0;
   MuonInfo_MuonType = 0;
   MuonInfo_Pt = 0;
   MuonInfo_Eta = 0;
   MuonInfo_Phi = 0;
   MuonInfo_Energy = 0;
   MuonInfo_BestTrackDz = 0;
   MuonInfo_BestTrackDxy = 0;
   MuonInfo_PFIsoDeltaBetaCorrR04 = 0;
   MuonInfo_TrackerBasedIsoR03 = 0;
   MuonInfo_CutBasedIdMedium = 0;
   MuonInfo_CutBasedIdTight = 0;
   MuonInfo_GenMatch = 0;
   MuonInfo_GenPdgID = 0;
   MuonInfo_GenPt = 0;
   MuonInfo_GenEta = 0;
   MuonInfo_GenPhi = 0;
   JetInfo_Pt = 0;
   JetInfo_Eta = 0;
   JetInfo_Phi = 0;
   JetInfo_Mass = 0;
   JetInfo_Energy = 0;
   JetInfo_PtRaw = 0;
   JetInfo_QGL = 0;
   JetInfo_RMS = 0;
   JetInfo_puJetIdMVA = 0;
   JetInfo_GenJetMatch = 0;
   JetInfo_pfCombinedInclusiveSecondaryVertexV2BJetTags = 0;
   JetInfo_pfCombinedMVAV2BJetTags = 0;
   JetInfo_pfDeepCSVJetTags_probb = 0;
   JetInfo_pfDeepCSVJetTags_probbb = 0;
   JetInfo_pfDeepCSVJetTags_probc = 0;
   JetInfo_pfDeepCSVJetTags_probudsg = 0;
   JetInfo_JECScale = 0;
   JetInfo_JERScale = 0;
   JetInfo_JECUnc = 0;
   JetInfo_JERUp = 0;
   JetInfo_JERDown = 0;
   JetInfo_GenPartonMatch = 0;
   JetInfo_GenPt = 0;
   JetInfo_GenEta = 0;
   JetInfo_GenPhi = 0;
   JetInfo_GenPdgID = 0;
   JetInfo_GenFlavor = 0;
   JetInfo_GenHadronFlavor = 0;
   GenPartInfo_Pt = 0;
   GenPartInfo_Eta = 0;
   GenPartInfo_Phi = 0;
   GenPartInfo_Mass = 0;
   GenPartInfo_PdgID = 0;
   GenPartInfo_Status = 0;
   GenPartInfo_nMo = 0;
   GenPartInfo_nDa = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EvtInfo.NPu", &EvtInfo_NPu, &b_EvtInfo_NPu);
   fChain->SetBranchAddress("EvtInfo.NVtx", &EvtInfo_NVtx, &b_EvtInfo_NVtx);
   fChain->SetBranchAddress("EvtInfo.passTrigger", &EvtInfo_passTrigger, &b_EvtInfo_passTrigger);
   fChain->SetBranchAddress("EvtInfo.genweight", &EvtInfo_genweight, &b_EvtInfo_genweight);
   fChain->SetBranchAddress("EvtInfo.Rho", &EvtInfo_Rho, &b_EvtInfo_Rho);
   fChain->SetBranchAddress("EvtInfo.PVz", &EvtInfo_PVz, &b_EvtInfo_PVz);
   fChain->SetBranchAddress("EvtInfo.BSsigmaz", &EvtInfo_BSsigmaz, &b_EvtInfo_BSsigmaz);
   fChain->SetBranchAddress("EvtInfo.Flag_HBHENoiseFilter", &EvtInfo_Flag_HBHENoiseFilter, &b_EvtInfo_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("EvtInfo.Flag_HBHENoiseIsoFilter", &EvtInfo_Flag_HBHENoiseIsoFilter, &b_EvtInfo_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("EvtInfo.Flag_EcalDeadCellTriggerPrimitiveFilter", &EvtInfo_Flag_EcalDeadCellTriggerPrimitiveFilter, &b_EvtInfo_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("EvtInfo.Flag_goodVertices", &EvtInfo_Flag_goodVertices, &b_EvtInfo_Flag_goodVertices);
   fChain->SetBranchAddress("EvtInfo.Flag_globalSuperTightHalo2016Filter", &EvtInfo_Flag_globalSuperTightHalo2016Filter, &b_EvtInfo_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("EvtInfo.Flag_BadPFMuonFilter", &EvtInfo_Flag_BadPFMuonFilter, &b_EvtInfo_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("EvtInfo.Flag_BadChargedCandidateFilter", &EvtInfo_Flag_BadChargedCandidateFilter, &b_EvtInfo_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("EvtInfo.Flag_ecalBadCalibFilter", &EvtInfo_Flag_ecalBadCalibFilter, &b_EvtInfo_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("EvtInfo.Flag_eeBadScFilter", &EvtInfo_Flag_eeBadScFilter, &b_EvtInfo_Flag_eeBadScFilter);
   fChain->SetBranchAddress("DiPhoInfo.mass", &DiPhoInfo_mass, &b_DiPhoInfo_mass);
   fChain->SetBranchAddress("DiPhoInfo.pt", &DiPhoInfo_pt, &b_DiPhoInfo_pt);
   fChain->SetBranchAddress("DiPhoInfo.leadPt", &DiPhoInfo_leadPt, &b_DiPhoInfo_leadPt);
   fChain->SetBranchAddress("DiPhoInfo.leadEta", &DiPhoInfo_leadEta, &b_DiPhoInfo_leadEta);
   fChain->SetBranchAddress("DiPhoInfo.leadPhi", &DiPhoInfo_leadPhi, &b_DiPhoInfo_leadPhi);
   fChain->SetBranchAddress("DiPhoInfo.leadE", &DiPhoInfo_leadE, &b_DiPhoInfo_leadE);
   fChain->SetBranchAddress("DiPhoInfo.leadsigEOverE", &DiPhoInfo_leadsigEOverE, &b_DiPhoInfo_leadsigEOverE);
   fChain->SetBranchAddress("DiPhoInfo.leadR9", &DiPhoInfo_leadR9, &b_DiPhoInfo_leadR9);
   fChain->SetBranchAddress("DiPhoInfo.leadsieie", &DiPhoInfo_leadsieie, &b_DiPhoInfo_leadsieie);
   fChain->SetBranchAddress("DiPhoInfo.leadhoe", &DiPhoInfo_leadhoe, &b_DiPhoInfo_leadhoe);
   fChain->SetBranchAddress("DiPhoInfo.leadIDMVA", &DiPhoInfo_leadIDMVA, &b_DiPhoInfo_leadIDMVA);
   fChain->SetBranchAddress("DiPhoInfo.leadhasPixelSeed", &DiPhoInfo_leadhasPixelSeed, &b_DiPhoInfo_leadhasPixelSeed);
   fChain->SetBranchAddress("DiPhoInfo.leadGenMatch", &DiPhoInfo_leadGenMatch, &b_DiPhoInfo_leadGenMatch);
   fChain->SetBranchAddress("DiPhoInfo.leadGenMatchType", &DiPhoInfo_leadGenMatchType, &b_DiPhoInfo_leadGenMatchType);
   fChain->SetBranchAddress("DiPhoInfo.subleadPt", &DiPhoInfo_subleadPt, &b_DiPhoInfo_subleadPt);
   fChain->SetBranchAddress("DiPhoInfo.subleadEta", &DiPhoInfo_subleadEta, &b_DiPhoInfo_subleadEta);
   fChain->SetBranchAddress("DiPhoInfo.subleadPhi", &DiPhoInfo_subleadPhi, &b_DiPhoInfo_subleadPhi);
   fChain->SetBranchAddress("DiPhoInfo.subleadE", &DiPhoInfo_subleadE, &b_DiPhoInfo_subleadE);
   fChain->SetBranchAddress("DiPhoInfo.subleadsigEOverE", &DiPhoInfo_subleadsigEOverE, &b_DiPhoInfo_subleadsigEOverE);
   fChain->SetBranchAddress("DiPhoInfo.subleadR9", &DiPhoInfo_subleadR9, &b_DiPhoInfo_subleadR9);
   fChain->SetBranchAddress("DiPhoInfo.subleadsieie", &DiPhoInfo_subleadsieie, &b_DiPhoInfo_subleadsieie);
   fChain->SetBranchAddress("DiPhoInfo.subleadhoe", &DiPhoInfo_subleadhoe, &b_DiPhoInfo_subleadhoe);
   fChain->SetBranchAddress("DiPhoInfo.subleadIDMVA", &DiPhoInfo_subleadIDMVA, &b_DiPhoInfo_subleadIDMVA);
   fChain->SetBranchAddress("DiPhoInfo.subleadhasPixelSeed", &DiPhoInfo_subleadhasPixelSeed, &b_DiPhoInfo_subleadhasPixelSeed);
   fChain->SetBranchAddress("DiPhoInfo.subleadGenMatch", &DiPhoInfo_subleadGenMatch, &b_DiPhoInfo_subleadGenMatch);
   fChain->SetBranchAddress("DiPhoInfo.subleadGenMatchType", &DiPhoInfo_subleadGenMatchType, &b_DiPhoInfo_subleadGenMatchType);
   fChain->SetBranchAddress("DiPhoInfo.SelectedVz", &DiPhoInfo_SelectedVz, &b_DiPhoInfo_SelectedVz);
   fChain->SetBranchAddress("DiPhoInfo.GenVz", &DiPhoInfo_GenVz, &b_DiPhoInfo_GenVz);
   fChain->SetBranchAddress("DiPhoInfo.centralWeight", &DiPhoInfo_centralWeight, &b_DiPhoInfo_centralWeight);
   fChain->SetBranchAddress("DiPhoInfo.MvaLinearSystUp", &DiPhoInfo_MvaLinearSystUp, &b_DiPhoInfo_MvaLinearSystUp);
   fChain->SetBranchAddress("DiPhoInfo.MvaLinearSystDown", &DiPhoInfo_MvaLinearSystDown, &b_DiPhoInfo_MvaLinearSystDown);
   fChain->SetBranchAddress("DiPhoInfo.LooseMvaSFUp", &DiPhoInfo_LooseMvaSFUp, &b_DiPhoInfo_LooseMvaSFUp);
   fChain->SetBranchAddress("DiPhoInfo.LooseMvaSFDown", &DiPhoInfo_LooseMvaSFDown, &b_DiPhoInfo_LooseMvaSFDown);
   fChain->SetBranchAddress("DiPhoInfo.PreselSFUp", &DiPhoInfo_PreselSFUp, &b_DiPhoInfo_PreselSFUp);
   fChain->SetBranchAddress("DiPhoInfo.PreselSFDown", &DiPhoInfo_PreselSFDown, &b_DiPhoInfo_PreselSFDown);
   fChain->SetBranchAddress("DiPhoInfo.electronVetoSFUp", &DiPhoInfo_electronVetoSFUp, &b_DiPhoInfo_electronVetoSFUp);
   fChain->SetBranchAddress("DiPhoInfo.electronVetoSFDown", &DiPhoInfo_electronVetoSFDown, &b_DiPhoInfo_electronVetoSFDown);
   fChain->SetBranchAddress("DiPhoInfo.TriggerWeightUp", &DiPhoInfo_TriggerWeightUp, &b_DiPhoInfo_TriggerWeightUp);
   fChain->SetBranchAddress("DiPhoInfo.TriggerWeightDown", &DiPhoInfo_TriggerWeightDown, &b_DiPhoInfo_TriggerWeightDown);
   fChain->SetBranchAddress("DiPhoInfo.FracRVWeightUp", &DiPhoInfo_FracRVWeightUp, &b_DiPhoInfo_FracRVWeightUp);
   fChain->SetBranchAddress("DiPhoInfo.FracRVWeightDown", &DiPhoInfo_FracRVWeightDown, &b_DiPhoInfo_FracRVWeightDown);
   fChain->SetBranchAddress("DiPhoInfo.FracRVNvtxWeightUp", &DiPhoInfo_FracRVNvtxWeightUp, &b_DiPhoInfo_FracRVNvtxWeightUp);
   fChain->SetBranchAddress("DiPhoInfo.FracRVNvtxWeightDown", &DiPhoInfo_FracRVNvtxWeightDown, &b_DiPhoInfo_FracRVNvtxWeightDown);
   fChain->SetBranchAddress("ElecInfo.Size", &ElecInfo_Size, &b_ElecInfo_Size);
   fChain->SetBranchAddress("ElecInfo.Charge", &ElecInfo_Charge, &b_ElecInfo_Charge);
   fChain->SetBranchAddress("ElecInfo.Pt", &ElecInfo_Pt, &b_ElecInfo_Pt);
   fChain->SetBranchAddress("ElecInfo.Eta", &ElecInfo_Eta, &b_ElecInfo_Eta);
   fChain->SetBranchAddress("ElecInfo.Phi", &ElecInfo_Phi, &b_ElecInfo_Phi);
   fChain->SetBranchAddress("ElecInfo.Energy", &ElecInfo_Energy, &b_ElecInfo_Energy);
   fChain->SetBranchAddress("ElecInfo.EtaSC", &ElecInfo_EtaSC, &b_ElecInfo_EtaSC);
   fChain->SetBranchAddress("ElecInfo.PhiSC", &ElecInfo_PhiSC, &b_ElecInfo_PhiSC);
   fChain->SetBranchAddress("ElecInfo.GsfTrackDz", &ElecInfo_GsfTrackDz, &b_ElecInfo_GsfTrackDz);
   fChain->SetBranchAddress("ElecInfo.GsfTrackDxy", &ElecInfo_GsfTrackDxy, &b_ElecInfo_GsfTrackDxy);
   fChain->SetBranchAddress("ElecInfo.EGMCutBasedIDVeto", &ElecInfo_EGMCutBasedIDVeto, &b_ElecInfo_EGMCutBasedIDVeto);
   fChain->SetBranchAddress("ElecInfo.EGMCutBasedIDLoose", &ElecInfo_EGMCutBasedIDLoose, &b_ElecInfo_EGMCutBasedIDLoose);
   fChain->SetBranchAddress("ElecInfo.EGMCutBasedIDMedium", &ElecInfo_EGMCutBasedIDMedium, &b_ElecInfo_EGMCutBasedIDMedium);
   fChain->SetBranchAddress("ElecInfo.EGMCutBasedIDTight", &ElecInfo_EGMCutBasedIDTight, &b_ElecInfo_EGMCutBasedIDTight);
   fChain->SetBranchAddress("ElecInfo.fggPhoVeto", &ElecInfo_fggPhoVeto, &b_ElecInfo_fggPhoVeto);
   fChain->SetBranchAddress("ElecInfo.tmpPhoVeto", &ElecInfo_tmpPhoVeto, &b_ElecInfo_tmpPhoVeto);
   fChain->SetBranchAddress("ElecInfo.EnergyCorrFactor", &ElecInfo_EnergyCorrFactor, &b_ElecInfo_EnergyCorrFactor);
   fChain->SetBranchAddress("ElecInfo.EnergyPostCorrErr", &ElecInfo_EnergyPostCorrErr, &b_ElecInfo_EnergyPostCorrErr);
   fChain->SetBranchAddress("ElecInfo.EnergyPostCorrScaleUp", &ElecInfo_EnergyPostCorrScaleUp, &b_ElecInfo_EnergyPostCorrScaleUp);
   fChain->SetBranchAddress("ElecInfo.EnergyPostCorrScaleDown", &ElecInfo_EnergyPostCorrScaleDown, &b_ElecInfo_EnergyPostCorrScaleDown);
   fChain->SetBranchAddress("ElecInfo.EnergyPostCorrSmearUp", &ElecInfo_EnergyPostCorrSmearUp, &b_ElecInfo_EnergyPostCorrSmearUp);
   fChain->SetBranchAddress("ElecInfo.EnergyPostCorrSmearDown", &ElecInfo_EnergyPostCorrSmearDown, &b_ElecInfo_EnergyPostCorrSmearDown);
   fChain->SetBranchAddress("ElecInfo.GenMatch", &ElecInfo_GenMatch, &b_ElecInfo_GenMatch);
   fChain->SetBranchAddress("ElecInfo.GenPdgID", &ElecInfo_GenPdgID, &b_ElecInfo_GenPdgID);
   fChain->SetBranchAddress("ElecInfo.GenPt", &ElecInfo_GenPt, &b_ElecInfo_GenPt);
   fChain->SetBranchAddress("ElecInfo.GenEta", &ElecInfo_GenEta, &b_ElecInfo_GenEta);
   fChain->SetBranchAddress("ElecInfo.GenPhi", &ElecInfo_GenPhi, &b_ElecInfo_GenPhi);
   fChain->SetBranchAddress("MuonInfo.Size", &MuonInfo_Size, &b_MuonInfo_Size);
   fChain->SetBranchAddress("MuonInfo.Charge", &MuonInfo_Charge, &b_MuonInfo_Charge);
   fChain->SetBranchAddress("MuonInfo.MuonType", &MuonInfo_MuonType, &b_MuonInfo_MuonType);
   fChain->SetBranchAddress("MuonInfo.Pt", &MuonInfo_Pt, &b_MuonInfo_Pt);
   fChain->SetBranchAddress("MuonInfo.Eta", &MuonInfo_Eta, &b_MuonInfo_Eta);
   fChain->SetBranchAddress("MuonInfo.Phi", &MuonInfo_Phi, &b_MuonInfo_Phi);
   fChain->SetBranchAddress("MuonInfo.Energy", &MuonInfo_Energy, &b_MuonInfo_Energy);
   fChain->SetBranchAddress("MuonInfo.BestTrackDz", &MuonInfo_BestTrackDz, &b_MuonInfo_BestTrackDz);
   fChain->SetBranchAddress("MuonInfo.BestTrackDxy", &MuonInfo_BestTrackDxy, &b_MuonInfo_BestTrackDxy);
   fChain->SetBranchAddress("MuonInfo.PFIsoDeltaBetaCorrR04", &MuonInfo_PFIsoDeltaBetaCorrR04, &b_MuonInfo_PFIsoDeltaBetaCorrR04);
   fChain->SetBranchAddress("MuonInfo.TrackerBasedIsoR03", &MuonInfo_TrackerBasedIsoR03, &b_MuonInfo_TrackerBasedIsoR03);
   fChain->SetBranchAddress("MuonInfo.CutBasedIdMedium", &MuonInfo_CutBasedIdMedium, &b_MuonInfo_CutBasedIdMedium);
   fChain->SetBranchAddress("MuonInfo.CutBasedIdTight", &MuonInfo_CutBasedIdTight, &b_MuonInfo_CutBasedIdTight);
   fChain->SetBranchAddress("MuonInfo.GenMatch", &MuonInfo_GenMatch, &b_MuonInfo_GenMatch);
   fChain->SetBranchAddress("MuonInfo.GenPdgID", &MuonInfo_GenPdgID, &b_MuonInfo_GenPdgID);
   fChain->SetBranchAddress("MuonInfo.GenPt", &MuonInfo_GenPt, &b_MuonInfo_GenPt);
   fChain->SetBranchAddress("MuonInfo.GenEta", &MuonInfo_GenEta, &b_MuonInfo_GenEta);
   fChain->SetBranchAddress("MuonInfo.GenPhi", &MuonInfo_GenPhi, &b_MuonInfo_GenPhi);
   fChain->SetBranchAddress("jets_size", &jets_size, &b_jets_size);
   fChain->SetBranchAddress("JetInfo.Pt", &JetInfo_Pt, &b_JetInfo_Pt);
   fChain->SetBranchAddress("JetInfo.Eta", &JetInfo_Eta, &b_JetInfo_Eta);
   fChain->SetBranchAddress("JetInfo.Phi", &JetInfo_Phi, &b_JetInfo_Phi);
   fChain->SetBranchAddress("JetInfo.Mass", &JetInfo_Mass, &b_JetInfo_Mass);
   fChain->SetBranchAddress("JetInfo.Energy", &JetInfo_Energy, &b_JetInfo_Energy);
   fChain->SetBranchAddress("JetInfo.PtRaw", &JetInfo_PtRaw, &b_JetInfo_PtRaw);
   fChain->SetBranchAddress("JetInfo.QGL", &JetInfo_QGL, &b_JetInfo_QGL);
   fChain->SetBranchAddress("JetInfo.RMS", &JetInfo_RMS, &b_JetInfo_RMS);
   fChain->SetBranchAddress("JetInfo.puJetIdMVA", &JetInfo_puJetIdMVA, &b_JetInfo_puJetIdMVA);
   fChain->SetBranchAddress("JetInfo.GenJetMatch", &JetInfo_GenJetMatch, &b_JetInfo_GenJetMatch);
   fChain->SetBranchAddress("JetInfo.pfCombinedInclusiveSecondaryVertexV2BJetTags", &JetInfo_pfCombinedInclusiveSecondaryVertexV2BJetTags, &b_JetInfo_pfCombinedInclusiveSecondaryVertexV2BJetTags);
   fChain->SetBranchAddress("JetInfo.pfCombinedMVAV2BJetTags", &JetInfo_pfCombinedMVAV2BJetTags, &b_JetInfo_pfCombinedMVAV2BJetTags);
   fChain->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probb", &JetInfo_pfDeepCSVJetTags_probb, &b_JetInfo_pfDeepCSVJetTags_probb);
   fChain->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probbb", &JetInfo_pfDeepCSVJetTags_probbb, &b_JetInfo_pfDeepCSVJetTags_probbb);
   fChain->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probc", &JetInfo_pfDeepCSVJetTags_probc, &b_JetInfo_pfDeepCSVJetTags_probc);
   fChain->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probudsg", &JetInfo_pfDeepCSVJetTags_probudsg, &b_JetInfo_pfDeepCSVJetTags_probudsg);
   fChain->SetBranchAddress("JetInfo.JECScale", &JetInfo_JECScale, &b_JetInfo_JECScale);
   fChain->SetBranchAddress("JetInfo.JERScale", &JetInfo_JERScale, &b_JetInfo_JERScale);
   fChain->SetBranchAddress("JetInfo.JECUnc", &JetInfo_JECUnc, &b_JetInfo_JECUnc);
   fChain->SetBranchAddress("JetInfo.JERUp", &JetInfo_JERUp, &b_JetInfo_JERUp);
   fChain->SetBranchAddress("JetInfo.JERDown", &JetInfo_JERDown, &b_JetInfo_JERDown);
   fChain->SetBranchAddress("JetInfo.GenPartonMatch", &JetInfo_GenPartonMatch, &b_JetInfo_GenPartonMatch);
   fChain->SetBranchAddress("JetInfo.GenPt", &JetInfo_GenPt, &b_JetInfo_GenPt);
   fChain->SetBranchAddress("JetInfo.GenEta", &JetInfo_GenEta, &b_JetInfo_GenEta);
   fChain->SetBranchAddress("JetInfo.GenPhi", &JetInfo_GenPhi, &b_JetInfo_GenPhi);
   fChain->SetBranchAddress("JetInfo.GenPdgID", &JetInfo_GenPdgID, &b_JetInfo_GenPdgID);
   fChain->SetBranchAddress("JetInfo.GenFlavor", &JetInfo_GenFlavor, &b_JetInfo_GenFlavor);
   fChain->SetBranchAddress("JetInfo.GenHadronFlavor", &JetInfo_GenHadronFlavor, &b_JetInfo_GenHadronFlavor);
   fChain->SetBranchAddress("MetInfo.Pt", &MetInfo_Pt, &b_MetInfo_Pt);
   fChain->SetBranchAddress("MetInfo.Phi", &MetInfo_Phi, &b_MetInfo_Phi);
   fChain->SetBranchAddress("MetInfo.Px", &MetInfo_Px, &b_MetInfo_Px);
   fChain->SetBranchAddress("MetInfo.Py", &MetInfo_Py, &b_MetInfo_Py);
   fChain->SetBranchAddress("MetInfo.SumET", &MetInfo_SumET, &b_MetInfo_SumET);
   fChain->SetBranchAddress("MetInfo.CorrPtShiftJetEnUp", &MetInfo_CorrPtShiftJetEnUp, &b_MetInfo_CorrPtShiftJetEnUp);
   fChain->SetBranchAddress("MetInfo.CorrPtShiftJetEnDown", &MetInfo_CorrPtShiftJetEnDown, &b_MetInfo_CorrPtShiftJetEnDown);
   fChain->SetBranchAddress("MetInfo.CorrPtShiftJetResUp", &MetInfo_CorrPtShiftJetResUp, &b_MetInfo_CorrPtShiftJetResUp);
   fChain->SetBranchAddress("MetInfo.CorrPtShiftJetResDown", &MetInfo_CorrPtShiftJetResDown, &b_MetInfo_CorrPtShiftJetResDown);
   fChain->SetBranchAddress("MetInfo.CorrPtShiftUncEnUp", &MetInfo_CorrPtShiftUncEnUp, &b_MetInfo_CorrPtShiftUncEnUp);
   fChain->SetBranchAddress("MetInfo.CorrPtShiftUncEnDown", &MetInfo_CorrPtShiftUncEnDown, &b_MetInfo_CorrPtShiftUncEnDown);
   fChain->SetBranchAddress("MetInfo.CorrPtShiftPhoEnUp", &MetInfo_CorrPtShiftPhoEnUp, &b_MetInfo_CorrPtShiftPhoEnUp);
   fChain->SetBranchAddress("MetInfo.CorrPtShiftPhoEnDown", &MetInfo_CorrPtShiftPhoEnDown, &b_MetInfo_CorrPtShiftPhoEnDown);
   fChain->SetBranchAddress("MetInfo.CorrPhiShiftJetEnUp", &MetInfo_CorrPhiShiftJetEnUp, &b_MetInfo_CorrPhiShiftJetEnUp);
   fChain->SetBranchAddress("MetInfo.CorrPhiShiftJetEnDown", &MetInfo_CorrPhiShiftJetEnDown, &b_MetInfo_CorrPhiShiftJetEnDown);
   fChain->SetBranchAddress("MetInfo.CorrPhiShiftJetResUp", &MetInfo_CorrPhiShiftJetResUp, &b_MetInfo_CorrPhiShiftJetResUp);
   fChain->SetBranchAddress("MetInfo.CorrPhiShiftJetResDown", &MetInfo_CorrPhiShiftJetResDown, &b_MetInfo_CorrPhiShiftJetResDown);
   fChain->SetBranchAddress("MetInfo.CorrPhiShiftUncEnUp", &MetInfo_CorrPhiShiftUncEnUp, &b_MetInfo_CorrPhiShiftUncEnUp);
   fChain->SetBranchAddress("MetInfo.CorrPhiShiftUncEnDown", &MetInfo_CorrPhiShiftUncEnDown, &b_MetInfo_CorrPhiShiftUncEnDown);
   fChain->SetBranchAddress("MetInfo.CorrPhiShiftPhoEnUp", &MetInfo_CorrPhiShiftPhoEnUp, &b_MetInfo_CorrPhiShiftPhoEnUp);
   fChain->SetBranchAddress("MetInfo.CorrPhiShiftPhoEnDown", &MetInfo_CorrPhiShiftPhoEnDown, &b_MetInfo_CorrPhiShiftPhoEnDown);
   fChain->SetBranchAddress("GenPartInfo.size", &GenPartInfo_size, &b_GenPartInfo_size);
   fChain->SetBranchAddress("GenPartInfo.Pt", &GenPartInfo_Pt, &b_GenPartInfo_Pt);
   fChain->SetBranchAddress("GenPartInfo.Eta", &GenPartInfo_Eta, &b_GenPartInfo_Eta);
   fChain->SetBranchAddress("GenPartInfo.Phi", &GenPartInfo_Phi, &b_GenPartInfo_Phi);
   fChain->SetBranchAddress("GenPartInfo.Mass", &GenPartInfo_Mass, &b_GenPartInfo_Mass);
   fChain->SetBranchAddress("GenPartInfo.PdgID", &GenPartInfo_PdgID, &b_GenPartInfo_PdgID);
   fChain->SetBranchAddress("GenPartInfo.Status", &GenPartInfo_Status, &b_GenPartInfo_Status);
   fChain->SetBranchAddress("GenPartInfo.nMo", &GenPartInfo_nMo, &b_GenPartInfo_nMo);
   fChain->SetBranchAddress("GenPartInfo.nDa", &GenPartInfo_nDa, &b_GenPartInfo_nDa);
   fChain->SetBranchAddress("HTXSstage0cat", &HTXSstage0cat, &b_HTXSstage0cat);
   fChain->SetBranchAddress("HTXSstage1cat", &HTXSstage1cat, &b_HTXSstage1cat);
   fChain->SetBranchAddress("HTXSnjets", &HTXSnjets, &b_HTXSnjets);
   fChain->SetBranchAddress("HTXSpTH", &HTXSpTH, &b_HTXSpTH);
   fChain->SetBranchAddress("HTXSpTV", &HTXSpTV, &b_HTXSpTV);
   Notify();
}

Bool_t flashggStdTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void flashggStdTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t flashggStdTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef flashggStdTree_cxx
