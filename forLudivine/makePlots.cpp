// vim: set fdm=marker: 
#include "PlotHelper.h"
#include <iomanip>
// histogram index to sample{{{
// 0: ttH
// 1: Drell-Yan
// 2: Gamma Gamma + Jets
// 3: Gamma + Jets
// 4: QCD
// 5: TT + Gamma Gamma 
// 6: TT + Gamma + Jets
// 7: V + Gamma
// 8: W + Jets
// 9: TT + Jets
// }}}
//std::map<TString, TString> mLabels {{{
std::map<TString, TString> mLabels = {
	{"DY", "Drell-Yan"}, 
	{"DiPhoton", "#gamma#gamma + Jets"},
    //{"DiPhoton", "Multi-jet + #gamma#gamma"},
    {"GammaJets", "#gamma + Jets"},
    //{"GammaJets", "Multi-jet + #gamma"},
    {"GammaJets_Madgraph", "#gamma + Jets (Madgraph)"},
	{"QCD", "QCD"},
	{"TTGG", "t#bar{t} + #gamma#gamma"},
	{"TTGJets", "t#bar{t} + #gamma"},
	{"VG", "V + #gamma"},
	{"WJets", "W+Jets"},
	{"TTJets", "t#bar{t} + jets"}, 
	{"THQ", "tHq"},
	{"TGamma", "t+#gamma+Jets"},
	{"QCD_GammaJets_imputed", "(#gamma) + Jets"},
    //{"QCD_GammaJets_imputed", "Multi-jet (+ #gamma)"},
    {"TTZ", "t#bar{t}Z"},
    {"TTW", "t#bar{t}W"},
        {"VV", "VV"},
        {"tV", "t + V"},
	{"ttH", "t#bar{t}H(125)"},
	{"TT_FCNC_hut", "t#bar{t} FCNC (Hut)"},
	{"TT_FCNC_hct", "t#bar{t} FCNC (Hct)"},
	{"ST_FCNC_hut", "t FCNC (Hut)"},
        {"ST_FCNC_hct", "t FCNC (Hct)"}
};
//}}}
//std::map<TString, int> mColors {{{
std::map<TString, int> mColors = {
        {"DY", kCyan-7},
        {"DiPhoton", kBlue - 4},
        {"GammaJets", kAzure + 1},
	{"GammaJets_Madgraph", kRed},
        {"QCD", kCyan-7},
        {"TTGG", kGreen-2},
        {"TTGJets", kGreen-7},
        {"VG", kViolet-9},
        {"WJets", kBlue+2},
        {"TTJets", kSpring+10},
	{"TGamma", kYellow-9},
	//{"QCD_GammaJets_imputed", kGray},
	{"QCD_GammaJets_imputed", kCyan-9},
	{"TTZ", kAzure+1},
    {"TTW", kRed},
	{"tV", kPink-6},
	{"VV", kPink+6},
    {"ttH", kCyan-9}
};
//}}}
//std::map<TString, TString> mLatex {{{
std::map<TString, TString> mLatex = {
        {"DY", "Drell-Yan"},
        {"DiPhoton", "$\\gamma\\gamma$ + Jets"},
        {"GammaJets", "$\\gamma$ + Jets"},
	{"GammaJets_Madgraph", "$\\gamma$ + Jets (Madgraph)"},
        {"QCD", "QCD"},
        {"TTGG", "$t\\bar{t}+\\gamma\\gamma$"},
        {"TTGJets", "$t\\bar{t}+\\gamma$ + Jets"},
        {"VG", "V + $\\gamma$"},
        {"WJets", "W + Jets"},
        {"TTJets", "$t\\bar{t}$ + Jets"},
        {"TGamma", "$t + \\gamma$"},
	{"QCD_GammaJets_imputed", "($\\gamma$) + Jets (Data Sideband)"},
	{"TTZ", "$t\\bar{t}+Z$"},
    {"TTW", "$t\\bar{t}+W$"},
        {"VV", "$VV$"},
        {"tV", "$tV$"},
	{"ttH", "$t\\bar{t}H$"},
	{"TT_FCNC_hut", "$t\\bar{t}$ FCNC (Hut)"},
        {"TT_FCNC_hct", "$t\\bar{t}$ FCNC (Hct)"},
        {"ST_FCNC_hut", "$t$ FCNC (Hut)"},
        {"ST_FCNC_hct", "$t$ FCNC (Hct)"}
};
//}}}
const int nGenPhotonCats = 3;
const int nGenLeptonCats = 6;
//std::map<int, TString> mPhotons {{{
std::map<int, TString> mPhotons = {
	{0, "(F/F)"},
	{1, "(F/P)"},
	{2, "(P/P)"}
};
//}}}
//std::map<int, TString> mLeptons {{{
std::map<int, TString> mLeptons = {
	{0, "(Had.)"},
	{1, "(e)"},
	{2, "(#mu)"},
	{3, "(e#mu)"},
	{4, "(ee)"},
	{5, "(#mu#mu)"}
};
//}}}
//std::map<int, TString> mLeptonsLatex {{{
std::map<int, TString> mLeptonsLatex = {
        {0, "Hadronic"},
        {1, "$e$"},
        {2, "$\\mu$"},
        {3, "$e\\mu$"},
        {4, "$ee$"},
        {5, "$\\mu\\mu$"}
};
//}}}
//const std::vector<TString> syst_ext {{{
const std::vector<TString> syst_ext = {
    /*
    "MvaShift",
    "SigmaEOverEShift",
    "MaterialCentralBarrel",
    "MaterialOuterBarrel",
    "MaterialForward",
    "FNUFEB",
    "FNUFEE",
    "MCScaleGain6EB",
    "MCScaleGain1EB",
    "ShowerShapeHighR9EB",
    "MCScaleHighR9EB",
    "MCSmearHighR9EBRho",
    "MCSmearHighR9EBPhi",
    "ShowerShapeHighR9EE",
    "MCScaleHighR9EE",
    "MCSmearHighR9EERho",
    "MCSmearHighR9EEPhi",
    "ShowerShapeLowR9EB",
    "MCScaleLowR9EB",
    "MCSmearLowR9EBRho",
    "MCSmearLowR9EBPhi",
    "ShowerShapeLowR9EE",
    "MCScaleLowR9EE",
    "MCSmearLowR9EERho",
    "MCSmearLowR9EEPhi",
    "JEC",
    "JER",
    "PUJIDShift",
    "metJecUncertainty",
    "metJerUncertainty",
    "metPhoUncertainty",
    "metUncUncertainty",
    "UnmatchedPUWeight",
    "MvaLinearSyst",
    "LooseMvaSF",
    "PreselSF",
    "electronVetoSF",
    "TriggerWeight",
    "FracRVWeight",
    "ElectronWeight",
    "JetBTagCutWeight",
    "JetBTagReshapeWeight",
    */

    "MvaShift",
    "JEC",
    "JetBTagReshapeWeight",
    "MCSmearHighR9EBPhi",
    "MCSmearLowR9EBPhi",
    "JER",
    "PreselSF",
    "electronVetoSF",
};
//}}}

void make_plot(TCanvas* c1, TFile* file, string output_name, TString hist_name, TString x_label, vector<TString> vBkgs, vector<TString> vSigs, int idx, TString type = "std", TString year = "2016", bool loose_mva_cut = false, TFile* file_ref = nullptr, vector<TString> vExtraInfo = {}, int yearIdx = -1, bool doSyst = false, bool doRatio = true) {
  // setup{{{
  TString extension = loose_mva_cut ? "MVACategories_1" : "";
  extension = yearIdx == -1 ? "" : (yearIdx == 0 ? "Year_0" : (yearIdx == 1 ? "Year_1" : (yearIdx == 2 ? "Year_2" : "")));

  cout << hist_name + "_Data" + extension << endl;

  TH1D* hData = (TH1D*)file->Get(hist_name + "_Data" + extension);
  TH1D* hData_ref;
  if (file_ref != nullptr) {
    hData_ref = (TH1D*)file_ref->Get(hist_name + "_Data" + extension);
    //hData_ref->Scale(45.966/41.5); 
  }

  /*
  TH1D* hSig_TTH = (TH1D*)file->Get(hist_name + "_ttH" + extension);
  TH1D* hSig_THQ = (TH1D*)file->Get(hist_name + "_THQ" + extension);
  TH1D* hSig_THW = (TH1D*)file->Get(hist_name + "_THW" + extension);
  TH1D* hSig_FCNC_hut = (TH1D*)file->Get(hist_name + "_FCNC_hut" + extension);
  TH1D* hSig_FCNC_hct = (TH1D*)file->Get(hist_name + "_FCNC_hct" + extension);
  vector<TH1D*> hSig;// = {hSig_TTH, hSig_THQ, hSig_THW};
  hSig = {hSig_TTH, hSig_FCNC_hut, hSig_FCNC_hct};
  */
  vector<TH1D*> hSig;

  vector<TH1D*> hBkg;
  vector<TH1D*> hBkgSystUp;
  vector<TH1D*> hBkgSystDown;

  Comparison* c;
  
  vector<TString> vLegendLabels; 
  vector<int> vColors;

  TString output = output_name;
  //}}}
  //if (type.Contains("std")) {{{
  if (type.Contains("std")) {
    if (vSigs.size() == 0) { 
      vLegendLabels = {year + "Data", "ttH (M125)", "FCNC_hut", "FCNC_hct"};
      TH1D* hSig_TTH = (TH1D*)file->Get(hist_name + "_ttH" + extension);
      TH1D* hSig_THQ = (TH1D*)file->Get(hist_name + "_THQ" + extension);
      TH1D* hSig_THW = (TH1D*)file->Get(hist_name + "_THW" + extension);
      TH1D* hSig_FCNC_hut = (TH1D*)file->Get(hist_name + "_FCNC_hut" + extension);
      TH1D* hSig_FCNC_hct = (TH1D*)file->Get(hist_name + "_FCNC_hct" + extension);
      hSig = {hSig_TTH, hSig_FCNC_hut, hSig_FCNC_hct};	
    }
    else {
      vLegendLabels = {year + " Data"};
      for (int i = 0; i < vSigs.size(); i++) {
	    vLegendLabels.push_back(mLabels.find(vSigs[i])->second);
        hSig.push_back((TH1D*)file->Get(hist_name + "_" + vSigs[i] + extension));
      }
            
    }

    for (int i = 0; i < vBkgs.size(); i++) {
      if (vBkgs[i] != "Other" && vBkgs[i] != "Other2") {
        hBkg.push_back((TH1D*)file->Get(hist_name + "_" + vBkgs[i] + extension));
        vLegendLabels.push_back(mLabels.find(vBkgs[i])->second);
        vColors.push_back(mColors.find(vBkgs[i])->second);
      }
      else if (vBkgs[i] == "Other") {
        TH1D* hOther = (TH1D*)file->Get(hist_name + "_" + "DY" + extension);
        vector<TString> rares = { "TGamma", "TTV", "VV", "tV"};
        //if (output.Contains("Hadronic"))
        //    rares.push_back("VG");
        for (unsigned int i = 0; i < rares.size(); i++) {
          hOther->Add((TH1D*)file->Get(hist_name + "_" + rares[i] + extension));
        }
        hBkg.push_back(hOther);
        vLegendLabels.push_back("Other");
        vColors.push_back(kRed);
      }
      else {
        TH1D* hOther = (TH1D*)file->Get(hist_name + "_" + "TGamma" + extension);
        vector<TString> rares = { "TTJets", "TTGG", "TTGJets", "DiPhoton", "VV", "tV", "TTW"};
        if (output.Contains("Hadronic"))
            rares.push_back("QCD_GammaJets_imputed");
        else
            rares.push_back("GammaJets");
        for (unsigned int i = 0; i < rares.size(); i++) {
          hOther->Add((TH1D*)file->Get(hist_name + "_" + rares[i] + extension));
        }
        hBkg.push_back(hOther);
        vLegendLabels.push_back("Other");
        vColors.push_back(kGreen+1);
      }
    }
 
    // Get systs
    vector<TString> syst_bkgs;
    for (unsigned int i = 0; i < vBkgs.size(); i++) {
      if (vBkgs[i] == "Other") {
        syst_bkgs.push_back("TGamma");
        syst_bkgs.push_back("TTW");
        syst_bkgs.push_back("TTZ");
        syst_bkgs.push_back("VV");
        syst_bkgs.push_back("tV");
        syst_bkgs.push_back("DY");
        //if (output.Contains("Hadronic"))
        //    syst_bkgs.push_back("VG");
      }
      else if (vBkgs[i] == "Other2") {
        syst_bkgs.push_back("TGamma");
        syst_bkgs.push_back("VV");
        syst_bkgs.push_back("tV");
        syst_bkgs.push_back("DiPhoton");
        if (output.Contains("Hadronic"))
            syst_bkgs.push_back("QCD_GammaJets_imputed");
        else
            syst_bkgs.push_back("GammaJets");
        syst_bkgs.push_back("TTGG");
        syst_bkgs.push_back("TTGJets");
        syst_bkgs.push_back("TTJets");
        syst_bkgs.push_back("VG");
      }
      else
        syst_bkgs.push_back(vBkgs[i]);
    }

    int n_systs = doSyst? syst_ext.size(): 0;
    for (int i = 0; i < n_systs; i++) {
      for (int j = 0; j < syst_bkgs.size(); j++) {
        TString syst_hist_name_up = "h" + syst_ext[i] + "Up01sigma" + hist_name(1,hist_name.Length());
        TString syst_hist_name_down = "h" + syst_ext[i] + "Down01sigma" + hist_name(1,hist_name.Length());
        TH1D* hSystUp;
        TH1D* hSystDown;
        if (syst_bkgs[j] != "QCD_GammaJets_imputed") {
            cout << syst_hist_name_up + "_" + syst_bkgs[j] + extension << endl;
            hSystUp = (TH1D*)(file->Get(syst_hist_name_up + "_" + syst_bkgs[j] + extension))->Clone("hSystUp" + hist_name + syst_bkgs[j]);
            hSystDown = (TH1D*)(file->Get(syst_hist_name_down + "_" + syst_bkgs[j] + extension))->Clone("hSystDown" + hist_name + syst_bkgs[j]);
        }
        else {
            hSystUp = (TH1D*)(file->Get(hist_name + "_" + syst_bkgs[j] + extension))->Clone("hSystUp" + hist_name + syst_bkgs[j]);
            hSystDown = (TH1D*)(file->Get(hist_name + "_" + syst_bkgs[j] + extension))->Clone("hSystDown" + hist_name + syst_bkgs[j]);
        }
        if (j == 0) {
          hBkgSystUp.push_back(hSystUp);
          hBkgSystDown.push_back(hSystDown);
        }
        else {
          hBkgSystUp[i]->Add(hSystUp);
          hBkgSystDown[i]->Add(hSystDown);
        }
        //delete hSystUp;
        //delete hSystDown;
      }
      /*
      // Insert dummy systematics
      for (int j = 0; j < hBkgSystUp[i]->GetSize(); j++) {
        hBkgSystUp[i]->SetBinContent(j, hBkgSystUp[i]->GetBinContent(j) * 1.0000001);
        hBkgSystDown[i]->SetBinContent(j, hBkgSystDown[i]->GetBinContent(j) * 0.9999999);
      } 
      */
    }

    // print latex table
    if (hist_name == "hNVtx") {
      int n_bins = hData->GetSize()-2;
      double yield_data = hData->Integral(0, n_bins+1);
      double yield_signal = hSig[0]->Integral(0, n_bins+1);
      double yield_bkg = 0;
      vector<double> yield_mc;
      for (int i = 0; i < hBkg.size(); i++) {
        yield_bkg += hBkg[i]->Integral(0, n_bins+1);
        yield_mc.push_back(hBkg[i]->Integral(0, n_bins+1));
      }
      cout.setf(ios::fixed);
      cout << std::setprecision(2) << endl;
      cout << "\\begin{center} \\Fontvi" << endl;
      cout << "\\begin{tabular}{| r || r | r|} \\hline" << endl;
      cout << "Process & Yield & Frac. of total bkg. \\\\ \\hline" << endl;
      for (int i = 0; i < vBkgs.size(); i++) {
        if (i == vBkgs.size() - 1)
          cout << mLatex.find(vBkgs[i])->second << " & " << yield_mc[i] << " & " << yield_mc[i]/yield_bkg << " \\\\ \\hline" << endl;
        else 
          cout << mLatex.find(vBkgs[i])->second << " & " << yield_mc[i] << " & " << yield_mc[i]/yield_bkg << " \\\\" << endl;
      }
      cout << "All bkg. & " << yield_bkg << " & " << yield_bkg/yield_bkg << " \\\\ \\hline" << endl;
      cout << "Data & " << yield_data << " & " << yield_data/yield_bkg << " \\\\ \\hline" << endl;
      cout << "\\end{tabular} \\end{center}" << endl;
    }
  }
  //}}}
  //else if (type == "shape") {{{
  else if (type == "shape") {
    vLegendLabels = {"ttH (M125)"};
    for (int i = 0; i < vBkgs.size(); i++) {
      hBkg.push_back((TH1D*)file->Get(hist_name + "_" + vBkgs[i]));
      vLegendLabels.push_back(mLabels.find(vBkgs[i])->second);
      vColors.push_back(mColors.find(vBkgs[i])->second);
    }
  }
  //}}}
  //else if (type == "GJet_shape") {{{
  else if (type == "GJet_shape") {
    vLegendLabels = {"#gamma + jets (Pythia + MadGraph Reweighted)", "#gamma + jets (Madgraph)"};
    for (int i = 0; i < vBkgs.size(); i++) {
      hBkg.push_back((TH1D*)file->Get(hist_name + "_" + vBkgs[i]));
      vColors.push_back(mColors.find(vBkgs[i])->second);
    }
  }
  //}}}
  //else if (type == "individual_shape") {{{
  else if (type == "individual_shape") {
    hBkg.push_back(hSig[0]);
    vColors = {kBlack};
    vLegendLabels = {"ttH (M125)"}; 
    for (int i = 0; i < vBkgs.size(); i++) {
      hBkg.push_back((TH1D*)file->Get(hist_name + "_" + vBkgs[i] + extension));
      vLegendLabels.push_back(mLabels.find(vBkgs[i])->second);
      vColors.push_back(mColors.find(vBkgs[i])->second);
    }
  }
  //}}}
  //else if (type =="genPhoton") {{{
  else if (type =="genPhoton") {
    vLegendLabels = {"ttH (M125)"};
    vColors = {kBlue+2, kAzure+1, kCyan-7, kYellow, kGreen -4, kTeal + 3, kRed+3, kRed, kMagenta-9};
    for (int i = 0; i < vBkgs.size(); i++) {
      for (int j = 0; j < nGenPhotonCats; j++) {
        hBkg.push_back((TH1D*)file->Get(hist_name + "_" + vBkgs[i] + "GenPhoton_" + to_string(j)));
        vLegendLabels.push_back(mLabels.find(vBkgs[i])->second + mPhotons.find(j)->second);
      }   
    }
    // print latex table
    if (hist_name == "hMass") {
      int n_bins = hData->GetSize()-2;
      int start_bin;
      if (!(output_name == "ttHLeptonic_ttbarCR_plotsgenPhoton.pdf" || output_name == "ttHLeptonic_ttbarCR_plotsstd.pdf"))
        start_bin = 21; // start at m_gg of 80 GeV
      else
	start_bin = 1;
      double yield_data_unc(0), yield_signal_unc(0), yield_bkg_unc(0);
      vector<double> yield_mc_unc;
      double yield_data = hData->IntegralAndError(start_bin, n_bins+1, yield_data_unc);
      double yield_signal = hSig[0]->IntegralAndError(start_bin, n_bins+1, yield_signal_unc);
      double yield_bkg = 0;
      vector<double> yield_mc;
    
      vector<TH1D*> hBkgNominal;
      for (int i = 0; i < vBkgs.size(); i++)
        hBkgNominal.push_back((TH1D*)file->Get(hist_name + "_" + vBkgs[i]));

      for (int i = 0; i < hBkg.size(); i++) {
        yield_bkg += hBkg[i]->IntegralAndError(start_bin, n_bins+1, yield_bkg_unc);
        yield_mc_unc.push_back(0);
        yield_mc.push_back(hBkg[i]->IntegralAndError(start_bin, n_bins+1, yield_mc_unc[i]));
      }
      cout.setf(ios::fixed);
      cout << std::setprecision(2) << endl;
      cout << "\\begin{center} \\Fontvi" << endl;
      cout << "\\begin{tabular}{| c | r || r |} \\hline" << endl;
      cout << "Process & Component & Yield \\\\ \\hline" << endl;
      for (int i = 0; i < vBkgs.size(); i++) {
        double yield_process_unc;
        double yield_process = hBkgNominal[i]->IntegralAndError(start_bin, n_bins+1, yield_process_unc);
        for (int j = 0; j < nGenPhotonCats; j++) {
          double yield_component = yield_mc[(i*nGenPhotonCats)+j];
          if (j == 0)
            cout << "\\multirow{" << nGenPhotonCats+1 << "}{*}{" << mLatex.find(vBkgs[i])->second << "} & " << mPhotons.find(j)->second << " & " << yield_component << " $\\pm$ " << yield_mc_unc[(i*nGenPhotonCats)+j] << " \\\\" << endl;
          else if (j == nGenPhotonCats-1)
            cout << " & " << mPhotons.find(j)->second<< " & " << yield_component << " $\\pm$ " << yield_mc_unc[(i*nGenPhotonCats)+j] << " \\\\ \\cline{2-3}" << endl;
          else
            cout << " & " << mPhotons.find(j)->second<< " & " << yield_component << " $\\pm$ " << yield_mc_unc[(i*nGenPhotonCats)+j] << " \\\\" << endl;
        }
        cout << " & All Components & " << yield_process << " $\\pm$ " << yield_process_unc << " \\\\ \\hline" << endl;

      }
      //cout << "\\multicolumn{2}{|c||}{All background} & " << yield_bkg << " & " << yield_bkg/yield_bkg << " \\\\ \\hline" << endl;
      //cout << "\\multicolumn{2}{|c||}{Signal} & " << yield_signal << " & " << yield_signal/yield_bkg << " \\\\ \\hline" << endl;
      cout << "\\end{tabular} \\end{center}" << endl;

    }
  }
  //}}}
  //else if (type =="genLepton") {{{
  else if (type =="genLepton") {
    vLegendLabels = {"ttH (M125)"};
    vColors = {kBlue+2, kAzure+1, kCyan-7, kYellow, kGreen -4, kTeal + 3, kRed+3, kRed, kRed - 7, kMagenta-9, kOrange, kViolet+1, kOrange-9};
    for (int i = 0; i < vBkgs.size(); i++) {
      for (int j = 0; j < nGenLeptonCats; j++) {
        hBkg.push_back((TH1D*)file->Get(hist_name + "_" + vBkgs[i] + "GenLepton_" + to_string(j)));
        vLegendLabels.push_back(mLabels.find(vBkgs[i])->second + mLeptons.find(j)->second);
      }
    }

    if (hist_name == "hMass") {
      int n_bins = hData->GetSize()-2;
      int start_bin = 21; // start at m_gg of 80 GeV
      double yield_data_unc(0), yield_signal_unc(0), yield_bkg_unc(0);
      vector<double> yield_mc_unc;
      double yield_data = hData->IntegralAndError(start_bin, n_bins+1, yield_data_unc);
      double yield_signal = hSig[0]->IntegralAndError(start_bin, n_bins+1, yield_signal_unc);
      double yield_bkg = 0;
      vector<double> yield_mc;

      vector<TH1D*> hBkgNominal;
      for (int i = 0; i < vBkgs.size(); i++)
        hBkgNominal.push_back((TH1D*)file->Get(hist_name + "_" + vBkgs[i]));   
 
      for (int i = 0; i < hBkg.size(); i++) {
        yield_bkg += hBkg[i]->IntegralAndError(start_bin, n_bins+1, yield_bkg_unc);
	yield_mc_unc.push_back(0);
        yield_mc.push_back(hBkg[i]->IntegralAndError(start_bin, n_bins+1, yield_mc_unc[i]));
      }
      cout.setf(ios::fixed);
      cout << std::setprecision(2) << endl;
      cout << "\\begin{center} \\Fontvi" << endl;
      cout << "\\begin{tabular}{| c | r || r |} \\hline" << endl;
      cout << "Process & Component & Yield \\\\ \\hline" << endl;
      for (int i = 0; i < vBkgs.size(); i++) {
        double yield_process_unc;
        double yield_process = hBkgNominal[i]->IntegralAndError(start_bin, n_bins+1, yield_process_unc);
	for (int j = 0; j < nGenLeptonCats; j++) {
          double yield_component = yield_mc[(i*nGenLeptonCats)+j];
          if (j == 0)
            cout << "\\multirow{" << nGenLeptonCats+1 << "}{*}{" << mLatex.find(vBkgs[i])->second << "} & " << mLeptonsLatex.find(j)->second << " & " << yield_component << " $\\pm$ " << yield_mc_unc[(i*nGenLeptonCats)+j] << " \\\\" << endl;
          else if (j == nGenLeptonCats-1)
            cout << " & " << mLeptonsLatex.find(j)->second<< " & " << yield_component << " $\\pm$ " << yield_mc_unc[(i*nGenLeptonCats)+j] << " \\\\ \\cline{2-3}" << endl;
          else
            cout << " & " << mLeptonsLatex.find(j)->second<< " & " << yield_component << " $\\pm$ " << yield_mc_unc[(i*nGenLeptonCats)+j] << " \\\\" << endl;
        }
        cout << " & All Components & " << yield_process << " $\\pm$ " << yield_process_unc << " \\\\ \\hline" << endl;
      }
      cout << "\\end{tabular} \\end{center}" << endl;
    }
  }
  //}}}
  // lumi{{{
  //TString output = output_name;
  std::map<TString, double> lumi_map = {
	{"2016", 35.9},
	{"2017", 41.5},
	{"2018", 59.76},
	{"RunII", 137.16},
    {"", 137.16}
  };

  double lumi = lumi_map[year];
  cout << "Lumi is " << lumi << endl;
  //double lumi = year == "All" ? 77.4 : (year == "2018" ? 45.996 : ((year == "2017" ? 41.5 : 35.9))); 
  //}}}
  //if (type.Contains("std")) {{{
  if (type.Contains("std")) {
    if (type.Contains("shape")) {
      hSig[0]->Scale(1./hSig[0]->Integral(0, hSig[0]->GetNbinsX()+1));
      if (type.Contains("sig_vs_data")) {
        c = new Comparison(c1, hSig[0], hData);
        c->set_both_data();
        c->set_rat_label("#frac{Signal}{Data}");
      }
      else {
        c = new Comparison(c1, hSig[0], hBkg);
        c->set_rat_label("#frac{Signal}{Background}");
      }
      c->set_scale(-1);
      c->set_y_label("Fraction of events");
    }
    else {
      if (file_ref == nullptr && !doSyst) {
        c = new Comparison(c1, {hData}, hSig, hBkg);
      }
      else if (file_ref == nullptr && doSyst)
        c = new Comparison(c1, {hData}, hSig, hBkg, hBkgSystUp, hBkgSystDown);
      else {
        c = new Comparison(c1, {hData}, hSig, {hData_ref});
        //c = new Comparison(c1, {hData, hData_ref}, hSig, hBkg);
      }
      c->set_data_drawOpt("E");
      c->set_rat_label("#frac{Data}{MC}");
      c->set_y_label("Events");
      if (hist_name.Contains("tthMVA_RunII") && !output.Contains("ttZ"))
          c->set_stack_order({6,5,4,3,2,1,0});
      //else if (hist_name.Contains("tthMVA_RunII"))
      //    c->set_stack_order({3,2,1,0});
      //if (hist_name.Contains("tthMVA_RunII") && output.Contains("Hadronic"))
      //    c->set_stack_order({5,4,3,1,2,0});
      //else if (hist_name.Contains("tthMVA_RunII") && output.Contains("Leptonic"))
      //    c->set_stack_order({6,1,4,5,0,3,2});
    }
    if (hist_name.Contains("htthMVA_RunII_transf") && !output.Contains("ttZ")) {
        vector<double> vlines;
        vector<double> cp_lines;
        if (output.Contains("Leptonic")) {
            vlines = { 0.8997816, 0.95635754, 0.9725133, 0.9870608 }; 
            cp_lines = { 0.9597816 };
            c->set_y_lim_range({0.333, 2*pow(10,5)});
        }
        else if (output.Contains("Hadronic")) {
            vlines = { 0.986025, 0.9948537, 0.9983046, 0.9990729 }; 
            c->set_y_lim_range({0.5, pow(10,6)});
            cp_lines = { 0.99722563 };
        }
        for (unsigned int i = 0; i < vlines.size(); i++) {
            vlines[i] = -log(1-vlines[i]);
        }
        for (unsigned int i = 0; i < cp_lines.size(); i++)
            cp_lines[i] = -log(1-cp_lines[i]);
        c->give_vlines(vlines);
        c->give_vlines_dotted(cp_lines);
        c->give_vshade({0.,vlines[0]});
        c->add_paper_info(output.Contains("Hadronic") ? "Had" : "Lep");
        c->skip_cp(); // uncomment to remove cp lines
    }

    else if (hist_name.Contains("htthMVA_RunII_transf") && hist_name.Contains("ttZ")) {
        vector<double> vlines;
        if (output.Contains("Leptonic")) {
            vlines = { 0.8997816, 0.95635754, 0.9725133, 0.9870608 };
        }
        else if (output.Contains("Hadronic")) {
            vlines = { 0.986025, 0.9948537, 0.9983046, 0.9990729 };
        }
        for (unsigned int i = 0; i < vlines.size(); i++) {
            vlines[i] = -log(1-vlines[i]);
        }
        c->give_vshade({0.,vlines[0]});
        c->skip_signal();
    }

    if (hist_name.Contains("htthMVA_RunII_transf_ttZ") && output.Contains("Leptonic")) {
        c->set_x_bin_range({2,8});
        c->set_y_lim_range({0.5, pow(10,3)});
        c->set_stack_order({2,1,0});
    }
    else if (hist_name.Contains("htthMVA_RunII_transf_ttZ_") && output.Contains("Hadronic") && !(hist_name.Contains("v3") || hist_name.Contains("v4"))) {
        c->set_x_bin_range({3,14});
        c->set_y_lim_range({0.5, pow(10,4)});
        c->set_stack_order({2,1,0});
    }

    if (!doRatio)
        c->set_no_ratio();
    c->set_lumi(lumi);
    //c->set_log_rat();
    c->set_rat_lim_range({0.0, 2.0});
    //if (hist_name.Contains("SigmaEOverE") || hist_name.Contains("DiphotonMassResolution")) {
    //  c->set_no_log();
    //}
    if (type.Contains("linear"))
      c->set_no_log();

  }
  //}}}
  //else if (type == "shape") {{{
  else if (type == "shape") {
    c = new Comparison(c1, hSig[0], hBkg);
    c->set_data_drawOpt("HIST");
    c->set_rat_label("#frac{Signal}{Background}");
    c->set_scale(-1);
    c->set_log_rat();
    c->set_rat_lim_range({0.1, 10.0});
  }
  //}}}
  //else if (type == "GJet_shape") {{{
  else if (type == "GJet_shape") {
    hBkg[0]->Scale(1/hBkg[0]->Integral(0, hBkg[0]->GetSize()-1));
    c = new Comparison(c1, hBkg[0], hBkg[1]);
    c->set_x_label(x_label);
    c->set_y_label("Fraction of Events");
    c->set_no_lumi();
    c->set_scale(-1);
    c->set_rat_label("#frac{Pythia+MadGraph}{MadGraph}");
    c->set_both_data();
    if (output.Contains("wWeights"))
      c->give_info("Pythia + MadGraph Reweighted");
    c->set_y_lim_range({0.005, 3.0});
  }
  //}}}
  //else if (type == "individual_shape") {{{
  else if (type == "individual_shape") {
    c = new Comparison(c1, hBkg);
    c->set_data_drawOpt("HIST");
    c->set_x_label(x_label);
    c->set_y_label("Fraction of Events");
    c->set_scale(-1);
    c->set_no_lumi();
    c->set_no_log();
    c->set_y_lim_range({0.0, 1.0}); 
  }
  //}}}
  //else {{{
  else {
    c = new Comparison(c1, hSig[0], hBkg);
    c->set_data_drawOpt("HIST");
    c->set_rat_label("#frac{Signal}{Background}");
  }
  //}}}
  // set hist{{{
  if (type.Contains("std") && type.Contains("shape"))
      vLegendLabels.erase(vLegendLabels.begin());
  if (type.Contains("std") && type.Contains("shape") && type.Contains("sig_vs_data")) {
      vLegendLabels = {"ttH (M125)", "Data"};
      vColors = {kBlack, kBlack};
  }
  c->set_legend_labels(vLegendLabels);
  c->set_colors(vColors);

  c->set_filename(output_name);
  c->set_x_label(x_label);
  //c->set_y_label("Events");
  //double lumi = year == "All" ? 77.4 : (year == "2018" ? 45.996 : ((year == "2017" ? 41.5 : 35.9)));
  //c->set_lumi(lumi);

  if ((hist_name == "hNJets" || hist_name == "hNbLoose") && !output.Contains("GJet_Reweight")) {
    if (output.Contains("ttHHadronic_2017_Presel")) {
      c->set_no_log();
      c->set_y_lim_range({0.0, 90000});
    }
    //else if (output.Contains("Leptonic"))
    //  c->set_y_lim_range({0.01, pow(10,4)});
    //else if (output.Contains("Hadronic"))
    //  c->set_y_lim_range({0.01, pow(10,4)});
  }

  //if (hist_name.Contains("IDMVA"))
  //  c->set_no_log();

  if (hist_name == "hDiphotonMassResolution")
    c->set_x_bin_range({1, 100});

  if (hist_name.Contains("LeadPToM") || hist_name.Contains("SubleadPToM"))
    c->set_x_bin_range({1,10});

  if (hist_name == "hMassAN") {
    c->set_no_flow();
    c->set_no_log();
    if (output.Contains("SR3") && output.Contains("Hadronic"))
      c->set_y_lim_range({0, 5});
    if (output.Contains("SR2") && output.Contains("Hadronic"))
      c->set_y_lim_range({0, 8});
    if (output.Contains("SR1") && output.Contains("Hadronic"))
      c->set_y_lim_range({0, 25});
    if (output.Contains("SR1") && output.Contains("Leptonic"))
      c->set_y_lim_range({0, 10});
    if (output.Contains("SR2") && output.Contains("Leptonic"))
      c->set_y_lim_range({0, 3.5});
    if (output.Contains("dilep"))
      c->set_y_lim_range({0, 3.5});
    else if (output.Contains("Leptonic") && output.Contains("2016"))
      c->set_y_lim_range({0, 10});
    c->set_x_bin_range({1,80});
    cout << "Data yield in [100,120], [130,180]: " << hData->Integral() << endl;
    cout << "Data yield in [115, 120], [130, 135]: " << hData->Integral(16,20) + hData->Integral(31,35) << endl;
    cout << "Signal yield in [120, 130]: " << hSig[0]->Integral(21,30) << endl;
  }
  

  if (hist_name == "hMass") {
   if (!(output_name == "ttHLeptonic_ttbarCR_plotsgenPhoton.pdf" || output_name == "ttHLeptonic_ttbarCR_plotsstd.pdf") ) {
    c->set_no_underflow();
    //c->set_x_bin_range({21,50});
    }
    if (type == "std" || type == "std_linear")
      c->set_verbose();
  }

  if (hist_name.Contains("SigmaIEtaIEta"))      c->set_x_bin_range({1,50});

  if (hist_name.Contains("htthMVA_RunII_transf_ttZ"))
      c->set_lower_lim(0.1);

  for (int i = 0; i < vExtraInfo.size(); i++)
    c->give_info(vExtraInfo[i]);


  if (type == "individual_shape")
    c->plot(idx, false);
  else
    c->plot(idx);
  delete hData;
  if (type != "individual_shape") {
    //delete hSig;
    for (int i = 0; i < hSig.size(); i++)
      delete hSig[i];
  }
  for (int i = 0; i < hBkg.size(); i++)
    delete hBkg[i];
  delete c;
  //}}}
  cout << "Hello World!" << endl;
}
 

int main(int argc, char* argv[])
{
  // Parse args{{{
  if (argc < 3) {
    cout << "Please provide two arguments: type of plot to make (e.g. 'std') and input file (e.g. '../ttHHadronicLoose_histograms.root')" << endl;
    return 0;
  }

  TString type = argv[1];
  TString type_s = argv[1];

  int yearIdx = -1;
  if (type.Contains("2016"))
    yearIdx = 0;
  else if (type.Contains("2017"))
    yearIdx = 1;
  else if (type.Contains("2018"))
    yearIdx = 2;

  TString file_path = argv[2];
  bool impute_gjets = file_path.Contains("impute") && !(file_path.Contains("presel") || file_path.Contains("sideband")); 
  TString year = file_path.Contains("All") ? "All" : file_path.Contains("2018") ? "2018" : ((file_path.Contains("2017") ? "2017" : "2016"));
  if (file_path.Contains("RunII"))
    year = "";

  if (yearIdx == 0)
    year = "2016";
  else if (yearIdx == 1)
    year = "2017";
  else if (yearIdx == 2)
    year = "2018";
  cout << "Year: " << year << endl;

  TString tag = file_path.Contains("Hadronic") ? "Hadronic" : "Leptonic";

  TString info = argv[3]; 
  TObjArray *tx = info.Tokenize("|");
  vector<TString> vInfo = {};
  for (int i = 0; i < tx->GetEntries(); i++)
    vInfo.push_back(((TObjString *)(tx->At(i)))->String());

  //TString file_path_ref = argv[4];
  //TString year_ref = file_path_ref.Contains("RunII") ? "2017" : file_path_ref.Contains("2018") ? "2018" : ((file_path_ref.Contains("2017") ? "2017" : "2016"));

  bool doSyst = false;
  bool doRatio = true;
  bool loose_mva_cut = false; //argc > 4;
  TString mva_ext = loose_mva_cut ? "_looseMVACut" : "";

  TFile* f = new TFile(file_path);
  vector<TFile*> vFiles = {f};
  string output = (file_path.ReplaceAll("../", "")).ReplaceAll(".root", type + mva_ext + ".pdf").Data();

  TFile* f_ref; 
  f_ref = nullptr;
  //if (argc > 4) f_ref = new TFile(file_path_ref);
  //else f_ref = nullptr;
  
  vector<string> vNames = {output};

  // Decide which backgrounds you want to plot
  //if (type.Contains("std_linear"))
  //  type = "std_linear";
  //else if (type.Contains("std"))
  //  type = "std";

  vector<TString> vSigs = {"ttH"};
  vector<TString> vBkgs;
  if (type == "std" || type == "shape" || type == "std_linear") { 
    vBkgs = {"DiPhoton", "GammaJets", "TTGG", "TTGJets", "TTJets", "DY", "VG", "TGamma", "TTV", "VV", "tV"};
    if (file_path.Contains("impute"))
      vBkgs = {"DiPhoton", "QCD_GammaJets_imputed", "TTGG", "TTGJets", "TTJets", "DY", "VG", "TGamma", "TTV", "VV", "tV"};
  }
  else if (type == "individual_shape") {
    vBkgs = {"DiPhoton", "GammaJets", "TTGG", "TTGJets"};
  }
  else if (type == "genPhoton") {
    //vBkgs = {"DiPhoton", "GammaJets", "QCD"};
    //vBkgs = {"TTGG", "TTGJets", "TTJets"};
    vBkgs = {"TTGG", "TTGJets"};
  } 
  else if (type == "genLepton") {
    //vBkgs = {"TTGG", "TTGJets", "TTJets"};
    //vBkgs = {"TTGJets"};
    //vBkgs = {"TTGG"};
    vBkgs = {"TTJets"};
    //vBkgs = {"TTGG", "TTGJets"};
  }

  else if (type == "GJet_shape") {
    vBkgs = {"GammaJets", "GammaJets_Madgraph"};
  }

  

  if (argc > 4) {
    TString bkg_list = argv[4];
    TObjArray *t_bkg = bkg_list.Tokenize("|");
    vBkgs = {};
    for (int i = 0; i < t_bkg->GetEntries(); i++)
      vBkgs.push_back(((TObjString *)(t_bkg->At(i)))->String());
    cout << "Plotting the following as bkgs: ";
    for (int i = 0; i < vBkgs.size(); i++)
      cout << vBkgs[i] << ", ";
    cout << endl;
  }

  if (argc > 5) {
    TString sig_list = argv[5];
    TObjArray *t_sig = sig_list.Tokenize("|");
    vSigs = {};
    for (int i = 0; i < t_sig->GetEntries(); i++)
      vSigs.push_back(((TObjString *)(t_sig->At(i)))->String());
    cout << "Plotting the following as signals: "; 
    for (int i = 0; i < vSigs.size(); i++)
      cout << vSigs[i] << ", ";
    cout << endl;
  }
  //}}}
  // Style options{{{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow,0);
  gStyle->SetPaintTextFormat(".2f");
  gStyle->SetTickLength(0.01);
  //gStyle->SetErrorX(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTitleAlign(33);
  gStyle->SetTitleX(.93);
  gStyle->SetTitleFontSize(0.03); 

  TCanvas* c1 = new TCanvas("c1", "histos", 600, 800);
  TCanvas* c2 = new TCanvas("c2", "histos2", 800, 800);
  TCanvas* c3 = new TCanvas("c3", "histos3", 1000, 800);
  //}}}
  // Make plots{{{
  for (int i = 0; i < vFiles.size(); i++) {
    make_plot(c1, vFiles[i], vNames[i], "hMass", "m_{#gamma#gamma} [GeV]", vBkgs, vSigs, 0,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //others{{{
    make_plot(c1, vFiles[i], vNames[i], "hMass_v2", "m_{#gamma#gamma} [GeV]", vBkgs, vSigs, 0,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hRapidity", "Y_{#gamma#gamma} [GeV^{1/2}]", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hDiphotonSumPt", "p_{T}(#gamma_{1}) + p_{T}(#gamma_{2}) [GeV]", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hDiphotonCosPhi", "|cos(#Delta #phi_{#gamma 1, #gamma 2})|", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);


    //make_plot(c1, vFiles[i], vNames[i], "hMassTop1", "m_{#gamma#gammaj} [GeV]", vBkgs, vSigs, 0,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hMassTop2", "m_{jjj} [GeV]", vBkgs, vSigs, 0,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    make_plot(c1, vFiles[i], vNames[i], "hNJets", "N_{jets}", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hNbLoose", "N_{b-jets} (Loose)", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hNbMedium", "N_{b-jets} (Medium)", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);   
    make_plot(c1, vFiles[i], vNames[i], "hNbTight", "N_{b-jets} (Tight)", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    make_plot(c1, vFiles[i], vNames[i], "hJet1pT", "Jet1 p_{T} [GeV]", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hJet2pT", "Jet2 p_{T} [GeV]", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hJet3pT", "Jet3 p_{T} [GeV]", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hJet4pT", "Jet4 p_{T} [GeV]", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hbJet1pT", "bJet1 p_{T} [GeV]", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hbJet2pT", "bJet2 p_{T} [GeV]", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    make_plot(c1, vFiles[i], vNames[i], "hPhotonLeadPt", "p_{T}(#gamma_{1}) [GeV]", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonLeadEta", "#eta(#gamma_{1})", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonLeadPhi", "#phi(#gamma_{1})", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonLeadSigmaIEtaIEta", "#sigma_{i#eta,i#eta}(#gamma_{1})", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonLeadHOverE", "h/E(#gamma_{1})", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonLeadR9", "R9(#gamma_{1})", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonLeadIDMVA", "Photon ID MVA(#gamma_{1})", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonLeadPToM", "p_{T}/m_{#gamma#gamma} (#gamma_{1})", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonLeadSigmaEOverE", "#sigma_{E}/E (#gamma_{1})", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    make_plot(c1, vFiles[i], vNames[i], "hPhotonSubleadPt", "p_{T}(#gamma_{2}) [GeV]", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonSubleadEta", "#eta(#gamma_{2})", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonSubleadPhi", "#phi(#gamma_{2})", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonSubleadSigmaIEtaIEta", "#sigma_{i#eta,i#eta}(#gamma_{2})", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonSubleadHOverE", "h/E(#gamma_{2})", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonSubleadR9", "R9(#gamma_{2})", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonSubleadIDMVA", "Photon ID MVA(#gamma_{2})", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonSubleadPToM", "p_{T}/m_{#gamma#gamma} (#gamma_{2})", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonSubleadSigmaEOverE", "#sigma_{E}/E (#gamma_{2})", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    make_plot(c1, vFiles[i], vNames[i], "htthMVA", "tth MVA", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hMaxBTag", "max b-tag response", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hSecondMaxBTag", "2nd max b-tag response", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    make_plot(c1, vFiles[i], vNames[i], "hJet1Eta", "Jet1 #eta", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hJet2Eta", "Jet2 #eta", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hJet3Eta", "Jet3 #eta", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hJet4Eta", "Jet4 #eta", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    make_plot(c1, vFiles[i], vNames[i], "hHT", "HT [GeV]", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hMetPt", "E_{T}^{miss} [GeV]", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    if (vNames[i] == "ttHLeptonic_ttbarCR_plots" + type_s + ".pdf") {
      make_plot(c1, vFiles[i], vNames[i], "hPhotonMaxIDMVA", "Max #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hPhotonMinIDMVA", "Min #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hPhotonMinIDMVA_coarse", "Min #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hDiphoMVA", "Diphoton MVA", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    }
    //make_plot(c2, vFiles[i], vNames[i], "hMassAN", "m_{#gamma#gamma} [GeV]", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hMassAN", "m_{#gamma#gamma} [GeV]", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hPhotonIDMVA_prompt", "#gamma ID (Prompt)", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hPhotonIDMVA_elec", "#gamma ID (Elec)", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hPhotonIDMVA_fake", "#gamma ID (Fake)", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    make_plot(c1, vFiles[i], vNames[i], "hMT", "m_{T}(E_{T}^{miss}, lep) [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonMaxIDMVA_fine", "Max #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonMinIDMVA_fine", "Min #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonMinIDMVA_coarse", "Min #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonMaxIDMVA_coarse", "Max #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hDiphoMVA", "Diphoton MVA", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    //make_plot(c1, vFiles[i], vNames[i], "h" + tag + "MVA", tag + " MVA", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    if (tag == "Leptonic") {
      make_plot(c1, vFiles[i], vNames[i], "hLeptonPt", "Lepton p_{T} [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hLeptonEta", "Lepton #eta", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hMuonMiniIsolation", "Muon Mini-Iso", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    }

    make_plot(c1, vFiles[i], vNames[i], "hPtHiggs", "DiPhoton p_{T} [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hMinDrDiphoJet", "Min #Delta R(p_{#gamma#gamma}, jet)", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    if (tag == "Leptonic") {
      /*
      make_plot(c1, vFiles[i], vNames[i], "hDijetClosestWMass", "Min |m_{jj} - m_{W}| [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hDijetMass", "m_{jj} (all pairs)", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hDeltaRDiphoW", "#Delta R(p_{#gamma#gamma}, p_{W})", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hDeltaRDiphoLep", "#Delta R(p_{#gamma#gamma}, p_{lep})", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hTopPt", "Hadronic Top p_{T} [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hTopMass", "Hadronic Top Mass [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hTopEta", "Hadronic Top Eta [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hDeltaRDiphoTop", "#Delta R(p_{#gamma#gamma}, p_{top (had.)})", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      */

      make_plot(c1, vFiles[i], vNames[i], "hPhotonDeltaR", "#Delta R(#gamma_{1}, #gamma_{2})", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

      make_plot(c1, vFiles[i], vNames[i], "hMaxBTag", "Max b-tag", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hSecondMaxBTag", "2nd Max b-tag", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

      make_plot(c1, vFiles[i], vNames[i], "hJet5pT", "Jet5 p_{T} [GeV]", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);     
      make_plot(c1, vFiles[i], vNames[i], "hJet5Eta", "Jet5 #eta", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

      make_plot(c1, vFiles[i], vNames[i], "hNLepLoose", "N_{lep} (loose ID)", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hNLepMedium", "N_{lep} (medium ID)", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hNLepTight", "N_{lep} (tight ID)", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    }

    //if (vNames[i] == "ttHHadronicLoose_plots_" + type_s + ".pdf" || vNames[i] == "ttHHadronicCustom_plots_" + type_s + ".pdf")
    //  make_plot(c1, vFiles[i], vNames[i], "hHadronicMVA", "Hadronic MVA Score", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio); 

    //if (vNames[i] == "ttHLeptonicLoose_plots_" + type_s + ".pdf" || vNames[i] == "ttHLeptonicCustom_plots_" + type_s + ".pdf")
    //make_plot(c1, vFiles[i], vNames[i], "hLeptonicMVA", "Leptonic MVA Score", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    //make_plot(c1, vFiles[i], vNames[i], "hLeadMinDr", "Min #Delta R(#gamma_1, leps/jets)", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hSubleadMinDr", "Min #Delta R(#gamma_2, leps/jets)", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    make_plot(c1, vFiles[i], vNames[i], "hAbsCosHelicity", "|cos(#theta_{helicity})| (p_{#gamma #gamma})", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    
    if (tag == "Leptonic") {
      //make_plot(c1, vFiles[i], vNames[i], "hPhotonMinIDMVA_passPSV", "Min #gamma ID (pass PSV)", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      //make_plot(c1, vFiles[i], vNames[i], "hPhotonMinIDMVA_failPSV", "Min #gamma ID (fail PSV)", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    }

    //make_plot(c1, vFiles[i], vNames[i], "hPixelSeed", "Pixel Seed Veto", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hPixelSeedEB", "Pixel Seed Veto (EB)", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hPixelSeedEE", "Pixel Seed Veto (EE)", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    //make_plot(c1, vFiles[i], vNames[i], "hDiphotonMassResolution", "#sigma_{m_{#gamma#gamma}} / m_{#gamma#gamma}", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    if (tag == "Hadronic") {
      /*  
      make_plot(c1, vFiles[i], vNames[i], "hDiphotonMassResolutionHighMVA", "#sigma_{m_{#gamma#gamma}} / m_{#gamma#gamma}", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hDiphotonMassResolutionMedMVA", "#sigma_{m_{#gamma#gamma}} / m_{#gamma#gamma}", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
      make_plot(c1, vFiles[i], vNames[i], "hDiphotonMassResolutionLowMVA", "#sigma_{m_{#gamma#gamma}} / m_{#gamma#gamma}", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

      make_plot(c1, vFiles[i], vNames[i], "hPhotonMaxIDMVA_NJets2", "Max #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, {"N_{jets} == 2"});
      make_plot(c1, vFiles[i], vNames[i], "hPhotonMinIDMVA_NJets2", "Min #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, {"N_{jets} == 2"});     
      make_plot(c1, vFiles[i], vNames[i], "hPhotonMaxIDMVA_NJets3", "Max #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, {"N_{jets} == 3"});
      make_plot(c1, vFiles[i], vNames[i], "hPhotonMinIDMVA_NJets3", "Min #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, {"N_{jets} == 3"});
      make_plot(c1, vFiles[i], vNames[i], "hPhotonMaxIDMVA_NJets4+", "Max #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, {"N_{jets} #geq 4"});
      make_plot(c1, vFiles[i], vNames[i], "hPhotonMinIDMVA_NJets4+", "Min #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, {"N_{jets} #geq 4"});
      */
    }

    make_plot(c2, vFiles[i], vNames[i], "hNJets", "N_{jets}", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    if (file_path.Contains("GJet_Reweight_Preselection"))
      make_plot(c1, vFiles[i], vNames[i], "hGJet_BDT", "#gamma + jets BDT Score", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    //make_plot(c1, vFiles[i], vNames[i], "hPhotonMinIDMVA_coarse_0b", "Min #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hPhotonMaxIDMVA_coarse_0b", "Max #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    make_plot(c1, vFiles[i], vNames[i], "hTopTagger_score", "Top tagger score", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    //make_plot(c1, vFiles[i], vNames[i], "hMinIDPhotonPt", "Min. ID #gamma p_{T} [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hMinIDPhotonEta", "Min. ID #gamma #eta", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hMaxIDPhotonPt", "Max. ID #gamma p_{T} [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hMaxIDPhotonEta", "Max. ID #gamma #eta", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    make_plot(c1, vFiles[i], vNames[i], "hPhotonDeltaR", "#Delta R(#gamma_{1}, #gamma_{2})", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hDiphotonPtOverMass", "p_{T}^{#gamma#gamma} / m_{#gamma#gamma}", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    make_plot(c1, vFiles[i], vNames[i], "hPhotonLeadPixelSeed", "Lead #gamma PSV", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonSubleadPixelSeed", "Sublead #gamma PSV", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    make_plot(c1, vFiles[i], vNames[i], "hJet1BTag", "Jet 1 b-tag score", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hJet2BTag", "Jet 2 b-tag score", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hJet3BTag", "Jet 3 b-tag score", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hJet4BTag", "Jet 4 b-tag score", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
 
    //make_plot(c1, vFiles[i], vNames[i], "hHadronicMVA_coarse", "Hadronic MVA Score", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hHadronicMVA_fine", "Hadronic MVA Score", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio); 

    //make_plot(c1, vFiles[i], vNames[i], "htthMVA_RunII", "MVA Score", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    /*
    vInfo.push_back("N_{jets} #geq 5");
    make_plot(c1, vFiles[i], vNames[i], "hPhotonMaxIDMVA_NJets5+", "Max #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonMinIDMVA_NJets5+", "Min #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio); 
    vInfo[vInfo.size()-1] = "N_{jets} #geq 7"; 
    make_plot(c1, vFiles[i], vNames[i], "hPhotonMaxIDMVA_NJets7+", "Max #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hPhotonMinIDMVA_NJets7+", "Min #gamma ID", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio); 

    vInfo[vInfo.size()-1] = "Evts Passing p_{T}/m_{#gamma#gamma} Cuts";
    make_plot(c1, vFiles[i], vNames[i], "hMass_PassPtToM", "m_{#gamma#gamma} [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    vInfo[vInfo.size()-1] = "Evts Failing p_{T}/m_{#gamma#gamma} Cuts";
    make_plot(c1, vFiles[i], vNames[i], "hMass_FailPtToM", "m_{#gamma#gamma} [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    vInfo[vInfo.size()-1] = "After Cut on MVA";
    vInfo.push_back("");    
    vInfo[vInfo.size()-1] = "Evts Passing p_{T}/m_{#gamma#gamma} Cuts";
    make_plot(c1, vFiles[i], vNames[i], "hMass_PassPtToM_AfterBDTCut", "m_{#gamma#gamma} [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    vInfo[vInfo.size()-1] = "Evts Failing p_{T}/m_{#gamma#gamma} Cuts";
    make_plot(c1, vFiles[i], vNames[i], "hMass_FailPtToM_AfterBDTCut", "m_{#gamma#gamma} [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);


    vInfo[vInfo.size()-1] = "";

    make_plot(c1, vFiles[i], vNames[i], "hMassTop_Hq_1", "m_{#gamma#gammaq,1} [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hMassTop_Hq_2", "m_{#gamma#gammaq,2} [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hMassTop_Hq_3", "m_{#gamma#gammaq,3} [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    make_plot(c1, vFiles[i], vNames[i], "hMassTop_qqq_1", "m_{bqq,1} [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hMassTop_qqq_2", "m_{bqq,2} [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //make_plot(c1, vFiles[i], vNames[i], "hMassTop_qqq_3", "m_{bqq,3} [GeV]", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);

    vInfo.erase(vInfo.end() - 2, vInfo.end());
    */


    if (!(file_path.Contains("hct") || file_path.Contains("hut"))) {
        make_plot(c1, vFiles[i], vNames[i], "htthMVA_RunII_transf", "BDT-bkg", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
        
        make_plot(c1, vFiles[i], vNames[i], "hDNNScore_ttH_vs_ttGG", "DNN Score (t#bar{t}H vs. t#bar{t} + #gamma#gamma)", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
        if (tag == "Hadronic")
           make_plot(c1, vFiles[i], vNames[i], "hDNNScore_ttH_vs_dipho", "DNN Score (t#bar{t}H vs. #gamma#gamma + jets)", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio); 

        make_plot(c1, vFiles[i], vNames[i], "htthMVA_RunII_transf_ttZ", "BDT-bkg", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
        make_plot(c1, vFiles[i], vNames[i], "htthMVA_RunII_transf_ttZ_v2", "BDT-bkg", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
        make_plot(c1, vFiles[i], vNames[i], "htthMVA_RunII_transf_ttZ_v3", "BDT-bkg", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
        make_plot(c1, vFiles[i], vNames[i], "htthMVA_RunII_transf_ttZ_v4", "BDT-bkg", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
        make_plot(c1, vFiles[i], vNames[i], "htthMVA_RunII_transf_ttZ_v5", "BDT-bkg", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
        make_plot(c3, vFiles[i], vNames[i], "htthMVA_RunII_transf", "BDT-bkg", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, false);
        //make_plot(c3, vFiles[i], vNames[i], "htthMVA_RunII_transf_bounded_v2", "BDT-bkg", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, false);
    }
    if (file_path.Contains("hct") || file_path.Contains("hut")) 
        make_plot(c1, vFiles[i], vNames[i], "hMVA_transf", "MVA Score", vBkgs, vSigs, 1, type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst);
    make_plot(c1, vFiles[i], vNames[i], "hRho", "Rho", vBkgs, vSigs, 1,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    make_plot(c1, vFiles[i], vNames[i], "hNVtx", "# Vertices", vBkgs, vSigs, 2,type, year, loose_mva_cut, f_ref, vInfo, yearIdx, doSyst, doRatio);
    //}}}
  }
  //}}}
}
