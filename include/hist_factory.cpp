#include "hist_factory.h"
#include <string>
using namespace std;
//#include "hist_components.cpp"
        //component_plots chist_deltaR_wboson_topKinFit;
        //chist_deltaR_gen_reco_wboson_topKinFit.Init("ntuples_skimmed", "chist_deltaR_gen_reco_wboson_topKinFit", "regDisc", "irrDisc", 40, 0, 6);
    //    chist_deltaR_gen_reco_wboson_topKinFit.Fill_hist(canFoundSolution_topKinFit, deltaR);
    //chist_deltaR_gen_reco_wboson_topKinFit.Report_hist("regular", "irregular");
    //chist_deltaR_gen_reco_wboson_topKinFit.Draw_hist(c1, "regular", "irregular");

hist_factory::hist_factory(const char* dirName, const char* objName, const char* tag, double _mass_scope_upper_bound)
{
    sprintf(_dirName, "%s", dirName);
    sprintf(_objName, "%s", objName);
    double _mass = _mass_scope_upper_bound;

    legend = new TLegend(0.64,0.66,0.86,0.86);

    char _w[4] = "W";
    char _t[4] = "t";
    char _particle[16] = "";
    char _titles[128] = "";
    if( (string) _objName == "wboson") sprintf(_particle, "%s", _w);
    if( (string) _objName == "top")    sprintf(_particle, "%s", _t);
    if( (string) _objName == "tbw")    sprintf(_particle, "%s", _t);
    if( (string) _objName == "tqh")    sprintf(_particle, "%s", _t);
    //printf("[check] %s\n", _particle);

    // basic quantities
    hist_pz      = new TH1D(Form("hist_%s_pz_%s"      , objName , tag) , ";P_{z} [GeV/c];Entries"      , 40 , -200    , 200);
    hist_pt      = new TH1D(Form("hist_%s_pt_%s"      , objName , tag) , ";P_{T} [GeV/c];Entries"      , 20 , 0    , 200);
    hist_eta     = new TH1D(Form("hist_%s_eta_%s"     , objName , tag) , ";#eta;Entries"               , 20 , -2.5 , 2.5);
    hist_phi     = new TH1D(Form("hist_%s_phi_%s"     , objName , tag) , ";#phi;Entries"               , 20 , -3.0 , 3.0);
    hist_mass    = new TH1D(Form("hist_%s_mass_%s"    , objName , tag) , ";Mass [GeV/c^{2}];Entries"   , 100 , 0    , _mass);
    hist_gen_mass    = new TH1D(Form("gen_hist_%s_mass_%s"    , objName , tag) , ";Mass [GeV/c^{2}];Entries"   , 100 , 0    , _mass);

    // check reco-gen quantities
    sprintf(_titles, "%s", Form(";M(%s)^{gen} - M(%s)^{reco} [GeV/c^{2}];Entries", _particle, _particle));
    hist_deltaM  = new TH1D(Form("hist_deltaM_%s_%s"  , objName , tag) , _titles , 40 , -200 , 200);
    hist_deltaR  = new TH1D(Form("hist_deltaR_%s_%s"  , objName , tag) , ";deltaR;Entries"             , 40 , 0    , 6);
    hist_deltaPT = new TH1D(Form("hist_deltaPT_%s_%s" , objName , tag) , ";P_{T}^{gen} - P_{T}^{reco} [GeV/c];Entries"    , 40 , -200 , 200);
    hist_deltaPz = new TH1D(Form("hist_deltaPz_%s_%s" , objName , tag) , ";Pz^{gen} - Pz^{reco} [GeV/c];Entries"    , 40 , -200 , 200);

    sprintf(_titles, "%s", Form(";(M(%s)^{gen} - M(%s)^{reco}) / M(%s)^{gen};Entries", _particle, _particle, _particle));
    hist_deltaM_ratio  = new TH1D(Form("hist_ratio_deltaM_%s_%s"  , objName , tag) , _titles , 40 , -4 , 4);
    hist_deltaPT_ratio = new TH1D(Form("hist_ratio_deltaPT_%s_%s" , objName , tag) , ";(P_{T}^{gen} - P_{T}^{reco}) / P_{T}^{gen};Entries" , 40 , -4 , 4);
    hist_deltaPz_ratio = new TH1D(Form("hist_ratio_deltaPz_%s_%s" , objName , tag) , ";(Pz^{gen} - Pz^{reco}) / Pz^{gen};Entries" , 40 , -4 , 4);

    sprintf(_titles, "%s", Form(";M(%s)^{gen} [GeV/c^{2}];(M(%s)^{gen} - M(%s)^{reco}) [GeV/c^{2}]", _particle, _particle, _particle));
    hist_deltaM_genM      = new TH2D(Form("hist_2D_deltaM_genM_%s_%s"      , objName , tag) , _titles , 40 , 0 , _mass , 40 , -100 , 100);
    sprintf(_titles, "%s", Form(";M(%s)_mean^{gen} [GeV/c^{2}];(M(%s)^{gen} - M(%s)^{reco}) [GeV/c^{2}]", _particle, _particle, _particle));
    hist_deltaM_genM_mean = new TH2D(Form("hist_2D_deltaM_genM_mean_%s_%s" , objName , tag) , _titles , 40 , 0 , _mass , 40 , -100 , 100);
    sprintf(_titles, "%s", Form(";M(%s)_rms^{gen} [GeV/c^{2}];(M(%s)^{gen} - M(%s)^{reco}) [GeV/c^{2}]", _particle, _particle, _particle));
    hist_deltaM_genM_rms  = new TH2D(Form("hist_2D_deltaM_genM_rms_%s_%s"  , objName , tag) , _titles , 40 , 0 , _mass , 40 , -100 , 100);

    hist_recPz_genPz      = new TH2D(Form("hist_2D_recPz_genPz_%s_%s"      , objName , tag) , ";Pz^{gen} [GeV/c];Pz^{reco} [GeV/c]" , 40 , -200 , 200 , 40 , -200 , 200);
}

hist_factory::~hist_factory()
{
}

void hist_factory::Fill_hist(TLorentzVector obj_reco, TLorentzVector obj_gen)
{
    double _pz      = obj_reco.Pz();
    double _pt      = obj_reco.Pt();
    double _eta     = obj_reco.Eta();
    double _phi     = obj_reco.Phi();
    double _mass    = obj_reco.M();
    double _gen_mass= obj_gen.M();
    double _deltaR  = obj_reco.DeltaR(obj_gen);
    double _deltaM  = obj_gen.M() - obj_reco.M();
    double _deltaPT = obj_gen.Pt() - obj_reco.Pt();
    double _deltaPz = obj_gen.Pz() - obj_reco.Pz();
    hist_pz      -> Fill(_pz);
    hist_pt      -> Fill(_pt);
    hist_eta     -> Fill(_eta);
    hist_phi     -> Fill(_phi);
    hist_mass    -> Fill(_mass);
    hist_gen_mass-> Fill(_gen_mass);
    hist_deltaR  -> Fill(_deltaR);
    hist_deltaM  -> Fill(_deltaM);
    hist_deltaPT -> Fill(_deltaPT);
    hist_deltaPz -> Fill(_deltaPz);

    double _genM = obj_gen.M();
    double _genPT = obj_gen.Pt();
    double _genPz = obj_gen.Pz();
    double _deltaM_ratio = _deltaM / _genM;
    double _deltaPT_ratio = _deltaPT / _genPT;
    double _deltaPz_ratio = _deltaPz / _genPz;
    hist_deltaM_ratio  -> Fill(_deltaM_ratio);
    hist_deltaPT_ratio -> Fill(_deltaPT_ratio);
    hist_deltaPz_ratio -> Fill(_deltaPz_ratio);

    double _recPz = _pz;
    double _genM_mean = obj_gen.M();
    double _genM_rms = obj_gen.M();
    hist_recPz_genPz      -> Fill(_genPz, _recPz);
    hist_deltaM_genM      -> Fill(_genM, _deltaM);
    hist_deltaM_genM_mean -> Fill(_genM_mean, _deltaM);
    hist_deltaM_genM_rms  -> Fill(_genM_rms, _deltaM);
}

int hist_factory::Get_histEntries()
{
    int entries = hist_pt->GetEntries();
    return entries;
}

void hist_factory::Draw_individual_hist(TCanvas *c1, TH1D *hist)
{
    hist -> Draw("hist");
    hist -> SetMinimum(0.);
    hist -> SetLineWidth(2);
    c1->SaveAs(Form("%s/%s.png", _dirName, hist->GetName()));
}

void hist_factory::Draw_individual_hist(TCanvas *c1, TH2D *hist, const char* option)
{
    gPad -> SetRightMargin(0.10);
    if(option == "colz")
    {
        gPad -> SetRightMargin(0.15);
        hist -> SetStats(0); // remove stats for better display
    }
    hist -> GetYaxis() -> SetTitleOffset(1.05);
    hist -> Draw(option);
    c1->SaveAs(Form("%s/%s.png", _dirName, hist->GetName()));

    //printf("[INFO-hist_factory::Draw_individual_hist] Correlation Factor of %s: %f\n", hist->GetName(), hist->GetCorrelationFactor());
}

void hist_factory::Draw_gen_reco(TCanvas *c1, TH1D *h1, TH1D *h2)
{
    double max1 = h1 -> GetMaximum();
    double max2 = h2 -> GetMaximum();
    double max = (max1>max2) ? max1 : max2;
    h1 -> SetStats(0);
    h1 -> SetMaximum(max*1.2);
    //h1 -> SetMinimum(1);
    //h2 -> SetMinimum(1);
    h1 -> Draw("hist");
    h2 -> Draw("hist, same");
    h1 -> SetLineWidth(2);
    h2 -> SetLineWidth(2);
    h2 -> SetLineColor(2);

    legend->Clear();
    legend->AddEntry(h1, "Reco", "l");
    legend->AddEntry(h2, "Gen-level", "l");
    legend->SetLineColor(0);
    legend->Draw("same");

    //gPad->SetLogy(1);

    c1->SaveAs(Form("%s/gen_reco_%s.png", _dirName, h1->GetName()));
}

void hist_factory::Draw_all_hist(TCanvas *c1)
{
    Draw_individual_hist(c1, hist_pz);
    Draw_individual_hist(c1, hist_pt);
    Draw_individual_hist(c1, hist_eta);
    Draw_individual_hist(c1, hist_phi);
    Draw_individual_hist(c1, hist_mass);
    Draw_individual_hist(c1, hist_gen_mass);
    Draw_individual_hist(c1, hist_deltaR);
    Draw_individual_hist(c1, hist_deltaM);
    Draw_individual_hist(c1, hist_deltaPT);
    Draw_individual_hist(c1, hist_deltaPz);

    Draw_individual_hist(c1, hist_deltaM_ratio);
    Draw_individual_hist(c1, hist_deltaPT_ratio);
    Draw_individual_hist(c1, hist_deltaPz_ratio);

    Draw_individual_hist(c1, hist_recPz_genPz, "colz");
    Draw_individual_hist(c1, hist_deltaM_genM, "box");
    Draw_individual_hist(c1, hist_deltaM_genM_mean, "box");
    Draw_individual_hist(c1, hist_deltaM_genM_rms, "box");

    Draw_gen_reco(c1, hist_mass, hist_gen_mass);
}

double GetMaxScope(TH1D* h1, TH1D* h2){
    double max1 = h1 -> GetMaximum();
    double max2 = h2 -> GetMaximum();
    double max = (max1>max2) ? max1 : max2;
    return max;
}

void PlotTogether(TCanvas* c1, TLegend *legend, TH1D *h1, TH1D *h2, const char* path)
{
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.15);

    double scale = 1.2;
    double max = GetMaxScope(h1, h2);
    h1->SetStats(0);
    h1->SetLineWidth(2);
    h1->SetLineColor(kBlue);
    h1->SetMaximum(max*scale);
    h1->Draw();

    h2->SetStats(0);
    h2->SetLineWidth(2);
    h2->SetLineStyle(2);
    h2->SetLineColor(kRed);
    h2->Draw("same");

    legend->Clear();
    legend->AddEntry(h1, "quadratic", "l");
    legend->AddEntry(h2, "topKinFit", "l");
    legend->SetLineColor(0);
    legend->Draw("same");

    c1->SaveAs(Form("%s.png", path));

}

void MakeComparisonPlots(const char* _dir, TCanvas *c1, hist_factory hf1, hist_factory hf2)
{
    char _file_suffix[32] = "hist_comparison";
    PlotTogether(c1, hf1.legend, hf1.hist_deltaR, hf2.hist_deltaR, Form("%s/%s_deltaR_%s", _dir, _file_suffix, hf2._objName) );
    PlotTogether(c1, hf1.legend, hf1.hist_deltaM, hf2.hist_deltaM, Form("%s/%s_deltaM_%s", _dir, _file_suffix, hf2._objName) );
    PlotTogether(c1, hf1.legend, hf1.hist_deltaPT, hf2.hist_deltaPT, Form("%s/%s_deltaPT_%s", _dir, _file_suffix, hf2._objName) );
    PlotTogether(c1, hf1.legend, hf1.hist_deltaPz, hf2.hist_deltaPz, Form("%s/%s_deltaPz_%s", _dir, _file_suffix, hf2._objName) );

    PlotTogether(c1, hf1.legend, hf1.hist_deltaM_ratio, hf2.hist_deltaM_ratio, Form("%s/%s_deltaM_ratio_%s", _dir, _file_suffix, hf2._objName) );
    PlotTogether(c1, hf1.legend, hf1.hist_deltaPT_ratio, hf2.hist_deltaPT_ratio, Form("%s/%s_deltaPT_ratio_%s", _dir, _file_suffix, hf2._objName) );
    PlotTogether(c1, hf1.legend, hf1.hist_deltaPz_ratio, hf2.hist_deltaPz_ratio, Form("%s/%s_deltaPz_ratio_%s", _dir, _file_suffix, hf2._objName) );
}

void PrintCorrelationFactors(const char* title, TH2D* h_quadratic, TH2D* h_topKinFit)
{
    printf("[INFO] %s: %9.6f %9.6f\n", title, h_quadratic->GetCorrelationFactor(), h_topKinFit->GetCorrelationFactor());
}
void MakeComparison_CorrelationFactors(hist_factory hf1, hist_factory hf2)
{
    //printf("[INFO] Correlation Factor of quadratic vs. topKinFit\n");
    PrintCorrelationFactors(Form("%9s_recPz_genPz", hf2._objName), hf1.hist_recPz_genPz, hf2.hist_recPz_genPz);
}
