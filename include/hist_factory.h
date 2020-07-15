#ifndef HIST_FACTORY_H
#define HIST_FACTORY_H

#include <stdio.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLorentzVector.h>

class hist_factory
{
    public:
        hist_factory();
        hist_factory(const char* , const char* , const char* , double);
        ~hist_factory();

        void Fill_hist(TLorentzVector obj_reco, TLorentzVector obj_gen);
        int  Get_histEntries();
        void Draw_individual_hist(TCanvas *c1, TH1D *hist);
        void Draw_individual_hist(TCanvas *c1, TH2D *hist, const char* option);
        void Draw_all_hist(TCanvas *c1);
        void Draw_gen_reco(TCanvas *c1, TH1D *h1, TH1D *h2);

        friend double GetMaxScope(TH1D* h1, TH1D* h2);
        friend void PlotTogether(TCanvas* c1, TLegend *legend, TH1D *h1, TH1D *h2, const char* path);
        friend void MakeComparisonPlots(const char* _dir, TCanvas *c1, hist_factory hf1, hist_factory hf2);
        friend void PrintCorrelationFactors(const char* title, TH2D* h_quadratic, TH2D* h_topKinFit);
        friend void MakeComparison_CorrelationFactors(hist_factory hf1, hist_factory hf2);


    private:
        char _objName[128];
        char _dirName[64];

        TLegend *legend;

        // basik quantities
        TH1D* hist_pz;
        TH1D* hist_pt;
        TH1D* hist_eta;
        TH1D* hist_phi;
        TH1D* hist_mass;
        TH1D* hist_gen_mass;

        // check reco-gen quantities
        TH1D* hist_deltaR;
        TH1D* hist_deltaM;
        TH1D* hist_deltaPT;
        TH1D* hist_deltaPz;

        TH1D* hist_deltaM_ratio;
        TH1D* hist_deltaPT_ratio;
        TH1D* hist_deltaPz_ratio;

        TH2D* hist_recPz_genPz; // check neutrino Pz
        TH2D* hist_deltaM_genM; // check W and top mass resolution
        TH2D* hist_deltaM_genM_mean;
        TH2D* hist_deltaM_genM_rms;

};
#endif
