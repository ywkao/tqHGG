#ifndef HIST_COMPONENTS_H
#define HIST_COMPONENTS_H

#include <stdio.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

class component_plots
{
    public:
        component_plots();
        virtual ~component_plots();

        void Init(const char* dirName, const char* histName, const char* ext_01, const char* ext_02, int nbins, double lower_bound, double upper_bound);
        void Fill_hist(bool canFindSolution, double value);
        void Draw_hist(const char* label1, const char* label2);
        void Report_hist(const char* label1, const char* label2);

    private:
        char histName_[128];
        char dirName_[64];

        TCanvas *c1;
        TLegend *legend;

        TH1D* hist;
        TH1D* hist_positive;
        TH1D* hist_negative;
};
#endif
