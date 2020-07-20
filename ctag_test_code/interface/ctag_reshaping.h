#ifndef _CTAG_RESHAPING_H_
#define _CTAG_RESHAPING_H_

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <TFile.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TString.h>

class retrieve_scale_factor{
    public:
        retrieve_scale_factor();
        ~retrieve_scale_factor();

        void set_type_sys_uncertainty(TString);
        void set_cvsl_cvsb(double, double);
        void convert_disciminants_into_bin_numbers();
        double get_scale_factor(TString, double, double);
        void debug_mode();

    private:
        TFile *file;
        TH2D *h;
        int bin_cvsl;
        int bin_cvsb;
        double scale_factor;

        TString type_sys_uncertainty;
        double cvsl;
        double cvsb;

        bool debug_;
};

#endif
