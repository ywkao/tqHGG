#ifndef _HELLO_H_
#define _HELLO_H_

#include <iostream>
#include <stdio.h>
#include <TFile.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TString.h>

void test_run();
void print_th2d_content(TH2D* h);
//double get_scale_factor(double cvsl, double cvsb, TString type_sys_uncertainty);
//int convert_cvsb_into_bin_number(double cvsb, TH2D* h);
//int convert_cvsl_into_bin_number(double cvsl, TH2D* h);
//void check_region(TH2D* h);

class retrieve_scale_factor{
    public:
        retrieve_scale_factor();
        ~retrieve_scale_factor();

        void set_type_sys_uncertainty(TString);
        void set_cvsl_cvsb(double, double);
        void convert_disciminants_into_bin_numbers();

        double get_scale_factor(TString, double, double);

        void print_th2d_content(TString);

    private:
        TFile *file;
        TH2D *h;
        int bin_cvsl;
        int bin_cvsb;
        double scale_factor;

        TString type_sys_uncertainty;
        double cvsl;
        double cvsb;
};

#endif
