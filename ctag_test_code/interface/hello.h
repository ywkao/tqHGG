#ifndef _HELLO_H_
#define _HELLO_H_

#include <stdio.h>
#include <TFile.h>
#include <TH2D.h>
#include <TRandom3.h>

void test_run(TH2D* h);
double get_scale_factor(double cvsl, double cvsb, TH2D* h);
int convert_cvsb_into_bin_number(double cvsb, TH2D* h);
int convert_cvsl_into_bin_number(double cvsl, TH2D* h);
void check_region(TH2D* h);
void print_th2d_content(TH2D* h);

#endif
