#ifndef __PUTOOLS_H__
#define __PUTOOLS_H__
#include <TH1D.h>

void Rebin(TH1D *h_mc_input, TH1D* &h_mc){
    for(int i=0; i<75; ++i){//rebin from 100 to 75
        double content = h_mc_input->GetBinContent(i+1);
        double error = h_mc_input->GetBinError(i+1);
        h_mc->SetBinContent(i+1, content);
        h_mc->SetBinError(i+1, error);
    }
}

double SumContents(TH1D* &h){
    double Ntotal = 0.;
    //for(int i=0; i<99; ++i){
    //    Ntotal += h->GetBinContent(i+1);
    //}
    Ntotal = h->Integral();
    return Ntotal;
}

void Normalization(TH1D* &h, double Ntotal){
    int num = h->GetNbinsX();
    //printf("num = %d\n", num);
    for(int i=0; i<num; ++i){
        double content = h->GetBinContent(i+1);
        double error = h->GetBinError(i+1);
        h->SetBinContent(i+1, content/Ntotal);
        h->SetBinError(i+1, error/Ntotal);
    }
}

#endif
