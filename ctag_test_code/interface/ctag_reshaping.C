#include "ctag_reshaping.h"

retrieve_scale_factor::retrieve_scale_factor()
{
    file = TFile::Open("DeepCSV_ctagSF_MiniAOD94X_2017_pTincl.root");
    debug_ = false;
}

retrieve_scale_factor::~retrieve_scale_factor()
{
    file->Close();
    std::cout << "This is the end of destructor" << std::endl;
}

void retrieve_scale_factor::set_type_sys_uncertainty(TString input)
{
    type_sys_uncertainty = input;
    h = (TH2D*) file->Get(type_sys_uncertainty);
    if(debug_) printf("type:%s, ", type_sys_uncertainty.Data());
}

void retrieve_scale_factor::set_cvsl_cvsb(double cvsl_, double cvsb_)
{
    cvsl = cvsl_;
    cvsb = cvsb_;
    if(debug_) printf("cvsl = %.3f, cvsb = %.3f, ", cvsl, cvsb);
}

double retrieve_scale_factor::get_scale_factor(TString input, double cvsl, double cvsb)
{
    set_type_sys_uncertainty(input);
    set_cvsl_cvsb(cvsl, cvsb);
    convert_disciminants_into_bin_numbers();
    
    scale_factor = h->GetBinContent(bin_cvsl, bin_cvsb);
    if(debug_) printf("sf = %.2f\n", scale_factor);
    return scale_factor;
}

void retrieve_scale_factor::convert_disciminants_into_bin_numbers()
{
    int nbinx = h->GetNbinsX();
    int nbiny = h->GetNbinsY();
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();
    double ymin = h->GetYaxis()->GetXmin();
    double ymax = h->GetYaxis()->GetXmax();
    double cvsl_ = cvsl;
    double cvsb_ = cvsb;
    double step;
    int counter;

    step = (ymax - ymin) / (double) nbiny;
    counter = 0;
    while(cvsb_ >= step)
    {
        cvsb_ -= step;
        counter += 1;
    }
    bin_cvsb = counter + 1;

    step = (xmax - xmin) / (double) nbinx;
    counter = 0;
    while(cvsl_ >= step)
    {
        cvsl_ -= step;
        counter += 1;
    }
    bin_cvsl = counter + 1;
}

void retrieve_scale_factor::debug_mode()
{
    debug_ = true;
}
