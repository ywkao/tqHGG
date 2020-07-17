#include "hello.h"

void test_run()
{
    retrieve_scale_factor sf;
    TRandom3 random(1234);

    TString type_flavour[3] = {"b", "c", "l"}; 
    TString type_uncertainty[14] = {"Stat", "EleIDSF", "LHEScaleWeight_muF", "LHEScaleWeight_muR", "MuIDSF", "PSWeightFSR", "PSWeightISR", "PUWeight", "XSec_DYJets", "XSec_ST", "XSec_WJets", "XSec_ttbar", "jer", "jesTotal"};

    for(int i=0; i<3; ++i)
    {
        for(int j=0; j<14; ++j)
        {
            if(!(i==0 and j==5)) continue;
    //------------------------------------------------------------------------------------------//
            TString name = "SF" + type_flavour[i] + "_hist" + "_" + type_uncertainty[j] + "Up";
            std::cout << name << std::endl;
            sf.print_th2d_content(name);

            for(int i=0; i<100; ++i)
            {
                double cvsl = random.Rndm();
                double cvsb = random.Rndm();
                double value = sf.get_scale_factor(name, cvsl, cvsb);
                printf("sf = %.3f\n", value);
            }
    //------------------------------------------------------------------------------------------//
        }
    }
    std::cout << "This is the end of test run" << std::endl;
}

//------------------------------------------ class functions ------------------------------------------//
retrieve_scale_factor::retrieve_scale_factor()
{
    file = TFile::Open("DeepCSV_ctagSF_MiniAOD94X_2017_pTincl.root");
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
}

void retrieve_scale_factor::set_cvsl_cvsb(double _cvsl, double _cvsb)
{
    cvsl = _cvsl;
    cvsb = _cvsb;
}

double retrieve_scale_factor::get_scale_factor(TString input, double cvsl, double cvsb)
{
    set_type_sys_uncertainty(input);
    set_cvsl_cvsb(cvsl, cvsb);
    convert_disciminants_into_bin_numbers();
    
    printf("(cvsl, cvsb) = (%.2f, %.2f), bins = (%d, %d), ", cvsl, cvsb, bin_cvsl, bin_cvsb);

    scale_factor = h->GetBinContent(bin_cvsl, bin_cvsb);
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

void retrieve_scale_factor::print_th2d_content(TString input)
{
    set_type_sys_uncertainty(input);
    int nbinx = h->GetNbinsX();
    int nbiny = h->GetNbinsY();

    for(int j=0; j < nbiny; ++j)
    {
        for(int i=0; i < nbinx; ++i)
        {
            double sf = h->GetBinContent(i+1, nbiny - j);
            printf("%3.1f ", sf);
        }
        printf("\n");
    }
}

// legacy code{{{
//double get_scale_factor(double cvsl, double cvsb, TString type_sys_uncertainty)
//{
//    TFile *file = TFile::Open("DeepCSV_ctagSF_MiniAOD94X_2017_pTincl.root");
//    TH2D *h = (TH2D*) file->Get(type_sys_uncertainty);
//
//    int bin_cvsl = convert_cvsl_into_bin_number(cvsl, h);
//    int bin_cvsb = convert_cvsb_into_bin_number(cvsb, h);
//    double scale_factor = h->GetBinContent(bin_cvsl, bin_cvsb);
//
//    printf("(cvsl, cvsb) = (%.2f, %.2f), bins = (%d, %d), ", cvsl, cvsb, bin_cvsl, bin_cvsb);
//
//    file->Close();
//    return scale_factor;
//}
//
//int convert_cvsb_into_bin_number(double cvsb, TH2D* h)
//{
//    int nbiny = h->GetNbinsY();
//    double ymin = h->GetYaxis()->GetXmin();
//    double ymax = h->GetYaxis()->GetXmax();
//
//    double step = (ymax - ymin) / (double) nbiny;
//
//    int counter = 0;
//    while(cvsb >= step)
//    {
//        cvsb -= step;
//        counter += 1;
//    }
//
//    return counter + 1;
//}
//
//int convert_cvsl_into_bin_number(double cvsl, TH2D* h)
//{
//    int nbinx = h->GetNbinsX();
//    double xmin = h->GetXaxis()->GetXmin();
//    double xmax = h->GetXaxis()->GetXmax();
//
//    double step = (xmax - xmin) / (double) nbinx;
//
//    int counter = 0;
//    while(cvsl >= step)
//    {
//        cvsl -= step;
//        counter += 1;
//    }
//
//    return counter + 1;
//}
//
//void check_region(TH2D* h)
//{
//    double xmin = h->GetXaxis()->GetXmin();
//    double xmax = h->GetXaxis()->GetXmax();
//    double ymin = h->GetYaxis()->GetXmin();
//    double ymax = h->GetYaxis()->GetXmax();
//
//    printf("[xmin, xmax] = [%f, %f]; ", xmin, xmax);
//    printf("[ymin, ymax] = [%f, %f]\n", ymin, ymax);
//}
//}}}
