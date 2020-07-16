#include "hello.h"

void test_run(TH2D* h)
{
    TRandom3 random(1234);

    for(int i=0; i<100; ++i)
    {
        double cvsl = random.Rndm();
        double cvsb = random.Rndm();
        double sf = get_scale_factor(cvsl, cvsb, h);
        printf("sf = %.3f\n", sf);
    }

}

double get_scale_factor(double cvsl, double cvsb, TH2D* h)
{
    int bin_cvsl = convert_cvsl_into_bin_number(cvsl, h);
    int bin_cvsb = convert_cvsb_into_bin_number(cvsb, h);
    double scale_factor = h->GetBinContent(bin_cvsl, bin_cvsb);

    printf("(cvsl, cvsb) = (%.2f, %.2f), bins = (%d, %d), ", cvsl, cvsb, bin_cvsl, bin_cvsb);

    return scale_factor;
}

int convert_cvsb_into_bin_number(double cvsb, TH2D* h)
{
    int nbiny = h->GetNbinsY();
    double ymin = h->GetYaxis()->GetXmin();
    double ymax = h->GetYaxis()->GetXmax();

    double step = (ymax - ymin) / (double) nbiny;

    int counter = 0;
    while(cvsb >= step)
    {
        cvsb -= step;
        counter += 1;
    }

    return counter + 1;
}

int convert_cvsl_into_bin_number(double cvsl, TH2D* h)
{
    int nbinx = h->GetNbinsX();
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();

    double step = (xmax - xmin) / (double) nbinx;

    int counter = 0;
    while(cvsl >= step)
    {
        cvsl -= step;
        counter += 1;
    }

    return counter + 1;
}

void check_region(TH2D* h)
{
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();
    double ymin = h->GetYaxis()->GetXmin();
    double ymax = h->GetYaxis()->GetXmax();

    printf("[xmin, xmax] = [%f, %f]; ", xmin, xmax);
    printf("[ymin, ymax] = [%f, %f]\n", ymin, ymax);
}

void print_th2d_content(TH2D* h)
{
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
