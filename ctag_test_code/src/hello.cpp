#include "../interface/hello.C"

int main()
{
    printf("Hello World!\n");

    TFile *file = TFile::Open("DeepCSV_ctagSF_MiniAOD94X_2017_pTincl.root");
    TH2D *h = (TH2D*) file->Get("SFb_hist_StatUp");

    test_run(h);
    //double sf = get_scale_factor(cvsl, cvsb, h);
    //print_th2d_content(h);

    file->Close();
    return 0;
}
