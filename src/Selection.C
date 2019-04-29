#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>

void Selection(){
    char input_file[256]; sprintf(input_file, "%s", "ntuples_skimmed/ntuple_DiPhotonJetsBox_M40_80-Sherpa.root"); printf("[INFO] input_file  = %s\n", input_file);
    char output_file[256]; sprintf(output_file, "%s", "plots/test.root"); printf("[INFO] output_file = %s\n", output_file);

    TFile *fout = new TFile(output_file, "RECREATE");
    TChain *treeReader = new TChain("mytree");
    treeReader->Add(input_file);

    TCanvas *c1 = new TCanvas("c1", "c1", 700, 800);
    treeReader->Draw("num_jets");
    c1->SaveAs("plots/hist_num_jets_test.png");
}
