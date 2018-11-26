#include <stdio.h>
#include <TFile.h>
#include <TTree.h>
#include <vector>

//#define PATH "/home/xiaokao/Desktop/play/analysis/tthTest.root"

int main(){
    TFile *fin = TFile::Open("/home/xiaokao/Desktop/play/analysis/tthTest.root");
    //TTree *flashggStdTree = (TTree*)fin->Get("flashggNtuple/flashggStdTree");

    Int_t jets_size=0;
    //flashggStdTree->SetBranchAddress("jets_size", &jets_size);
    //flashggStdTree->GetEntry(100);
    printf("HelloWorld\n");
    printf("Jet_size = %d\n", jets_size);
    return 1;
}
