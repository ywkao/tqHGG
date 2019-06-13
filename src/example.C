{
    char directory[512];
    //=== NOTE: replace with your directory ===//
    sprintf(directory, "%s", "/wk_cms2/youying/public/RunII2016_8028_flashgg_ntuples/DoubleEG_sethzenz-LegacyReReco-07Aug2017-2_6_1-2_6_1-v1-Run2016B-07Aug17_ver2-v2/");

    TChain *treeReader;
    //=== NOTE: replace with your treeName ===//
    treeReader = new TChain("flashggNtuples/flashggStdTree");
    treeReader->Add(Form("%s/*.root", directory));

    int jets_size;
    treeReader->SetBranchAddress("jets_size", &jets_size);

    int nentries = treeReader->GetEntries(); printf("[INFO] N_entries = %d\n", nentries);

    for(int ientry=0; ientry<nentries; ientry++){
        if((ientry+1)%1000==0 || (ientry+1)==nentries) printf("ientry = %d\r", ientry);
        treeReader->GetEntry(ientry);
        if(jets_size > 0) printf("jets_size = %d\r", jets_size);
    }
}
