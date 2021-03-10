// vim: set fdm=marker:
//***************************************************************************
//
// FileName    : covarianceMatrixStudy.cc
// Purpose     : Study covariance matrix of Mjj and Mbjj.
// Description : 
// Author      : Yu-Wei Kao [ykao@cern.ch]
// OutputFile  : include/covMatrix.C
//
//***************************************************************************
#include "../include/covarianceMatrixStudy.C"
double w_boson_mass = 80.379;//GeV
double top_quark_mass = 172.5;//GeV

//int covarianceMatrixStudy(std::vector<char*> input_file, string tag){
int covarianceMatrixStudy(string input_file, string tag){
    int CONSTRAINT_NUM_JETS = get_constraint_on_num_jets(tag); //TT=4, ST=3
    printf("check CONSTRAINT_NUM_JETS = %d\n", CONSTRAINT_NUM_JETS);
    
    // ### I/O {{{
    char output_file[512]; sprintf(output_file, "result_top_reco_study/covariance_%s_hadronic.root", tag.c_str());
    TFile *fout = new TFile(output_file, "RECREATE");
    //==============================//
    //----- Read input file(s) -----//
    //==============================//
    flashggStdTreeReader treeReader;
    treeReader.InitChain("flashggNtuples/flashggStdTree");
    ////treeReader.AddSingleRootFile(input_file);
    //if(input_file.size() > 1) for(size_t i=0; i<input_file.size(); ++i) treeReader.AddSingleRootFile(input_file[i]);
    //else if(input_file.size() == 1)                                     treeReader.AddSingleRootFile(input_file[0]);
    //else                                                               {printf("[Warnning] number of input files is unexpected\n"); return 1;}
    size_t _pos1_ = input_file.find(",");
    size_t _pos2_ = input_file.length();
    if(_pos1_ != -1)
    {
        string sub_file_1 = input_file.substr(0, _pos1_);
        string sub_file_2 = input_file.substr(_pos1_+2, _pos2_);
        //printf("input_file = %s\n", sub_file_1.c_str());
        //printf("input_file = %s\n", sub_file_2.c_str());
        treeReader.AddSingleRootFile(sub_file_1);
        treeReader.AddSingleRootFile(sub_file_2);
    } else
        treeReader.AddSingleRootFile(input_file);

    treeReader.SetBranchAddresses();
    //===============================//
    //----- Prepare output file -----//
    //===============================//
    myTreeClass mytree;
    mytree.InitTree();
    mytree.MakeNewBranchAddresses();
    //------------------------------
    TNtuple *nt_mass_w   = new TNtuple("nt_mass_w", "", "mass");
    TNtuple *nt_mass_tbw = new TNtuple("nt_mass_tbw", "", "mass");
    TNtuple *nt_mass_tqh = new TNtuple("nt_mass_tqh", "", "mass");
    //TH1D *hist_deltaR_gen_reco  = new TH1D("hist_deltaR_gen_reco"    , "hist_deltaR_gen_reco; #Delta R; Entries"      , 10 , 0 , 1.);
    TH1D *hist_mass_gen_w         = new TH1D("hist_mass_gen_w"         , "hist_mass_gen_w; Mass [GeV]; Entries"         , 40 , 0 , 160);
    TH1D *hist_mass_gen_sm_top    = new TH1D("hist_mass_gen_sm_top"    , "hist_mass_gen_sm_top; Mass [GeV]; Entries"    , 40 , 0 , 320);
    TH1D *hist_mass_gen_fcnc_top  = new TH1D("hist_mass_gen_fcnc_top"  , "hist_mass_gen_fcnc_top; Mass [GeV]; Entries"  , 40 , 0 , 320);
    TH1D *hist_mass_reco_w        = new TH1D("hist_mass_reco_w"        , "hist_mass_reco_w; Mass [GeV]; Entries"        , 40 , 0 , 160);
    TH1D *hist_mass_reco_sm_top   = new TH1D("hist_mass_reco_sm_top"   , "hist_mass_reco_sm_top; Mass [GeV]; Entries"   , 40 , 0 , 320);
    TH1D *hist_mass_reco_fcnc_top = new TH1D("hist_mass_reco_fcnc_top" , "hist_mass_reco_fcnc_top; Mass [GeV]; Entries" , 40 , 0 , 320);
    //}}}
    //Declare matrix & counters{{{
    double covarianceMatrix[3][3] = {0};
    double mean[3] = {0};
    int Nevents_pass_selection = 0, Counter_before_selection_on_genParticle = 0;
    int nentries = treeReader.GetEntries(); printf("[INFO] N_entries = %d\n", nentries);
    //}}}
    
    //for(int ientry=0; ientry<1000; ientry++){
    for(int ientry=0; ientry<nentries; ientry++){
        treeReader.flashggStdTree->GetEntry(ientry);//load data
        // ### Basic Selections{{{
        //==================================================//
        //-------------   Reset Parameters   ---------------//
        //==================================================//
        mytree.Clear();
        //==================================================//
        //--------------   Basic Selectoin   ---------------//
        //==================================================//
        //require MC events pass trigger.
        if(!treeReader.EvtInfo_passTrigger) continue;// require MC events pass trigger.
        //control region
        if(treeReader.DiPhoInfo_mass<100 || treeReader.DiPhoInfo_mass>180) continue;
        //if(treeReader.DiPhoInfo_mass>120 && treeReader.DiPhoInfo_mass<130) continue;
        //require the quality of photons. (PT requirement)
        //bool pass_leadingPhotonPT = treeReader.DiPhoInfo_leadPt > treeReader.DiPhoInfo_mass / 2.;
        //bool pass_subleadingPhotonPT = treeReader.DiPhoInfo_subleadPt > treeReader.DiPhoInfo_mass / 4.;
        bool pass_leadingPhotonPT = treeReader.DiPhoInfo_leadPt > 35.;
        bool pass_subleadingPhotonPT = treeReader.DiPhoInfo_subleadPt > 25.;
        bool pass_photon_criteria_pt = pass_leadingPhotonPT && pass_subleadingPhotonPT;

        bool pass_leadingPhotonEta =  (treeReader.DiPhoInfo_leadEta < 1.4442) || (treeReader.DiPhoInfo_leadEta > 1.566 && treeReader.DiPhoInfo_leadEta < 2.4);
        bool pass_subleadingPhotonEta = (treeReader.DiPhoInfo_subleadEta < 1.4442) || (treeReader.DiPhoInfo_subleadEta > 1.566 && treeReader.DiPhoInfo_subleadEta < 2.4);
        bool pass_photon_criteria_eta = pass_leadingPhotonEta && pass_subleadingPhotonEta;

        bool pass_photon_criteria = pass_photon_criteria_pt && pass_photon_criteria_eta;
        if(!pass_photon_criteria) continue;

        bool pass_photon_IDMVA = treeReader.DiPhoInfo_leadIDMVA>-0.7 && treeReader.DiPhoInfo_subleadIDMVA>-0.7;
        if(!pass_photon_IDMVA) continue;
        //====== Hadronic Channel Criteria ======//
        if(treeReader.jets_size<3) continue;//quick skimmed
        //---  Store Diphoton Inf{{{
        TLorentzVector leading_photon, subleading_photon, diphoton;
        leading_photon.SetPtEtaPhiE(treeReader.DiPhoInfo_leadPt, treeReader.DiPhoInfo_leadEta, treeReader.DiPhoInfo_leadPhi, treeReader.DiPhoInfo_leadE);
        subleading_photon.SetPtEtaPhiE(treeReader.DiPhoInfo_subleadPt, treeReader.DiPhoInfo_subleadEta, treeReader.DiPhoInfo_subleadPhi, treeReader.DiPhoInfo_subleadE);
        diphoton = leading_photon + subleading_photon;
        //}}}
        //---  Select Electron{{{
        int num_electrons = 0;
        std::vector<TLorentzVector> Electrons;
        bool bool_AtLeastOneElectron = treeReader.ElecInfo_Size>0 ? true : false;//treeReader.ElecInfo_Size = -999 => event without diphoton candidate
        if(bool_AtLeastOneElectron){
            for(int i=0; i<treeReader.ElecInfo_Size; i++){
                if( !treeReader.ElecInfo_EGMCutBasedIDMedium->at(i) ) continue;
                if( fabs(treeReader.ElecInfo_Eta->at(i)) > 2.4 ) continue;
                if( fabs(treeReader.ElecInfo_Eta->at(i)) > 1.4442 && fabs(treeReader.ElecInfo_Eta->at(i)) < 1.566 ) continue;
                if( fabs(treeReader.ElecInfo_Pt->at(i))  < 10.  ) continue;
                //--- check deltaR(electron,photon) ---//
                TLorentzVector electron; 
                electron.SetPtEtaPhiE(treeReader.ElecInfo_Pt->at(i), treeReader.ElecInfo_Eta->at(i), treeReader.ElecInfo_Phi->at(i), treeReader.ElecInfo_Energy->at(i));
                double delta_R = electron.DeltaR(leading_photon);
                if( delta_R<0.2 ) continue;
                delta_R = electron.DeltaR(subleading_photon);
                if( delta_R<0.2 ) continue;
                num_electrons+=1;
                Electrons.push_back(electron);
            }
        }
        else{
                num_electrons=treeReader.ElecInfo_Size;
        }
        bool bool_AtLeastOneSelectedElectron = mytree.num_electrons>0 ? true : false;//for calculation of deltaR(e,j).
        //}}}
        //---  Select Muon{{{
        int num_muons = 0, num_leptons = 0;
        std::vector<TLorentzVector> Muons;
        bool bool_AtLeastOneMuon = treeReader.MuonInfo_Size>0 ? true : false;//treeReader.MuonInfo_Size = -999 => event without diphoton candidate
        if(bool_AtLeastOneMuon){
            for(int i=0; i<treeReader.MuonInfo_Size; i++){
                if( !treeReader.MuonInfo_CutBasedIdTight->at(i) ) continue;
                if( fabs(treeReader.MuonInfo_Eta->at(i)) > 2.4 ) continue;
                if( fabs(treeReader.MuonInfo_Pt->at(i))  < 5.  ) continue;
                if( treeReader.MuonInfo_PFIsoDeltaBetaCorrR04->at(i) > 0.25  ) continue;
                //--- check deltaR(muon,photon) ---//
                TLorentzVector muon; 
                muon.SetPtEtaPhiE(treeReader.MuonInfo_Pt->at(i), treeReader.MuonInfo_Eta->at(i), treeReader.MuonInfo_Phi->at(i), treeReader.MuonInfo_Energy->at(i));
                double delta_R = muon.DeltaR(leading_photon);
                if( delta_R<0.2 ) continue;
                delta_R = muon.DeltaR(subleading_photon);
                if( delta_R<0.2 ) continue;
                num_muons+=1;
                Muons.push_back(muon);
            }
        }
        else{
                num_muons=treeReader.MuonInfo_Size;
        }
        num_leptons = num_electrons + num_muons;
        bool bool_AtLeastOneSelectedMuon = mytree.num_muons>0 ? true : false;//for calculation of deltaR(mu,j).
        //}}}
        if(num_leptons > 0) continue;
        //}}}
        //### Start GenMatching for Jets{{{
        //========================================//
        //-----  Start GenMatching for Jets  -----//
        //========================================//
        //### GenMatching: find the gen particle (MC truth) for each jet (reconstructed). 
        //### pdgID: (1, 2, 3, 4, 5, 6) = (d, u, s, c, b, t)
        //### This is the simplest version. Identify the corresponding gen particle by selecting smallest deltaR(gen, jet).
        std::vector<int> index_GenParticles, GenParticles_PdgID, GenParticles_MomPdgID;
        std::vector<TLorentzVector> Jets, GenParticles;
        //--- GenMatching for each reco jet ---//
        for(int i=0; i<treeReader.jets_size; i++){
            //--- basic jet selections{{{
            if( fabs(treeReader.JetInfo_Eta->at(i)) > 2.4 ) continue;
            if( fabs(treeReader.JetInfo_Pt->at(i))  < 25.  ) continue;
            TLorentzVector jet; jet.SetPtEtaPhiE(treeReader.JetInfo_Pt->at(i), treeReader.JetInfo_Eta->at(i), treeReader.JetInfo_Phi->at(i), treeReader.JetInfo_Energy->at(i));
            double delta_R = jet.DeltaR(leading_photon);
            if( delta_R<0.4 ) continue;
            delta_R = jet.DeltaR(subleading_photon);
            if( delta_R<0.4 ) continue;
            //--------------------------------------------------//
            bool bool_passJetLeptonSeparation = true;//if no leptons selected, the jet pass the delta_R criterion automatically.
            //--- check deltaR(jet,e) ---//
            if(bool_AtLeastOneSelectedElectron){
                for(int i=0; i<mytree.num_electrons; i++){
                    delta_R = jet.DeltaR(Electrons.at(i));
                    if( delta_R<0.4 ) bool_passJetLeptonSeparation = false;
                }
            }
            if(!bool_passJetLeptonSeparation) continue;//if not pass, reject the jet.
            //--- check deltaR(jet,mu) ---//
            if(bool_AtLeastOneSelectedMuon){
                for(int i=0; i<mytree.num_muons; i++){
                    delta_R = jet.DeltaR(Muons.at(i));
                    if( delta_R<0.4 ) bool_passJetLeptonSeparation = false;
                }
            }
            if(!bool_passJetLeptonSeparation) continue;//if not pass, reject the jet.
            //--------------------------------------------------//
            //}}}
            //--- gen matching{{{
            TLorentzVector truelove;//we are looking for the right genParticle to match the jet.
            int index = -1, truelove_PdgID = -999, truelove_MomPdgID = -999; double delta_R_min = 999.;
            for(int j=0; j<treeReader.GenPartInfo_size; j++){
                if( abs(treeReader.GenPartInfo_Eta->at(j)) > 10000 ) continue;//remove incoming particle
                if( abs(treeReader.GenPartInfo_PdgID->at(j)) == 6 ) continue;//exclude top quark
                bool isAvailable = checkAvailability(j, index_GenParticles); if(!isAvailable) continue;//the truelove of previous jets shall not become another one's truelove.
                //--------------------
                TLorentzVector genParticle;
                genParticle.SetPtEtaPhiM(treeReader.GenPartInfo_Pt->at(j), treeReader.GenPartInfo_Eta->at(j), treeReader.GenPartInfo_Phi->at(j), treeReader.GenPartInfo_Mass->at(j));
                double delta_R = jet.DeltaR(genParticle);
                //select quark & require min(deltaR)
                if( abs(treeReader.GenPartInfo_PdgID->at(j)) < 7 && delta_R < 0.4 && delta_R < delta_R_min){
                    index = j;//record the matched genParticle
                    delta_R_min = delta_R;
                    truelove = genParticle;
                    truelove_PdgID = treeReader.GenPartInfo_PdgID->at(j);
                    truelove_MomPdgID = treeReader.GenPartInfo_MomPdgID->at(j);
                }
            }//end of gen loop
            //}}}
            index_GenParticles.push_back(index);
            Jets.push_back(jet);
            GenParticles.push_back(truelove);
            GenParticles_PdgID.push_back(truelove_PdgID);
            GenParticles_MomPdgID.push_back(truelove_MomPdgID);
        }//end of jet loop
        //}}}
        //check geninfo(skipped){{{
//        printf("\n");
//        printf("[GenCheck] GenPartInfo_size = %d\n", treeReader.GenPartInfo_size);
//        for(int i=0; i<treeReader.GenPartInfo_size; i++){
//            int pdgID = treeReader.GenPartInfo_PdgID->at(i);
//            bool isPromptFinalState = treeReader.GenPartInfo_isPromptFinalState->at(i);
//            bool isNeutrino = (abs(pdgID) == 12 || abs(pdgID) == 14 || abs(pdgID) == 16) && isPromptFinalState;
//            bool isChargedLepton = (abs(pdgID) == 11 || abs(pdgID) == 13 || abs(pdgID) == 15) && isPromptFinalState;
//            bool isWboson = (abs(pdgID) == 24);
//            bool isTop = (abs(pdgID) == 6);
//
//            printf("(%d) ", i);
//            printf("Status = %3d, ", treeReader.GenPartInfo_Status->at(i));
//            printf("PdgID = %3d, ", treeReader.GenPartInfo_PdgID->at(i));
//            printf("Pt = %6.2f, ", treeReader.GenPartInfo_Pt->at(i));
//            printf("Eta = %9.2f, ", treeReader.GenPartInfo_Eta->at(i));
//            printf("Phi = %6.2f, ", treeReader.GenPartInfo_Phi->at(i));
//            printf("Mass = %6.2f, ", treeReader.GenPartInfo_Mass->at(i));
//            printf("isHardProcess = %3d, ", treeReader.GenPartInfo_isHardProcess->at(i) ? 1 : 0);
//            printf("isPromptFinalState = %3d, ", treeReader.GenPartInfo_isPromptFinalState->at(i) ? 1 : 0);
//            printf("MomPdgID = %5d, ", treeReader.GenPartInfo_MomPdgID->at(i));
//            printf("MomStatus = %3d\n", treeReader.GenPartInfo_MomStatus->at(i));
//        }
        //}}}
        Counter_before_selection_on_genParticle += 1;
        //### Selection on gen-matched quarks: at least 3(4) quarks && at least 1 b quark{{{
        //--- require at least 3(4) quarks ---//
        int  count_quarks = 0;
        for(std::size_t i=0; i<GenParticles_PdgID.size(); ++i){
            //printf("(%d) ", i);
            //printf("PdgID = %d, ", GenParticles_PdgID[i]);
            //printf("MomPdgID = %d\n", GenParticles_MomPdgID[i]);
            if( abs(GenParticles_PdgID.at(i)) < 7){ count_quarks += 1; }
        }
        if(count_quarks<CONSTRAINT_NUM_JETS) continue;

        //--- reject event without bottom quark ---//
        bool has_bottom_quark;
        for(std::size_t i=0; i<GenParticles_PdgID.size(); ++i){
            if( abs(GenParticles_PdgID[i]) == 5 ){ has_bottom_quark = true; break; }
            else has_bottom_quark = false;
        }
        if(!has_bottom_quark) continue;
        //}}}
        //===============================//
        //------  Reconstruction   ------//
        //===============================//
        TLorentzVector bjet, wjets[2], tqh_jet;
        TLorentzVector bquark, wquarks[2], tqh_quark, Higgs;
        //### Check GenInfo (skipped){{{
        //printf("\n");
        //for(std::size_t i=0; i<GenParticles_PdgID.size(); ++i){
        //    printf("(%2d) ", i);
        //    printf("PdgID = %2d, ", GenParticles_PdgID[i]);
        //    printf("MomPdgID = %2d\n", GenParticles_MomPdgID[i]);
        //}
        //}}}
        //### Identified particles according to mom info && filter out events with -999 index{{{
        int index_bjet = -999;
        int index_tqh_quark = -999;
        int index_counter = 0, index_wjets[2]; std::fill(std::begin(index_wjets), std::end(index_wjets), -999);
        for(std::size_t i=0; i<GenParticles_PdgID.size(); ++i){
            if( abs(GenParticles_PdgID[i]) == 5 && abs(GenParticles_MomPdgID[i]) == 6 ){
                index_bjet = i;
                bjet = Jets[i];
                bquark = GenParticles[i];
                continue;
            }
            //### CHECK!!!
            if( abs(GenParticles_PdgID[i]) != 5 && abs(GenParticles_MomPdgID[i]) == 6 ){
                index_tqh_quark = i;
                tqh_jet = Jets[i];
                tqh_quark = GenParticles[i];
                continue;
            }
            if( abs(GenParticles_MomPdgID[i]) == 24 ){
                index_wjets[index_counter] = i;
                wjets[index_counter] = Jets[i];
                wquarks[index_counter] = GenParticles[i];
                index_counter += 1;
            }
        }
        //}}}

        vector<int> indices_reco_jets = {index_bjet, index_tqh_quark, index_wjets[0], index_wjets[1]};
        bool pass_particle_non_void_check = get_particle_non_void_check_result( tag, indices_reco_jets );
        if(!pass_particle_non_void_check) continue;

        //### gen/reco W/tops reconstructon{{{
        TLorentzVector w_candidate_reco, sm_top_candidate_reco, fcnc_top_candidate_reco;
        w_candidate_reco = wjets[0] + wjets[1];
        sm_top_candidate_reco = bjet + w_candidate_reco;
        fcnc_top_candidate_reco = tqh_jet + diphoton;
        hist_mass_reco_w -> Fill(w_candidate_reco.M());
        hist_mass_reco_sm_top -> Fill(sm_top_candidate_reco.M());
        hist_mass_reco_fcnc_top -> Fill(fcnc_top_candidate_reco.M());
        nt_mass_w -> Fill(w_candidate_reco.M());
        nt_mass_tbw -> Fill(sm_top_candidate_reco.M());
        nt_mass_tqh -> Fill(fcnc_top_candidate_reco.M());

        //NOTE: the reconstructed w boson and top quark from GEN-PARTICLES are for double check purpose only.
        //### Higgs gen info{{{
        for(int i=0; i<treeReader.GenPartInfo_size; i++){
            int pdgID = treeReader.GenPartInfo_PdgID->at(i);
            if(abs(pdgID) == 25)
                Higgs.SetPtEtaPhiM(treeReader.GenPartInfo_Pt->at(i), treeReader.GenPartInfo_Eta->at(i), treeReader.GenPartInfo_Phi->at(i), treeReader.GenPartInfo_Mass->at(i));
        }
        //}}}
        TLorentzVector w_candidate_gen, sm_top_candidate_gen, fcnc_top_candidate_gen;
        w_candidate_gen = wquarks[0] + wquarks[1];
        sm_top_candidate_gen = bquark + w_candidate_gen;
        fcnc_top_candidate_gen = tqh_quark + Higgs;
        hist_mass_gen_w -> Fill(w_candidate_gen.M());
        hist_mass_gen_sm_top -> Fill(sm_top_candidate_gen.M());
        hist_mass_gen_fcnc_top -> Fill(fcnc_top_candidate_gen.M());
        //}}}
        //### Covariant Matrix{{{
        //==================================================//
        //-------------   Covariant Matrix   ---------------//
        //==================================================//
        mean[0] += w_candidate_reco.M();
        mean[1] += sm_top_candidate_reco.M();
        mean[2] += fcnc_top_candidate_reco.M();

        double element_0 = w_candidate_reco.M() - w_boson_mass;
        double element_1 = sm_top_candidate_reco.M() - top_quark_mass;
        double element_2 = fcnc_top_candidate_reco.M() - top_quark_mass;
        covarianceMatrix[0][0] += element_0 * element_0;
        covarianceMatrix[0][1] += element_0 * element_1;
        covarianceMatrix[0][2] += element_0 * element_2;
        covarianceMatrix[1][0] += element_1 * element_0;
        covarianceMatrix[1][1] += element_1 * element_1;
        covarianceMatrix[1][2] += element_1 * element_2;
        covarianceMatrix[2][0] += element_2 * element_0;
        covarianceMatrix[2][1] += element_2 * element_1;
        covarianceMatrix[2][2] += element_2 * element_2;
        //}}}
        //==================================================//
        //-------------   Event Counting     ---------------//
        //==================================================//
        Nevents_pass_selection += 1;
        mytree.Fill();
    }// End of event loop.

    //### Report{{{
    //==================================================//
    //---------------------  Report  -------------------//
    //==================================================//
    double eff = 100. * (double) Nevents_pass_selection / (double) Counter_before_selection_on_genParticle;
    printf("[INFO] Nevents_pass_selection = %d\n", Nevents_pass_selection);
    printf("[INFO] before/after gen selection: %d/%d (%.3f%%)\n", Nevents_pass_selection, Counter_before_selection_on_genParticle, eff);

    //TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    //MakePlots(c1 , hist_mass_gen_w         , "" , Form("cov_hist_mass_%s_gen_w.png", tag.c_str()) );
    //MakePlots(c1 , hist_mass_gen_sm_top    , "" , Form("cov_hist_mass_%s_gen_sm_top.png", tag.c_str()) );
    //MakePlots(c1 , hist_mass_gen_fcnc_top  , "" , Form("cov_hist_mass_%s_gen_fcnc_top.png", tag.c_str()) );
    //MakePlots(c1 , hist_mass_reco_w        , "" , Form("cov_hist_mass_%s_reco_w.png", tag.c_str()) );
    //MakePlots(c1 , hist_mass_reco_sm_top   , "" , Form("cov_hist_mass_%s_reco_sm_top.png", tag.c_str()) );
    //MakePlots(c1 , hist_mass_reco_fcnc_top , "" , Form("cov_hist_mass_%s_reco_fcnc_top.png", tag.c_str()) );

    //}}}
    //### Covariant Matrix{{{
    // calculation{{{
    mean[0] /= Nevents_pass_selection;
    mean[1] /= Nevents_pass_selection;
    mean[2] /= Nevents_pass_selection;
    covarianceMatrix[0][0] /= Nevents_pass_selection;
    covarianceMatrix[0][1] /= Nevents_pass_selection;
    covarianceMatrix[0][2] /= Nevents_pass_selection;
    covarianceMatrix[1][0] /= Nevents_pass_selection;
    covarianceMatrix[1][1] /= Nevents_pass_selection;
    covarianceMatrix[1][2] /= Nevents_pass_selection;
    covarianceMatrix[2][0] /= Nevents_pass_selection;
    covarianceMatrix[2][1] /= Nevents_pass_selection;
    covarianceMatrix[2][2] /= Nevents_pass_selection;
    //}}}
    // printf{{{
    printf("mean = %6.2f, sigma_prime_jj  = %6.2f\n", mean[0], sqrt(covarianceMatrix[0][0]));
    printf("mean = %6.2f, sigma_prime_bjj = %6.2f\n", mean[1], sqrt(covarianceMatrix[1][1]));
    printf("mean = %6.2f, sigma_prime_qgg = %6.2f\n", mean[2], sqrt(covarianceMatrix[2][2]));
    printf("---------------\n");
    printf("matrix(0,0) = %6.2f; matrix(0,1) = %6.2f;\n", covarianceMatrix[0][0], 0.);
    printf("matrix(1,0) = %6.2f; matrix(1,1) = %6.2f;\n", 0., covarianceMatrix[1][1]);
    printf("---------------\n");
    printf("matrix(0,0) = %6.2f; matrix(0,1) = %6.2f;\n", covarianceMatrix[0][0], covarianceMatrix[0][1]);
    printf("matrix(1,0) = %6.2f; matrix(1,1) = %6.2f;\n", covarianceMatrix[1][0], covarianceMatrix[1][1]);
    printf("---------------\n");
    printf("matrix(0,0) = %6.2f; matrix(0,1) = %6.2f; matrix(0,2) = %6.2f;\n", covarianceMatrix[0][0], covarianceMatrix[0][1], covarianceMatrix[0][2]);
    printf("matrix(1,0) = %6.2f; matrix(1,1) = %6.2f; matrix(1,2) = %6.2f;\n", covarianceMatrix[1][0], covarianceMatrix[1][1], covarianceMatrix[1][2]);
    printf("matrix(2,0) = %6.2f; matrix(2,1) = %6.2f; matrix(2,2) = %6.2f;\n", covarianceMatrix[2][0], covarianceMatrix[2][1], covarianceMatrix[2][2]);
    printf("\n");
    //}}}
    //}}}


    char* message = new char[1024];
    ofstream myfile;
    // output file{{{
    myfile.open(Form("include/covMatrix_%s.C", tag.c_str()));
    myfile << "//### AUTOMATICALLY GENERATED BY src/covarianceMatrixStudy.cpp ###//\n";
    sprintf( message, Form("//### DATASET: %s\n", input_file.c_str()) ); myfile << message;
    myfile << "#include <TMatrixD.h>\n";
    myfile << "\n\n";
    myfile << "double Chi2_calculator_simple(double w_mass, double t_mass){\n";
    myfile << "    TVectorD vec_mass(2);\n";
    myfile << "    vec_mass(0) = w_mass - w_boson_mass;\n";
    myfile << "    vec_mass(1) = t_mass - top_quark_mass;\n";
    myfile << "    TMatrixD matrix(2,2);\n";
    sprintf( message, Form("    matrix(0,0) = %6.2f; matrix(0,1) = %6.2f;\n", covarianceMatrix[0][0], 0.) ); myfile << message;
    sprintf( message, Form("    matrix(1,0) = %6.2f; matrix(1,1) = %6.2f;\n", 0., covarianceMatrix[1][1]) ); myfile << message;
    //myfile << "    matrix(0,0) = 345.84; matrix(0,1) =   0.00;\n";
    //myfile << "    matrix(1,0) =   0.00; matrix(1,1) = 961.79;\n";
    myfile << "    return matrix.Invert()*vec_mass*vec_mass;\n";
    myfile << "}\n";
    myfile << "\n";
    myfile << "double Chi2_calculator_modified(double w_mass, double t_mass){\n";
    myfile << "    TVectorD vec_mass(2);\n";
    myfile << "    vec_mass(0) = w_mass - w_boson_mass;\n";
    myfile << "    vec_mass(1) = t_mass - top_quark_mass;\n";
    myfile << "    TMatrixD matrix(2,2);\n";
    sprintf( message, Form("    matrix(0,0) = %6.2f; matrix(0,1) = %6.2f;\n", covarianceMatrix[0][0], covarianceMatrix[0][1]) ); myfile << message;
    sprintf( message, Form("    matrix(1,0) = %6.2f; matrix(1,1) = %6.2f;\n", covarianceMatrix[1][0], covarianceMatrix[1][1]) ); myfile << message;
    //myfile << "    matrix(0,0) = 345.84; matrix(0,1) = 378.21;\n";
    //myfile << "    matrix(1,0) = 378.21; matrix(1,1) = 961.79;\n";
    myfile << "    return matrix.Invert()*vec_mass*vec_mass;\n";
    myfile << "}\n";
    myfile << "\n";
    myfile << "double Chi2_calculator_improved(double w_mass, double t_mass, double fcnc_top_mass){\n";
    myfile << "    TVectorD vec_mass(3);\n";
    myfile << "    vec_mass(0) = w_mass - w_boson_mass;\n";
    myfile << "    vec_mass(1) = t_mass - top_quark_mass;\n";
    myfile << "    vec_mass(2) = fcnc_top_mass - top_quark_mass;\n";
    myfile << "    TMatrixD matrix(3,3);\n";
    sprintf( message, Form("    matrix(0,0) = %6.2f; matrix(0,1) = %6.2f; matrix(0,2) = %6.2f;\n", covarianceMatrix[0][0], covarianceMatrix[0][1], covarianceMatrix[0][2]) ); myfile << message;
    sprintf( message, Form("    matrix(1,0) = %6.2f; matrix(1,1) = %6.2f; matrix(1,2) = %6.2f;\n", covarianceMatrix[1][0], covarianceMatrix[1][1], covarianceMatrix[1][2]) ); myfile << message;
    sprintf( message, Form("    matrix(2,0) = %6.2f; matrix(2,1) = %6.2f; matrix(2,2) = %6.2f;\n", covarianceMatrix[2][0], covarianceMatrix[2][1], covarianceMatrix[2][2]) ); myfile << message;
    //myfile << "    matrix(0,0) = 345.84; matrix(0,1) = 378.21; matrix(0,2) =   3.26;\n";
    //myfile << "    matrix(1,0) = 378.21; matrix(1,1) = 961.79; matrix(1,2) =   3.78;\n";
    //myfile << "    matrix(2,0) =   3.26; matrix(2,1) =   3.78; matrix(2,2) = 272.53;\n";
    myfile << "    return matrix.Invert()*vec_mass*vec_mass;\n";
    myfile << "}\n";
    myfile.close();
    //}}}
    // another output file{{{
    myfile.open(Form("include/covMatrix_%s.txt", tag.c_str()));
    myfile << "//### AUTOMATICALLY GENERATED BY src/covarianceMatrixStudy.cpp ###//\n";
    sprintf( message, Form("//### DATASET: %s\n", input_file.c_str()) ); myfile << message;
    sprintf( message, Form("%6.2f, %6.2f, %6.2f\n", covarianceMatrix[0][0], covarianceMatrix[0][1], covarianceMatrix[0][2]) ); myfile << message;
    sprintf( message, Form("%6.2f, %6.2f, %6.2f\n", covarianceMatrix[1][0], covarianceMatrix[1][1], covarianceMatrix[1][2]) ); myfile << message;
    sprintf( message, Form("%6.2f, %6.2f, %6.2f\n", covarianceMatrix[2][0], covarianceMatrix[2][1], covarianceMatrix[2][2]) ); myfile << message;
    myfile.close();
    //}}}

    fout->Write();
    fout->Close();

    delete message;

    return 0;
}

int main(int argc, char *argv[]){
    vector<string> tags = {"st_hct", "tt_hct", "st_hut", "tt_hut"};

    char dir[512] = "/wk_cms2/youying/public/tH_FCNC/Era2017_RR-31Mar2018_v2";
    vector<string> vec_rootfiles =
    {
        Form("%s/ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root"     , dir),
        Form("%s/TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root, %s/TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8.root" , dir, dir),
        Form("%s/ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root"     , dir),
        Form("%s/TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root, %s/TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8.root" , dir, dir)
    };

    //vector< vector<char*> > vec_rootfiles = {
    //    {Form("%s/ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root"     , dir)},
    //    {Form("%s/TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8.root" , dir), Form("%s/TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8.root" , dir)},
    //    {Form("%s/ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root"     , dir)},
    //    {Form("%s/TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8.root" , dir), Form("%s/TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8.root" , dir)}
    //};

    for(size_t i=0; i<vec_rootfiles.size(); ++i)
    {
        printf("check tag: %s\n", tags[i].c_str());
        printf("check rootfiles = %s\n", vec_rootfiles[i].c_str());
        covarianceMatrixStudy( vec_rootfiles[i], tags[i] );
    }

    return 0;
}
