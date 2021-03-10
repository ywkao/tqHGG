#include "generalChiSquareStudy.h"
#include "covMatrix.C"
//chi-2 method related{{{
//chi-2 with 3x3 covariance matrix{{{
std::vector<int> get_bjjq_indices_min_chi2_3x3(std::vector<TLorentzVector> Jets, std::vector<int> indices_bjet, TLorentzVector diphoton)
{
    vector<double> vec_chi2;
    vector<vector<int>> vec_indices_jet;
    for(std::size_t i=0; i!=indices_bjet.size(); ++i){
        double chi2_min = 99999;
        std::vector<int> index_jet_chi2_simple = get_bjj_indices_min_chi2_3x3(Jets, indices_bjet[i], diphoton, chi2_min); // indices of b j j q
        vec_chi2.push_back(chi2_min);
        vec_indices_jet.push_back(index_jet_chi2_simple);
        //printf("[check-imp] [%d] (b, j, j, q) = ", i);
        //printf("(%d, ", index_jet_chi2_simple[0]);
        //printf("%d, " , index_jet_chi2_simple[1]);
        //printf("%d, " , index_jet_chi2_simple[2]);
        //printf("%d), ", index_jet_chi2_simple[3]);
        //printf("chi2 = %7.3f\n", chi2_min);
    }

    int min_index =  std::min_element(vec_chi2.begin(), vec_chi2.end()) - vec_chi2.begin();
    double min    = *std::min_element(vec_chi2.begin(), vec_chi2.end());
    //printf("[curiosity] begin, end, diff = %d, %d, %d\n", vec_chi2.begin(), vec_chi2.end(), vec_chi2.end() - vec_chi2.begin());
    //printf("[curiosity] *begin, *end     = %f, %f\n", *vec_chi2.begin(), *(vec_chi2.end()-1));
    //printf("[check-imp] min: [%d] (b, j, j, q) = ", min_index);
    //printf("(%d, ", vec_indices_jet[min_index][0]);
    //printf("%d, " , vec_indices_jet[min_index][1]);
    //printf("%d, " , vec_indices_jet[min_index][2]);
    //printf("%d), ", vec_indices_jet[min_index][3]);
    //printf("chi2 = %7.3f\n\n", min);
    return vec_indices_jet[min_index];
}


std::vector<int> get_bjj_indices_min_chi2_3x3(std::vector<TLorentzVector> Jets, int index_bjet, TLorentzVector diphoton, double &chi2_min)
{
    // 1. pick up 3 jets
    // 2. q-jj three combinations 
    int index_tqh_qjet;
    std::vector<int> index_wjets(2);
    TLorentzVector wjets[2], tqh_qjet;
    TLorentzVector w_boson, sm_top, fcnc_top; 
    TLorentzVector bjet = Jets[index_bjet];
    std::size_t num_jets = Jets.size();
    for(std::size_t i=0; i!=num_jets; ++i){
        if(i==index_bjet) continue;//bypass bjet
        for(std::size_t j=i+1; j!=num_jets; ++j){
            if(j==index_bjet) continue;//bypass bjet
            for(std::size_t k=j+1; k!=num_jets; ++k){
                if(k==index_bjet) continue;//bypass bjet
                //--- combinations{{{
                TLorentzVector jets_chosen[3];
                jets_chosen[0] = Jets[i];
                jets_chosen[1] = Jets[j];
                jets_chosen[2] = Jets[k];

                TLorentzVector tqh_q_chosen[3];
                tqh_q_chosen[0] = jets_chosen[0];
                tqh_q_chosen[1] = jets_chosen[1];
                tqh_q_chosen[2] = jets_chosen[2];

                TLorentzVector w_candidate[3];
                w_candidate[0] = jets_chosen[1] + jets_chosen[2];
                w_candidate[1] = jets_chosen[0] + jets_chosen[2];
                w_candidate[2] = jets_chosen[0] + jets_chosen[1];

                TLorentzVector top_candidate[3];
                TLorentzVector fcnc_top_candidate[3];
                //}}}
                //--- calculation{{{
                std::vector<double> chi2;
                double w_mass[3], t_mass[3], fcnc_top_mass[3];
                for(int x = 0; x<3; ++x){
                    w_mass[x] =  w_candidate[x].M();

                    top_candidate[x] = w_candidate[x] + bjet;
                    t_mass[x] =  top_candidate[x].M();

                    fcnc_top_candidate[x] = tqh_q_chosen[x] + diphoton;
                    fcnc_top_mass[x] = fcnc_top_candidate[x].M();

                    chi2.push_back( Chi2_calculator_improved(w_mass[x], t_mass[x], fcnc_top_mass[x]) );
                }
                //}}}
                //--- sorting{{{
                int smallest_chi2_index = std::min_element(chi2.begin(),chi2.end()) - chi2.begin();
                double smallest_chi2 = *std::min_element(chi2.begin(),chi2.end());

                if(smallest_chi2 < chi2_min){
                    if(smallest_chi2_index == 0){
                        tqh_qjet = tqh_q_chosen[0];
                        wjets[0] = jets_chosen[1];
                        wjets[1] = jets_chosen[2];
                        index_tqh_qjet = i;
                        index_wjets[0] = j;
                        index_wjets[1] = k;
                    } else if(smallest_chi2_index == 1){
                        tqh_qjet = tqh_q_chosen[1];
                        wjets[0] = jets_chosen[0];
                        wjets[1] = jets_chosen[2];
                        index_tqh_qjet = j;
                        index_wjets[0] = i;
                        index_wjets[1] = k;
                    } else{
                        tqh_qjet = tqh_q_chosen[2];
                        wjets[0] = jets_chosen[0];
                        wjets[1] = jets_chosen[1];
                        index_tqh_qjet = k;
                        index_wjets[0] = i;
                        index_wjets[1] = j;
                    }

                    chi2_min = smallest_chi2;
                }
                //}}}
            }
        }
    }//end of looping jets
    w_boson  = wjets[0] + wjets[1];
    sm_top   = w_boson + bjet;
    fcnc_top = tqh_qjet + diphoton;

    std::vector<int> indices = {index_bjet, index_wjets[0], index_wjets[1], index_tqh_qjet};
    return indices;
}
//}}}
// chi-2 with 2x2 covariance matrix{{{
vector<int> get_bjjq_indices_min_chi2(std::vector<TLorentzVector> Jets, std::vector<int> indices_bjet, TLorentzVector diphoton, bool is_chi2_modified)
{
    std::vector<int> indices = get_bjj_indices_min_chi2(Jets, indices_bjet, is_chi2_modified);
    int index_q = get_q_index_min_chi2(Jets, indices, diphoton);
    indices.push_back(index_q);
    return indices;
}


int get_q_index_min_chi2(std::vector<TLorentzVector> Jets, std::vector<int> indices_bjj, TLorentzVector diphoton)
{
    std::vector<int> indices;
    std::vector<double> top_fcnh_chi2;
    for(std::size_t i=0; i!=Jets.size(); ++i){
        if(i==indices_bjj[0] || i==indices_bjj[1] || i==indices_bjj[2]) continue; //skip the selected jets for bjj
        TLorentzVector top_fcnh_tmp = diphoton + Jets[i];
        double chi2 = (top_fcnh_tmp.M() - top_quark_mass) * (top_fcnh_tmp.M() - top_quark_mass);
        indices.push_back(i);
        top_fcnh_chi2.push_back(chi2);
        //printf("[check-ywk] q = ");
        //printf("%d, " , i);
        //printf("chi2 = %7.3f\n", chi2);
    }

    int min_index =  std::min_element(top_fcnh_chi2.begin(), top_fcnh_chi2.end()) - top_fcnh_chi2.begin();
    double min    = *std::min_element(top_fcnh_chi2.begin(), top_fcnh_chi2.end());
    //printf("[check-ywk] min: q = ");
    //printf("%d, " , indices[min_index]);
    //printf("chi2 = %7.3f\n", top_fcnh_chi2[min_index]);

    return indices[min_index];
}


vector<int> get_bjj_indices_min_chi2(std::vector<TLorentzVector> Jets, std::vector<int> indices_bjet, bool is_chi2_modified)
{
    vector<double> vec_chi2;
    vector<vector<int>> vec_indices_jet;
    for(std::size_t i=0; i!=indices_bjet.size(); ++i){
        double chi2_min = 99999;
        std::vector<int> index_jet_chi2_simple = get_indices_chi2(Jets, indices_bjet[i], chi2_min, is_chi2_modified); // indices of b j j
        vec_chi2.push_back(chi2_min);
        vec_indices_jet.push_back(index_jet_chi2_simple);
        //printf("[check-ywk] [%d] (b, j, j) = ", i);
        //printf("(%d, ", index_jet_chi2_simple[0]);
        //printf("%d, " , index_jet_chi2_simple[1]);
        //printf("%d), ", index_jet_chi2_simple[2]);
        //printf("chi2 = %7.3f\n", chi2_min);
    }

    int min_index =  std::min_element(vec_chi2.begin(), vec_chi2.end()) - vec_chi2.begin();
    double min    = *std::min_element(vec_chi2.begin(), vec_chi2.end());
    //printf("[curiosity] begin, end, diff = %d, %d, %d\n", vec_chi2.begin(), vec_chi2.end(), vec_chi2.end() - vec_chi2.begin());
    //printf("[curiosity] *begin, *end     = %f, %f\n", *vec_chi2.begin(), *(vec_chi2.end()-1));
    //printf("[check-ywk] min: [%d] (b, j, j) = ", min_index);
    //printf("(%d, ", vec_indices_jet[min_index][0]);
    //printf("%d, " , vec_indices_jet[min_index][1]);
    //printf("%d), ", vec_indices_jet[min_index][2]);
    //printf("chi2 = %7.3f\n\n", min);
    return vec_indices_jet[min_index];
}


std::vector<int> get_indices_chi2(std::vector<TLorentzVector> Jets, int index_bjet, double &chi2_min, bool is_chi2_modified)
{
    std::vector<int> indices_selected_jet(3, -999);
    indices_selected_jet[0] = index_bjet;
    TLorentzVector bjet = Jets[index_bjet];

    std::size_t num_jets = Jets.size();
    for(std::size_t i=0; i<num_jets; ++i){
        if(i==index_bjet) continue;//bypass bjet
        for(std::size_t j=i+1; j<num_jets; ++j){
            if(j==index_bjet) continue;//bypass bjet
            TLorentzVector w_candidate = Jets[i] + Jets[j];
            double w_mass = w_candidate.M();
            TLorentzVector top_candidate = w_candidate + bjet;
            double t_mass = top_candidate.M();
            double chi2 = !is_chi2_modified ? Chi2_calculator_simple(w_mass, t_mass) : Chi2_calculator_modified(w_mass, t_mass);
            if(chi2 < chi2_min){
                indices_selected_jet[1] = i;
                indices_selected_jet[2] = j;
                chi2_min = chi2;
            }
        }
    }//end of looping jets

    return indices_selected_jet; //indices of b j j
}
//}}}
//others{{{
bool checkAvailability(int index, std::vector<int> ID_IsChosen)
{
    bool result = true;//if pass the following for loop, the genParticle of the index is available.
    for(std::size_t i=0; i<ID_IsChosen.size(); ++i){
        if(index == ID_IsChosen[i]){ result = false; break; } 
    }
    return result;
}


TLorentzVector GetBestM1(double &M1, int num_jets, int index_bjet, std::vector<int> index_jet, TLorentzVector diphoton, std::vector<TLorentzVector> Jets, int &index_q, TLorentzVector &jet_q)
{
    //record all the possible combinations
    std::vector<int> indices;
    std::vector<double> top_fcnh_chi2;
    std::vector<TLorentzVector> top_fcnh_candidates;
    for(int i=0; i<num_jets; ++i){
        if(i==index_bjet || i==index_jet[0] || i==index_jet[1]) continue;//skip the jets for bjj
        TLorentzVector top_fcnh_tmp = diphoton + Jets[i];
        double chi2 = (top_fcnh_tmp.M() - top_quark_mass) * (top_fcnh_tmp.M() - top_quark_mass);
        indices.push_back(i);
        top_fcnh_chi2.push_back(chi2);
        top_fcnh_candidates.push_back(top_fcnh_tmp);
    }
    //choose the candidate with Mass closest to top mass
    int index_M1_min_chi2 = -999; // meta index
    double chi2_min_M1 = 40000;//200^2 > 178^2 ~ (x-top_mass)^2 //NOTE: will bound the range of M1!!!
    //pick the wanted one (with index)
    for(int i=0; i<top_fcnh_candidates.size(); ++i){ if(top_fcnh_chi2[i]<chi2_min_M1){ index_M1_min_chi2 = i; chi2_min_M1 = top_fcnh_chi2[i];} }
    TLorentzVector null;
    if(index_M1_min_chi2 == -999){ M1 = -999; index_q = -999; jet_q.SetPtEtaPhiM(0, 0, 0, 0); return null;}//No suitable candidate
    else{
        index_q = indices[index_M1_min_chi2];
        //jet_q = Jets[index_M1_min_chi2]; // bug: you should not use the meta-index. Please use the  real index for Jets!!
        jet_q = Jets[index_q];
        M1 = top_fcnh_candidates[index_M1_min_chi2].M();
        return top_fcnh_candidates[index_M1_min_chi2];
    }
}
//}}}
//}}}
//leading jet method related{{{
double GetM1M2_ratio(double M1, double M2){
    double ratio = abs(M1/M2-1) + abs(M2/M1-1);
    return ratio;
}
//}}}
//get_matching_result{{{
std::vector<bool> get_matching_result(std::vector<int> reco, std::vector<int> gen)
{
    bool tbwIsCorrectlyMatched = check_tbwIsCorrectlyMatched(reco, gen); 
    bool tqhIsCorrectlyMatched = check_tqhIsCorrectlyMatched(reco[3], gen[3]);
    bool sigIsCorrectlyReco    = tbwIsCorrectlyMatched && tqhIsCorrectlyMatched;
    bool bjet_is_matched       = check_bjet_matching(reco[0], gen[0]);
    std::vector<bool> matching_result = {tbwIsCorrectlyMatched, tqhIsCorrectlyMatched, sigIsCorrectlyReco, bjet_is_matched};
    return matching_result;
}


bool check_bjet_matching(int index_qjet, int index_gen_matched){
    bool tqhIsCorrectlyMatched = (index_qjet == index_gen_matched) && (index_gen_matched != -999);
    return tqhIsCorrectlyMatched;
}

bool check_tqhIsCorrectlyMatched(int index_qjet, int index_gen_matched){
    bool tqhIsCorrectlyMatched = (index_qjet == index_gen_matched) && (index_gen_matched != -999);
    return tqhIsCorrectlyMatched;
}


bool check_tbwIsCorrectlyMatched(std::vector<int> index_jet_chi2, std::vector<int> index_gen_matched)
{
    bool bjet_is_from_tbw = index_jet_chi2[0] == index_gen_matched[0];
    bool are_From_W_boson = (  ((index_jet_chi2[1] == index_gen_matched[1]) && (index_jet_chi2[2] == index_gen_matched[2]))\
                            || ((index_jet_chi2[1] == index_gen_matched[2]) && (index_jet_chi2[2] == index_gen_matched[1])) );
    bool tbwIsCorrectlyMatched = bjet_is_from_tbw && are_From_W_boson;
    return tbwIsCorrectlyMatched;
}


bool check_tbwIsCorrectlyMatched(std::vector<int> index_jet_chi2, std::vector<int> jetIndex_momPdgID_is_wboson, int index_bjet, int jetIndex_is_bquarkFromSMtop) // old version
{
    bool are_From_W_boson;

    if(jetIndex_momPdgID_is_wboson[0] == -999 || jetIndex_momPdgID_is_wboson[1] == -999)
        are_From_W_boson = false;
    else
        are_From_W_boson = ( ((index_jet_chi2[0] == jetIndex_momPdgID_is_wboson[0]) && (index_jet_chi2[1] == jetIndex_momPdgID_is_wboson[1]))\
                            || ((index_jet_chi2[0] == jetIndex_momPdgID_is_wboson[1]) && (index_jet_chi2[1] == jetIndex_momPdgID_is_wboson[0])) );
    bool tbwIsCorrectlyMatched = (index_bjet == jetIndex_is_bquarkFromSMtop) && are_From_W_boson;
    return tbwIsCorrectlyMatched;
}
//}}}
//check matching{{{
// ### Alternative Mother info
bool isMatched_with_Gen_tbw(std::vector<int> *GenPartInfo_PdgID, std::vector<int> *GenPartInfo_MomPdgID, int index_bquark, int index_quark1, int index_quark2){
    bool isCorrect = \
    index_bquark > 0 &&\
    index_quark1 > 0 &&\
    index_quark2 > 0 &&\
    abs(GenPartInfo_PdgID->at(index_bquark)) == 5 &&\
    abs(GenPartInfo_MomPdgID->at(index_bquark)) == 6 &&\
    abs(GenPartInfo_MomPdgID->at(index_quark1)) == 24 &&\
    abs(GenPartInfo_MomPdgID->at(index_quark2)) == 24 &&\
    (index_quark1 != index_quark2);

    return isCorrect;
}
//### genMatching
bool is_this_tqh_quark(TLorentzVector jet, Int_t GenPartInfo_size, std::vector<int> *GenPartInfo_MomPdgID, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_Status, std::vector<int> *GenPartInfo_PdgID){
    TLorentzVector truelove;//we are looking for the right genParticle to match the jet.
    int index = -999, truelove_PdgID = -999; double delta_R_min = 999.;
    for(int i=0; i<GenPartInfo_size; i++){
        if( abs(GenPartInfo_Status->at(i)) != 23 ) continue;//remove incoming/intermediate particles
        if( abs(GenPartInfo_PdgID->at(i)) > 6 ) continue;//exclude top quark & other particles
        //--------------------
        TLorentzVector genParticle;
        genParticle.SetPtEtaPhiM(GenPartInfo_Pt->at(i), GenPartInfo_Eta->at(i), GenPartInfo_Phi->at(i), GenPartInfo_Mass->at(i));
        double delta_R = jet.DeltaR(genParticle);
        //select quark & require min(deltaR)
        if( delta_R < 0.4 && delta_R < delta_R_min){
            index = i;//record the matched genParticle
            delta_R_min = delta_R;
            truelove = genParticle;
            truelove_PdgID = GenPartInfo_PdgID->at(i);
        }
    }//end of gen loop

    bool is_tqh_quark = ( (index != -999) && (abs(GenPartInfo_PdgID->at(index))!=5) && (abs(GenPartInfo_MomPdgID->at(index))==6) );
    //if(index!=-999 && !is_tqh_quark) printf("[CHECK-is_this_tqh_quark] pdgID = %d, Mom_PdgID = %d\n", GenPartInfo_PdgID->at(index), GenPartInfo_MomPdgID->at(index));
    //else if(index!=-999)             printf("[CHECK-is_this_tqh_quark] Matched! pdgID = %d, Mom_PdgID = %d\n", GenPartInfo_PdgID->at(index), GenPartInfo_MomPdgID->at(index));
    //else                             printf(" No matched particles\n");
    return is_tqh_quark;

}
bool CheckBJetID(TLorentzVector jet, Int_t GenPartInfo_size, std::vector<int> *GenPartInfo_MomPdgID, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_Status, std::vector<int> *GenPartInfo_PdgID){
    TLorentzVector truelove;//we are looking for the right genParticle to match the jet.
    int index = -999, truelove_PdgID = -999; double delta_R_min = 999.;
    for(int i=0; i<GenPartInfo_size; i++){
        if( abs(GenPartInfo_Status->at(i)) != 23 ) continue;//remove incoming/intermediate particles
        if( abs(GenPartInfo_PdgID->at(i)) == 6 ) continue;//exclude top quark
        //--------------------
        TLorentzVector genParticle;
        genParticle.SetPtEtaPhiM(GenPartInfo_Pt->at(i), GenPartInfo_Eta->at(i), GenPartInfo_Phi->at(i), GenPartInfo_Mass->at(i));
        double delta_R = jet.DeltaR(genParticle);
        //select quark & require min(deltaR)
        if( abs(GenPartInfo_PdgID->at(i)) < 7 && delta_R < 0.4 && delta_R < delta_R_min){
        //if( abs(GenPartInfo_PdgID->at(i)) < 7 && delta_R < delta_R_min){
            index = i;//record the matched genParticle
            delta_R_min = delta_R;
            truelove = genParticle;
            truelove_PdgID = GenPartInfo_PdgID->at(i);
        }
    }//end of gen loop
    // bjet is bquark coming from top
    bool bjet_is_bquark = ( abs(truelove_PdgID) == 5 && abs(GenPartInfo_MomPdgID->at(index)) == 6);
    return bjet_is_bquark;
}
TLorentzVector GetGenParticle(TLorentzVector jet, Int_t GenPartInfo_size, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_Status, std::vector<int> *GenPartInfo_PdgID, std::vector<int> &index_GenParticles, int &genParticle_PdgID){
    // GenMatching: find the gen particle (MC truth) for each jet (reconstructed). 
    // pdgID: (1, 2, 3, 4, 5, 6) = (d, u, s, c, b, t)
    TLorentzVector truelove;//we are looking for the right genParticle to match the jet.
    int index = -1, truelove_PdgID = -999; double delta_R_min = 999.;
    for(int i=0; i<GenPartInfo_size; i++){
        if( abs(GenPartInfo_Status->at(i)) != 23 ) continue;//remove incoming/intermediate particles
        if( abs(GenPartInfo_PdgID->at(i)) == 6 ) continue;//exclude top quark
        //bool isAvailable = checkAvailability(i, index_GenParticles); if(!isAvailable) continue;//the truelove of previous jets shall not become another one's truelove.
        //--------------------
        TLorentzVector genParticle;
        genParticle.SetPtEtaPhiM(GenPartInfo_Pt->at(i), GenPartInfo_Eta->at(i), GenPartInfo_Phi->at(i), GenPartInfo_Mass->at(i));
        double delta_R = jet.DeltaR(genParticle);
        //select quark & require min(deltaR)
        if( abs(GenPartInfo_PdgID->at(i)) < 7 && delta_R < 0.4 && delta_R < delta_R_min){
            index = i;//record the matched genParticle
            delta_R_min = delta_R;
            truelove = genParticle;
            truelove_PdgID = GenPartInfo_PdgID->at(i);
        }
    }//end of gen loop
    //---------- store particle info ----------//
    index_GenParticles.push_back(index);
    genParticle_PdgID = truelove_PdgID;
    return truelove;
}
//}}}
//obtain_deltaR{{{
double obtain_deltaR(bool bool_print, TLorentzVector jet, int &Matched_PdgID, int &Matched_MomPdgID, Int_t GenPartInfo_size, std::vector<int> *GenPartInfo_MomPdgID, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_Status, std::vector<int> *GenPartInfo_PdgID){
    TLorentzVector truelove;//we are looking for the right genParticle to match the jet.
    int index = -999, truelove_PdgID = -999; double delta_R_min = 999.;
    for(int i=0; i<GenPartInfo_size; i++){
        if( abs(GenPartInfo_Status->at(i)) != 23 ) continue;//remove incoming/intermediate particles
        if( abs(GenPartInfo_PdgID->at(i)) > 6 ) continue;//exclude top quark & other particles
        //--------------------
        TLorentzVector genParticle;
        genParticle.SetPtEtaPhiM(GenPartInfo_Pt->at(i), GenPartInfo_Eta->at(i), GenPartInfo_Phi->at(i), GenPartInfo_Mass->at(i));
        double delta_R = jet.DeltaR(genParticle);
        //select quark & require min(deltaR)
        if( delta_R < 0.4 && delta_R < delta_R_min){
            index = i;//record the matched genParticle
            delta_R_min = delta_R;
            truelove = genParticle;
            truelove_PdgID = GenPartInfo_PdgID->at(i);
        }
    }//end of gen loop

    if(index>0){
        Matched_PdgID = GenPartInfo_PdgID->at(index);
        Matched_MomPdgID = GenPartInfo_MomPdgID->at(index);
        if(bool_print){
            printf("deltaR = %6.3f, ", delta_R_min);
            printf("MomPdgID = %5d, ", GenPartInfo_MomPdgID->at(index));
            printf("PdgID = %3d\n", GenPartInfo_PdgID->at(index));
        }
    } else{
        Matched_PdgID = -999;
        Matched_MomPdgID = -999;
        if(bool_print) printf("[check-obtain-deltaR] Not found.\n");
    }

    return delta_R_min;
}
//}}}
// print gen info{{{
void obtain_gen_matched_ID(bool bool_print, TLorentzVector jet, int &index_gen, int &Matched_PdgID, int &Matched_MomPdgID, Int_t GenPartInfo_size, std::vector<int> *GenPartInfo_MomPdgID, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_Status, std::vector<int> *GenPartInfo_PdgID){
    TLorentzVector truelove;//we are looking for the right genParticle to match the jet.
    int index = -999, truelove_PdgID = -999; double delta_R_min = 999.; index_gen = -1;
    for(int i=0; i<GenPartInfo_size; i++){
        if( abs(GenPartInfo_Status->at(i)) != 23 ) continue;//remove incoming/intermediate particles
        //if( abs(GenPartInfo_PdgID->at(i)) > 6 ) continue;//exclude top quark & other particles
        if( abs(GenPartInfo_PdgID->at(i)) == 6 || abs(GenPartInfo_PdgID->at(i)) > 21) continue;//exclude top quark & photon, Z, W, H
        //--------------------
        TLorentzVector genParticle;
        genParticle.SetPtEtaPhiM(GenPartInfo_Pt->at(i), GenPartInfo_Eta->at(i), GenPartInfo_Phi->at(i), GenPartInfo_Mass->at(i));
        double delta_R = jet.DeltaR(genParticle);
        //select quark & require min(deltaR)
        if( delta_R < 0.4 && delta_R < delta_R_min){
            index = i;//record the matched genParticle
            delta_R_min = delta_R;
            truelove = genParticle;
            truelove_PdgID = GenPartInfo_PdgID->at(i);
            index_gen = i;
        }
    }//end of gen loop


    if(index>0){ // index = 0 or 1 corresponds to incoming particles
        Matched_PdgID = GenPartInfo_PdgID->at(index);
        Matched_MomPdgID = GenPartInfo_MomPdgID->at(index);
        //if(bool_print){{{
        if(bool_print){
            //printf("Status = %3d, ", GenPartInfo_Status->at(index));
            printf("Pt = %6.2f, ", GenPartInfo_Pt->at(index));
            printf("Eta = %9.2f, ", GenPartInfo_Eta->at(index));
            printf("Phi = %6.2f, ", GenPartInfo_Phi->at(index));
            printf("Mass = %6.2f, ", GenPartInfo_Mass->at(index));
            printf("deltaR = %6.3f, ", delta_R_min);
            printf("MomPdgID = %5d, ", GenPartInfo_MomPdgID->at(index));
            printf("PdgID = %3d\n", GenPartInfo_PdgID->at(index));
        }
        //}}}
    } else{
        Matched_PdgID = -999;
        Matched_MomPdgID = -999;
        if(bool_print) printf("[check-obtain] Not found.\n");
    }
}
//}}}
//### report{{{
void report_rate(const char* name, int counter, int N){
    double eff = (double) counter / (double) N;
    double n = (double) N;
    double p = (double) counter / (double) N;
    double q = 1 - p;
    double error = sqrt(p * q / n);
    printf("[INFO] report_rate::%-60s = %7d / %7d ( %6.2f%% +/- %6.2f%% )\n", name, counter, N, 100.*eff, 100.*error);
}
void kinematics_info(int index, TLorentzVector Particle, TLorentzVector gen){
    double deltaR = Particle.DeltaR(gen);
    printf("(%2d)                            ", index);
    printf("Pt = %6.2f, ", Particle.Pt());
    printf("Eta = %9.2f, ", Particle.Eta());
    printf("Phi = %6.2f, ", Particle.Phi());
    printf("Mass = %6.2f, ", Particle.M());
    printf("deltaR = %6.3f\n", deltaR);
}
void kinematics_info(const char* Title, TLorentzVector Particle, int index){
    printf("Pt = %6.2f, ", Particle.Pt());
    printf("Eta = %9.2f, ", Particle.Eta());
    printf("Phi = %6.2f, ", Particle.Phi());
    printf("Mass = %6.2f, ", Particle.M());
    printf("index = %5d\n", index);
}
void kinematics_report(const char* recoTitle, TLorentzVector recoParticle, int id_recoParticle, const char* genTitle, TLorentzVector genParticle, int genParticle_PdgID){
        double delta_R = genParticle.DeltaR(recoParticle);
        printf("(%s) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f, Energy = %6.2f, Mass = %6.2f, id = %d\n", recoTitle, recoParticle.Pt(), recoParticle.Eta(), recoParticle.Phi(), recoParticle.Energy(), recoParticle.M(), id_recoParticle);
        printf("(%s) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f, Energy = %6.2f, Mass = %6.2f, id = %d, delta_R = %6.2f\n", genTitle, genParticle.Pt(), genParticle.Eta(), genParticle.Phi(), genParticle.Energy(), genParticle.M(), genParticle_PdgID, delta_R);
}
void hist_bin_fraction(TH1D *hist, const char* title, int entries_matched_in_bin2){
    int total_entries = hist->GetEntries();
    printf("[INFO-hist] %s (%d)\n", title, total_entries);
    printf("[INFO-hist] fractions = ");
    for(int i=0; i<3; ++i){
        if(i<2) printf("%6.4f, ", (double) hist->GetBinContent(i+1) / (double)total_entries);
        else{
            int accumulation = 0;
            for(int j=i; j<10; ++j) accumulation += (double) hist->GetBinContent(j+1);
            printf("%6.4f; ", i, (double) accumulation / (double)total_entries);
            printf("%6.4f, ", i, (double) entries_matched_in_bin2 / (double)total_entries);
            printf("%6.4f\n", i, (double) entries_matched_in_bin2 / (double)hist->GetBinContent(2));
        }
    }
}
void hist_report(TH1D *hist, const char* chi2_type){
    double mean = hist->GetMean();
    double sigma = hist->GetMeanError();
    double sigma2 = hist->GetRMS();
    int bin_maximum = hist->GetMaximumBin();
    double maxbin_lower_edge = hist->GetBinLowEdge(bin_maximum);
    double bin_width = hist->GetBinWidth(bin_maximum);
    double maxbin_upper_edge = maxbin_lower_edge + bin_width;

    //printf("[INFO] chi2 %s: mean = %6.2f, sigma = %6.2f, %6.2f\n", chi2_type, mean, sigma, sigma2);
    printf("[INFO] chi2 %s: mean = %6.2f, sigma = %6.2f, mode lies in [%6.2f, %6.2f]\n", chi2_type, mean, sigma2, maxbin_lower_edge, maxbin_upper_edge);
}
//}}}
//### mkplots{{{
void MakeTwoPlots(TCanvas *c1, TH1D* hist_gen, TH1D* hist_reco, TLegend *legend, const char* name){
    hist_reco->SetStats(0);
    hist_reco->SetLineWidth(2);
    hist_reco->SetLineColor(kBlue);
    hist_reco->Draw();
    hist_gen->SetLineWidth(2);
    hist_gen->SetLineColor(kRed);
    hist_gen->Draw("same;hist");

    legend->Clear();
    legend->AddEntry(hist_gen, "gen-level", "l");
    legend->AddEntry(hist_reco, "reco", "l");
    legend->SetLineColor(0);
    legend->Draw("same");

    c1->SaveAs(Form("%s/%s", output_dir, name));
}
void MakePlots(TCanvas *c1, TH1D* hist, const char* title, const char* outputFile){
    hist->Draw();
    hist->SetTitle(title);
    hist->SetXTitle(title);
    hist->SetYTitle("Entries");
    hist->GetYaxis()->SetTitleOffset(1.4);
    hist->Write();
    c1->SaveAs(outputFile);
}
void MakeFinalPlots(TCanvas *c1, TH1D* hist_simple, TH1D* hist_modified, TH1D*hist_yfyj, TLegend *legend, const char* name){
    double max;
    double scale = 1.2;
    double max_simple = hist_simple -> GetMaximum();
    double max_modified = hist_modified -> GetMaximum();
    double max_yfyj = hist_yfyj -> GetMaximum();
    max = (max_modified > max_simple) ? max_modified : max_simple;
    max = (max_yfyj > max_modified) ? max_yfyj : max_modified;
    double overflow_simple = hist_simple->GetBinContent(hist_simple->GetNbinsX() + 1);
    double overflow_modified = hist_modified->GetBinContent(hist_modified->GetNbinsX() + 1);
    double overflow_yfyj = hist_yfyj->GetBinContent(hist_yfyj->GetNbinsX() + 1);
    //max = (max > overflow_simple) ? max : overflow_simple;
    //max = (max > overflow_modified) ? max : overflow_modified;
    //max = (max > overflow_yfyj) ? max : overflow_yfyj;

    hist_simple->SetStats(0);
    hist_simple->SetMaximum(max*scale);
    hist_simple->SetLineWidth(2);
    hist_simple->SetLineColor(kBlue);
    hist_simple->GetXaxis()->SetRange(1, hist_simple->GetNbinsX() + 1);
    hist_simple->Draw();
    hist_modified->SetLineWidth(2);
    hist_modified->SetLineColor(kRed);
    hist_modified->GetXaxis()->SetRange(1, hist_modified->GetNbinsX() + 1);
    hist_modified->Draw("same");
    hist_yfyj->SetLineWidth(2);
    hist_yfyj->SetLineColor(kGreen+4);
    hist_yfyj->GetXaxis()->SetRange(1, hist_yfyj->GetNbinsX() + 1);
    hist_yfyj->Draw("same");

    legend->Clear();
    legend->SetTextSize(0.03);
    legend->AddEntry(hist_yfyj, "leading-jets method", "l");
    legend->AddEntry(hist_simple, "simple #chi^{2} method", "l");
    legend->AddEntry(hist_modified, "modified #chi^{2} method", "l");
    legend->SetLineColor(0);
    legend->Draw("same");

    c1->SaveAs(Form("%s/%s", output_dir, name));
}
//###}}}
//### bool functions{{{
bool isThisMCsignal(char* dataset){
    if((string)dataset == "ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-T2HJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    return false;
}
bool isThisDataOrNot(char* dataset){
    if((string)dataset == "DoubleEG_B") return true;
    if((string)dataset == "DoubleEG_C") return true;
    if((string)dataset == "DoubleEG_D") return true;
    if((string)dataset == "DoubleEG_E") return true;
    if((string)dataset == "DoubleEG_F") return true;
    return false;
}
bool isThisMultiFile(char* dataset){
    if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return true;
    if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return true;
    return false;
}
//###}}}
//### flashggStdTree {{{
flashggStdTreeParameters::flashggStdTreeParameters(){
    GenPartInfo_Pt = new std::vector<float>;
    GenPartInfo_Eta = new std::vector<float>;
    GenPartInfo_Phi = new std::vector<float>;
    GenPartInfo_Mass = new std::vector<float>;
    GenPartInfo_PdgID = new std::vector<int>;
    GenPartInfo_Status = new std::vector<int>;
    GenPartInfo_nMo = new std::vector<int>;
    GenPartInfo_nDa = new std::vector<int>;
    GenPartInfo_isHardProcess = new std::vector<bool>;
    GenPartInfo_fromHardProcessFinalState = new std::vector<bool>;
    GenPartInfo_isPromptFinalState = new std::vector<bool>;
    GenPartInfo_isDirectPromptTauDecayProductFinalState = new std::vector<bool>;
    GenPartInfo_MomPdgID = new std::vector<int>;
    GenPartInfo_MomStatus = new std::vector<int>;
    GenPartInfo_MomPt = new std::vector<float>;
    GenPartInfo_MomEta = new std::vector<float>;
    GenPartInfo_MomPhi = new std::vector<float>;
    GenPartInfo_MomMass = new std::vector<float>;
    //------------------------
    JetInfo_Pt = new std::vector<float>;
    JetInfo_Eta = new std::vector<float>;
    JetInfo_Phi = new std::vector<float>;
    JetInfo_Mass = new std::vector<float>;
    JetInfo_Energy = new std::vector<float>;
    JetInfo_GenFlavor = new std::vector<int>;
    JetInfo_pfDeepCSVJetTags_probb = new std::vector<float>;
    JetInfo_pfDeepCSVJetTags_probbb = new std::vector<float>;
    //------------------------
    ElecInfo_Charge = new std::vector<int>;
    ElecInfo_Pt = new std::vector<float>;
    ElecInfo_Eta = new std::vector<float>;
    ElecInfo_Phi = new std::vector<float>;
    ElecInfo_Energy = new std::vector<float>;
    ElecInfo_EtaSC = new std::vector<float>;
    ElecInfo_PhiSC = new std::vector<float>;
    ElecInfo_GsfTrackDz = new std::vector<float>;
    ElecInfo_GsfTrackDxy = new std::vector<float>;
    ElecInfo_EGMCutBasedIDVeto = new std::vector<bool>;
    ElecInfo_EGMCutBasedIDLoose = new std::vector<bool>;
    ElecInfo_EGMCutBasedIDMedium = new std::vector<bool>;
    ElecInfo_EGMCutBasedIDTight = new std::vector<bool>;
    ElecInfo_fggPhoVeto = new std::vector<bool>;
    ElecInfo_EnergyCorrFactor = new std::vector<float>;
    ElecInfo_EnergyPostCorrErr = new std::vector<float>;
    ElecInfo_EnergyPostCorrScaleUp = new std::vector<float>;
    ElecInfo_EnergyPostCorrScaleDown = new std::vector<float>;
    ElecInfo_EnergyPostCorrSmearUp = new std::vector<float>;
    ElecInfo_EnergyPostCorrSmearDown = new std::vector<float>;
    ElecInfo_GenMatch = new std::vector<bool>;
    ElecInfo_GenPdgID = new std::vector<int>;
    ElecInfo_GenPt = new std::vector<float>;
    ElecInfo_GenEta = new std::vector<float>;
    ElecInfo_GenPhi = new std::vector<float>;
    //------------------------
    MuonInfo_Charge = new std::vector<int>;
    MuonInfo_MuonType = new std::vector<float>;
    MuonInfo_Pt = new std::vector<float>;
    MuonInfo_Eta = new std::vector<float>;
    MuonInfo_Phi = new std::vector<float>;
    MuonInfo_Energy = new std::vector<float>;
    MuonInfo_BestTrackDz = new std::vector<float>;
    MuonInfo_BestTrackDxy = new std::vector<float>;
    MuonInfo_PFIsoDeltaBetaCorrR04 = new std::vector<float>;
    MuonInfo_TrackerBasedIsoR03 = new std::vector<float>;
    MuonInfo_CutBasedIdMedium = new std::vector<bool>;
    MuonInfo_CutBasedIdTight = new std::vector<bool>;
    MuonInfo_GenMatch = new std::vector<bool>;
    MuonInfo_GenPdgID = new std::vector<int>;
    MuonInfo_GenPt = new std::vector<float>;
    MuonInfo_GenEta = new std::vector<float>;
    MuonInfo_GenPhi = new std::vector<float>;
    //------------------------
}
flashggStdTreeParameters::~flashggStdTreeParameters(){
    delete GenPartInfo_Pt;
    delete GenPartInfo_Eta;
    delete GenPartInfo_Phi;
    delete GenPartInfo_Mass;
    delete GenPartInfo_PdgID;
    delete GenPartInfo_Status;
    delete GenPartInfo_nMo;
    delete GenPartInfo_nDa;
    delete GenPartInfo_isHardProcess;
    delete GenPartInfo_fromHardProcessFinalState;
    delete GenPartInfo_isPromptFinalState;
    delete GenPartInfo_isDirectPromptTauDecayProductFinalState;
    delete GenPartInfo_MomPdgID;
    delete GenPartInfo_MomStatus;
    delete GenPartInfo_MomPt;
    delete GenPartInfo_MomEta;
    delete GenPartInfo_MomPhi;
    delete GenPartInfo_MomMass;
    //------------------------
    delete JetInfo_Pt;
    delete JetInfo_Eta;
    delete JetInfo_Phi;
    delete JetInfo_Mass;
    delete JetInfo_Energy;
    delete JetInfo_GenFlavor;
    delete JetInfo_pfDeepCSVJetTags_probb;
    delete JetInfo_pfDeepCSVJetTags_probbb;
    //------------------------
    delete ElecInfo_Charge;
    delete ElecInfo_Pt;
    delete ElecInfo_Eta;
    delete ElecInfo_Phi;
    delete ElecInfo_Energy;
    delete ElecInfo_EtaSC;
    delete ElecInfo_PhiSC;
    delete ElecInfo_GsfTrackDz;
    delete ElecInfo_GsfTrackDxy;
    delete ElecInfo_EGMCutBasedIDVeto;
    delete ElecInfo_EGMCutBasedIDLoose;
    delete ElecInfo_EGMCutBasedIDMedium;
    delete ElecInfo_EGMCutBasedIDTight;
    delete ElecInfo_fggPhoVeto;
    delete ElecInfo_EnergyCorrFactor;
    delete ElecInfo_EnergyPostCorrErr;
    delete ElecInfo_EnergyPostCorrScaleUp;
    delete ElecInfo_EnergyPostCorrScaleDown;
    delete ElecInfo_EnergyPostCorrSmearUp;
    delete ElecInfo_EnergyPostCorrSmearDown;
    delete ElecInfo_GenMatch;
    delete ElecInfo_GenPdgID;
    delete ElecInfo_GenPt;
    delete ElecInfo_GenEta;
    delete ElecInfo_GenPhi;
    //------------------------
    delete MuonInfo_Charge;
    delete MuonInfo_MuonType;
    delete MuonInfo_Pt;
    delete MuonInfo_Eta;
    delete MuonInfo_Phi;
    delete MuonInfo_Energy;
    delete MuonInfo_BestTrackDz;
    delete MuonInfo_BestTrackDxy;
    delete MuonInfo_PFIsoDeltaBetaCorrR04;
    delete MuonInfo_TrackerBasedIsoR03;
    delete MuonInfo_CutBasedIdMedium;
    delete MuonInfo_CutBasedIdTight;
    delete MuonInfo_GenMatch;
    delete MuonInfo_GenPdgID;
    delete MuonInfo_GenPt;
    delete MuonInfo_GenEta;
    delete MuonInfo_GenPhi;
    //------------------------
}
flashggStdTreeReader::flashggStdTreeReader(void){
    printf("[INFO] Reading data...\n");
}
void flashggStdTreeReader::InitChain(const char* treeName){
    flashggStdTree = new TChain(treeName);
    printf("[INFO] flashggStdTreeReader::InitChain : Finished!\n");
}
void flashggStdTreeReader::AddSingleRootFile(char* input_file){
    flashggStdTree->Add(input_file);
    printf(Form("[INFO] flashggStdTreeReader::AddSingleRootFile : %s\n", input_file));
}
void flashggStdTreeReader::AddMultiRootFile(char* input_file){
    flashggStdTree->Add(Form("%s/*.root", input_file));
    printf("[INFO] flashggStdTreeReader::AddMultiRootFile : Finished!\n");
}
int flashggStdTreeReader::GetEntries(void){
    printf("[INFO] flashggStdTreeReader::GetEntries : %d\n", flashggStdTree->GetEntries());
    return flashggStdTree->GetEntries();
}
double flashggStdTreeReader::GetGenWeight(void){
    return EvtInfo_genweight;
}
TChain* flashggStdTreeReader::GetTChain(void){
    return flashggStdTree;
}
void flashggStdTreeReader::SetBranchAddresses(){
    flashggStdTree->SetBranchAddress("GenPartInfo.size", &GenPartInfo_size);
    flashggStdTree->SetBranchAddress("GenPartInfo.Pt", &GenPartInfo_Pt);
    flashggStdTree->SetBranchAddress("GenPartInfo.Eta", &GenPartInfo_Eta);
    flashggStdTree->SetBranchAddress("GenPartInfo.Phi", &GenPartInfo_Phi);
    flashggStdTree->SetBranchAddress("GenPartInfo.Mass", &GenPartInfo_Mass);
    flashggStdTree->SetBranchAddress("GenPartInfo.PdgID", &GenPartInfo_PdgID);
    flashggStdTree->SetBranchAddress("GenPartInfo.Status", &GenPartInfo_Status);
    flashggStdTree->SetBranchAddress("GenPartInfo.nMo", &GenPartInfo_nMo);
    flashggStdTree->SetBranchAddress("GenPartInfo.nDa", &GenPartInfo_nDa);
    flashggStdTree->SetBranchAddress("GenPartInfo.isHardProcess", &GenPartInfo_isHardProcess);
    flashggStdTree->SetBranchAddress("GenPartInfo.fromHardProcessFinalState", &GenPartInfo_fromHardProcessFinalState);
    flashggStdTree->SetBranchAddress("GenPartInfo.isPromptFinalState", &GenPartInfo_isPromptFinalState);
    flashggStdTree->SetBranchAddress("GenPartInfo.isDirectPromptTauDecayProductFinalState", &GenPartInfo_isDirectPromptTauDecayProductFinalState);
    flashggStdTree->SetBranchAddress("GenPartInfo.MomPdgID", &GenPartInfo_MomPdgID);
    flashggStdTree->SetBranchAddress("GenPartInfo.MomStatus", &GenPartInfo_MomStatus);
    flashggStdTree->SetBranchAddress("GenPartInfo.MomPt", &GenPartInfo_MomPt);
    flashggStdTree->SetBranchAddress("GenPartInfo.MomEta", &GenPartInfo_MomEta);
    flashggStdTree->SetBranchAddress("GenPartInfo.MomPhi", &GenPartInfo_MomPhi);
    flashggStdTree->SetBranchAddress("GenPartInfo.MomMass", &GenPartInfo_MomMass);
    //------------------------
    flashggStdTree->SetBranchAddress("EvtInfo.passTrigger", &EvtInfo_passTrigger);
    flashggStdTree->SetBranchAddress("EvtInfo.NPu", &EvtInfo_NPu);
    flashggStdTree->SetBranchAddress("EvtInfo.NVtx", &EvtInfo_NVtx);
    flashggStdTree->SetBranchAddress("EvtInfo.Rho", &EvtInfo_Rho);
    flashggStdTree->SetBranchAddress("EvtInfo.genweight", &EvtInfo_genweight);
    //------------------------
    flashggStdTree->SetBranchAddress("DiPhoInfo.mass", &DiPhoInfo_mass);
    flashggStdTree->SetBranchAddress("DiPhoInfo.pt", &DiPhoInfo_pt);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadPt", &DiPhoInfo_leadPt);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadEta", &DiPhoInfo_leadEta);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadPhi", &DiPhoInfo_leadPhi);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadE", &DiPhoInfo_leadE);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadhoe", &DiPhoInfo_leadhoe);
    flashggStdTree->SetBranchAddress("DiPhoInfo.leadIDMVA", &DiPhoInfo_leadIDMVA);
    //------------------------
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadPt", &DiPhoInfo_subleadPt);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadEta", &DiPhoInfo_subleadEta);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadPhi", &DiPhoInfo_subleadPhi);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadE", &DiPhoInfo_subleadE);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadhoe", &DiPhoInfo_subleadhoe);
    flashggStdTree->SetBranchAddress("DiPhoInfo.subleadIDMVA", &DiPhoInfo_subleadIDMVA);
    //------------------------
    flashggStdTree->SetBranchAddress("jets_size", &jets_size);
    flashggStdTree->SetBranchAddress("JetInfo.Pt", &JetInfo_Pt);
    flashggStdTree->SetBranchAddress("JetInfo.Eta", &JetInfo_Eta);
    flashggStdTree->SetBranchAddress("JetInfo.Phi", &JetInfo_Phi);
    flashggStdTree->SetBranchAddress("JetInfo.Mass", &JetInfo_Mass);
    flashggStdTree->SetBranchAddress("JetInfo.Energy", &JetInfo_Energy);
    flashggStdTree->SetBranchAddress("JetInfo.GenFlavor", &JetInfo_GenFlavor);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probb", &JetInfo_pfDeepCSVJetTags_probb);
    flashggStdTree->SetBranchAddress("JetInfo.pfDeepCSVJetTags_probbb", &JetInfo_pfDeepCSVJetTags_probbb);
    //------------------------
    flashggStdTree->SetBranchAddress("ElecInfo.Size", &ElecInfo_Size);
    flashggStdTree->SetBranchAddress("ElecInfo.Charge", &ElecInfo_Charge);
    flashggStdTree->SetBranchAddress("ElecInfo.Pt", &ElecInfo_Pt);
    flashggStdTree->SetBranchAddress("ElecInfo.Eta", &ElecInfo_Eta);
    flashggStdTree->SetBranchAddress("ElecInfo.Phi", &ElecInfo_Phi);
    flashggStdTree->SetBranchAddress("ElecInfo.Energy", &ElecInfo_Energy);
    flashggStdTree->SetBranchAddress("ElecInfo.EtaSC", &ElecInfo_EtaSC);
    flashggStdTree->SetBranchAddress("ElecInfo.PhiSC", &ElecInfo_PhiSC);
    flashggStdTree->SetBranchAddress("ElecInfo.GsfTrackDz", &ElecInfo_GsfTrackDz);
    flashggStdTree->SetBranchAddress("ElecInfo.GsfTrackDxy", &ElecInfo_GsfTrackDxy);
    flashggStdTree->SetBranchAddress("ElecInfo.EGMCutBasedIDVeto", &ElecInfo_EGMCutBasedIDVeto);
    flashggStdTree->SetBranchAddress("ElecInfo.EGMCutBasedIDLoose", &ElecInfo_EGMCutBasedIDLoose);
    flashggStdTree->SetBranchAddress("ElecInfo.EGMCutBasedIDMedium", &ElecInfo_EGMCutBasedIDMedium);
    flashggStdTree->SetBranchAddress("ElecInfo.EGMCutBasedIDTight", &ElecInfo_EGMCutBasedIDTight);
    flashggStdTree->SetBranchAddress("ElecInfo.fggPhoVeto", &ElecInfo_fggPhoVeto);
    flashggStdTree->SetBranchAddress("ElecInfo.EnergyCorrFactor", &ElecInfo_EnergyCorrFactor);
    flashggStdTree->SetBranchAddress("ElecInfo.EnergyPostCorrErr", &ElecInfo_EnergyPostCorrErr);
    flashggStdTree->SetBranchAddress("ElecInfo.EnergyPostCorrScaleUp", &ElecInfo_EnergyPostCorrScaleUp);
    flashggStdTree->SetBranchAddress("ElecInfo.EnergyPostCorrScaleDown", &ElecInfo_EnergyPostCorrScaleDown);
    flashggStdTree->SetBranchAddress("ElecInfo.EnergyPostCorrSmearUp", &ElecInfo_EnergyPostCorrSmearUp);
    flashggStdTree->SetBranchAddress("ElecInfo.EnergyPostCorrSmearDown", &ElecInfo_EnergyPostCorrSmearDown);
    flashggStdTree->SetBranchAddress("ElecInfo.GenMatch", &ElecInfo_GenMatch);
    flashggStdTree->SetBranchAddress("ElecInfo.GenPdgID", &ElecInfo_GenPdgID);
    flashggStdTree->SetBranchAddress("ElecInfo.GenPt", &ElecInfo_GenPt);
    flashggStdTree->SetBranchAddress("ElecInfo.GenEta", &ElecInfo_GenEta);
    flashggStdTree->SetBranchAddress("ElecInfo.GenPhi", &ElecInfo_GenPhi);
    //------------------------
    flashggStdTree->SetBranchAddress("MuonInfo.Size", &MuonInfo_Size);
    flashggStdTree->SetBranchAddress("MuonInfo.Charge", &MuonInfo_Charge);
    flashggStdTree->SetBranchAddress("MuonInfo.MuonType", &MuonInfo_MuonType);
    flashggStdTree->SetBranchAddress("MuonInfo.Pt", &MuonInfo_Pt);
    flashggStdTree->SetBranchAddress("MuonInfo.Eta", &MuonInfo_Eta);
    flashggStdTree->SetBranchAddress("MuonInfo.Phi", &MuonInfo_Phi);
    flashggStdTree->SetBranchAddress("MuonInfo.Energy", &MuonInfo_Energy);
    flashggStdTree->SetBranchAddress("MuonInfo.BestTrackDz", &MuonInfo_BestTrackDz);
    flashggStdTree->SetBranchAddress("MuonInfo.BestTrackDxy", &MuonInfo_BestTrackDxy);
    flashggStdTree->SetBranchAddress("MuonInfo.PFIsoDeltaBetaCorrR04", &MuonInfo_PFIsoDeltaBetaCorrR04);
    flashggStdTree->SetBranchAddress("MuonInfo.TrackerBasedIsoR03", &MuonInfo_TrackerBasedIsoR03);
    flashggStdTree->SetBranchAddress("MuonInfo.CutBasedIdMedium", &MuonInfo_CutBasedIdMedium);
    flashggStdTree->SetBranchAddress("MuonInfo.CutBasedIdTight", &MuonInfo_CutBasedIdTight);
    flashggStdTree->SetBranchAddress("MuonInfo.GenMatch", &MuonInfo_GenMatch);
    flashggStdTree->SetBranchAddress("MuonInfo.GenPdgID", &MuonInfo_GenPdgID);
    flashggStdTree->SetBranchAddress("MuonInfo.GenPt", &MuonInfo_GenPt);
    flashggStdTree->SetBranchAddress("MuonInfo.GenEta", &MuonInfo_GenEta);
    flashggStdTree->SetBranchAddress("MuonInfo.GenPhi", &MuonInfo_GenPhi);
    //------------------------
    flashggStdTree->SetBranchAddress("MetInfo.Pt", &MetInfo_Pt);
    flashggStdTree->SetBranchAddress("MetInfo.Phi", &MetInfo_Phi);
    flashggStdTree->SetBranchAddress("MetInfo.Px", &MetInfo_Px);
    flashggStdTree->SetBranchAddress("MetInfo.Py", &MetInfo_Py);
    flashggStdTree->SetBranchAddress("MetInfo.SumET", &MetInfo_SumET);
    //------------------------
    printf("[INFO] flashggStdTreeReader::SetBranchAddresses : Finished!\n");
}
///### }}}
//### myTreeClass {{{
myParameters::myParameters(){
    JetInfo_jet_pt_selection = new std::vector<float>;
    JetInfo_jet_eta_selection = new std::vector<float>;
    JetInfo_jet_phi_selection = new std::vector<float>;
    JetInfo_jet_energy_selection = new std::vector<float>;
    JetInfo_jet_diphoton_deltaR_selection = new std::vector<float>;
    JetInfo_jet_leadingPhoton_deltaR_selection = new std::vector<float>;
    JetInfo_jet_subleadingPhoton_deltaR_selection = new std::vector<float>;
    JetInfo_jet_pfDeepCSVJetTags_probb_selection = new std::vector<float>;
    JetInfo_jet_pfDeepCSVJetTags_probbb_selection = new std::vector<float>;
    ElecInfo_electron_pt_selection = new std::vector<float>;
    ElecInfo_electron_eta_selection = new std::vector<float>;
    ElecInfo_electron_phi_selection = new std::vector<float>;
    ElecInfo_electron_energy_selection = new std::vector<float>;
    ElecInfo_electron_diphoton_deltaR_selection = new std::vector<float>;
    ElecInfo_electron_leadingPhoton_deltaR_selection = new std::vector<float>;
    ElecInfo_electron_subleadingPhoton_deltaR_selection = new std::vector<float>;
    MuonInfo_muon_pt_selection = new std::vector<float>;
    MuonInfo_muon_eta_selection = new std::vector<float>;
    MuonInfo_muon_phi_selection = new std::vector<float>;
    MuonInfo_muon_energy_selection = new std::vector<float>;
    MuonInfo_muon_diphoton_deltaR_selection = new std::vector<float>;
    MuonInfo_muon_leadingPhoton_deltaR_selection = new std::vector<float>;
    MuonInfo_muon_subleadingPhoton_deltaR_selection = new std::vector<float>;
}
myParameters::~myParameters(){
    delete JetInfo_jet_pt_selection;
    delete JetInfo_jet_eta_selection;
    delete JetInfo_jet_phi_selection;
    delete JetInfo_jet_energy_selection;
    delete JetInfo_jet_diphoton_deltaR_selection;
    delete JetInfo_jet_leadingPhoton_deltaR_selection;
    delete JetInfo_jet_subleadingPhoton_deltaR_selection;
    delete JetInfo_jet_pfDeepCSVJetTags_probb_selection;
    delete JetInfo_jet_pfDeepCSVJetTags_probbb_selection;
    delete ElecInfo_electron_pt_selection;
    delete ElecInfo_electron_eta_selection;
    delete ElecInfo_electron_phi_selection;
    delete ElecInfo_electron_energy_selection;
    delete ElecInfo_electron_diphoton_deltaR_selection;
    delete ElecInfo_electron_leadingPhoton_deltaR_selection;
    delete ElecInfo_electron_subleadingPhoton_deltaR_selection;
    delete MuonInfo_muon_pt_selection;
    delete MuonInfo_muon_eta_selection;
    delete MuonInfo_muon_phi_selection;
    delete MuonInfo_muon_energy_selection;
    delete MuonInfo_muon_diphoton_deltaR_selection;
    delete MuonInfo_muon_leadingPhoton_deltaR_selection;
    delete MuonInfo_muon_subleadingPhoton_deltaR_selection;
}

void myTreeClass::InitTree(){
    mytree =  new TTree("mytree", "mytree");
}
void myTreeClass::MakeNewBranchAddresses(){
    mytree -> Branch("Mass_w_candidate_chi2_simple", &Mass_w_candidate_chi2_simple, "Mass_w_candidate_chi2_simple/F");
    mytree -> Branch("Mass_top_candidate_chi2_simple", &Mass_top_candidate_chi2_simple, "Mass_top_candidate_chi2_simple/F");
    mytree -> Branch("Mass_w_candidate_chi2_modified", &Mass_w_candidate_chi2_modified, "Mass_w_candidate_chi2_modified/F");
    mytree -> Branch("Mass_top_candidate_chi2_modified", &Mass_top_candidate_chi2_modified, "Mass_top_candidate_chi2_modified/F");
    mytree -> Branch("Mass_gen_w_candidate_chi2_simple", &Mass_gen_w_candidate_chi2_simple, "Mass_gen_w_candidate_chi2_simple/F");
    mytree -> Branch("Mass_gen_top_candidate_chi2_simple", &Mass_gen_top_candidate_chi2_simple, "Mass_gen_top_candidate_chi2_simple/F");
    mytree -> Branch("Mass_gen_w_candidate_chi2_modified", &Mass_gen_w_candidate_chi2_modified, "Mass_gen_w_candidate_chi2_modified/F");
    mytree -> Branch("Mass_gen_top_candidate_chi2_modified", &Mass_gen_top_candidate_chi2_modified, "Mass_gen_top_candidate_chi2_modified/F");
    //------------------------
    mytree -> Branch("EvtInfo_totalEntry_before_preselection", &EvtInfo_totalEntry_before_preselection, "EvtInfo_totalEntry_before_preselection/I");
    mytree -> Branch("EvtInfo_NormalizationFactor_lumi", &EvtInfo_NormalizationFactor_lumi, "EvtInfo_NormalizationFactor_lumi/F");
    mytree -> Branch("EvtInfo_NPu", &EvtInfo_NPu, "EvtInfo_NPu/I");
    mytree -> Branch("EvtInfo_Rho", &EvtInfo_Rho, "EvtInfo_Rho/F");
    mytree -> Branch("EvtInfo_NVtx", &EvtInfo_NVtx, "EvtInfo_NVtx/I");
    mytree -> Branch("EvtInfo_genweight", &EvtInfo_genweight, "EvtInfo_genweight/F");
    //------------------------
    mytree -> Branch("DiPhoInfo_mass", &DiPhoInfo_mass, "DiPhoInfo_mass/F");
    mytree -> Branch("DiPhoInfo_pt", &DiPhoInfo_pt, "DiPhoInfo_pt/F");
    mytree -> Branch("DiPhoInfo_eta", &DiPhoInfo_eta, "DiPhoInfo_eta/F");
    mytree -> Branch("DiPhoInfo_phi", &DiPhoInfo_phi, "DiPhoInfo_phi/F");
    mytree -> Branch("DiPhoInfo_energy", &DiPhoInfo_energy, "DiPhoInfo_energy/F");
    mytree -> Branch("DiPhoInfo_leadPt", &DiPhoInfo_leadPt, "DiPhoInfo_leadPt/F");
    mytree -> Branch("DiPhoInfo_leadEta", &DiPhoInfo_leadEta, "DiPhoInfo_leadEta/F");
    mytree -> Branch("DiPhoInfo_leadPhi", &DiPhoInfo_leadPhi, "DiPhoInfo_leadPhi/F");
    mytree -> Branch("DiPhoInfo_leadE", &DiPhoInfo_leadE, "DiPhoInfo_leadE/F");
    mytree -> Branch("DiPhoInfo_leadhoe", &DiPhoInfo_leadhoe, "DiPhoInfo_leadhoe/F");
    mytree -> Branch("DiPhoInfo_leadIDMVA", &DiPhoInfo_leadIDMVA, "DiPhoInfo_leadIDMVA/F");
    mytree -> Branch("DiPhoInfo_subleadPt", &DiPhoInfo_subleadPt, "DiPhoInfo_subleadPt/F");
    mytree -> Branch("DiPhoInfo_subleadEta", &DiPhoInfo_subleadEta, "DiPhoInfo_subleadEta/F");
    mytree -> Branch("DiPhoInfo_subleadPhi", &DiPhoInfo_subleadPhi, "DiPhoInfo_subleadPhi/F");
    mytree -> Branch("DiPhoInfo_subleadE", &DiPhoInfo_subleadE, "DiPhoInfo_subleadE/F");
    mytree -> Branch("DiPhoInfo_subleadhoe", &DiPhoInfo_subleadhoe, "DiPhoInfo_subleadhoe/F");
    mytree -> Branch("DiPhoInfo_subleadIDMVA", &DiPhoInfo_subleadIDMVA, "DiPhoInfo_subleadIDMVA/F");
    //------------------------
    mytree -> Branch("ElecInfo_Size", &ElecInfo_Size, "ElecInfo_Size/I");
    mytree -> Branch("MuonInfo_Size", &MuonInfo_Size, "MuonInfo_Size/I");
    mytree -> Branch("num_leptons", &num_leptons, "num_leptons/I");// # of selected objects.
    mytree -> Branch("num_electrons", &num_electrons, "num_electrons/I");// # of selected objects.
    mytree -> Branch("num_muons", &num_muons, "num_muons/I");// # of selected objects.
    mytree -> Branch("ElecInfo_electron_pt", &ElecInfo_electron_pt);
    mytree -> Branch("ElecInfo_electron_eta", &ElecInfo_electron_eta);
    mytree -> Branch("ElecInfo_electron_phi", &ElecInfo_electron_phi);
    mytree -> Branch("ElecInfo_electron_energy", &ElecInfo_electron_energy);
    mytree -> Branch("ElecInfo_electron_diphoton_deltaR", &ElecInfo_electron_diphoton_deltaR);
    mytree -> Branch("ElecInfo_electron_leadingPhoton_deltaR", &ElecInfo_electron_leadingPhoton_deltaR);
    mytree -> Branch("ElecInfo_electron_subleadingPhoton_deltaR", &ElecInfo_electron_subleadingPhoton_deltaR);
    mytree -> Branch("MuonInfo_muon_pt", &MuonInfo_muon_pt);
    mytree -> Branch("MuonInfo_muon_eta", &MuonInfo_muon_eta);
    mytree -> Branch("MuonInfo_muon_phi", &MuonInfo_muon_phi);
    mytree -> Branch("MuonInfo_muon_energy", &MuonInfo_muon_energy);
    mytree -> Branch("MuonInfo_muon_diphoton_deltaR", &MuonInfo_muon_diphoton_deltaR);
    mytree -> Branch("MuonInfo_muon_leadingPhoton_deltaR", &MuonInfo_muon_leadingPhoton_deltaR);
    mytree -> Branch("MuonInfo_muon_subleadingPhoton_deltaR", &MuonInfo_muon_subleadingPhoton_deltaR);
    //------------------------
    mytree -> Branch("jets_size", &jets_size, "jets_size/I");
    mytree -> Branch("num_jets", &num_jets, "num_jets/I");
    mytree -> Branch("JetInfo_jet_pt", &JetInfo_jet_pt);
    mytree -> Branch("JetInfo_jet_eta", &JetInfo_jet_eta);
    mytree -> Branch("JetInfo_jet_phi", &JetInfo_jet_phi);
    mytree -> Branch("JetInfo_jet_energy", &JetInfo_jet_energy);
    mytree -> Branch("JetInfo_jet_diphoton_deltaR", &JetInfo_jet_diphoton_deltaR);
    mytree -> Branch("JetInfo_jet_leadingPhoton_deltaR", &JetInfo_jet_leadingPhoton_deltaR);
    mytree -> Branch("JetInfo_jet_subleadingPhoton_deltaR", &JetInfo_jet_subleadingPhoton_deltaR);
    mytree -> Branch("JetInfo_jet_pfDeepCSVJetTags_probb", &JetInfo_jet_pfDeepCSVJetTags_probb);
    mytree -> Branch("JetInfo_jet_pfDeepCSVJetTags_probbb", &JetInfo_jet_pfDeepCSVJetTags_probbb);
    //mytree -> Branch("num_btagged_jets", &num_btagged_jets, "num_btagged_jets/I");
    //mytree -> Branch("num_nonbtagged_jets", &num_nonbtagged_jets, "num_nonbtagged_jets/I");
    //------------------------
    //mytree -> Branch("inv_mass_dijet", &inv_mass_dijet, "inv_mass_dijet/F");
    //mytree -> Branch("inv_mass_diphoton", &inv_mass_diphoton, "inv_mass_diphoton/F");
    //mytree -> Branch("inv_mass_tbw", &inv_mass_tbw, "inv_mass_tbw/F");
}
void myTreeClass::Fill(){
    mytree -> Fill();
}
void myParameters::Clear(){
    EvtInfo_totalEntry_before_preselection = 0;
    EvtInfo_NormalizationFactor_lumi = 0;
    //------------------------
    DiPhoInfo_eta = 0;
    DiPhoInfo_phi = 0;
    DiPhoInfo_energy = 0;
    //------------------------
    num_jets = 0;
    JetInfo_jet_pt.clear();
    JetInfo_jet_eta.clear();
    JetInfo_jet_phi.clear();
    JetInfo_jet_energy.clear();
    JetInfo_jet_diphoton_deltaR.clear();
    JetInfo_jet_leadingPhoton_deltaR.clear();
    JetInfo_jet_subleadingPhoton_deltaR.clear();
    JetInfo_jet_pfDeepCSVJetTags_probb.clear();
    JetInfo_jet_pfDeepCSVJetTags_probbb.clear();
    //------------------------
    num_leptons = 0;
    num_electrons = 0;
    num_muons = 0;
    ElecInfo_electron_pt.clear();
    ElecInfo_electron_eta.clear();
    ElecInfo_electron_phi.clear();
    ElecInfo_electron_energy.clear();
    ElecInfo_electron_diphoton_deltaR.clear();
    ElecInfo_electron_leadingPhoton_deltaR.clear();
    ElecInfo_electron_subleadingPhoton_deltaR.clear();
    MuonInfo_muon_pt.clear();
    MuonInfo_muon_eta.clear();
    MuonInfo_muon_phi.clear();
    MuonInfo_muon_energy.clear();
    MuonInfo_muon_diphoton_deltaR.clear();
    MuonInfo_muon_leadingPhoton_deltaR.clear();
    MuonInfo_muon_subleadingPhoton_deltaR.clear();
    //------------------------
    //Not used in preselection stage
    //------------------------
    num_btagged_jets = 0;
    num_nonbtagged_jets = 0;
    //------------------------
    //Chi-2 sorting related
    //------------------------
    Mass_w_candidate_chi2_simple = 0;
    Mass_top_candidate_chi2_simple = 0;
    Mass_w_candidate_chi2_modified = 0;
    Mass_top_candidate_chi2_modified = 0;
    Mass_gen_w_candidate_chi2_simple = 0;
    Mass_gen_top_candidate_chi2_simple = 0;
    Mass_gen_w_candidate_chi2_modified = 0;
    Mass_gen_top_candidate_chi2_modified = 0;

    inv_mass_dijet = 0;
    inv_mass_diphoton = 0;
    inv_mass_tbw = 0;
    //------------------------
    JetInfo_dijet_delta_eta = 0;
    JetInfo_dijet_delta_phi = 0;
    JetInfo_dijet_delta_angle = 0;
    //------------------------
}
//### }}}

//=== leptonic ===//
// print and mkplots{{{
void PrintCountsAndRatio(const char* title, int a, int b){
    double n = (double) b;
    double p = (double) a / (double) b;
    double q = 1 - p;
    double error = sqrt(p * q / n);
    printf("[INFO] %-50s = %6d / %6d (%6.2f%% +/- %6.2f%%)\n", title, a, b, 100. * (double)a / (double)b, 100. * error);
}
void plots_comparison_two_methods(const char* dirName, const char* plotName, TCanvas *c1, TLegend *legend, TH1D* hist_quadratic, TH1D* hist_topKinFit)
{
    double scale = 1.2;
    double max_01 = hist_quadratic -> GetMaximum();
    double max_02 = hist_topKinFit -> GetMaximum();
    double max = (max_01 > max_02) ? max_01 : max_02;

    hist_quadratic -> SetStats(0);
    hist_quadratic -> SetMaximum(max*scale);
    hist_quadratic->Draw();

    hist_topKinFit -> SetStats(0);
    hist_topKinFit->SetLineColor(kRed);
    hist_topKinFit->Draw("same");

    legend->Clear();
    legend->AddEntry(hist_quadratic,  "quardratic", "l");
    legend->AddEntry(hist_topKinFit,  "topKinFit", "l");
    legend->SetLineColor(0);
    legend->Draw("same");
    c1->SaveAs( Form("%s/%s.png", dirName, plotName) )    ;
}
void MakePlots_coeffD(TCanvas* c1, TLegend* legend_ratio, const char* histName, const char* label1, const char* label2, TH1D* hist, TH1D* hist_positiveD, TH1D* hist_negativeD){
    hist -> Draw("hist");
    hist -> SetLineWidth(2);
    hist_positiveD -> Draw("hist;same");
    hist_positiveD -> SetLineColor(kRed);
    hist_positiveD -> SetLineWidth(2);
    hist_positiveD -> SetLineStyle(2);
    hist_negativeD -> Draw("hist;same");
    hist_negativeD -> SetLineColor(kGreen);
    hist_negativeD -> SetLineWidth(2);
    legend_ratio->Clear();
    legend_ratio->AddEntry(hist, "All", "l");
    legend_ratio->AddEntry(hist_positiveD, label1, "l");
    legend_ratio->AddEntry(hist_negativeD, label2, "l");
    legend_ratio->SetLineColor(0);
    legend_ratio->Draw("same");
    c1->SaveAs(Form("ntuples_skimmed/%s.png", histName));
}
void Set2DPlot(TH2D *&h){
    h->GetXaxis()->SetTitleSize(40);//25
    h->GetXaxis()->SetTitleFont(43);
    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetTitleSize(32);//20
    h->GetYaxis()->SetTitleFont(43);
    h->GetYaxis()->SetTitleOffset(0.9);//1.2
}
///}}}
/*
// ### Alternative Mother info{{{
bool isMatched_with_Gen_tbw(std::vector<int> *GenPartInfo_MomPdgID, int index_bjet, int index_jet1, int index_jet2){
    bool isCorrect = \
    index_bjet > 0 &&\
    index_jet1 > 0 &&\
    index_jet2 > 0 &&\
    abs(GenPartInfo_MomPdgID->at(index_bjet)) == 6 &&\
    abs(GenPartInfo_MomPdgID->at(index_jet1)) == 24 &&\
    abs(GenPartInfo_MomPdgID->at(index_jet2)) == 24 &&\
    (index_jet1 != index_jet2);

    return isCorrect;
}
//old code{{{
//bool isMatched_with_Gen_W_Boson(TLorentzVector gen_w_sel, TH1D *&hist, Int_t GenPartInfo_size, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_PdgID){
//    double delta_R_gen_w = 999;
//    for(int j=0; j<GenPartInfo_size; j++){
//        if( abs(GenPartInfo_PdgID->at(j)) == 24 ){
//            TLorentzVector genParticle;
//            genParticle.SetPtEtaPhiM(GenPartInfo_Pt->at(j), GenPartInfo_Eta->at(j), GenPartInfo_Phi->at(j), GenPartInfo_Mass->at(j));
//            delta_R_gen_w = gen_w_sel.DeltaR(genParticle);
//        }
//    }//end of gen loop
//    //--- remove event with bad combination ---//
//    //if(gen_w_sel.M() < 20) continue; 
//    hist->Fill(delta_R_gen_w);
//    if(delta_R_gen_w > 0.005 || gen_w_sel.M() < 20) return false;
//    else return true;
//}
//}}}
//}}}
// ##### print gen info{{{
void print_matched_gen_info(TLorentzVector jet, Int_t GenPartInfo_size, std::vector<int> *GenPartInfo_MomPdgID, std::vector<float> *GenPartInfo_Pt, std::vector<float> *GenPartInfo_Eta, std::vector<float> *GenPartInfo_Phi, std::vector<float> *GenPartInfo_Mass, std::vector<int> *GenPartInfo_Status, std::vector<int> *GenPartInfo_PdgID){
    TLorentzVector truelove;//we are looking for the right genParticle to match the jet.
    int index = -999, truelove_PdgID = -999; double delta_R_min = 999.;
    for(int i=0; i<GenPartInfo_size; i++){
        if( abs(GenPartInfo_Status->at(i)) != 23 ) continue;//remove incoming/intermediate particles
        if( abs(GenPartInfo_PdgID->at(i)) > 6 ) continue;//exclude top quark & other particles
        //--------------------
        TLorentzVector genParticle;
        genParticle.SetPtEtaPhiM(GenPartInfo_Pt->at(i), GenPartInfo_Eta->at(i), GenPartInfo_Phi->at(i), GenPartInfo_Mass->at(i));
        double delta_R = jet.DeltaR(genParticle);
        //select quark & require min(deltaR)
        if( delta_R < 0.4 && delta_R < delta_R_min){
            index = i;//record the matched genParticle
            delta_R_min = delta_R;
            truelove = genParticle;
            truelove_PdgID = GenPartInfo_PdgID->at(i);
        }
    }//end of gen loop

    if(index>0){
        printf("[check-tqh] \n");
            printf("Status = %3d, ", GenPartInfo_Status->at(index));
            printf("PdgID = %3d, ", GenPartInfo_PdgID->at(index));
            printf("Pt = %6.2f, ", GenPartInfo_Pt->at(index));
            printf("Eta = %9.2f, ", GenPartInfo_Eta->at(index));
            printf("Phi = %6.2f, ", GenPartInfo_Phi->at(index));
            printf("Mass = %6.2f, ", GenPartInfo_Mass->at(index));
            printf("MomPdgID = %5d\n", GenPartInfo_MomPdgID->at(index));
    } else{
        printf("[check-tqh] Not found.\n");
    }
}
//}}}
//### report{{{
//void kinematics_info(const char* Title, TLorentzVector Particle, int index){
//        printf("(%s) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f, Energy = %6.2f, Mass = %6.2f\n", Title, Particle.Pt(), Particle.Eta(), Particle.Phi(), Particle.Energy(), Particle.M());
//}
void kinematics_report(const char* recoTitle, TLorentzVector recoParticle, int id_recoParticle, const char* genTitle, TLorentzVector genParticle, int genParticle_PdgID){
        double delta_R = genParticle.DeltaR(recoParticle);
        printf("(%s) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f, Energy = %6.2f, Mass = %6.2f, id = %d\n", recoTitle, recoParticle.Pt(), recoParticle.Eta(), recoParticle.Phi(), recoParticle.Energy(), recoParticle.M(), id_recoParticle);
        printf("(%s) Pt = %6.2f, Eta = %6.2f, Phi = %6.2f, Energy = %6.2f, Mass = %6.2f, id = %d, delta_R = %6.2f\n", genTitle, genParticle.Pt(), genParticle.Eta(), genParticle.Phi(), genParticle.Energy(), genParticle.M(), genParticle_PdgID, delta_R);
}
void hist_bin_fraction(TH1D *hist, const char* title, int entries_matched_in_bin2){
    int total_entries = hist->GetEntries();
    printf("[INFO-hist] %s (%d)\n", title, total_entries);
    printf("[INFO-hist] fractions = ");
    for(int i=0; i<3; ++i){
        if(i<2) printf("%6.4f, ", (double) hist->GetBinContent(i+1) / (double)total_entries);
        else{
            int accumulation = 0;
            for(int j=i; j<10; ++j) accumulation += (double) hist->GetBinContent(j+1);
            printf("%6.4f; ", i, (double) accumulation / (double)total_entries);
            printf("%6.4f, ", i, (double) entries_matched_in_bin2 / (double)total_entries);
            printf("%6.4f\n", i, (double) entries_matched_in_bin2 / (double)hist->GetBinContent(2));
        }
    }
}
void hist_report(TH1D *hist, const char* chi2_type){
    double mean = hist->GetMean();
    double sigma = hist->GetMeanError();
    double sigma2 = hist->GetRMS();
    int bin_maximum = hist->GetMaximumBin();
    double maxbin_lower_edge = hist->GetBinLowEdge(bin_maximum);
    double bin_width = hist->GetBinWidth(bin_maximum);
    double maxbin_upper_edge = maxbin_lower_edge + bin_width;

    //printf("[INFO] chi2 %s: mean = %6.2f, sigma = %6.2f, %6.2f\n", chi2_type, mean, sigma, sigma2);
    printf("[INFO] chi2 %s: mean = %6.2f, sigma = %6.2f, mode lies in [%6.2f, %6.2f]\n", chi2_type, mean, sigma2, maxbin_lower_edge, maxbin_upper_edge);
}
//}}}
//### mkplots{{{
void MakeTwoPlots(TCanvas *c1, TH1D* hist_gen, TH1D* hist_reco, TLegend *legend, const char* name){
    hist_reco->SetStats(0);
    hist_reco->SetLineWidth(2);
    hist_reco->SetLineColor(kBlue);
    hist_reco->Draw();
    hist_gen->SetLineWidth(2);
    hist_gen->SetLineColor(kRed);
    hist_gen->Draw("same;hist");

    legend->Clear();
    legend->AddEntry(hist_gen, "gen-level", "l");
    legend->AddEntry(hist_reco, "reco", "l");
    legend->SetLineColor(0);
    legend->Draw("same");

    c1->SaveAs(Form("ntuples_skimmed/%s", name));
}
void MakePlots(TCanvas *c1, TH1D* hist, const char* title, const char* outputFile){
    hist->Draw();
    hist->SetTitle(title);
    hist->SetXTitle(title);
    hist->SetYTitle("Entries");
    hist->GetYaxis()->SetTitleOffset(1.4);
    hist->Write();
    c1->SaveAs(outputFile);
}
void MakeFinalPlots(TCanvas *c1, TH1D* hist_simple, TH1D* hist_modified, TH1D*hist_yfyj, TLegend *legend, const char* name){
    double max;
    double scale = 1.2;
    double max_simple = hist_simple -> GetMaximum();
    double max_modified = hist_modified -> GetMaximum();
    double max_yfyj = hist_yfyj -> GetMaximum();
    max = (max_modified > max_simple) ? max_modified : max_simple;
    max = (max_yfyj > max_modified) ? max_yfyj : max_modified;
    double overflow_simple = hist_simple->GetBinContent(hist_simple->GetNbinsX() + 1);
    double overflow_modified = hist_modified->GetBinContent(hist_modified->GetNbinsX() + 1);
    double overflow_yfyj = hist_yfyj->GetBinContent(hist_yfyj->GetNbinsX() + 1);
    //max = (max > overflow_simple) ? max : overflow_simple;
    //max = (max > overflow_modified) ? max : overflow_modified;
    //max = (max > overflow_yfyj) ? max : overflow_yfyj;

    hist_simple->SetStats(0);
    hist_simple->SetMaximum(max*scale);
    hist_simple->SetLineWidth(2);
    hist_simple->SetLineColor(kBlue);
    hist_simple->GetXaxis()->SetRange(1, hist_simple->GetNbinsX() + 1);
    hist_simple->Draw();
    hist_modified->SetLineWidth(2);
    hist_modified->SetLineColor(kRed);
    hist_modified->GetXaxis()->SetRange(1, hist_modified->GetNbinsX() + 1);
    hist_modified->Draw("same");
    hist_yfyj->SetLineWidth(2);
    hist_yfyj->SetLineColor(kGreen+4);
    hist_yfyj->GetXaxis()->SetRange(1, hist_yfyj->GetNbinsX() + 1);
    hist_yfyj->Draw("same");

    legend->Clear();
    legend->SetTextSize(0.03);
    legend->AddEntry(hist_yfyj, "leading-jets method", "l");
    legend->AddEntry(hist_simple, "simple #chi^{2} method", "l");
    legend->AddEntry(hist_modified, "modified #chi^{2} method", "l");
    legend->SetLineColor(0);
    legend->Draw("same");

    c1->SaveAs(Form("ntuples_skimmed/%s", name));
}
//###}}}
//### bool functions{{{
bool isThisMCsignal(char* dataset){
    if((string)dataset == "ST_FCNC-TH_Thadronic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "ST_FCNC-TH_Thadronic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "ST_FCNC-TH_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "ST_FCNC-TH_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-T2HJ_aTleptonic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-TtoHJ_aThadronic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-TtoHJ_aTleptonic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-aTtoHJ_Thadronic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hct-MadGraph5-pythia8") return true;
    if((string)dataset == "TT_FCNC-aTtoHJ_Tleptonic_HToaa_eta_hut-MadGraph5-pythia8") return true;
    return false;
}
bool isThisDataOrNot(char* dataset){
    if((string)dataset == "DoubleEG_B") return true;
    if((string)dataset == "DoubleEG_C") return true;
    if((string)dataset == "DoubleEG_D") return true;
    if((string)dataset == "DoubleEG_E") return true;
    if((string)dataset == "DoubleEG_F") return true;
    return false;
}
bool isThisMultiFile(char* dataset){
    if((string)dataset == "DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa") return true;
    if((string)dataset == "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8") return true;
    return false;
}
//###}}}
*/
