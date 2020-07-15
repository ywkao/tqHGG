#ifndef _generalChiSquareStudy_leptonic_h_
#define _generalChiSquareStudy_leptonic_h_

#include <vector>
#include <TH1D.h>
#include <TLorentzVector.h>
using namespace std;

void initHists();
double evaluate_neutrino_pz(TLorentzVector lepton, vector<double> met_info);
TLorentzVector derive_reco_neutrino(TLorentzVector lepton, vector<double> met_info);
TLorentzVector derive_reco_wboson(TLorentzVector lepton, TLorentzVector reco_neutrino);
TLorentzVector derive_reco_tbw(TLorentzVector reco_wboson, TLorentzVector bjet);
int get_q_index_min_chi2(std::vector<TLorentzVector> Jets, int index_bjet, TLorentzVector diphoton);
bool isMatched(int a, int b);

#endif
