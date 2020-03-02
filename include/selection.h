#ifndef __SELECTION_H__
#define __SELECTION_H__

void Selection(char* input_file, char* output_file, char* output_tree, char* dataset, char* output_dir, char* channel);
TLorentzVector GetBestM1(float &M1, int num_jets, int index_bjet, std::vector<int> index_jet, TLorentzVector diphoton, std::vector<TLorentzVector> Jets, int &index_q, TLorentzVector &jet_q);
void MakePlots(TCanvas *c1, TH1D* hist, const char* outputFile);
float Chi2_calculator_modified(float w_mass, float t_mass);
#endif
