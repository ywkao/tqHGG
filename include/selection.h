#ifndef __SELECTION_H__
#define __SELECTION_H__

void Selection(char* input_file, char* output_file, char* dataset, char* output_dir, char* channel);
void MakePlots(TCanvas *c1, TH1D* hist, const char* outputFile);
#endif
