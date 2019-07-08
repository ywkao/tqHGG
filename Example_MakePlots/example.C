#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPad.h>

void example(){
    TCanvas *c1 = new TCanvas("c1", "c1", 700, 800);
    //--------------------
    TFile *file_data = TFile::Open("data.root");
    TFile *file_mc  = TFile::Open("mc.root");
    //--------------------
    TH1D *hist_data  = (TH1D*)file_data->Get("hist_DiPhoInfo_mass");
    TH1D *hist_mc    = (TH1D*)file_mc->Get("hist_DiPhoInfo_mass");   hist_mc->Scale(120.);  //Artificial scalling
    TH1D *hist_err   = (TH1D*)hist_mc->Clone();
    TH1D *hist_ratio = (TH1D*)hist_data->Clone();//will be ploted in lower pad
    TH1D *hist_ratio_mc_uncertainty = (TH1D*)hist_mc->Clone();//will be ploted in lower pad
    hist_ratio->Divide(hist_mc);
    //--------------------
    hist_ratio->SetTitle("");
    hist_ratio->SetStats(0); //No statistics on lower plot
    hist_ratio->SetMarkerStyle(20);
    hist_ratio->SetMarkerSize(1.2);
    //--------------------
    hist_data->SetTitle("");
    hist_data->SetStats(0); //No statistics on lower plot
    hist_data->SetMarkerStyle(20);
    hist_data->SetMarkerSize(1.2);
    //--------------------
    hist_mc->SetTitle("");
    hist_mc->SetStats(0); //No statistics on lower plot
    hist_mc->SetFillColor(kBlue);
    hist_mc->SetLineColor(kBlue);
    //--------------------
    hist_err->SetMarkerColor(kGray+2);
    hist_err->SetFillColor(kGray+2);
    hist_err->SetLineColor(kGray+2);
    hist_err->SetFillStyle(3001);


    //=============================================//
    //===== Calculate relative MC uncertainty =====//
    //=============================================//
    int nbins = hist_mc->GetNbinsX();
    //printf("nbins=%d\n", nbins);
    for(int i=0; i<nbins; i++){
        double mean = hist_mc->GetBinContent(i+1);
        double error = hist_mc->GetBinError(i+1);
        double upper_rel_error = (mean==0.) ? 0. : (mean+error)/mean;
        double lower_rel_error = (mean==0.) ? 0. : (mean-error)/mean;
        double error_mean = (upper_rel_error + lower_rel_error)/2.;
        double error_error = (upper_rel_error - lower_rel_error)/2.;
        hist_ratio_mc_uncertainty->SetBinContent(i+1, error_mean);
        hist_ratio_mc_uncertainty->SetBinError(i+1, error_error);
        //printf("bin = %d, mean=%f, error=%f, error_error=%f\n", i+1, mean, error, error_error);
    }



    //============================//
    //===== Draw upper plots =====//
    //============================//
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
    pad1->SetBottomMargin(0); //Upper and lower pads are joined
    pad1->SetAttLinePS(kBlack,1,2);
    pad1->Draw();
    pad1->cd(); //pad1 becomes current pad
    gPad->SetTicks(1,1);
    //--------------------
    //Determine which hist to plot first, and set scale automatically.
    //--------------------
    double scale = 2.5;
    double max_data = hist_data->GetMaximum();
    double max_mc  = hist_mc->GetMaximum();
    if(max_mc > max_data){
        hist_mc->Draw("hist");//histogram
        hist_err->Draw("E2,same");//plot the error of mc
        hist_data->Draw("p,E1,same");//point, error bar, same pad
        hist_mc->SetMaximum(max_mc*scale);
    } else{
        hist_data->Draw("p,E1");//point, error bar
        hist_err->Draw("E2,same");//plot the error of mc
        hist_mc->Draw("hist,same");//histogram
        hist_data->Draw("p,E1,same");//Make data points over the histogram
        hist_data->SetMaximum(max_data*scale);
    } 
    //--------------------
    //TLegend setting
    //--------------------
    TLegend *legend = new TLegend(0.50,0.55,0.85,0.85);
    legend->AddEntry(hist_data, "Observed", "lep");
    legend->AddEntry(hist_mc, "Background", "f");
    legend->AddEntry(hist_err, "Uncertainty", "f");
    legend->SetLineColor(0);
    legend->Draw("same");
    //--------------------
    //Add some texts
    //--------------------
    TLatex latex, latex_lumi;
    latex.SetNDC(kTRUE);
    latex.SetTextFont(43);
    latex.SetTextSize(22);
    latex.SetTextAlign(13);
    latex.DrawLatex(0.14, 0.84, "#bf{CMS}");
    latex.DrawLatex(0.14, 0.78, "#it{Preliminary}");//NOTE: change to "Open Data""
    //--------------------
    latex_lumi.SetNDC(kTRUE);
    latex_lumi.SetTextFont(43);
    latex_lumi.SetTextSize(22);
    latex_lumi.SetTextAlign(31);
    latex_lumi.DrawLatex(0.89, 0.92, "41.5 fb^{-1} (2017, 13 TeV)");//NOTE: change to "?? fb^{-1} (2012, 8 TeV)"

    //============================//
    //===== Draw lower plots =====//
    //============================//
    c1->cd(); //Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
    pad2->SetTopMargin(0.01); 
    pad2->SetBottomMargin(0.4);
    pad2->SetAttLinePS(kBlack,1,2);
    pad2->SetGridx(1); //Upper and lower pads are joined
    pad2->Draw();
    pad2->cd(); //pad2 becomes current pad
    gPad->SetTicks(1,1);
    //--------------------
    hist_ratio->SetMaximum(1.75);
    hist_ratio->SetMinimum(0.5);
    hist_ratio->Draw("p,E1");
    hist_ratio_mc_uncertainty->SetFillColor(kGray);
    hist_ratio_mc_uncertainty->SetLineColor(kGray);
    hist_ratio_mc_uncertainty->Draw("E2,same");
    hist_ratio->Draw("p,E1,same");
    //--------------------
    hist_ratio->GetYaxis()->SetTitle("Obs/Exp");
    hist_ratio->GetYaxis()->SetNdivisions(5);
    hist_ratio->GetYaxis()->SetTitleSize(20);
    hist_ratio->GetYaxis()->SetTitleFont(43);
    hist_ratio->GetYaxis()->SetTitleOffset(1.2);
    hist_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hist_ratio->GetYaxis()->SetLabelSize(15);
    //--------------------
    hist_ratio->GetXaxis()->SetTitle("Diphoton Invariant Mass [GeV/c^2]");
    hist_ratio->GetXaxis()->SetTitleSize(20);
    hist_ratio->GetXaxis()->SetTitleFont(43);
    hist_ratio->GetXaxis()->SetTitleOffset(4.);
    hist_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hist_ratio->GetXaxis()->SetLabelSize(15);
    //--------------------
    c1->Update();//update the value of pad2->GetUxmax().
    TLine line;
    line.SetLineStyle(2);
    line.DrawLine(pad2->GetUxmin(),0.5,pad2->GetUxmax(),0.5);
    line.DrawLine(pad2->GetUxmin(),1.0,pad2->GetUxmax(),1.0);
    line.DrawLine(pad2->GetUxmin(),1.5,pad2->GetUxmax(),1.5);
    line.DrawLine(pad2->GetUxmin(),2.0,pad2->GetUxmax(),2.0);
    //--------------------
    c1->SaveAs("example.png");
}
