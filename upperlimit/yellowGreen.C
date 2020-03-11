{
    double signal_strength_25  = 0.0032;
    double signal_strength_160 = 0.0051;
    double signal_strength_500 = 0.0093;
    double signal_strength_840 = 0.0184;
    double signal_strength_975 = 0.0279;

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    TH1D* hist = new TH1D("hist", ";;signal strength", 2, 0, 2);
    hist->SetStats(0);
    hist->SetMaximum(signal_strength_975*1.5);
    hist->GetXaxis()->SetBinLabel(1, "Leptonic");
    hist->GetXaxis()->SetBinLabel(2, "Hadronic");
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.08);
    hist->GetXaxis()->SetTitleOffset(0.5);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetTitleOffset(1.4);
    hist->SetLineWidth(2);
    hist->SetLineStyle(2);
    hist->SetBinContent(1, signal_strength_500);
    hist->Draw();
    //--------------------------------------------------
    TH1D* h_yellow = new TH1D("h_yellow", "", 2, 0, 2);
    double mean_yellow  = (signal_strength_975 + signal_strength_25) / 2. ;
    double error_yellow = (signal_strength_975 - signal_strength_25) / 2. ;
    h_yellow->SetBinContent(1, mean_yellow);
    h_yellow->SetBinError(1, error_yellow);
    h_yellow->SetFillColor(kYellow);
    h_yellow->Draw("E2,same");
    //--------------------------------------------------
    TH1D* h_green = new TH1D("h_green", "", 2, 0, 2);
    double mean_green  = (signal_strength_840 + signal_strength_160) / 2. ;
    double error_green = (signal_strength_840 - signal_strength_160) / 2. ;
    h_green->SetBinContent(1, mean_green);
    h_green->SetBinError(1, error_green);
    h_green->SetFillColor(kGreen);
    h_green->Draw("E2,same");
    //--------------------------------------------------
    TLegend *legend = new TLegend(0.68,0.60,0.88,0.88);
    legend->SetLineColor(0);
    legend->SetTextFont(42);
    legend->AddEntry(hist,     "Expected limit", "l");
    legend->AddEntry(h_green,  "1#sigma stat. error", "f");
    legend->AddEntry(h_yellow, "2#sigma stat. error", "f");
    hist->Draw("same");
    legend->Draw("same");
    //--------------------------------------------------
    TLatex latex;
    latex.SetNDC(kTRUE);
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.SetTextAlign(11);
    latex.DrawLatex(0.10, 0.91, "#bf{CMS} #it{Simulation}");
    //--------------------------------------------------
    gPad->SetLeftMargin(0.10);
    gPad->SetTicks(0,1);
    c1->SaveAs("plot_upperlimit_hct.png");
}
