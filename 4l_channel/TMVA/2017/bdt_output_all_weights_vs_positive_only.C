void bdt_outut_all_weights_vs_positive_only()
{
    gStyle->SetOptStat(0);

    TFile *f2017_all = new TFile("./all/TMVA.root");
    TFile *f2017_positive = new TFile("./positive_only/TMVA.root");

    f2017_all->cd("dataset/Method_BDT/BDT");
    TH1D *hBDT_all = (TH1D*)gDirectory->Get("MVA_BDT_S");
    hBDT_all->SetLineColor(kRed);

    f2017_positive->cd("dataset/Method_BDT/BDT");
    TH1D *hBDT_positive = (TH1D*)gDirectory->Get("MVA_BDT_S");
    hBDT_positive->SetLineColor(kGreen);

    TCanvas *c1 = new TCanvas();
    c1->cd();
    hBDT_all->SetTitle("2017 (signal)");
    hBDT_all->GetYaxis()->SetRangeUser(0,4);
    hBDT_all->Draw();
    hBDT_positive->Draw("same");

    TLegend *legend = new TLegend(0.15,0.65,0.25,0.88);
	legend->SetTextFont(132);
	legend->SetTextSize(0.05);
	legend->AddEntry(hBDT_all,"all weights","l");
    legend->AddEntry(hBDT_positive,"positive weights only","l");
	legend->SetBorderSize(0);
    legend->Draw();

    //--------------------------------------------------------

    f2017_all->cd("dataset/Method_BDT/BDT");
    hBDT_all = (TH1D*)gDirectory->Get("MVA_BDT_B");
    hBDT_all->SetLineColor(kRed);

    f2017_positive->cd("dataset/Method_BDT/BDT");
    hBDT_positive = (TH1D*)gDirectory->Get("MVA_BDT_B");
    hBDT_positive->SetLineColor(kGreen);

    TCanvas *c2 = new TCanvas();
    c2->cd();
    hBDT_all->SetTitle("2017 (background)");
    hBDT_all->GetYaxis()->SetRangeUser(0,10);
    hBDT_all->Draw();
    hBDT_positive->Draw("same");

    legend = new TLegend(0.15,0.65,0.25,0.88);
	legend->SetTextFont(132);
	legend->SetTextSize(0.05);
	legend->AddEntry(hBDT_all,"all weights","l");
    legend->AddEntry(hBDT_positive,"positive weights only","l");
	legend->SetBorderSize(0);
    legend->Draw();
}