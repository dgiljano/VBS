void roc_all_weights_vs_positive_only()
{
    gStyle->SetOptStat(0);

    TFile *f2017_all = new TFile("./all/TMVA.root");
    TFile *f2017_positive = new TFile("./positive_only/TMVA.root");

    f2017_all->cd("dataset/Method_BDT/BDT");
    TH1D *hROC_all = (TH1D*)gDirectory->Get("MVA_BDT_rejBvsS");
    hROC_all->SetLineColor(kRed);

    f2017_positive->cd("dataset/Method_BDT/BDT");
    TH1D *hROC_positive = (TH1D*)gDirectory->Get("MVA_BDT_rejBvsS");
    hROC_positive->SetLineColor(kGreen);

    hROC_all->SetTitle("2017");
    hROC_all->Draw();
    hROC_positive->Draw("same");

    TLegend *legend = new TLegend(0.15,0.45,0.25,0.68);
	legend->SetTextFont(132);
	legend->SetTextSize(0.05);
	legend->AddEntry(hROC_all,"all weights","l");
    legend->AddEntry(hROC_positive,"positive weights only","l");
	legend->SetBorderSize(0);
    legend->Draw();
}