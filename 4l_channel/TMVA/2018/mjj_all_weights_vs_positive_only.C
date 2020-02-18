void mjj_all_weights_vs_positive_only()
{
    gStyle->SetOptStat(0);

    bool normalise = true;

    TFile *f2018_all = new TFile("./all/TMVA.root");
    TFile *f2018_positive = new TFile("./positive_only/TMVA.root");

    f2018_all->cd("dataset/Method_BDT/BDT");
    TH1D *hmjj_all = (TH1D*)gDirectory->Get("mjj__Signal");
    hmjj_all->SetLineColor(kRed);

    f2018_positive->cd("dataset/Method_BDT/BDT");
    TH1D *hmjj_positive = (TH1D*)gDirectory->Get("mjj__Signal");
    hmjj_positive->SetLineColor(kGreen);

    TCanvas *c1 = new TCanvas();
    c1->cd();
    hmjj_all->SetTitle("mjj 2018 (signal)");
    hmjj_all->GetYaxis()->SetRangeUser(0,0.5);
    if (normalise)
    {
        hmjj_all->Scale(1/hmjj_all->Integral());
        hmjj_positive->Scale(1/hmjj_positive->Integral());
    }
    hmjj_all->GetYaxis()->SetRangeUser(0,0.12);
    hmjj_all->Draw();
    hmjj_positive->Draw("same");

    TLegend *legend = new TLegend(0.15,0.65,0.25,0.88);
	legend->SetTextFont(132);
	legend->SetTextSize(0.05);
	legend->AddEntry(hmjj_all,"all weights","l");
    legend->AddEntry(hmjj_positive,"positive weights only","l");
	legend->SetBorderSize(0);
    legend->Draw();

    //--------------------------------------------------------

    f2018_all->cd("dataset/Method_BDT/BDT");
    hmjj_all = (TH1D*)gDirectory->Get("mjj__Background");
    hmjj_all->SetLineColor(kRed);

    f2018_positive->cd("dataset/Method_BDT/BDT");
    hmjj_positive = (TH1D*)gDirectory->Get("mjj__Background");
    hmjj_positive->SetLineColor(kGreen);

    TCanvas *c2 = new TCanvas();
    c2->cd();
    hmjj_all->SetTitle("mjj 2018 (background)");
    hmjj_all->GetYaxis()->SetRangeUser(0,10);
    if (normalise)
    {
        hmjj_all->Scale(1/hmjj_all->Integral());
        hmjj_positive->Scale(1/hmjj_positive->Integral());
    }
    hmjj_all->Draw();
    hmjj_positive->Draw("same");

    legend = new TLegend(0.15,0.65,0.25,0.88);
	legend->SetTextFont(132);
	legend->SetTextSize(0.05);
	legend->AddEntry(hmjj_all,"all weights","l");
    legend->AddEntry(hmjj_positive,"positive weights only","l");
	legend->SetBorderSize(0);
    legend->Draw();
}