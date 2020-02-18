void mjj_all_years_signal()
{
    gStyle->SetOptStat(0);
    
    bool normalise = true;

    TFile *f2016_all = new TFile("./2016/all/TMVA.root");
    TFile *f2016_positive = new TFile("./2016/positive_only/TMVA.root");
    TFile *f2017_all = new TFile("./2017/all/TMVA.root");
    TFile *f2017_positive = new TFile("./2017/positive_only/TMVA.root");
    TFile *f2018_all = new TFile("./2018/all/TMVA.root");
    TFile *f2018_positive = new TFile("./2018/positive_only/TMVA.root");

    f2016_all->cd("dataset/Method_BDT/BDT");
    TH1D *hmjj_all_2016 = (TH1D*)gDirectory->Get("mjj__Signal");
    hmjj_all_2016->SetLineColor(kRed);

    f2017_all->cd("dataset/Method_BDT/BDT");
    TH1D *hmjj_all_2017 = (TH1D*)gDirectory->Get("mjj__Signal");
    hmjj_all_2017->SetLineColor(kGreen);

    f2018_all->cd("dataset/Method_BDT/BDT");
    TH1D *hmjj_all_2018 = (TH1D*)gDirectory->Get("mjj__Signal");
    hmjj_all_2018->SetLineColor(kBlue);

    TCanvas *c1 = new TCanvas();
    c1->cd();
    hmjj_all_2016->SetTitle("mjj signal (all weights)");
    if (normalise)
    {
        hmjj_all_2016->Scale(1/hmjj_all_2016->Integral());
        hmjj_all_2017->Scale(1/hmjj_all_2017->Integral());
        hmjj_all_2018->Scale(1/hmjj_all_2018->Integral());
    }
    hmjj_all_2016->GetYaxis()->SetRangeUser(0,0.12);
    hmjj_all_2016->Draw();
    hmjj_all_2017->Draw("same");
    hmjj_all_2018->Draw("same");

    TLegend *legend = new TLegend(0.15,0.65,0.25,0.88);
	legend->SetTextFont(132);
	legend->SetTextSize(0.05);
	legend->AddEntry(hmjj_all_2016,"2016","l");
    legend->AddEntry(hmjj_all_2017,"2017","l");
    legend->AddEntry(hmjj_all_2018,"2018","l");
	legend->SetBorderSize(0);
    legend->Draw();

    //---------------------------------------------------------------------

    f2016_positive->cd("dataset/Method_BDT/BDT");
    TH1D *hmjj_positive_2016 = (TH1D*)gDirectory->Get("mjj__Signal");
    hmjj_positive_2016->SetLineColor(kRed);

    f2017_positive->cd("dataset/Method_BDT/BDT");
    TH1D *hmjj_positive_2017 = (TH1D*)gDirectory->Get("mjj__Signal");
    hmjj_positive_2017->SetLineColor(kGreen);

    f2018_positive->cd("dataset/Method_BDT/BDT");
    TH1D *hmjj_positive_2018 = (TH1D*)gDirectory->Get("mjj__Signal");
    hmjj_positive_2018->SetLineColor(kBlue);

    TCanvas *c2 = new TCanvas();
    c2->cd();
    hmjj_positive_2016->SetTitle("mjj signal (positive weights only)");
    if (normalise)
    {
        hmjj_positive_2016->Scale(1/hmjj_positive_2016->Integral());
        hmjj_positive_2017->Scale(1/hmjj_positive_2017->Integral());
        hmjj_positive_2018->Scale(1/hmjj_positive_2018->Integral());
    }
    hmjj_positive_2016->GetYaxis()->SetRangeUser(0,0.12);
    hmjj_positive_2016->Draw();
    hmjj_positive_2017->Draw("same");
    hmjj_positive_2018->Draw("same");

    legend = new TLegend(0.15,0.65,0.25,0.88);
	legend->SetTextFont(132);
	legend->SetTextSize(0.05);
	legend->AddEntry(hmjj_positive_2016,"2016","l");
    legend->AddEntry(hmjj_positive_2017,"2017","l");
    legend->AddEntry(hmjj_positive_2018,"2018","l");
	legend->SetBorderSize(0);
    legend->Draw();
}
