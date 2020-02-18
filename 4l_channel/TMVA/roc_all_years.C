void roc_all_years()
{
    gStyle->SetOptStat(0);


    TFile *f2016_all = new TFile("./2016/all/TMVA.root");
    TFile *f2016_positive = new TFile("./2016/positive_only/TMVA.root");
    TFile *f2017_all = new TFile("./2017/all/TMVA.root");
    TFile *f2017_positive = new TFile("./2017/positive_only/TMVA.root");
    TFile *f2018_all = new TFile("./2018/all/TMVA.root");
    TFile *f2018_positive = new TFile("./2018/positive_only/TMVA.root");

    f2016_all->cd("dataset/Method_BDT/BDT");
    TH1D *hROC_all_2016 = (TH1D*)gDirectory->Get("MVA_BDT_rejBvsS");
    hROC_all_2016->SetLineColor(kRed);

    f2017_all->cd("dataset/Method_BDT/BDT");
    TH1D *hROC_all_2017 = (TH1D*)gDirectory->Get("MVA_BDT_rejBvsS");
    hROC_all_2017->SetLineColor(kGreen);

    f2018_all->cd("dataset/Method_BDT/BDT");
    TH1D *hROC_all_2018 = (TH1D*)gDirectory->Get("MVA_BDT_rejBvsS");
    hROC_all_2018->SetLineColor(kBlue);

    TCanvas *c1 = new TCanvas();
    c1->cd();
    hROC_all_2016->SetTitle("All years (all weights)");
    hROC_all_2016->Draw();
    hROC_all_2017->Draw("same");
    hROC_all_2018->Draw("same");

    TLegend *legend = new TLegend(0.15,0.45,0.25,0.68);
	legend->SetTextFont(132);
	legend->SetTextSize(0.05);
	legend->AddEntry(hROC_all_2016,"2016","l");
    legend->AddEntry(hROC_all_2017,"2017","l");
    legend->AddEntry(hROC_all_2018,"2018","l");
	legend->SetBorderSize(0);
    legend->Draw();

    //---------------------------------------------------------------------

    f2016_positive->cd("dataset/Method_BDT/BDT");
    TH1D *hROC_positive_2016 = (TH1D*)gDirectory->Get("MVA_BDT_rejBvsS");
    hROC_positive_2016->SetLineColor(kRed);

    f2017_positive->cd("dataset/Method_BDT/BDT");
    TH1D *hROC_positive_2017 = (TH1D*)gDirectory->Get("MVA_BDT_rejBvsS");
    hROC_positive_2017->SetLineColor(kGreen);

    f2018_positive->cd("dataset/Method_BDT/BDT");
    TH1D *hROC_positive_2018 = (TH1D*)gDirectory->Get("MVA_BDT_rejBvsS");
    hROC_positive_2018->SetLineColor(kBlue);

    TCanvas *c2 = new TCanvas();
    c2->cd();
    hROC_positive_2016->SetTitle("All years (positive weights only)");
    hROC_positive_2016->Draw();
    hROC_positive_2017->Draw("same");
    hROC_positive_2018->Draw("same");

    legend = new TLegend(0.15,0.45,0.25,0.68);
	legend->SetTextFont(132);
	legend->SetTextSize(0.05);
	legend->AddEntry(hROC_positive_2016,"2016","l");
    legend->AddEntry(hROC_positive_2017,"2017","l");
    legend->AddEntry(hROC_positive_2018,"2018","l");
	legend->SetBorderSize(0);
    legend->Draw();
}