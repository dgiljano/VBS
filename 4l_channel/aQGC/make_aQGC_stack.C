void make_aQGC_stack()
{
    Float_t bins_FT[] = { 0, 200, 300, 400, 500, 600, 800, 1000, 1200, 1400 };

    // ---------------------------------------------- 2016 histograms ------------------------------------------------

    TFile *f_2016 = new TFile("../onlymjjCut_jet_pt_gt_30/m4l_histos_2016.root");

    TH1F *hsum2_2016 = new TH1F("hsum2_2016","hsum2_2016", 9, bins_FT);
    hsum2_2016 = (TH1F*)f_2016->Get("hsum2_1");

    TH1F *hsum1_2016 = new TH1F("hsum1_2016","hsum1_2016", 9, bins_FT);
    hsum1_2016 = (TH1F*)f_2016->Get("hsum1_1");

    TH1F *hqqzz_2016 = new TH1F("hqqzz_2016","hqqzz_2016", 9, bins_FT);
    hqqzz_2016 = (TH1F*)f_2016->Get("hqqzz_1");

    TH1F *httzwwz_2016 = new TH1F("httzwwz_2016","httzwwz_2016", 9, bins_FT);
    httzwwz_2016 = (TH1F*)f_2016->Get("httzwzz_1");

    TH1F *hzx_2016 = new TH1F("hzx_2016","hzx_2016", 9, bins_FT);
    hzx_2016 = (TH1F*)f_2016->Get("hzx1");

    TH1F *hdata_2016 = new TH1F("hdata_2016","hdata_2016", 9, bins_FT);
    hdata_2016 = (TH1F*)f_2016->Get("hdata_1");

    TH1F *hft8_2016 = new TH1F("hft8_2016","hft8_2016", 9, bins_FT);
    hft8_2016 = (TH1F*)f_2016->Get("BLS_hvbs_1_rescaled_FT8_1");

    TH1F *hft9_2016 = new TH1F("hft9_2016","hft9_2016", 9, bins_FT);
    hft9_2016 = (TH1F*)f_2016->Get("BLS_hvbs_1_rescaled_FT9_2");

    // ---------------------------------------------- 2017 histograms ------------------------------------------------

    TFile *f_2017 = new TFile("../onlymjjCut_jet_pt_gt_30/m4l_histos_2017.root");

    TH1F *hsum2_2017 = new TH1F("hsum2_2017","hsum2_2017", 9, bins_FT);
    hsum2_2017 = (TH1F*)f_2017->Get("hsum2_1");

    TH1F *hsum1_2017 = new TH1F("hsum1_2017","hsum1_2017", 9, bins_FT);
    hsum1_2017 = (TH1F*)f_2017->Get("hsum1_1");

    TH1F *hqqzz_2017 = new TH1F("hqqzz_2017","hqqzz_2017", 9, bins_FT);
    hqqzz_2017 = (TH1F*)f_2017->Get("hqqzz_1");

    TH1F *httzwwz_2017 = new TH1F("httzwwz_2017","httzwwz_2017", 9, bins_FT);
    httzwwz_2017 = (TH1F*)f_2017->Get("httzwzz_1");

    TH1F *hzx_2017 = new TH1F("hzx_2017","hzx_2017", 9, bins_FT);
    hzx_2017 = (TH1F*)f_2017->Get("hzx1");

    TH1F *hdata_2017 = new TH1F("hdata_2017","hdata_2017", 9, bins_FT);
    hdata_2017 = (TH1F*)f_2017->Get("hdata_1");

    TH1F *hft8_2017 = new TH1F("hft8_2017","hft8_2017", 9, bins_FT);
    hft8_2017 = (TH1F*)f_2017->Get("BLS_diboson_rescaled_FT8_1");

    TH1F *hft9_2017 = new TH1F("hft9_2017","hft9_2017", 9, bins_FT);
    hft9_2017 = (TH1F*)f_2017->Get("BLS_diboson_rescaled_FT9_2");

    // ---------------------------------------------- 2018 histograms ------------------------------------------------

    TFile *f_2018 = new TFile("../onlymjjCut_jet_pt_gt_30/m4l_histos_2018.root");

    TH1F *hsum2_2018 = new TH1F("hsum2_2018","hsum2_2018", 9, bins_FT);
    hsum2_2018 = (TH1F*)f_2018->Get("hsum2_1");

    TH1F *hsum1_2018 = new TH1F("hsum1_2018","hsum1_2018", 9, bins_FT);
    hsum1_2018 = (TH1F*)f_2018->Get("hsum1_1");

    TH1F *hqqzz_2018 = new TH1F("hqqzz_2018","hqqzz_2018", 9, bins_FT);
    hqqzz_2018 = (TH1F*)f_2018->Get("hqqzz_1");

    TH1F *httzwwz_2018 = new TH1F("httzwwz_2018","httzwwz_2018", 9, bins_FT);
    httzwwz_2018 = (TH1F*)f_2018->Get("httzwzz_1");

    TH1F *hzx_2018 = new TH1F("hzx_2018","hzx_2018", 9, bins_FT);
    hzx_2018 = (TH1F*)f_2018->Get("hzx1");

    TH1F *hdata_2018 = new TH1F("hdata_2018","hdata_2018", 9, bins_FT);
    hdata_2018 = (TH1F*)f_2018->Get("hdata_1");

    TH1F *hft8_2018 = new TH1F("hft8_2018","hft8_2018", 9, bins_FT);
    hft8_2018 = (TH1F*)f_2018->Get("BLS_diboson_rescaled_FT8_1");

    TH1F *hft9_2018 = new TH1F("hft9_2018","hft9_2018", 9, bins_FT);
    hft9_2018 = (TH1F*)f_2018->Get("BLS_diboson_rescaled_FT9_2");


    // -------------------------------------------- summing years -----------------------------------------------------

    TH1F *hsum2 = new TH1F("hsum2","hsum2", 9, bins_FT);
    TH1F *hsum1 = new TH1F("hsum1","hsum1", 9, bins_FT);
    TH1F *hqqzz = new TH1F("hqqzz","hqqzz", 9, bins_FT);
    TH1F *httzwwz = new TH1F("httzwwz","httzwwz", 9, bins_FT);
    TH1F *hzx = new TH1F("hzx","hzx", 9, bins_FT);
    TH1F *hdata = new TH1F("hdata","hdata", 9, bins_FT);
    TH1F *hft8 = new TH1F("hft8","hft8", 9, bins_FT);
    TH1F *hft9 = new TH1F("hft9","hft9", 9, bins_FT);

    for (int i = 1; i < 10; i++)
    {
        hsum2->SetBinContent(i, hsum2_2016->GetBinContent(i) + hsum2_2017->GetBinContent(i) + hsum2_2018->GetBinContent(i));
        hsum1->SetBinContent(i, hsum1_2016->GetBinContent(i) + hsum1_2017->GetBinContent(i) + hsum1_2018->GetBinContent(i));
        hqqzz->SetBinContent(i, hqqzz_2016->GetBinContent(i) + hqqzz_2017->GetBinContent(i) + hqqzz_2018->GetBinContent(i));
        httzwwz->SetBinContent(i, httzwwz_2016->GetBinContent(i) + httzwwz_2017->GetBinContent(i) + httzwwz_2018->GetBinContent(i));
        hzx->SetBinContent(i, hzx_2016->GetBinContent(i) + hzx_2017->GetBinContent(i) + hzx_2018->GetBinContent(i));
        hdata->SetBinContent(i, hdata_2016->GetBinContent(i) + hdata_2017->GetBinContent(i) + hdata_2018->GetBinContent(i));
        hft8->SetBinContent(i, hft8_2016->GetBinContent(i) + hft8_2017->GetBinContent(i) + hft8_2018->GetBinContent(i));
        hft9->SetBinContent(i, hft9_2016->GetBinContent(i) + hft9_2017->GetBinContent(i) + hft9_2018->GetBinContent(i));
    }

    hzx->SetFillColor(kGreen);
    httzwwz->SetFillColor(kYellow);
    hqqzz->SetFillColor(kCyan);
    hsum1->SetFillColor(kBlue);
    hsum2->SetFillColor(kMagenta);
    hdata->SetMarkerStyle(20);

    hft8->SetLineColor(kYellow+1);
	hft8->SetLineWidth(3);
	hft8->SetLineStyle(9);

    hft9->SetLineColor(kRed+1);
	hft9->SetLineWidth(3);
	hft9->SetLineStyle(9);

    // ------------------------------------ stacking histograms --------------------------------------------------------

    //THStack *hs = new THStack();
    char filetitle[300];
	sprintf(filetitle,"CMS Preliminary                                                               137.1 fb^{-1}");  
    THStack *hs = new THStack("",filetitle);

    hs->Add(hsum2,"hist");
	hs->Add(hsum1,"hist");
	hs->Add(hqqzz,"hist");
    hs->Add(httzwwz,"hist");
	hs->Add(hzx,"hist");
    hs->Add(hft8);
	hs->Add(hft9);
    hs->Add(hdata,"E1");

    TCanvas *c1 = new TCanvas("c1","example",800,1000);
    gPad->SetLogy();
    c1->cd();

    hs->Draw("nostack");

    TLegend *legend = new TLegend(0.7,0.65,0.895,0.88);
	legend->SetTextFont(132);
	legend->SetTextSize(0.02);
	legend->AddEntry(hzx,"Z+X","f");
    legend->AddEntry(httzwwz,"t#bar{t}Z, WWZ","f");
	legend->AddEntry(hqqzz,"q#bar{q}#rightarrowZZ","f");
	legend->AddEntry(hsum1,"gg#rightarrowZZ","f");
	legend->AddEntry(hsum2,"VBS","f");
	legend->AddEntry(hdata,"Data","lep");
    legend->AddEntry(hft8,"f_{T8} /#Lambda^{4} = 1 TeV^{-4}","f");
	legend->AddEntry(hft9,"f_{T9} /#Lambda^{4} = 2 TeV^{-4}","f");
	legend->SetBorderSize(0);
    legend->Draw();

    hs->GetXaxis()->SetTitle("M_{4l} [GeV]");
    hs->GetYaxis()->SetTitle("Event/bin");
    hs->GetYaxis()->SetTitleOffset(1.4);
    hs->GetXaxis()->SetTitleOffset(1.4);

    hs->SetMinimum(0.2);

    c1->SaveAs("m4l_plot_allMCatNLO.pdf");
}