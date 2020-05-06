void make_aTGC_stack()
{
    Float_t bins_FT[] = { 0, 200, 300, 400, 500, 600, 800, 1000, 1200, 1400 };

    // ---------------------------------------------- 2016 histograms ------------------------------------------------

    TFile *f_2016 = new TFile("../onlymjjCut_jet_pt_gt_30/2016/m4l_histos_2016.root");

    TH1F *hsum2_2016 = new TH1F("hsum2_2016","hsum2_2016", 9, bins_FT); //qq
    hsum2_2016 = (TH1F*)f_2016->Get("hsum2_1");

    TH1F *hsum1_2016 = new TH1F("hsum1_2016","hsum1_2016", 9, bins_FT); //gg
    hsum1_2016 = (TH1F*)f_2016->Get("hsum1_1");

    TH1F *httzwwz_2016 = new TH1F("httzwwz_2016","httzwwz_2016", 9, bins_FT); //ttZWWzWZZ
    httzwwz_2016 = (TH1F*)f_2016->Get("httzwzz_1");

    TH1F *hzx_2016 = new TH1F("hzx_2016","hzx_2016", 9, bins_FT);   //ZX
    hzx_2016 = (TH1F*)f_2016->Get("hzx1");

    TH1F *hsum2_atgc_2016 = new TH1F("hsum2_atgc_2016","hsum2_atgc_2016", 9, bins_FT);  //aTGC
    hsum2_atgc_2016 = (TH1F*)f_2016->Get("hsum2_atgc_1");

    // ---------------------------------------------- 2017 histograms ------------------------------------------------

    TFile *f_2017 = new TFile("../onlymjjCut_jet_pt_gt_30/2017/m4l_histos_2017.root");

    TH1F *hsum2_2017 = new TH1F("hsum2_2017","hsum2_2017", 9, bins_FT); //qq
    hsum2_2017 = (TH1F*)f_2017->Get("hsum2_1");

    TH1F *hsum1_2017 = new TH1F("hsum1_2017","hsum1_2017", 9, bins_FT); //gg
    hsum1_2017 = (TH1F*)f_2017->Get("hsum1_1");

    TH1F *httzwwz_2017 = new TH1F("httzwwz_2017","httzwwz_2017", 9, bins_FT);   //ttZWWzWZZ
    httzwwz_2017 = (TH1F*)f_2017->Get("httzwzz_1");

    TH1F *hzx_2017 = new TH1F("hzx_2017","hzx_2017", 9, bins_FT);   //ZX
    hzx_2017 = (TH1F*)f_2017->Get("hzx1");

    // ---------------------------------------------- 2018 histograms ------------------------------------------------

    TFile *f_2018 = new TFile("../onlymjjCut_jet_pt_gt_30/2018/m4l_histos_2018.root");

    TH1F *hsum2_2018 = new TH1F("hsum2_2018","hsum2_2018", 9, bins_FT); //qq
    hsum2_2018 = (TH1F*)f_2018->Get("hsum2_1");

    TH1F *hsum1_2018 = new TH1F("hsum1_2018","hsum1_2018", 9, bins_FT); //gg
    hsum1_2018 = (TH1F*)f_2018->Get("hsum1_1");

    TH1F *httzwwz_2018 = new TH1F("httzwwz_2018","httzwwz_2018", 9, bins_FT);   //ttZWWzWZZ
    httzwwz_2018 = (TH1F*)f_2018->Get("httzwzz_1");

    TH1F *hzx_2018 = new TH1F("hzx_2018","hzx_2018", 9, bins_FT);   //ZX
    hzx_2018 = (TH1F*)f_2018->Get("hzx1");


    // -------------------------------------------- summing years -----------------------------------------------------

    TH1F *hsum2 = new TH1F("hsum2","hsum2", 9, bins_FT);
    TH1F *hsum1 = new TH1F("hsum1","hsum1", 9, bins_FT);
    TH1F *httzwwz = new TH1F("httzwwz","httzwwz", 9, bins_FT);
    TH1F *hzx = new TH1F("hzx","hzx", 9, bins_FT);

    TH1F *hsum2_atgc = new TH1F("hsum2_atgc","hsum2_atgc", 9, bins_FT);

    for (int i = 1; i < 10; i++)
    {
        hsum2->SetBinContent(i, hsum2_2016->GetBinContent(i) + hsum2_2017->GetBinContent(i) + hsum2_2018->GetBinContent(i));
        hsum1->SetBinContent(i, hsum1_2016->GetBinContent(i) + hsum1_2017->GetBinContent(i) + hsum1_2018->GetBinContent(i));
        httzwwz->SetBinContent(i, httzwwz_2016->GetBinContent(i) + httzwwz_2017->GetBinContent(i) + httzwwz_2018->GetBinContent(i));
        hzx->SetBinContent(i, hzx_2016->GetBinContent(i) + hzx_2017->GetBinContent(i) + hzx_2018->GetBinContent(i));

        hsum2_atgc->SetBinContent(i, hsum2_atgc_2016->GetBinContent(i) * 3.8189);
    }

    hzx->SetFillColor(kGreen);
    httzwwz->SetFillColor(kYellow);
    hsum1->SetFillColor(kBlue);
    hsum2->SetFillColor(kCyan);

    //hsum2_atgc->SetFillColor(kGray);
    hsum2_atgc->SetLineColor(kRed);
	hsum2_atgc->SetLineWidth(3);
	hsum2_atgc->SetLineStyle(9);

    // ------------------------------------ stacking histograms --------------------------------------------------------

    //THStack *hs = new THStack();
    char filetitle[300];
	sprintf(filetitle,"CMS Preliminary                                                               137.1 fb^{-1}");  
    THStack *hs = new THStack("",filetitle);

    hs->Add(hsum2,"hist");
    hs->Add(hsum2_atgc,"hist");
    //hs->Add(hsum2,"hist");
	hs->Add(hsum1,"hist");
    hs->Add(httzwwz,"hist");
	hs->Add(hzx,"hist");


    TCanvas *c1 = new TCanvas("c1","example",800,1000);
    gPad->SetLogy();
    c1->cd();

    hs->Draw("nostack");

    TLegend *legend=new TLegend(0.6,0.55,0.75,0.88);
	legend->SetTextFont(72);
	legend->SetTextSize(0.04);
	legend->AddEntry(hzx,"Z+X","f");
	legend->AddEntry(httzwwz,"t#bar{t}Z, WWZ, WZZ","f");
	legend->AddEntry(hsum1,"gg#rightarrowZZ","f");
	legend->AddEntry(hsum2,"qq#rightarrowZZ","f");
    legend->AddEntry(hsum2_atgc,"aTGC","f");
	legend->SetBorderSize(0);
    legend->Draw("same");

    hs->GetXaxis()->SetTitle("M_{4l} [GeV]");
    hs->GetYaxis()->SetTitle("Event/bin");
    hs->GetYaxis()->SetTitleOffset(1.4);
    hs->GetXaxis()->SetTitleOffset(1.4);

    hs->SetMinimum(0.2);

    c1->SaveAs("m4l_plot_aTGC.pdf");
}