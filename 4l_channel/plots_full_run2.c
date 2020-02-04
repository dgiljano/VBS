void plots_full_run2()
{
    // ---------------------------------------------- 2016 ------------------------------------------------

    TFile *fmjj_2016 = new TFile("./onlymjjCut_jet_pt_gt_30/2016/mjj_2016.root");

    TH1F *hvbs_mjj_2016 = new TH1F("vbs_mjj_2016","vbs_mjj_2016", 20, 100, 1000);
    hvbs_mjj_2016 = (TH1F*)fmjj_2016->Get("vbs");

    TH1F *hqq_mjj_2016 = new TH1F("QCD_qq_mjj_2016","QCD_qq_mjj_2016", 20, 100, 1000);
    hqq_mjj_2016 = (TH1F*)fmjj_2016->Get("QCD_qq");

    TH1F *hgg_mjj_2016 = new TH1F("QCD_gg_mjj_2016","QCD_gg_mjj_2016", 20, 100, 1000);
    hgg_mjj_2016 = (TH1F*)fmjj_2016->Get("QCD_gg");

    TH1F *hzx_mjj_2016 = new TH1F("zx_mjj_2016","zx_mjj_2016", 20, 100, 1000);
    hzx_mjj_2016 = (TH1F*)fmjj_2016->Get("ZX");

    TH1F *httz_mjj_2016 = new TH1F("ttz_mjj_2016","ttz_mjj_2016", 20, 100, 1000);
    httz_mjj_2016 = (TH1F*)fmjj_2016->Get("ttZ");

    TH1F *hwwz_mjj_2016 = new TH1F("wwz_mjj_2016","wwz_mjj_2016", 20, 100, 1000);
    hwwz_mjj_2016 = (TH1F*)fmjj_2016->Get("WWZ");

    TH1F *hdata_mjj_2016 = new TH1F("data_mjj_2016","data_mjj_2016", 20, 100, 1000);
    hdata_mjj_2016 = (TH1F*)fmjj_2016->Get("data");

    TH1F *httzwwz_mjj_2016 = new TH1F("ttzwwz_mjj_2016","ttzwwz_mjj_2016", 20, 100, 1000);
    httzwwz_mjj_2016 = (TH1F*)httz_mjj_2016->Clone();
    httzwwz_mjj_2016->Add(hwwz_mjj_2016);

    
    TFile *fdetajj_2016 = new TFile("./onlymjjCut_jet_pt_gt_30/2016/detajj_2016.root");

    TH1F *hvbs_detajj_2016 = new TH1F("vbs_detajj_2016","vbs_detajj_2016", 20, 100, 1000);
    hvbs_detajj_2016 = (TH1F*)fdetajj_2016->Get("vbs");

    TH1F *hqq_detajj_2016 = new TH1F("QCD_qq_detajj_2016","QCD_qq_detajj_2016", 20, 100, 1000);
    hqq_detajj_2016 = (TH1F*)fdetajj_2016->Get("QCD_qq");

    TH1F *hgg_detajj_2016 = new TH1F("QCD_gg_detajj_2016","QCD_gg_detajj_2016", 20, 100, 1000);
    hgg_detajj_2016 = (TH1F*)fdetajj_2016->Get("QCD_gg");

    TH1F *hzx_detajj_2016 = new TH1F("zx_detajj_2016","zx_detajj_2016", 20, 100, 1000);
    hzx_detajj_2016 = (TH1F*)fdetajj_2016->Get("ZX");

    TH1F *httz_detajj_2016 = new TH1F("ttz_detajj_2016","ttz_detajj_2016", 20, 100, 1000);
    httz_detajj_2016 = (TH1F*)fdetajj_2016->Get("ttZ");

    TH1F *hwwz_detajj_2016 = new TH1F("wwz_detajj_2016","wwz_detajj_2016", 20, 100, 1000);
    hwwz_detajj_2016 = (TH1F*)fdetajj_2016->Get("WWZ");

    TH1F *hdata_detajj_2016 = new TH1F("data_detajj_2016","data_detajj_2016", 20, 100, 1000);
    hdata_detajj_2016 = (TH1F*)fdetajj_2016->Get("data");

    TH1F *httzwwz_detajj_2016 = new TH1F("ttzwwz_detajj_2016","ttzwwz_detajj_2016", 20, 100, 1000);
    httzwwz_detajj_2016 = (TH1F*)httz_detajj_2016->Clone();
    httzwwz_detajj_2016->Add(hwwz_detajj_2016);


    // ---------------------------------------------- 2017 ------------------------------------------------

    TFile *fmjj_2017 = new TFile("./onlymjjCut_jet_pt_gt_30/2017/mjj_2017.root");

    TH1F *hvbs_mjj_2017 = new TH1F("vbs_mjj_2017","vbs_mjj_2017", 20, 100, 1000);
    hvbs_mjj_2017 = (TH1F*)fmjj_2017->Get("vbs");

    TH1F *hqq_mjj_2017 = new TH1F("QCD_qq_mjj_2017","QCD_qq_mjj_2017", 20, 100, 1000);
    hqq_mjj_2017 = (TH1F*)fmjj_2017->Get("QCD_qq");

    TH1F *hgg_mjj_2017 = new TH1F("QCD_gg_mjj_2017","QCD_gg_mjj_2017", 20, 100, 1000);
    hgg_mjj_2017 = (TH1F*)fmjj_2017->Get("QCD_gg");

    TH1F *hzx_mjj_2017 = new TH1F("zx_mjj_2017","zx_mjj_2017", 20, 100, 1000);
    hzx_mjj_2017 = (TH1F*)fmjj_2017->Get("ZX");

    TH1F *httz_mjj_2017 = new TH1F("ttz_mjj_2017","ttz_mjj_2017", 20, 100, 1000);
    httz_mjj_2017 = (TH1F*)fmjj_2017->Get("ttZ");

    TH1F *hwwz_mjj_2017 = new TH1F("wwz_mjj_2017","wwz_mjj_2017", 20, 100, 1000);
    hwwz_mjj_2017 = (TH1F*)fmjj_2017->Get("WWZ");

    TH1F *hdata_mjj_2017 = new TH1F("data_mjj_2017","data_mjj_2017", 20, 100, 1000);
    hdata_mjj_2017 = (TH1F*)fmjj_2017->Get("data");

    TH1F *httzwwz_mjj_2017 = new TH1F("ttzwwz_mjj_2017","ttzwwz_mjj_2017", 20, 100, 1000);
    httzwwz_mjj_2017 = (TH1F*)httz_mjj_2017->Clone();
    httzwwz_mjj_2017->Add(hwwz_mjj_2017);

    
    TFile *fdetajj_2017 = new TFile("./onlymjjCut_jet_pt_gt_30/2017/detajj_2017.root");

    TH1F *hvbs_detajj_2017 = new TH1F("vbs_detajj_2017","vbs_detajj_2017", 20, 100, 1000);
    hvbs_detajj_2017 = (TH1F*)fdetajj_2017->Get("vbs");

    TH1F *hqq_detajj_2017 = new TH1F("QCD_qq_detajj_2017","QCD_qq_detajj_2017", 20, 100, 1000);
    hqq_detajj_2017 = (TH1F*)fdetajj_2017->Get("QCD_qq");

    TH1F *hgg_detajj_2017 = new TH1F("QCD_gg_detajj_2017","QCD_gg_detajj_2017", 20, 100, 1000);
    hgg_detajj_2017 = (TH1F*)fdetajj_2017->Get("QCD_gg");

    TH1F *hzx_detajj_2017 = new TH1F("zx_detajj_2017","zx_detajj_2017", 20, 100, 1000);
    hzx_detajj_2017 = (TH1F*)fdetajj_2017->Get("ZX");

    TH1F *httz_detajj_2017 = new TH1F("ttz_detajj_2017","ttz_detajj_2017", 20, 100, 1000);
    httz_detajj_2017 = (TH1F*)fdetajj_2017->Get("ttZ");

    TH1F *hwwz_detajj_2017 = new TH1F("wwz_detajj_2017","wwz_detajj_2017", 20, 100, 1000);
    hwwz_detajj_2017 = (TH1F*)fdetajj_2017->Get("WWZ");

    TH1F *hdata_detajj_2017 = new TH1F("data_detajj_2017","data_detajj_2017", 20, 100, 1000);
    hdata_detajj_2017 = (TH1F*)fdetajj_2017->Get("data");

    TH1F *httzwwz_detajj_2017 = new TH1F("ttzwwz_detajj_2017","ttzwwz_detajj_2017", 20, 100, 1000);
    httzwwz_detajj_2017 = (TH1F*)httz_detajj_2017->Clone();
    httzwwz_detajj_2017->Add(hwwz_detajj_2017);


    // ---------------------------------------------- 2018 ------------------------------------------------

    TFile *fmjj_2018 = new TFile("./onlymjjCut_jet_pt_gt_30/2018/mjj_2018.root");

    TH1F *hvbs_mjj_2018 = new TH1F("vbs_mjj_2018","vbs_mjj_2018", 20, 100, 1000);
    hvbs_mjj_2018 = (TH1F*)fmjj_2018->Get("vbs");

    TH1F *hqq_mjj_2018 = new TH1F("QCD_qq_mjj_2018","QCD_qq_mjj_2018", 20, 100, 1000);
    hqq_mjj_2018 = (TH1F*)fmjj_2018->Get("QCD_qq");

    TH1F *hgg_mjj_2018 = new TH1F("QCD_gg_mjj_2018","QCD_gg_mjj_2018", 20, 100, 1000);
    hgg_mjj_2018 = (TH1F*)fmjj_2018->Get("QCD_gg");

    TH1F *hzx_mjj_2018 = new TH1F("zx_mjj_2018","zx_mjj_2018", 20, 100, 1000);
    hzx_mjj_2018 = (TH1F*)fmjj_2018->Get("ZX");

    TH1F *httz_mjj_2018 = new TH1F("ttz_mjj_2018","ttz_mjj_2018", 20, 100, 1000);
    httz_mjj_2018 = (TH1F*)fmjj_2018->Get("ttZ");

    TH1F *hwwz_mjj_2018 = new TH1F("wwz_mjj_2018","wwz_mjj_2018", 20, 100, 1000);
    hwwz_mjj_2018 = (TH1F*)fmjj_2018->Get("WWZ");

    TH1F *hdata_mjj_2018 = new TH1F("data_mjj_2018","data_mjj_2018", 20, 100, 1000);
    hdata_mjj_2018 = (TH1F*)fmjj_2018->Get("data");

    TH1F *httzwwz_mjj_2018 = new TH1F("ttzwwz_mjj_2018","ttzwwz_mjj_2018", 20, 100, 1000);
    httzwwz_mjj_2018 = (TH1F*)httz_mjj_2018->Clone();
    httzwwz_mjj_2018->Add(hwwz_mjj_2018);

    
    TFile *fdetajj_2018 = new TFile("./onlymjjCut_jet_pt_gt_30/2018/detajj_2018.root");

    TH1F *hvbs_detajj_2018 = new TH1F("vbs_detajj_2018","vbs_detajj_2018", 20, 100, 1000);
    hvbs_detajj_2018 = (TH1F*)fdetajj_2018->Get("vbs");

    TH1F *hqq_detajj_2018 = new TH1F("QCD_qq_detajj_2018","QCD_qq_detajj_2018", 20, 100, 1000);
    hqq_detajj_2018 = (TH1F*)fdetajj_2018->Get("QCD_qq");

    TH1F *hgg_detajj_2018 = new TH1F("QCD_gg_detajj_2018","QCD_gg_detajj_2018", 20, 100, 1000);
    hgg_detajj_2018 = (TH1F*)fdetajj_2018->Get("QCD_gg");

    TH1F *hzx_detajj_2018 = new TH1F("zx_detajj_2018","zx_detajj_2018", 20, 100, 1000);
    hzx_detajj_2018 = (TH1F*)fdetajj_2018->Get("ZX");

    TH1F *httz_detajj_2018 = new TH1F("ttz_detajj_2018","ttz_detajj_2018", 20, 100, 1000);
    httz_detajj_2018 = (TH1F*)fdetajj_2018->Get("ttZ");

    TH1F *hwwz_detajj_2018 = new TH1F("wwz_detajj_2018","wwz_detajj_2018", 20, 100, 1000);
    hwwz_detajj_2018 = (TH1F*)fdetajj_2018->Get("WWZ");

    TH1F *hdata_detajj_2018 = new TH1F("data_detajj_2018","data_detajj_2018", 20, 100, 1000);
    hdata_detajj_2018 = (TH1F*)fdetajj_2018->Get("data");

    TH1F *httzwwz_detajj_2018 = new TH1F("ttzwwz_detajj_2018","ttzwwz_detajj_2018", 20, 100, 1000);
    httzwwz_detajj_2018 = (TH1F*)httz_detajj_2018->Clone();
    httzwwz_detajj_2018->Add(hwwz_detajj_2018);


    // ---------------------------------------- combined -------------------------------------------------

    TH1F *hvbs_mjj = new TH1F("vbs_mjj","vbs_mjj", 20, 100, 1000);
    hvbs_mjj = (TH1F*)hvbs_mjj_2016->Clone();
    hvbs_mjj->Add(hvbs_mjj_2017);
    hvbs_mjj->Add(hvbs_mjj_2018);

    TH1F *hqq_mjj = new TH1F("qq_mjj","qq_mjj", 20, 100, 1000);
    hqq_mjj = (TH1F*)hqq_mjj_2016->Clone();
    hqq_mjj->Add(hqq_mjj_2017);
    hqq_mjj->Add(hqq_mjj_2018);

    TH1F *hgg_mjj = new TH1F("gg_mjj","gg_mjj", 20, 100, 1000);
    hgg_mjj = (TH1F*)hgg_mjj_2016->Clone();
    hgg_mjj->Add(hgg_mjj_2017);
    hgg_mjj->Add(hgg_mjj_2018);

    TH1F *hzx_mjj = new TH1F("zx_mjj","zx_mjj", 20, 100, 1000);
    hzx_mjj = (TH1F*)hzx_mjj_2016->Clone();
    hzx_mjj->Add(hzx_mjj_2017);
    hzx_mjj->Add(hzx_mjj_2018);

    TH1F *httzwwz_mjj = new TH1F("ttzwwz_mjj","ttzwwz_mjj", 20, 100, 1000);
    httzwwz_mjj = (TH1F*)httzwwz_mjj_2016->Clone();
    httzwwz_mjj->Add(httzwwz_mjj_2017);
    httzwwz_mjj->Add(httzwwz_mjj_2018);

    TH1F *hdata_mjj = new TH1F("data_mjj","data_mjj", 20, 100, 1000);
    hdata_mjj = (TH1F*)hdata_mjj_2016->Clone();
    hdata_mjj->Add(hdata_mjj_2017);
    hdata_mjj->Add(hdata_mjj_2018);
    
    
    TH1F *hvbs_detajj = new TH1F("vbs_detajj","vbs_detajj", 20, 100, 1000);
    hvbs_detajj = (TH1F*)hvbs_detajj_2016->Clone();
    hvbs_detajj->Add(hvbs_detajj_2017);
    hvbs_detajj->Add(hvbs_detajj_2018);

    TH1F *hqq_detajj = new TH1F("qq_detajj","qq_detajj", 20, 100, 1000);
    hqq_detajj = (TH1F*)hqq_detajj_2016->Clone();
    hqq_detajj->Add(hqq_detajj_2017);
    hqq_detajj->Add(hqq_detajj_2018);

    TH1F *hgg_detajj = new TH1F("gg_detajj","gg_detajj", 20, 100, 1000);
    hgg_detajj = (TH1F*)hgg_detajj_2016->Clone();
    hgg_detajj->Add(hgg_detajj_2017);
    hgg_detajj->Add(hgg_detajj_2018);

    TH1F *hzx_detajj = new TH1F("zx_detajj","zx_detajj", 20, 100, 1000);
    hzx_detajj = (TH1F*)hzx_detajj_2016->Clone();
    hzx_detajj->Add(hzx_detajj_2017);
    hzx_detajj->Add(hzx_detajj_2018);

    TH1F *httzwwz_detajj = new TH1F("ttzwwz_detajj","ttzwwz_detajj", 20, 100, 1000);
    httzwwz_detajj = (TH1F*)httzwwz_detajj_2016->Clone();
    httzwwz_detajj->Add(httzwwz_detajj_2017);
    httzwwz_detajj->Add(httzwwz_detajj_2018);

    TH1F *hdata_detajj = new TH1F("data_detajj","data_detajj", 20, 100, 1000);
    hdata_detajj = (TH1F*)hdata_detajj_2016->Clone();
    hdata_detajj->Add(hdata_detajj_2017);
    hdata_detajj->Add(hdata_detajj_2018);


    //------------------------------------------------ making stack -----------------------------------------

    hzx_detajj->SetFillColor(kGreen);
    httzwwz_detajj->SetFillColor(kYellow);
    hqq_detajj->SetFillColor(kCyan);
    hgg_detajj->SetFillColor(kBlue);
    hvbs_detajj->SetFillColor(kMagenta);
    hdata_detajj->SetMarkerStyle(20);

    hzx_detajj->SetFillColor(kGreen);
    httzwwz_detajj->Add(httzwwz_detajj,hzx_detajj,1,1);    //tt
	httzwwz_detajj->SetFillColor(kYellow);
	hqq_detajj->Add(httzwwz_detajj,hqq_detajj,1,1); //real ew
	hqq_detajj->SetFillColor(kCyan);
	hgg_detajj->Add(hqq_detajj,hgg_detajj,1,1); //real gg
	hgg_detajj->SetFillColor(kBlue);
	hvbs_detajj->Add(hgg_detajj,hvbs_detajj,1,1);    //real vbs
	hvbs_detajj->SetFillColor(kMagenta);

    char filetitle[300];
	sprintf(filetitle,"CMS Preliminary                                                               137.1 fb^{-1}");  
    THStack *hs_detajj = new THStack("",filetitle);

    hs_detajj->Add(hvbs_detajj,"hist");
	hs_detajj->Add(hgg_detajj,"hist");
	hs_detajj->Add(hqq_detajj,"hist");
    hs_detajj->Add(httzwwz_detajj,"hist");
	hs_detajj->Add(hzx_detajj,"hist");
    hs_detajj->Add(hdata_detajj,"E1");

    TCanvas *c1 = new TCanvas("c1","example",800,1000);
    //gPad->SetLogy();
    c1->cd();

    hs_detajj->Draw("nostack");

    TLegend *legend = new TLegend(0.7,0.65,0.895,0.88);
	legend->SetTextFont(132);
	legend->SetTextSize(0.02);
	legend->AddEntry(hzx_detajj,"Z+X","f");
    legend->AddEntry(httzwwz_detajj,"t#bar{t}Z, WWZ","f");
	legend->AddEntry(hqq_detajj,"q#bar{q}#rightarrowZZ","f");
	legend->AddEntry(hgg_detajj,"gg#rightarrowZZ","f");
	legend->AddEntry(hvbs_detajj,"VBS","f");
	legend->AddEntry(hdata_detajj,"Data","lep");
	legend->SetBorderSize(0);
    legend->Draw();

    hs_detajj->GetXaxis()->SetTitle("#Delta#eta_{jj}");
    hs_detajj->GetYaxis()->SetTitle("Event/bin");
    hs_detajj->GetYaxis()->SetTitleOffset(1.4);
    hs_detajj->GetXaxis()->SetTitleOffset(1.4);

    hs_detajj->SetMinimum(0);




    hzx_mjj->SetFillColor(kGreen);
    httzwwz_mjj->SetFillColor(kYellow);
    hqq_mjj->SetFillColor(kCyan);
    hgg_mjj->SetFillColor(kBlue);
    hvbs_mjj->SetFillColor(kMagenta);
    hdata_mjj->SetMarkerStyle(20);

    hzx_mjj->SetFillColor(kGreen);
    httzwwz_mjj->Add(httzwwz_mjj,hzx_mjj,1,1);    //tt
	httzwwz_mjj->SetFillColor(kYellow);
	hqq_mjj->Add(httzwwz_mjj,hqq_mjj,1,1); //real ew
	hqq_mjj->SetFillColor(kCyan);
	hgg_mjj->Add(hqq_mjj,hgg_mjj,1,1); //real gg
	hgg_mjj->SetFillColor(kBlue);
	hvbs_mjj->Add(hgg_mjj,hvbs_mjj,1,1);    //real vbs
	hvbs_mjj->SetFillColor(kMagenta);

    filetitle[300];
	sprintf(filetitle,"CMS Preliminary                                                               137.1 fb^{-1}");  
    THStack *hs_mjj = new THStack("",filetitle);

    hs_mjj->Add(hvbs_mjj,"hist");
	hs_mjj->Add(hgg_mjj,"hist");
	hs_mjj->Add(hqq_mjj,"hist");
    hs_mjj->Add(httzwwz_mjj,"hist");
	hs_mjj->Add(hzx_mjj,"hist");
    hs_mjj->Add(hdata_mjj,"E1");

    TCanvas *c2 = new TCanvas("c2","example",800,1000);
    //gPad->SetLogy();
    c2->cd();

    hs_mjj->Draw("nostack");

    legend = new TLegend(0.7,0.65,0.895,0.88);
	legend->SetTextFont(132);
	legend->SetTextSize(0.02);
	legend->AddEntry(hzx_mjj,"Z+X","f");
    legend->AddEntry(httzwwz_mjj,"t#bar{t}Z, WWZ","f");
	legend->AddEntry(hqq_mjj,"q#bar{q}#rightarrowZZ","f");
	legend->AddEntry(hgg_mjj,"gg#rightarrowZZ","f");
	legend->AddEntry(hvbs_mjj,"VBS","f");
	legend->AddEntry(hdata_mjj,"Data","lep");
	legend->SetBorderSize(0);
    legend->Draw();

    hs_mjj->GetXaxis()->SetTitle("m_{jj} [GeV]");
    hs_mjj->GetYaxis()->SetTitle("Event/bin");
    hs_mjj->GetYaxis()->SetTitleOffset(1.4);
    hs_mjj->GetXaxis()->SetTitleOffset(1.4);

    hs_mjj->SetMinimum(0);
    
}