#include "external_cConstants.h"
void plotter(){

	//histogram stack
	THStack *hs = new THStack("hs","Kinematic variable distribution;K_{D};Events/0.05");

	//histograms
	TH1F *h_complete_data = new TH1F("","",20,0.,1.); //all data
	TH1F *h00 = new TH1F("h00","",20,0.,1.); //bkg_kin<0.7 cut full background plot
	TH1F *h0 = new TH1F("h0","",20,0.,1.); //data, because we are hiding higher energies in this phase
	TH1F *h1 = new TH1F("h1","h1;xtitle;ytitle",20,0.,1.); //ew
	TH1F *h1bis = new TH1F("h1bis","",20,0.,1.); //ew+zx -> real ew
	TH1F *h2 = new TH1F("h2","",20,0.,1.); //gg
        TH1F *h3 = new TH1F("h3","",20,0.,1.); //vbs
	TH1F *h4 = new TH1F("h4","",20,0.,1.); //gg+ew -> real gg
        TH1F *h5 = new TH1F("h5","",20,0.,1.); //gg+ew+vbs
	TH1F *h6 = new TH1F("h6","",20,0.,1.); //gg+ew+vbs+zx ->real vbs
	TH1F *h7 = new TH1F("h7","",20,0.,1.); //qqzz e
	TH1F *h8 = new TH1F("h8","",20,0.,1.); //qqzz mu
	TH1F *h9 = new TH1F("h9","",20,0.,1.); //qqzz e mu
	TH1F *h10 = new TH1F("h10","",20,0.,1.);//ggzz e
	TH1F *h11 = new TH1F("h11","",20,0.,1.);//ggzz mu
	TH1F *h12 = new TH1F("h12","",20,0.,1.); //ggzz e mu
	TH1F *h13 = new TH1F("h13","",20,0.,1.);//vbs e
	TH1F *h14 = new TH1F("h14","",20,0.,1.);//vbs mu
	TH1F *h15 = new TH1F("h15","",20,0.,1.);//vbs e mu
	TH1F *h0_ee = new TH1F("h16","",20,0.,1.);//vbs e
        TH1F *h0_mm = new TH1F("h17","",20,0.,1.);//vbs mu
        TH1F *h0_em = new TH1F("h18","",20,0.,1.);//vbs e mu

	gStyle->SetPalette(1);
	TFile *input_file;
	TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
	//float lumi = 1;
	float lumi = 35.9E03;
	//for loop for different samples
	for(int j = 0; j < 4; j++){

		//print cycle
		std::cout << endl << j << endl;
		//open and close input file to get gen_su_weights (sum of all tree events weights)
		if (j==0)
		input_file= TFile::Open("root://lxcms03://data3/Higgs/170222/ZZTo4l/ZZ4lAnalysis.root");
		if (j==1) input_file= TFile::Open("root://lxcms03://data3/Higgs/170222/VBFTo4eJJ_0PMH125Contin_phantom128/ZZ4lAnalysis.root");
		if (j==2) input_file= TFile::Open("root://lxcms03://data3/Higgs/170222/ggTo4tau_Contin_MCFM701/ZZ4lAnalysis.root");
		if (j==3) input_file= TFile::Open("root://lxcms03//data3/Higgs/170222/AllData/ZZ4lAnalysis.root");


		TH1F *hCounters= (TH1F*)input_file->Get("ZZTree/Counters");
		float gen_sum_weights = hCounters->GetBinContent(40);
                std::cout<<endl<<j<<"  "<< gen_sum_weights<<endl;
		input_file->Close();

		//tchain and add function for multiple input files
		TChain *tqqzz= new TChain("ZZTree/candTree");
		if (j==0)
		tqqzz->Add("root://lxcms03://data3/Higgs/170222/ZZTo4l/ZZ4lAnalysis.root");
		if (j==1) {
                                tqqzz->Add("root://lxcms03://data3/Higgs/170222/ggTo2e2mu_Contin_MCFM701/ZZ4lAnalysis.root");
                                tqqzz->Add("root://lxcms03://data3/Higgs/170222/ggTo2e2tau_Contin_MCFM701/ZZ4lAnalysis.root");
                                tqqzz->Add("root://lxcms03://data3/Higgs/170222/ggTo2mu2tau_Contin_MCFM701/ZZ4lAnalysis.root");
                                tqqzz->Add("root://lxcms03://data3/Higgs/170222/ggTo4e_Contin_MCFM701/ZZ4lAnalysis.root");
                                tqqzz->Add("root://lxcms03://data3/Higgs/170222/ggTo4mu_Contin_MCFM701/ZZ4lAnalysis.root");
                                tqqzz->Add("root://lxcms03://data3/Higgs/170222/ggTo4tau_Contin_MCFM701/ZZ4lAnalysis.root");
                }
		if (j==2) {
		          	tqqzz->Add("root://lxcms03://data3/Higgs/170222/VBFTo4eJJ_0PMH125Contin_phantom128/ZZ4lAnalysis.root");
		                tqqzz->Add("root://lxcms03://data3/Higgs/170222/VBFTo4muJJ_0PMH125Contin_phantom128/ZZ4lAnalysis.root");
                		tqqzz->Add("root://lxcms03://data3/Higgs/170222/VBFTo2e2muJJ_0PMH125Contin_phantom128/ZZ4lAnalysis.root");
		}
		if (j==3)
                tqqzz->Add("root://lxcms03://data3/Higgs/170222/AllData/ZZ4lAnalysis.root");

		//histogram declaration
		//TH1F *kin_zz = new TH1F("kin_zz","",20,0.,1.); //was 100 bins

		//original variable declarations
		float ZZPt,ZZMass,DiJetMass;
		float xsec,KFactorEWKqqZZ,overallEventWeight,KFactorQCDqqZZ_M;
		vector<float> *LepPt=new vector<float>;
		short Z1Flav,Z2Flav;
		short nCleanedJetsPt30;
		float pvbf_VAJHU_old;
		float phjj_VAJHU_old;
		float bkg_VAMCFM,p0plus_VAJHU;
		short ZZsel;
		vector<short> *LepLepId=0;
		short nExtraLep;
		short nCleanedJetsPt30BTagged_bTagSF;

		//new variable declarations
		float p_JJEW_BKG_MCFM_JECNominal; //not in use
		float p_JJQCD_BKG_MCFM_JECNominal;
		float p_JJVBF_BKG_MCFM_JECNominal;
		float KFactorQCDggzz_Nominal;
	//	float lumi = 35.9E03;
		float my_sum = 0;
		float mc_integral;
		float data_integral;
		//original brach addresses
		tqqzz->SetBranchAddress("nExtraLep",&nExtraLep);
		tqqzz->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF",&nCleanedJetsPt30BTagged_bTagSF);
		tqqzz->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal",&pvbf_VAJHU_old);
		tqqzz->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal",&phjj_VAJHU_old);
		tqqzz->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1_1_JHUGen",&p0plus_VAJHU);
		tqqzz->SetBranchAddress("p_QQB_BKG_MCFM",&bkg_VAMCFM);
		tqqzz->SetBranchAddress("ZZPt",&ZZPt);
		tqqzz->SetBranchAddress("ZZMass",&ZZMass);
		tqqzz->SetBranchAddress("Z1Flav",&Z1Flav);
		tqqzz->SetBranchAddress("Z2Flav",&Z2Flav);
		tqqzz->SetBranchAddress("KFactor_EW_qqZZ",&KFactorEWKqqZZ);
		tqqzz->SetBranchAddress("xsec",&xsec);
		tqqzz->SetBranchAddress("ZZsel",&ZZsel);
		tqqzz->SetBranchAddress("LepLepId",&LepLepId);
		tqqzz->SetBranchAddress("LepPt",&LepPt);
		tqqzz->SetBranchAddress("overallEventWeight",&overallEventWeight);
		tqqzz->SetBranchAddress("KFactor_QCD_qqZZ_M",&KFactorQCDqqZZ_M);
		tqqzz->SetBranchAddress("nCleanedJetsPt30",&nCleanedJetsPt30);

		//new branch addresses
		tqqzz->SetBranchAddress("p_JJEW_BKG_MCFM_JECNominal",&p_JJEW_BKG_MCFM_JECNominal);
		tqqzz->SetBranchAddress("p_JJVBF_BKG_MCFM_JECNominal",&p_JJVBF_BKG_MCFM_JECNominal);
                tqqzz->SetBranchAddress("p_JJQCD_BKG_MCFM_JECNominal",&p_JJQCD_BKG_MCFM_JECNominal);
		tqqzz->SetBranchAddress("DiJetMass",&DiJetMass);
		tqqzz->SetBranchAddress("KFactor_QCD_ggZZ_Nominal",&KFactorQCDggzz_Nominal);
		//additional and output variable declarations
		float weight, weight_up, weight_dn;
		float weight_vbf, weight_vbf_up, weight_vbf_dn;
		int chan;
		int vbfcate = 0;
		float dbkg_kin;

		//output branches
		TTree *tnew = new TTree("SelectedTree","SelectedTree");
		tnew->Branch("mreco",&ZZMass,"mreco/F");
		tnew->Branch("dbkg_kin",&dbkg_kin,"dbkg_kin/F");
		tnew->Branch("weight",&weight,"weight/F");
		tnew->Branch("weight_up",&weight_up,"weight_up/F");
		tnew->Branch("weight_dn",&weight_dn,"weight_dn/F");
		tnew->Branch("weight_vbf",&weight_vbf,"weight_vbf/F");
		tnew->Branch("weight_vbf_up",&weight_vbf_up,"weight_vbf_up/F");
		tnew->Branch("weight_vbf_dn",&weight_vbf_dn,"weight_vbf_dn/F");
		tnew->Branch("chan",&chan,"chan/I");
		tnew->Branch("vbfcate",&vbfcate,"vbfcate/I");

		//loop on entries
		for(int i=0;i<tqqzz->GetEntries();i++){
			tqqzz->GetEntry(i);

			//unique selection condition (see paper page 8) & DiJetMass condition
			if(DiJetMass>100 && nExtraLep==0 && ZZMass > 160 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0))){

			//unique selection
		//	if(nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0))){

			//set vbf_category
                        	 vbfcate=1;
			//weight definition
			//KFactorEWKqqZZ = 1;
			//KFactorQCDqqZZ_M = 1;
			//weight=1;
			weight= (xsec*KFactorEWKqqZZ * overallEventWeight*KFactorQCDqqZZ_M*lumi)/(gen_sum_weights);
			//if (j==1) weight= (xsec*KFactorQCDggzz_Nominal *lumi)/(gen_sum_weights);
			//weight= (xsec*KFactorEWKqqZZ * overallEventWeight*KFactorQCDqqZZ_M*lumi)/(gen_sum_weights);
			//dbkg_kin=overallEventWeight;
			//kin variable
			dbkg_kin =p_JJVBF_BKG_MCFM_JECNominal/(p_JJVBF_BKG_MCFM_JECNominal+ p_JJQCD_BKG_MCFM_JECNominal*0.02);

			//division in channels
			 if(abs(Z1Flav)==abs(Z2Flav) && abs(Z1Flav)==121){
                                 chan=2;
                        }
                        else if (abs(Z1Flav)==abs(Z2Flav) && abs(Z1Flav)!=121){
                                 chan=1;
                        }
                        else{
                                 chan=3;
                        }


			//1D kin var hist fill
			//this is the normalization histogram
			if (j!=3){
			           if (dbkg_kin <= 1.1) h00->Fill(dbkg_kin,weight); 				  				 h_complete_data->Fill(dbkg_kin,weight);
			}

			if (j==3){
				   if (dbkg_kin <= 0.7)   h0->Fill(dbkg_kin);
				   h0->SetMarkerStyle(20);
				   if (chan == 2 ) h0_ee->Fill(dbkg_kin);
				   if (chan == 1) h0_mm->Fill(dbkg_kin);
			           if (chan == 3) h0_em->Fill(dbkg_kin);
			}
			if (j==0){                              //qqzz
				    h1->Fill(dbkg_kin,weight);

				    if (chan == 1) h8->Fill(dbkg_kin,weight); //mu
				    if (chan == 2) h7->Fill(dbkg_kin,weight); //e
				    if (chan == 3)  h9->Fill(dbkg_kin,weight);          //mu+e
			}
 			if (j==1){                              //gg
                                    h2->Fill(dbkg_kin,weight);
                                    h2->SetFillColor(kBlue);

			            if (chan == 1) h11->Fill(dbkg_kin,weight); //mu
                                    if (chan == 2) h10->Fill(dbkg_kin,weight); //e
                                    if (chan == 3)  h12->Fill(dbkg_kin,weight);          //mu+e

                        }
			if (j==2){                              //vbs
                                    h3->Fill(dbkg_kin,weight);
                                    h3->SetFillColor(kMagenta);

				    if (chan == 1) h14->Fill(dbkg_kin,weight); //mu
                                    if (chan == 2) h13->Fill(dbkg_kin,weight); //e
                                    if (chan == 3) h15->Fill(dbkg_kin,weight);          //mu+e

			}
			//kin_zz->GetYaxis()->SetTitle("Events/0.05");
			//add histogram to stack
			//for cycle ends here
			}

		}//entries loop  end

	}//file loop  end

	//ZX CONTRIBUTION

       	 	TChain *tqqzz_zx= new TChain("candTree");
        	tqqzz_zx->Add("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/CMSSW_8_0_24_patch1/src/HZZ4l-plotter/ZXinput_mjj.root");

                //histogram declaration
                TH1F *kin_zz_zx = new TH1F("kin_zz","",20,0.,1.); //full histogram

		//new variable declarations
                float dbkg_kin_zx;
                float weight_zx;

                //original brach addresses
                tqqzz_zx->SetBranchAddress("dbkg_kin",&dbkg_kin_zx);
                tqqzz_zx->SetBranchAddress("weight",&weight_zx);

                //entries loop
                for(int i=0;i<tqqzz_zx->GetEntries();i++){
                tqqzz_zx->GetEntry(i);

                //1D kin var hist fill
                kin_zz_zx->Fill(dbkg_kin_zx,weight_zx); //full histogram filling
		if (dbkg_kin_zx<=0.7) h00->Fill(dbkg_kin_zx,weight_zx);
	        }

		       	 //INTEGRAL CHECK

			 cout << "Integral check" << endl;
			 cout <<"qqzz,         integral is " << h1->Integral() << endl;
			 cout <<"qqzz (4e),    integral is " << h7->Integral() << endl;
			 cout <<"qqzz (4mu),   integral is " << h8->Integral() << endl;
			 cout <<"qqzz (2e2mu), integral is " << h9->Integral() << endl;
			 cout <<"ggzz,         integral is " << h2->Integral() << endl;
			 cout <<"ggzz (4e),    integral is " << h10->Integral() << endl;
                         cout <<"ggzz (4mu),   integral is " << h11->Integral() << endl;
                         cout <<"ggzz (2e2mu), integral is " << h12->Integral() << endl;
			 cout <<"vbs,          integral is " << h3->Integral() << endl;
			 cout <<"vbs (4e),     integral is " << h13->Integral() << endl;
                         cout <<"vbs (4mu),    integral is " << h14->Integral() << endl;
                         cout <<"vbs (2e2mu),  integral is " << h15->Integral() << endl;
                         cout <<"data (4e),     integral is " << h0_ee->Integral() << endl;
                         cout <<"data (4mu),    integral is " << h0_mm->Integral() << endl;
                         cout <<"data (2e2mu),  integral is " << h0_em->Integral() << endl;

			 float workspace_integral = (2.00933+1.97468+3.96839);
                         kin_zz_zx->Scale(workspace_integral/kin_zz_zx->Integral());

			 //HISTOGRAMS ADDED TO STACK
			 kin_zz_zx->SetFillColor(kGreen);
			 h1bis->Add(kin_zz_zx,h1,1,1); //real ew
		         h1bis->SetFillColor(kCyan);
			 h4->Add(h1bis,h2,1,1); //real gg
               		 h4->SetFillColor(kBlue);
               		 h5->Add(h4,h3,1,1);    //real vbs
               		 h5->SetFillColor(kMagenta);

			 //data normalisation
			 data_integral = 1;// h5->Integral();
			 mc_integral = 1; //h_complete_data->Integral();
			 kin_zz_zx->Scale(data_integral/mc_integral);
			 h1bis->Scale(data_integral/mc_integral);
			 h4->Scale(data_integral/mc_integral);
			 h5->Scale(data_integral/mc_integral);
			 //add histograms to stack
               		 hs->Add(h5,"hist");
               		 hs->Add(h4,"hist");
               		 hs->Add(h1bis,"hist");
			 hs->Add(kin_zz_zx,"hist");
			 TH1F *h0divide = (TH1F*)h0->Clone();
			 hs->Add(h0,"E1");

                         // draw the legend
                         legend->SetTextFont(72);
                         legend->SetTextSize(0.04);
                         legend->AddEntry(kin_zz_zx,"Z+X","f");
			 legend->AddEntry(h1bis,"q#bar{q}#rightarrowZZ","f");
			 legend->AddEntry(h4,"gg#rightarrowZZ","f");
			 legend->AddEntry(h5,"VBS","f");
			 legend->AddEntry(h0,"Data","lep");
		         legend->SetBorderSize(0);

	//GRAPHICS
	//hs->GetXaxis()->SetTitle("K_{D}");
	//hs->GetXaxis()->SetTitleSize(.15);
        //hs->GetYaxis()->SetTitle("Events/0.05");
        //hs->GetYaxis()->SetTitleSize(.15);

        //RATIO PLOT
        //canvas
	TCanvas *c1 = new TCanvas("c1","example",800,1000);
        //pad1
        float eps =0;// 0.006;
	TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1,0);
	pad1->SetBottomMargin(0);
	pad1->Draw();
	pad1->cd();
	//top plot
	hs->Draw("nostack"); //old
	hs->GetXaxis()->SetTitle("K_{D}");
	legend->Draw("same");
	//switch?
	c1->cd();
	//pad2
	TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.3+eps,0); //old position
	pad2->SetTopMargin(0);
	//make bottom pad transparent
	//pad2->SetFrameFillColor(0);
	//pad2->SetFrameBorderMode(0);
	//pad2->SetFrameFillColor(0);
	//pad2->SetFrameBorderMode(0);

	pad2->Draw();
	pad2->cd();
	//bottom plot
	TH1F *h0copy = (TH1F*) h0divide->Clone();
	//axis labels
	h0copy->GetXaxis()->SetLabelFont(59);//change this for font type
	h0copy->GetXaxis()->SetLabelSize(22);
	h0copy->GetYaxis()->SetLabelFont(59);//change this for font type
	h0copy->GetYaxis()->SetLabelSize(22);
	//axis titles
	h0copy->GetXaxis()->SetTitleFont(59); //change this for font type
        h0copy->GetXaxis()->SetTitleSize(22);
        h0copy->GetYaxis()->SetTitleFont(59); //change this for font type
        h0copy->GetYaxis()->SetTitleSize(22);
	h0copy->GetXaxis()->SetTitleOffset(4.5);
	h0copy->GetYaxis()->SetTitleOffset(1.7);
	h0copy->GetYaxis()->SetRangeUser(-1.,3.);

	h0copy->Sumw2();
	h0copy->SetStats(0); //clear stat box
	h0copy->Divide(h5); //invert divide
	h0copy->SetMarkerStyle(20);
	h0copy->SetTitle("; K_{D}; Data/MC");
	//h0copy->GetXaxis()->SetTitleSize(50);
	//h0copy->GetXaxis()->SetLabelSize(35);
	//h0copy->GetYaxis()->SetTitleSize(15);
	//h0copy->GetYaxis()->SetLabelSize(15);

	//gStyle->SetLabelSize(2,"x").
	h0copy->Draw("ep");

	//Orizontal line
  	TLine *line = new TLine(0,1,1,1);
  	line->SetLineColor(kRed);
	line->SetLineStyle(2);
  	line->Draw("same");

        //close and print on file
	c1->cd();
	gPad->Print("template/plots/kin_variable_plot_new.png");
}
