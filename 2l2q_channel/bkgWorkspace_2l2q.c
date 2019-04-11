#include "ProcessNormalization.cc"

using namespace RooFit;

// highmass = 0 onshell
// highmass = 1 offshell
// highmass = 2 search
//
//
// cate_vbf =0 ggH
// cate_vbf =1 VBF
// cate_vbf =2 RSE
// cate_vbf =3 TLE

void dosomething(TString chan="all", int cate_vbf=1, int highmass=0,int is2D=0){ //was cate_vbf=1


	//luminosity
	double lumi = 35.9;

  //constants
	double parzx_all[5][6]={
		//		ee
		1,141.9,21.3,0,0,0,
		//mm
		1,130.4,15.6,0,0,0,
		//all
		0.45,131.1,18.1,0.55,133.8,18.9,
		//ee RSE
		0.298623315141,238.804039513,38.5525989967,4.83114145422,-0.0097489713697,0,
		//all RSE
		0.000171310428786,209.221006175,26.5346636174,11.193766044,-0.00296426709129,0
	};
	double parzx_rse_all[2][8]{
		//ee
		2.18216e+01  ,2.99737e+02   ,-2.48222e+02   ,1.38710e+01   ,2.76729e-03   ,2.02956e+02   ,2.54751e+01   ,4.84386e-05 ,
		//all RSE
		2.19408e+00   ,-1.25292e+02  ,-8.72399e+01   ,1.50523e+01   ,8.52961e-03   ,2.40365e+02   ,4.26689e+01   ,6.19292e-05
	};
	double parzx[6]={0.};
	double parzx_rse[8]={0.};

  //deal with constants sepaeretly for each channel
	if(cate_vbf!=2){
		if (chan=="ee") 	{
			for (int i=0;i<6;i++){parzx[i]=parzx_all[0][i];}
		}
		if (chan=="mm")  {
			for (int i=0;i<6;i++){parzx[i]=parzx_all[1][i];}
		}
		if (chan=="all") {
			for (int i=0;i<6;i++){parzx[i]=parzx_all[2][i];}
		}
	}
	else{
		if (chan=="all") {
			for (int i=0;i<8;i++){parzx_rse[i]=parzx_rse_all[1][i];}
		}
		if (chan=="ee") {
			for (int i=0;i<8;i++){parzx_rse[i]=parzx_rse_all[0][i];}
		}
		}
		//write yields file
		ofstream yields("yields.txt",std::fstream::app);

    //RooWorkspace declaration
		RooWorkspace w("w");

    //binning choices
  	double recolowarr[3]={500,500,500};
		double recohigharr[3]={3500.,3500.,3500.};
		const int reconbinsarr[3]={1500,1500,1500};


		const double low_reco=recolowarr[cate_vbf];
		const double high_reco=recohigharr[cate_vbf];
		const int nbins_reco=reconbinsarr[cate_vbf];




    //Roo variable mreco
		RooRealVar* TT_M_mreco= new RooRealVar("mreco","M_{ZZ} (GeV)",1500,low_reco,high_reco);//working with 180
		if(cate_vbf==2){
			TT_M_mreco->SetName("mreco_rse");
			TT_M_mreco->SetTitle("mreco_rse");
		}

    //Roo variable dbkg
		//RooRealVar* dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",0.5,0.,1.);
		RooRealVar* TT_M_dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",30,0.,1.);
		RooPlot* TT_M_frame= TT_M_mreco->frame(low_reco,high_reco) ;
		//RooPlot* frame1= dbkg->frame(0,1) ; //???

		if(highmass==2){
			TT_M_mreco->setBins(nbins_reco);
			TT_M_dbkg->setBins(30);
		}

		//tree management
		TString TT_M_treename[4]={"","","_rse","_tle"};
		TChain *TT_M_tTT_merged= new TChain("SelectedTree"+TT_M_treename[cate_vbf]);

		TH1F *TT_M_hTT_merged= new TH1F ("TT_M_hTT_merged","",nbins_reco,low_reco,high_reco);
		TH1F *TT_M_hTT_merged_up= new TH1F ("TT_M_hTT_merged_up","",nbins_reco,low_reco,high_reco);
		TH1F *TT_M_hTT_merged_dn= new TH1F ("TT_M_hTT_merged_dn","",nbins_reco,low_reco,high_reco);

		TT_M_tTT_merged->Add("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_TT_merged.root");

    //a code is assigned to each channel
		int TT_M_channum=1;
		if(chan=="all")
			TT_M_channum=1;
		else if(chan=="ee")
			TT_M_channum=0;
		else
			TT_M_channum=2;

    //extra Roo variables and normal variables
		// RooRealVar* wt= new RooRealVar("weight","wt",1500,0.000,1000.); //125
		// RooRealVar* wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,1000.);//125
		// RooRealVar* wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,1000.); //125
		RooRealVar* TT_M_wt= new RooRealVar("weight","wt",1500,0.000,5000.); //125
		RooRealVar* TT_M_wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,5000.);//125
		RooRealVar* TT_M_wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,5000.); //125

		float TT_M_ZZMass;
		float TT_M_weight, TT_M_weight_up, TT_M_weight_dn;
		int TT_M_channel;
		TT_M_tTT_merged->SetBranchAddress("mreco",&TT_M_ZZMass);
		//tTT_merged->SetBranchAddress("dbkg",&ZZMass); //???

		if(cate_vbf!=1){
			TT_M_tTT_merged->SetBranchAddress("weight",&TT_M_weight);
			TT_M_tTT_merged->SetBranchAddress("weight_up",&TT_M_weight_up);
			TT_M_tTT_merged->SetBranchAddress("weight_dn",&TT_M_weight_dn);
			TT_M_tTT_merged->Draw("mreco>>hTT_merged",Form("weight*(chan==%d)",TT_M_channum));
			TT_M_tTT_merged->Draw("mreco>>hTT_merged_dn",Form("weight_dn*(chan==%d)",TT_M_channum));
			TT_M_tTT_merged->Draw("mreco>>hTT_merged_up",Form("weight_up*(chan==%d)",TT_M_channum));
		}
		else{
			TT_M_tTT_merged->SetBranchAddress("weight_vbf",&TT_M_weight);
			TT_M_tTT_merged->SetBranchAddress("weight_vbf_up",&TT_M_weight_up);
			TT_M_tTT_merged->SetBranchAddress("weight_vbf_dn",&TT_M_weight_dn);
			TT_M_tTT_merged->Draw("mreco>>hTT_merged",Form("weight_vbf*(chan==%d)",TT_M_channum));
			TT_M_tTT_merged->Draw("mreco>>hTT_merged_dn",Form("weight_vbf_dn*(chan==%d)",TT_M_channum));
			TT_M_tTT_merged->Draw("mreco>>hTT_merged_up",Form("weight_vbf_up*(chan==%d)",TT_M_channum));
		}

		TT_M_tTT_merged->SetBranchAddress("chan",&TT_M_channel);

    // //integral output
		// if(cate_vbf!=2){
		// 	cout<<"test"<<endl;
		// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<TT_M_hTT_merged->Integral()<<endl;
		// 	cout<<"integral is " << hTT_merged->Integral()<<endl;
		// 	cout<<"lumi is" << lumi << endl;
		// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hTT_merged_up->Integral()<<endl;
		// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hTT_merged_dn->Integral()<<endl;
		// }
		// else{
		// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hTT_merged->Integral()<<endl;  /// Remember why 1.6?
		// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hTT_merged_up->Integral()<<endl;
		// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hTT_merged_dn->Integral()<<endl;
		// }

    //new tree
		TTree *TT_M_cuttree =new TTree("SelectedTree","SelectedTree");
		TString TT_M_mrecob = Form("%s/F",TT_M_mreco->GetName());
		TT_M_cuttree->Branch(TT_M_mreco->GetName(),&TT_M_ZZMass,TT_M_mrecob);
		TT_M_cuttree->Branch("weight",&TT_M_weight,"weight/F");
		TT_M_cuttree->Branch("weight_up",&TT_M_weight_up,"weight_up/F");
		TT_M_cuttree->Branch("weight_dn",&TT_M_weight_dn,"weight_dn/F");
		//std::cout<< endl << "PROBLEM" << endl;
		for(int i =0;i<TT_M_tTT_merged->GetEntries();i++){
			TT_M_tTT_merged->GetEntry(i);
			if(TT_M_channel==TT_M_channum){
				TT_M_cuttree->Fill();
			}
		}

		//Roo objects declared and filled
		double TT_M_rho =1;
		if(highmass==0)
			TT_M_rho=4;
		RooDataSet TT_M_bkgdata ("TT_M_bkgdata"+chan+Form("_%d",cate_vbf),"",TT_M_cuttree,RooArgSet(*TT_M_mreco,*TT_M_wt),"weight");
		RooKeysPdf TT_M_TT_mergedpdf_1d("TT_M_merged_1d","",*TT_M_mreco,TT_M_bkgdata,RooKeysPdf::MirrorBoth,TT_M_rho);

    //file input
		TFile *TT_M_ff = new TFile("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_TT_merged.root");
    //2D templated created and check for empty or NaN bins
		TH2F *TT_M_temp_zz=(TH2F*)TT_M_ff->Get("temp_rebin_"+chan+"_mergedSRvbf");

		//Preparazione METODO RICERCA ERRORI
		Double_t error_up_sum = 0;
		Double_t error_down_sum = 0;
		Double_t integral = 0;
		for(int bx=0;bx<=TT_M_temp_zz->GetNbinsX();bx++){
			for(int by=0;by<=TT_M_temp_zz->GetNbinsY();by++){

				if(!(TT_M_temp_zz->GetBinContent(bx,by) >0 || TT_M_temp_zz->GetBinContent(bx,by)<0)	){
					TT_M_temp_zz->SetBinContent(bx,by,1.0e-10);
				}

				//IMPLEMENTAZIONE METODO RICERCA ERRORI UP E DOWN
				integral = integral + TT_M_temp_zz->GetBinContent(bx,by);
				if (TT_M_temp_zz->GetBinContent(bx,by) > TT_M_temp_zz->GetBinErrorLow(bx,by)){
					  error_down_sum = error_down_sum + TT_M_temp_zz->GetBinErrorLow(bx,by);
				}
				error_up_sum = error_up_sum + TT_M_temp_zz->GetBinErrorUp(bx,by);

					//std::cout << endl << temp_zz->GetBinContent(bx,by) << " at " << bx << " "<< by <<endl;
			}
		}

		//ULTIMAZIONE METODO RICERCA ERRORI
		Double_t ratio_down = (-error_down_sum+integral)/integral;
		Double_t ratio_up = (error_up_sum+integral)/integral;
		std::cout << endl << "ratio down = " << ratio_down << endl << " ratio up = " << ratio_up << endl;

		 if (TT_M_channum == 1) {
		 	TT_M_temp_zz->Draw("colz");
		 	gPad->Print("workspace/template_check.png");
		 }

		//other Roo objects
		RooDataHist* TT_M_template_sig= new RooDataHist("TT_M_temp_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*TT_M_mreco,*TT_M_dbkg),TT_M_temp_zz);
		RooHistPdf* TT_M_pdf_2d_sig = new RooHistPdf("TT_M_pdf_2d_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*TT_M_mreco,*TT_M_dbkg),*TT_M_template_sig);
		RooProdPdf *TT_M_TT_mergedpdf_2d= new RooProdPdf("TT_M_merged_2d_"+chan+Form("_%d",cate_vbf),"",TT_M_TT_mergedpdf_1d,Conditional(*TT_M_pdf_2d_sig,*TT_M_dbkg));
// *********************************************************************************************************************************************************


//Roo variable mreco
RooRealVar* DY_M_mreco= new RooRealVar("mreco","M_{ZZ} (GeV)",1500,low_reco,high_reco);//working with 180
if(cate_vbf==2){
	DY_M_mreco->SetName("mreco_rse");
	DY_M_mreco->SetTitle("mreco_rse");
}

//Roo variable dbkg
//RooRealVar* dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",0.5,0.,1.);
RooRealVar* DY_M_dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",30,0.,1.);
RooPlot* DY_M_frame= DY_M_mreco->frame(low_reco,high_reco) ;
//RooPlot* frame1= dbkg->frame(0,1) ; //???

if(highmass==2){
	DY_M_mreco->setBins(nbins_reco);
	DY_M_dbkg->setBins(30);
}

//tree management
TString DY_M_treename[4]={"","","_rse","_tle"};
TChain *DY_M_tDY_merged= new TChain("SelectedTree"+DY_M_treename[cate_vbf]);

TH1F *DY_M_hDY_merged= new TH1F ("DY_M_hDY_merged","",nbins_reco,low_reco,high_reco);
TH1F *DY_M_hDY_merged_up= new TH1F ("DY_M_hDY_merged_up","",nbins_reco,low_reco,high_reco);
TH1F *DY_M_hDY_merged_dn= new TH1F ("DY_M_hDY_merged_dn","",nbins_reco,low_reco,high_reco);

DY_M_tDY_merged->Add("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_DY_merged.root");

//a code is assigned to each channel
int DY_M_channum=1;
if(chan=="all")
	DY_M_channum=1;
else if(chan=="ee")
	DY_M_channum=0;
else
	DY_M_channum=2;

//extra Roo variables and normal variables
// RooRealVar* wt= new RooRealVar("weight","wt",1500,0.000,1000.); //125
// RooRealVar* wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,1000.);//125
// RooRealVar* wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,1000.); //125
RooRealVar* DY_M_wt= new RooRealVar("weight","wt",1500,0.000,5000.); //125
RooRealVar* DY_M_wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,5000.);//125
RooRealVar* DY_M_wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,5000.); //125

float DY_M_ZZMass;
float DY_M_weight, DY_M_weight_up, DY_M_weight_dn;
int DY_M_channel;
DY_M_tDY_merged->SetBranchAddress("mreco",&DY_M_ZZMass);
//tDY_merged->SetBranchAddress("dbkg",&ZZMass); //???

if(cate_vbf!=1){
	DY_M_tDY_merged->SetBranchAddress("weight",&DY_M_weight);
	DY_M_tDY_merged->SetBranchAddress("weight_up",&DY_M_weight_up);
	DY_M_tDY_merged->SetBranchAddress("weight_dn",&DY_M_weight_dn);
	DY_M_tDY_merged->Draw("mreco>>hDY_merged",Form("weight*(chan==%d)",DY_M_channum));
	DY_M_tDY_merged->Draw("mreco>>hDY_merged_dn",Form("weight_dn*(chan==%d)",DY_M_channum));
	DY_M_tDY_merged->Draw("mreco>>hDY_merged_up",Form("weight_up*(chan==%d)",DY_M_channum));
}
else{
	DY_M_tDY_merged->SetBranchAddress("weight_vbf",&DY_M_weight);
	DY_M_tDY_merged->SetBranchAddress("weight_vbf_up",&DY_M_weight_up);
	DY_M_tDY_merged->SetBranchAddress("weight_vbf_dn",&DY_M_weight_dn);
	DY_M_tDY_merged->Draw("mreco>>hDY_merged",Form("weight_vbf*(chan==%d)",DY_M_channum));
	DY_M_tDY_merged->Draw("mreco>>hDY_merged_dn",Form("weight_vbf_dn*(chan==%d)",DY_M_channum));
	DY_M_tDY_merged->Draw("mreco>>hDY_merged_up",Form("weight_vbf_up*(chan==%d)",DY_M_channum));
}

DY_M_tDY_merged->SetBranchAddress("chan",&DY_M_channel);

// //integral output
// if(cate_vbf!=2){
// 	cout<<"test"<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<DY_M_hDY_merged->Integral()<<endl;
// 	cout<<"integral is " << hDY_merged->Integral()<<endl;
// 	cout<<"lumi is" << lumi << endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hDY_merged_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hDY_merged_dn->Integral()<<endl;
// }
// else{
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hDY_merged->Integral()<<endl;  /// Remember why 1.6?
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hDY_merged_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hDY_merged_dn->Integral()<<endl;
// }

//new tree
TTree *DY_M_cuttree =new TTree("SelectedTree","SelectedTree");
TString DY_M_mrecob = Form("%s/F",DY_M_mreco->GetName());
DY_M_cuttree->Branch(DY_M_mreco->GetName(),&DY_M_ZZMass,DY_M_mrecob);
DY_M_cuttree->Branch("weight",&DY_M_weight,"weight/F");
DY_M_cuttree->Branch("weight_up",&DY_M_weight_up,"weight_up/F");
DY_M_cuttree->Branch("weight_dn",&DY_M_weight_dn,"weight_dn/F");
//std::cout<< endl << "PROBLEM" << endl;
for(int i =0;i<DY_M_tDY_merged->GetEntries();i++){
	DY_M_tDY_merged->GetEntry(i);
	if(DY_M_channel==DY_M_channum){
		DY_M_cuttree->Fill();
	}
}

//Roo objects declared and filled
double DY_M_rho =1;
if(highmass==0)
	DY_M_rho=4;
RooDataSet DY_M_bkgdata ("DY_M_bkgdata"+chan+Form("_%d",cate_vbf),"",DY_M_cuttree,RooArgSet(*DY_M_mreco,*DY_M_wt),"weight");
RooKeysPdf DY_M_DY_mergedpdf_1d("DY_M_merged_1d","",*DY_M_mreco,DY_M_bkgdata,RooKeysPdf::MirrorBoth,DY_M_rho);

//file input
TFile *DY_M_ff = new TFile("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_DY_merged.root");
//2D templated created and check for empty or NaN bins
TH2F *DY_M_temp_zz=(TH2F*)DY_M_ff->Get("temp_rebin_"+chan+"_mergedSRvbf");

for(int bx=0;bx<=DY_M_temp_zz->GetNbinsX();bx++){
	for(int by=0;by<=DY_M_temp_zz->GetNbinsY();by++){

		if(!(DY_M_temp_zz->GetBinContent(bx,by) >0 || DY_M_temp_zz->GetBinContent(bx,by)<0)	){
			DY_M_temp_zz->SetBinContent(bx,by,1.0e-10);

		}
			//std::cout << endl << temp_zz->GetBinContent(bx,by) << " at " << bx << " "<< by <<endl;
	}
}
 if (DY_M_channum == 1) {
	gStyle->SetOptStat(0);
	TCanvas *c = new TCanvas();
  c->SetLeftMargin(0.2);
	c->SetRightMargin(0.2);
	c->SetBottomMargin(0.2);
	DY_M_temp_zz->SetZTitle("weight");
	DY_M_temp_zz->SetXTitle("m_{ZZ}");
	DY_M_temp_zz->SetYTitle("D_{k}");
	DY_M_temp_zz->Draw("colz");
	gPad->Print("workspace/DY_merged.png");
 }

//other Roo objects
RooDataHist* DY_M_template_sig= new RooDataHist("DY_M_temp_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*DY_M_mreco,*DY_M_dbkg),DY_M_temp_zz);
RooHistPdf* DY_M_pdf_2d_sig = new RooHistPdf("DY_M_pdf_2d_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*DY_M_mreco,*DY_M_dbkg),*DY_M_template_sig);
RooProdPdf *DY_M_DY_mergedpdf_2d= new RooProdPdf("DY_M_merged_2d_"+chan+Form("_%d",cate_vbf),"",DY_M_DY_mergedpdf_1d,Conditional(*DY_M_pdf_2d_sig,*DY_M_dbkg));
// *********************************************************************************************************************************************************


		//Data
		TChain *tdata = new TChain("SelectedTree"+TT_M_treename[cate_vbf]);
		tdata->Add("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/RooFitInput_Data_R.root");
		TTree* reducetree= tdata->CopyTree(Form("chan==%d&&vbfcate==%d",TT_M_channum,cate_vbf));
		RooDataSet* data_obs_1d= new RooDataSet("data_obs_1d","data_obs_1d",reducetree,*TT_M_mreco);
		RooDataSet* data_obs_2d = new RooDataSet("data_obs_2d","data_obs_2d",reducetree,RooArgSet(*TT_M_mreco,*TT_M_dbkg));

// *********************************************************************************************************************************************************

//Roo variable mreco
RooRealVar* WW_ZZ_M_mreco= new RooRealVar("mreco","M_{ZZ} (GeV)",1500,low_reco,high_reco);//working with 180
if(cate_vbf==2){
	WW_ZZ_M_mreco->SetName("mreco_rse");
	WW_ZZ_M_mreco->SetTitle("mreco_rse");
}

//Roo variable dbkg
//RooRealVar* dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",0.5,0.,1.);
RooRealVar* WW_ZZ_M_dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",30,0.,1.);
RooPlot* WW_ZZ_M_frame= WW_ZZ_M_mreco->frame(low_reco,high_reco) ;
//RooPlot* frame1= dbkg->frame(0,1) ; //???

if(highmass==2){
	WW_ZZ_M_mreco->setBins(nbins_reco);
	WW_ZZ_M_dbkg->setBins(30);
}

//tree management
TString WW_ZZ_M_treename[4]={"","","_rse","_tle"};
TChain *WW_ZZ_M_tWW_ZZ_merged= new TChain("SelectedTree"+WW_ZZ_M_treename[cate_vbf]);

TH1F *WW_ZZ_M_hWW_ZZ_merged= new TH1F ("WW_ZZ_M_hWW_ZZ_merged","",nbins_reco,low_reco,high_reco);
TH1F *WW_ZZ_M_hWW_ZZ_merged_up= new TH1F ("WW_ZZ_M_hWW_ZZ_merged_up","",nbins_reco,low_reco,high_reco);
TH1F *WW_ZZ_M_hWW_ZZ_merged_dn= new TH1F ("WW_ZZ_M_hWW_ZZ_merged_dn","",nbins_reco,low_reco,high_reco);

WW_ZZ_M_tWW_ZZ_merged->Add("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_WW_ZZ_merged.root");

//a code is assigned to each channel
int WW_ZZ_M_channum=1;
if(chan=="all")
	WW_ZZ_M_channum=1;
else if(chan=="ee")
	WW_ZZ_M_channum=0;
else
	WW_ZZ_M_channum=2;

//extra Roo variables and normal variables
// RooRealVar* wt= new RooRealVar("weight","wt",1500,0.000,1000.); //125
// RooRealVar* wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,1000.);//125
// RooRealVar* wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,1000.); //125
RooRealVar* WW_ZZ_M_wt= new RooRealVar("weight","wt",1500,0.000,5000.); //125
RooRealVar* WW_ZZ_M_wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,5000.);//125
RooRealVar* WW_ZZ_M_wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,5000.); //125

float WW_ZZ_M_ZZMass;
float WW_ZZ_M_weight, WW_ZZ_M_weight_up, WW_ZZ_M_weight_dn;
int WW_ZZ_M_channel;
WW_ZZ_M_tWW_ZZ_merged->SetBranchAddress("mreco",&WW_ZZ_M_ZZMass);
//tWW_ZZ_merged->SetBranchAddress("dbkg",&ZZMass); //???

if(cate_vbf!=1){
	WW_ZZ_M_tWW_ZZ_merged->SetBranchAddress("weight",&WW_ZZ_M_weight);
	WW_ZZ_M_tWW_ZZ_merged->SetBranchAddress("weight_up",&WW_ZZ_M_weight_up);
	WW_ZZ_M_tWW_ZZ_merged->SetBranchAddress("weight_dn",&WW_ZZ_M_weight_dn);
	WW_ZZ_M_tWW_ZZ_merged->Draw("mreco>>hWW_ZZ_merged",Form("weight*(chan==%d)",WW_ZZ_M_channum));
	WW_ZZ_M_tWW_ZZ_merged->Draw("mreco>>hWW_ZZ_merged_dn",Form("weight_dn*(chan==%d)",WW_ZZ_M_channum));
	WW_ZZ_M_tWW_ZZ_merged->Draw("mreco>>hWW_ZZ_merged_up",Form("weight_up*(chan==%d)",WW_ZZ_M_channum));
}
else{
	WW_ZZ_M_tWW_ZZ_merged->SetBranchAddress("weight_vbf",&WW_ZZ_M_weight);
	WW_ZZ_M_tWW_ZZ_merged->SetBranchAddress("weight_vbf_up",&WW_ZZ_M_weight_up);
	WW_ZZ_M_tWW_ZZ_merged->SetBranchAddress("weight_vbf_dn",&WW_ZZ_M_weight_dn);
	WW_ZZ_M_tWW_ZZ_merged->Draw("mreco>>hWW_ZZ_merged",Form("weight_vbf*(chan==%d)",WW_ZZ_M_channum));
	WW_ZZ_M_tWW_ZZ_merged->Draw("mreco>>hWW_ZZ_merged_dn",Form("weight_vbf_dn*(chan==%d)",WW_ZZ_M_channum));
	WW_ZZ_M_tWW_ZZ_merged->Draw("mreco>>hWW_ZZ_merged_up",Form("weight_vbf_up*(chan==%d)",WW_ZZ_M_channum));
}

WW_ZZ_M_tWW_ZZ_merged->SetBranchAddress("chan",&WW_ZZ_M_channel);

// //integral output
// if(cate_vbf!=2){
// 	cout<<"test"<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<WW_ZZ_M_hWW_ZZ_merged->Integral()<<endl;
// 	cout<<"integral is " << hWW_ZZ_merged->Integral()<<endl;
// 	cout<<"lumi is" << lumi << endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hWW_ZZ_merged_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hWW_ZZ_merged_dn->Integral()<<endl;
// }
// else{
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hWW_ZZ_merged->Integral()<<endl;  /// Remember why 1.6?
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hWW_ZZ_merged_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hWW_ZZ_merged_dn->Integral()<<endl;
// }

//new tree
TTree *WW_ZZ_M_cuttree =new TTree("SelectedTree","SelectedTree");
TString WW_ZZ_M_mrecob = Form("%s/F",WW_ZZ_M_mreco->GetName());
WW_ZZ_M_cuttree->Branch(WW_ZZ_M_mreco->GetName(),&WW_ZZ_M_ZZMass,WW_ZZ_M_mrecob);
WW_ZZ_M_cuttree->Branch("weight",&WW_ZZ_M_weight,"weight/F");
WW_ZZ_M_cuttree->Branch("weight_up",&WW_ZZ_M_weight_up,"weight_up/F");
WW_ZZ_M_cuttree->Branch("weight_dn",&WW_ZZ_M_weight_dn,"weight_dn/F");
//std::cout<< endl << "PROBLEM" << endl;
for(int i =0;i<WW_ZZ_M_tWW_ZZ_merged->GetEntries();i++){
	WW_ZZ_M_tWW_ZZ_merged->GetEntry(i);
	if(WW_ZZ_M_channel==WW_ZZ_M_channum){
		WW_ZZ_M_cuttree->Fill();
	}
}

//Roo objects declared and filled
double WW_ZZ_M_rho =1;
if(highmass==0)
	WW_ZZ_M_rho=4;
RooDataSet WW_ZZ_M_bkgdata ("WW_ZZ_M_bkgdata"+chan+Form("_%d",cate_vbf),"",WW_ZZ_M_cuttree,RooArgSet(*WW_ZZ_M_mreco,*WW_ZZ_M_wt),"weight");
RooKeysPdf WW_ZZ_M_WW_ZZ_mergedpdf_1d("WW_ZZ_M_merged_1d","",*WW_ZZ_M_mreco,WW_ZZ_M_bkgdata,RooKeysPdf::MirrorBoth,WW_ZZ_M_rho);

//file input
TFile *WW_ZZ_M_ff = new TFile("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_WW_ZZ_merged.root");
//2D templated created and check for empty or NaN bins
TH2F *WW_ZZ_M_temp_zz=(TH2F*)WW_ZZ_M_ff->Get("temp_rebin_"+chan+"_mergedSRvbf");

for(int bx=0;bx<=WW_ZZ_M_temp_zz->GetNbinsX();bx++){
	for(int by=0;by<=WW_ZZ_M_temp_zz->GetNbinsY();by++){

		if(!(WW_ZZ_M_temp_zz->GetBinContent(bx,by) >0 || WW_ZZ_M_temp_zz->GetBinContent(bx,by)<0)	){
			WW_ZZ_M_temp_zz->SetBinContent(bx,by,1.0e-10);

		}
			//std::cout << endl << temp_zz->GetBinContent(bx,by) << " at " << bx << " "<< by <<endl;
	}
}
 if (WW_ZZ_M_channum == 1) {
	WW_ZZ_M_temp_zz->Draw("colz");
	gPad->Print("workspace/template_check.png");
 }

//other Roo objects
RooDataHist* WW_ZZ_M_template_sig= new RooDataHist("WW_ZZ_M_temp_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*WW_ZZ_M_mreco,*WW_ZZ_M_dbkg),WW_ZZ_M_temp_zz);
RooHistPdf* WW_ZZ_M_pdf_2d_sig = new RooHistPdf("WW_ZZ_M_pdf_2d_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*WW_ZZ_M_mreco,*WW_ZZ_M_dbkg),*WW_ZZ_M_template_sig);
RooProdPdf *WW_ZZ_M_WW_ZZ_mergedpdf_2d= new RooProdPdf("WW_ZZ_M_merged_2d_"+chan+Form("_%d",cate_vbf),"",WW_ZZ_M_WW_ZZ_mergedpdf_1d,Conditional(*WW_ZZ_M_pdf_2d_sig,*WW_ZZ_M_dbkg));

//**************************************************************

//Roo variable mreco
RooRealVar* sig_WZ_M_mreco= new RooRealVar("mreco","M_{ZZ} (GeV)",1500,low_reco,high_reco);//working with 180
if(cate_vbf==2){
	sig_WZ_M_mreco->SetName("mreco_rse");
	sig_WZ_M_mreco->SetTitle("mreco_rse");
}

//Roo variable dbkg
//RooRealVar* dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",0.5,0.,1.);
RooRealVar* sig_WZ_M_dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",30,0.,1.);
RooPlot* sig_WZ_M_frame= sig_WZ_M_mreco->frame(low_reco,high_reco) ;
//RooPlot* frame1= dbkg->frame(0,1) ; //???

if(highmass==2){
	sig_WZ_M_mreco->setBins(nbins_reco);
	sig_WZ_M_dbkg->setBins(30);
}

//tree management
TString sig_WZ_M_treename[4]={"","","_rse","_tle"};
TChain *sig_WZ_M_tsig_WZ_merged= new TChain("SelectedTree"+sig_WZ_M_treename[cate_vbf]);

TH1F *sig_WZ_M_hsig_WZ_merged= new TH1F ("sig_WZ_M_hsig_WZ_merged","",nbins_reco,low_reco,high_reco);
TH1F *sig_WZ_M_hsig_WZ_merged_up= new TH1F ("sig_WZ_M_hsig_WZ_merged_up","",nbins_reco,low_reco,high_reco);
TH1F *sig_WZ_M_hsig_WZ_merged_dn= new TH1F ("sig_WZ_M_hsig_WZ_merged_dn","",nbins_reco,low_reco,high_reco);

sig_WZ_M_tsig_WZ_merged->Add("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_sig_WZ_merged.root");

//a code is assigned to each channel
int sig_WZ_M_channum=1;
if(chan=="all")
	sig_WZ_M_channum=1;
else if(chan=="ee")
	sig_WZ_M_channum=0;
else
	sig_WZ_M_channum=2;

//extra Roo variables and normal variables
// RooRealVar* wt= new RooRealVar("weight","wt",1500,0.000,1000.); //125
// RooRealVar* wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,1000.);//125
// RooRealVar* wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,1000.); //125
RooRealVar* sig_WZ_M_wt= new RooRealVar("weight","wt",1500,0.000,5000.); //125
RooRealVar* sig_WZ_M_wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,5000.);//125
RooRealVar* sig_WZ_M_wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,5000.); //125

float sig_WZ_M_ZZMass;
float sig_WZ_M_weight, sig_WZ_M_weight_up, sig_WZ_M_weight_dn;
int sig_WZ_M_channel;
sig_WZ_M_tsig_WZ_merged->SetBranchAddress("mreco",&sig_WZ_M_ZZMass);
//tsig_WZ_merged->SetBranchAddress("dbkg",&ZZMass); //???

if(cate_vbf!=1){
	sig_WZ_M_tsig_WZ_merged->SetBranchAddress("weight",&sig_WZ_M_weight);
	sig_WZ_M_tsig_WZ_merged->SetBranchAddress("weight_up",&sig_WZ_M_weight_up);
	sig_WZ_M_tsig_WZ_merged->SetBranchAddress("weight_dn",&sig_WZ_M_weight_dn);
	sig_WZ_M_tsig_WZ_merged->Draw("mreco>>hsig_WZ_merged",Form("weight*(chan==%d)",sig_WZ_M_channum));
	sig_WZ_M_tsig_WZ_merged->Draw("mreco>>hsig_WZ_merged_dn",Form("weight_dn*(chan==%d)",sig_WZ_M_channum));
	sig_WZ_M_tsig_WZ_merged->Draw("mreco>>hsig_WZ_merged_up",Form("weight_up*(chan==%d)",sig_WZ_M_channum));
}
else{
	sig_WZ_M_tsig_WZ_merged->SetBranchAddress("weight_vbf",&sig_WZ_M_weight);
	sig_WZ_M_tsig_WZ_merged->SetBranchAddress("weight_vbf_up",&sig_WZ_M_weight_up);
	sig_WZ_M_tsig_WZ_merged->SetBranchAddress("weight_vbf_dn",&sig_WZ_M_weight_dn);
	sig_WZ_M_tsig_WZ_merged->Draw("mreco>>hsig_WZ_merged",Form("weight_vbf*(chan==%d)",sig_WZ_M_channum));
	sig_WZ_M_tsig_WZ_merged->Draw("mreco>>hsig_WZ_merged_dn",Form("weight_vbf_dn*(chan==%d)",sig_WZ_M_channum));
	sig_WZ_M_tsig_WZ_merged->Draw("mreco>>hsig_WZ_merged_up",Form("weight_vbf_up*(chan==%d)",sig_WZ_M_channum));
}

sig_WZ_M_tsig_WZ_merged->SetBranchAddress("chan",&sig_WZ_M_channel);

// //integral output
// if(cate_vbf!=2){
// 	cout<<"test"<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<sig_WZ_M_hsig_WZ_merged->Integral()<<endl;
// 	cout<<"integral is " << hsig_WZ_merged->Integral()<<endl;
// 	cout<<"lumi is" << lumi << endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsig_WZ_merged_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsig_WZ_merged_dn->Integral()<<endl;
// }
// else{
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsig_WZ_merged->Integral()<<endl;  /// Remember why 1.6?
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsig_WZ_merged_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsig_WZ_merged_dn->Integral()<<endl;
// }

//new tree
TTree *sig_WZ_M_cuttree =new TTree("SelectedTree","SelectedTree");
TString sig_WZ_M_mrecob = Form("%s/F",sig_WZ_M_mreco->GetName());
sig_WZ_M_cuttree->Branch(sig_WZ_M_mreco->GetName(),&sig_WZ_M_ZZMass,sig_WZ_M_mrecob);
sig_WZ_M_cuttree->Branch("weight",&sig_WZ_M_weight,"weight/F");
sig_WZ_M_cuttree->Branch("weight_up",&sig_WZ_M_weight_up,"weight_up/F");
sig_WZ_M_cuttree->Branch("weight_dn",&sig_WZ_M_weight_dn,"weight_dn/F");
//std::cout<< endl << "PROBLEM" << endl;
for(int i =0;i<sig_WZ_M_tsig_WZ_merged->GetEntries();i++){
	sig_WZ_M_tsig_WZ_merged->GetEntry(i);
	if(sig_WZ_M_channel==sig_WZ_M_channum){
		sig_WZ_M_cuttree->Fill();
	}
}

//Roo objects declared and filled
double sig_WZ_M_rho =1;
if(highmass==0)
	sig_WZ_M_rho=4;
RooDataSet sig_WZ_M_bkgdata ("sig_WZ_M_bkgdata"+chan+Form("_%d",cate_vbf),"",sig_WZ_M_cuttree,RooArgSet(*sig_WZ_M_mreco,*sig_WZ_M_wt),"weight");
RooKeysPdf sig_WZ_M_sig_WZ_mergedpdf_1d("sig_WZ_M_merged_1d","",*sig_WZ_M_mreco,sig_WZ_M_bkgdata,RooKeysPdf::MirrorBoth,sig_WZ_M_rho);

//file input
TFile *sig_WZ_M_ff = new TFile("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_sig_WZ_merged.root");
//2D templated created and check for empty or NaN bins
TH2F *sig_WZ_M_temp_zz=(TH2F*)sig_WZ_M_ff->Get("temp_rebin_"+chan+"_mergedSRvbf");

for(int bx=0;bx<=sig_WZ_M_temp_zz->GetNbinsX();bx++){
	for(int by=0;by<=sig_WZ_M_temp_zz->GetNbinsY();by++){

		if(!(sig_WZ_M_temp_zz->GetBinContent(bx,by) >0 || sig_WZ_M_temp_zz->GetBinContent(bx,by)<0)	){
			sig_WZ_M_temp_zz->SetBinContent(bx,by,1.0e-10);

		}
			//std::cout << endl << temp_zz->GetBinContent(bx,by) << " at " << bx << " "<< by <<endl;
	}
}
 if (sig_WZ_M_channum == 1) {
	sig_WZ_M_temp_zz->Draw("colz");
	gPad->Print("workspace/template_check.png");
 }

//other Roo objects
RooDataHist* sig_WZ_M_template_sig= new RooDataHist("sig_WZ_M_temp_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*sig_WZ_M_mreco,*sig_WZ_M_dbkg),sig_WZ_M_temp_zz);
RooHistPdf* sig_WZ_M_pdf_2d_sig = new RooHistPdf("sig_WZ_M_pdf_2d_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*sig_WZ_M_mreco,*sig_WZ_M_dbkg),*sig_WZ_M_template_sig);
RooProdPdf *sig_WZ_M_sig_WZ_mergedpdf_2d= new RooProdPdf("sig_WZ_M_merged_2d_"+chan+Form("_%d",cate_vbf),"",sig_WZ_M_sig_WZ_mergedpdf_1d,Conditional(*sig_WZ_M_pdf_2d_sig,*sig_WZ_M_dbkg));

//***************************************************************

//Roo variable mreco
RooRealVar* sig_ZZ_M_mreco= new RooRealVar("mreco","M_{ZZ} (GeV)",1500,low_reco,high_reco);//working with 180
if(cate_vbf==2){
	sig_ZZ_M_mreco->SetName("mreco_rse");
	sig_ZZ_M_mreco->SetTitle("mreco_rse");
}

//Roo variable dbkg
//RooRealVar* dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",0.5,0.,1.);
RooRealVar* sig_ZZ_M_dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",30,0.,1.);
RooPlot* sig_ZZ_M_frame= sig_ZZ_M_mreco->frame(low_reco,high_reco) ;
//RooPlot* frame1= dbkg->frame(0,1) ; //???

if(highmass==2){
	sig_ZZ_M_mreco->setBins(nbins_reco);
	sig_ZZ_M_dbkg->setBins(30);
}

//tree management
TString sig_ZZ_M_treename[4]={"","","_rse","_tle"};
TChain *sig_ZZ_M_tsig_ZZ_merged= new TChain("SelectedTree"+sig_ZZ_M_treename[cate_vbf]);

TH1F *sig_ZZ_M_hsig_ZZ_merged= new TH1F ("sig_ZZ_M_hsig_ZZ_merged","",nbins_reco,low_reco,high_reco);
TH1F *sig_ZZ_M_hsig_ZZ_merged_up= new TH1F ("sig_ZZ_M_hsig_ZZ_merged_up","",nbins_reco,low_reco,high_reco);
TH1F *sig_ZZ_M_hsig_ZZ_merged_dn= new TH1F ("sig_ZZ_M_hsig_ZZ_merged_dn","",nbins_reco,low_reco,high_reco);

sig_ZZ_M_tsig_ZZ_merged->Add("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_sig_ZZ_merged.root");

//a code is assigned to each channel
int sig_ZZ_M_channum=1;
if(chan=="all")
	sig_ZZ_M_channum=1;
else if(chan=="ee")
	sig_ZZ_M_channum=0;
else
	sig_ZZ_M_channum=2;

//extra Roo variables and normal variables
// RooRealVar* wt= new RooRealVar("weight","wt",1500,0.000,1000.); //125
// RooRealVar* wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,1000.);//125
// RooRealVar* wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,1000.); //125
RooRealVar* sig_ZZ_M_wt= new RooRealVar("weight","wt",1500,0.000,5000.); //125
RooRealVar* sig_ZZ_M_wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,5000.);//125
RooRealVar* sig_ZZ_M_wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,5000.); //125

float sig_ZZ_M_ZZMass;
float sig_ZZ_M_weight, sig_ZZ_M_weight_up, sig_ZZ_M_weight_dn;
int sig_ZZ_M_channel;
sig_ZZ_M_tsig_ZZ_merged->SetBranchAddress("mreco",&sig_ZZ_M_ZZMass);
//tsig_ZZ_merged->SetBranchAddress("dbkg",&ZZMass); //???

if(cate_vbf!=1){
	sig_ZZ_M_tsig_ZZ_merged->SetBranchAddress("weight",&sig_ZZ_M_weight);
	sig_ZZ_M_tsig_ZZ_merged->SetBranchAddress("weight_up",&sig_ZZ_M_weight_up);
	sig_ZZ_M_tsig_ZZ_merged->SetBranchAddress("weight_dn",&sig_ZZ_M_weight_dn);
	sig_ZZ_M_tsig_ZZ_merged->Draw("mreco>>hsig_ZZ_merged",Form("weight*(chan==%d)",sig_ZZ_M_channum));
	sig_ZZ_M_tsig_ZZ_merged->Draw("mreco>>hsig_ZZ_merged_dn",Form("weight_dn*(chan==%d)",sig_ZZ_M_channum));
	sig_ZZ_M_tsig_ZZ_merged->Draw("mreco>>hsig_ZZ_merged_up",Form("weight_up*(chan==%d)",sig_ZZ_M_channum));
}
else{
	sig_ZZ_M_tsig_ZZ_merged->SetBranchAddress("weight_vbf",&sig_ZZ_M_weight);
	sig_ZZ_M_tsig_ZZ_merged->SetBranchAddress("weight_vbf_up",&sig_ZZ_M_weight_up);
	sig_ZZ_M_tsig_ZZ_merged->SetBranchAddress("weight_vbf_dn",&sig_ZZ_M_weight_dn);
	sig_ZZ_M_tsig_ZZ_merged->Draw("mreco>>hsig_ZZ_merged",Form("weight_vbf*(chan==%d)",sig_ZZ_M_channum));
	sig_ZZ_M_tsig_ZZ_merged->Draw("mreco>>hsig_ZZ_merged_dn",Form("weight_vbf_dn*(chan==%d)",sig_ZZ_M_channum));
	sig_ZZ_M_tsig_ZZ_merged->Draw("mreco>>hsig_ZZ_merged_up",Form("weight_vbf_up*(chan==%d)",sig_ZZ_M_channum));
}

sig_ZZ_M_tsig_ZZ_merged->SetBranchAddress("chan",&sig_ZZ_M_channel);

// //integral output
// if(cate_vbf!=2){
// 	cout<<"test"<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<sig_ZZ_M_hsig_ZZ_merged->Integral()<<endl;
// 	cout<<"integral is " << hsig_ZZ_merged->Integral()<<endl;
// 	cout<<"lumi is" << lumi << endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsig_ZZ_merged_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsig_ZZ_merged_dn->Integral()<<endl;
// }
// else{
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsig_ZZ_merged->Integral()<<endl;  /// Remember why 1.6?
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsig_ZZ_merged_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsig_ZZ_merged_dn->Integral()<<endl;
// }

//new tree
TTree *sig_ZZ_M_cuttree =new TTree("SelectedTree","SelectedTree");
TString sig_ZZ_M_mrecob = Form("%s/F",sig_ZZ_M_mreco->GetName());
sig_ZZ_M_cuttree->Branch(sig_ZZ_M_mreco->GetName(),&sig_ZZ_M_ZZMass,sig_ZZ_M_mrecob);
sig_ZZ_M_cuttree->Branch("weight",&sig_ZZ_M_weight,"weight/F");
sig_ZZ_M_cuttree->Branch("weight_up",&sig_ZZ_M_weight_up,"weight_up/F");
sig_ZZ_M_cuttree->Branch("weight_dn",&sig_ZZ_M_weight_dn,"weight_dn/F");
//std::cout<< endl << "PROBLEM" << endl;
for(int i =0;i<sig_ZZ_M_tsig_ZZ_merged->GetEntries();i++){
	sig_ZZ_M_tsig_ZZ_merged->GetEntry(i);
	if(sig_ZZ_M_channel==sig_ZZ_M_channum){
		sig_ZZ_M_cuttree->Fill();
	}
}

//Roo objects declared and filled
double sig_ZZ_M_rho =1;
if(highmass==0)
	sig_ZZ_M_rho=4;
RooDataSet sig_ZZ_M_bkgdata ("sig_ZZ_M_bkgdata"+chan+Form("_%d",cate_vbf),"",sig_ZZ_M_cuttree,RooArgSet(*sig_ZZ_M_mreco,*sig_ZZ_M_wt),"weight");
RooKeysPdf sig_ZZ_M_sig_ZZ_mergedpdf_1d("sig_ZZ_M_merged_1d","",*sig_ZZ_M_mreco,sig_ZZ_M_bkgdata,RooKeysPdf::MirrorBoth,sig_ZZ_M_rho);

//file input
TFile *sig_ZZ_M_ff = new TFile("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_sig_ZZ_merged.root");
//2D templated created and check for empty or NaN bins
TH2F *sig_ZZ_M_temp_zz=(TH2F*)sig_ZZ_M_ff->Get("temp_rebin_"+chan+"_mergedSRvbf");

for(int bx=0;bx<=sig_ZZ_M_temp_zz->GetNbinsX();bx++){
	for(int by=0;by<=sig_ZZ_M_temp_zz->GetNbinsY();by++){

		if(!(sig_ZZ_M_temp_zz->GetBinContent(bx,by) >0 || sig_ZZ_M_temp_zz->GetBinContent(bx,by)<0)	){
			sig_ZZ_M_temp_zz->SetBinContent(bx,by,1.0e-10);

		}
			//std::cout << endl << temp_zz->GetBinContent(bx,by) << " at " << bx << " "<< by <<endl;
	}
}
 if (sig_ZZ_M_channum == 1) {
	gStyle->SetOptStat(0);
	sig_ZZ_M_temp_zz->Draw("colz");
	gPad->Print("workspace/sig_zz.png");
 }

//other Roo objects
RooDataHist* sig_ZZ_M_template_sig= new RooDataHist("sig_ZZ_M_temp_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*sig_ZZ_M_mreco,*sig_ZZ_M_dbkg),sig_ZZ_M_temp_zz);
RooHistPdf* sig_ZZ_M_pdf_2d_sig = new RooHistPdf("sig_ZZ_M_pdf_2d_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*sig_ZZ_M_mreco,*sig_ZZ_M_dbkg),*sig_ZZ_M_template_sig);
RooProdPdf *sig_ZZ_M_sig_ZZ_mergedpdf_2d= new RooProdPdf("sig_ZZ_M_merged_2d_"+chan+Form("_%d",cate_vbf),"",sig_ZZ_M_sig_ZZ_mergedpdf_1d,Conditional(*sig_ZZ_M_pdf_2d_sig,*sig_ZZ_M_dbkg));

//*************************************************************+

RooRealVar* TT_R_mreco= new RooRealVar("mreco","M_{ZZ} (GeV)",1500,low_reco,high_reco);//working with 180
if(cate_vbf==2){
  TT_R_mreco->SetName("mreco_rse");
  TT_R_mreco->SetTitle("mreco_rse");
}

//Roo variable dbkg
//RooRealVar* dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",0.5,0.,1.);
RooRealVar* TT_R_dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",30,0.,1.);
RooPlot* TT_R_frame= TT_R_mreco->frame(low_reco,high_reco) ;
//RooPlot* frame1= dbkg->frame(0,1) ; //???

if(highmass==2){
  TT_R_mreco->setBins(nbins_reco);
  TT_R_dbkg->setBins(30);
}

//tree management
TString TT_R_treename[4]={"","","_rse","_tle"};
TChain *TT_R_tTT_resolved= new TChain("SelectedTree"+TT_R_treename[cate_vbf]);

TH1F *TT_R_hTT_resolved= new TH1F ("TT_R_hTT_resolved","",nbins_reco,low_reco,high_reco);
TH1F *TT_R_hTT_resolved_up= new TH1F ("TT_R_hTT_resolved_up","",nbins_reco,low_reco,high_reco);
TH1F *TT_R_hTT_resolved_dn= new TH1F ("TT_R_hTT_resolved_dn","",nbins_reco,low_reco,high_reco);

TT_R_tTT_resolved->Add("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_TT_resolved.root");

//a code is assigned to each channel
int TT_R_channum=1;
if(chan=="all")
  TT_R_channum=1;
else if(chan=="ee")
  TT_R_channum=0;
else
  TT_R_channum=2;

//extra Roo variables and normal variables
// RooRealVar* wt= new RooRealVar("weight","wt",1500,0.000,1000.); //125
// RooRealVar* wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,1000.);//125
// RooRealVar* wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,1000.); //125
RooRealVar* TT_R_wt= new RooRealVar("weight","wt",1500,0.000,5000.); //125
RooRealVar* TT_R_wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,5000.);//125
RooRealVar* TT_R_wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,5000.); //125

float TT_R_ZZMass;
float TT_R_weight, TT_R_weight_up, TT_R_weight_dn;
int TT_R_channel;
TT_R_tTT_resolved->SetBranchAddress("mreco",&TT_R_ZZMass);
//tTT_resolved->SetBranchAddress("dbkg",&ZZMass); //???

if(cate_vbf!=1){
  TT_R_tTT_resolved->SetBranchAddress("weight",&TT_R_weight);
  TT_R_tTT_resolved->SetBranchAddress("weight_up",&TT_R_weight_up);
  TT_R_tTT_resolved->SetBranchAddress("weight_dn",&TT_R_weight_dn);
  TT_R_tTT_resolved->Draw("mreco>>hTT_resolved",Form("weight*(chan==%d)",TT_R_channum));
  TT_R_tTT_resolved->Draw("mreco>>hTT_resolved_dn",Form("weight_dn*(chan==%d)",TT_R_channum));
  TT_R_tTT_resolved->Draw("mreco>>hTT_resolved_up",Form("weight_up*(chan==%d)",TT_R_channum));
}
else{
  TT_R_tTT_resolved->SetBranchAddress("weight_vbf",&TT_R_weight);
  TT_R_tTT_resolved->SetBranchAddress("weight_vbf_up",&TT_R_weight_up);
  TT_R_tTT_resolved->SetBranchAddress("weight_vbf_dn",&TT_R_weight_dn);
  TT_R_tTT_resolved->Draw("mreco>>hTT_resolved",Form("weight_vbf*(chan==%d)",TT_R_channum));
  TT_R_tTT_resolved->Draw("mreco>>hTT_resolved_dn",Form("weight_vbf_dn*(chan==%d)",TT_R_channum));
  TT_R_tTT_resolved->Draw("mreco>>hTT_resolved_up",Form("weight_vbf_up*(chan==%d)",TT_R_channum));
}

TT_R_tTT_resolved->SetBranchAddress("chan",&TT_R_channel);

// //integral output
// if(cate_vbf!=2){
// 	cout<<"test"<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<TT_R_hTT_resolved->Integral()<<endl;
// 	cout<<"integral is " << hTT_resolved->Integral()<<endl;
// 	cout<<"lumi is" << lumi << endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hTT_resolved_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hTT_resolved_dn->Integral()<<endl;
// }
// else{
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hTT_resolved->Integral()<<endl;  /// Remember why 1.6?
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hTT_resolved_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hTT_resolved_dn->Integral()<<endl;
// }

//new tree
TTree *TT_R_cuttree =new TTree("SelectedTree","SelectedTree");
TString TT_R_mrecob = Form("%s/F",TT_R_mreco->GetName());
TT_R_cuttree->Branch(TT_R_mreco->GetName(),&TT_R_ZZMass,TT_R_mrecob);
TT_R_cuttree->Branch("weight",&TT_R_weight,"weight/F");
TT_R_cuttree->Branch("weight_up",&TT_R_weight_up,"weight_up/F");
TT_R_cuttree->Branch("weight_dn",&TT_R_weight_dn,"weight_dn/F");
//std::cout<< endl << "PROBLEM" << endl;
for(int i =0;i<TT_R_tTT_resolved->GetEntries();i++){
  TT_R_tTT_resolved->GetEntry(i);
  if(TT_R_channel==TT_R_channum){
    TT_R_cuttree->Fill();
  }
}

//Roo objects declared and filled
double TT_R_rho =1;
if(highmass==0)
  TT_R_rho=4;
RooDataSet TT_R_bkgdata ("TT_R_bkgdata"+chan+Form("_%d",cate_vbf),"",TT_R_cuttree,RooArgSet(*TT_R_mreco,*TT_R_wt),"weight");
RooKeysPdf TT_R_TT_resolvedpdf_1d("TT_R_resolved_1d","",*TT_R_mreco,TT_R_bkgdata,RooKeysPdf::MirrorBoth,TT_R_rho);

//file input
TFile *TT_R_ff = new TFile("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_TT_resolved.root");
//2D templated created and check for empty or NaN bins
TH2F *TT_R_temp_zz=(TH2F*)TT_R_ff->Get("temp_rebin_"+chan+"_resolvedSRvbf");

for(int bx=0;bx<=TT_R_temp_zz->GetNbinsX();bx++){
  for(int by=0;by<=TT_R_temp_zz->GetNbinsY();by++){

    if(!(TT_R_temp_zz->GetBinContent(bx,by) >0 || TT_R_temp_zz->GetBinContent(bx,by)<0)	){
      TT_R_temp_zz->SetBinContent(bx,by,1.0e-10);

    }
      //std::cout << endl << temp_zz->GetBinContent(bx,by) << " at " << bx << " "<< by <<endl;
  }
}
 if (TT_R_channum == 1) {
  TT_R_temp_zz->Draw("colz");
  gPad->Print("workspace/template_check.png");
 }

//other Roo objects
RooDataHist* TT_R_template_sig= new RooDataHist("TT_R_temp_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*TT_R_mreco,*TT_R_dbkg),TT_R_temp_zz);
RooHistPdf* TT_R_pdf_2d_sig = new RooHistPdf("TT_R_pdf_2d_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*TT_R_mreco,*TT_R_dbkg),*TT_R_template_sig);
RooProdPdf *TT_R_TT_resolvedpdf_2d= new RooProdPdf("TT_R_resolved_2d_"+chan+Form("_%d",cate_vbf),"",TT_R_TT_resolvedpdf_1d,Conditional(*TT_R_pdf_2d_sig,*TT_R_dbkg));
// *********************************************************************************************************************************************************


//Roo variable mreco
RooRealVar* DY_R_mreco= new RooRealVar("mreco","M_{ZZ} (GeV)",1500,low_reco,high_reco);//working with 180
if(cate_vbf==2){
DY_R_mreco->SetName("mreco_rse");
DY_R_mreco->SetTitle("mreco_rse");
}

//Roo variable dbkg
//RooRealVar* dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",0.5,0.,1.);
RooRealVar* DY_R_dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",30,0.,1.);
RooPlot* DY_R_frame= DY_R_mreco->frame(low_reco,high_reco) ;
//RooPlot* frame1= dbkg->frame(0,1) ; //???

if(highmass==2){
DY_R_mreco->setBins(nbins_reco);
DY_R_dbkg->setBins(30);
}

//tree management
TString DY_R_treename[4]={"","","_rse","_tle"};
TChain *DY_R_tDY_resolved= new TChain("SelectedTree"+DY_R_treename[cate_vbf]);

TH1F *DY_R_hDY_resolved= new TH1F ("DY_R_hDY_resolved","",nbins_reco,low_reco,high_reco);
TH1F *DY_R_hDY_resolved_up= new TH1F ("DY_R_hDY_resolved_up","",nbins_reco,low_reco,high_reco);
TH1F *DY_R_hDY_resolved_dn= new TH1F ("DY_R_hDY_resolved_dn","",nbins_reco,low_reco,high_reco);

DY_R_tDY_resolved->Add("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_DY_resolved.root");

//a code is assigned to each channel
int DY_R_channum=1;
if(chan=="all")
DY_R_channum=1;
else if(chan=="ee")
DY_R_channum=0;
else
DY_R_channum=2;

//extra Roo variables and normal variables
// RooRealVar* wt= new RooRealVar("weight","wt",1500,0.000,1000.); //125
// RooRealVar* wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,1000.);//125
// RooRealVar* wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,1000.); //125
RooRealVar* DY_R_wt= new RooRealVar("weight","wt",1500,0.000,5000.); //125
RooRealVar* DY_R_wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,5000.);//125
RooRealVar* DY_R_wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,5000.); //125

float DY_R_ZZMass;
float DY_R_weight, DY_R_weight_up, DY_R_weight_dn;
int DY_R_channel;
DY_R_tDY_resolved->SetBranchAddress("mreco",&DY_R_ZZMass);
//tDY_resolved->SetBranchAddress("dbkg",&ZZMass); //???

if(cate_vbf!=1){
DY_R_tDY_resolved->SetBranchAddress("weight",&DY_R_weight);
DY_R_tDY_resolved->SetBranchAddress("weight_up",&DY_R_weight_up);
DY_R_tDY_resolved->SetBranchAddress("weight_dn",&DY_R_weight_dn);
DY_R_tDY_resolved->Draw("mreco>>hDY_resolved",Form("weight*(chan==%d)",DY_R_channum));
DY_R_tDY_resolved->Draw("mreco>>hDY_resolved_dn",Form("weight_dn*(chan==%d)",DY_R_channum));
DY_R_tDY_resolved->Draw("mreco>>hDY_resolved_up",Form("weight_up*(chan==%d)",DY_R_channum));
}
else{
DY_R_tDY_resolved->SetBranchAddress("weight_vbf",&DY_R_weight);
DY_R_tDY_resolved->SetBranchAddress("weight_vbf_up",&DY_R_weight_up);
DY_R_tDY_resolved->SetBranchAddress("weight_vbf_dn",&DY_R_weight_dn);
DY_R_tDY_resolved->Draw("mreco>>hDY_resolved",Form("weight_vbf*(chan==%d)",DY_R_channum));
DY_R_tDY_resolved->Draw("mreco>>hDY_resolved_dn",Form("weight_vbf_dn*(chan==%d)",DY_R_channum));
DY_R_tDY_resolved->Draw("mreco>>hDY_resolved_up",Form("weight_vbf_up*(chan==%d)",DY_R_channum));
}

DY_R_tDY_resolved->SetBranchAddress("chan",&DY_R_channel);

// //integral output
// if(cate_vbf!=2){
// 	cout<<"test"<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<DY_R_hDY_resolved->Integral()<<endl;
// 	cout<<"integral is " << hDY_resolved->Integral()<<endl;
// 	cout<<"lumi is" << lumi << endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hDY_resolved_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hDY_resolved_dn->Integral()<<endl;
// }
// else{
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hDY_resolved->Integral()<<endl;  /// Remember why 1.6?
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hDY_resolved_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hDY_resolved_dn->Integral()<<endl;
// }

//new tree
TTree *DY_R_cuttree =new TTree("SelectedTree","SelectedTree");
TString DY_R_mrecob = Form("%s/F",DY_R_mreco->GetName());
DY_R_cuttree->Branch(DY_R_mreco->GetName(),&DY_R_ZZMass,DY_R_mrecob);
DY_R_cuttree->Branch("weight",&DY_R_weight,"weight/F");
DY_R_cuttree->Branch("weight_up",&DY_R_weight_up,"weight_up/F");
DY_R_cuttree->Branch("weight_dn",&DY_R_weight_dn,"weight_dn/F");
//std::cout<< endl << "PROBLEM" << endl;
for(int i =0;i<DY_R_tDY_resolved->GetEntries();i++){
DY_R_tDY_resolved->GetEntry(i);
if(DY_R_channel==DY_R_channum){
DY_R_cuttree->Fill();
}
}

//Roo objects declared and filled
double DY_R_rho =1;
if(highmass==0)
DY_R_rho=4;
RooDataSet DY_R_bkgdata ("DY_R_bkgdata"+chan+Form("_%d",cate_vbf),"",DY_R_cuttree,RooArgSet(*DY_R_mreco,*DY_R_wt),"weight");
RooKeysPdf DY_R_DY_resolvedpdf_1d("DY_R_resolved_1d","",*DY_R_mreco,DY_R_bkgdata,RooKeysPdf::MirrorBoth,DY_R_rho);

//file input
TFile *DY_R_ff = new TFile("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_DY_resolved.root");
//2D templated created and check for empty or NaN bins
TH2F *DY_R_temp_zz=(TH2F*)DY_R_ff->Get("temp_rebin_"+chan+"_resolvedSRvbf");

for(int bx=0;bx<=DY_R_temp_zz->GetNbinsX();bx++){
for(int by=0;by<=DY_R_temp_zz->GetNbinsY();by++){

if(!(DY_R_temp_zz->GetBinContent(bx,by) >0 || DY_R_temp_zz->GetBinContent(bx,by)<0)	){
  DY_R_temp_zz->SetBinContent(bx,by,1.0e-10);

}
  //std::cout << endl << temp_zz->GetBinContent(bx,by) << " at " << bx << " "<< by <<endl;
}
}
if (DY_R_channum == 1) {
DY_R_temp_zz->Draw("colz");
gPad->Print("workspace/template_check.png");
}

//other Roo objects
RooDataHist* DY_R_template_sig= new RooDataHist("DY_R_temp_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*DY_R_mreco,*DY_R_dbkg),DY_R_temp_zz);
RooHistPdf* DY_R_pdf_2d_sig = new RooHistPdf("DY_R_pdf_2d_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*DY_R_mreco,*DY_R_dbkg),*DY_R_template_sig);
RooProdPdf *DY_R_DY_resolvedpdf_2d= new RooProdPdf("DY_R_resolved_2d_"+chan+Form("_%d",cate_vbf),"",DY_R_DY_resolvedpdf_1d,Conditional(*DY_R_pdf_2d_sig,*DY_R_dbkg));
// *********************************************************************************************************************************************************


// // //Data
// // TChain *tdata = new TChain("SelectedTree"+TT_R_treename[cate_vbf]);
// // tdata->Add("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/data_2l.root");
// // TTree* reducetree= tdata->CopyTree(Form("chan==%d&&vbfcate==%d",TT_R_channum,cate_vbf));
// // RooDataSet* data_obs_1d= new RooDataSet("data_obs_1d","data_obs_1d",reducetree,*TT_R_mreco);
// // RooDataSet* data_obs_2d = new RooDataSet("data_obs_2d","data_obs_2d",reducetree,RooArgSet(*TT_R_mreco,*TT_R_dbkg));
//
// // *********************************************************************************************************************************************************
//
//Roo variable mreco
RooRealVar* WW_ZZ_R_mreco= new RooRealVar("mreco","M_{ZZ} (GeV)",1500,low_reco,high_reco);//working with 180
if(cate_vbf==2){
WW_ZZ_R_mreco->SetName("mreco_rse");
WW_ZZ_R_mreco->SetTitle("mreco_rse");
}

//Roo variable dbkg
//RooRealVar* dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",0.5,0.,1.);
RooRealVar* WW_ZZ_R_dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",30,0.,1.);
RooPlot* WW_ZZ_R_frame= WW_ZZ_R_mreco->frame(low_reco,high_reco) ;
//RooPlot* frame1= dbkg->frame(0,1) ; //???

if(highmass==2){
WW_ZZ_R_mreco->setBins(nbins_reco);
WW_ZZ_R_dbkg->setBins(30);
}

//tree management
TString WW_ZZ_R_treename[4]={"","","_rse","_tle"};
TChain *WW_ZZ_R_tWW_ZZ_resolved= new TChain("SelectedTree"+WW_ZZ_R_treename[cate_vbf]);

TH1F *WW_ZZ_R_hWW_ZZ_resolved= new TH1F ("WW_ZZ_R_hWW_ZZ_resolved","",nbins_reco,low_reco,high_reco);
TH1F *WW_ZZ_R_hWW_ZZ_resolved_up= new TH1F ("WW_ZZ_R_hWW_ZZ_resolved_up","",nbins_reco,low_reco,high_reco);
TH1F *WW_ZZ_R_hWW_ZZ_resolved_dn= new TH1F ("WW_ZZ_R_hWW_ZZ_resolved_dn","",nbins_reco,low_reco,high_reco);

WW_ZZ_R_tWW_ZZ_resolved->Add("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_WW_ZZ_resolved.root");

//a code is assigned to each channel
int WW_ZZ_R_channum=1;
if(chan=="all")
WW_ZZ_R_channum=1;
else if(chan=="ee")
WW_ZZ_R_channum=0;
else
WW_ZZ_R_channum=2;

//extra Roo variables and normal variables
// RooRealVar* wt= new RooRealVar("weight","wt",1500,0.000,1000.); //125
// RooRealVar* wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,1000.);//125
// RooRealVar* wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,1000.); //125
RooRealVar* WW_ZZ_R_wt= new RooRealVar("weight","wt",1500,0.000,5000.); //125
RooRealVar* WW_ZZ_R_wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,5000.);//125
RooRealVar* WW_ZZ_R_wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,5000.); //125

float WW_ZZ_R_ZZMass;
float WW_ZZ_R_weight, WW_ZZ_R_weight_up, WW_ZZ_R_weight_dn;
int WW_ZZ_R_channel;
WW_ZZ_R_tWW_ZZ_resolved->SetBranchAddress("mreco",&WW_ZZ_R_ZZMass);
//tWW_ZZ_resolved->SetBranchAddress("dbkg",&ZZMass); //???

if(cate_vbf!=1){
WW_ZZ_R_tWW_ZZ_resolved->SetBranchAddress("weight",&WW_ZZ_R_weight);
WW_ZZ_R_tWW_ZZ_resolved->SetBranchAddress("weight_up",&WW_ZZ_R_weight_up);
WW_ZZ_R_tWW_ZZ_resolved->SetBranchAddress("weight_dn",&WW_ZZ_R_weight_dn);
WW_ZZ_R_tWW_ZZ_resolved->Draw("mreco>>hWW_ZZ_resolved",Form("weight*(chan==%d)",WW_ZZ_R_channum));
WW_ZZ_R_tWW_ZZ_resolved->Draw("mreco>>hWW_ZZ_resolved_dn",Form("weight_dn*(chan==%d)",WW_ZZ_R_channum));
WW_ZZ_R_tWW_ZZ_resolved->Draw("mreco>>hWW_ZZ_resolved_up",Form("weight_up*(chan==%d)",WW_ZZ_R_channum));
}
else{
WW_ZZ_R_tWW_ZZ_resolved->SetBranchAddress("weight_vbf",&WW_ZZ_R_weight);
WW_ZZ_R_tWW_ZZ_resolved->SetBranchAddress("weight_vbf_up",&WW_ZZ_R_weight_up);
WW_ZZ_R_tWW_ZZ_resolved->SetBranchAddress("weight_vbf_dn",&WW_ZZ_R_weight_dn);
WW_ZZ_R_tWW_ZZ_resolved->Draw("mreco>>hWW_ZZ_resolved",Form("weight_vbf*(chan==%d)",WW_ZZ_R_channum));
WW_ZZ_R_tWW_ZZ_resolved->Draw("mreco>>hWW_ZZ_resolved_dn",Form("weight_vbf_dn*(chan==%d)",WW_ZZ_R_channum));
WW_ZZ_R_tWW_ZZ_resolved->Draw("mreco>>hWW_ZZ_resolved_up",Form("weight_vbf_up*(chan==%d)",WW_ZZ_R_channum));
}

WW_ZZ_R_tWW_ZZ_resolved->SetBranchAddress("chan",&WW_ZZ_R_channel);

// //integral output
// if(cate_vbf!=2){
// 	cout<<"test"<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<WW_ZZ_R_hWW_ZZ_resolved->Integral()<<endl;
// 	cout<<"integral is " << hWW_ZZ_resolved->Integral()<<endl;
// 	cout<<"lumi is" << lumi << endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hWW_ZZ_resolved_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hWW_ZZ_resolved_dn->Integral()<<endl;
// }
// else{
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hWW_ZZ_resolved->Integral()<<endl;  /// Remember why 1.6?
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hWW_ZZ_resolved_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hWW_ZZ_resolved_dn->Integral()<<endl;
// }

//new tree
TTree *WW_ZZ_R_cuttree =new TTree("SelectedTree","SelectedTree");
TString WW_ZZ_R_mrecob = Form("%s/F",WW_ZZ_R_mreco->GetName());
WW_ZZ_R_cuttree->Branch(WW_ZZ_R_mreco->GetName(),&WW_ZZ_R_ZZMass,WW_ZZ_R_mrecob);
WW_ZZ_R_cuttree->Branch("weight",&WW_ZZ_R_weight,"weight/F");
WW_ZZ_R_cuttree->Branch("weight_up",&WW_ZZ_R_weight_up,"weight_up/F");
WW_ZZ_R_cuttree->Branch("weight_dn",&WW_ZZ_R_weight_dn,"weight_dn/F");
//std::cout<< endl << "PROBLEM" << endl;
for(int i =0;i<WW_ZZ_R_tWW_ZZ_resolved->GetEntries();i++){
WW_ZZ_R_tWW_ZZ_resolved->GetEntry(i);
if(WW_ZZ_R_channel==WW_ZZ_R_channum){
WW_ZZ_R_cuttree->Fill();
}
}

//Roo objects declared and filled
double WW_ZZ_R_rho =1;
if(highmass==0)
WW_ZZ_R_rho=4;
RooDataSet WW_ZZ_R_bkgdata ("WW_ZZ_R_bkgdata"+chan+Form("_%d",cate_vbf),"",WW_ZZ_R_cuttree,RooArgSet(*WW_ZZ_R_mreco,*WW_ZZ_R_wt),"weight");
RooKeysPdf WW_ZZ_R_WW_ZZ_resolvedpdf_1d("WW_ZZ_R_resolved_1d","",*WW_ZZ_R_mreco,WW_ZZ_R_bkgdata,RooKeysPdf::MirrorBoth,WW_ZZ_R_rho);

//file input
TFile *WW_ZZ_R_ff = new TFile("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_WW_ZZ_resolved.root");
//2D templated created and check for empty or NaN bins
TH2F *WW_ZZ_R_temp_zz=(TH2F*)WW_ZZ_R_ff->Get("temp_rebin_"+chan+"_resolvedSRvbf");

for(int bx=0;bx<=WW_ZZ_R_temp_zz->GetNbinsX();bx++){
for(int by=0;by<=WW_ZZ_R_temp_zz->GetNbinsY();by++){

if(!(WW_ZZ_R_temp_zz->GetBinContent(bx,by) >0 || WW_ZZ_R_temp_zz->GetBinContent(bx,by)<0)	){
  WW_ZZ_R_temp_zz->SetBinContent(bx,by,1.0e-10);

}
  //std::cout << endl << temp_zz->GetBinContent(bx,by) << " at " << bx << " "<< by <<endl;
}
}
if (WW_ZZ_R_channum == 1) {
WW_ZZ_R_temp_zz->Draw("colz");
gPad->Print("workspace/template_check.png");
}

//other Roo objects
RooDataHist* WW_ZZ_R_template_sig= new RooDataHist("WW_ZZ_R_temp_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*WW_ZZ_R_mreco,*WW_ZZ_R_dbkg),WW_ZZ_R_temp_zz);
RooHistPdf* WW_ZZ_R_pdf_2d_sig = new RooHistPdf("WW_ZZ_R_pdf_2d_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*WW_ZZ_R_mreco,*WW_ZZ_R_dbkg),*WW_ZZ_R_template_sig);
RooProdPdf *WW_ZZ_R_WW_ZZ_resolvedpdf_2d= new RooProdPdf("WW_ZZ_R_resolved_2d_"+chan+Form("_%d",cate_vbf),"",WW_ZZ_R_WW_ZZ_resolvedpdf_1d,Conditional(*WW_ZZ_R_pdf_2d_sig,*WW_ZZ_R_dbkg));

// //**************************************************************
//
// //Roo variable mreco
RooRealVar* sig_ZZ_R_mreco= new RooRealVar("mreco","M_{ZZ} (GeV)",1500,low_reco,high_reco);//working with 180
if(cate_vbf==2){
sig_ZZ_R_mreco->SetName("mreco_rse");
sig_ZZ_R_mreco->SetTitle("mreco_rse");
}

//Roo variable dbkg
//RooRealVar* dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",0.5,0.,1.);
RooRealVar* sig_ZZ_R_dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",30,0.,1.);
RooPlot* sig_ZZ_R_frame= sig_ZZ_R_mreco->frame(low_reco,high_reco) ;
//RooPlot* frame1= dbkg->frame(0,1) ; //???

if(highmass==2){
sig_ZZ_R_mreco->setBins(nbins_reco);
sig_ZZ_R_dbkg->setBins(30);
}

//tree management
TString sig_ZZ_R_treename[4]={"","","_rse","_tle"};
TChain *sig_ZZ_R_tsig_ZZ_resolved= new TChain("SelectedTree"+sig_ZZ_R_treename[cate_vbf]);

TH1F *sig_ZZ_R_hsig_ZZ_resolved= new TH1F ("sig_ZZ_R_hsig_ZZ_resolved","",nbins_reco,low_reco,high_reco);
TH1F *sig_ZZ_R_hsig_ZZ_resolved_up= new TH1F ("sig_ZZ_R_hsig_ZZ_resolved_up","",nbins_reco,low_reco,high_reco);
TH1F *sig_ZZ_R_hsig_ZZ_resolved_dn= new TH1F ("sig_ZZ_R_hsig_ZZ_resolved_dn","",nbins_reco,low_reco,high_reco);

sig_ZZ_R_tsig_ZZ_resolved->Add("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_sig_ZZ_resolved.root");

//a code is assigned to each channel
int sig_ZZ_R_channum=1;
if(chan=="all")
sig_ZZ_R_channum=1;
else if(chan=="ee")
sig_ZZ_R_channum=0;
else
sig_ZZ_R_channum=2;

//extra Roo variables and normal variables
// RooRealVar* wt= new RooRealVar("weight","wt",1500,0.000,1000.); //125
// RooRealVar* wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,1000.);//125
// RooRealVar* wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,1000.); //125
RooRealVar* sig_ZZ_R_wt= new RooRealVar("weight","wt",1500,0.000,5000.); //125
RooRealVar* sig_ZZ_R_wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,5000.);//125
RooRealVar* sig_ZZ_R_wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,5000.); //125

float sig_ZZ_R_ZZMass;
float sig_ZZ_R_weight, sig_ZZ_R_weight_up, sig_ZZ_R_weight_dn;
int sig_ZZ_R_channel;
sig_ZZ_R_tsig_ZZ_resolved->SetBranchAddress("mreco",&sig_ZZ_R_ZZMass);
//tsig_ZZ_resolved->SetBranchAddress("dbkg",&ZZMass); //???

if(cate_vbf!=1){
sig_ZZ_R_tsig_ZZ_resolved->SetBranchAddress("weight",&sig_ZZ_R_weight);
sig_ZZ_R_tsig_ZZ_resolved->SetBranchAddress("weight_up",&sig_ZZ_R_weight_up);
sig_ZZ_R_tsig_ZZ_resolved->SetBranchAddress("weight_dn",&sig_ZZ_R_weight_dn);
sig_ZZ_R_tsig_ZZ_resolved->Draw("mreco>>hsig_ZZ_resolved",Form("weight*(chan==%d)",sig_ZZ_R_channum));
sig_ZZ_R_tsig_ZZ_resolved->Draw("mreco>>hsig_ZZ_resolved_dn",Form("weight_dn*(chan==%d)",sig_ZZ_R_channum));
sig_ZZ_R_tsig_ZZ_resolved->Draw("mreco>>hsig_ZZ_resolved_up",Form("weight_up*(chan==%d)",sig_ZZ_R_channum));
}
else{
sig_ZZ_R_tsig_ZZ_resolved->SetBranchAddress("weight_vbf",&sig_ZZ_R_weight);
sig_ZZ_R_tsig_ZZ_resolved->SetBranchAddress("weight_vbf_up",&sig_ZZ_R_weight_up);
sig_ZZ_R_tsig_ZZ_resolved->SetBranchAddress("weight_vbf_dn",&sig_ZZ_R_weight_dn);
sig_ZZ_R_tsig_ZZ_resolved->Draw("mreco>>hsig_ZZ_resolved",Form("weight_vbf*(chan==%d)",sig_ZZ_R_channum));
sig_ZZ_R_tsig_ZZ_resolved->Draw("mreco>>hsig_ZZ_resolved_dn",Form("weight_vbf_dn*(chan==%d)",sig_ZZ_R_channum));
sig_ZZ_R_tsig_ZZ_resolved->Draw("mreco>>hsig_ZZ_resolved_up",Form("weight_vbf_up*(chan==%d)",sig_ZZ_R_channum));
}

sig_ZZ_R_tsig_ZZ_resolved->SetBranchAddress("chan",&sig_ZZ_R_channel);

// //integral output
// if(cate_vbf!=2){
// 	cout<<"test"<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<sig_ZZ_R_hsig_ZZ_resolved->Integral()<<endl;
// 	cout<<"integral is " << hsig_ZZ_resolved->Integral()<<endl;
// 	cout<<"lumi is" << lumi << endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsig_ZZ_resolved_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsig_ZZ_resolved_dn->Integral()<<endl;
// }
// else{
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsig_ZZ_resolved->Integral()<<endl;  /// Remember why 1.6?
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsig_ZZ_resolved_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsig_ZZ_resolved_dn->Integral()<<endl;
// }

//new tree
TTree *sig_ZZ_R_cuttree =new TTree("SelectedTree","SelectedTree");
TString sig_ZZ_R_mrecob = Form("%s/F",sig_ZZ_R_mreco->GetName());
sig_ZZ_R_cuttree->Branch(sig_ZZ_R_mreco->GetName(),&sig_ZZ_R_ZZMass,sig_ZZ_R_mrecob);
sig_ZZ_R_cuttree->Branch("weight",&sig_ZZ_R_weight,"weight/F");
sig_ZZ_R_cuttree->Branch("weight_up",&sig_ZZ_R_weight_up,"weight_up/F");
sig_ZZ_R_cuttree->Branch("weight_dn",&sig_ZZ_R_weight_dn,"weight_dn/F");
//std::cout<< endl << "PROBLEM" << endl;
for(int i =0;i<sig_ZZ_R_tsig_ZZ_resolved->GetEntries();i++){
sig_ZZ_R_tsig_ZZ_resolved->GetEntry(i);
if(sig_ZZ_R_channel==sig_ZZ_R_channum){
sig_ZZ_R_cuttree->Fill();
}
}

//Roo objects declared and filled
double sig_ZZ_R_rho =1;
if(highmass==0)
sig_ZZ_R_rho=4;
RooDataSet sig_ZZ_R_bkgdata ("sig_ZZ_R_bkgdata"+chan+Form("_%d",cate_vbf),"",sig_ZZ_R_cuttree,RooArgSet(*sig_ZZ_R_mreco,*sig_ZZ_R_wt),"weight");
RooKeysPdf sig_ZZ_R_sig_ZZ_resolvedpdf_1d("sig_ZZ_R_resolved_1d","",*sig_ZZ_R_mreco,sig_ZZ_R_bkgdata,RooKeysPdf::MirrorBoth,sig_ZZ_R_rho);

//file input
TFile *sig_ZZ_R_ff = new TFile("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_sig_ZZ_resolved.root");
//2D templated created and check for empty or NaN bins
TH2F *sig_ZZ_R_temp_zz=(TH2F*)sig_ZZ_R_ff->Get("temp_rebin_"+chan+"_resolvedSRvbf");

for(int bx=0;bx<=sig_ZZ_R_temp_zz->GetNbinsX();bx++){
for(int by=0;by<=sig_ZZ_R_temp_zz->GetNbinsY();by++){

if(!(sig_ZZ_R_temp_zz->GetBinContent(bx,by) >0 || sig_ZZ_R_temp_zz->GetBinContent(bx,by)<0)	){
  sig_ZZ_R_temp_zz->SetBinContent(bx,by,1.0e-10);

}
  //std::cout << endl << temp_zz->GetBinContent(bx,by) << " at " << bx << " "<< by <<endl;
}
}
if (sig_ZZ_R_channum == 1) {
	TCanvas *c = new TCanvas();
	c->SetRightMargin(0.2);
	c->SetLeftMargin(0.2);
	c->SetBottomMargin(0.2);
gStyle->SetOptStat(0);
sig_ZZ_R_temp_zz->SetZTitle("weight");
sig_ZZ_R_temp_zz->SetXTitle("m_{ZZ}");
sig_ZZ_R_temp_zz->SetYTitle("D_{k}");
sig_ZZ_R_temp_zz->Draw("colz");
gPad->Print("workspace/sig_zz_resolved.png");
}

//other Roo objects
RooDataHist* sig_ZZ_R_template_sig= new RooDataHist("sig_ZZ_R_temp_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*sig_ZZ_R_mreco,*sig_ZZ_R_dbkg),sig_ZZ_R_temp_zz);
RooHistPdf* sig_ZZ_R_pdf_2d_sig = new RooHistPdf("sig_ZZ_R_pdf_2d_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*sig_ZZ_R_mreco,*sig_ZZ_R_dbkg),*sig_ZZ_R_template_sig);
RooProdPdf *sig_ZZ_R_sig_ZZ_resolvedpdf_2d= new RooProdPdf("sig_ZZ_R_resolved_2d_"+chan+Form("_%d",cate_vbf),"",sig_ZZ_R_sig_ZZ_resolvedpdf_1d,Conditional(*sig_ZZ_R_pdf_2d_sig,*sig_ZZ_R_dbkg));

//***************************************************************


//Roo variable mreco
RooRealVar* segnaleWZ_R_mreco= new RooRealVar("mreco","M_{ZZ} (GeV)",1500,low_reco,high_reco);//working with 180
if(cate_vbf==2){
segnaleWZ_R_mreco->SetName("mreco_rse");
segnaleWZ_R_mreco->SetTitle("mreco_rse");
}

//Roo variable dbkg
//RooRealVar* dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",0.5,0.,1.);
RooRealVar* segnaleWZ_R_dbkg= new RooRealVar("dbkg","Dbkg_{kin} ",30,0.,1.);
RooPlot* segnaleWZ_R_frame= segnaleWZ_R_mreco->frame(low_reco,high_reco) ;
//RooPlot* frame1= dbkg->frame(0,1) ; //???

if(highmass==2){
segnaleWZ_R_mreco->setBins(nbins_reco);
segnaleWZ_R_dbkg->setBins(30);
}

//tree management
TString segnaleWZ_R_treename[4]={"","","_rse","_tle"};
TChain *segnaleWZ_R_tsegnaleWZ_resolved= new TChain("SelectedTree"+segnaleWZ_R_treename[cate_vbf]);

TH1F *segnaleWZ_R_hsegnaleWZ_resolved= new TH1F ("segnaleWZ_R_hsegnaleWZ_resolved","",nbins_reco,low_reco,high_reco);
TH1F *segnaleWZ_R_hsegnaleWZ_resolved_up= new TH1F ("segnaleWZ_R_hsegnaleWZ_resolved_up","",nbins_reco,low_reco,high_reco);
TH1F *segnaleWZ_R_hsegnaleWZ_resolved_dn= new TH1F ("segnaleWZ_R_hsegnaleWZ_resolved_dn","",nbins_reco,low_reco,high_reco);

segnaleWZ_R_tsegnaleWZ_resolved->Add("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_sig_WZ_resolved.root");

//a code is assigned to each channel
int segnaleWZ_R_channum=1;
if(chan=="all")
segnaleWZ_R_channum=1;
else if(chan=="ee")
segnaleWZ_R_channum=0;
else
segnaleWZ_R_channum=2;

//extra Roo variables and normal variables
// RooRealVar* wt= new RooRealVar("weight","wt",1500,0.000,1000.); //125
// RooRealVar* wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,1000.);//125
// RooRealVar* wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,1000.); //125
RooRealVar* segnaleWZ_R_wt= new RooRealVar("weight","wt",1500,0.000,5000.); //125
RooRealVar* segnaleWZ_R_wt_up= new RooRealVar("weight_up","wt_up",1500,0.000,5000.);//125
RooRealVar* segnaleWZ_R_wt_dn= new RooRealVar("weight_dn","wt_dn",1500,0.000,5000.); //125

float segnaleWZ_R_ZZMass;
float segnaleWZ_R_weight, segnaleWZ_R_weight_up, segnaleWZ_R_weight_dn;
int segnaleWZ_R_channel;
segnaleWZ_R_tsegnaleWZ_resolved->SetBranchAddress("mreco",&segnaleWZ_R_ZZMass);
//tsegnaleWZ_resolved->SetBranchAddress("dbkg",&ZZMass); //???

if(cate_vbf!=1){
segnaleWZ_R_tsegnaleWZ_resolved->SetBranchAddress("weight",&segnaleWZ_R_weight);
segnaleWZ_R_tsegnaleWZ_resolved->SetBranchAddress("weight_up",&segnaleWZ_R_weight_up);
segnaleWZ_R_tsegnaleWZ_resolved->SetBranchAddress("weight_dn",&segnaleWZ_R_weight_dn);
segnaleWZ_R_tsegnaleWZ_resolved->Draw("mreco>>hsegnaleWZ_resolved",Form("weight*(chan==%d)",segnaleWZ_R_channum));
segnaleWZ_R_tsegnaleWZ_resolved->Draw("mreco>>hsegnaleWZ_resolved_dn",Form("weight_dn*(chan==%d)",segnaleWZ_R_channum));
segnaleWZ_R_tsegnaleWZ_resolved->Draw("mreco>>hsegnaleWZ_resolved_up",Form("weight_up*(chan==%d)",segnaleWZ_R_channum));
}
else{
segnaleWZ_R_tsegnaleWZ_resolved->SetBranchAddress("weight_vbf",&segnaleWZ_R_weight);
segnaleWZ_R_tsegnaleWZ_resolved->SetBranchAddress("weight_vbf_up",&segnaleWZ_R_weight_up);
segnaleWZ_R_tsegnaleWZ_resolved->SetBranchAddress("weight_vbf_dn",&segnaleWZ_R_weight_dn);
segnaleWZ_R_tsegnaleWZ_resolved->Draw("mreco>>hsegnaleWZ_resolved",Form("weight_vbf*(chan==%d)",segnaleWZ_R_channum));
segnaleWZ_R_tsegnaleWZ_resolved->Draw("mreco>>hsegnaleWZ_resolved_dn",Form("weight_vbf_dn*(chan==%d)",segnaleWZ_R_channum));
segnaleWZ_R_tsegnaleWZ_resolved->Draw("mreco>>hsegnaleWZ_resolved_up",Form("weight_vbf_up*(chan==%d)",segnaleWZ_R_channum));
}

segnaleWZ_R_tsegnaleWZ_resolved->SetBranchAddress("chan",&segnaleWZ_R_channel);

// //integral output
// if(cate_vbf!=2){
// 	cout<<"test"<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<segnaleWZ_R_hsegnaleWZ_resolved->Integral()<<endl;
// 	cout<<"integral is " << hsegnaleWZ_resolved->Integral()<<endl;
// 	cout<<"lumi is" << lumi << endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsegnaleWZ_resolved_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsegnaleWZ_resolved_dn->Integral()<<endl;
// }
// else{
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsegnaleWZ_resolved->Integral()<<endl;  /// Remember why 1.6?
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsegnaleWZ_resolved_up->Integral()<<endl;
// 	cout<<chan<<"\t"<<cate_vbf<<"\t"<<hsegnaleWZ_resolved_dn->Integral()<<endl;
// }

//new tree
TTree *segnaleWZ_R_cuttree =new TTree("SelectedTree","SelectedTree");
TString segnaleWZ_R_mrecob = Form("%s/F",segnaleWZ_R_mreco->GetName());
segnaleWZ_R_cuttree->Branch(segnaleWZ_R_mreco->GetName(),&segnaleWZ_R_ZZMass,segnaleWZ_R_mrecob);
segnaleWZ_R_cuttree->Branch("weight",&segnaleWZ_R_weight,"weight/F");
segnaleWZ_R_cuttree->Branch("weight_up",&segnaleWZ_R_weight_up,"weight_up/F");
segnaleWZ_R_cuttree->Branch("weight_dn",&segnaleWZ_R_weight_dn,"weight_dn/F");
//std::cout<< endl << "PROBLEM" << endl;
for(int i =0;i<segnaleWZ_R_tsegnaleWZ_resolved->GetEntries();i++){
segnaleWZ_R_tsegnaleWZ_resolved->GetEntry(i);
if(segnaleWZ_R_channel==segnaleWZ_R_channum){
segnaleWZ_R_cuttree->Fill();
}
}

//Roo objects declared and filled
double segnaleWZ_R_rho =1;
if(highmass==0)
segnaleWZ_R_rho=4;
RooDataSet segnaleWZ_R_bkgdata ("segnaleWZ_R_bkgdata"+chan+Form("_%d",cate_vbf),"",segnaleWZ_R_cuttree,RooArgSet(*segnaleWZ_R_mreco,*segnaleWZ_R_wt),"weight");
RooKeysPdf segnaleWZ_R_segnaleWZ_resolvedpdf_1d("segnaleWZ_R_resolved_1d","",*segnaleWZ_R_mreco,segnaleWZ_R_bkgdata,RooKeysPdf::MirrorBoth,segnaleWZ_R_rho);

//file input
TFile *segnaleWZ_R_ff = new TFile("/afs/cern.ch/work/c/cthorbur/VBF_ANALYSIS/2l2j/CMSSW_8_0_26_patch1/src/ZZAnalysis/Plotter/RooFitInput/MC_sig_WZ_resolved.root");
//2D templated created and check for empty or NaN bins
TH2F *segnaleWZ_R_temp_zz=(TH2F*)segnaleWZ_R_ff->Get("temp_rebin_"+chan+"_resolvedSRvbf");

for(int bx=0;bx<=segnaleWZ_R_temp_zz->GetNbinsX();bx++){
for(int by=0;by<=segnaleWZ_R_temp_zz->GetNbinsY();by++){

if(!(segnaleWZ_R_temp_zz->GetBinContent(bx,by) >0 || segnaleWZ_R_temp_zz->GetBinContent(bx,by)<0)	){
  segnaleWZ_R_temp_zz->SetBinContent(bx,by,1.0e-10);

}
  std::cout << endl << segnaleWZ_R_temp_zz->GetBinContent(bx,by) << " at " << bx << " "<< by <<endl;
}
}
if (segnaleWZ_R_channum == 1) {
segnaleWZ_R_temp_zz->Draw("colz");
gPad->Print("workspace/template_check.png");
}

//other Roo objects
RooDataHist* segnaleWZ_R_template_sig= new RooDataHist("segnaleWZ_R_temp_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*segnaleWZ_R_mreco,*segnaleWZ_R_dbkg),segnaleWZ_R_temp_zz);
RooHistPdf* segnaleWZ_R_pdf_2d_sig = new RooHistPdf("segnaleWZ_R_pdf_2d_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*segnaleWZ_R_mreco,*segnaleWZ_R_dbkg),*segnaleWZ_R_template_sig);
RooProdPdf *segnaleWZ_R_segnaleWZ_resolvedpdf_2d= new RooProdPdf("segnaleWZ_R_resolved_2d_"+chan+Form("_%d",cate_vbf),"",segnaleWZ_R_segnaleWZ_resolvedpdf_1d,Conditional(*segnaleWZ_R_pdf_2d_sig,*segnaleWZ_R_dbkg));


//*************************************************************+

//Roo real var mreco

//Import to workspace and write output file
		if(is2D){
			data_obs_2d->SetNameTitle("data_obs","data_obs");
		// 	 WW_ZZ_M_WW_ZZ_mergedpdf_2d->SetNameTitle("bkg_vz","bkg_vz");
		// 	 TT_M_TT_mergedpdf_2d->SetNameTitle("bkg_ttbar","bkg_ttbar");
		// 	 DY_M_DY_mergedpdf_2d->SetNameTitle("bkg_zjets","bkg_zjets");
		// sig_WZ_M_sig_WZ_mergedpdf_2d->SetNameTitle("sig_wz","sig_wz");
		// 	sig_ZZ_M_sig_ZZ_mergedpdf_2d->SetNameTitle("sig_zz","sig_zz");
			WW_ZZ_R_WW_ZZ_resolvedpdf_2d->SetNameTitle("bkg_vz","bkg_vz");
		   	TT_R_TT_resolvedpdf_2d->SetNameTitle("bkg_ttbar","bkg_ttbar");
		   	DY_R_DY_resolvedpdf_2d->SetNameTitle("bkg_zjets","bkg_zjets");
		   	segnaleWZ_R_segnaleWZ_resolvedpdf_2d->SetNameTitle("sig_wz","sig_wz");
		   	sig_ZZ_R_sig_ZZ_resolvedpdf_2d->SetNameTitle("sig_zz","sig_zz");

			// w.import(*WW_ZZ_M_WW_ZZ_mergedpdf_2d);
			// w.import(*TT_M_TT_mergedpdf_2d);
			// w.import(*DY_M_DY_mergedpdf_2d);
			// w.import(*sig_WZ_M_sig_WZ_mergedpdf_2d);
			// w.import(*sig_ZZ_M_sig_ZZ_mergedpdf_2d);
			w.import(*WW_ZZ_R_WW_ZZ_resolvedpdf_2d);
		  w.import(*TT_R_TT_resolvedpdf_2d);
		  w.import(*DY_R_DY_resolvedpdf_2d);
		  w.import(*segnaleWZ_R_segnaleWZ_resolvedpdf_2d);
		   w.import(*sig_ZZ_R_sig_ZZ_resolvedpdf_2d);

			w.import(*data_obs_2d);
		}
		else{
			data_obs_1d->SetNameTitle("data_obs","data_obs");
			// WW_ZZ_M_WW_ZZ_mergedpdf_1d.SetNameTitle("bkg_vz","bkg_vz");
			// TT_M_TT_mergedpdf_1d.SetNameTitle("bkg_ttbar","bkg_ttbar");
			// DY_M_DY_mergedpdf_1d.SetNameTitle("bkg_zjets","bkg_zjets");
			// sig_WZ_M_sig_WZ_mergedpdf_1d.SetNameTitle("sig_wz","sig_wz");
			// sig_ZZ_M_sig_ZZ_mergedpdf_1d.SetNameTitle("sig_zz","sig_zz");
			WW_ZZ_R_WW_ZZ_resolvedpdf_1d.SetNameTitle("bkg_vz","bkg_vz");
			TT_R_TT_resolvedpdf_1d.SetNameTitle("bkg_ttbar","bkg_ttbar");
			DY_R_DY_resolvedpdf_1d.SetNameTitle("bkg_zjets","bkg_zjets");
			 segnaleWZ_R_segnaleWZ_resolvedpdf_1d.SetNameTitle("sig_wz","sig_wz");
			 sig_ZZ_R_sig_ZZ_resolvedpdf_1d.SetNameTitle("sig_zz","sig_zz");


			// w.import(WW_ZZ_M_WW_ZZ_mergedpdf_1d);
			// w.import(TT_M_TT_mergedpdf_1d);
			// w.import(DY_M_DY_mergedpdf_1d);
			// w.import(sig_WZ_M_sig_WZ_mergedpdf_1d);
			// w.import(sig_ZZ_M_sig_ZZ_mergedpdf_1d);
			w.import(WW_ZZ_R_WW_ZZ_resolvedpdf_1d);
		  w.import(TT_R_TT_resolvedpdf_1d);
		  w.import(DY_R_DY_resolvedpdf_1d);
		   w.import(segnaleWZ_R_segnaleWZ_resolvedpdf_1d);
		  w.import(sig_ZZ_R_sig_ZZ_resolvedpdf_1d);

			w.import(*data_obs_1d);
		}



TFile *fwork ;
if(highmass==0)
	fwork= new TFile("workspace/"+chan+Form("%dS_13TeV.input_func.root",cate_vbf),"recreate");
if(highmass==1)
	fwork= new TFile("workspace/"+chan+Form("%dS_13TeV.input_func.root",cate_vbf),"recreate");
	if(highmass==2){
		if(is2D)
			fwork= new TFile("workspace/"+chan+Form("%dS_13TeV.input_func.root",cate_vbf),"recreate");
		else
			fwork= new TFile("workspace/"+chan+Form("%dS_13TeV.input_func.root",cate_vbf),"recreate");
	}
fwork->cd();
w.Write();
fwork->Close();
}

void bkgWorkspace_2l2q(TString chan, int vbfcate, int highmass, int is2D){
	dosomething(chan,vbfcate,highmass,is2D);
}
