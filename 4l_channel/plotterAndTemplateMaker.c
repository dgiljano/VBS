// #include "external_cConstants.h"
#include <TSpline.h>
#include <TString.h>
#include <memory>
#include "TMath.h"
#include "helper_functions.h"

string jet_pt_cut = "jet_pt_gt_30";

TH2F* rebinTemplate(TH2F* orig, int year=2016, int itype=0) 
{
  	char filename[300];    
	char pname[30];
  	if (itype == 0) sprintf(pname,"qqzz");
  	if (itype == 1) sprintf(pname,"ggzz");
  	if (itype == 2) sprintf(pname,"vbs");

  	TH2F* tempt = (TH2F*)orig->Clone();
   	for(int binx=0;binx<tempt->GetXaxis()->GetNbins();binx++)
   	{
		double inttmp1 = tempt->Integral(binx+1,binx+1);
	   	for(int biny=0;biny<tempt->GetNbinsY();biny++)
	   	{
			if(inttmp1!=0 )
	  			tempt->SetBinContent(binx+1, biny+1, tempt->GetBinContent(binx+1,biny+1)/inttmp1);
    	}
   	}      
   	
	tempt->Smooth();
   	int nbinsfinal=230;
   	double xbinfinal[231] = {160,162,164,166,168,170,172,174,176,178,180,182,184,186,188,190,192,194,196,198,200,202,204,206,208,210,212,214,216,218,220,222,224,226,228,230,232,234,236,238,240,242,244,246,248,250,252,254,256,258,260,262,264,266,268,270,272,274,276,278,280,282,284,286,288,290,292,294,296,298,
	   						 300,302,304,306,308,310,312,314,316,318,320,322,324,326,328,330,332,334,336,338,340,342,344,346,348,350,352,354,356,358,360,362,364,366,368,370,372,374,376,378,380,382,384,386,388,390,392,394,396,398,400,410,420,430,440,450,460,470,480,490,500,
							 510,520,530,540,550,560,570,580,590,600,610,620,630,640,650,660,670,680,690,700,710,720,730,740,750,760,770,780,790,800,810,820,830,840,850,860,870,880,890,900,910,920,930,940,950,960,970,980,990,1000,1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,
							 1110,1120,1130,1140,1150,1160,1170,1180,1190,1200,1210,1220,1230,1240,1250,1260,1270,1280,1290,1300,1310,1320,1330,1340,1350,1360,1370,1380,1390,1400,1410,1420,1430,1440,1450,1460,1470,1480,1490,1500};
   
   	TH2F *result = new TH2F(orig->GetName(),orig->GetTitle(),nbinsfinal,xbinfinal,40,0.,1.);
    
   	for(int binx=0;binx<result->GetXaxis()->GetNbins();binx++)
	{
		int b = tempt->GetXaxis()->FindBin(result->GetXaxis()->GetBinCenter(binx+1));
      	for(int biny=0;biny<40;biny++)
		{
			float bc = tempt->GetBinContent(b,biny+1);
			result->SetBinContent(binx+1,biny+1,bc);
      	}
    }

   	result->Draw("colz");
   	sprintf(filename,"template/plots/%s_%s_%s_%d.png",orig->GetName(),pname,jet_pt_cut.c_str(),year);
   	gPad->Print(filename);
   	return result;    
}

void plotterAndTemplateMaker(int year = 2016, int useMCatNLO = 1)
{
    //useMCatNLO = 0 : use just POWHEG
    //useMCatNLO = 1 : use just aMCatNLO
    //useMCatNLO = 2 : use aMCatNLO for shape, POWHEG for integral

    float lumi = 35.9;
    if (year == 2017) lumi = 41.5;
	if (year == 2018) lumi = 59.7;

	static const int vars = 16;
    string titlex[vars] = {"K_{D}","M_{4l} [GeV]","M_{jj} [GeV]","#Delta #eta_{jj}","p_{T,j}","#eta_{j}","#eta(j_1)","#eta(j_2)","p_T(j_1)","p_T(j_2)","sum(#eta_j)","m_{jj}/#Delta#eta_{jj}","#eta(Z_{1}*)","#eta(Z_{2}*)","R(p_{T}^{hard})","R(p_{T}^{jet})"};        
    string titley[vars] = {"Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin"};     
	string namegif[vars] = {"Dbkgkin","m4l","mjj","detajj","ptj","etaj","eta_j1","eta_j2","pt_jet1","pt_jet2","eta_j_sum","mjj_over_detajj","eta_Z1_star","eta_Z2_star","R_pt_hard","R_pt_jet"};
    int bins[vars] = {20,20,20,20,30,20,20,20,30,30,20,30,20,20,30,30};
    float xmin[vars] = {0.,0.,100.,0.,0.,-5.,-5.,-5.,25.,25.,-8.,-5.,-6.,-6.,0.,0.};
	float xmax[vars] = {1.,1400.,1000.,8.,300.,5.,5.,5.,600.,600.,8.,400.,6.,6.,1.,1.};	
	bool drawSignal[vars] = {false,false,true,true,true,true,true,true,true,true,true,true,true,true,true,true};

	//histogram stack
    char filename[300]; char filetitle[300];
	THStack *hs[vars];
    for(int iv = 0; iv < vars; iv++)
	{
	  	sprintf(filename,"hs%d",iv);   
	  	sprintf(filetitle,"CMS Preliminary                %2.1f fb^{-1};%s;%s",lumi,titlex[iv].c_str(),titley[iv].c_str());  
	  	hs[iv] = new THStack(filename,filetitle);
	}  

    lumi *= 1000;

    //histograms


	TH1F *hdata[vars]; //data, because we are hiding higher energies in this phase
	TH1F *hqqzz_powheg[vars]; //ew
	TH1F *hqqzz[vars]; //ew+zx -> real ew
	TH1F *hggzz[vars]; //gg
	TH1F *hvbs[vars]; //vbs
	TH1F *hsum1[vars]; //gg+ew -> real gg
	TH1F *hsum2[vars]; //gg+ew+vbs
	TH1F *hqqzz_ee[vars]; //qqzz e
	TH1F *hqqzz_mm[vars]; //qqzz mu
	TH1F *hqqzz_em[vars]; //qqzz e mu
	TH1F *hggzz_ee[vars];//ggzz e
	TH1F *hggzz_mm[vars];//ggzz mu
	TH1F *hggzz_em[vars]; //ggzz e mu
	TH1F *hvbs_ee[vars];//vbs e
	TH1F *hvbs_mm[vars];//vbs mu
	TH1F *hvbs_em[vars];//vbs e mu
	TH1F *hdata_ee[vars];//vbs e
	TH1F *hdata_mm[vars];//vbs mu
	TH1F *hdata_em[vars];//vbs e mu
	TH1F *httzwwz[vars];//ttz + wwz
	TH1F *httzwwz_ee[vars];//ttz + wwz e
	TH1F *httzwwz_mm[vars];//ttz + wwz mu
	TH1F *httzwwz_em[vars];//ttz + wwz e mu

	//------------------------------------------------------------------------------- my histograms ------------------------------------------------------------------------------

	TH1F *hZ1Mass_duje = new TH1F ("Z1Mass_duje", "", 40, 50, 130);
	TH1F *hZ1Mass = new TH1F ("Z1Mass", "", 40, 50, 130);
	TH1F *hZ2Mass_duje = new TH1F ("Z2Mass_duje", "", 40, 50, 130);
	TH1F *hZ2Mass = new TH1F ("Z2Mass", "", 40, 50, 130);
	TH1F *hZMass_duje = new TH1F ("ZMass_duje", "", 40, 50, 130);
	TH1F *hZMass = new TH1F ("ZMass", "", 40, 50, 130);
	TH1F *hZ1M_difference = new TH1F ("hZ1M_difference", "", 100, -2, 2);
	TH1F *hZ2M_difference = new TH1F ("hZ2M_difference", "", 100, -2, 2);


	//--------------------------------------------------------------------------- end of my histograms ---------------------------------------------------------------------------

	//TF_8 and TF_9 files
	bool calculate_aQGC_limits = true;
	//TH1F *hewk_FT8, *hewk_FT9;

	//if (calculate_aQGC_limits)
	//{
		Float_t bins_FT[] = { 0, 200, 300, 400, 500, 600, 800, 1000, 1200, 1400 };

		TFile *f_ewk_FT8 = new TFile("aQGC/ewk_FT8.root");
		TH1F *hewk_FT8 = new TH1F("ewk_FT8","ewk_FT8", 9, bins_FT);
		hewk_FT8 = (TH1F*)f_ewk_FT8->Get("BLS_hvbs_1_rescaled_FT8_1");

		TFile *f_ewk_FT9 = new TFile("aQGC/ewk_FT9.root");
		TH1F *hewk_FT9 = new TH1F("ewk_FT9","ewk_FT9", 9, bins_FT);
		hewk_FT9 = (TH1F*)f_ewk_FT9->Get("BLS_hvbs_1_rescaled_FT9_2");

	//}
	

	for(int iv = 0; iv < vars; iv++)
	{
	  	sprintf(filename,"hdata_%d",iv);   hdata[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //data, because we are hiding higher energies in this phase
	  	sprintf(filename,"hqqzz_powheg_%d",iv);   hqqzz_powheg[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //ew
	  	sprintf(filename,"hqqzz_%d",iv);   hqqzz[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //ew+zx -> real ew
	  	sprintf(filename,"hggzz_%d",iv);   hggzz[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //gg
	  	sprintf(filename,"hvbs_%d",iv);   hvbs[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //vbs
	  	sprintf(filename,"hsum1_%d",iv);   hsum1[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //gg+ew -> real gg
	  	sprintf(filename,"hsum2_%d",iv);   hsum2[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //gg+ew+vbs
	  	sprintf(filename,"hqqzz_ee_%d",iv);   hqqzz_ee[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //qqzz e
	  	sprintf(filename,"hqqzz_mm_%d",iv);   hqqzz_mm[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //qqzz mu
	  	sprintf(filename,"hqqzz_em_%d",iv);   hqqzz_em[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //qqzz e mu
	  	sprintf(filename,"hggzz_ee_%d",iv);   hggzz_ee[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//ggzz e
	  	sprintf(filename,"hggzz_mm_%d",iv);   hggzz_mm[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//ggzz mu
	  	sprintf(filename,"hggzz_em_%d",iv);  hggzz_em[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //ggzz e mu
	  	sprintf(filename,"hvbs_ee_%d",iv); hvbs_ee[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//vbs e
	  	sprintf(filename,"hvbs_mm_%d",iv); hvbs_mm[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//vbs mu
	  	sprintf(filename,"hvbs_em_%d",iv); hvbs_em[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//vbs e mu
	  	sprintf(filename,"hdataee_%d",iv); hdata_ee[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//vbs e
	  	sprintf(filename,"hdatamm_%d",iv); hdata_mm[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//vbs mu
	  	sprintf(filename,"hdataem_%d",iv); hdata_em[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//vbs e mu
	  	sprintf(filename,"httzwwz_%d",iv); httzwwz[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//for 2018 rescale
	  	sprintf(filename,"httzwwz_ee_%d",iv); httzwwz_ee[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//ttzwwz e
	  	sprintf(filename,"httzwwz_mm_%d",iv); httzwwz_mm[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//ttzwwz mu
	  	sprintf(filename,"httzwwz_em_%d",iv); httzwwz_em[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//ttzwwz e mu

		if (iv == 1)
		{
			Float_t bins_m4l[] = { 0, 200, 300, 400, 500, 600, 800, 1000, 1200, 1400 };

	  		sprintf(filename,"hdata_%d",iv);   hdata[iv] = new TH1F(filename,"",9,bins_m4l); //data, because we are hiding higher energies in this phase
	  		sprintf(filename,"hqqzz_powheg_%d",iv);   hqqzz_powheg[iv] = new TH1F(filename,"",9,bins_m4l); //ew
	  		sprintf(filename,"hqqzz_%d",iv);   hqqzz[iv] = new TH1F(filename,"",9,bins_m4l); //ew+zx -> real ew
	  		sprintf(filename,"hggzz_%d",iv);   hggzz[iv] = new TH1F(filename,"",9,bins_m4l); //gg
	  		sprintf(filename,"hvbs_%d",iv);   hvbs[iv] = new TH1F(filename,"",9,bins_m4l); //vbs
	  		sprintf(filename,"hsum1_%d",iv);   hsum1[iv] = new TH1F(filename,"",9,bins_m4l); //gg+ew -> real gg
	  		sprintf(filename,"hsum2_%d",iv);   hsum2[iv] = new TH1F(filename,"",9,bins_m4l); //gg+ew+vbs
	  		sprintf(filename,"hqqzz_ee_%d",iv);   hqqzz_ee[iv] = new TH1F(filename,"",9,bins_m4l); //qqzz e
	  		sprintf(filename,"hqqzz_mm_%d",iv);   hqqzz_mm[iv] = new TH1F(filename,"",9,bins_m4l); //qqzz mu
	  		sprintf(filename,"hqqzz_em_%d",iv);   hqqzz_em[iv] = new TH1F(filename,"",9,bins_m4l); //qqzz e mu
	  		sprintf(filename,"hggzz_ee_%d",iv);   hggzz_ee[iv] = new TH1F(filename,"",9,bins_m4l);//ggzz e
	  		sprintf(filename,"hggzz_mm_%d",iv);   hggzz_mm[iv] = new TH1F(filename,"",9,bins_m4l);//ggzz mu
	  		sprintf(filename,"hggzz_em_%d",iv);  hggzz_em[iv] = new TH1F(filename,"",9,bins_m4l); //ggzz e mu
	  		sprintf(filename,"hvbs_ee_%d",iv); hvbs_ee[iv] = new TH1F(filename,"",9,bins_m4l);//vbs e
	  		sprintf(filename,"hvbs_mm_%d",iv); hvbs_mm[iv] = new TH1F(filename,"",9,bins_m4l);//vbs mu
	  		sprintf(filename,"hvbs_em_%d",iv); hvbs_em[iv] = new TH1F(filename,"",9,bins_m4l);//vbs e mu
	  		sprintf(filename,"hdataee_%d",iv); hdata_ee[iv] = new TH1F(filename,"",9,bins_m4l);//vbs e
	  		sprintf(filename,"hdatamm_%d",iv); hdata_mm[iv] = new TH1F(filename,"",9,bins_m4l);//vbs mu
	  		sprintf(filename,"hdataem_%d",iv); hdata_em[iv] = new TH1F(filename,"",9,bins_m4l);//vbs e mu
	  		sprintf(filename,"httzwwz_%d",iv); httzwwz[iv] = new TH1F(filename,"",9,bins_m4l);//for 2018 rescale
	  		sprintf(filename,"httzwwz_ee_%d",iv); httzwwz_ee[iv] = new TH1F(filename,"",9,bins_m4l);//ttzwwz e
	  		sprintf(filename,"httzwwz_mm_%d",iv); httzwwz_mm[iv] = new TH1F(filename,"",9,bins_m4l);//ttzwwz mu
	  		sprintf(filename,"httzwwz_em_%d",iv); httzwwz_em[iv] = new TH1F(filename,"",9,bins_m4l);//ttzwwz e mu
			  
		}	
	}   
	gStyle->SetPalette(1);
	TFile *input_file;
    int nbins=17;
	double xbin[18] = {160,166,170,176,182,188,194,200,208,216,224,234,246,260,278,302,338,1500}; 

	float c_constant = 8.5;
    // if (year == 2017) c_constant = 3.5;
	// if (year == 2018) c_constant = 3.5; 
	TFile* f_ = TFile::Open("/home/llr/cms/giljanovic/scratch/CMSSW_10_2_15/src/ZZAnalysis/AnalysisStep/data/cconstants/SmoothKDConstant_m4l_DjjVBF13TeV.root");
	TSpline3* ts = (TSpline3*)(f_->Get("sp_gr_varReco_Constant_Smooth")->Clone());
	f_->Close();

    // find available samples  
	int nSamp = 0;  
	TString rootname[40];
 	sprintf(filename,"newsamples%d_withMGggZZandVBS.txt",year);

	ifstream parInput(filename);
        
	if (parInput.is_open()) 
	{
	  	while ( parInput.good() ) 
		{
	  	  	parInput >> rootname[nSamp]; 
      	    cout << nSamp << " " <<  rootname[nSamp] << endl; 
	  	  	nSamp++;
	  	}
	  	parInput.close();
	} 

	//original variable declarations
	float ZZPt,ZZMass,Z1Mass,Z2Mass,DiJetMass,DiJetDEta;
	float xsec,KFactorEWKqqZZ,overallEventWeight,L1prefiringWeight,genHEPMCweight,KFactorQCDqqZZ_M;
	vector<float> *LepPt=new vector<float>;
	vector<float> *JetPt=new vector<float>;
	vector<float> *JetEta=new vector<float>;
	vector<float> *JetPhi=new vector<float>;
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

	//additional and output variable declarations
	float weight, weight_up, weight_dn;
	float weight_vbf, weight_vbf_up, weight_vbf_dn;
	int chan;
	int vbfcate = 0;
	float dbkg_kin, theVar; 
	
	//output branches
	TTree *tnew[4]; 
	TFile *fnew[4];
	//template declarations (1D) 
	TH1F *temp_1d_4e[4];
	TH1F *temp_1d_4mu[4];
	TH1F *temp_1d_2e2mu[4];
	//template declarations (2D) 
	TH2F *temp_zz_4e[4];
	TH2F *temp_zz_4mu[4];
	TH2F *temp_zz_2e2mu[4];

	for (int it=0; it < 4; it++) 
	{
	  	if (it==0) sprintf(filename,"template/root_output_files/qqzz_Moriond_%s_%d.root",jet_pt_cut.c_str(), year); 
	  	if (it==1) sprintf(filename,"template/root_output_files/ggzz_Moriond_%s_%d.root",jet_pt_cut.c_str(), year); 
	  	if (it==2) sprintf(filename,"template/root_output_files/vbs_Moriond_%s_%d.root",jet_pt_cut.c_str(), year); 
	  	if (it==3) sprintf(filename,"template/root_output_files/data_%s_%d.root",jet_pt_cut.c_str(), year); 
	  	fnew[it] = new TFile(filename,"recreate");
	  	tnew[it] = new TTree("SelectedTree","SelectedTree");
	  	tnew[it]->Branch("mreco",&ZZMass,"mreco/F");
	  	tnew[it]->Branch("dbkg_kin",&dbkg_kin,"dbkg_kin/F");
	  	tnew[it]->Branch("weight",&weight,"weight/F");
	  	tnew[it]->Branch("weight_up",&weight_up,"weight_up/F");
	  	tnew[it]->Branch("weight_dn",&weight_dn,"weight_dn/F");
	  	tnew[it]->Branch("chan",&chan,"chan/I");
	  	tnew[it]->Branch("vbfcate",&vbfcate,"vbfcate/I");
	  	temp_zz_4e[it] = new TH2F("temp_zz_4e","",nbins,xbin,40,0.,1.);
	  	temp_zz_4mu[it] = new TH2F("temp_zz_4mu","",nbins,xbin,40,0.,1.);
	  	temp_zz_2e2mu[it] = new TH2F("temp_zz_2e2mu","",nbins,xbin,40,0.,1.);
	  	temp_1d_4e[it] = new TH1F("temp_1d_4e","",50,0.,1.);
	  	temp_1d_4mu[it] = new TH1F("temp_1d_4mu","",50,0.,1.);
	  	temp_1d_2e2mu[it] = new TH1F("temp_1d_2e2mu","",50,0.,1.);
	}

	// ------------------------------------------------------------------------------------------------- my declarations ---------------------------------------------------------------------------------------------

	vector<short> *lepId = new vector<short>;
	vector<float> *lepPt = new vector<float>;
	vector<float> *lepEta = new vector<float>;
	vector<float> *lepPhi = new vector<float>;

	vector<TLorentzVector> electrons;
	vector<TLorentzVector> muons;
	vector<short> electrons_charge;
	vector<short> muons_charge;

	TLorentzVector Z1, Z2;

	float eta_Z1_star, eta_Z2_star, R_pt_hard, R_pt_jet;

	// --------------------------------------------------------------------------------------------- end of my declarations ------------------------------------------------------------------------------------------

	
	//for loop for different samples
	for(int is = 0; is < nSamp-1; is++)
	{

	  	//print cycle
	  	std::cout << endl << is << endl;
	  
        TFile* input_file = TFile::Open(rootname[is].Data()); 
	  	TH1F *hCounters= (TH1F*)input_file->Get("ZZTree/Counters");
	  	float gen_sum_weights = hCounters->GetBinContent(40);
	  	std::cout<<endl<<is<<"  "<< gen_sum_weights<<endl;
	  
	  	//tchain and add function for multiple input files
	  	TChain *tqqzz= new TChain("ZZTree/candTree");
	  	tqqzz->Add(rootname[is].Data());
          
        //process class
        int j = 0;   // qqzz powheg
	  	
		if (rootname[is].Contains("AllData")) j = 3;   
	  	if (rootname[is].Contains("ggTo") || rootname[is].Contains("ggZZnew")) j = 1;      
        if (rootname[is].Contains("VBFTo")) j = 2;
	  	if (rootname[is].Contains("amcatnlo")) j = 4;
	  	//if (rootname[is].Contains("WWZ") || rootname[is].Contains("TTZ")) j = 5;   
	  	//histogram declaration
	  	//TH1F *kin_zz = new TH1F("kin_zz","",bins,xmin,xmax); //was 100 bins
	  
	
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
	  	tqqzz->SetBranchAddress("Z1Mass",&Z1Mass);
	  	tqqzz->SetBranchAddress("Z2Mass",&Z2Mass);
	  	tqqzz->SetBranchAddress("KFactor_EW_qqZZ",&KFactorEWKqqZZ);
	  	tqqzz->SetBranchAddress("xsec",&xsec);
	  	tqqzz->SetBranchAddress("ZZsel",&ZZsel);
	  	tqqzz->SetBranchAddress("LepLepId",&LepLepId);
	  	tqqzz->SetBranchAddress("LepPt",&LepPt);
        tqqzz->SetBranchAddress("JetPt",&JetPt);
        tqqzz->SetBranchAddress("JetEta",&JetEta);
		tqqzz->SetBranchAddress("JetPhi",&JetPhi);
	  	tqqzz->SetBranchAddress("overallEventWeight",&overallEventWeight);
	  	tqqzz->SetBranchAddress("genHEPMCweight",&genHEPMCweight);
	  	tqqzz->SetBranchAddress("L1prefiringWeight",&L1prefiringWeight);
	  	tqqzz->SetBranchAddress("KFactor_QCD_qqZZ_M",&KFactorQCDqqZZ_M);
	  	tqqzz->SetBranchAddress("nCleanedJetsPt30",&nCleanedJetsPt30);
	  
	  	//new branch addresses
	  	tqqzz->SetBranchAddress("p_JJEW_BKG_MCFM_JECNominal",&p_JJEW_BKG_MCFM_JECNominal);
	  	tqqzz->SetBranchAddress("p_JJVBF_BKG_MCFM_JECNominal",&p_JJVBF_BKG_MCFM_JECNominal);
	  	tqqzz->SetBranchAddress("p_JJQCD_BKG_MCFM_JECNominal",&p_JJQCD_BKG_MCFM_JECNominal);
	  	tqqzz->SetBranchAddress("DiJetMass",&DiJetMass);
        tqqzz->SetBranchAddress("DiJetDEta",&DiJetDEta);
	  	tqqzz->SetBranchAddress("KFactor_QCD_ggZZ_Nominal",&KFactorQCDggzz_Nominal);


		// --------------------------------------------------------------------------------- my branches ---------------------------------------------------------------

		tqqzz->SetBranchAddress("LepLepId",&lepId);
		tqqzz->SetBranchAddress("LepPt",&lepPt);
		tqqzz->SetBranchAddress("LepEta",&lepEta);
		tqqzz->SetBranchAddress("LepPhi",&lepPhi);

		// ----------------------------------------------------------------------------- end of my branches ------------------------------------------------------------


	  	//loop on entries
	  	int enne = tqqzz->GetEntries();

	  	// preliminary loop to fix MG wrong weights (only VBS 2017-18)
	  	float resum = gen_sum_weights;
      	if (j==2 && year>2016) resum = hCounters->GetBinContent(1);
  
	  	for(int i=0;i<enne;i++)
		{
	    	tqqzz->GetEntry(i);
	    
	   	    //unique selection condition (see paper page 8) & DiJetMass condition
	    	// if(DiJetMass>100 && nExtraLep==0 && ZZMass > 160 &&(((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0))){
	    
	    	if(DiJetMass>100 && ZZMass > 180 && nCleanedJetsPt30>1 && Z1Mass < 120 && Z1Mass > 60 && Z2Mass < 120 && Z2Mass > 60)
			{
				
	      		// ------------------------------------------------------------ construc electron and muon objects -----------------------------------------
				
				clear_vectors(electrons, muons, electrons_charge, muons_charge);
				create_electron_and_muon_objects(electrons, muons, electrons_charge, muons_charge, lepId, lepPt, lepEta, lepPhi);
				build_ZZ_pair(electrons, muons, electrons_charge, muons_charge, Z1, Z2);

				if (is == 0)
				{
					hZ1Mass->Fill(Z1Mass);
					hZ1Mass_duje->Fill(Z1.M());
					hZ2Mass->Fill(Z2Mass);
					hZ2Mass_duje->Fill(Z2.M());

					hZMass->Fill(Z1Mass);
					hZMass->Fill(Z2Mass);
					hZMass_duje->Fill(Z1.M());
					hZMass_duje->Fill(Z2.M());

					hZ1M_difference->Fill(Z1Mass-Z1.M());
					hZ2M_difference->Fill(Z2Mass-Z2.M());
				}
				
				// ------------------------------------------------------- end of construc electron and muon objects ---------------------------------------

				// ------------------------------------------------------------ calculate aditional variables ----------------------------------------------

				calculate_Zeppenfeld_Z(Z1, Z2, JetEta, eta_Z1_star, eta_Z2_star);
				calculate_R_pt_hard(Z1, Z2, JetEta, JetPhi, JetPt, R_pt_hard);
				calculate_R_pt_jet(JetEta, JetPhi, JetPt, R_pt_jet);

				// -------------------------------------------------------- end of calculate aditional variables -------------------------------------------
				  
	      		//set vbf_category
	      		vbfcate=1;

	      		weight= (xsec*KFactorEWKqqZZ*overallEventWeight*KFactorQCDqqZZ_M*L1prefiringWeight*lumi)/(resum);
	      		// correct k-factor for NNLO/NLO?
	      		// if (j==1) weight= (xsec*overallEventWeight*KFactorQCDggzz_Nominal*L1prefiringWeight*lumi)/(resum);
	      		// if (j==1 && useMCatNLO==1) weight /=1.7;
	      		if (j==1) weight= (xsec*overallEventWeight*1.3*L1prefiringWeight*lumi)/(resum);
          		if (j==2 && year==2016) weight= (xsec*overallEventWeight*L1prefiringWeight*lumi)/(resum);
          		if (j==2 && year>2016) weight= (xsec*overallEventWeight*L1prefiringWeight*lumi)/(genHEPMCweight*resum);             
	      		if (j==5) weight= (xsec*overallEventWeight*L1prefiringWeight*lumi)/(resum);
	      		if (j==3) weight=1.; 

	      		//TEMPORARY FOR MISSING 2e2mu SAMPLE
	      		//if (j==2 && year==2017) weight *= 2.;

	      		//division in channels
	      		if(abs(Z1Flav)==abs(Z2Flav) && abs(Z1Flav)==121) chan=2;
	      		else if (abs(Z1Flav)==abs(Z2Flav) && abs(Z1Flav)!=121) chan=1;
	      		else chan=3;

	      		//kin variable
            	float c_mzz = c_constant*ts->Eval(ZZMass);
	      		dbkg_kin = p_JJVBF_BKG_MCFM_JECNominal/(p_JJVBF_BKG_MCFM_JECNominal+ p_JJQCD_BKG_MCFM_JECNominal*c_mzz);
	      		if (dbkg_kin < 0.00 || dbkg_kin > 1.00) continue;

	      		// fill templates
	      		if ((useMCatNLO > 0 && j==4) || j==5) 
				{
					tnew[0]->Fill();
					if(chan==1) { temp_zz_4mu[0]->Fill(ZZMass,dbkg_kin,weight); temp_1d_4mu[0]->Fill(dbkg_kin,weight); } 
					else if(chan==2) { temp_zz_4e[0]->Fill(ZZMass,dbkg_kin,weight); temp_1d_4e[0]->Fill(dbkg_kin,weight); }
	        		else { temp_zz_2e2mu[0]->Fill(ZZMass,dbkg_kin,weight); temp_1d_2e2mu[0]->Fill(dbkg_kin,weight); }
	      		}
				else if ((useMCatNLO == 0 && j==0) || j==5) 
				{
					tnew[0]->Fill();
					if(chan==1) { temp_zz_4mu[0]->Fill(ZZMass,dbkg_kin,weight); temp_1d_4mu[0]->Fill(dbkg_kin,weight); }
					else if(chan==2) { temp_zz_4e[0]->Fill(ZZMass,dbkg_kin,weight);  temp_1d_4e[0]->Fill(dbkg_kin,weight); }
	        		else { temp_zz_2e2mu[0]->Fill(ZZMass,dbkg_kin,weight); temp_1d_2e2mu[0]->Fill(dbkg_kin,weight); }
	      		}
		  		else if (j==1) 
				{
					tnew[1]->Fill();
					temp_zz_4mu[1]->Fill(ZZMass,dbkg_kin,weight);  temp_1d_4mu[1]->Fill(dbkg_kin,weight);
					temp_zz_4e[1]->Fill(ZZMass,dbkg_kin,weight);  temp_1d_4e[1]->Fill(dbkg_kin,weight);
	        		temp_zz_2e2mu[1]->Fill(ZZMass,dbkg_kin,weight);  temp_1d_2e2mu[1]->Fill(dbkg_kin,weight);
	      		} 
				else 
				{
					tnew[j]->Fill();
					if(chan==1) { temp_zz_4mu[j]->Fill(ZZMass,dbkg_kin,weight);  temp_1d_4mu[j]->Fill(dbkg_kin,weight); }
					else if(chan==2) { temp_zz_4e[j]->Fill(ZZMass,dbkg_kin,weight); temp_1d_4e[j]->Fill(dbkg_kin,weight); } 
	        		else { temp_zz_2e2mu[j]->Fill(ZZMass,dbkg_kin,weight);  temp_1d_2e2mu[j]->Fill(dbkg_kin,weight); }
	    		}

            	// fill plots
	      		for(int il = 0; il < vars+2; il++)
				{
            	    int iv = il; 
					if (il == 0) theVar = dbkg_kin;
					if (il == 1) theVar = ZZMass;
					if (il == 2) theVar = DiJetMass;
					if (il == 3) theVar = fabs(DiJetDEta);
            	    if (il == 4) theVar = JetPt->at(0);
            	    if (il == 5) theVar = JetEta->at(0);
            	    if (il == 6) 
					{
						theVar = JetPt->at(1);   
						iv = 4;
					}
					if (il == 7) 
					{
						theVar = JetEta->at(1);   
						iv = 5;
					}
					if (il == 8)
					{
						theVar = JetEta->at(0);
						iv = 6;
					}
					if (il == 9)
					{
						theVar = JetEta->at(1);
						iv = 7;
					}
					if (il == 10)
					{
						theVar = JetPt->at(0);
						iv = 8;
					}
					if (il == 11)
					{
						theVar = JetPt->at(1);
						iv = 9;
					}
					if (il == 12)
					{
						theVar = JetEta->at(0) + JetEta->at(1);
						iv = 10;
					}
					if (il == 13)
					{
						theVar = ZZMass/fabs(DiJetDEta);
						iv = 11;
					}
					if (il == 14)
					{
						theVar = eta_Z1_star;
						iv = 12;
					}
					if (il == 15)
					{
						theVar = eta_Z2_star;
						iv = 13;
					}
					if (il == 16)
					{
						theVar = R_pt_hard;
						iv = 14;
					}
					if (il == 17)
					{
						theVar = R_pt_jet;
						iv = 15;
					}
	
					//1D kin var hist fill
					//this is the normalization histogram
					// if (j<3){
					//  if (theVar <= 1.1) h00[iv]->Fill(theVar,weight);    
					// h_complete_data[iv]->Fill(theVar,weight);
					//}

					if (j==3)
					{
					  	if (iv > 0 || theVar < 0.75 || year == 2016) hdata[iv]->Fill(theVar);
					  	else cout << "Blinded event!" << endl;
					  	hdata[iv]->SetMarkerStyle(20);
					  	if (chan == 2) hdata_ee[iv]->Fill(theVar);
					  	if (chan == 1) hdata_mm[iv]->Fill(theVar);
					  	if (chan == 3) hdata_em[iv]->Fill(theVar);
					}
					if (j==0) //qqzz
					{                              
					  hqqzz_powheg[iv]->Fill(theVar,weight); 

					  /* if (chan == 1) hqqzz_mm[iv]->Fill(theVar,weight); //mu
					  if (chan == 2) hqqzz_ee[iv]->Fill(theVar,weight); //e
					  if (chan == 3)  hqqzz_em[iv]->Fill(theVar,weight);  */        //mu+e
					}
					if (j==1) //gg
					{                              
					  	hggzz[iv]->Fill(theVar,weight);
					  	hggzz[iv]->SetFillColor(kBlue);

					  	if (chan == 1) hggzz_mm[iv]->Fill(theVar,weight); //mu
					  	if (chan == 2) hggzz_ee[iv]->Fill(theVar,weight); //e
					  	if (chan == 3)  hggzz_em[iv]->Fill(theVar,weight);          //mu+e

					}
					if (j==2) //vbs
					{                              
					  	hvbs[iv]->Fill(theVar,weight);
					  	hvbs[iv]->SetFillColor(kMagenta);

					  	if (chan == 1) hvbs_mm[iv]->Fill(theVar,weight); //mu
					  	if (chan == 2) hvbs_ee[iv]->Fill(theVar,weight); //e
					  	if (chan == 3) hvbs_em[iv]->Fill(theVar,weight);          //mu+e

					}
					if (j==4) 
					{
					  	hqqzz[iv]->Fill(theVar,weight);

					  	if (chan == 1) hqqzz_mm[iv]->Fill(theVar,weight); //mu
					  	if (chan == 2) hqqzz_ee[iv]->Fill(theVar,weight); //e
					  	if (chan == 3)  hqqzz_em[iv]->Fill(theVar,weight);
					}
					if (j==5) 
					{
					  	httzwwz[iv]->Fill(theVar,weight);

					  	if (chan == 1) httzwwz_mm[iv]->Fill(theVar,weight); //mu
					  	if (chan == 2) httzwwz_ee[iv]->Fill(theVar,weight); //e
					  	if (chan == 3) httzwwz_em[iv]->Fill(theVar,weight);          //mu+e
					}

					//kin_zz->GetYaxis()->SetTitle("Events/0.05");
					//add histogram to stack
					//for cycle ends here
	    		}
	    	}
		}//entries loop  end
	}//file loop  end
       
	//ZX CONTRIBUTION
	  
	TChain *tqqzz_zx= new TChain("candTree");
	sprintf(filename,"/home/llr/cms/giljanovic/scratch/VBS/CMSSW_10_2_15/src/vbs_analysis/4l_channel/data_driven_MC/ZX%d_%s.root",year, jet_pt_cut.c_str()); 
	tqqzz_zx->Add(filename);
	
	//histogram declaration
	TH1F *hzx[vars]; 
	TH1F *hzx_ee[vars]; 
	TH1F *hzx_em[vars]; 
	TH1F *hzx_mm[vars]; 
	for(int iv = 0; iv < vars; iv++)
	{
	  	sprintf(filename,"hzx%d",iv);   	hzx[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //full histogram
	  	sprintf(filename,"hzx_ee%d",iv);   	hzx_ee[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //full histogram
	  	sprintf(filename,"hzx_mm%d",iv);   	hzx_mm[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //full histogram
	  	sprintf(filename,"hzx_em%d",iv);   	hzx_em[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //full histogram

		if (iv == 1)
		{
			Float_t bins_m4l[] = { 0, 200, 300, 400, 500, 600, 800, 1000, 1200, 1400 };

			sprintf(filename,"hzx%d",iv);   	hzx[iv] = new TH1F(filename,"",9,bins_m4l); //full histogram
	  		sprintf(filename,"hzx_ee%d",iv);   	hzx_ee[iv] = new TH1F(filename,"",9,bins_m4l); //full histogram
	  		sprintf(filename,"hzx_mm%d",iv);   	hzx_mm[iv] = new TH1F(filename,"",9,bins_m4l); //full histogram
	  		sprintf(filename,"hzx_em%d",iv);   	hzx_em[iv] = new TH1F(filename,"",9,bins_m4l); //full histogram
		}
	}
	
	//new variable declarations
	float var_zx, dbkg_kin_zx, ZZMass_zx, DiJetMass_zx, DiJetDEta_zx;
	float ptjet1_zx,ptjet2_zx,etajet1_zx,etajet2_zx; 
	float weight_zx;
	int chan_zx;
	
	//original branch addresses
	tqqzz_zx->SetBranchAddress("dbkg_kin",&dbkg_kin_zx);
	tqqzz_zx->SetBranchAddress("ZZMass",&ZZMass_zx);
	tqqzz_zx->SetBranchAddress("DiJetMass",&DiJetMass_zx);
	tqqzz_zx->SetBranchAddress("DiJetDEta",&DiJetDEta_zx);
	tqqzz_zx->SetBranchAddress("weight",&weight_zx);
	tqqzz_zx->SetBranchAddress("ptjet1",&ptjet1_zx); 
	tqqzz_zx->SetBranchAddress("ptjet2",&ptjet2_zx);   
	tqqzz_zx->SetBranchAddress("etajet1",&etajet1_zx); 
	tqqzz_zx->SetBranchAddress("etajet2",&etajet2_zx);
	tqqzz_zx->SetBranchAddress("chan",&chan_zx);

	//entries loop
	for(int i=0;i<tqqzz_zx->GetEntries();i++)
	{
	  	tqqzz_zx->GetEntry(i);
	  
	  	//1D kin var hist fill
	  	for(int il = 0; il < vars+2; il++)
		{
            int iv = il; 
	    	if (iv == 0) var_zx = dbkg_kin_zx;
	    	if (iv == 1) var_zx = ZZMass_zx;
	    	if (iv == 2) var_zx = DiJetMass_zx;
	    	if (iv == 3) var_zx = fabs(DiJetDEta_zx);
	    	if (il == 4) var_zx = ptjet1_zx;
	    	if (il == 5) var_zx = etajet1_zx;
	    	if (il == 6) 
			{
				var_zx = ptjet2_zx;   
				iv = 4;
			}
	    	if (il == 7) 
			{
				var_zx = etajet2_zx;   
				iv = 5;
			}
			if (il == 8)
			{
				var_zx = etajet1_zx;
				iv = 6;
			}
			if (il == 9)
			{
				var_zx = etajet1_zx;
				iv = 7;
			}
			if (il == 10)
			{
				var_zx = ptjet1_zx;
				iv = 8;
			}
			if (il == 11)
			{
				var_zx = ptjet2_zx;
				iv = 9;
			}
			if (il == 12)
			{
				var_zx = etajet1_zx + etajet2_zx;
				iv = 10;
			}
			if (il == 13)
			{
				theVar = ZZMass_zx/fabs(DiJetDEta_zx);
				iv = 11;
			}			

	    	if (fabs(weight_zx) < 100000. && iv < 12) 
			{
	      		hzx[iv]->Fill(var_zx,weight_zx);
	      		if (chan_zx == 2) hzx_ee[iv]->Fill(var_zx,weight_zx);
	      		if (chan_zx == 1) hzx_mm[iv]->Fill(var_zx,weight_zx);
	      		if (chan_zx == 3) hzx_em[iv]->Fill(var_zx,weight_zx);
	    	}
	    	// if (var_zx<=0.7) h00[iv]->Fill(var_zx,weight_zx);
	  	}
	}
	
	for (int it=0; it < 4; it++) 
	{
	  	fnew[it]->cd();
	  	if (it < 3) 
		{
	    	TH2F* final_4e = rebinTemplate(temp_zz_4e[it],year,it);
	    	final_4e->Write();
	    	temp_1d_4e[it]->Write();
	    	TH2F* final_4mu = rebinTemplate(temp_zz_4mu[it],year,it);
	    	final_4mu->Write();
	    	temp_1d_4mu[it]->Write();
	    	TH2F* final_2e2mu = rebinTemplate(temp_zz_2e2mu[it],year,it);
	    	final_2e2mu->Write();
	    	temp_1d_2e2mu[it]->Write();
	  	}
	  	for(int il = 0; il < vars; il++)
		{
	    	if (it==0) {hqqzz[il]->Write();    httzwwz[il]->Write();}
	    	if (it==1) hggzz[il]->Write();
	    	if (it==2) hvbs[il]->Write();
	    	if (it==3) {hdata[il]->Write();    hzx[il]->Write();}
	  	}
	  	tnew[it]->Write();
	  	fnew[it]->Close();
	}
	
	//INTEGRAL CHECK
	sprintf(filename,"MCyields_%s_%d.txt",jet_pt_cut.c_str(), year);
	ofstream yields(filename,std::fstream::app);
    sprintf(filename,"datayields_%s_%d.txt",jet_pt_cut.c_str(), year);
	ofstream yields2(filename,std::fstream::app);
	sprintf(filename,"MCyields_highMELA_%s_%d.txt",jet_pt_cut.c_str(), year);
	ofstream yields3(filename,std::fstream::app);
    sprintf(filename,"datayields_highMELA_%s_%d.txt",jet_pt_cut.c_str(), year);
	ofstream yields4(filename,std::fstream::app);

	for(int iv = 0; iv < vars; iv++)
	{
	  	cout << "Integral check" << endl;
	  	cout <<"qqzz POWHEG,         integral is " << hqqzz_powheg[iv]->Integral() << endl;
        cout <<"qqzz,         integral is " << hqqzz[iv]->Integral() << endl;
	  	cout <<"qqzz (4e),    integral is " << hqqzz_ee[iv]->Integral() << endl;
	  	cout <<"qqzz (4mu),   integral is " << hqqzz_mm[iv]->Integral() << endl;
	  	cout <<"qqzz (2e2mu), integral is " << hqqzz_em[iv]->Integral() << endl;
	  	cout <<"ggzz,         integral is " << hggzz[iv]->Integral() << endl;
	  	cout <<"ggzz (4e),    integral is " << hggzz_ee[iv]->Integral() << endl;
	  	cout <<"ggzz (4mu),   integral is " << hggzz_mm[iv]->Integral() << endl;
	  	cout <<"ggzz (2e2mu), integral is " << hggzz_em[iv]->Integral() << endl;
	  	cout <<"vbs,          integral is " << hvbs[iv]->Integral() << endl;
	  	cout <<"vbs (4e),     integral is " << hvbs_ee[iv]->Integral() << endl;
	  	cout <<"vbs (4mu),    integral is " << hvbs_mm[iv]->Integral() << endl;
	  	cout <<"vbs (2e2mu),  integral is " << hvbs_em[iv]->Integral() << endl;
	  	cout <<"zx,          integral is " << hzx[iv]->Integral() << endl;
	  	cout <<"zx (4e),     integral is " << hzx_ee[iv]->Integral() << endl;
	  	cout <<"zx (4mu),    integral is " << hzx_mm[iv]->Integral() << endl;
	  	cout <<"zx (2e2mu),  integral is " << hzx_em[iv]->Integral() << endl;
	  	cout <<"data,          integral is " << hdata[iv]->Integral() << endl;
	  	cout <<"data (4e),     integral is " << hdata_ee[iv]->Integral() << endl;
	  	cout <<"data (4mu),    integral is " << hdata_mm[iv]->Integral() << endl;
	  	cout <<"data (2e2mu),  integral is " << hdata_em[iv]->Integral() << endl;
	  
	
	  	if (iv == 0) 
		{ 
	    	float vbsee  = hvbs_ee[iv]->Integral();
	    	float vbsemu = hvbs_em[iv]->Integral();
            float vbsmumu = hvbs_mm[iv]->Integral();
	    	float vbsee_high  = hvbs_ee[iv]->Integral(19,20);
	    	float vbsemu_high = hvbs_em[iv]->Integral(19,20);
            float vbsmumu_high = hvbs_mm[iv]->Integral(19,20);

	    	// TEMPORARY MISSING SAMPLE
	    	/* if (year == 2017) {
	    	  vbsemu_high = (vbsee_high+vbsmumu_high)/2.;
	    	  vbsee_high /= 2.;       vbsmumu_high /= 2.;
	    	  vbsemu = (vbsee+vbsmumu)/2.;
	    	  vbsee /= 2.;       vbsmumu /= 2.;
	    	  } */
	    
	    	yields << vbsmumu << "\t" << hzx_mm[iv]->Integral() << "\t" << hqqzz_mm[iv]->Integral()+httzwwz_mm[iv]->Integral() << "\t" << hggzz_mm[iv]->Integral() << "\t" << vbsee << "\t" << hzx_ee[iv]->Integral() << "\t" <<  hqqzz_ee[iv]->Integral()+httzwwz_ee[iv]->Integral() << "\t" << hggzz_ee[iv]->Integral() << "\t" << vbsemu << "\t" << hzx_em[iv]->Integral() << "\t" << hqqzz_em[iv]->Integral()+httzwwz_em[iv]->Integral() << "\t" << hggzz_em[iv]->Integral() << endl;
	    	yields2 << hdata_mm[iv]->Integral() << "\t" << hdata_ee[iv]->Integral() << "\t" << hdata_em[iv]->Integral() << endl;
	    	yields3 << vbsmumu_high << "\t" << hzx_mm[iv]->Integral(19,20) << "\t" << hqqzz_mm[iv]->Integral(19,20)+ httzwwz_mm[iv]->Integral(19,20) << "\t" << hggzz_mm[iv]->Integral(19,20) << "\t" << vbsee_high << "\t" << hzx_ee[iv]->Integral(19,20) << "\t" <<  hqqzz_ee[iv]->Integral(19,20)+httzwwz[iv]->Integral(19,20) << "\t" << hggzz_ee[iv]->Integral(19,20) << "\t" << vbsemu_high << "\t" << hzx_em[iv]->Integral(19,20) << "\t" << hqqzz_em[iv]->Integral(19,20)+httzwwz_em[iv]->Integral(19,20) << "\t" << hggzz_em[iv]->Integral(19,20) << endl;
	    	yields4 << hdata_mm[iv]->Integral(19,20) << "\t" << hdata_ee[iv]->Integral(19,20) << "\t" << hdata_em[iv]->Integral(19,20) << endl;
	  	}
	  
	  	float zx_integral = (2.00933+1.97468+3.96839)*lumi/35.9E03;
        const float powheg_integral = hqqzz_powheg[iv]->Integral();
	  
	  	/* if (year == 2018 && useMCatNLO > 0) {
	  	  hqqzz_powheg[iv]->Multiply(hqqzz_powheg[iv],hqqzz[iv]);
	  	  hqqzz_powheg[iv]->Divide(hqqzz_powheg[iv],httzwwz[iv]);
	  	  }  */
	  	//hzx[iv]->Scale(zx_integral/hzx[iv]->Integral());

	  	//data normalisation
	  	if (useMCatNLO == 2) 
		{
	    	// if (year == 2018) hqqzz_powheg[iv]->Scale(powheg_integral/hqqzz_powheg[iv]->Integral());
	    	// else 
	    	hqqzz[iv]->Scale(powheg_integral/hqqzz[iv]->Integral());
	  	}
	  
	  	//HISTOGRAMS ADDED TO STACK
	  	hzx[iv]->SetFillColor(kGreen);
        httzwwz[iv]->Add(httzwwz[iv],hzx[iv],1,1);    //tt
	  	httzwwz[iv]->SetFillColor(kYellow);
	  	if (useMCatNLO == 0) hqqzz[iv]->Add(httzwwz[iv],hqqzz_powheg[iv],1,1); //real ew
	  	else hqqzz[iv]->Add(httzwwz[iv],hqqzz[iv],1,1); //real ew
	  	hqqzz[iv]->SetFillColor(kCyan);
	  	hsum1[iv]->Add(hqqzz[iv],hggzz[iv],1,1); //real gg
	  	hsum1[iv]->SetFillColor(kBlue);
	  	hsum2[iv]->Add(hsum1[iv],hvbs[iv],1,1);    //real vbs
	  	hsum2[iv]->SetFillColor(kMagenta);





		//INCLUDING OVERFLOW INTO THE LAST BIN (for m4l plot)
		if (iv == 1)
		{
			hsum2[iv]->SetBinContent(hsum2[iv]->GetNbinsX(), hsum2[iv]->GetBinContent(hsum2[iv]->GetNbinsX()) + hsum2[iv]->GetBinContent(hsum2[iv]->GetNbinsX() + 1));	//gg+ew+vbs (magenta -> signal)
			hsum1[iv]->SetBinContent(hsum1[iv]->GetNbinsX(), hsum1[iv]->GetBinContent(hsum1[iv]->GetNbinsX()) + hsum1[iv]->GetBinContent(hsum1[iv]->GetNbinsX() + 1));	//gg+ew ->real gg (blue)
			hqqzz[iv]->SetBinContent(hqqzz[iv]->GetNbinsX(), hqqzz[iv]->GetBinContent(hqqzz[iv]->GetNbinsX()) + hqqzz[iv]->GetBinContent(hqqzz[iv]->GetNbinsX() + 1));	//real ew (cyan)
			hzx[iv]->SetBinContent(hzx[iv]->GetNbinsX(), hzx[iv]->GetBinContent(hzx[iv]->GetNbinsX()) + hzx[iv]->GetBinContent(hzx[iv]->GetNbinsX() + 1));	//full Z+X (green)
			hdata[iv]->SetBinContent(hdata[iv]->GetNbinsX(), hdata[iv]->GetBinContent(hdata[iv]->GetNbinsX()) + hdata[iv]->GetBinContent(hdata[iv]->GetNbinsX() + 1));	//data
		}


		//saving EWK, qq, gg and data histogram in root file for aQGC part
		if (iv == 1 && calculate_aQGC_limits == false)
		{

			TFile *ewk_hist = new TFile("ewk.root","recreate");
			hvbs[iv]->Write();
			ewk_hist->Write();
			ewk_hist->Close();

			TFile *qq_hist = new TFile("qq.root","recreate");
			hqqzz[iv]->Write();
			qq_hist->Write();
			qq_hist->Close();

			TFile *gg_hist = new TFile("gg.root","recreate");
			hggzz[iv]->Write();
			gg_hist->Write();
			gg_hist->Close();

			TFile *data_hist = new TFile("data.root","recreate");
			hdata[iv]->Write();
			data_hist->Write();
			data_hist->Close();

			TFile *zx_hist = new TFile("zx.root","recreate");
			hzx[iv]->Write();
			zx_hist->Write();
			zx_hist->Close();


			TFile *h_all_contributions = new TFile("./aQGC/raw_histos/all_contributions.root", "recreate");
			TH1F *h_ewk = (TH1F*) hvbs[iv]->Clone();
			h_ewk->SetName("diboson");
			h_ewk->Write();

			TH1F *h_qq = (TH1F*) hqqzz[iv]->Clone();
			h_qq->SetName("QCD_qq");
			h_qq->Write();

			TH1F *h_gg = (TH1F*) hggzz[iv]->Clone();
			h_gg->SetName("QCD_gg");
			h_gg->Write();

			TH1F *h_data = (TH1F*) hdata[iv]->Clone();
			h_data->SetName("data_obs");
			h_data->Write();

			TH1F *h_zx = (TH1F*) hzx[iv]->Clone();
			h_zx->SetName("ZpX");
			h_zx->Write();
			
			h_all_contributions->Close();
		}

		//FT8 and FT9 ADDED TO STACK
		if (iv == 1 && calculate_aQGC_limits)
		{
			hewk_FT8->Add(hsum1[iv], 1);
			hewk_FT8->SetLineColor(kYellow+1);
			hewk_FT8->SetLineWidth(3);
			hewk_FT8->SetLineStyle(9);
			//hewk_FT8->SetFillColor(kYellow+1);
			
			hewk_FT9->Add(hsum1[iv], 1);
			hewk_FT9->SetLineColor(kRed+1);
			hewk_FT9->SetLineWidth(3);
			hewk_FT9->SetLineStyle(9);
			//hewk_FT9->SetFillColor(kRed+1);
		}




	  
	  	//add histograms to stack
	  	hs[iv]->Add(hsum2[iv],"hist");
	  	hs[iv]->Add(hsum1[iv],"hist");
	  	hs[iv]->Add(hqqzz[iv],"hist");
      	hs[iv]->Add(httzwwz[iv],"hist");
	  	hs[iv]->Add(hzx[iv],"hist");
	  	TH1F *hdatadivide = (TH1F*)hdata[iv]->Clone();
	  	hs[iv]->Add(hdata[iv],"E1");
		if (iv == 1 && calculate_aQGC_limits)
		{
			hs[iv]->Add(hewk_FT8);
			hs[iv]->Add(hewk_FT9);
		}
	  
	  	// draw the legend
	  	TLegend *legend=new TLegend(0.6,0.55,0.85,0.88);
	  	legend->SetTextFont(72);
	  	legend->SetTextSize(0.04);
	  	legend->AddEntry(hzx[iv],"Z+X","f");
        legend->AddEntry(httzwwz[iv],"t#bar{t}Z, WWZ","f");
	  	legend->AddEntry(hqqzz[iv],"q#bar{q}#rightarrowZZ","f");
	  	legend->AddEntry(hsum1[iv],"gg#rightarrowZZ","f");
	  	legend->AddEntry(hsum2[iv],"VBS","f");
	  	legend->AddEntry(hdata[iv],"Data","lep");
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
		if (iv == 1) pad1->SetLogy();
	  	pad1->cd();

	  	//top plot
	  	if (iv == 1) hs[iv]->SetMaximum(100.*lumi/35.9E3);	// for m4l
		if (iv == 3) hs[iv]->SetMaximum(20.*lumi/35.9E3);  
        if (iv == 0 || iv > 3) hs[iv]->SetMaximum(45.*lumi/35.9E3);
		if (iv == 6 || iv == 7) hs[iv]->SetMaximum(30.*lumi/35.9E3); // for eta_jet1 and eta_jet2
		if (iv == 8 || iv == 9) hs[iv]->SetMaximum(50.*lumi/35.9E3); // for pt_jet1 and pt_jet2
		if (iv == 12) hs[iv]->SetMaximum(60);

	  	hs[iv]->Draw("nostack"); //old
	  	if (drawSignal[iv]) 
		{
	    	TH1F* h77 = (TH1F*)hvbs[iv]->Clone();
	    	h77->SetLineWidth(3);
	    	h77->SetLineColor(kMagenta);
	    	h77->SetFillStyle(0);
	    	h77->Scale(30);   
	    	legend->AddEntry(h77,"VBS (x30)","l");
	    	h77->Draw("histsame");
	  	}
	  	hs[iv]->GetXaxis()->SetTitle(titlex[iv].c_str());
	  	legend->Draw("same");
	  	//switch?
	  	c1->cd();
	  	//pad2
	  	TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.3+eps,0); //old position
	  	pad2->SetTopMargin(0);
	  	pad2->SetBottomMargin(0.25);
	  	//make bottom pad transparent
	  	//pad2->SetFrameFillColor(0);
	  	//pad2->SetFrameBorderMode(0);
	  	//pad2->SetFrameFillColor(0);
	  	//pad2->SetFrameBorderMode(0);
	  
	  	pad2->Draw();
	  	pad2->cd();
	  	//bottom plot
	  	TH1F *hdatacopy = (TH1F*) hdatadivide->Clone();
	  	//axis labels
	  	hdatacopy->GetXaxis()->SetLabelFont(59);//change this for font type
	  	hdatacopy->GetXaxis()->SetLabelSize(22);
	  	hdatacopy->GetYaxis()->SetLabelFont(59);//change this for font type
	  	hdatacopy->GetYaxis()->SetLabelSize(22);
	  	//axis titles
	  	hdatacopy->GetXaxis()->SetTitleFont(59); //change this for font type
	  	hdatacopy->GetXaxis()->SetTitleSize(22);
	  	hdatacopy->GetYaxis()->SetTitleFont(59); //change this for font type
	  	hdatacopy->GetYaxis()->SetTitleSize(22);
	  	hdatacopy->GetXaxis()->SetTitleOffset(4.5);
	  	hdatacopy->GetYaxis()->SetTitleOffset(1.7);
	  	hdatacopy->GetYaxis()->SetRangeUser(-1.,3.);

	  	hdatacopy->Sumw2();
	  	hdatacopy->SetStats(0); //clear stat box
	  	hdatacopy->Divide(hsum2[iv]); //invert divide
	  	hdatacopy->SetMarkerStyle(20);
	  	string title = " ;"+titlex[iv]+"; Data/MC";
	  	hdatacopy->SetTitle(title.c_str());
	  	//hdatacopy->GetXaxis()->SetTitleSize(50);
	  	//hdatacopy->GetXaxis()->SetLabelSize(35);
	  	//hdatacopy->GetYaxis()->SetTitleSize(15);
	  	//hdatacopy->GetYaxis()->SetLabelSize(15);

	  	//gStyle->SetLabelSize(2,"x").
	  	hdatacopy->Draw("ep");
	  
	  	//Orizontal line
	  	TLine *line = new TLine(xmin[iv],1,xmax[iv],1);
	  	line->SetLineColor(kRed);
	  	line->SetLineStyle(2);
	  	line->Draw("same");
	  
	  	//close and print on file
	  	c1->cd();
	  	if (useMCatNLO == 0) sprintf(filename,"onlymjjCut_%s/%s_plot_allPOWHEG_%d.png",jet_pt_cut.c_str(), namegif[iv].c_str(),year);      
	  	if (useMCatNLO == 1) sprintf(filename,"onlymjjCut_%s/%s_plot_allMCatNLO_%d.png",jet_pt_cut.c_str(), namegif[iv].c_str(),year);     
	  	if (useMCatNLO == 2) sprintf(filename,"onlymjjCut_%s/%s_plot_MCatNLOshape_POWHEGint_%d.png",jet_pt_cut.c_str(), namegif[iv].c_str(),year);
	  	gPad->Print(filename);

		cout << filename << endl;
      	c1->SaveAs(filename);
	}
	TCanvas *c2 = new TCanvas();
	c2->cd();
	hZ1Mass->SetMaximum(4500);
	hZ1Mass->SetLineColor(kRed);
	hZ1Mass->SetLineWidth(3);
	hZ1Mass->Draw();

	hZ1Mass_duje->SetLineColor(kBlue);
	hZ1Mass_duje->SetLineWidth(3);
	hZ1Mass_duje->Draw("same");
	c2->SaveAs("Z1Mass.png");

	TCanvas *c3 = new TCanvas();
	c3->cd();
	hZ2Mass->SetMaximum(2600);
	hZ2Mass->SetLineColor(kRed);
	hZ2Mass->SetLineWidth(3);
	hZ2Mass->Draw();

	hZ2Mass_duje->SetLineColor(kBlue);
	hZ2Mass_duje->SetLineWidth(3);
	hZ2Mass_duje->Draw("same");
	c3->SaveAs("Z2Mass.png");

	TCanvas *c4 = new TCanvas();
	c4->cd();
	hZMass->SetMaximum(2600);
	hZMass->SetLineColor(kRed);
	hZMass->SetLineWidth(3);
	hZMass->Draw();

	hZMass_duje->SetLineColor(kBlue);
	hZMass_duje->SetLineWidth(3);
	hZMass_duje->Draw("same");
	c4->SaveAs("ZMass.png");


	TCanvas *c5 = new TCanvas();
	c5->cd();
	hZ1M_difference->SetMaximum(2600);
	hZ1M_difference->SetLineColor(kRed);
	hZ1M_difference->SetLineWidth(3);
	hZ1M_difference->Draw();
	c5->SaveAs("Z1Mass_difference.png");

	TCanvas *c6 = new TCanvas();
	c6->cd();
	hZ2M_difference->SetMaximum(2600);
	hZ2M_difference->SetLineColor(kRed);
	hZ2M_difference->SetLineWidth(3);
	hZ2M_difference->Draw();
	c6->SaveAs("Z2Mass_difference.png");

}