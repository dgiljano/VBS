// #include "external_cConstants.h"
#include <TSpline.h>
#include <TString.h>
#include <memory>
#include "TMath.h"
#include "helper_functions.h"

string jet_pt_cut = "jet_pt_gt_30";

TH2F* rebinTemplate(TH2F* orig, int year=2016, int itype=0) {

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
   	double xbinfinal[231 ]= {160,162,164,166,168,170,172,174,176,178,180,182,184,186,188,190,192,194,196,198,200,202,204,206,208,210,212,214,216,218,220,222,224,226,228,230,232,234,236,238,240,242,244,246,248,250,252,254,256,258,260,262,264,266,268,270,272,274,276,278,280,282,284,286,288,290,292,294,296,298,
   							 300,302,304,306,308,310,312,314,316,318,320,322,324,326,328,330,332,334,336,338,340,342,344,346,348,350,352,354,356,358,360,362,364,366,368,370,372,374,376,378,380,382,384,386,388,390,392,394,396,398,400,410,420,430,440,450,460,470,480,490,
							 500,510,520,530,540,550,560,570,580,590,600,610,620,630,640,650,660,670,680,690,700,710,720,730,740,750,760,770,780,790,800,810,820,830,840,850,860,870,880,890,900,910,920,930,940,950,960,970,980,990,1000,
							 1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,1110,1120,1130,1140,1150,1160,1170,1180,1190,1200,1210,1220,1230,1240,1250,1260,1270,1280,1290,1300,1310,1320,1330,1340,1350,1360,1370,1380,1390,1400,1410,1420,1430,1440,1450,1460,1470,1480,1490,1500};
   
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
   	sprintf(filename,"template/plots/%s_%s_%d.png",orig->GetName(),pname,year);
   	gPad->Print(filename);
   	return result;    
}

void plotterAndTemplateMaker_MCStat_test(int year = 2016, int useMCatNLO = 1, int enriched = 0)
{
    //useMCatNLO = 0 : use just POWHEG
    //useMCatNLO = 1 : use just aMCatNLO
    //useMCatNLO = 2 : use aMCatNLO for shape, POWHEG for integral

    float lumi = 35.9;
    if (year == 2017) lumi = 41.5;
	if (year == 2018) lumi = 59.7;

	string theExtra = "";
    if (enriched == 1) theExtra = "_VBSenr";
	if (enriched == 2) theExtra = "_superVBSenr";
	if (enriched == 3) theExtra = "_bkgdEnr";
	if (enriched == 4) theExtra = "_ptjet50";

	static const int vars = 34;
    string titlex[vars] = {"K_{D}","M_{4l} [GeV]","M_{jj} [GeV]","#Delta #eta_{jj}","p_{T}(j)","#eta_{j}","#eta(j_{1})","#eta(j_{2})","p_{T}(j_{1})","p_{T}(j_{2})","sum(#eta_{j})","m_{jj}/#Delta#eta_{jj}","#eta*(Z_{1})","#eta*(Z_{2})","R(p_{T}^{hard})","R(p_{T}^{jet})","|#eta_{min}(j)|","|#eta_{max}(j)|","|#eta_{min}(l)|","|#eta_{max}(l)|","#Delta#Phi(Z_{1}Z_{2})","y(Z_{1})","y(Z_{2})","y(j_{1})","y(j_{2})","p_{T}(Z_{1})","p_{T}(Z_{2})","p_{T}(l_{3})","qg tagger(j_{1})","qg tagger(j_{2})","M_{4l} [GeV]","qg tagger(j_{1})","qg tagger(j_{2})","sum(|#eta(j)|)"};        
    string titley[vars] = {"Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin","Events/bin"};     
	string namegif[vars] = {"Dbkgkin","m4l","mjj","detajj","ptj","etaj","eta_j1","eta_j2","pt_jet1","pt_jet2","eta_j_sum","mjj_over_detajj","eta_Z1_star","eta_Z2_star","R_pt_hard","R_pt_jet","abs_etajet_min","abs_etajet_max","abs_etalep_min","abs_etalep_max","delta_phi_ZZ","rapidity_Z1","rapidity_Z2","rapidity_j1","rapidity_j2","pt_Z1","pt_Z2","pt_l3","j1_qg_tagger_check","j2_qg_tagger_check","m4l_original","j1_qg_tagger","j2_qg_tagger","abs_etajet_sum"};
    int bins[vars] = {20,20,20,20,30,20,20,20,30,30,20,30,20,20,30,30,30,30,30,30,70,30,30,30,30,30,30,30,50,50,20,50,50,30};
    float xmin[vars] = {0.,0.,100.,0.,0.,-5.,-5.,-5.,25.,25.,-8.,-5.,-6.,-6.,0.,0.,0.,0.,0.,0.,0.,-2.5,-2.5,-2.5,-2.5,0.,0.,0.,-1.3,-1.3,0.,0.,0.,0.};
	float xmax[vars] = {1.,1400.,1000.,8.,300.,5.,5.,5.,600.,600.,8.,400.,6.,6.,1.,1.,3.,3.,3.,3.,3.15,2.5,2.5,2.5,2.5,600.,600.,600.,1.3,1.3,1400.,1.,1.,7.};	
	bool drawSignal[vars] = {false,false,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true};

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
	TH1F *httzwzz[vars];//ttz + wwz
	TH1F *httzwzz_ee[vars];//ttz + wwz e
	TH1F *httzwzz_mm[vars];//ttz + wwz mu
	TH1F *httzwzz_em[vars];//ttz + wwz e mu

	//------------------------------------------------------------------------------- my histograms ------------------------------------------------------------------------------

	TH1F *hZ1Mass_duje = new TH1F ("Z1Mass_duje", "", 40, 50, 130);
	TH1F *hZ1Mass = new TH1F ("Z1Mass", "", 40, 50, 130);
	TH1F *hZ2Mass_duje = new TH1F ("Z2Mass_duje", "", 40, 50, 130);
	TH1F *hZ2Mass = new TH1F ("Z2Mass", "", 40, 50, 130);
	TH1F *hZMass_duje = new TH1F ("ZMass_duje", "", 40, 50, 130);
	TH1F *hZMass = new TH1F ("ZMass", "", 40, 50, 130);
	TH1F *hZ1M_difference = new TH1F ("hZ1M_difference", "", 100, -2, 2);
	TH1F *hZ2M_difference = new TH1F ("hZ2M_difference", "", 100, -2, 2);

	TH1F *httz[vars];//ttz only
	TH1F *hwwz[vars];//wwz only
	TH1F *hwzz[vars];//wzz only


	//--------------------------------------------------------------------------- end of my histograms ---------------------------------------------------------------------------

	//TF_8 and TF_9 files
	bool calculate_aQGC_limits = false;
	TString aQGC_filename = "onlymjjCut_jet_pt_gt_30/m4l_histos_" + to_string(year) + ".root";
	TFile *aQGC_histos_file = new TFile(aQGC_filename, "RECREATE");
	//TH1F *hewk_FT8, *hewk_FT9;

	//if (calculate_aQGC_limits)
	//{
		Float_t bins_FT[] = { 0, 200, 300, 400, 500, 600, 800, 1000, 1200, 1400 };

		TFile *f_ewk_FT8 = new TFile("aQGC/ewk_FT8.root");
		TH1F *hewk_FT8 = new TH1F("ewk_FT8","ewk_FT8", 9, bins_FT);
		hewk_FT8 = (TH1F*)f_ewk_FT8->Get("BLS_diboson_rescaled_FT8_1");

		TFile *f_ewk_FT9 = new TFile("aQGC/ewk_FT9.root");
		TH1F *hewk_FT9 = new TH1F("ewk_FT9","ewk_FT9", 9, bins_FT);
		hewk_FT9 = (TH1F*)f_ewk_FT9->Get("BLS_diboson_rescaled_FT9_2");

	//}
	

	for(int iv = 0; iv < vars; iv++)
	{
		if (enriched == 1 || enriched == 2) bins[iv] /= 2;
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
	  	sprintf(filename,"httzwzz_%d",iv); httzwzz[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//for 2018 rescale
	  	sprintf(filename,"httzwzz_ee_%d",iv); httzwzz_ee[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//ttzwzz e
	  	sprintf(filename,"httzwzz_mm_%d",iv); httzwzz_mm[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//ttzwzz mu
	  	sprintf(filename,"httzwzz_em_%d",iv); httzwzz_em[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//ttzwzz e mu
		sprintf(filename,"httz_%d",iv); httz[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);
		sprintf(filename,"hwwz_%d",iv); hwwz[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);
		sprintf(filename,"hwzz_%d",iv); hwzz[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);

		if (iv == 1)
		{
			Float_t bins_m4l[] = { 0, 1 };

	  		sprintf(filename,"hdata_%d",iv);   hdata[iv] = new TH1F(filename,"",1,bins_m4l); //data, because we are hiding higher energies in this phase
	  		sprintf(filename,"hqqzz_powheg_%d",iv);   hqqzz_powheg[iv] = new TH1F(filename,"",1,bins_m4l); //ew
	  		sprintf(filename,"hqqzz_%d",iv);   hqqzz[iv] = new TH1F(filename,"",1,bins_m4l); //ew+zx -> real ew
	  		sprintf(filename,"hggzz_%d",iv);   hggzz[iv] = new TH1F(filename,"",1,bins_m4l); //gg
	  		sprintf(filename,"hvbs_%d",iv);   hvbs[iv] = new TH1F(filename,"",1,bins_m4l); //vbs
	  		sprintf(filename,"hsum1_%d",iv);   hsum1[iv] = new TH1F(filename,"",1,bins_m4l); //gg+ew -> real gg
	  		sprintf(filename,"hsum2_%d",iv);   hsum2[iv] = new TH1F(filename,"",1,bins_m4l); //gg+ew+vbs
	  		sprintf(filename,"hqqzz_ee_%d",iv);   hqqzz_ee[iv] = new TH1F(filename,"",1,bins_m4l); //qqzz e
	  		sprintf(filename,"hqqzz_mm_%d",iv);   hqqzz_mm[iv] = new TH1F(filename,"",1,bins_m4l); //qqzz mu
	  		sprintf(filename,"hqqzz_em_%d",iv);   hqqzz_em[iv] = new TH1F(filename,"",1,bins_m4l); //qqzz e mu
	  		sprintf(filename,"hggzz_ee_%d",iv);   hggzz_ee[iv] = new TH1F(filename,"",1,bins_m4l);//ggzz e
	  		sprintf(filename,"hggzz_mm_%d",iv);   hggzz_mm[iv] = new TH1F(filename,"",1,bins_m4l);//ggzz mu
	  		sprintf(filename,"hggzz_em_%d",iv);  hggzz_em[iv] = new TH1F(filename,"",1,bins_m4l); //ggzz e mu
	  		sprintf(filename,"hvbs_ee_%d",iv); hvbs_ee[iv] = new TH1F(filename,"",1,bins_m4l);//vbs e
	  		sprintf(filename,"hvbs_mm_%d",iv); hvbs_mm[iv] = new TH1F(filename,"",1,bins_m4l);//vbs mu
	  		sprintf(filename,"hvbs_em_%d",iv); hvbs_em[iv] = new TH1F(filename,"",1,bins_m4l);//vbs e mu
	  		sprintf(filename,"hdataee_%d",iv); hdata_ee[iv] = new TH1F(filename,"",1,bins_m4l);//vbs e
	  		sprintf(filename,"hdatamm_%d",iv); hdata_mm[iv] = new TH1F(filename,"",1,bins_m4l);//vbs mu
	  		sprintf(filename,"hdataem_%d",iv); hdata_em[iv] = new TH1F(filename,"",1,bins_m4l);//vbs e mu
	  		sprintf(filename,"httzwzz_%d",iv); httzwzz[iv] = new TH1F(filename,"",1,bins_m4l);//for 2018 rescale
	  		sprintf(filename,"httzwzz_ee_%d",iv); httzwzz_ee[iv] = new TH1F(filename,"",1,bins_m4l);//ttzwwz e
	  		sprintf(filename,"httzwzz_mm_%d",iv); httzwzz_mm[iv] = new TH1F(filename,"",1,bins_m4l);//ttzwwz mu
	  		sprintf(filename,"httzwzz_em_%d",iv); httzwzz_em[iv] = new TH1F(filename,"",1,bins_m4l);//ttzwwz e mu
			sprintf(filename,"httz_%d",iv); httz[iv] = new TH1F(filename,"",1,bins_m4l);
			sprintf(filename,"hwwz_%d",iv); hwwz[iv] = new TH1F(filename,"",1,bins_m4l);
			sprintf(filename,"hwzz_%d",iv); hwzz[iv] = new TH1F(filename,"",1,bins_m4l);	  
		}
	}   
	gStyle->SetPalette(1);
	TFile *input_file;
    int nbins=17;
	double xbin[18]={160,166,170,176,182,188,194,200,208,216,224,234,246,260,278,302,338,1500}; 

	float c_constant = 14.0;
    // if (year == 2017) c_constant = 3.5;
	// if (year == 2018) c_constant = 3.5; 
	TFile* f_ = TFile::Open("/home/llr/cms/giljanovic/scratch/CMSSW_10_2_15/src/ZZAnalysis/AnalysisStep/data/cconstants/SmoothKDConstant_m4l_DjjVBF13TeV.root");
	TSpline3* ts = (TSpline3*)(f_->Get("sp_gr_varReco_Constant_Smooth")->Clone());
	f_->Close();

    // find available samples  
	int nSamp = 0;  
	TString rootname[40];
 	sprintf(filename,"cutbased_samples%d.txt",year);

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
	TTree *tnew[5]; 
	TFile *fnew[5];
	//template declarations (1D) 
	TH1F *temp_1d_4e[5];
	TH1F *temp_1d_4mu[5];
	TH1F *temp_1d_2e2mu[5];
	//template declarations (2D) 
	TH2F *temp_zz_4e[5];
	TH2F *temp_zz_4mu[5];
	TH2F *temp_zz_2e2mu[5];

	for (int it=0; it < 5; it++) 
	{
	  	if (it==0) sprintf(filename,"template/root_output_files/qqzz_Moriond_%d%s.root",year,theExtra.c_str()); 
	  	if (it==1) sprintf(filename,"template/root_output_files/ggzz_Moriond_%d%s.root",year,theExtra.c_str()); 
	  	if (it==2) sprintf(filename,"template/root_output_files/vbs_Moriond_%d%s.root",year,theExtra.c_str()); 
	  	if (it==3) sprintf(filename,"template/root_output_files/data_%d%s.root",year,theExtra.c_str()); 
	  	if (it==4) sprintf(filename,"template/root_output_files/ttzwzz_Moriond_%d%s.root",year,theExtra.c_str()); 
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
	float abs_etajet_min, abs_etajet_max, abs_etalep_min, abs_etalep_max;
	float delta_phi_ZZ;

	float rapidity_Z1, rapidity_Z2, rapidity_j1, rapidity_j2;
	float pt_Z1, pt_Z2, pt_l3;

	vector<float> *jet_qg_tagger = new vector<float>;

	// --------------------------------------------------------------------------------------------- end of my declarations ------------------------------------------------------------------------------------------

	// ----------------------------------------------------------------------------------------------- preparing MVA trees -------------------------------------------------------------------------------------------

	TFile *mva_file = new TFile("TMVA/MVA.root","RECREATE");

	TTree *tree_sig = new TTree("signal","signal");

	float mva_sig_mjj, mva_sig_deta_jj, mva_sig_m4l, mva_sig_eta_star_Z1, mva_sig_eta_star_Z2, mva_sig_R_pt_hard, mva_sig_R_pt_jet;

	float mva_sig_abs_etajet_min, mva_sig_abs_etajet_max, mva_sig_abs_etalep_min, mva_sig_abs_etalep_max;
	float mva_sig_delta_phi_ZZ;
	float mva_sig_rapidity_Z1, mva_sig_rapidity_Z2, mva_sig_rapidity_j1, mva_sig_rapidity_j2;
	float mva_sig_pt_Z1, mva_sig_pt_Z2, mva_sig_pt_l3;
	float mva_sig_jet1_qg_tagger, mva_sig_jet2_qg_tagger;
	float mva_sig_dbkg_kin, mva_sig_eta_j1, mva_sig_eta_j2, mva_sig_pt_jet1, mva_sig_pt_jet2, mva_sig_eta_j_sum;
	float mva_sig_mjj_over_detajj, mva_sig_abs_etajet_sum;
	float mva_sig_weight;

	tree_sig->Branch("mjj", &mva_sig_mjj, "mva_sig_mjj/F");
	tree_sig->Branch("deta_jj", &mva_sig_deta_jj, "mva_sig_deta_jj/F");
	tree_sig->Branch("m4l", &mva_sig_m4l, "mva_sig_m4l/F");
	tree_sig->Branch("eta_star_Z1", &mva_sig_eta_star_Z1, "mva_sig_eta_star_Z1/F");
	tree_sig->Branch("eta_star_Z2", &mva_sig_eta_star_Z2, "mva_sig_eta_star_Z2/F");
	tree_sig->Branch("R_pt_hard", &mva_sig_R_pt_hard, "mva_sig_R_pt_hard/F");
	tree_sig->Branch("R_pt_jet", &mva_sig_R_pt_jet, "mva_sig_R_pt_jet/F");
	
	tree_sig->Branch("abs_etajet_min", &mva_sig_abs_etajet_min, "mva_sig_abs_etajet_min/F");
	tree_sig->Branch("abs_etajet_max", &mva_sig_abs_etajet_max, "mva_sig_abs_etajet_max/F");
	tree_sig->Branch("abs_etalep_min", &mva_sig_abs_etalep_min, "mva_sig_abs_etalep_min/F");
	tree_sig->Branch("abs_etalep_max", &mva_sig_abs_etalep_max, "mva_sig_abs_etalep_max/F");
	tree_sig->Branch("delta_phi_ZZ", &mva_sig_delta_phi_ZZ, "mva_sig_delta_phi_ZZ/F");
	tree_sig->Branch("rapidity_Z1", &mva_sig_rapidity_Z1, "mva_sig_rapidity_Z1/F");
	tree_sig->Branch("rapidity_Z2", &mva_sig_rapidity_Z2, "mva_sig_rapidity_Z2/F");
	tree_sig->Branch("rapidity_j1", &mva_sig_rapidity_j1, "mva_sig_rapidity_j1/F");
	tree_sig->Branch("rapidity_j2", &mva_sig_rapidity_j2, "mva_sig_rapidity_j2/F");
	tree_sig->Branch("pt_Z1", &mva_sig_pt_Z1, "mva_sig_pt_Z1/F");
	tree_sig->Branch("pt_Z2", &mva_sig_pt_Z2, "mva_sig_pt_Z2/F");
	tree_sig->Branch("pt_l3", &mva_sig_pt_l3, "mva_sig_pt_l3/F");
	tree_sig->Branch("jet1_qg_tagger", &mva_sig_jet1_qg_tagger, "mva_sig_jet1_qg_tagger/F");
	tree_sig->Branch("jet2_qg_tagger", &mva_sig_jet2_qg_tagger, "mva_sig_jet2_qg_tagger/F");
	tree_sig->Branch("dbkg_kin", &mva_sig_dbkg_kin, "mva_sig_dbkg_kin/F");
	tree_sig->Branch("eta_j1", &mva_sig_eta_j1, "mva_sig_eta_j1/F");
	tree_sig->Branch("eta_j2", &mva_sig_eta_j2, "mva_sig_eta_j2/F");
	tree_sig->Branch("pt_jet1", &mva_sig_pt_jet1, "mva_sig_pt_jet1/F");
	tree_sig->Branch("pt_jet2", &mva_sig_pt_jet2, "mva_sig_pt_jet2/F");
	tree_sig->Branch("eta_j_sum", &mva_sig_eta_j_sum, "mva_sig_eta_j_sum/F");
	tree_sig->Branch("mjj_over_detajj", &mva_sig_mjj_over_detajj, "mva_sig_mjj_over_detajj/F");
	tree_sig->Branch("abs_etajet_sum", &mva_sig_abs_etajet_sum, "mva_sig_abs_etajet_sum/F");
	tree_sig->Branch("weight", &mva_sig_weight, "mva_sig_weight/F");

	TTree *tree_qq = new TTree("qq","qq");

	float mva_qq_mjj, mva_qq_deta_jj, mva_qq_m4l, mva_qq_eta_star_Z1, mva_qq_eta_star_Z2, mva_qq_R_pt_hard, mva_qq_R_pt_jet;
	float mva_qq_abs_etajet_min, mva_qq_abs_etajet_max, mva_qq_abs_etalep_min, mva_qq_abs_etalep_max;
	float mva_qq_delta_phi_ZZ;
	float mva_qq_rapidity_Z1, mva_qq_rapidity_Z2, mva_qq_rapidity_j1, mva_qq_rapidity_j2;
	float mva_qq_pt_Z1, mva_qq_pt_Z2, mva_qq_pt_l3;
	float mva_qq_jet1_qg_tagger, mva_qq_jet2_qg_tagger;
	float mva_qq_dbkg_kin, mva_qq_eta_j1, mva_qq_eta_j2, mva_qq_pt_jet1, mva_qq_pt_jet2, mva_qq_eta_j_sum;
	float mva_qq_mjj_over_detajj, mva_qq_abs_etajet_sum;
	float mva_qq_weight;

	tree_qq->Branch("mjj", &mva_qq_mjj, "mva_qq_mjj/F");
	tree_qq->Branch("deta_jj", &mva_qq_deta_jj, "mva_qq_deta_jj/F");
	tree_qq->Branch("m4l", &mva_qq_m4l, "mva_qq_m4l/F");
	tree_qq->Branch("eta_star_Z1", &mva_qq_eta_star_Z1, "mva_qq_eta_star_Z1/F");
	tree_qq->Branch("eta_star_Z2", &mva_qq_eta_star_Z2, "mva_qq_eta_star_Z2/F");
	tree_qq->Branch("R_pt_hard", &mva_qq_R_pt_hard, "mva_qq_R_pt_hard/F");
	tree_qq->Branch("R_pt_jet", &mva_qq_R_pt_jet, "mva_qq_R_pt_jet/F");

	tree_qq->Branch("abs_etajet_min", &mva_qq_abs_etajet_min, "mva_qq_abs_etajet_min/F");
	tree_qq->Branch("abs_etajet_max", &mva_qq_abs_etajet_max, "mva_qq_abs_etajet_max/F");
	tree_qq->Branch("abs_etalep_min", &mva_qq_abs_etalep_min, "mva_qq_abs_etalep_min/F");
	tree_qq->Branch("abs_etalep_max", &mva_qq_abs_etalep_max, "mva_qq_abs_etalep_max/F");
	tree_qq->Branch("delta_phi_ZZ", &mva_qq_delta_phi_ZZ, "mva_qq_delta_phi_ZZ/F");
	tree_qq->Branch("rapidity_Z1", &mva_qq_rapidity_Z1, "mva_qq_rapidity_Z1/F");
	tree_qq->Branch("rapidity_Z2", &mva_qq_rapidity_Z2, "mva_qq_rapidity_Z2/F");
	tree_qq->Branch("rapidity_j1", &mva_qq_rapidity_j1, "mva_qq_rapidity_j1/F");
	tree_qq->Branch("rapidity_j2", &mva_qq_rapidity_j2, "mva_qq_rapidity_j2/F");
	tree_qq->Branch("pt_Z1", &mva_qq_pt_Z1, "mva_qq_pt_Z1/F");
	tree_qq->Branch("pt_Z2", &mva_qq_pt_Z2, "mva_qq_pt_Z2/F");
	tree_qq->Branch("pt_l3", &mva_qq_pt_l3, "mva_qq_pt_l3/F");
	tree_qq->Branch("jet1_qg_tagger", &mva_qq_jet1_qg_tagger, "mva_qq_jet1_qg_tagger/F");
	tree_qq->Branch("jet2_qg_tagger", &mva_qq_jet2_qg_tagger, "mva_qq_jet2_qg_tagger/F");
	tree_qq->Branch("dbkg_kin", &mva_qq_dbkg_kin, "mva_qq_dbkg_kin/F");
	tree_qq->Branch("eta_j1", &mva_qq_eta_j1, "mva_qq_eta_j1/F");
	tree_qq->Branch("eta_j2", &mva_qq_eta_j2, "mva_qq_eta_j2/F");
	tree_qq->Branch("pt_jet1", &mva_qq_pt_jet1, "mva_qq_pt_jet1/F");
	tree_qq->Branch("pt_jet2", &mva_qq_pt_jet2, "mva_qq_pt_jet2/F");
	tree_qq->Branch("eta_j_sum", &mva_qq_eta_j_sum, "mva_qq_eta_j_sum/F");
	tree_qq->Branch("mjj_over_detajj", &mva_qq_mjj_over_detajj, "mva_qq_mjj_over_detajj/F");
	tree_qq->Branch("abs_etajet_sum", &mva_qq_abs_etajet_sum, "mva_qq_abs_etajet_sum/F");
	tree_qq->Branch("weight", &mva_qq_weight, "mva_qq_weight/F");


	TTree *tree_gg = new TTree("gg","gg");

	float mva_gg_mjj, mva_gg_deta_jj, mva_gg_m4l, mva_gg_eta_star_Z1, mva_gg_eta_star_Z2, mva_gg_R_pt_hard, mva_gg_R_pt_jet;
	float mva_gg_abs_etajet_min, mva_gg_abs_etajet_max, mva_gg_abs_etalep_min, mva_gg_abs_etalep_max;
	float mva_gg_delta_phi_ZZ;
	float mva_gg_rapidity_Z1, mva_gg_rapidity_Z2, mva_gg_rapidity_j1, mva_gg_rapidity_j2;
	float mva_gg_pt_Z1, mva_gg_pt_Z2, mva_gg_pt_l3;
	float mva_gg_jet1_qg_tagger, mva_gg_jet2_qg_tagger;
	float mva_gg_dbkg_kin, mva_gg_eta_j1, mva_gg_eta_j2, mva_gg_pt_jet1, mva_gg_pt_jet2, mva_gg_eta_j_sum;
	float mva_gg_mjj_over_detajj, mva_gg_abs_etajet_sum;
	float mva_gg_weight;

	tree_gg->Branch("mjj", &mva_gg_mjj, "mva_gg_mjj/F");
	tree_gg->Branch("deta_jj", &mva_gg_deta_jj, "mva_gg_deta_jj/F");
	tree_gg->Branch("m4l", &mva_gg_m4l, "mva_gg_m4l/F");
	tree_gg->Branch("eta_star_Z1", &mva_gg_eta_star_Z1, "mva_gg_eta_star_Z1/F");
	tree_gg->Branch("eta_star_Z2", &mva_gg_eta_star_Z2, "mva_gg_eta_star_Z2/F");
	tree_gg->Branch("R_pt_hard", &mva_gg_R_pt_hard, "mva_gg_R_pt_hard/F");
	tree_gg->Branch("R_pt_jet", &mva_gg_R_pt_jet, "mva_gg_R_pt_jet/F");

	tree_gg->Branch("abs_etajet_min", &mva_gg_abs_etajet_min, "mva_gg_abs_etajet_min/F");
	tree_gg->Branch("abs_etajet_max", &mva_gg_abs_etajet_max, "mva_gg_abs_etajet_max/F");
	tree_gg->Branch("abs_etalep_min", &mva_gg_abs_etalep_min, "mva_gg_abs_etalep_min/F");
	tree_gg->Branch("abs_etalep_max", &mva_gg_abs_etalep_max, "mva_gg_abs_etalep_max/F");
	tree_gg->Branch("delta_phi_ZZ", &mva_gg_delta_phi_ZZ, "mva_gg_delta_phi_ZZ/F");
	tree_gg->Branch("rapidity_Z1", &mva_gg_rapidity_Z1, "mva_gg_rapidity_Z1/F");
	tree_gg->Branch("rapidity_Z2", &mva_gg_rapidity_Z2, "mva_gg_rapidity_Z2/F");
	tree_gg->Branch("rapidity_j1", &mva_gg_rapidity_j1, "mva_gg_rapidity_j1/F");
	tree_gg->Branch("rapidity_j2", &mva_gg_rapidity_j2, "mva_gg_rapidity_j2/F");
	tree_gg->Branch("pt_Z1", &mva_gg_pt_Z1, "mva_gg_pt_Z1/F");
	tree_gg->Branch("pt_Z2", &mva_gg_pt_Z2, "mva_gg_pt_Z2/F");
	tree_gg->Branch("pt_l3", &mva_gg_pt_l3, "mva_gg_pt_l3/F");
	tree_gg->Branch("jet1_qg_tagger", &mva_gg_jet1_qg_tagger, "mva_gg_jet1_qg_tagger/F");
	tree_gg->Branch("jet2_qg_tagger", &mva_gg_jet2_qg_tagger, "mva_gg_jet2_qg_tagger/F");
	tree_gg->Branch("dbkg_kin", &mva_gg_dbkg_kin, "mva_gg_dbkg_kin/F");
	tree_gg->Branch("eta_j1", &mva_gg_eta_j1, "mva_gg_eta_j1/F");
	tree_gg->Branch("eta_j2", &mva_gg_eta_j2, "mva_gg_eta_j2/F");
	tree_gg->Branch("pt_jet1", &mva_gg_pt_jet1, "mva_gg_pt_jet1/F");
	tree_gg->Branch("pt_jet2", &mva_gg_pt_jet2, "mva_gg_pt_jet2/F");
	tree_gg->Branch("eta_j_sum", &mva_gg_eta_j_sum, "mva_gg_eta_j_sum/F");
	tree_gg->Branch("mjj_over_detajj", &mva_gg_mjj_over_detajj, "mva_gg_mjj_over_detajj/F");
	tree_gg->Branch("abs_etajet_sum", &mva_gg_abs_etajet_sum, "mva_gg_abs_etajet_sum/F");
	tree_gg->Branch("weight", &mva_gg_weight, "mva_gg_weight/F");


	TTree *tree_zx = new TTree("zx","zx");

	float mva_zx_mjj, mva_zx_deta_jj, mva_zx_m4l, mva_zx_eta_star_Z1, mva_zx_eta_star_Z2, mva_zx_R_pt_hard, mva_zx_R_pt_jet;

	float mva_zx_abs_etajet_min, mva_zx_abs_etajet_max, mva_zx_abs_etalep_min, mva_zx_abs_etalep_max;
	float mva_zx_delta_phi_ZZ;
	float mva_zx_rapidity_Z1, mva_zx_rapidity_Z2, mva_zx_rapidity_j1, mva_zx_rapidity_j2;
	float mva_zx_pt_Z1, mva_zx_pt_Z2, mva_zx_pt_l3;
	float mva_zx_jet1_qg_tagger, mva_zx_jet2_qg_tagger;
	float mva_zx_dbkg_kin, mva_zx_eta_j1, mva_zx_eta_j2, mva_zx_pt_jet1, mva_zx_pt_jet2, mva_zx_eta_j_sum;
	float mva_zx_mjj_over_detajj, mva_zx_abs_etajet_sum;
	float mva_zx_weight;

	tree_zx->Branch("mjj", &mva_zx_mjj, "mva_zx_mjj/F");
	tree_zx->Branch("deta_jj", &mva_zx_deta_jj, "mva_zx_deta_jj/F");
	tree_zx->Branch("m4l", &mva_zx_m4l, "mva_zx_m4l/F");
	tree_zx->Branch("eta_star_Z1", &mva_zx_eta_star_Z1, "mva_zx_eta_star_Z1/F");
	tree_zx->Branch("eta_star_Z2", &mva_zx_eta_star_Z2, "mva_zx_eta_star_Z2/F");
	tree_zx->Branch("R_pt_hard", &mva_zx_R_pt_hard, "mva_zx_R_pt_hard/F");
	tree_zx->Branch("R_pt_jet", &mva_zx_R_pt_jet, "mva_zx_R_pt_jet/F");

	tree_zx->Branch("abs_etajet_min", &mva_zx_abs_etajet_min, "mva_zx_abs_etajet_min/F");
	tree_zx->Branch("abs_etajet_max", &mva_zx_abs_etajet_max, "mva_zx_abs_etajet_max/F");
	tree_zx->Branch("abs_etalep_min", &mva_zx_abs_etalep_min, "mva_zx_abs_etalep_min/F");
	tree_zx->Branch("abs_etalep_max", &mva_zx_abs_etalep_max, "mva_zx_abs_etalep_max/F");
	tree_zx->Branch("delta_phi_ZZ", &mva_zx_delta_phi_ZZ, "mva_zx_delta_phi_ZZ/F");
	tree_zx->Branch("rapidity_Z1", &mva_zx_rapidity_Z1, "mva_zx_rapidity_Z1/F");
	tree_zx->Branch("rapidity_Z2", &mva_zx_rapidity_Z2, "mva_zx_rapidity_Z2/F");
	tree_zx->Branch("rapidity_j1", &mva_zx_rapidity_j1, "mva_zx_rapidity_j1/F");
	tree_zx->Branch("rapidity_j2", &mva_zx_rapidity_j2, "mva_zx_rapidity_j2/F");
	tree_zx->Branch("pt_Z1", &mva_zx_pt_Z1, "mva_zx_pt_Z1/F");
	tree_zx->Branch("pt_Z2", &mva_zx_pt_Z2, "mva_zx_pt_Z2/F");
	tree_zx->Branch("pt_l3", &mva_zx_pt_l3, "mva_zx_pt_l3/F");
	tree_zx->Branch("jet1_qg_tagger", &mva_zx_jet1_qg_tagger, "mva_zx_jet1_qg_tagger/F");
	tree_zx->Branch("jet2_qg_tagger", &mva_zx_jet2_qg_tagger, "mva_zx_jet2_qg_tagger/F");
	tree_zx->Branch("dbkg_kin", &mva_zx_dbkg_kin, "mva_zx_dbkg_kin/F");
	tree_zx->Branch("eta_j1", &mva_zx_eta_j1, "mva_zx_eta_j1/F");
	tree_zx->Branch("eta_j2", &mva_zx_eta_j2, "mva_zx_eta_j2/F");
	tree_zx->Branch("pt_jet1", &mva_zx_pt_jet1, "mva_zx_pt_jet1/F");
	tree_zx->Branch("pt_jet2", &mva_zx_pt_jet2, "mva_zx_pt_jet2/F");
	tree_zx->Branch("eta_j_sum", &mva_zx_eta_j_sum, "mva_zx_eta_j_sum/F");
	tree_zx->Branch("mjj_over_detajj", &mva_zx_mjj_over_detajj, "mva_zx_mjj_over_detajj/F");
	tree_zx->Branch("abs_etajet_sum", &mva_zx_abs_etajet_sum, "mva_zx_abs_etajet_sum/F");
	tree_zx->Branch("weight", &mva_zx_weight, "mva_zx_weight/F");


	TTree *tree_ttzwwz = new TTree("ttzwwz","ttzwwz");

	float mva_ttzwwz_mjj, mva_ttzwwz_deta_jj, mva_ttzwwz_m4l, mva_ttzwwz_eta_star_Z1, mva_ttzwwz_eta_star_Z2, mva_ttzwwz_R_pt_hard, mva_ttzwwz_R_pt_jet;
	float mva_ttzwwz_abs_etajet_min, mva_ttzwwz_abs_etajet_max, mva_ttzwwz_abs_etalep_min, mva_ttzwwz_abs_etalep_max;
	float mva_ttzwwz_delta_phi_ZZ;
	float mva_ttzwwz_rapidity_Z1, mva_ttzwwz_rapidity_Z2, mva_ttzwwz_rapidity_j1, mva_ttzwwz_rapidity_j2;
	float mva_ttzwwz_pt_Z1, mva_ttzwwz_pt_Z2, mva_ttzwwz_pt_l3;
	float mva_ttzwwz_jet1_qg_tagger, mva_ttzwwz_jet2_qg_tagger;
	float mva_ttzwwz_dbkg_kin, mva_ttzwwz_eta_j1, mva_ttzwwz_eta_j2, mva_ttzwwz_pt_jet1, mva_ttzwwz_pt_jet2, mva_ttzwwz_eta_j_sum;
	float mva_ttzwwz_mjj_over_detajj, mva_ttzwwz_abs_etajet_sum;
	float mva_ttzwwz_weight;

	tree_ttzwwz->Branch("mjj", &mva_ttzwwz_mjj, "mva_ttzwwz_mjj/F");
	tree_ttzwwz->Branch("deta_jj", &mva_ttzwwz_deta_jj, "mva_ttzwwz_deta_jj/F");
	tree_ttzwwz->Branch("m4l", &mva_ttzwwz_m4l, "mva_ttzwwz_m4l/F");
	tree_ttzwwz->Branch("eta_star_Z1", &mva_ttzwwz_eta_star_Z1, "mva_ttzwwz_eta_star_Z1/F");
	tree_ttzwwz->Branch("eta_star_Z2", &mva_ttzwwz_eta_star_Z2, "mva_ttzwwz_eta_star_Z2/F");
	tree_ttzwwz->Branch("R_pt_hard", &mva_ttzwwz_R_pt_hard, "mva_ttzwwz_R_pt_hard/F");
	tree_ttzwwz->Branch("R_pt_jet", &mva_ttzwwz_R_pt_jet, "mva_ttzwwz_R_pt_jet/F");

	tree_ttzwwz->Branch("abs_etajet_min", &mva_ttzwwz_abs_etajet_min, "mva_ttzwwz_abs_etajet_min/F");
	tree_ttzwwz->Branch("abs_etajet_max", &mva_ttzwwz_abs_etajet_max, "mva_ttzwwz_abs_etajet_max/F");
	tree_ttzwwz->Branch("abs_etalep_min", &mva_ttzwwz_abs_etalep_min, "mva_ttzwwz_abs_etalep_min/F");
	tree_ttzwwz->Branch("abs_etalep_max", &mva_ttzwwz_abs_etalep_max, "mva_ttzwwz_abs_etalep_max/F");
	tree_ttzwwz->Branch("delta_phi_ZZ", &mva_ttzwwz_delta_phi_ZZ, "mva_ttzwwz_delta_phi_ZZ/F");
	tree_ttzwwz->Branch("rapidity_Z1", &mva_ttzwwz_rapidity_Z1, "mva_ttzwwz_rapidity_Z1/F");
	tree_ttzwwz->Branch("rapidity_Z2", &mva_ttzwwz_rapidity_Z2, "mva_ttzwwz_rapidity_Z2/F");
	tree_ttzwwz->Branch("rapidity_j1", &mva_ttzwwz_rapidity_j1, "mva_ttzwwz_rapidity_j1/F");
	tree_ttzwwz->Branch("rapidity_j2", &mva_ttzwwz_rapidity_j2, "mva_ttzwwz_rapidity_j2/F");
	tree_ttzwwz->Branch("pt_Z1", &mva_ttzwwz_pt_Z1, "mva_ttzwwz_pt_Z1/F");
	tree_ttzwwz->Branch("pt_Z2", &mva_ttzwwz_pt_Z2, "mva_ttzwwz_pt_Z2/F");
	tree_ttzwwz->Branch("pt_l3", &mva_ttzwwz_pt_l3, "mva_ttzwwz_pt_l3/F");
	tree_ttzwwz->Branch("jet1_qg_tagger", &mva_ttzwwz_jet1_qg_tagger, "mva_ttzwwz_jet1_qg_tagger/F");
	tree_ttzwwz->Branch("jet2_qg_tagger", &mva_ttzwwz_jet2_qg_tagger, "mva_ttzwwz_jet2_qg_tagger/F");
	tree_ttzwwz->Branch("dbkg_kin", &mva_ttzwwz_dbkg_kin, "mva_ttzwwz_dbkg_kin/F");
	tree_ttzwwz->Branch("eta_j1", &mva_ttzwwz_eta_j1, "mva_ttzwwz_eta_j1/F");
	tree_ttzwwz->Branch("eta_j2", &mva_ttzwwz_eta_j2, "mva_ttzwwz_eta_j2/F");
	tree_ttzwwz->Branch("pt_jet1", &mva_ttzwwz_pt_jet1, "mva_ttzwwz_pt_jet1/F");
	tree_ttzwwz->Branch("pt_jet2", &mva_ttzwwz_pt_jet2, "mva_ttzwwz_pt_jet2/F");
	tree_ttzwwz->Branch("eta_j_sum", &mva_ttzwwz_eta_j_sum, "mva_ttzwwz_eta_j_sum/F");
	tree_ttzwwz->Branch("mjj_over_detajj", &mva_ttzwwz_mjj_over_detajj, "mva_ttzwwz_mjj_over_detajj/F");
	tree_ttzwwz->Branch("abs_etajet_sum", &mva_ttzwwz_abs_etajet_sum, "mva_ttzwwz_abs_etajet_sum/F");
	tree_ttzwwz->Branch("weight", &mva_ttzwwz_weight, "mva_ttzwwz_weight/F");


	TTree *tree_data = new TTree("data","data");

	float mva_data_mjj, mva_data_deta_jj, mva_data_m4l, mva_data_eta_star_Z1, mva_data_eta_star_Z2, mva_data_R_pt_hard, mva_data_R_pt_jet;
	float mva_data_abs_etajet_min, mva_data_abs_etajet_max, mva_data_abs_etalep_min, mva_data_abs_etalep_max;
	float mva_data_delta_phi_ZZ;
	float mva_data_rapidity_Z1, mva_data_rapidity_Z2, mva_data_rapidity_j1, mva_data_rapidity_j2;
	float mva_data_pt_Z1, mva_data_pt_Z2, mva_data_pt_l3;
	float mva_data_jet1_qg_tagger, mva_data_jet2_qg_tagger;
	float mva_data_dbkg_kin, mva_data_eta_j1, mva_data_eta_j2, mva_data_pt_jet1, mva_data_pt_jet2, mva_data_eta_j_sum;
	float mva_data_mjj_over_detajj, mva_data_abs_etajet_sum;
	float mva_data_weight;

	tree_data->Branch("mjj", &mva_data_mjj, "mva_data_mjj/F");
	tree_data->Branch("deta_jj", &mva_data_deta_jj, "mva_data_deta_jj/F");
	tree_data->Branch("m4l", &mva_data_m4l, "mva_data_m4l/F");
	tree_data->Branch("eta_star_Z1", &mva_data_eta_star_Z1, "mva_data_eta_star_Z1/F");
	tree_data->Branch("eta_star_Z2", &mva_data_eta_star_Z2, "mva_data_eta_star_Z2/F");
	tree_data->Branch("R_pt_hard", &mva_data_R_pt_hard, "mva_data_R_pt_hard/F");
	tree_data->Branch("R_pt_jet", &mva_data_R_pt_jet, "mva_data_R_pt_jet/F");

	tree_data->Branch("abs_etajet_min", &mva_data_abs_etajet_min, "mva_data_abs_etajet_min/F");
	tree_data->Branch("abs_etajet_max", &mva_data_abs_etajet_max, "mva_data_abs_etajet_max/F");
	tree_data->Branch("abs_etalep_min", &mva_data_abs_etalep_min, "mva_data_abs_etalep_min/F");
	tree_data->Branch("abs_etalep_max", &mva_data_abs_etalep_max, "mva_data_abs_etalep_max/F");
	tree_data->Branch("delta_phi_ZZ", &mva_data_delta_phi_ZZ, "mva_data_delta_phi_ZZ/F");
	tree_data->Branch("rapidity_Z1", &mva_data_rapidity_Z1, "mva_data_rapidity_Z1/F");
	tree_data->Branch("rapidity_Z2", &mva_data_rapidity_Z2, "mva_data_rapidity_Z2/F");
	tree_data->Branch("rapidity_j1", &mva_data_rapidity_j1, "mva_data_rapidity_j1/F");
	tree_data->Branch("rapidity_j2", &mva_data_rapidity_j2, "mva_data_rapidity_j2/F");
	tree_data->Branch("pt_Z1", &mva_data_pt_Z1, "mva_data_pt_Z1/F");
	tree_data->Branch("pt_Z2", &mva_data_pt_Z2, "mva_data_pt_Z2/F");
	tree_data->Branch("pt_l3", &mva_data_pt_l3, "mva_data_pt_l3/F");
	tree_data->Branch("jet1_qg_tagger", &mva_data_jet1_qg_tagger, "mva_data_jet1_qg_tagger/F");
	tree_data->Branch("jet2_qg_tagger", &mva_data_jet2_qg_tagger, "mva_data_jet2_qg_tagger/F");
	tree_data->Branch("dbkg_kin", &mva_data_dbkg_kin, "mva_data_dbkg_kin/F");
	tree_data->Branch("eta_j1", &mva_data_eta_j1, "mva_data_eta_j1/F");
	tree_data->Branch("eta_j2", &mva_data_eta_j2, "mva_data_eta_j2/F");
	tree_data->Branch("pt_jet1", &mva_data_pt_jet1, "mva_data_pt_jet1/F");
	tree_data->Branch("pt_jet2", &mva_data_pt_jet2, "mva_data_pt_jet2/F");
	tree_data->Branch("eta_j_sum", &mva_data_eta_j_sum, "mva_data_eta_j_sum/F");
	tree_data->Branch("mjj_over_detajj", &mva_data_mjj_over_detajj, "mva_data_mjj_over_detajj/F");
	tree_data->Branch("abs_etajet_sum", &mva_data_abs_etajet_sum, "mva_data_abs_etajet_sum/F");
	tree_data->Branch("weight", &mva_data_weight, "mva_data_weight/F");


	// ------------------------------------------------------------------------------------------- end of preparing MVA trees ----------------------------------------------------------------------------------------
	
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
		bool is_ttZ = false;
		bool is_wwZ = false;
		bool is_wZZ = false;
		if (rootname[is].Contains("AllData")) j = 3;   
		if (rootname[is].Contains("ggTo") || rootname[is].Contains("ggZZnew")) j = 1;      
    	if (rootname[is].Contains("VBFTo")) j = 2;
		if (rootname[is].Contains("amcatnlo")) j = 5;
		if (rootname[is].Contains("WWZ") || rootname[is].Contains("TTZ") || rootname[is].Contains("WZZ")) j = 4;   	  
		if (rootname[is].Contains("WWZ")) is_wwZ = true;
		if (rootname[is].Contains("TTZ")) is_ttZ = true;
		if (rootname[is].Contains("WZZ")) is_wZZ = true;
	
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
		tqqzz->SetBranchAddress("JetQGLikelihood",&jet_qg_tagger);

		// ----------------------------------------------------------------------------- end of my branches ------------------------------------------------------------

		//loop on entries
		int enne = tqqzz->GetEntries();

		// preliminary loop to fix MG wrong weights (only VBS 2017-18)
		float resum = gen_sum_weights;
    	if (j==2 && year>2016) resum = hCounters->GetBinContent(1);
	
		for(int i=0;i<enne;i++)
		{
	    	tqqzz->GetEntry(i);

	    	if(DiJetMass > 100 && ZZMass > 180 && nCleanedJetsPt30 > 1 && Z1Mass < 120 && Z1Mass > 60 && Z2Mass < 120 && Z2Mass > 60)
			{
				if (enriched == 1 && (DiJetMass < 400 || fabs(DiJetDEta) < 2.4)) continue;
	      		if (enriched == 2 && (DiJetMass < 400 || fabs(DiJetDEta) < 5.0)) continue;
	      		if (enriched == 3 && DiJetMass > 400 && fabs(DiJetDEta) > 2.4) continue;
	      		if (enriched == 4 && (JetPt->at(0) < 50 || JetPt->at(1) < 50)) continue;
				// ------------------------------------------------------------ construc electron and muon objects -----------------------------------------
				
				clear_vectors(electrons, muons, electrons_charge, muons_charge);
				create_electron_and_muon_objects(electrons, muons, electrons_charge, muons_charge, lepId, lepPt, lepEta, lepPhi);
				build_ZZ_pair(electrons, muons, electrons_charge, muons_charge, Z1, Z2);

				if (is == 3)
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
				calculate_min_max_jet_eta(JetEta, abs_etajet_min, abs_etajet_max);
				calculate_min_max_lepton_eta(electrons, muons, abs_etalep_min, abs_etalep_max);
				calculate_dphi_ZZ(Z1, Z2, delta_phi_ZZ);
				
				calculate_rapidity_Z1_Z2(Z1, Z2, rapidity_Z1, rapidity_Z2);
				calculate_rapidity_j1_j2(JetEta, JetPhi, JetPt, rapidity_j1, rapidity_j2);
				calculate_pt_Z1_Z2_l3(Z1, Z2, lepPt, pt_Z1, pt_Z2, pt_l3);


				// -------------------------------------------------------- end of calculate aditional variables -------------------------------------------
				  
	      		//set vbf_category
	      		vbfcate=1;
				  
				// make sure prefiring weight is 1 for real data
	      		float prefiringWeight = L1prefiringWeight;
	      		if (j==2) prefiringWeight = 1.;
	      		
				weight= (xsec*KFactorEWKqqZZ*overallEventWeight*KFactorQCDqqZZ_M*prefiringWeight*lumi)/(resum);
				if (j==1) weight= (xsec*overallEventWeight*1.53*1.64*prefiringWeight*lumi)/(resum);
              	if (j==2 && year==2016) weight= (xsec*overallEventWeight*prefiringWeight*lumi)/(resum);
              	if (j==2 && year>2016) weight= (xsec*overallEventWeight*prefiringWeight*lumi)/(genHEPMCweight*resum);             
	      		if (j==5) weight= (xsec*overallEventWeight*prefiringWeight*lumi)/(resum);
	      		if (j==3) weight=prefiringWeight; 

	      		//division in channels
	      		if(abs(Z1Flav)==abs(Z2Flav) && abs(Z1Flav)==121) chan=2;
	      		else if (abs(Z1Flav)==abs(Z2Flav) && abs(Z1Flav)!=121) chan=1;
	      		else chan=3;

	      		//kin variable
            	float c_mzz = c_constant*ts->Eval(ZZMass);
	      		dbkg_kin = p_JJVBF_BKG_MCFM_JECNominal/(p_JJVBF_BKG_MCFM_JECNominal+ p_JJQCD_BKG_MCFM_JECNominal*c_mzz);
	      		if (dbkg_kin < 0.00 || dbkg_kin > 1.00) continue;

	      		// fill templates
	      		if (useMCatNLO > 0 && j==5) 
				{
					tnew[0]->Fill();
					temp_zz_4mu[0]->Fill(ZZMass,dbkg_kin,weight); temp_1d_4mu[0]->Fill(dbkg_kin,weight);  
					temp_zz_4e[0]->Fill(ZZMass,dbkg_kin,weight); temp_1d_4e[0]->Fill(dbkg_kin,weight); 
					temp_zz_2e2mu[0]->Fill(ZZMass,dbkg_kin,weight); temp_1d_2e2mu[0]->Fill(dbkg_kin,weight); 
	      		} 
				else if (useMCatNLO == 0 && j==0) 
				{
					tnew[0]->Fill();
					temp_zz_4mu[0]->Fill(ZZMass,dbkg_kin,weight); temp_1d_4mu[0]->Fill(dbkg_kin,weight); 
					temp_zz_4e[0]->Fill(ZZMass,dbkg_kin,weight);  temp_1d_4e[0]->Fill(dbkg_kin,weight); 
					temp_zz_2e2mu[0]->Fill(ZZMass,dbkg_kin,weight); temp_1d_2e2mu[0]->Fill(dbkg_kin,weight); 
	      		} 
				else if (j==1 || j==4) 
				{
					tnew[j]->Fill();
					temp_zz_4mu[j]->Fill(ZZMass,dbkg_kin,weight);  temp_1d_4mu[j]->Fill(dbkg_kin,weight);
					temp_zz_4e[j]->Fill(ZZMass,dbkg_kin,weight);  temp_1d_4e[j]->Fill(dbkg_kin,weight);
	        		temp_zz_2e2mu[j]->Fill(ZZMass,dbkg_kin,weight);  temp_1d_2e2mu[j]->Fill(dbkg_kin,weight);
	      		} 
				else 
				{
					tnew[j]->Fill();
					if(chan==1) 
					{ 
						temp_zz_4mu[j]->Fill(ZZMass,dbkg_kin,weight);
						temp_1d_4mu[j]->Fill(dbkg_kin,weight);
					}
					else if(chan==2)
					{
						temp_zz_4e[j]->Fill(ZZMass,dbkg_kin,weight);
						temp_1d_4e[j]->Fill(dbkg_kin,weight);
					} 
	    			else
					{
						temp_zz_2e2mu[j]->Fill(ZZMass,dbkg_kin,weight);
						temp_1d_2e2mu[j]->Fill(dbkg_kin,weight);
					}
	      		}

            	// fill plots
	      		for(int il = 0; il < vars+2; il++)
				{
            	    int iv = il; 
					if (il == 0) theVar = dbkg_kin;
					if (il == 1) theVar = ZZMass;
					//if (il == 1) theVar = 0.5;
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
						theVar = DiJetMass/fabs(DiJetDEta);
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
					if (il == 18)
					{
						theVar = abs_etajet_min;
						iv = 16;
					}
					if (il == 19)
					{
						theVar = abs_etajet_max;
						iv = 17;
					}
					if (il == 20)
					{
						theVar = abs_etalep_min;
						iv = 18;
					}
					if (il == 21)
					{
						theVar = abs_etalep_max;
						iv = 19;
					}
					if (il == 22)
					{
						theVar = delta_phi_ZZ;
						iv = 20;
					}
					if (il == 23)
					{
						theVar = rapidity_Z1;
						iv = 21;
					}
					if (il == 24)
					{
						theVar = rapidity_Z2;
						iv = 22;
					}
					if (il == 25)
					{
						theVar = rapidity_j1;
						iv = 23;
					}
					if (il == 26)
					{
						theVar = rapidity_j2;
						iv = 24;
					}
					if (il == 27)
					{
						theVar = pt_Z1;
						iv = 25;
					}
					if (il == 28)
					{
						theVar = pt_Z2;
						iv = 26;
					}
					if (il == 29)
					{
						theVar = pt_l3;
						iv = 27;
					}
					if (il == 30)
					{
						theVar = jet_qg_tagger->at(0);
						iv = 28;
					}
					if (il == 31)
					{
						theVar = jet_qg_tagger->at(1);
						iv = 29;
					}
					if (il == 32)
					{
						theVar = ZZMass;
						iv = 30;
					}
					if (il == 33)
					{
						theVar = jet_qg_tagger->at(0);
						iv = 31;
					}
					if (il == 34)
					{
						theVar = jet_qg_tagger->at(1);
						iv = 32;
					}
					if (il == 35)
					{
						theVar = fabs(JetEta->at(0)) + fabs(JetEta->at(1));
						iv = 33;
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
					if (j==0)	//qqzz
					{ 
					  	if (il == 1 && theVar >= 1200 && theVar < 1400) hqqzz_powheg[iv]->Fill(0.5,weight);

					  /* if (chan == 1) hqqzz_mm[iv]->Fill(theVar,weight); //mu
					  if (chan == 2) hqqzz_ee[iv]->Fill(theVar,weight); //e
					  if (chan == 3)  hqqzz_em[iv]->Fill(theVar,weight);  */        //mu+e
					}
					if (j==1)	//gg
					{ 	                             
					  	if (il == 1 && theVar >= 1200 && theVar < 1400) hggzz[iv]->Fill(0.5,weight);
					  	hggzz[iv]->SetFillColor(kBlue);

					  	if (chan == 1) hggzz_mm[iv]->Fill(theVar,weight); //mu
					  	if (chan == 2) hggzz_ee[iv]->Fill(theVar,weight); //e
					  	if (chan == 3)  hggzz_em[iv]->Fill(theVar,weight);          //mu+e

					}
					if (j==2)	//vbs
					{                             
					  	if (il == 1 && theVar >= 1200 && theVar < 1400) hvbs[iv]->Fill(0.5,weight);
					  	hvbs[iv]->SetFillColor(kMagenta);

					  	if (chan == 1) hvbs_mm[iv]->Fill(theVar,weight); //mu
					  	if (chan == 2) hvbs_ee[iv]->Fill(theVar,weight); //e
					  	if (chan == 3) hvbs_em[iv]->Fill(theVar,weight);          //mu+e

					}
					if (j==5) 
					{
					  	if (il == 1 && theVar >= 1200 && theVar < 1400) hqqzz[iv]->Fill(0.5,weight);

					  	if (chan == 1) hqqzz_mm[iv]->Fill(theVar,weight); //mu
					  	if (chan == 2) hqqzz_ee[iv]->Fill(theVar,weight); //e
					  	if (chan == 3) hqqzz_em[iv]->Fill(theVar,weight);
					}
					if (j==4) 
					{
					  	if (il == 1 && theVar >= 1200 && theVar < 1400) httzwzz[iv]->Fill(0.5,weight);

					  	if (chan == 1) httzwzz_mm[iv]->Fill(theVar,weight); //mu
					  	if (chan == 2) httzwzz_ee[iv]->Fill(theVar,weight); //e
					  	if (chan == 3) httzwzz_em[iv]->Fill(theVar,weight);          //mu+e
					}
					if (is_ttZ) httz[iv]->Fill(theVar,weight);
					if (is_wwZ) hwwz[iv]->Fill(theVar,weight);
					if (is_wZZ) hwzz[iv]->Fill(theVar,weight);
	      		}
                // filling trees for TMVA

				if (j==2)	//signal
				{
					mva_sig_mjj = DiJetMass;
					mva_sig_deta_jj = fabs(DiJetDEta);
					mva_sig_m4l = ZZMass;
					mva_sig_eta_star_Z1 = eta_Z1_star;
					mva_sig_eta_star_Z2 = eta_Z2_star;
					mva_sig_R_pt_hard = R_pt_hard;
					mva_sig_R_pt_jet = R_pt_jet;

					mva_sig_abs_etajet_min = abs_etajet_min;
					mva_sig_abs_etajet_max = abs_etajet_max;
					mva_sig_abs_etalep_min = abs_etalep_min;
					mva_sig_abs_etalep_max = abs_etalep_max;
					mva_sig_delta_phi_ZZ = delta_phi_ZZ;
					mva_sig_rapidity_Z1 = rapidity_Z1;
					mva_sig_rapidity_Z2 = rapidity_Z2;
					mva_sig_rapidity_j1 = rapidity_j1;
					mva_sig_rapidity_j2 = rapidity_j2;
					mva_sig_pt_Z1 = pt_Z1;
					mva_sig_pt_Z2 = pt_Z2;
					mva_sig_pt_l3 = pt_l3;
					mva_sig_jet1_qg_tagger = jet_qg_tagger->at(0);
					mva_sig_jet2_qg_tagger = jet_qg_tagger->at(1);
					mva_sig_dbkg_kin = dbkg_kin;
					mva_sig_eta_j1 = JetEta->at(0);
					mva_sig_eta_j2 = JetEta->at(1);
					mva_sig_pt_jet1 = JetPt->at(0);
					mva_sig_pt_jet2 = JetPt->at(1);
					mva_sig_eta_j_sum = JetEta->at(0) + JetEta->at(1);
					mva_sig_mjj_over_detajj = DiJetMass/fabs(DiJetDEta);
					mva_sig_abs_etajet_sum = fabs(JetEta->at(0)) + fabs(JetEta->at(1));
					mva_sig_weight = weight;

					tree_sig->Fill();
				}
				if (j==5)	//qq
				{
					mva_qq_mjj = DiJetMass;
					mva_qq_deta_jj = fabs(DiJetDEta);
					mva_qq_m4l = ZZMass;
					mva_qq_eta_star_Z1 = eta_Z1_star;
					mva_qq_eta_star_Z2 = eta_Z2_star;
					mva_qq_R_pt_hard = R_pt_hard;
					mva_qq_R_pt_jet = R_pt_jet;

					mva_qq_abs_etajet_min = abs_etajet_min;
					mva_qq_abs_etajet_max = abs_etajet_max;
					mva_qq_abs_etalep_min = abs_etalep_min;
					mva_qq_abs_etalep_max = abs_etalep_max;
					mva_qq_delta_phi_ZZ = delta_phi_ZZ;
					mva_qq_rapidity_Z1 = rapidity_Z1;
					mva_qq_rapidity_Z2 = rapidity_Z2;
					mva_qq_rapidity_j1 = rapidity_j1;
					mva_qq_rapidity_j2 = rapidity_j2;
					mva_qq_pt_Z1 = pt_Z1;
					mva_qq_pt_Z2 = pt_Z2;
					mva_qq_pt_l3 = pt_l3;
					mva_qq_jet1_qg_tagger = jet_qg_tagger->at(0);
					mva_qq_jet2_qg_tagger = jet_qg_tagger->at(1);
					mva_qq_dbkg_kin = dbkg_kin;
					mva_qq_eta_j1 = JetEta->at(0);
					mva_qq_eta_j2 = JetEta->at(1);
					mva_qq_pt_jet1 = JetPt->at(0);
					mva_qq_pt_jet2 = JetPt->at(1);
					mva_qq_eta_j_sum = JetEta->at(0) + JetEta->at(1);
					mva_qq_mjj_over_detajj = DiJetMass/fabs(DiJetDEta);
					mva_qq_abs_etajet_sum = fabs(JetEta->at(0)) + fabs(JetEta->at(1));
					mva_qq_weight = weight;

					tree_qq->Fill();
				}
				if (j==1)	//gg
				{
					mva_gg_mjj = DiJetMass;
					mva_gg_deta_jj = fabs(DiJetDEta);
					mva_gg_m4l = ZZMass;
					mva_gg_eta_star_Z1 = eta_Z1_star;
					mva_gg_eta_star_Z2 = eta_Z2_star;
					mva_gg_R_pt_hard = R_pt_hard;
					mva_gg_R_pt_jet = R_pt_jet;

					mva_gg_abs_etajet_min = abs_etajet_min;
					mva_gg_abs_etajet_max = abs_etajet_max;
					mva_gg_abs_etalep_min = abs_etalep_min;
					mva_gg_abs_etalep_max = abs_etalep_max;
					mva_gg_delta_phi_ZZ = delta_phi_ZZ;
					mva_gg_rapidity_Z1 = rapidity_Z1;
					mva_gg_rapidity_Z2 = rapidity_Z2;
					mva_gg_rapidity_j1 = rapidity_j1;
					mva_gg_rapidity_j2 = rapidity_j2;
					mva_gg_pt_Z1 = pt_Z1;
					mva_gg_pt_Z2 = pt_Z2;
					mva_gg_pt_l3 = pt_l3;
					mva_gg_jet1_qg_tagger = jet_qg_tagger->at(0);
					mva_gg_jet2_qg_tagger = jet_qg_tagger->at(1);
					mva_gg_dbkg_kin = dbkg_kin;
					mva_gg_eta_j1 = JetEta->at(0);
					mva_gg_eta_j2 = JetEta->at(1);
					mva_gg_pt_jet1 = JetPt->at(0);
					mva_gg_pt_jet2 = JetPt->at(1);
					mva_gg_eta_j_sum = JetEta->at(0) + JetEta->at(1);
					mva_gg_mjj_over_detajj = DiJetMass/fabs(DiJetDEta);
					mva_gg_abs_etajet_sum = fabs(JetEta->at(0)) + fabs(JetEta->at(1));
					mva_gg_weight = weight;

					tree_gg->Fill();
				}
				if (j==4)	//ttZ,WWZ
				{
					mva_ttzwwz_mjj = DiJetMass;
					mva_ttzwwz_deta_jj = fabs(DiJetDEta);
					mva_ttzwwz_m4l = ZZMass;
					mva_ttzwwz_eta_star_Z1 = eta_Z1_star;
					mva_ttzwwz_eta_star_Z2 = eta_Z2_star;
					mva_ttzwwz_R_pt_hard = R_pt_hard;
					mva_ttzwwz_R_pt_jet = R_pt_jet;

					mva_ttzwwz_abs_etajet_min = abs_etajet_min;
					mva_ttzwwz_abs_etajet_max = abs_etajet_max;
					mva_ttzwwz_abs_etalep_min = abs_etalep_min;
					mva_ttzwwz_abs_etalep_max = abs_etalep_max;
					mva_ttzwwz_delta_phi_ZZ = delta_phi_ZZ;
					mva_ttzwwz_rapidity_Z1 = rapidity_Z1;
					mva_ttzwwz_rapidity_Z2 = rapidity_Z2;
					mva_ttzwwz_rapidity_j1 = rapidity_j1;
					mva_ttzwwz_rapidity_j2 = rapidity_j2;
					mva_ttzwwz_pt_Z1 = pt_Z1;
					mva_ttzwwz_pt_Z2 = pt_Z2;
					mva_ttzwwz_pt_l3 = pt_l3;
					mva_ttzwwz_jet1_qg_tagger = jet_qg_tagger->at(0);
					mva_ttzwwz_jet2_qg_tagger = jet_qg_tagger->at(1);
					mva_ttzwwz_dbkg_kin = dbkg_kin;
					mva_ttzwwz_eta_j1 = JetEta->at(0);
					mva_ttzwwz_eta_j2 = JetEta->at(1);
					mva_ttzwwz_pt_jet1 = JetPt->at(0);
					mva_ttzwwz_pt_jet2 = JetPt->at(1);
					mva_ttzwwz_eta_j_sum = JetEta->at(0) + JetEta->at(1);
					mva_ttzwwz_mjj_over_detajj = DiJetMass/fabs(DiJetDEta);
					mva_ttzwwz_abs_etajet_sum = fabs(JetEta->at(0)) + fabs(JetEta->at(1));
					mva_ttzwwz_weight = weight;

					tree_ttzwwz->Fill();
				}
				if (j==3)	//data
				{
					mva_data_mjj = DiJetMass;
					mva_data_deta_jj = fabs(DiJetDEta);
					mva_data_m4l = ZZMass;
					mva_data_eta_star_Z1 = eta_Z1_star;
					mva_data_eta_star_Z2 = eta_Z2_star;
					mva_data_R_pt_hard = R_pt_hard;
					mva_data_R_pt_jet = R_pt_jet;

					mva_data_abs_etajet_min = abs_etajet_min;
					mva_data_abs_etajet_max = abs_etajet_max;
					mva_data_abs_etalep_min = abs_etalep_min;
					mva_data_abs_etalep_max = abs_etalep_max;
					mva_data_delta_phi_ZZ = delta_phi_ZZ;
					mva_data_rapidity_Z1 = rapidity_Z1;
					mva_data_rapidity_Z2 = rapidity_Z2;
					mva_data_rapidity_j1 = rapidity_j1;
					mva_data_rapidity_j2 = rapidity_j2;
					mva_data_pt_Z1 = pt_Z1;
					mva_data_pt_Z2 = pt_Z2;
					mva_data_pt_l3 = pt_l3;
					mva_data_jet1_qg_tagger = jet_qg_tagger->at(0);
					mva_data_jet2_qg_tagger = jet_qg_tagger->at(1);
					mva_data_dbkg_kin = dbkg_kin;
					mva_data_eta_j1 = JetEta->at(0);
					mva_data_eta_j2 = JetEta->at(1);
					mva_data_pt_jet1 = JetPt->at(0);
					mva_data_pt_jet2 = JetPt->at(1);
					mva_data_eta_j_sum = JetEta->at(0) + JetEta->at(1);
					mva_data_mjj_over_detajj = DiJetMass/fabs(DiJetDEta);
					mva_data_abs_etajet_sum = fabs(JetEta->at(0)) + fabs(JetEta->at(1));
					mva_data_weight = weight;

					tree_data->Fill();
				}
	    	}
		}//entries loop  end
	}//file loop  end
       
	//ZX CONTRIBUTION
	  
	TChain *tqqzz_zx= new TChain("candTree");
	sprintf(filename,"/home/llr/cms/giljanovic/scratch/VBS/CMSSW_10_2_15/src/vbs_analysis/4l_channel/data_driven_MC/ZX%d%s_%s.root",year,theExtra.c_str(), jet_pt_cut.c_str());
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
			Float_t bins_m4l[] = {0, 1 };

			sprintf(filename,"hzx%d",iv);   	hzx[iv] = new TH1F(filename,"",1,bins_m4l); //full histogram
	  		sprintf(filename,"hzx_ee%d",iv);   	hzx_ee[iv] = new TH1F(filename,"",1,bins_m4l); //full histogram
	  		sprintf(filename,"hzx_mm%d",iv);   	hzx_mm[iv] = new TH1F(filename,"",1,bins_m4l); //full histogram
	  		sprintf(filename,"hzx_em%d",iv);   	hzx_em[iv] = new TH1F(filename,"",1,bins_m4l); //full histogram
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

	// -------------------------------------------- my variables -----------------------------------------

	float eta_Z1_star_zx, eta_Z2_star_zx, R_pt_hard_zx, R_pt_jet_zx;
	float abs_etajet_min_zx, abs_etajet_max_zx, abs_etalep_min_zx, abs_etalep_max_zx;
	float delta_phi_ZZ_zx;

	float rapidity_Z1_zx, rapidity_Z2_zx, rapidity_j1_zx, rapidity_j2_zx;
	float pt_Z1_zx, pt_Z2_zx, pt_l3_zx;

	float j1_qg_tagger_zx, j2_qg_tagger_zx;

	// ---------------------------------------------------------------------------------------------------

	// ---------------------------------------- my branches ----------------------------------------------

	tqqzz_zx->SetBranchAddress("eta_Z1_star",&eta_Z1_star_zx);
	tqqzz_zx->SetBranchAddress("eta_Z2_star",&eta_Z2_star_zx);
	tqqzz_zx->SetBranchAddress("R_pt_hard",&R_pt_hard_zx);
	tqqzz_zx->SetBranchAddress("R_pt_jet",&R_pt_jet_zx);

	tqqzz_zx->SetBranchAddress("abs_etajet_min",&abs_etajet_min_zx);
	tqqzz_zx->SetBranchAddress("abs_etajet_max",&abs_etajet_max_zx);
	tqqzz_zx->SetBranchAddress("abs_etalep_min",&abs_etalep_min_zx);
	tqqzz_zx->SetBranchAddress("abs_etalep_max",&abs_etalep_max_zx);

	tqqzz_zx->SetBranchAddress("delta_phi_ZZ",&delta_phi_ZZ_zx);

	tqqzz_zx->SetBranchAddress("rapidity_Z1",&rapidity_Z1_zx);
	tqqzz_zx->SetBranchAddress("rapidity_Z2",&rapidity_Z2_zx);
	tqqzz_zx->SetBranchAddress("rapidity_j1",&rapidity_j1_zx);
	tqqzz_zx->SetBranchAddress("rapidity_j2",&rapidity_j2_zx);

	tqqzz_zx->SetBranchAddress("pt_Z1",&pt_Z1_zx);
	tqqzz_zx->SetBranchAddress("pt_Z2",&pt_Z2_zx);
	tqqzz_zx->SetBranchAddress("pt_l3",&pt_l3_zx);

	tqqzz_zx->SetBranchAddress("j1_qg_tagger",&j1_qg_tagger_zx);
	tqqzz_zx->SetBranchAddress("j2_qg_tagger",&j2_qg_tagger_zx);

	// ---------------------------------------------------------------------------------------------------

	//entries loop
	for(int i=0;i<tqqzz_zx->GetEntries();i++)
	{
	  	tqqzz_zx->GetEntry(i);
	  
	  	//1D kin var hist fill
	  	for(int il = 0; il < vars+2; il++)
		{
            int iv = il; 
	    	if (iv == 0) var_zx = dbkg_kin_zx;
	    	//if (iv == 1) var_zx = ZZMass_zx;
			if (iv == 1) var_zx = 0.5;
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
				var_zx = etajet2_zx;
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
				var_zx = DiJetMass_zx/fabs(DiJetDEta_zx);
				iv = 11;
			}	
			if (il == 14)
			{
				var_zx = eta_Z1_star_zx;
				iv = 12;
			}
			if (il == 15)
			{
				var_zx = eta_Z2_star_zx;
				iv = 13;
			}
			if (il == 16)
			{
				var_zx = R_pt_hard_zx;
				iv = 14;
			}
			if (il == 17)
			{
				var_zx = R_pt_jet_zx;
				iv = 15;
			}
			if (il == 18)
			{
				var_zx = abs_etajet_min_zx;
				iv = 16;
			}
			if (il == 19)
			{
				var_zx = abs_etajet_max_zx;
				iv = 17;
			}
			if (il == 20)
			{
				var_zx = abs_etalep_min_zx;
				iv = 18;
			}
			if (il == 21)
			{
				var_zx = abs_etalep_max_zx;
				iv = 19;
			}
			if (il == 22)
			{
				var_zx = delta_phi_ZZ_zx;
				iv = 20;
			}
			if (il == 23)
			{
				var_zx = rapidity_Z1_zx;
				iv = 21;
			}
			if (il == 24)
			{
				var_zx = rapidity_Z2_zx;
				iv = 22;
			}
			if (il == 25)
			{
				var_zx = rapidity_j1_zx;
				iv = 23;
			}
			if (il == 26)
			{
				var_zx = rapidity_j2_zx;
				iv = 24;
			}
			if (il == 27)
			{
				var_zx = pt_Z1_zx;
				iv = 25;
			}
			if (il == 28)
			{
				var_zx = pt_Z2_zx;
				iv = 26;
			}
			if (il == 29)
			{
				var_zx = pt_l3_zx;
				iv = 27;
			}
			if (il == 30)
			{
				var_zx = j1_qg_tagger_zx;
				iv = 28;
			}
			if (il == 31)
			{
				var_zx = j2_qg_tagger_zx;
				iv = 29;
			}
			if (il == 32)
			{
				var_zx = ZZMass_zx;
				iv = 30;
			}
			if (il == 33)
			{
				var_zx = j1_qg_tagger_zx;
				iv = 31;
			}
			if (il == 34)
			{
				var_zx = j2_qg_tagger_zx;
				iv = 32;
			}
			if (il == 35)
			{
				var_zx = fabs(etajet1_zx) + fabs(etajet2_zx);
				iv = 33;
			}	

	    	if (fabs(weight_zx) < 100000.) 
			{
				if (il == 1)
					if (var_zx < 1200 || var_zx > 1400) continue;
	      		hzx[iv]->Fill(var_zx,weight_zx);
	      		if (chan_zx == 2) hzx_ee[iv]->Fill(var_zx,weight_zx);
	      		if (chan_zx == 1) hzx_mm[iv]->Fill(var_zx,weight_zx);
	      		if (chan_zx == 3) hzx_em[iv]->Fill(var_zx,weight_zx);
	    	}	    
	  	}

		mva_zx_mjj = DiJetMass_zx;
		mva_zx_deta_jj = fabs(DiJetDEta_zx);
		mva_zx_m4l = ZZMass_zx;
		mva_zx_eta_star_Z1 = eta_Z1_star_zx;
		mva_zx_eta_star_Z2 = eta_Z2_star_zx;
		mva_zx_R_pt_hard = R_pt_hard_zx;
		mva_zx_R_pt_jet = R_pt_jet_zx;

		mva_zx_abs_etajet_min = abs_etajet_min_zx;
		mva_zx_abs_etajet_max = abs_etajet_max_zx;
		mva_zx_abs_etalep_min = abs_etalep_min_zx;
		mva_zx_abs_etalep_max = abs_etalep_max_zx;
		mva_zx_delta_phi_ZZ = delta_phi_ZZ_zx;
		mva_zx_rapidity_Z1 = rapidity_Z1_zx;
		mva_zx_rapidity_Z2 = rapidity_Z2_zx;
		mva_zx_rapidity_j1 = rapidity_j1_zx;
		mva_zx_rapidity_j2 = rapidity_j2_zx;
		mva_zx_pt_Z1 = pt_Z1_zx;
		mva_zx_pt_Z2 = pt_Z2_zx;
		mva_zx_pt_l3 = pt_l3_zx;
		mva_zx_jet1_qg_tagger = j1_qg_tagger_zx;
		mva_zx_jet2_qg_tagger = j2_qg_tagger_zx;
		mva_zx_dbkg_kin = dbkg_kin_zx;
		mva_zx_eta_j1 = etajet1_zx;
		mva_zx_eta_j2 = etajet2_zx;
		mva_zx_pt_jet1 = ptjet1_zx;
		mva_zx_pt_jet2 = ptjet2_zx;
		mva_zx_eta_j_sum = etajet1_zx + etajet2_zx;
		mva_zx_mjj_over_detajj = DiJetMass_zx/fabs(DiJetDEta_zx);
		mva_zx_abs_etajet_sum = fabs(etajet1_zx) + fabs(etajet2_zx);
		mva_zx_weight = weight_zx;

		tree_zx->Fill();
	}
	
	for (int it=0; it < 5; it++) 
	{
	  	fnew[it]->cd();
	  	if (it < 3 || it == 4) 
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
	    	if (it==0) hqqzz[il]->Write();  
	    	if (it==1) hggzz[il]->Write();
	    	if (it==2) hvbs[il]->Write();
	    	if (it==3) {hdata[il]->Write();    hzx[il]->Write();}
	    	if (it==4) httzwzz[il]->Write();
	  	}
	  	tnew[it]->Write();
	  	fnew[it]->Close();
	}

	
	//INTEGRAL CHECK
	sprintf(filename,"MCyields_%d.txt",year);
	ofstream yields(filename,std::fstream::app);
    sprintf(filename,"datayields_%d.txt",year);
	ofstream yields2(filename,std::fstream::app);
	sprintf(filename,"MCyields_highMELA_%d.txt",year);
	ofstream yields3(filename,std::fstream::app);
    sprintf(filename,"datayields_highMELA_%d.txt",year);
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
	    
	    	yields << vbsmumu << "\t" << hzx_mm[iv]->Integral() << "\t" << hqqzz_mm[iv]->Integral()  << "\t" << hggzz_mm[iv]->Integral() << "\t" << httzwzz_mm[iv]->Integral() << "\t" << vbsee << "\t" << hzx_ee[iv]->Integral() << "\t" <<  hqqzz_ee[iv]->Integral()  << "\t" << hggzz_ee[iv]->Integral() << "\t" << httzwzz_ee[iv]->Integral() << "\t" << vbsemu << "\t" << hzx_em[iv]->Integral() << "\t" << hqqzz_em[iv]->Integral()  << "\t" << hggzz_em[iv]->Integral() << "\t" << httzwzz_em[iv]->Integral() << endl;
	    	yields2 << hdata_mm[iv]->Integral() << "\t" << hdata_ee[iv]->Integral() << "\t" << hdata_em[iv]->Integral() << endl;
	    	yields3 << vbsmumu_high << "\t" << hzx_mm[iv]->Integral(19,20) << "\t" << hqqzz_mm[iv]->Integral(19,20)  << "\t" << hggzz_mm[iv]->Integral(19,20) << "\t" << "\t" <<  httzwzz_mm[iv]->Integral(19,20) << vbsee_high << "\t" << hzx_ee[iv]->Integral(19,20) << "\t" <<  hqqzz_ee[iv]->Integral(19,20)  << "\t" << hggzz_ee[iv]->Integral(19,20) << "\t" << httzwzz[iv]->Integral(19,20) << "\t" << vbsemu_high << "\t" << hzx_em[iv]->Integral(19,20) << "\t" << hqqzz_em[iv]->Integral(19,20) << "\t" << hggzz_em[iv]->Integral(19,20) << "\t" <<  httzwzz_em[iv]->Integral(19,20) << endl;
	    	yields4 << hdata_mm[iv]->Integral(19,20) << "\t" << hdata_ee[iv]->Integral(19,20) << "\t" << hdata_em[iv]->Integral(19,20) << endl;
	  	}
	  
	  	float zx_integral = (2.00933+1.97468+3.96839)*lumi/35.9E03;
        const float powheg_integral = hqqzz_powheg[iv]->Integral();
	  
	  	/* if (year == 2018 && useMCatNLO > 0) {
	  	  hqqzz_powheg[iv]->Multiply(hqqzz_powheg[iv],hqqzz[iv]);
	  	  hqqzz_powheg[iv]->Divide(hqqzz_powheg[iv],httzwzz[iv]);
	  	  }  */
	  	//hzx[iv]->Scale(zx_integral/hzx[iv]->Integral());
	  
	  	//data normalisation
	  	if (useMCatNLO == 2) 
		{
	    	// if (year == 2018) hqqzz_powheg[iv]->Scale(powheg_integral/hqqzz_powheg[iv]->Integral());
	    	// else 
	    	hqqzz[iv]->Scale(powheg_integral/hqqzz[iv]->Integral());
	  	}


		// ---------------------------------------------- saving mjj, detajj and kD to file ------------------------------------

		if (iv == 2)
		{
			TString name = "./onlymjjCut_jet_pt_gt_30/mjj_" + to_string(year) + ".root";
			TFile *h_all_contributions = new TFile(name, "recreate");
			TH1F *h_ewk = (TH1F*) hvbs[iv]->Clone();
			h_ewk->SetName("vbs");
			h_ewk->Write();

			TH1F *h_qq = (TH1F*) hqqzz[iv]->Clone();
			h_qq->SetName("QCD_qq");
			h_qq->Write();

			TH1F *h_gg = (TH1F*) hggzz[iv]->Clone();
			h_gg->SetName("QCD_gg");
			h_gg->Write();

			TH1F *h_data = (TH1F*) hdata[iv]->Clone();
			h_data->SetName("data");
			h_data->Write();

			TH1F *h_zx = (TH1F*) hzx[iv]->Clone();
			h_zx->SetName("ZX");
			h_zx->Write();

			TH1F *h_ttz = (TH1F*) httz[iv]->Clone();
			h_ttz->SetName("ttZ");
			h_ttz->Write();

			TH1F *h_wwz = (TH1F*) hwwz[iv]->Clone();
			h_wwz->SetName("WWZ");
			h_wwz->Write();

			TH1F *h_wzz = (TH1F*) hwzz[iv]->Clone();
			h_wzz->SetName("WZZ");
			h_wzz->Write();
			
			h_all_contributions->Close();
		}

		if (iv == 3)
		{
			TString name = "./onlymjjCut_jet_pt_gt_30/detajj_" + to_string(year) + ".root";
			TFile *h_all_contributions = new TFile(name, "recreate");
			TH1F *h_ewk = (TH1F*) hvbs[iv]->Clone();
			h_ewk->SetName("vbs");
			h_ewk->Write();

			TH1F *h_qq = (TH1F*) hqqzz[iv]->Clone();
			h_qq->SetName("QCD_qq");
			h_qq->Write();

			TH1F *h_gg = (TH1F*) hggzz[iv]->Clone();
			h_gg->SetName("QCD_gg");
			h_gg->Write();

			TH1F *h_data = (TH1F*) hdata[iv]->Clone();
			h_data->SetName("data");
			h_data->Write();

			TH1F *h_zx = (TH1F*) hzx[iv]->Clone();
			h_zx->SetName("ZX");
			h_zx->Write();

			TH1F *h_ttz = (TH1F*) httz[iv]->Clone();
			h_ttz->SetName("ttZ");
			h_ttz->Write();

			TH1F *h_wwz = (TH1F*) hwwz[iv]->Clone();
			h_wwz->SetName("WWZ");
			h_wwz->Write();

			TH1F *h_wzz = (TH1F*) hwzz[iv]->Clone();
			h_wzz->SetName("WZZ");
			h_wzz->Write();
			
			h_all_contributions->Close();
		}

		/*if (iv == 0)
		{
			TString name = "./onlymjjCut_jet_pt_gt_30/kd_" + to_string(year) + ".root";
			TFile *h_all_contributions = new TFile(name, "recreate");
			TH1F *h_ewk = (TH1F*) hvbs[iv]->Clone();
			h_ewk->SetName("bkg_vbs");
			h_ewk->Write();

			TH1F *h_qq = (TH1F*) hqqzz[iv]->Clone();
			h_qq->SetName("bkg_qqzz");
			h_qq->Write();

			TH1F *h_gg = (TH1F*) hggzz[iv]->Clone();
			h_gg->SetName("bkg_ggzz");
			h_gg->Write();

			TH1F *h_data = (TH1F*) hdata[iv]->Clone();
			h_data->SetName("data_obs");
			h_data->Write();

			TH1F *h_zx = (TH1F*) hzx[iv]->Clone();
			h_zx->SetName("bkg_zjet");
			h_zx->Write();

			TH1F *h_ttz = (TH1F*) httz[iv]->Clone();
			TH1F *h_wwz = (TH1F*) hwwz[iv]->Clone();
			TH1F *h_ttzwwz = (TH1F*) httz[iv]->Clone();
			h_ttzwwz->Add(hwwz[iv]);
			h_ttzwwz->SetName("bkg_ttzwzz");
			h_ttzwwz->Write();
			
			h_all_contributions->Close();
		}*/

		// ------------------------------------------ end of saving mjj, detajj and kD to file ---------------------------------

		//saving EWK, qq, gg and data histogram in root file for aQGC part
		if (iv == 1 && calculate_aQGC_limits == false)
		{
			TString name = "./aQGC/raw_histos/all_contributions_" + to_string(year) + ".root";
			TFile *h_all_contributions = new TFile(name, "recreate");
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

			TH1F *h_ttz = (TH1F*) httz[iv]->Clone();
			h_ttz->SetName("ttZ");
			h_ttz->Write();

			TH1F *h_wwz = (TH1F*) hwwz[iv]->Clone();
			h_wwz->SetName("WWZ");
			h_wwz->Write();

			TH1F *h_wzz = (TH1F*) hwzz[iv]->Clone();
			h_wzz->SetName("WZZ");
			h_wzz->Write();
			
			h_all_contributions->Close();
		}

	  
	  	//HISTOGRAMS ADDED TO STACK
	  	hzx[iv]->SetFillColor(kGreen);
        httzwzz[iv]->Add(httzwzz[iv],hzx[iv],1,1);    //tt
	  	httzwzz[iv]->SetFillColor(kYellow);
	  	if (useMCatNLO == 0) hqqzz[iv]->Add(httzwzz[iv],hqqzz_powheg[iv],1,1); //real ew
	  	else hqqzz[iv]->Add(httzwzz[iv],hqqzz[iv],1,1); //real ew
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
			httzwzz[iv]->SetBinContent(httzwzz[iv]->GetNbinsX(), httzwzz[iv]->GetBinContent(httzwzz[iv]->GetNbinsX()) + httzwzz[iv]->GetBinContent(httzwzz[iv]->GetNbinsX() + 1));	//ttzwwz
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

			//aQGC_histos_file->cd();
			//hewk_FT8->Write();
			//aQGC_histos_file->Close();
		}
	  
	  	//add histograms to stack
	  	hs[iv]->Add(hsum2[iv],"hist");
	  	hs[iv]->Add(hsum1[iv],"hist");
	  	hs[iv]->Add(hqqzz[iv],"hist");
      	hs[iv]->Add(httzwzz[iv],"hist");
	  	hs[iv]->Add(hzx[iv],"hist");
	  	TH1F *hdatadivide = (TH1F*)hdata[iv]->Clone();
	  	hs[iv]->Add(hdata[iv],"E1");
		if (iv == 1 && calculate_aQGC_limits)
		{
			aQGC_histos_file->cd();
			hsum2[iv]->Write();
			hsum1[iv]->Write();
			hqqzz[iv]->Write();
			httzwzz[iv]->Write();
			hzx[iv]->Write();
			hdata[iv]->Write();
			hewk_FT8->Write();
			hewk_FT9->Write();
			aQGC_histos_file->Close();

			hs[iv]->Add(hewk_FT8);
			hs[iv]->Add(hewk_FT9);
		}
	  
	  	// draw the legend
	  	TLegend *legend=new TLegend(0.6,0.55,0.85,0.88);
	  	legend->SetTextFont(72);
	  	legend->SetTextSize(0.04);
	  	legend->AddEntry(hzx[iv],"Z+X","f");
		legend->AddEntry(httzwzz[iv],"t#bar{t}Z, WWZ, WZZ","f");
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
		set_y_axis_scale(hs[iv], iv, year);
		if (iv ==0 && year == 2016) hs[iv]->SetMaximum(48);
		if (iv ==0 && year == 2017) hs[iv]->SetMaximum(55);
		if (iv ==0 && year == 2018) hs[iv]->SetMaximum(75);
		

	  	hs[iv]->Draw("nostack"); //old
	  	if (drawSignal[iv] && !enriched) 
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
	  	if (useMCatNLO == 0) sprintf(filename,"onlymjjCut_%s/%s_plot_allPOWHEG_%d_%s.png",jet_pt_cut.c_str(), namegif[iv].c_str(),year, theExtra.c_str());      
	  	if (useMCatNLO == 1) sprintf(filename,"onlymjjCut_%s/%s_plot_allMCatNLO_%d_%s.png",jet_pt_cut.c_str(), namegif[iv].c_str(),year, theExtra.c_str());     
	  	if (useMCatNLO == 2) sprintf(filename,"onlymjjCut_%s/%s_plot_MCatNLOshape_POWHEGint_%d_%s.png",jet_pt_cut.c_str(), namegif[iv].c_str(),year, theExtra.c_str());
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
	//hZ1Mass->SetMaximum(350000);

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
	//hZ2Mass->SetMaximum(150000);

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
	//hZMass->SetMaximum(350000);

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
	//hZ1M_difference->SetMaximum(350000);
	c5->SaveAs("Z1Mass_difference.png");

	TCanvas *c6 = new TCanvas();
	c6->cd();
	hZ2M_difference->SetMaximum(2600);
	hZ2M_difference->SetLineColor(kRed);
	hZ2M_difference->SetLineWidth(3);
	hZ2M_difference->Draw();
	//hZ2M_difference->SetMaximum(350000);
	c6->SaveAs("Z2Mass_difference.png");

	// ------------------------------------- Write TMVA trees ----------------------------------------

	mva_file->cd();
	tree_sig->Write();
	tree_qq->Write();
	tree_gg->Write();
	tree_ttzwwz->Write();
	tree_data->Write();
	tree_zx->Write();

	mva_file->Close();
}
