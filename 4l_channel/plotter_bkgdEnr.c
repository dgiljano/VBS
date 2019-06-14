// #include "external_cConstants.h"
void plotter_bkgdEnr(int year = 2018, int useMCatNLO = 1){
  
       //useMCatNLO = 0 : use just POWHEG
       //useMCatNLO = 1 : use just aMCatNLO
       //useMCatNLO = 2 : use aMCatNLO for shape, POWHEG for integral

        float lumi = 35.9;
        if (year == 2017) lumi = 41.5;
	if (year == 2018) lumi = 59.7;

	static const int vars = 4;
        string titlex[vars] = {"K_{D}","M_{4l} [GeV]","M_{jj} [GeV]","#Delta #eta_{jj}"};        
        string titley[vars] = {"Events/0.025","Events/16 GeV","Events/22.5 GeV","Events/0.175"};     
	string namegif[vars] = {"Dbkgkin","m4l","mjj","detajj"};
        int bins[vars] = {40,40,40,40};
        float xmin[vars] = {0.,160.,100.,0.};
	float xmax[vars] = {1.,800.,1000.,8.};	
	bool drawSignal[vars] = {false,false,true,true};

	//histogram stack
        char filename[300]; char filetitle[300];
	THStack *hs[vars];
        for(int iv = 0; iv < vars; iv++){
	  sprintf(filename,"hs%d",iv);   
	  sprintf(filetitle,"CMS Preliminary                %2.1f fb^{-1};%s;%s",lumi,titlex[iv].c_str(),titley[iv].c_str());  
	  hs[iv] = new THStack(filename,filetitle);
	}  

        lumi *= 1000;

        //histograms

	TH1F *h_complete_data[vars]; //all data
	TH1F *h00[vars]; //bkg_kin<0.7 cut full background plot
	TH1F *h0[vars]; //data, because we are hiding higher energies in this phase
	TH1F *h1[vars]; //ew
	TH1F *h1bis[vars]; //ew+zx -> real ew
	TH1F *h2[vars]; //gg
	TH1F *h3[vars]; //vbs
	TH1F *h4[vars]; //gg+ew -> real gg
	TH1F *h5[vars]; //gg+ew+vbs
	TH1F *h6[vars]; //gg+ew+vbs+zx ->real vbs
	TH1F *h7[vars]; //qqzz e
	TH1F *h8[vars]; //qqzz mu
	TH1F *h9[vars]; //qqzz e mu
	TH1F *h10[vars];//ggzz e
	TH1F *h11[vars];//ggzz mu
	TH1F *h12[vars]; //ggzz e mu
	TH1F *h13[vars];//vbs e
	TH1F *h14[vars];//vbs mu
	TH1F *h15[vars];//vbs e mu
	TH1F *h0_ee[vars];//vbs e
	TH1F *h0_mm[vars];//vbs mu
	TH1F *h0_em[vars];//vbs e mu
	TH1F *hnum[vars];//for 2018 rescale
	TH1F *hden[vars];//for 2018 rescale

	for(int iv = 0; iv < vars; iv++){
	  sprintf(filename,"hcd%d",iv);   h_complete_data[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //all data
	  sprintf(filename,"h00_%d",iv);   h00[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //bkg_kin<0.7 cut full background plot
	  sprintf(filename,"h0_%d",iv);   h0[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //data, because we are hiding higher energies in this phase
	  sprintf(filename,"h1_%d",iv);   h1[iv] = new TH1F(filename,"h1;xtitle;ytitle",bins[iv],xmin[iv],xmax[iv]); //ew
	  sprintf(filename,"h1bis_%d",iv);   h1bis[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //ew+zx -> real ew
	  sprintf(filename,"h2_%d",iv);   h2[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //gg
	  sprintf(filename,"h3_%d",iv);   h3[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //vbs
	  sprintf(filename,"h4_%d",iv);   h4[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //gg+ew -> real gg
	  sprintf(filename,"h5_%d",iv);   h5[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //gg+ew+vbs
	  sprintf(filename,"h6_%d",iv);   h6[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //gg+ew+vbs+zx ->real vbs
	  sprintf(filename,"h7_%d",iv);   h7[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //qqzz e
	  sprintf(filename,"h8_%d",iv);   h8[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //qqzz mu
	  sprintf(filename,"h9_%d",iv);   h9[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //qqzz e mu
	  sprintf(filename,"h10_%d",iv);   h10[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//ggzz e
	  sprintf(filename,"h11_%d",iv);   h11[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//ggzz mu
	  sprintf(filename,"h12_%d",iv);  h12[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //ggzz e mu
	  sprintf(filename,"h13_%d",iv); h13[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//vbs e
	  sprintf(filename,"h14_%d",iv); h14[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//vbs mu
	  sprintf(filename,"h15_%d",iv); h15[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//vbs e mu
	  sprintf(filename,"h0ee_%d",iv); h0_ee[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//vbs e
	  sprintf(filename,"h0mm_%d",iv); h0_mm[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//vbs mu
	  sprintf(filename,"h0em_%d",iv); h0_em[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//vbs e mu
	  sprintf(filename,"hnum_%d",iv); hnum[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//for 2018 rescale
	  sprintf(filename,"hden_%d",iv); hden[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]);//for 2018 rescale
	}   
	gStyle->SetPalette(1);
	TFile *input_file;

	float c_constant = 0.05;
        if (year == 2017) c_constant = 2.3;
	if (year == 2018) c_constant = 2.3;   

        // find available samples  
	int nSamp = 0;  
	TString rootname[40];
 	sprintf(filename,"samples%d.txt",year);

	ifstream parInput(filename);
        
	if (parInput.is_open()) {
	  while ( parInput.good() ) {
	    parInput >> rootname[nSamp]; 
            cout << nSamp << " " <<  rootname[nSamp] << endl; 
	    nSamp++;
	  }
	  parInput.close();
	} 

	//for loop for different samples
	for(int is = 0; is < nSamp-1; is++){

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
	  if (rootname[is].Contains("ggTo")) j = 1;      
          if (rootname[is].Contains("VBFTo")) j = 2;
	  if (rootname[is].Contains("amcatnlo")) j = 4;
	  if (year == 2018 && rootname[is].Contains("180721_2017")) j = 5;

	  //histogram declaration
	  //TH1F *kin_zz = new TH1F("kin_zz","",bins,xmin,xmax); //was 100 bins
	  
	  //original variable declarations
	  float ZZPt,ZZMass,DiJetMass,DiJetDEta;
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
          tqqzz->SetBranchAddress("DiJetDEta",&DiJetDEta);
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
	  int enne = tqqzz->GetEntries();
  
	  for(int i=0;i<enne;i++){
	    tqqzz->GetEntry(i);
	    
	    //unique selection condition (see paper page 8) & DiJetMass condition
	    if(((fabs(DiJetDEta) > 1 && fabs(DiJetDEta) < 2.4) || (DiJetMass>100 && DiJetMass<400)) && nExtraLep==0 && ZZMass > 160 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0))){
	      
	      //unique selection
	      //	if(nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0))){
	      
	      //set vbf_category
	      vbfcate=1;
	      //weight definition
	      //KFactorEWKqqZZ = 1;
	      //KFactorQCDqqZZ_M = 1;
	      //weight=1;
	      weight= (xsec*KFactorEWKqqZZ * overallEventWeight*KFactorQCDqqZZ_M*lumi)/(gen_sum_weights);
	      if (j==1) weight= (xsec* overallEventWeight *KFactorQCDggzz_Nominal *lumi)/(gen_sum_weights);
	      if (j==3) weight=1.; 
	      //weight= (xsec*KFactorEWKqqZZ * overallEventWeight*KFactorQCDqqZZ_M*lumi)/(gen_sum_weights);
	      //dbkg_kin=overallEventWeight;

 //division in channels
	      if(abs(Z1Flav)==abs(Z2Flav) && abs(Z1Flav)==121){
		chan=2;
	      }
	      else if (abs(Z1Flav)==abs(Z2Flav) && abs(Z1Flav)!=121){
		chan=1;
	      }
	      else chan=3;
	      
	      //kin variable
	      for(int iv = 0; iv < vars; iv++){
		if (iv == 0) dbkg_kin = p_JJVBF_BKG_MCFM_JECNominal/(p_JJVBF_BKG_MCFM_JECNominal+ p_JJQCD_BKG_MCFM_JECNominal*c_constant);
		if (iv == 1) dbkg_kin = ZZMass;
		if (iv == 2) dbkg_kin = DiJetMass;
		if (iv == 3) dbkg_kin = fabs(DiJetDEta);
   	      
		//1D kin var hist fill
		//this is the normalization histogram
		if (j<3){
		  if (dbkg_kin <= 1.1) h00[iv]->Fill(dbkg_kin,weight);    
		  h_complete_data[iv]->Fill(dbkg_kin,weight);
		}
		
		if (j==3){
		  /*if (dbkg_kin <= 0.7)*/   h0[iv]->Fill(dbkg_kin);
		  h0[iv]->SetMarkerStyle(20);
		  if (chan == 2 ) h0_ee[iv]->Fill(dbkg_kin);
		  if (chan == 1) h0_mm[iv]->Fill(dbkg_kin);
		  if (chan == 3) h0_em[iv]->Fill(dbkg_kin);
		}
		if (j==0){                              //qqzz
		  h1[iv]->Fill(dbkg_kin,weight);
		  
		  if (chan == 1) h8[iv]->Fill(dbkg_kin,weight); //mu
		  if (chan == 2) h7[iv]->Fill(dbkg_kin,weight); //e
		  if (chan == 3)  h9[iv]->Fill(dbkg_kin,weight);          //mu+e
		}
		if (j==1){                              //gg
		  h2[iv]->Fill(dbkg_kin,weight);
		  h2[iv]->SetFillColor(kBlue);
		  
		  if (chan == 1) h11[iv]->Fill(dbkg_kin,weight); //mu
		  if (chan == 2) h10[iv]->Fill(dbkg_kin,weight); //e
		  if (chan == 3)  h12[iv]->Fill(dbkg_kin,weight);          //mu+e
		  
		}
		if (j==2){                              //vbs
		  h3[iv]->Fill(dbkg_kin,weight);
		  h3[iv]->SetFillColor(kMagenta);
		  
		  if (chan == 1) h14[iv]->Fill(dbkg_kin,weight); //mu
		  if (chan == 2) h13[iv]->Fill(dbkg_kin,weight); //e
		  if (chan == 3) h15[iv]->Fill(dbkg_kin,weight);          //mu+e
		  
		}
		if (j==4) hnum[iv]->Fill(dbkg_kin,weight);
		if (j==5) hden[iv]->Fill(dbkg_kin,weight);
		//kin_zz->GetYaxis()->SetTitle("Events/0.05");
		//add histogram to stack
		//for cycle ends here
	      }    
	    }
	      
	  }//entries loop  end
	  
	}//file loop  end
	//ZX CONTRIBUTION
	  
	TChain *tqqzz_zx= new TChain("candTree");
	sprintf(filename,"/afs/cern.ch/work/c/covarell/vbs2017/CMSSW_8_0_26_patch1/src/data_driven_MC/ZX%d_bkgdEnr.root",year); 
	tqqzz_zx->Add(filename);
	
	//histogram declaration
	TH1F *kin_zz_zx[vars]; 
	for(int iv = 0; iv < vars; iv++){
	  sprintf(filename,"kinzz%d",iv);   kin_zz_zx[iv] = new TH1F(filename,"",bins[iv],xmin[iv],xmax[iv]); //full histogram
	}
	
	//new variable declarations
	float var_zx, dbkg_kin_zx, ZZMass_zx, DiJetMass_zx, DiJetDEta_zx;
	float weight_zx;
	
	//original branch addresses
	tqqzz_zx->SetBranchAddress("dbkg_kin",&dbkg_kin_zx);
	tqqzz_zx->SetBranchAddress("ZZMass",&ZZMass_zx);
	tqqzz_zx->SetBranchAddress("DiJetMass",&DiJetMass_zx);
	tqqzz_zx->SetBranchAddress("DiJetDEta",&DiJetDEta_zx);
	tqqzz_zx->SetBranchAddress("weight",&weight_zx);
	
	//entries loop
	for(int i=0;i<tqqzz_zx->GetEntries();i++){
	  tqqzz_zx->GetEntry(i);
	  
	  //1D kin var hist fill
	  for(int iv = 0; iv < vars; iv++){
	    if (iv == 0) var_zx = dbkg_kin_zx;
	    if (iv == 1) var_zx = ZZMass_zx;
	    if (iv == 2) var_zx = DiJetMass_zx;
	    if (iv == 3) var_zx = fabs(DiJetDEta_zx);
	    if (year == 2016 || weight_zx < 1000000.) {
	      kin_zz_zx[iv]->Fill(var_zx,weight_zx);
	    }
	    if (var_zx<=0.7) h00[iv]->Fill(var_zx,weight_zx);
	    
	  }
	}
	
	//INTEGRAL CHECK
	
	for(int iv = 0; iv < vars; iv++){
	  cout << "Integral check" << endl;
	  cout <<"qqzz,         integral is " << h1[iv]->Integral() << endl;
	  cout <<"qqzz (4e),    integral is " << h7[iv]->Integral() << endl;
	  cout <<"qqzz (4mu),   integral is " << h8[iv]->Integral() << endl;
	  cout <<"qqzz (2e2mu), integral is " << h9[iv]->Integral() << endl;
	  cout <<"ggzz,         integral is " << h2[iv]->Integral() << endl;
	  cout <<"ggzz (4e),    integral is " << h10[iv]->Integral() << endl;
	  cout <<"ggzz (4mu),   integral is " << h11[iv]->Integral() << endl;
	  cout <<"ggzz (2e2mu), integral is " << h12[iv]->Integral() << endl;
	  cout <<"vbs,          integral is " << h3[iv]->Integral() << endl;
	  cout <<"vbs (4e),     integral is " << h13[iv]->Integral() << endl;
	  cout <<"vbs (4mu),    integral is " << h14[iv]->Integral() << endl;
	  cout <<"vbs (2e2mu),  integral is " << h15[iv]->Integral() << endl;
	  cout <<"data (4e),     integral is " << h0_ee[iv]->Integral() << endl;
	  cout <<"data (4mu),    integral is " << h0_mm[iv]->Integral() << endl;
	  cout <<"data (2e2mu),  integral is " << h0_em[iv]->Integral() << endl;
	  
	  float zx_integral = (2.00933+1.97468+3.96839)*lumi/35.9E03;
          const float powheg_integral = h1[iv]->Integral();
	  
	  if (year == 2018 && useMCatNLO > 0) {
	    h1[iv]->Multiply(h1[iv],hnum[iv]);
	    h1[iv]->Divide(h1[iv],hden[iv]);
	  }
	  kin_zz_zx[iv]->Scale(zx_integral/kin_zz_zx[iv]->Integral());
	  
	  //data normalisation
	  float data_integral = 1.;// h5->Integral();
	  float mc_integral = 1.; //h_complete_data->Integral();
	  kin_zz_zx[iv]->Scale(data_integral/mc_integral);
	  h1[iv]->Scale(data_integral/mc_integral);
	  h4[iv]->Scale(data_integral/mc_integral);
	  h5[iv]->Scale(data_integral/mc_integral);
	  hnum[iv]->Scale(data_integral/mc_integral);
	  if (useMCatNLO == 2) {
	    if (year == 2018) h1[iv]->Scale(powheg_integral/h1[iv]->Integral());
	    else hnum[iv]->Scale(powheg_integral/hnum[iv]->Integral());
	  }
	  
	  //HISTOGRAMS ADDED TO STACK
	  kin_zz_zx[iv]->SetFillColor(kGreen);
	  if (useMCatNLO == 0 || year == 2018) h1bis[iv]->Add(kin_zz_zx[iv],h1[iv],1,1); //real ew
	  else h1bis[iv]->Add(kin_zz_zx[iv],hnum[iv],1,1); //real ew
	  h1bis[iv]->SetFillColor(kCyan);
	  h4[iv]->Add(h1bis[iv],h2[iv],1,1); //real gg
	  h4[iv]->SetFillColor(kBlue);
	  h5[iv]->Add(h4[iv],h3[iv],1,1);    //real vbs
	  h5[iv]->SetFillColor(kMagenta);
	  
	  //add histograms to stack
	  hs[iv]->Add(h5[iv],"hist");
	  hs[iv]->Add(h4[iv],"hist");
	  hs[iv]->Add(h1bis[iv],"hist");
	  hs[iv]->Add(kin_zz_zx[iv],"hist");
	  TH1F *h0divide = (TH1F*)h0[iv]->Clone();
	  hs[iv]->Add(h0[iv],"E1");
	  
	  // draw the legend
	  TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
	  legend->SetTextFont(72);
	  legend->SetTextSize(0.04);
	  legend->AddEntry(kin_zz_zx[iv],"Z+X","f");
	  legend->AddEntry(h1bis[iv],"q#bar{q}#rightarrowZZ","f");
	  legend->AddEntry(h4[iv],"gg#rightarrowZZ","f");
	  legend->AddEntry(h5[iv],"VBS","f");
	  legend->AddEntry(h0[iv],"Data","lep");
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
	  hs[iv]->SetMaximum(12.*lumi/35.9E3); 
	  hs[iv]->Draw("nostack"); //old
	  if (drawSignal[iv]) {
	    TH1F* h77 = (TH1F*)h3[iv]->Clone();
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
	  pad2->SetBottomMargin(0.15);
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
	  h0copy->Divide(h5[iv]); //invert divide
	  h0copy->SetMarkerStyle(20);
	  string title = " ;"+titlex[iv]+"; Data/MC";
	  h0copy->SetTitle(title.c_str());
	  //h0copy->GetXaxis()->SetTitleSize(50);
	  //h0copy->GetXaxis()->SetLabelSize(35);
	  //h0copy->GetYaxis()->SetTitleSize(15);
	  //h0copy->GetYaxis()->SetLabelSize(15);
	  
	  //gStyle->SetLabelSize(2,"x").
	  h0copy->Draw("ep");
	  
	  //Orizontal line
	  TLine *line = new TLine(xmin[iv],1,xmax[iv],1);
	  line->SetLineColor(kRed);
	  line->SetLineStyle(2);
	  line->Draw("same");
	  
	  //close and print on file
	  c1->cd();
	  if (useMCatNLO == 0) sprintf(filename,"bkgdEnriched/%s_plot_allPOWHEG_%d.png",namegif[iv].c_str(),year);	  
	  if (useMCatNLO == 1) sprintf(filename,"bkgdEnriched/%s_plot_allMCatNLO_%d.png",namegif[iv].c_str(),year);	  
	  if (useMCatNLO == 2) sprintf(filename,"bkgdEnriched/%s_plot_MCatNLOshape_POWHEGint_%d.png",namegif[iv].c_str(),year);
	  gPad->Print(filename);
	}
}
