#include "include/candTree.h"
#include "include/cConstants.h"
//#include "external_cConstants.h"
#include "include/bitops.h"
#include "include/FakeRates.h"
#include "include/FinalStates.h"
#include "include/Settings.h"
#include <TSpline.h>
#include <vector>
#include "TMath.h"
#include "TLorentzVector.h"
#include "../helper_functions.h"

using namespace std;

int FindFinalStateZX(short Z1Flav, short Z2Flav);

int main( int argc, char *argv[] )
{	
    int year = 2018;
	int enriched = 0;
	//int enriched = atoi(argv[2]);
	string pt_cut = "jet_pt_gt_30";

	string theExtra = "";
    if (enriched == 1) theExtra = "_VBSenr";
    if (enriched == 2) theExtra = "_superVBSenr";
    if (enriched == 3) theExtra = "_bkgdEnr";
    if (enriched == 4) theExtra = "_ptjet50";

    int passtot = 0;
    int n2e2mu = 0;
    int n2mu2e = 0;
    int n4e = 0;
    int n4mu = 0;
    int njet2 = 0;
    int njet3 = 0;
    int njet4 = 0;
    int njet5p = 0;

	vector<float> _fs_ROS_SS;
    _fs_ROS_SS.push_back(1.04);//4mu  
    _fs_ROS_SS.push_back(1.01);//4e                                                                         
    _fs_ROS_SS.push_back(1.04);//2e2mu                                                                            
    _fs_ROS_SS.push_back(1.00);//2mu2e   	
 
    char name[200];
	if (year == 2016) sprintf(name,"/data_cms_upgrade/giljanovic/VBS/cut_based/2016/FakeRates/FakeRates_SS_2016.root");
	if (year == 2017) sprintf(name,"/data_cms_upgrade/giljanovic/VBS/cut_based/2017/FakeRates/FakeRates_SS_2017.root"); 
    if (year == 2018) sprintf(name,"/data_cms_upgrade/giljanovic/VBS/cut_based/2018/FakeRates/FakeRates_SS_2018.root");

	FakeRates *FR = new FakeRates(name);
    
	TChain *t = new TChain("CRZLLTree/candTree");
	if (year == 2016) t->Add("/data_cms_upgrade/giljanovic/VBS/cut_based/2016/AllData/ZZ4lAnalysis.root");
    if (year == 2017) t->Add("/data_cms_upgrade/giljanovic/VBS/cut_based/2017/AllData/ZZ4lAnalysis.root");
    if (year == 2018) t->Add("/data_cms_upgrade/giljanovic/VBS/cut_based/2018/AllData/ZZ4lAnalysis.root");

	float c_constant = 14.0;
    //if (year == 2017) c_constant = 2.3;
	//if (year == 2018) c_constant = 2.3; 
        
	TFile* f_ = TFile::Open("/home/llr/cms/giljanovic/scratch/CMSSW_10_2_15/src/ZZAnalysis/AnalysisStep/data/cconstants/SmoothKDConstant_m4l_DjjVBF13TeV.root");
	TSpline3* ts = (TSpline3*)(f_->Get("sp_gr_varReco_Constant_Smooth")->Clone());
	f_->Close();

	candTree data(t);
	Long64_t nentries = data.fChain->GetEntries();
	cout<< nentries<<endl;
	data.fChain->SetBranchStatus("*", 0);
	data.fChain->SetBranchStatus("DiJetMass", 1);
    data.fChain->SetBranchStatus("DiJetDEta", 1);
	data.fChain->SetBranchStatus("ZZMass", 1);
	data.fChain->SetBranchStatus("ZZMassErrCorr", 1);
	data.fChain->SetBranchStatus("p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal", 1);
	data.fChain->SetBranchStatus("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal", 1);
	data.fChain->SetBranchStatus("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal", 1);
	data.fChain->SetBranchStatus("p_m4l_SIG", 1);
	data.fChain->SetBranchStatus("p_m4l_BKG", 1);
	data.fChain->SetBranchStatus("CRflag", 1);
	data.fChain->SetBranchStatus("nCleanedJetsPt30", 1);
	data.fChain->SetBranchStatus("Z1Flav", 1);
	data.fChain->SetBranchStatus("Z2Flav", 1);
	data.fChain->SetBranchStatus("Z1Mass", 1);
	data.fChain->SetBranchStatus("Z2Mass", 1);
	data.fChain->SetBranchStatus("LepEta", 1);
	data.fChain->SetBranchStatus("LepPt", 1);
    data.fChain->SetBranchStatus("JetEta", 1);
	data.fChain->SetBranchStatus("JetPhi", 1);
	data.fChain->SetBranchStatus("JetPt", 1);      
	data.fChain->SetBranchStatus("LepLepId", 1);
	data.fChain->SetBranchStatus("p_QQB_BKG_MCFM", 1);
	data.fChain->SetBranchStatus("p_GG_SIG_ghg2_1_ghz1_1_JHUGen", 1);
	data.fChain->SetBranchStatus("nExtraLep", 1);
	data.fChain->SetBranchStatus("nCleanedJetsPt30BTagged_bTagSF", 1);

	// --------------------------------------------------------------------------------- my branches ---------------------------------------------------------------

	data.fChain->SetBranchStatus("LepLepId", 1);
	data.fChain->SetBranchStatus("LepPt", 1);
	data.fChain->SetBranchStatus("LepEta", 1);
	data.fChain->SetBranchStatus("LepPhi", 1);
	data.fChain->SetBranchStatus("JetQGLikelihood", 1);

	data.fChain->SetBranchStatus("Nvtx", 1);


	// ----------------------------------------------------------------------------- end of my branches ------------------------------------------------------------

	

	//my branches	
	data.fChain->SetBranchStatus("p_JJVBF_BKG_MCFM_JECNominal",1); //sign
	data.fChain->SetBranchStatus("p_JJQCD_BKG_MCFM_JECNominal",1); //bkg
	
	float dbkg_kin;
	float dbkg;
	float ZZMass_new;
    float dijmass_new;
    float dijeta_new;
    float ptjet1;
    float ptjet2;
    float etajet1;
    float etajet2;   
	int chan;
	int vbfcate;
	float weight;
	float d2jet;
	float d_2j;
	float ZZMassErrCorr_new;
	short njet;

	// --------------------------------------------------- my declarations ------------------------------------------------------------

	vector<short> *lepId = data.LepLepId;
	vector<float> *lepPt = data.LepPt;
	vector<float> *lepEta = data.LepEta;
	vector<float> *lepPhi = data.LepPhi;
	vector<float> *jet_qg_tagger = data.JetQGLikelihood;


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

	float j1_qg_tagger, j2_qg_tagger;

	int Nvtx = data.Nvtx;

	// ------------------------------------------------ end of my declarations --------------------------------------------------------

    sprintf(name,"ZX%d_%s.root",year,pt_cut.c_str());
	TFile *f = new TFile(name,"recreate");

	TTree *tnew =new TTree("candTree","");
	tnew->Branch("dbkg_kin",&dbkg_kin,"dbkg_kin/F");
	tnew->Branch("dbkg",&dbkg,"dbkg/F");
	tnew->Branch("ZZMass",&ZZMass_new,"ZZMass/F");
    tnew->Branch("DiJetMass",&dijmass_new,"DiJetMass/F");
	tnew->Branch("DiJetDEta",&dijeta_new,"DiJetDEta/F");
    tnew->Branch("ptjet1",&ptjet1,"ptjet1/F");
    tnew->Branch("ptjet2",&ptjet2,"ptjet2/F");
    tnew->Branch("etajet1",&etajet1,"etajet1/F");
    tnew->Branch("etajet2",&etajet2,"etajet2/F");
	tnew->Branch("ZZMassErrCorr",&ZZMassErrCorr_new,"ZZMassErrCorr/F");
	tnew->Branch("weight",&weight,"weight/F");
	tnew->Branch("chan",&chan,"chan/I");
	tnew->Branch("vbfcate",&vbfcate,"vbfcate/I");
	tnew->Branch("vbfMela",&d2jet,"vbfMela/F");
	tnew->Branch("d_2j",&d_2j,"d_2j/F");
	tnew->Branch("nCleanedJetsPt30",&njet,"nCleanedJetsPt30/S");

	// ----------------------- my variables ------------------------------
	tnew->Branch("eta_Z1_star",&eta_Z1_star,"eta_Z1_star/F");
	tnew->Branch("eta_Z2_star",&eta_Z2_star,"eta_Z2_star/F");
	tnew->Branch("R_pt_hard",&R_pt_hard,"R_pt_hard/F");
	tnew->Branch("R_pt_jet",&R_pt_jet,"R_pt_jet/F");
	tnew->Branch("abs_etajet_min",&abs_etajet_min,"abs_etajet_min/F");
	tnew->Branch("abs_etajet_max",&abs_etajet_max,"abs_etajet_max/F");
	tnew->Branch("abs_etalep_min",&abs_etalep_min,"abs_etalep_min/F");
	tnew->Branch("abs_etalep_max",&abs_etalep_max,"abs_etalep_max/F");
	tnew->Branch("delta_phi_ZZ",&delta_phi_ZZ,"delta_phi_ZZ/F");
	tnew->Branch("rapidity_Z1",&rapidity_Z1,"rapidity_Z1/F");
	tnew->Branch("rapidity_Z2",&rapidity_Z2,"rapidity_Z2/F");
	tnew->Branch("rapidity_j1",&rapidity_j1,"rapidity_j1/F");
	tnew->Branch("rapidity_j2",&rapidity_j2,"rapidity_j2/F");
	tnew->Branch("pt_Z1",&pt_Z1,"pt_Z1/F");
	tnew->Branch("pt_Z2",&pt_Z2,"pt_Z2/F");
	tnew->Branch("pt_l3",&pt_l3,"pt_l3/F");
	tnew->Branch("j1_qg_tagger",&j1_qg_tagger,"j1_qg_tagger/F");
	tnew->Branch("j2_qg_tagger",&j2_qg_tagger,"j2_qg_tagger/F");

	tnew->Branch("Nvtx",&Nvtx,"Nvtx/I");
	// --------------------------------------------------------------------

	for(Long64_t jentry=0; jentry<nentries;jentry++)
	{
		data.fChain->GetEntry(jentry);
		if(jentry%1000==0)
			cout<< jentry<<endl;

		if ( !data.CRflag ) continue;
		if ( !test_bit(data.CRflag, CRZLLss) ) continue;

		int _current_final_state = FindFinalStateZX(data.Z1Flav,data.Z2Flav);
		switch (_current_final_state){ 
			case Settings::fs2e2mu : chan=3; break;
			case Settings::fs2mu2e : chan=3; break;
			case Settings::fs4mu:    chan=1; break;
			case Settings::fs4e:     chan=2; break;
			default: cerr<<"ERROR! No final state";
		}
		njet = data.nCleanedJetsPt30;

		float c_Mela2j = getDVBF2jetsConstant(data.ZZMass);
		d2jet = 1./(1.+ c_Mela2j*data.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/data.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
		d_2j= data.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/(data.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal+data.p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal); 
		float WP_VBF2j = getDVBF2jetsWP(data.ZZMass, 0);
		
		
        if(data.DiJetMass>100 && data.ZZMass > 180 && data.nCleanedJetsPt30>1 && data.Z1Mass < 120 && data.Z1Mass > 60 && data.Z2Mass < 120 && data.Z2Mass > 60)
		{
			if (enriched == 1 && (data.DiJetMass < 400 || fabs(data.DiJetDEta) < 2.4)) continue;
		  	if (enriched == 2 && (data.DiJetMass < 400 || fabs(data.DiJetDEta) < 5.0)) continue;
		  	if (enriched == 3 && data.DiJetMass > 400 && fabs(data.DiJetDEta) > 2.4) continue;
		  	if (enriched == 4 && (data.JetPt->at(0) < 50 || data.JetPt->at(1) < 50)) continue;
		  	passtot++;
		  
		  	vbfcate=1;
		  	switch (_current_final_state)
			{ 
		     	case Settings::fs2e2mu : n2e2mu++; break;
		     	case Settings::fs2mu2e : n2mu2e++; break;
		     	case Settings::fs4mu:    n4mu++; break;
		     	case Settings::fs4e:     n4e++; break;
		     	default: cerr<<"ERROR! No final state";
		  	} 
		  	if (njet == 2) njet2++;
            else if (njet == 3) njet3++;		
		  	else if (njet == 4) njet4++;
            else njet5p++;
		  
		  	weight = _fs_ROS_SS.at(_current_final_state)*FR->GetFakeRate(data.LepPt->at(2),data.LepEta->at(2),data.LepLepId->at(2))*FR->GetFakeRate(data.LepPt->at(3),data.LepEta->at(3),data.LepLepId->at(3));

			//kinematic variable
		  	float c_mzz = c_constant*ts->Eval(data.ZZMass);
		  	dbkg_kin = data.p_JJVBF_BKG_MCFM_JECNominal / ( data.p_JJVBF_BKG_MCFM_JECNominal + data.p_JJQCD_BKG_MCFM_JECNominal*c_mzz );
		  
		  
		  	dbkg = data.p_GG_SIG_ghg2_1_ghz1_1_JHUGen*data.p_m4l_SIG / ( data.p_m4l_SIG*data.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + data.p_m4l_BKG*data.p_QQB_BKG_MCFM*getDbkgkinConstant(data.Z1Flav*data.Z2Flav,data.ZZMass) );
		  	ZZMass_new= data.ZZMass;
		  	dijmass_new= data.DiJetMass;
		  	dijeta_new= data.DiJetDEta;
		  	ptjet1=data.JetPt->at(0);
		  	ptjet2=data.JetPt->at(1);
		  	etajet1=data.JetEta->at(0);
		  	etajet2=data.JetEta->at(1);
		  	ZZMassErrCorr_new= data.ZZMassErrCorr;

			// ------------------------------------------------- calculate additional variables ------------------------------------------------------
			bool ZX = true;

			clear_vectors(electrons, muons, electrons_charge, muons_charge);
			create_electron_and_muon_objects(electrons, muons, electrons_charge, muons_charge, lepId, lepPt, lepEta, lepPhi);
			build_ZZ_pair(electrons, muons, electrons_charge, muons_charge, Z1, Z2, ZX);

			calculate_Zeppenfeld_Z(Z1, Z2, data.JetEta, eta_Z1_star, eta_Z2_star);
			calculate_R_pt_hard(Z1, Z2, data.JetEta, data.JetPhi, data.JetPt, R_pt_hard);
			calculate_R_pt_jet(data.JetEta, data.JetPhi, data.JetPt, R_pt_jet);
			calculate_min_max_jet_eta(data.JetEta, abs_etajet_min, abs_etajet_max);
			calculate_min_max_lepton_eta(electrons, muons, abs_etalep_min, abs_etalep_max);
			calculate_dphi_ZZ(Z1, Z2, delta_phi_ZZ);			
			calculate_rapidity_Z1_Z2(Z1, Z2, rapidity_Z1, rapidity_Z2);
			calculate_rapidity_j1_j2(data.JetEta, data.JetPhi, data.JetPt, rapidity_j1, rapidity_j2);
			calculate_pt_Z1_Z2_l3(Z1, Z2, lepPt, pt_Z1, pt_Z2, pt_l3);
			
			j1_qg_tagger = data.JetQGLikelihood->at(0);
			j2_qg_tagger = data.JetQGLikelihood->at(1);

			Nvtx = data.Nvtx;
	
		  	tnew->Fill();
		}	
	}
        
	cout << "n2mu2e : " << n2mu2e << "/" << passtot << " = " << (float)n2mu2e/(float)passtot << endl;
	cout << "n2e2mu : " << n2e2mu << "/" << passtot << " = " << (float)n2e2mu/(float)passtot << endl;
	cout << "n4e : " << n4e << "/" << passtot << " = " << (float)n4e/(float)passtot << endl;
	cout << "n4mu : " << n4mu << "/" << passtot << " = " << (float)n4mu/(float)passtot << endl;
    cout << "njet = 2 : " << njet2 << "/" << passtot << " = " << (float)njet2/(float)passtot << endl;      
	cout << "njet = 3 : " << njet3 << "/" << passtot << " = " << (float)njet3/(float)passtot << endl;      
	cout << "njet = 4 : " << njet4 << "/" << passtot << " = " << (float)njet4/(float)passtot << endl;      
	cout << "njet >= 5 : " << njet5p << "/" << passtot << " = " << (float)njet5p/(float)passtot << endl;
	f->cd();
	tnew->Write();
	f->Close();
}

int FindFinalStateZX(short Z1Flav, short Z2Flav)
{
	int final_state = -999;

	if ( Z1Flav == -121 )
	{
		if ( Z2Flav == +121 )
			final_state = Settings::fs4e;
		else if ( Z2Flav == +169 )
			final_state = Settings::fs2e2mu;
		else
			cerr << "[ERROR] in event " ;//<< RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z2Flav = " << Z2Flav << endl;
	} 
	else if ( Z1Flav == -169 )
	{
		if ( Z2Flav == +121 )
			final_state = Settings::fs2mu2e;
		else if ( Z2Flav == +169 )
			final_state = Settings::fs4mu;
		else
			cerr << "[ERROR] in event " << endl;//RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z2Flav = " << Z2Flav << endl;
	}
	else
	{
		cerr << "[ERROR] in event " << endl;//RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z1Flav = " << Z1Flav << endl;
	}

	return final_state;
}