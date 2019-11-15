#include "include/candTree.h"
#include "include/cConstants.h"
//#include "external_cConstants.h"
#include "include/bitops.h"
#include "include/FakeRates.h"
#include "include/FinalStates.h"
#include "include/Settings.h"
#include <TSpline.h>
#include <vector>

using namespace std;

int FindFinalStateZX(short Z1Flav, short Z2Flav);

int main( int argc, char *argv[] ){
	
    int year = 2016;
    string pt_cut = "jet_pt_gt_30";

	vector<float> _fs_ROS_SS;
	_fs_ROS_SS.push_back(1.22);//4mu
	_fs_ROS_SS.push_back(0.27);//4e
	_fs_ROS_SS.push_back(1.30);//2e2mu
	_fs_ROS_SS.push_back(0.28);//2mu2e
	// 2016
 
    char name[200];
	if (year == 2016) sprintf(name,"/data_cms_upgrade/giljanovic/VBS/MC_2016/FakeRates/FakeRates_SS_2016_Legacy.root");
	if (year == 2017) sprintf(name,"/data_cms_upgrade/giljanovic/VBS/MC_2017/FakeRates/FakeRates_SS_2017_Legacy.root"); 
    if (year == 2018) sprintf(name,"/data_cms_upgrade/giljanovic/VBS/MC_2018/FakeRates/FakeRates_SS_2018_Legacy.root");

	FakeRates *FR = new FakeRates(name);
    
	TChain *t = new TChain("CRZLLTree/candTree");
    //2016
	if (year == 2016) t->Add("/data_cms_upgrade/giljanovic/VBS/Data_2016/AllData/ZZ4lAnalysis.root");
	//2017
    if (year == 2017) t->Add("/data_cms_upgrade/giljanovic/VBS/Data_2017/AllData/ZZ4lAnalysis.root");
	//2018
    if (year == 2018) t->Add("/data_cms_upgrade/giljanovic/VBS/Data_2018/AllData/ZZ4lAnalysis.root");

	float c_constant = 8.5;
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
	data.fChain->SetBranchStatus("JetPt", 1);      
	data.fChain->SetBranchStatus("LepLepId", 1);
	data.fChain->SetBranchStatus("p_QQB_BKG_MCFM", 1);
	data.fChain->SetBranchStatus("p_GG_SIG_ghg2_1_ghz1_1_JHUGen", 1);
	data.fChain->SetBranchStatus("nExtraLep", 1);
	data.fChain->SetBranchStatus("nCleanedJetsPt30BTagged_bTagSF", 1);

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

	for(Long64_t jentry=0; jentry<nentries;jentry++){
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
			case Settings::fs4e:    chan=2; break;
			default: cerr<<"ERROR! No final state";
		}
		njet = data.nCleanedJetsPt30;

		float c_Mela2j = getDVBF2jetsConstant(data.ZZMass);
		d2jet = 1./(1.+ c_Mela2j*data.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/data.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
		d_2j= data.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/(data.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal+data.p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal); 
		float WP_VBF2j = getDVBF2jetsWP(data.ZZMass, 0);
		
		//deta cut
		// if( data.ZZMass > 160  && fabs(data.DiJetDEta) > 1 && data.DiJetMass > 100 && data.nExtraLep==0 && (((data.nCleanedJetsPt30==2||data.nCleanedJetsPt30==3)&&data.nCleanedJetsPt30BTagged_bTagSF<=1) ||(data.nCleanedJetsPt30>=4&&data.nCleanedJetsPt30BTagged_bTagSF==0)) ){
		// no deta cut
		// if( data.ZZMass > 160  && data.DiJetMass > 100 && data.nExtraLep==0 && (((data.nCleanedJetsPt30==2||data.nCleanedJetsPt30==3)&&data.nCleanedJetsPt30BTagged_bTagSF<=1) ||(data.nCleanedJetsPt30>=4&&data.nCleanedJetsPt30BTagged_bTagSF==0)) ){
		// old fiducial region
        if(data.DiJetMass>100 && data.ZZMass > 180 && data.nCleanedJetsPt30>1 && data.Z1Mass < 120 && data.Z1Mass > 60 && data.Z2Mass < 120 && data.Z2Mass > 60){
                //if(data.DiJetMass>100 && data.ZZMass > 180 && data.nCleanedJetsPt30>1 && data.Z1Mass < 120 && data.Z1Mass > 60 && data.Z2Mass < 120 && data.Z2Mass > 60 && data.JetPt->at(0) > 50 && data.JetPt->at(1) > 50){
		// bkgd enriched
		// if( data.ZZMass > 160  && data.DiJetMass > 100 && (fabs(data.DiJetDEta) < 2.4 || data.DiJetMass < 400) && data.nExtraLep==0 && (((data.nCleanedJetsPt30==2||data.nCleanedJetsPt30==3)&&data.nCleanedJetsPt30BTagged_bTagSF<=1) ||(data.nCleanedJetsPt30>=4&&data.nCleanedJetsPt30BTagged_bTagSF==0)) ){
		  vbfcate=1;
			//}
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
		}
		tnew->Fill();

	}
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

