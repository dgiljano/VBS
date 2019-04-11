
/*
 * usage:
 * -specify parameters at the end of this file
 * -run with:
 *   root -l alpha.C++
 */

//Initial Include
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Math/DistFunc.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TF1.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TLine.h"
#include "TSpline.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TEfficiency.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "tdrstyle.C"
// #include "CMS_lumi.C"
#include "plotUtils.C"
#include "ZZAnalysis/AnalysisStep/src/kFactors.C"

#include <ZZAnalysis/AnalysisStep/src/Category.cc>
#include <ZZAnalysis/AnalysisStep/src/bitops.cc>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter/fit_functions.C>

using namespace std;
//Variable declaration set 1
int useHTBinned = 2;         // 0 - use simple DY inclusive
                             // 1 - use ht binned
                             // 2 - use jet binned + b-enricchement

bool enforceNarrowWidth = false;
bool unblind = true;
bool wantFinalLogPlots = true;
bool includeHiggsMass = true;

int onlyOneLep = 1;          // 0 - ee
                             // 1 - all leptons
                             // 2 - mumu

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//Variable declaration set 2
const int nVariables = 27;
string varName[nVariables] = {
  "ZZMass",
  "ZZPt",
  "ZZMassRefit",
  "Z1Mass",
  "Z2Mass",
  "Z1Pt",
  "Z2Pt",
  "Z2Flav",
  "Jet1Pt",
  "Jet2Pt",
  "JetQG",
  "JetBtagger",
  "Lep1Pt",
  "Lep2Pt",
  "AbsCosTheta1",
  "CosTheta2",
  "CosThetaStar",
  "PhiStar",
  "Phi",
  "NExtraJets",
  "MET",
  "Z1Tau21",
  "JetQGProduct",
  "vbsMELA",
  "vbsMELA2",
  "ZZMasshighMELA",
  "vbsMELA3"
};
string varXLabel[nVariables] = {
  "m_{2#font[12]{l}2q} (GeV)",
  "p_{T2#font[12]{l}2q} (GeV)",
  "m_{2#font[12]{l}2q} (GeV)",
  "m_{jj} (GeV)",
  "m_{#font[12]{l}#font[12]{l}} (GeV)",
  "p_{T,jj} (GeV)",
  "p_{T,#font[12]{l}#font[12]{l}} (GeV)",
  "flavor",
  "p_{Tj1} (GeV)",
  "p_{Tj2} (GeV)",
  "Quark gluon likelihood",
  "B tagger",
  "p_{T#font[12]{l}1} (GeV)",
  "p_{T#font[12]{l}2} (GeV)",
  "|cos(#theta_{1})|",
  "cos(#theta_{2})",
  "cos(#theta^{*})",
  "#Phi^{*}",
  "#Phi",
  "N_{extra-jets}",
  "MET (GeV)",
  "#tau_{21} (J)",
  "Quark gluon likelihood product",
  "vbsMELA",
  "vbsMELA2",
  "m_{2#font[12]{l}2q} (GeV) for vbsMELA > 0.5",
  "vbsMELA3"
};
string varYLabel[nVariables] = {
  "Events / 25 GeV",
  "Events / 10 GeV",
  "Events / 25 GeV",
  "Events / 2.5 GeV",
  "Events / 2.5 GeV",
  "Events / 10 GeV",
  "Events / 10 GeV",
  "Events",
  "Events / 10 GeV",
  "Events / 10 GeV",
  "Events / 0.028",
  "Events / 0.056",
  "Events / 20 GeV",
  "Events / 20 GeV",
  "Events",
  "Events",
  "Events",
  "Events",
  "Events",
  "Events",
  "Events / 6 GeV",
  "Events / 0.04",
  "Events / 0.025",
  "Events / 0.025",
  "Events / 0.025",
  "Events / 25 GeV",
  "Events / 0.025"
};
Int_t  varNbin[nVariables] = { 70, 50, 70,  56,  44, 50,50, 400,  50,  50,  50,  50,  50,  50,  50, 50, 50, 25, 25, 12, 50, 78, 50, 40, 40, 50, 40};
Float_t varMin[nVariables] = {  250,  0,  250,  40,  40,  90, 90, -200,  0, 0, -0.2, -0.2, 0,  0, -0.2, -1.2, -1.2, 0., 0., 0.5, 0., -0.05,-0.2,0.,0., 250,0.};
Float_t varMax[nVariables] = { 2000, 500, 2000, 180, 150, 800, 800, 0, 500, 500, 1.2, 1.2, 500, 500, 1.2, 1.2, 1.2 , 3.15, 3.15, 6.5, 300., 1.05, 1.2, 1.,1., 2000.,1.};
Bool_t varLogx[nVariables] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
Bool_t varLogy[nVariables] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0};

const int nMasses = 14;
string signalMasses[nMasses] = {"200","250","300","350","400","450","500","550","600","700","900","1000","1500","2000"};

enum Process {Data=0, ggSpin0=1, VBFSpin0=2, DYjets=3, TTBar=4, Diboson=5}; // Spin2=6};
const int nProcesses = 6;
string sProcess[nProcesses] = {"Data", "Spin0900", "Spin01500", "DY", "TT", "VV"}; // "Spin2800"};
string processLabel[nProcesses] = {"Data", "VBS WZ (x100)", "VBS ZZ (x100)", "Z + jets", "t#bar{t},WW", "ZZ, WZ"};

// WITH NNLO k-FACTOR FOR Z+JET
Float_t scaleF[nProcesses] = {1.,100.,100.,1.231,1.,1.12864505708};

const int nFS = 3;
string sFS[nFS] = {"ee","all","mm"};

const int nType = 12;
string typeS[nType] = {"resolvedSB","mergedSB","mergedSR","resolvedSR","resolvedSBbtag","mergedSBbtag","mergedSRbtag","resolvedSRbtag","resolvedSBvbf","mergedSBvbf","mergedSRvbf","resolvedSRvbf"};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//Variable declaration set 3 + various functions
bool passAdditionalCuts(float tmvaZ2Mass, float tmvaZ1Pt, float tmvaZ2Pt, float tmvaZ1tau21, float zzmass, bool merged, bool syncevents, int draw) {
  // if (syncevents) cout << "*** Enter cuts: mZ2 " << tmvaZ2Mass << endl;
  if (tmvaZ2Mass < 60.) return false;
  // if (syncevents) cout << "*** Passed mZ2 " << tmvaZ2Pt << " " << tmvaZ1Pt << endl;
  if (tmvaZ1Pt < 100. || tmvaZ2Pt < 100.) return false;   // TEST!
  // if (syncevents) cout << "*** Passed pTZ1 Z2 " << tmvaZ1tau21 << endl;
  if (merged && tmvaZ1tau21 > 0.6) return false;
  // if (syncevents) cout << "*** Passed Z1tau21 " << zzmass << endl;
  // if (zzmass < 500.) return false;
  if (zzmass < 400.) return false;
  // if (syncevents) cout << "*** Passed " << endl;
  return true;
}

TGraphAsymmErrors* doBkgEstGraph(int n, TH1F* central, TH1F* up, TH1F* down, bool isRatio, bool useUpDown = true) {
  float theX[n];
  float theY[n];
  float theEX[n];
  float theEYup[n];
  float theEYdown[n];
  for(int ibin=1; ibin<=central->GetNbinsX(); ibin++){
    theX[ibin-1] = central->GetXaxis()->GetBinCenter(ibin);
    theY[ibin-1] = central->GetBinContent(ibin);
    if (isRatio) theY[ibin-1] = 0.;
    theEX[ibin-1] = central->GetXaxis()->GetBinCenter(ibin) - central->GetXaxis()->GetBinLowEdge(ibin);
    theEYup[ibin-1] = up->GetBinContent(ibin)-central->GetBinContent(ibin);
    if (!useUpDown) theEYup[ibin-1] = central->GetBinError(ibin);
    if (isRatio) theEYup[ibin-1] /= central->GetBinContent(ibin);
    theEYdown[ibin-1] = central->GetBinContent(ibin)-down->GetBinContent(ibin);
    if (!useUpDown) theEYdown[ibin-1] = central->GetBinError(ibin);
    if (isRatio) theEYdown[ibin-1] /= central->GetBinContent(ibin);
  }
  TGraphAsymmErrors* tg = new TGraphAsymmErrors(n,theX,theY,theEX,theEX,theEYup,theEYdown);
  tg->SetLineColor(kBlue+2);
  tg->SetMarkerColor(kBlue+2);
  tg->SetFillColor(kBlue+2);
  if (!useUpDown) tg->SetLineColor(kGreen+2);
  if (!useUpDown) tg->SetMarkerColor(kGreen+2);
  if (!useUpDown) tg->SetFillColor(kGreen+2);
  tg->SetFillStyle(3004);
  tg->SetLineWidth(3);
  tg->SetMarkerStyle(25);
  return tg;
}

/* float DoBinCenter(TH1F* hmass, int ib) {
  float minimX = hmass->GetBinLowEdge(ib);
  float centrX = hmass->GetBinCenter(ib);
  if (centrX < 900.) return centrX-0.12*(centrX-minimX);
  else return centrX+0.12*(centrX-minimX);
  }

TF1* extender(TF1* basis, float xMin, float xMax, int type=0)
{
  TF1* func;
  if (type==0) func = new TF1("pippo",myfunction,xMin,xMax,4);
  else if (type==1) func = new TF1("pippo",myfunctionErrUp,xMin,xMax,20);
  else func = new TF1("pippo",myfunctionErrDown,xMin,xMax,20);
  func->SetName(basis->GetName());
  func->SetLineColor(basis->GetLineColor());
  func->SetLineStyle(basis->GetLineStyle());
  func->SetLineWidth(basis->GetLineWidth());
  for (int i=0; i<basis->GetNpar(); i++) {
    func->SetParameter(i,basis->GetParameter(i));
  }
  return func;
  }
*/

float deltaPhi(float phi1, float phi2)
{
  float delt = phi1-phi2;
  while (delt >= 3.1416) delt -= 6.2832;
  while (delt < -3.1416) delt += 6.2832;
  return delt;
}

void densityHist(TH1F* hist)
{
  for (int i=1; i<=hist->GetNbinsX(); i++) {
    float binsize = hist->GetXaxis()->GetBinUpEdge(i) - hist->GetXaxis()->GetBinLowEdge(i);
    hist->SetBinContent(i,hist->GetBinContent(i)*50./binsize);
    hist->SetBinError(i,hist->GetBinError(i)*50./binsize);
  }
  hist->GetXaxis()->SetTitle("m_{ZZ} (GeV)");
  hist->GetYaxis()->SetTitle("Events / 50 GeV");
}

float getDVBF2jetsConstant(float ZZMass){
  float par[9]={
    1.876,
    -55.488,
    403.32,
    0.3906,
    80.8,
    27.7,
    -0.06,
    54.97,
    309.96
  };
  float kappa =
    pow(1.-atan((ZZMass-par[1])/par[2])*2./TMath::Pi(), par[0])
    + par[3]*exp(-pow((ZZMass-par[4])/par[5], 2))
    + par[6]*exp(-pow((ZZMass-par[7])/par[8], 2));
  float constant = kappa/(1.-kappa);
  return constant;
}

float getDZjjspin0Constant(float ZZMass){
  float constant = 0.0035*(3.05+(0.0005*ZZMass)-(2.73/(1.+exp((ZZMass-2500.)/258.26)))) ;
  return constant;
}

float getDZjjspin2Constant(float ZZMass){
  return 0.14;
}

float getDVBSConstant(float ZZMass){
  return 0.014;
}

void alpha(string dirout = "test13TeV", string theNtupleFile = "s.txt")
//void plotDataVsMC_2l2q(string dirout = "test13TeV", string theNtupleFile = "goodDatasets_oldProd_VBS.txt")
{
  bool sync = false;
  bool CR = false;
  int draw = 2;
  bool weightMCtrig = true;
  bool weighttau21 = false;

  // draw = -1: do not draw, but save inputs for fits (fine bins)
  //      = 0 : do not draw, but save inputs for fits (large bins)
  //      = 1 : draw only PAS-style plots
  //      = 2 : draw all
  // draw = 0;
//Variable delcaration set 4 (including lumi)
  float lumin = 35.8;   // ICHEP total
  setTDRStyle();

  const int nDatasets = 53;          // Moriond: 11
  const int nDatasetsMC = 13;         // Moriond: 9

  TFile* inputFile[nDatasets];
  TChain* inputTree[nDatasets];
  TH1F* hCounters[nDatasets];
  Long64_t NGenEvt[nDatasets];
  Float_t NEvtNarrow[nDatasets];
  Float_t sumWeights[nDatasets];
  Int_t mass[nDatasets];

  string Dsname[nDatasets] = {"VBSZWminus_mjj100","VBSZWplus_mjj100","VBSZZ_mjj100"/* ,"DY0Jet" */, "DY1Jet","DY2Jet","DY3Jet","DY4Jet","DYBJet","DYBFiltJet","TTBar","WW2l2nDib","WZDib","ZZDib","DoubleEG2016Bv1","DoubleMu2016Bv1","DoubleMu2016Bv2","DoubleEG2016Bv2","DoubleEG2016C","DoubleMu2016C","DoubleEG2016D","DoubleMu2016D","DoubleEG2016E","DoubleMu2016E","DoubleEG2016F","DoubleMu2016F","DoubleEG2016G","DoubleMu2016G","DoubleEG2016Hv1","DoubleMu2016Hv1","DoubleEG2016Hv2","DoubleMu2016Hv2","DoubleEG2016Hv3","DoubleMu2016Hv3","SingleEl2016Bv1","SingleMu2016Bv1","SingleEl2016Bv2","SingleMu2016Bv2","SingleEl2016C","SingleMu2016C","SingleEl2016D","SingleMu2016D","SingleEl2016E","SingleMu2016E","SingleEl2016F","SingleMu2016F","SingleEl2016G","SingleMu2016G","SingleEl2016Hv1","SingleMu2016Hv1","SingleEl2016Hv2","SingleMu2016Hv2","SingleEl2016Hv3","SingleMu2016Hv3"};

  // ZPT reweighting
  TFile* fKzpt = new TFile("ZPTCorrection/Kfac_pTZ.root");
  TH1D* NLO = (TH1D*)fKzpt->Get("NLO");// This NLO is 2015 measured results
  TH1D* LO = (TH1D*)fKzpt->Get("LO");// This LO is the LO*1.238 where 1.238 in inclusive NLO/LO K factor
  TH1D* Kzpt_ratio = (TH1D*)NLO->Clone("Kzpt_ratio");
  Kzpt_ratio->Divide(LO);
  TGraphErrors* gKzpt_ratio = new TGraphErrors(Kzpt_ratio);

  // EWK NOT AVAILABLE
  /* TFile* fewk = new TFile("kfactor.root");
  TH1D* hEWK = (TH1D*)fewk->Get("kfactor_ewk");
  TGraphErrors* gKewk = new TGraphErrors(hEWK);  */

  if (useHTBinned == 0) {
    Dsname[2] = "DYJetsToLL";
  } else if (useHTBinned == 1) {
    Dsname[2] = "DYHT100";
    Dsname[3] = "DYHT200";
    Dsname[4] = "DYHT400";
    Dsname[5] = "DYHT600";
    processLabel[3] = "Z + jets (HT > 100 GeV)";
  }

  // trigger weights
  TFile *f_data_se  = new TFile("trigeff/eff_data_2d_se.root","READ");
  TFile *f_data_de1 = new TFile("trigeff/eff_data_2d_de1.root","READ");
  TFile *f_data_de2 = new TFile("trigeff/eff_data_2d_de2.root","READ");
  TFile *f_data_sm  = new TFile("trigeff/eff_data_2d_sm.root","READ");
  TFile *f_data_dm1 = new TFile("trigeff/eff_data_2d_dm1.root","READ");
  TFile *f_data_dm2 = new TFile("trigeff/eff_data_2d_dm2.root","READ");

  TCanvas *c_data_se = (TCanvas*)f_data_se->Get("c1");
  TEfficiency *eff_data_se = (TEfficiency*)c_data_se->GetPrimitive("den_se_2d_clone");
  TCanvas *c_data_de1 = (TCanvas*)f_data_de1->Get("c1");
  TEfficiency *eff_data_de1 = (TEfficiency*)c_data_de1->GetPrimitive("den_de1_2d_clone");
  TCanvas *c_data_de2 = (TCanvas*)f_data_de2->Get("c1");
  TEfficiency *eff_data_de2 = (TEfficiency*)c_data_de2->GetPrimitive("den_de2_2d_clone");
  TCanvas *c_data_sm = (TCanvas*)f_data_sm->Get("c1");
  TEfficiency *eff_data_sm = (TEfficiency*)c_data_sm->GetPrimitive("den_sm_2d_clone");
  TCanvas *c_data_dm1 = (TCanvas*)f_data_dm1->Get("c1");
  TEfficiency *eff_data_dm1 = (TEfficiency*)c_data_dm1->GetPrimitive("den_dm1_2d_clone");
  TCanvas *c_data_dm2 = (TCanvas*)f_data_dm2->Get("c1");
  TEfficiency *eff_data_dm2 = (TEfficiency*)c_data_dm2->GetPrimitive("den_dm2_2d_clone");

  TFile *f_mc_se  = new TFile("trigeff/eff_mc_2d_se.root","READ");
  TFile *f_mc_de1 = new TFile("trigeff/eff_mc_2d_de1.root","READ");
  TFile *f_mc_de2 = new TFile("trigeff/eff_mc_2d_de2.root","READ");
  TFile *f_mc_sm  = new TFile("trigeff/eff_mc_2d_sm.root","READ");
  TFile *f_mc_dm1 = new TFile("trigeff/eff_mc_2d_dm1.root","READ");
  TFile *f_mc_dm2 = new TFile("trigeff/eff_mc_2d_dm2.root","READ");

  TCanvas *c_mc_se = (TCanvas*)f_mc_se->Get("c2");
  TEfficiency *eff_mc_se = (TEfficiency*)c_mc_se->GetPrimitive("den_se_2d_clone");
  TCanvas *c_mc_de1 = (TCanvas*)f_mc_de1->Get("c2");
  TEfficiency *eff_mc_de1 = (TEfficiency*)c_mc_de1->GetPrimitive("den_de1_2d_clone");
  TCanvas *c_mc_de2 = (TCanvas*)f_mc_de2->Get("c2");
  TEfficiency *eff_mc_de2 = (TEfficiency*)c_mc_de2->GetPrimitive("den_de2_2d_clone");
  TCanvas *c_mc_sm = (TCanvas*)f_mc_sm->Get("c2");
  TEfficiency *eff_mc_sm = (TEfficiency*)c_mc_sm->GetPrimitive("den_sm_2d_clone");
  TCanvas *c_mc_dm1 = (TCanvas*)f_mc_dm1->Get("c2");
  TEfficiency *eff_mc_dm1 = (TEfficiency*)c_mc_dm1->GetPrimitive("den_dm1_2d_clone");
  TCanvas *c_mc_dm2 = (TCanvas*)f_mc_dm2->Get("c2");
  TEfficiency *eff_mc_dm2 = (TEfficiency*)c_mc_dm2->GetPrimitive("den_dm2_2d_clone");
  // end trigger weights

  float tau21bin[39] = {
    -0.05,
    -0.021794,
    0.0064102,
    0.0346154,
    0.0628205,
    0.0910256,
    0.119231,
    0.147436,
    0.175641,
    0.203846,
    0.232051,
    0.260256,
    0.288462,
    0.316667,
    0.344872,
    0.373077,
    0.401282,
    0.429487,
    0.457692,
    0.485897,
    0.514103,
    0.542308,
    0.570513,
    0.598718,
    0.626923,
    0.655128,
    0.683333,
    0.711538,
    0.739744,
    0.767949,
    0.796154,
    0.824359,
    0.852564,
    0.880769,
    0.908974,
    0.937179,
    0.965385,
    0.99359,
    1.02179
  };

  float tau21corr[38] = {
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.853534,
    -0.488688,
    0.220111,
    -0.416235,
    0.322711,
    0.283294,
    -0.0625665,
    -0.217968,
    -0.239136,
    0.132369,
    -0.0898708,
    -0.132034,
    -0.155736,
    -0.0835754,
    -0.106685,
    -0.0437803,
    -0.424335,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
  };


  /// I/O to TMVA
  float tmvaZZPt, tmvaZ2Mass, tmvaZ1Pt, tmvaZ1tau21, tmvaZ2Pt, tmvaLepPt1, tmvaLepPt2, tmvaJetPt1, tmvaJetPt2;
  float tmvaJetQGLikelihood1, tmvaJetQGLikelihood2, tmvaabshelcosthetaZ1, tmvahelcosthetaZ2, tmvacosthetastar, tmvahelphi, tmvaphistarZ1;
//2l2q PRELIMINARY TREE MANAGEMENT

  //String name
  TString outfileName_test( "RooFitInput.root" );

  //Output file creation
  TFile* outputFile = TFile::Open( outfileName_test, "RECREATE" );

  //Tree creation
  outputFile->cd();
  TTree *tnew = new TTree("SelectedTree","SelectedTree");

  //Variable declaration
  float weight, weight_up, weight_dn;
  float weight_vbf, weight_vbf_up, weight_vbf_dn;
  int chan = 0;
  int vbfcate = 0;
  float dbkg_kin = 0;
  float ZZMass1 = 0;

  //other declarations
  //output file template binning
                int nbins=32;
                double xbin[33]={
                                 110,120,130,140,150,160,170,180,190,
                                 200,210,220,230,240,250,270,290,310,
                                 330,350,370,390,420,460,500,540,580,
                                 620,660,700,750,800,3500
                                };


  //Branch declaration
  tnew->Branch("mreco",&ZZMass1,"mreco/F");
  tnew->Branch("dbkg_kin",&dbkg_kin,"dbkg_kin/F");
  tnew->Branch("chan",&chan,"chan/I");
  tnew->Branch("vbfcate",&vbfcate,"vbfcate/I");



//  TH2F *temp_mm_mergedSRvbf = new TH2F("temp_mm_mergedSRvbf","",nbins,xbin,30,0.,1.);
  //2D histogram declarations
  TH2F* temp[nFS][nType];
  TH2F* temp_rebin[nFS][nType];

  for(int rs=0; rs<nFS; rs++){  					 //ee, mumu, or all

	for(int nt=0; nt<nType; nt++){                                     //resolved,merged,SB,SR

		temp[rs][nt] = new TH2F(Form("h1_%s_%s",sFS[rs].c_str(),typeS[nt].c_str()),"",nbins,xbin,30,0.,1.);
		temp_rebin[rs][nt] = new TH2F(Form("temp_rebin_%s_%s",sFS[rs].c_str(),typeS[nt].c_str()),"",1695,110,3500,30,0.,1.);
	 }
  }

  //Tree Writing
  //tnew->Write();

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Histogram declaration and initialisation
 // TString outfileName( "RooFitInput.root" );
 // if (onlyOneLep == 0) outfileName = "TMVAAndRoofitInputs_ee.root";
 // if (onlyOneLep == 2) outfileName = "TMVAAndRoofitInputs_mumu.root";

 // TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  // ofstream mengch;
  // mengch.open("mengch.txt");

  char filestring[400];
  float xMin[nType],xMax[nType];

  Int_t RunNumber;
  Long64_t EventNumber;
  Int_t LumiNumber;
  Float_t genEventWeight;
  Float_t overallEventWeight;
  Float_t Nvtx;
  Float_t genHMass;
  vector<Short_t> *ZZsel = 0;
  vector<Float_t> *ZZMass = 0;
  vector<Float_t> *ZZMass_up = 0;
  vector<Float_t> *ZZMass_dn = 0;
  vector<Float_t> *ZZMassRefit = 0;
  vector<Float_t> *ZZPt = 0;
  vector<Float_t> *Z1Mass = 0;
  vector<Float_t> *Z2Mass = 0;
  vector<Float_t> *Z1Pt = 0;
  vector<Float_t> *Z1tau21 = 0;
  vector<Float_t> *Z2Pt = 0;
  vector<Short_t> *Z2Flav = 0;
  vector<Short_t> *ZZCandType = 0;
  vector<Float_t> *LepPt = 0;
  vector<Float_t> *LepEta = 0;
  vector<Float_t> *LepPhi = 0;
  Short_t nJets;
  Short_t nJetsBTagged;
  vector<Float_t> *JetPt = 0;
  vector<Float_t> *JetPhi = 0;
  vector<Float_t> *JetEta = 0;
  vector<bool> *JetIsInZZCand = 0;
  vector<Float_t> *JetBTagger = 0;
  vector<Float_t> *JetIsBtaggedWithSF = 0;
  vector<Float_t> *JetQGLikelihood = 0;
  vector<Float_t> *helcosthetaZ1 = 0;
  vector<Float_t> *helcosthetaZ2 = 0;
  vector<Float_t> *costhetastar = 0;
  vector<Float_t> *helphi = 0;
  vector<Float_t> *phistarZ1 = 0;
  Float_t DiJetMass;
  vector<Float_t> *p0plus_VAJHU = 0;
  vector<Float_t> *p2bplus_VAJHU = 0;
  vector<Float_t> *pqqZJJ_VAMCFM = 0;
  vector<Float_t> *pvbs_VAMCFM_highestPTJets = 0;
  vector<Float_t> *pzzjj_VAMCFM_highestPTJets = 0;
  vector<Float_t> *pvbf_VAJHU_highestPTJets = 0;
  vector<Float_t> *phjj_VAJHU_highestPTJets = 0;
  Float_t xsec;
  Float_t Met;
  Float_t GenLep1Pt;
  Float_t GenLep1Phi;
  Float_t GenLep2Pt;
  Float_t GenLep2Phi;

  TH1F* h1[nVariables][nProcesses][nFS][nType];
  TH1F* h_alpha[nVariables][nProcesses][nFS][nType];
  for(int rs=0; rs<nFS; rs++){   //ee, mumu, or all
    for(int pr=0; pr<nProcesses; pr++){
      for(int nt=0; nt<nType; nt++){
	for(int v=0; v<nVariables; v++){
	  if (nt == 0 || ( nt == 3 && (v!=0 && v!=2)) )
	    h1[v][pr][rs][nt] = new TH1F(Form("h1_%s_%s_%s_%s",varName[v].c_str(),sFS[rs].c_str(),typeS[nt].c_str(),sProcess[pr].c_str()),
					 Form("h1_%s_%s_%s_%s",varName[v].c_str(),sFS[rs].c_str(),typeS[nt].c_str(),sProcess[pr].c_str()),
					 varNbin[v],varMin[v],varMax[v]);
	  else
	    h1[v][pr][rs][nt] = new TH1F(Form("h1_%s_%s_%s_%s",varName[v].c_str(),sFS[rs].c_str(),typeS[nt].c_str(),sProcess[pr].c_str()),
					 Form("h1_%s_%s_%s_%s",varName[v].c_str(),sFS[rs].c_str(),typeS[nt].c_str(),sProcess[pr].c_str()),
					 varNbin[v]/2,varMin[v],varMax[v]);
          h1[v][pr][rs][nt]->Sumw2();
	}
      }
    }
  }

  TH1F* hmass[nProcesses][nType+7];
  TH1F* hmass_up[nProcesses][nType+4];
  TH1F* hmass_down[nProcesses][nType+4];

  Float_t binsincl[] = { 500, 510, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1300, 1400, 1500, 1600, 2000, 2400, 3000};
  Int_t  binnumincl = sizeof(binsincl)/sizeof(Float_t) - 1;
  Float_t binsbtag[] = { 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1400, 1600, 2000, 2400, 3000};
  Int_t  binnumbtag = sizeof(binsbtag)/sizeof(Float_t) - 1;
  TString hmasstags[nType+7] = {"","","untagged, merged jet","untagged, resolved jets",
				"","","b-tagged, merged jet","b-tagged, resolved jets",
				"","","VBF-tagged, merged jet","VBF-tagged, resolved jets",
				"","","","","all categories, resolved jets","all categories, resolved jets","all categories, resolved jets"};

  for(int pr=0; pr<nProcesses; pr++){
    for(int nt=0; nt<nType+7; nt++){
      if (nt == 0 || nt == 3) {        // mZZ high stats
	if (draw >= 0) hmass[pr][nt] = new TH1F(Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
					       Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
					       binnumincl,binsincl);
	else hmass[pr][nt] = new TH1F(Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
				      Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
				      600, 500, 3500);
	if (draw >= 0) hmass_up[pr][nt] = new TH1F(Form("hmass_up_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
						  Form("hmass_up_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
						  binnumincl,binsincl);
	else hmass_up[pr][nt] = new TH1F(Form("hmass_up_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
					 Form("hmass_up_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
					 600, 500, 3500);
	if (draw >= 0) hmass_down[pr][nt] = new TH1F(Form("hmass_down_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
						    Form("hmass_down_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
						    binnumincl,binsincl);
	else  hmass_down[pr][nt] = new TH1F(Form("hmass_down_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
					    Form("hmass_down_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
					    600,500,3500);
      } else if (nt<nType) { // mZZ low stats
	 if (draw >= 0)  hmass[pr][nt] = new TH1F(Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
						 Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
						 binnumbtag,binsbtag);
	 else hmass[pr][nt] = new TH1F(Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
				       Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
				       600, 500, 3500);
	 if (draw >= 0) hmass_up[pr][nt] = new TH1F(Form("hmass_up_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
						    Form("hmass_up_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
						    binnumbtag,binsbtag);
	 else hmass_up[pr][nt] = new TH1F(Form("hmass_up_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
					  Form("hmass_up_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
					  600, 500, 3500);
	 if (draw >= 0 ) hmass_down[pr][nt] = new TH1F(Form("hmass_down_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
						      Form("hmass_down_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
						      binnumbtag,binsbtag);
	 else  hmass_down[pr][nt] = new TH1F(Form("hmass_down_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
					     Form("hmass_down_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
					     600,500,3500);
      } else if (nt<nType+4) { // mZZ non b-tagged
	hmass[pr][nt] = new TH1F(Form("hmass_%snobtag_%s",typeS[nt-nType].c_str(),sProcess[pr].c_str()),
				 Form("hmass_%snobtag_%s",typeS[nt-nType].c_str(),sProcess[pr].c_str()),
				 binnumincl,binsincl);
        hmass_up[pr][nt] = new TH1F(Form("hmass_up_%snobtag_%s",typeS[nt-nType].c_str(),sProcess[pr].c_str()),
				    Form("hmass_up_%snobtag_%s",typeS[nt-nType].c_str(),sProcess[pr].c_str()),
				    binnumincl,binsincl);
	hmass_down[pr][nt] = new TH1F(Form("hmass_down_%snobtag_%s",typeS[nt-nType].c_str(),sProcess[pr].c_str()),
				      Form("hmass_down_%snobtag_%s",typeS[nt-nType].c_str(),sProcess[pr].c_str()),
				      binnumincl,binsincl);
      } else if (nt == nType+4)          // MELA spin-0
	hmass[pr][nt] = new TH1F(Form("hmelaspin0_resolvedSR_%s",sProcess[pr].c_str()),
				 Form("hmelaspin0_resolvedSR_%s",sProcess[pr].c_str()),
				 20,0.,1.);
      else if (nt == nType+5)         // MELA spin-2
	hmass[pr][nt] = new TH1F(Form("hmelaspin2_resolvedSR_%s",sProcess[pr].c_str()),
				 Form("hmelaspin2_resolvedSR_%s",sProcess[pr].c_str()),
				 20,0.,1.);
      else if (nt == nType+6)         // VBF-tagged
	hmass[pr][nt] = new TH1F(Form("hmelavbf_allSR_%s",sProcess[pr].c_str()),
				 Form("hmelavbf_allSR_%s",sProcess[pr].c_str()),
				 42,-1.05,1.05);

      hmass[pr][nt]->Sumw2();
    }
  }

  // background functions (if unblind)
  // TFile ffunc("results.root","READ");
  TF1* ffit[nType];
  TF1* ffitup[nType];
  TF1* ffitdown[nType];
  // TF1* ffit_temp[nType];
  // TF1* ffitup_temp[nType];
  // TF1* ffitdown_temp[nType];

  TH1F* hbkg[nType];
  TH1F* hbkg_up[nType];
  TH1F* hbkg_down[nType];

  for(int nt=0; nt<nType; nt++){
    hbkg[nt] = (TH1F*)hmass[0][nt]->Clone();
    hbkg_up[nt] = (TH1F*)hmass[0][nt]->Clone();
    hbkg_down[nt] = (TH1F*)hmass[0][nt]->Clone();
  }

  if (unblind) {

    TFile f1("/afs/cern.ch/work/x/xiaomeng/public/mlfit.root");
    for(int nt=0; nt<nType; nt++){
      // ifstream parInput;
      // float parValue[20];  string parName;
      // if (string(typeS[nt]).find("merged") != std::string::npos) xMin[nt] = 700.;
      xMin[nt] = 500.;
      xMax[nt] = 3500.;

      if (nt == 2 || nt == 3 || nt == 6 || nt == 10 || nt == 7 || nt == 11) {   // MODIFY HERE
        unsigned int icode, jcode;   // channel code in Combine output
        if (nt == 2) { icode = 2; jcode = 1; }
        if (nt == 3) { icode = 1; jcode = 1; }
        if (nt == 6) { icode = 2; jcode = 2; }
        if (nt == 7) { icode = 1; jcode = 2; }
        if (nt == 10) { icode = 2; jcode = 3; }
        if (nt == 11) { icode = 1; jcode = 3; }

	Float_t binarray[binnumincl+1];
        if (nt == 3) {
	  for (int aaa=0; aaa<binnumincl;aaa++) {
	    binarray[aaa]=binsincl[aaa];
	  }
	  binarray[binnumincl]=3500.;
	} else {
	  for (int aaa=0; aaa<binnumbtag;aaa++) {
	    binarray[aaa]=binsbtag[aaa];
	  }
          binarray[binnumbtag]=3500.;
	}



        // sprintf(filestring,"shapes_fit_b/ch%d_ch1_ch%d/total_background",icode,jcode);
        sprintf(filestring,"shapes_prefit/ch%d_ch1_ch%d/total_background",icode,jcode);
	TH1F* htemp = (TH1F*)f1.Get(filestring);
	// sprintf(filestring,"shapes_fit_b/ch%d_ch2_ch%d/total_background",icode,jcode);
        sprintf(filestring,"shapes_prefit/ch%d_ch2_ch%d/total_background",icode,jcode);
	TH1F* htemp2 = (TH1F*)f1.Get(filestring);
	cout << "Reading fit file " << nt << endl;
	// rebin before density
	int ipast = 1;  float summa = 0.;  float summaErr = 0.;
	for(int ibin=1; ibin<=htemp->GetNbinsX(); ibin++){
          if (htemp->GetXaxis()->GetBinCenter(ibin) < xMin[nt]) continue;
	  if (htemp->GetXaxis()->GetBinCenter(ibin) > xMax[nt]) break;
	  if (htemp->GetXaxis()->GetBinCenter(ibin) < binarray[ipast]) {
	    summa += 10.*htemp->GetBinContent(ibin);
	    summaErr += 100.*htemp->GetBinError(ibin)*htemp->GetBinError(ibin);
	    summa += 10.*htemp2->GetBinContent(ibin);
	    summaErr += 100.*htemp2->GetBinError(ibin)*htemp2->GetBinError(ibin);
	    if (htemp->GetXaxis()->GetBinCenter(ibin+1) > binarray[ipast]) {
	      hbkg[nt]->SetBinContent(ipast,summa);
	      hbkg_up[nt]->SetBinContent(ipast,summa+sqrt(summaErr));
	      hbkg_down[nt]->SetBinContent(ipast,summa-sqrt(summaErr));
	      ipast++;         summa = 0.;      summaErr = 0.;
	    }
	  }
	}
        if (nt==3) hbkg[nt]->Scale(0.95);
	if (nt==3) hbkg_up[nt]->Scale(0.95);
	if (nt==3) hbkg_down[nt]->Scale(0.95);
	densityHist(hbkg[nt]);
	densityHist(hbkg_up[nt]);
	densityHist(hbkg_down[nt]);
      }
    }


  }
//Input data management and preliminary tree read settings
  // end background functions

  //---------- Will loop over all datasets
  for (int d=0; d<nDatasets; d++) {
    NGenEvt[d] = 0;      NEvtNarrow[d] = 0.;
    sumWeights[d] = 0.;  mass[d] = 0;
    inputTree[d] = new TChain("ZZTree/candTree");
  }
std::cout << "TEST" << endl;
  ifstream list(theNtupleFile.c_str());
  char fileName[400];
  while (list >> fileName) {
    if (string(fileName).find("store") != std::string::npos) sprintf(filestring,"root://eoscms//eos/cms/%s",fileName);
    sprintf(filestring,"%s",fileName);
    std::cout<<fileName << endl;
    TFile* ftemp = TFile::Open(filestring);
    cout << filestring << endl;
    for (int d=0; d<nDatasets; d++) {
      if (string(fileName).find(Dsname[d].c_str()) != std::string::npos) {
	inputTree[d]->Add(filestring);
	if (d<nDatasetsMC) {
	  hCounters[d] = (TH1F*)ftemp->Get("ZZTree/Counters");
	  NGenEvt[d] += hCounters[d]->GetBinContent(1);
          sumWeights[d] += hCounters[d]->GetBinContent(40);  // all weights
	}
	if (d<3) {
	  for (int m=0; m<nMasses; m++) {
	    if (string(fileName).find(signalMasses[m].c_str()) != std::string::npos) mass[d] = atoi(signalMasses[m].c_str());
	  }
	}
      }
    }
    ftemp->Close();
    if (string(fileName).find("160601") != std::string::npos) {
      scaleF[1] = 30.;      scaleF[2] = 30.;   // forgot 2tau2q
    }
  }


  for (int d=0; d<nDatasets; d++) {

    if (useHTBinned == 0 && d>2 && d<8) continue;   // in this case there is just one DY
    if (useHTBinned == 1 && d>5 && d<8) continue;   // in this case there are just four DY

    inputTree[d]->SetBranchAddress("Nvtx", &Nvtx);
    inputTree[d]->SetBranchAddress("RunNumber", &RunNumber);
    inputTree[d]->SetBranchAddress("EventNumber", &EventNumber);
    inputTree[d]->SetBranchAddress("LumiNumber", &LumiNumber);
    inputTree[d]->SetBranchAddress("genHEPMCweight", &genEventWeight);
    inputTree[d]->SetBranchAddress("overallEventWeight", &overallEventWeight);
    inputTree[d]->SetBranchAddress("ZZsel", &ZZsel);
    inputTree[d]->SetBranchAddress("ZZMass", &ZZMass);
    inputTree[d]->SetBranchAddress("ZZMass_up", &ZZMass_up);
    inputTree[d]->SetBranchAddress("ZZMass_dn", &ZZMass_dn);
    inputTree[d]->SetBranchAddress("ZZMassRefit", &ZZMassRefit);
    inputTree[d]->SetBranchAddress("ZZPt", &ZZPt);
    inputTree[d]->SetBranchAddress("Z1Mass", &Z1Mass);
    inputTree[d]->SetBranchAddress("Z1tau21", &Z1tau21);
    inputTree[d]->SetBranchAddress("Z2Mass", &Z2Mass);
    inputTree[d]->SetBranchAddress("Z1Pt", &Z1Pt);
    inputTree[d]->SetBranchAddress("Z2Pt", &Z2Pt);
    inputTree[d]->SetBranchAddress("Z2Flav", &Z2Flav);
    inputTree[d]->SetBranchAddress("ZZCandType", &ZZCandType);
    inputTree[d]->SetBranchAddress("LepPt", &LepPt);
    inputTree[d]->SetBranchAddress("LepPhi", &LepPhi);
    inputTree[d]->SetBranchAddress("LepEta", &LepEta);
    inputTree[d]->SetBranchAddress("JetPt", &JetPt);
    inputTree[d]->SetBranchAddress("JetEta", &JetEta);
    inputTree[d]->SetBranchAddress("JetPhi", &JetPhi);
    inputTree[d]->SetBranchAddress("JetIsInZZCand", &JetIsInZZCand);
    inputTree[d]->SetBranchAddress("JetBTagger", &JetBTagger);
    inputTree[d]->SetBranchAddress("JetIsBtaggedWithSF", &JetIsBtaggedWithSF);
    inputTree[d]->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood);
    inputTree[d]->SetBranchAddress("helcosthetaZ1",&helcosthetaZ1);
    inputTree[d]->SetBranchAddress("helcosthetaZ2",&helcosthetaZ2);
    inputTree[d]->SetBranchAddress("costhetastar",&costhetastar);
    inputTree[d]->SetBranchAddress("helphi",	  &helphi);
    inputTree[d]->SetBranchAddress("phistarZ1",   &phistarZ1);
    inputTree[d]->SetBranchAddress("DiJetMass",   &DiJetMass);
    inputTree[d]->SetBranchAddress("xsec", &xsec);
    inputTree[d]->SetBranchAddress("GenHMass", &genHMass);
    inputTree[d]->SetBranchAddress("PFMET", &Met);
    inputTree[d]->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU );
    inputTree[d]->SetBranchAddress("p2bplus_VAJHU", &p2bplus_VAJHU );
    inputTree[d]->SetBranchAddress("pqqZJJ_VAMCFM", &pqqZJJ_VAMCFM );
    inputTree[d]->SetBranchAddress("pvbs_VAMCFM_highestPTJets", &pvbs_VAMCFM_highestPTJets );
    inputTree[d]->SetBranchAddress("pzzjj_VAMCFM_highestPTJets", &pzzjj_VAMCFM_highestPTJets );
    inputTree[d]->SetBranchAddress("pvbf_VAJHU_highestPTJets", &pvbf_VAJHU_highestPTJets );
    inputTree[d]->SetBranchAddress("phjj_VAJHU_highestPTJets", &phjj_VAJHU_highestPTJets );
    inputTree[d]->SetBranchAddress("GenLep1Pt", &GenLep1Pt);
    inputTree[d]->SetBranchAddress("GenLep1Phi", &GenLep1Phi);
    inputTree[d]->SetBranchAddress("GenLep2Pt", &GenLep2Pt);
    inputTree[d]->SetBranchAddress("GenLep2Phi", &GenLep2Phi);

//Process tree and hist fill 1
    //---------- Process tree

    Long64_t entries = inputTree[d]->GetEntries();
 // categorize datasets into process
    int process;
    if (d==0 || d==1) process=1;                  //signal WZ
    else if (d==2) process=2;             //signal ZZ
    else if (d>2 && d<9) process=3;       //DY
    else if (d==9 || d==10) process=4;     //TT+WW
    else if (d>=11 && d<13) process=5;    //WZ,ZZ
    else process=0;                       //data

    // for synchronization
    ofstream myfile;
    int nPassMerged = 0 ;
    int nPassResol = 0 ;
    if (process == 0 && sync) myfile.open("synchronization.txt");

    if (process>0) {
    float eff = float(entries)/float(NGenEvt[d]);
    cout<<"Processing dataset "<<Dsname[d]<<" ("<<entries<<" entries of "<< NGenEvt[d] << " = " << eff*100. << "%) ..."<<endl;
    cout<<"process "<<process<<" d "<<d<<" "<<Dsname[d]<<endl;

    } else {
      cout<<"Processing dataset "<<d<<" ("<<entries<<" entries)"<<endl;
    }

    for (Long64_t z=0; z<entries; ++z){

      //if (z%10000 == 0) cout<<"Processing entry "<<z<<endl;

      inputTree[d]->GetEntry(z);

      bool syncevents = process == 0 && ((RunNumber==273158 && EventNumber==402295935 && LumiNumber==256) ||
					(RunNumber==274200 && EventNumber==806347507 && LumiNumber==499) ||
					(RunNumber==278167 && EventNumber==1293600593  && LumiNumber==723));

      if (ZZsel->size() < 1) continue;

      bool writeThis = (process==0 && sync);

      Double_t eventWeight = 1. ;

      if (process > 0) eventWeight = ( lumin * 1000 * scaleF[process] * overallEventWeight  / sumWeights[d] ) * xsec ;
      if (process > 0 && z == 0) cout << "cross-section = " << xsec << " pb; eventweight = " << eventWeight << endl;

      /// Apply extra kFactors: LO to NLO in QCD, plus EWK corrections for di-boson
      /*
      if(process==5 && string(Dsname[d]).find("ZZ") != std::string::npos) {
         double Kzpt = gKzpt_ratio->Eval(ZPt);
         eventWeight *= Kewk;
      }
      */
      if(string(Dsname[d]).find("DY") != std::string::npos) //DY but only for LO jet-binned
      {
          double ZPt = sqrt(GenLep1Pt*GenLep1Pt + GenLep2Pt*GenLep2Pt + 2*GenLep1Pt*GenLep2Pt*cos(GenLep1Phi-GenLep2Phi));

          //cout<<"ZPT reweighting comparing GEN zpt "<<ZPt<<" and RECO zpt "<<double(Z2Pt->at(0))<<endl;
          double Kzpt = gKzpt_ratio->Eval(double(Z2Pt->at(0)));
	  eventWeight = eventWeight * Kzpt;

      }

      // Apply extra kFactors: LO to NLO in QCD, plus EWK corrections for di-boson
      if(string(Dsname[d]).find("TT") != std::string::npos) {
          eventWeight = eventWeight * 87.31/57.35;
      }

      // keep only events around nominal mass!!!
      if (enforceNarrowWidth && mass[d] > 0 && fabs(genHMass-(float)mass[d]) > 0.05*mass[d]) continue;
      NEvtNarrow[d] = NEvtNarrow[d] + 1.;

      // find leading jets (notice this vector also includes subjets (identified by a btagger value of -1) which must be treated apart
      float pt1stJet = 0.0001;
      float pt2ndJet = 0.0001;
      float btag1stJet = 0.;
      float btag2ndJet = 0.;
      bool isB1stJet = false;
      bool isB2ndJet = false;
      float qglik1stJet = 0.;
      float qglik2ndJet = 0.;

      float pt1stSubjet = 0.0001;
      float pt2ndSubjet = 0.0001;
      float btag1stSubjet = 0.;
      float btag2ndSubjet = 0.;

      int nInJets = 0;
      int nExtraJets = 0;

      for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
	if (JetQGLikelihood->at(nJet) > -800.) {            // real jets
 	  if (JetIsInZZCand->at(nJet) ) {
	    if (pt1stJet < JetPt->at(nJet)) {
	      pt2ndJet = pt1stJet;
	      pt1stJet = JetPt->at(nJet);
	      qglik2ndJet = qglik1stJet;
	      qglik1stJet = JetQGLikelihood->at(nJet);
	    } else if (pt2ndJet < JetPt->at(nJet)) {
	      pt2ndJet = JetPt->at(nJet);
	      qglik2ndJet = JetQGLikelihood->at(nJet);
	    }
	    if (btag1stJet < JetBTagger->at(nJet)) {
	      btag2ndJet = btag1stJet;
	      btag1stJet = JetBTagger->at(nJet);
              isB1stJet = (JetIsBtaggedWithSF->at(nJet) > 0);
	    } else if (btag2ndJet < JetBTagger->at(nJet)) {
              btag2ndJet = JetBTagger->at(nJet);
              isB2ndJet = (JetIsBtaggedWithSF->at(nJet) > 0);
	    }
	    nInJets++;
	  } else {
            nExtraJets++;
	  }
        }
      }

      for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
	if (JetQGLikelihood->at(nJet) < -800.) {            // subjets
	  if (JetIsInZZCand->at(nJet) ) {
	    if (pt1stSubjet < JetPt->at(nJet)) {
	      pt2ndSubjet = pt1stSubjet;
	      pt1stSubjet = JetPt->at(nJet);
	    } else if (pt2ndSubjet < JetPt->at(nJet)) {
	      pt2ndSubjet = JetPt->at(nJet);
	    }
	    if (btag1stSubjet < JetBTagger->at(nJet)) {
	      btag2ndSubjet = btag1stSubjet;
	      btag1stSubjet = JetBTagger->at(nJet);
	    } else if (btag2ndSubjet < JetBTagger->at(nJet)) {
              btag2ndSubjet = JetBTagger->at(nJet);
	    }
	  }
	}
      }

      // find leading leptons
      float pt1stLep = 0.0001;
      float pt2ndLep = 0.0001;
      float eta1stLep = 0.0001;
      float eta2ndLep = 0.0001;
      float phi1stLep = 0.0001;
      float phi2ndLep = 0.0001;
      for (unsigned int nLep=0; nLep<LepPt->size(); nLep++) {
	if (pt1stLep < LepPt->at(nLep)) {
	  pt2ndLep = pt1stLep;
          eta2ndLep = eta1stLep;
          phi2ndLep = phi1stLep;
	  pt1stLep = LepPt->at(nLep);
          eta1stLep = LepEta->at(nLep);
          phi1stLep = LepPhi->at(nLep);
	} else if (pt2ndLep < LepPt->at(nLep)) {
	  pt2ndLep = LepPt->at(nLep);
          eta2ndLep = LepEta->at(nLep);
          phi2ndLep = LepPhi->at(nLep);
	}
      }

      //----- find lepton flavor

      int fsstart,fsend;
      if (abs(Z2Flav->at(0))==121) {
	fsstart=0;
	fsend=2;
      } else {
	fsstart=1;
	fsend=3;
      }

      int preferType = 0;    // choose merged or resolved
                             // and dump for synchronization
      if (ZZMass->size() == 1 && abs(ZZCandType->at(0)) == 1) {
	if (ZZMass->at(0) > 400) {
	  nPassMerged++;
	}
	preferType = 1;
      }
      else if (ZZMass->size() == 1 && abs(ZZCandType->at(0)) == 2) {
	if (ZZMass->at(0) > 400) {
	  nPassResol++;
	}
        preferType = 2;
      }
      else if (ZZMass->size() == 2 && abs(ZZCandType->at(0)) == 1) {
	if ( ZZMass->at(0) > 400 ) {
	  float mela[2] = {0.,0.};
	  float vbfmela[2] = {0.,0.};
	  for(unsigned int theCand=0; theCand<2; theCand++){
	    mela[theCand] = 1./(1.+getDZjjspin0Constant(ZZMass->at(theCand))*(pqqZJJ_VAMCFM->at(theCand)/p0plus_VAJHU->at(theCand)));
	    vbfmela[theCand] = ((pqqZJJ_VAMCFM->at(theCand) > 0. && pvbs_VAMCFM_highestPTJets->at(theCand) > 0. && nExtraJets > 1) ? 1./(1.+getDZjjspin0Constant(ZZMass->at(theCand))*(pqqZJJ_VAMCFM->at(theCand)/pvbs_VAMCFM_highestPTJets->at(theCand))) : -1.);
	  }
	  nPassMerged++; nPassResol++;

	}
	  // merged -> resolved (but ask for J quality)
        if ((Z2Pt->at(0) > 200. && Z1Pt->at(0) > 300. && Z1tau21->at(0) < 0.6)  || !(passAdditionalCuts(Z2Mass->at(1),Z1Pt->at(1),Z2Pt->at(1),0.,ZZMassRefit->at(1),false,false,draw)) ) preferType = 1;
        else preferType = 2;
      }

      // end dump for synchronization and choice

      // apply trigger weights
      if (weightMCtrig) {
	if (process > 0) {  //mc
	  if (abs(Z2Flav->at(0))==121) {
	    int bin1 = eff_mc_se->FindFixBin(eta1stLep,TMath::Min(pt1stLep,(float)199.0));
	    int bin2 = eff_mc_se->FindFixBin(eta2ndLep,TMath::Min(pt2ndLep,(float)199.0));
	    // mc eff
	    double mceff_se = eff_mc_se->GetEfficiency(bin1)+eff_mc_se->GetEfficiency(bin2)-eff_mc_se->GetEfficiency(bin1)*eff_mc_se->GetEfficiency(bin2);
	    double mceff_de = eff_mc_de1->GetEfficiency(bin1)*eff_mc_de2->GetEfficiency(bin2);
	    eventWeight *= mceff_se+mceff_de-mceff_se*mceff_de;
	  } else {
	    int bin1 = eff_mc_sm->FindFixBin(eta1stLep,TMath::Min(pt1stLep,(float)199.0));
	    int bin2 = eff_mc_sm->FindFixBin(eta2ndLep,TMath::Min(pt2ndLep,(float)199.0));
	    // mc eff
	    double mceff_sm = eff_mc_sm->GetEfficiency(bin1)+eff_mc_sm->GetEfficiency(bin2)-eff_mc_sm->GetEfficiency(bin1)*eff_mc_sm->GetEfficiency(bin2);
	    double mceff_dm = eff_mc_dm1->GetEfficiency(bin1)*eff_mc_dm2->GetEfficiency(bin2);
	    eventWeight *= mceff_sm+mceff_dm-mceff_sm*mceff_dm;
	  }
	}

      }
      // ----- fill histos
// TString outfileName_test( "RooFitInput.root" );

 //Output file creation
 // TFile* outputFile = TFile::Open( outfileName_test, "RECREATE" );
      for(int rs=fsstart; rs<fsend; rs++){

	for(unsigned int theCand=0; theCand<Z1Mass->size(); theCand++){

	  if (!CR) {   // regular events
            if (abs(ZZCandType->at(theCand)) != preferType) continue;

	    // redefine type
	    int typ = -1;
	    if (abs(ZZCandType->at(theCand)) == 1) {
	      if (Z1Mass->at(theCand) > 65. && Z1Mass->at(theCand) < 105.) typ = 2;
	      else typ = 1;
	    } else {
	      if (Z1Mass->at(theCand) > 65. && Z1Mass->at(theCand) < 105.) typ = 3;
	      else typ = 0;
	    }

	    // blind higgs region
	    if (!includeHiggsMass && Z1Mass->at(theCand) > 105. && Z1Mass->at(theCand) < 135.) continue;

	    int whichTmvaTree = -1;
	    if (typ==2 && process == 2) whichTmvaTree = 0;
	    if (typ==2 && process == 3) whichTmvaTree = 1;
	    if (typ==3 && process == 2) whichTmvaTree = 2;
 	    if (typ==3 && process == 3) whichTmvaTree = 3;

            float mela = ((pzzjj_VAMCFM_highestPTJets->at(theCand) > 0. && pvbs_VAMCFM_highestPTJets->at(theCand) > 0. && nExtraJets > 1) ? 1./(1.+1000.*getDZjjspin0Constant(ZZMass->at(theCand))*(pzzjj_VAMCFM_highestPTJets->at(theCand)/pvbs_VAMCFM_highestPTJets->at(theCand))) : -1.);
            float mela2 =  ((pqqZJJ_VAMCFM->at(theCand) > 0. && pvbf_VAJHU_highestPTJets->at(theCand) > 0. && nExtraJets > 1) ? 1./(1.+getDZjjspin0Constant(ZZMass->at(theCand))*(pqqZJJ_VAMCFM->at(theCand)/pvbf_VAJHU_highestPTJets->at(theCand))) : -1.);
	    float mela3 = ((pqqZJJ_VAMCFM->at(theCand) > 0. && pvbs_VAMCFM_highestPTJets->at(theCand) > 0. && nExtraJets > 1) ? 1./(1.+getDZjjspin0Constant(ZZMass->at(theCand))*(pqqZJJ_VAMCFM->at(theCand)/pvbs_VAMCFM_highestPTJets->at(theCand))) : -1.);

	    if ((typ==0 || typ==3) && nExtraJets > 1 && mela > 0. /*1.043-460./(ZZMass->at(theCand)+634.)*/ && DiJetMass > 100.) typ=typ+8;
	    if ((typ==1 || typ==2) && nExtraJets > 1 && mela > 0. /*1.043-460./(ZZMass->at(theCand)+634.)*/ && DiJetMass > 100.) typ=typ+8;

	    // if (typ > 7) cout << " VBS: " << pqqZJJ_VAMCFM->at(theCand) << " " << pvbs_VAMCFM_highestPTJets->at(theCand) << " " << vbfmela << endl;

	    if ((typ==0 || typ==3) && btag1stJet>0.5426 && btag2ndJet>0.5426)
               //isB1stJet && isB2ndJet)
               typ=typ+4;
	    if ((typ==1 || typ==2) && btag1stSubjet > 0.5426 && btag2ndSubjet > 0.5426) typ=typ+4;

	    // correct tau21
	    float t12weight = 1.;

            if (process == 3 && (typ == 1 || typ == 2) && weighttau21) {
              for (int itau = 0; itau < 39; itau++) {
                if (Z1tau21->at(theCand) > tau21bin[itau] && Z1tau21->at(theCand) < tau21bin[itau+1]) t12weight = 1.+tau21corr[itau];
              }
	    }


	    tmvaZZPt = (float)ZZPt->at(theCand);
	    tmvaZ2Mass = (float)Z2Mass->at(theCand);
	    tmvaZ1Pt = (float)Z1Pt->at(theCand);
	    tmvaZ2Pt = (float)Z2Pt->at(theCand);
	    tmvaLepPt1 = (float)pt1stLep;
	    tmvaLepPt2 = (float)pt2ndLep;
	    tmvaJetQGLikelihood1 = (float)qglik1stJet;
	    tmvaJetQGLikelihood2 = (float)qglik2ndJet;
	    tmvaabshelcosthetaZ1 = (float)abs(helcosthetaZ1->at(theCand));
	    tmvahelcosthetaZ2 = (float)helcosthetaZ2->at(theCand);
	    tmvacosthetastar = (float)costhetastar->at(theCand);
	    tmvahelphi = (float)helphi->at(theCand);
	    tmvaphistarZ1 = (float)phistarZ1->at(theCand);
	    tmvaZ1tau21 = (float)Z1tau21->at(theCand);
	    if (typ==0 || typ==3 || typ==4 || typ==7 || typ==8 || typ==11) {   // only resolved
	      tmvaJetPt1 = (float)pt1stJet;
	      tmvaJetPt2 = (float)pt2ndJet;
	    } else {
	      tmvaJetPt1 = (float)pt1stSubjet;
	      tmvaJetPt2 = (float)pt2ndSubjet;
	    }

            if ((typ==1 || typ==2 || typ==5 || typ==6 || typ==9 || typ==10) && !(passAdditionalCuts(tmvaZ2Mass,tmvaZ1Pt,tmvaZ2Pt,tmvaZ1tau21,ZZMass->at(theCand),true,syncevents,draw))) continue;
            if (!(typ==1 || typ==2 || typ==5 || typ==6 || typ==9 || typ==10) && !(passAdditionalCuts(tmvaZ2Mass,tmvaZ1Pt,tmvaZ2Pt,tmvaZ1tau21,ZZMassRefit->at(theCand),false,syncevents,draw))) continue;

	    h1[0][process][rs][typ]->Fill(ZZMass->at(theCand),eventWeight*t12weight);
	    h1[1][process][rs][typ]->Fill(ZZPt->at(theCand),eventWeight*t12weight);

	    h1[2][process][rs][typ]->Fill(ZZMassRefit->at(theCand),eventWeight*t12weight);
	    if (typ==0 || typ==3 || typ==4 || typ==7 || typ==8 || typ==11) {
	      if (rs == onlyOneLep) {
		hmass[process][typ]->Fill(ZZMassRefit->at(theCand),eventWeight*t12weight);
                hmass_up[process][typ]->Fill(ZZMassRefit->at(theCand)*ZZMass_up->at(theCand)/ZZMass->at(theCand),eventWeight*t12weight);
		hmass_down[process][typ]->Fill(ZZMassRefit->at(theCand)*ZZMass_dn->at(theCand)/ZZMass->at(theCand),eventWeight*t12weight);
		if ( typ == 0 || typ==8 ) {
		  hmass[process][12]->Fill(ZZMassRefit->at(theCand),eventWeight*t12weight);
                  hmass_up[process][12]->Fill(ZZMassRefit->at(theCand)*ZZMass_up->at(theCand)/ZZMass->at(theCand),eventWeight*t12weight);
		  hmass_down[process][12]->Fill(ZZMassRefit->at(theCand)*ZZMass_dn->at(theCand)/ZZMass->at(theCand),eventWeight*t12weight);
		}
		if ( typ == 3 || typ==11 ) {
		  hmass[process][15]->Fill(ZZMassRefit->at(theCand),eventWeight*t12weight);
                  hmass_up[process][15]->Fill(ZZMassRefit->at(theCand)*ZZMass_up->at(theCand)/ZZMass->at(theCand),eventWeight*t12weight);
		  hmass_down[process][15]->Fill(ZZMassRefit->at(theCand)*ZZMass_dn->at(theCand)/ZZMass->at(theCand),eventWeight*t12weight);
		}
	      }
	      if (mela > 0.5) h1[25][process][rs][typ]->Fill(ZZMassRefit->at(theCand),eventWeight*t12weight);
	    } else {
              if (rs == onlyOneLep) {
		hmass[process][typ]->Fill(ZZMass->at(theCand),eventWeight*t12weight);
                hmass_up[process][typ]->Fill(ZZMass_up->at(theCand),eventWeight*t12weight);
		hmass_down[process][typ]->Fill(ZZMass_dn->at(theCand),eventWeight*t12weight);
		if ( typ == 1 || typ== 9 ) {
		  hmass[process][13]->Fill(ZZMass->at(theCand),eventWeight*t12weight);
		  hmass_up[process][13]->Fill(ZZMass_up->at(theCand),eventWeight*t12weight);
		  hmass_down[process][13]->Fill(ZZMass_dn->at(theCand),eventWeight*t12weight);
		}
		if ( typ == 2 || typ== 10 ) {
		  hmass[process][14]->Fill(ZZMass->at(theCand),eventWeight*t12weight);
		  hmass_up[process][14]->Fill(ZZMass_up->at(theCand),eventWeight*t12weight);
		  hmass_down[process][14]->Fill(ZZMass_dn->at(theCand),eventWeight*t12weight);
		}
	      }
	      if (mela > 0.5) h1[25][process][rs][typ]->Fill(ZZMass->at(theCand),eventWeight*t12weight);
	    }

	    h1[3][process][rs][typ]->Fill(Z1Mass->at(theCand),eventWeight*t12weight);

	    h1[4][process][rs][typ]->Fill(Z2Mass->at(theCand),eventWeight*t12weight);
	    h1[5][process][rs][typ]->Fill(Z1Pt->at(theCand),eventWeight*t12weight);
	    h1[6][process][rs][typ]->Fill(Z2Pt->at(theCand),eventWeight*t12weight);
	    h1[7][process][rs][typ]->Fill(Z2Flav->at(theCand),eventWeight*t12weight);

	    if (typ==0 || typ==3 || typ==4 || typ==7 || typ==8 || typ==11) {   // only resolved
	      h1[8][process][rs][typ]->Fill(pt1stJet,eventWeight*t12weight);
	      h1[9][process][rs][typ]->Fill(pt2ndJet,eventWeight*t12weight);
	    } else {
	      h1[8][process][rs][typ]->Fill(pt1stSubjet,eventWeight*t12weight);
	      h1[9][process][rs][typ]->Fill(pt2ndSubjet,eventWeight*t12weight);
	    }

	    if (typ==0 || typ==3 || typ==4 || typ==7  || typ==8 || typ==11) {   // only resolved
	      for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
		if (JetIsInZZCand->at(nJet) && JetQGLikelihood->at(nJet) > -800. /*remove subjets of fat jet also included in this collection! */) {
		  h1[11][process][rs][typ]->Fill(JetBTagger->at(nJet),eventWeight*t12weight);
		  h1[10][process][rs][typ]->Fill(JetQGLikelihood->at(nJet),eventWeight*t12weight);
		}
	      }
	    } else {
	      for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
		if (JetIsInZZCand->at(nJet) && JetQGLikelihood->at(nJet) < -800. /*subjets of fat jet also included in this collection! */) {
		  h1[11][process][rs][typ]->Fill(JetBTagger->at(nJet),eventWeight*t12weight);
		}
	      }
	    }

	    h1[12][process][rs][typ]->Fill(pt1stLep,eventWeight*t12weight);
	    h1[13][process][rs][typ]->Fill(pt2ndLep,eventWeight*t12weight);

	    h1[14][process][rs][typ]->Fill(abs(helcosthetaZ1->at(theCand)),eventWeight*t12weight);
	    h1[15][process][rs][typ]->Fill(helcosthetaZ2->at(theCand),eventWeight*t12weight);
	    h1[16][process][rs][typ]->Fill(costhetastar->at(theCand),eventWeight*t12weight);
	    h1[17][process][rs][typ]->Fill(phistarZ1->at(theCand),eventWeight*t12weight);
	    h1[18][process][rs][typ]->Fill(fabs(helphi->at(theCand)),eventWeight*t12weight);

	    h1[19][process][rs][typ]->Fill(nExtraJets,eventWeight*t12weight);
	    h1[20][process][rs][typ]->Fill(Met,eventWeight*t12weight);

	    if (typ==1 || typ==2 || typ==5 || typ==6  || typ==9 || typ==10)
	      h1[21][process][rs][typ]->Fill(Z1tau21->at(theCand),eventWeight*t12weight);   // only merged
	    else
	      h1[22][process][rs][typ]->Fill(tmvaJetQGLikelihood1*tmvaJetQGLikelihood2,eventWeight*t12weight);     // only resolved
            if (process > 0 || mela < 0.75) h1[23][process][rs][typ]->Fill(mela,eventWeight*t12weight);
            if (process > 0 || mela2 < 0.75) h1[24][process][rs][typ]->Fill(mela2,eventWeight*t12weight);
            if (process > 0 || mela3 < 0.75) h1[26][process][rs][typ]->Fill(mela3,eventWeight*t12weight);
	    if (rs == onlyOneLep && (typ == 3 || typ == 7 || typ == 11) ) {
	      if (process > 0 || mela < 0.75) hmass[process][16]->Fill(mela,eventWeight*t12weight);
	      if (process > 0 || mela2 < 0.75) hmass[process][17]->Fill(mela2,eventWeight*t12weight);
	      if (process > 0 || mela3 < 0.75)  hmass[process][18]->Fill(mela3,eventWeight*t12weight);
		}

//2l2q 2D HISTOGRAM FILL

   //alpha method implementation
  // for(int rs=0; rs<nFS; rs++){
   //   AsResolved = h1[2][3][rs][11]
   //   AbResolved = h1[2][3][rs][8]
   //   BbResolved = h1[2][4][rs][8]
   //   CbResolved = h1[2][5][rs][8]
   //   DbResolved = h1[2][0][rs][8]
   //   AsMerged = h1[0][3][rs][10]
   //   AbMerged = h1[0][3][rs][9]
   // BbMerged = h1[0][4][rs][9]
   // CbMerged = h1[0][5][rs][9]
   // DbMerged = h1[0][0][rs][9]
   //}
//h_alpha = (As/Ab)*(Db-Bb-Cb)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Hist fill 2
	    // if (rs == 1 && whichTmvaTree > -1) outputTree[whichTmvaTree]->Fill();
	  } else { // control region for QG (only fille some variables)

	    for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
	      if (JetIsInZZCand->at(nJet) && JetQGLikelihood->at(nJet) > -800. /*remove subjets of fat jet also included in this collection! */) {
		if (deltaPhi(phi1stLep, JetPhi->at(nJet)) < 2.2 || deltaPhi(phi2ndLep, JetPhi->at(nJet)) < 2.2) continue;
	      }
	    }

	    h1[3][process][rs][3]->Fill(Z1Mass->at(theCand),eventWeight);
	    h1[4][process][rs][3]->Fill(Z2Mass->at(theCand),eventWeight);
	    h1[5][process][rs][3]->Fill(Z1Pt->at(theCand),eventWeight);
	    h1[6][process][rs][3]->Fill(Z2Pt->at(theCand),eventWeight);
	    h1[7][process][rs][3]->Fill(Z2Flav->at(theCand),eventWeight);

	    h1[8][process][rs][3]->Fill(pt1stJet,eventWeight);

	    for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
	      if (JetIsInZZCand->at(nJet) && JetQGLikelihood->at(nJet) > -800. /*remove subjets of fat jet also included in this collection! */) {
		h1[11][process][rs][3]->Fill(JetBTagger->at(nJet),eventWeight);
		h1[10][process][rs][3]->Fill(JetQGLikelihood->at(nJet),eventWeight);
	      }
	    }

	    h1[12][process][rs][3]->Fill(pt1stLep,eventWeight);

	    h1[13][process][rs][3]->Fill(pt2ndLep,eventWeight);
	 }

	}
/////////////
   //cout << endl << "first" << endl;
  // for(int typ = 0; typ < 5 ; typ++){
 // temp[rs][typ]->Write();
//	}
///////////////
 } }
   // }

    if (mass[d] > 0 && enforceNarrowWidth) NEvtNarrow[d] = NEvtNarrow[d] / entries;
    if (process==0 && sync) {
      myfile.close();
      float eff = float(nPassMerged)/float(NGenEvt[d]);
      cout<<"Pass merged analysis = "<<nPassMerged<<"/"<< NGenEvt[d] << " (" << eff*100. << "%)"<<endl;
      eff = float(nPassResol)/float(NGenEvt[d]);
      cout<<"Pass resolved analysis = "<<nPassResol<<"/"<< NGenEvt[d] << " (" << eff*100. << "%) ..."<<endl;
      break;
    }
  }


  if (draw > 0) {
    TCanvas c1;
    c1.cd();

    TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
    pad1->Draw();
    TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.02,0.95,0.35);
    pad2->Draw();

    for(int rs=0; rs<nFS; rs++){   //ee, mumu, or all
      for(int v=0; v<nVariables; v++){
	for(int nt=8; nt<nType; nt++){    // CHANGED

	  if (enforceNarrowWidth) {
	    cout << "Only selecting " << NEvtNarrow[0]*100. << "% events close to nominal mass (mimic narrow width)" << endl;
	    cout << "Only selecting " << NEvtNarrow[1]*100. << "% events close to nominal mass (mimic narrow width)" << endl;
	    h1[v][1][rs][nt]->Scale(1./NEvtNarrow[0]);
	    h1[v][2][rs][nt]->Scale(1./NEvtNarrow[1]);
	  }

	  h1[v][4][rs][nt]->Add(h1[v][5][rs][nt]);
	  h1[v][3][rs][nt]->Add(h1[v][4][rs][nt]);
          // h1[v][1][rs][nt]->Add(h1[v][2][rs][nt]);

	  h1[v][3][rs][nt]->GetXaxis()->SetTitle(varXLabel[v].c_str());
	  h1[v][3][rs][nt]->GetYaxis()->SetTitle(varYLabel[v].c_str());
	  h1[v][0][rs][nt]->GetXaxis()->SetTitle(varXLabel[v].c_str());
	  h1[v][0][rs][nt]->GetYaxis()->SetTitle(varYLabel[v].c_str());
	  h1[v][3][rs][nt]->SetFillStyle(1);
	  h1[v][3][rs][nt]->SetMinimum(0.1);
	  h1[v][3][rs][nt]->SetLineColor(kGreen+2);
	  h1[v][3][rs][nt]->SetFillColor(kGreen+2);
	  h1[v][4][rs][nt]->SetFillStyle(1);
	  h1[v][4][rs][nt]->SetMinimum(0.1);
	  h1[v][4][rs][nt]->SetLineColor(kYellow+2);
	  h1[v][4][rs][nt]->SetFillColor(kYellow+2);
	  h1[v][5][rs][nt]->SetFillStyle(1);
	  h1[v][5][rs][nt]->SetMinimum(0.1);
	  h1[v][5][rs][nt]->SetLineColor(kMagenta-3);
	  h1[v][5][rs][nt]->SetFillColor(kMagenta-3);
	  h1[v][2][rs][nt]->SetLineColor(kBlue-1);
	  h1[v][2][rs][nt]->SetLineWidth(3);
	  h1[v][1][rs][nt]->SetLineColor(kRed-1);
	  h1[v][1][rs][nt]->SetLineWidth(3);
	  h1[v][0][rs][nt]->SetMarkerStyle(20);
	  h1[v][0][rs][nt]->SetMinimum(0.1);

	}
      }
    }

    for(int nt=8; nt<nType+7; nt++){        // CHANGED
      if (enforceNarrowWidth) {
	cout << "Only selecting " << NEvtNarrow[0]*100. << "% events close to nominal mass (mimic narrow width)" << endl;
	cout << "Only selecting " << NEvtNarrow[1]*100. << "% events close to nominal mass (mimic narrow width)" << endl;
	// cout << "Only selecting " << NEvtNarrow[25]*100. << "% events close to nominal mass (mimic narrow width)" << endl;
	hmass[1][nt]->Scale(1./NEvtNarrow[0]);
	hmass[2][nt]->Scale(1./NEvtNarrow[1]);
        // hmass[6][nt]->Scale(1./NEvtNarrow[25]);
      }

      hmass[4][nt]->Add(hmass[5][nt]);
      hmass[3][nt]->Add(hmass[4][nt]);

      if (nt<nType+4) {
	densityHist(hmass[0][nt]);
	densityHist(hmass[1][nt]);
	densityHist(hmass[2][nt]);
	densityHist(hmass[3][nt]);
	densityHist(hmass[4][nt]);
	densityHist(hmass[5][nt]);
        if (nt<nType) {
	  hbkg[nt]->Add(hmass[4][nt]);
	  hbkg_up[nt]->Add(hmass[4][nt]);
	  hbkg_down[nt]->Add(hmass[4][nt]);
	}
      } else {
	if (nt==nType+4) hmass[3][nt]->GetXaxis()->SetTitle("D^{Zjj}_{bkg}");
	else if (nt==nType+5) hmass[3][nt]->GetXaxis()->SetTitle("D^{Zjj} (spin 2)");
        else hmass[3][nt]->GetXaxis()->SetTitle("D^{VBF}_{2jet}");
	hmass[3][nt]->GetYaxis()->SetTitle("Events / bin");
	float normer = hmass[3][nt]->Integral();
        hmass[0][nt]->Scale(1./hmass[0][nt]->Integral());
	hmass[3][nt]->Scale(1./normer);
        hmass[4][nt]->Scale(1./normer);
	hmass[5][nt]->Scale(1./normer);
	hmass[1][nt]->Scale(1./hmass[1][nt]->Integral());
	hmass[2][nt]->Scale(1./hmass[2][nt]->Integral());
	// hmass[6][nt]->Scale(1./hmass[6][nt]->Integral());
      }

      hmass[3][nt]->SetFillStyle(3001);
      if (nt<4) hmass[3][nt]->SetMinimum(0.11);
      else if (nt<16) hmass[3][nt]->SetMinimum(0.0011);
      else if (nt==18) hmass[3][nt]->SetMinimum(0.0001);
      else hmass[3][nt]->SetMinimum(0.);
      hmass[3][nt]->SetLineColor(kGreen+2);
      hmass[3][nt]->SetFillColor(kGreen+2);
      if (nt<12) hmass[3][nt]->GetXaxis()->SetLabelColor(0);
      hmass[4][nt]->SetFillStyle(3001);
      hmass[4][nt]->SetLineColor(kYellow+2);
      hmass[4][nt]->SetFillColor(kYellow+2);
      hmass[5][nt]->SetFillStyle(3001);
      hmass[5][nt]->SetLineColor(kMagenta-3);
      hmass[5][nt]->SetFillColor(kMagenta-3);
      hmass[2][nt]->SetLineColor(kRed-1);
      hmass[2][nt]->SetLineWidth(4);
      hmass[2][nt]->SetLineStyle(kDashed);
      hmass[1][nt]->SetLineColor(kRed-1);
      hmass[1][nt]->SetLineWidth(3);
      // hmass[6][nt]->SetLineColor(kOrange+1);
      // hmass[6][nt]->SetLineWidth(3);
      hmass[0][nt]->SetMarkerStyle(20);

    }

    //2l2q alpha method
    std::cout << endl << " GET READY: ALPHA METHOD! " << endl;
    TH1F* alpha_ee; //(As/Ab)*(Db-Bb-Cb) Resolved
    TH1F* temp_1; //AsResolved, then
    TH1F* temp_2; //AbResolved
    TH1F* temp_3; //AsResolved/AbResolved
    TH1F* temp_4; //DbResolved
    TH1F* temp_5; //BbResolved
    TH1F* temp_6; //CbResolved
    TH1F* temp_7; //Db-Bb Resolved
    TH1F* temp_8; //Db-Bb-Cb Resolved
    //for(int rs=0; rs<nFS; rs++){
    temp_1 = h1[2][3][1][11];//AsResolved
    temp_7 = (TH1F*) temp_1->Clone();

    //temp_7->Draw("hist");
    //gPad->Print("ALPHA/output.png");

    temp_2 = h1[2][3][1][8];  //AbResolved
    temp_1->Divide(temp_2);   //AsResolved/AbResolved
    temp_4 = h1[2][0][1][8];  //DbResolved
    temp_5 = h1[2][4][1][8];  //BbResolved
    temp_6 = h1[2][5][1][8];  //CbResolved
    temp_4->Add(temp_5,-1);   //Db-Bb Resolved
    temp_4->Add(temp_6,-1);   //Db-Bb-Cb Resolved
    //temp_1->Multiply(temp_4); //(As/Ab)*(Db-Bb-Cb) Resolved
    temp_1 -> SetLineColor(kRed);   //new
    temp_7 -> SetLineColor(kGreen);  //old

    temp_1 -> Draw();

    //define interpolating function and fit data
    TF1 *g1    = new TF1("g1","pol1",400,temp_1->GetBinCenter(20));
    temp_1->Fit(g1,"R");

    //cycle on histogram, cut oscillating values at high masses,
    //interpolate remaining values and extrapolate at high masses
    for (Int_t ii = 21; ii < varNbin[0]; ii++) {
        temp_1 -> SetBinContent(ii,g1->Eval(temp_1->GetBinCenter(ii)));
    }

    g1->SetLineColor(kGreen);
    g1->Draw("same");
    gPad->Print("ALPHA/RATIO_PLOT.png");
    std::cout << endl << " THAT'S ALL FOLKS! " << endl;
      // }
    //Plot for alpha METHOD

    // //RATIO PLOT
    //     //canvas
    //     TCanvas *ratio_plot = new TCanvas("ratio_plot","ratio_plot",800,1000);
    //     //pad1
    //     float eps =0;// 0.006;
    //     TPad *first_pad = new TPad("first_pad","first_pad",0,0.3,1,1,0);
    //     first_pad ->SetBottomMargin(0);
    //     first_pad ->Draw();
    //     first_pad ->cd();
    //     //top plot
    //     temp_1 -> SetLineColor(kRed);   //new
    //     temp_7 -> SetLineColor(kGreen);  //old
    //     temp_1 -> Draw("hist");
    //     temp_7->Draw("hist same");
    //     //legend->Draw("same");
    //     //switch?
    //     ratio_plot->cd();
    //     //pad2
    //     TPad *second_pad  = new TPad("second_pad","second_pad",0,0.1,1,0.3+eps,0); //old position
    //     second_pad->SetTopMargin(0);
    //     //make bottom pad transparent
    //     //pad2->SetFrameFillColor(0);
    //     //pad2->SetFrameBorderMode(0);
    //     //pad2->SetFrameFillColor(0);
    //     //pad2->SetFrameBorderMode(0);
    //
    //     second_pad->Draw();
    //     second_pad->cd();
    //     //bottom plot
    //     TH1F *h0copy = (TH1F*) temp_1->Clone();
    //     //axis labels
    //     h0copy->GetXaxis()->SetLabelFont(59);//change this for font type
    //     h0copy->GetXaxis()->SetLabelSize(22);
    //     h0copy->GetYaxis()->SetLabelFont(59);//change this for font type
    //     h0copy->GetYaxis()->SetLabelSize(22);
    //     //axis titles
    //     h0copy->GetXaxis()->SetTitleFont(59); //change this for font type
    //     h0copy->GetXaxis()->SetTitleSize(22);
    //     h0copy->GetYaxis()->SetTitleFont(59); //change this for font type
    //     h0copy->GetYaxis()->SetTitleSize(22);
    //     h0copy->GetXaxis()->SetTitleOffset(4.5);
    //     h0copy->GetYaxis()->SetTitleOffset(1.7);
    //     h0copy->GetYaxis()->SetRangeUser(-1.,3.);
    //
    //     h0copy->Sumw2();
    //     h0copy->SetStats(0); //clear stat box
    //     h0copy->Divide(temp_7); //invert divide
    //     h0copy->SetMarkerStyle(20);
    //     //h0copy->SetTitle("; m_ZZ; new/old");
    //     //h0copy->GetXaxis()->SetTitleSize(50);
    //     //h0copy->GetXaxis()->SetLabelSize(35);
    //     //h0copy->GetYaxis()->SetTitleSize(15);
    //     //h0copy->GetYaxis()->SetLabelSize(15);
    //
    //     //gStyle->SetLabelSize(2,"x").
    //     h0copy->Draw("ep");
    //
    //     //Orizontal line
    //     TLine *line = new TLine(0,1,1,1);
    //     line->SetLineColor(kRed);
    //     line->SetLineStyle(2);
    //     line->Draw("same");
    //
    //     //close and print on file
    //     ratio_plot->cd();
    //    gPad->Print("ALPHA/RATIO_PLOT.png");


    //2l2q hist rebinning and saving

    //String name
    //TString outfileName_test( "RooFitInput.root" );

    //Output file creation
    //TFile* outputFile = TFile::Open( outfileName_test, "RECREATE" );

    //std::cout<< endl << "Output file correctly opened" << endl;
    // outputFile->cd();
    // std::cout<< endl << "Start writing on tree..." << endl;
    // //TTree writing
    // //tnew->Write();
    // //Cycle on hist
    //     for(int rs=0; rs<nFS; rs++){   // ee, mumu, or all
    //
    // 	 for(int nt=0; nt<12; nt++){      //Changed
    //
    //                 for(int binx=0;binx<temp[rs][nt]->GetXaxis()->GetNbins();binx++){
    //                         double inttmp1 = temp[rs][nt]->Integral(binx+1,binx+1);
    //                         for(int biny=0;biny<temp[rs][nt]->GetNbinsY();biny++){
    //                                 temp[rs][nt]->SetBinContent(binx+1, biny+1, temp[rs][nt]->GetBinContent(binx+1,biny+1)/inttmp1);
    //                         }
    //
    //                 }
    //                 temp[rs][nt]->Smooth();
    //                 int nbinsfinal=455;
    //                 double xbinfinal[456]={
    //                         110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,150,152,154,156,158,160,162,164,166,168,170,172,174,176,178,180,182,184,186,188,190,192,194,196,198,200,202,204,206,208,210,212,214,216,218,220,222,224,226,228,230,232,234,236,238,240,242,244,246,248,250,252,254,256,258,260,262,264,266,268,270,272,274,276,278,280,282,284,286,288,290,292,294,296,298,300,302,304,306,308,310,312,314,316,318,320,322,324,326,328,330,332,334,336,338,340,342,344,346,348,350,352,354,356,358,360,362,364,366,368,370,372,374,376,378,380,382,384,386,388,390,392,394,396,398,400,410,420,430,440,450,460,470,480,490,500,510,520,530,540,550,560,570,580,590,600,610,620,630,640,650,660,670,680,690,700,710,720,730,740,750,760,770,780,790,800,810,820,830,840,850,860,870,880,890,900,910,920,930,940,950,960,970,980,990,1000,1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,1110,1120,1130,1140,1150,1160,1170,1180,1190,1200,1210,1220,1230,1240,1250,1260,1270,1280,1290,1300,1310,1320,1330,1340,1350,1360,1370,1380,1390,1400,1410,1420,1430,1440,1450,1460,1470,1480,1490,1500,1510,1520,1530,1540,1550,1560,1570,1580,1590,1600,1610,1620,1630,1640,1650,1660,1670,1680,1690,1700,1710,1720,1730,1740,1750,1760,1770,1780,1790,1800,1810,1820,1830,1840,1850,1860,1870,1880,1890,1900,1910,1920,1930,1940,1950,1960,1970,1980,1990,2000,2010,2020,2030,2040,2050,2060,2070,2080,2090,2100,2110,2120,2130,2140,2150,2160,2170,2180,2190,2200,2210,2220,2230,2240,2250,2260,2270,2280,2290,2300,2310,2320,2330,2340,2350,2360,2370,2380,2390,2400,2410,2420,2430,2440,2450,2460,2470,2480,2490,2500,2510,2520,2530,2540,2550,2560,2570,2580,2590,2600,2610,2620,2630,2640,2650,2660,2670,2680,2690,2700,2710,2720,2730,2740,2750,2760,2770,2780,2790,2800,2810,2820,2830,2840,2850,2860,2870,2880,2890,2900,2910,2920,2930,2940,2950,2960,2970,2980,2990,3000,3010,3020,3030,3040,3050,3060,3070,3080,3090,3100,3110,3120,3130,3140,3150,3160,3170,3180,3190,3200,3210,3220,3230,3240,3250,3260,3270,3280,3290,3300,3310,3320,3330,3340,3350,3360,3370,3380,3390,3400,3410,3420,3430,3440,3450,3460,3470,3480,3490,3500
    //                         };
    //
    //                 for(int binx=0;binx<temp_rebin[rs][nt]->GetXaxis()->GetNbins();binx++){
    //                         int b = temp[rs][nt]->GetXaxis()->FindBin(temp_rebin[rs][nt]->GetXaxis()->GetBinCenter(binx+1));
    //                         for(int biny=0;biny<30;biny++){
    //                                 float bc_4e = temp[rs][nt]->GetBinContent(b,biny+1);
    //                                 temp_rebin[rs][nt]->SetBinContent(binx+1,biny+1,bc_4e);
    //                         }
    //                 }
    //
    //                 //draw and save hist
    //                 temp_rebin[rs][nt]->Draw("colz");
    //                 gPad->Print( Form("templates/template_%s_%s.png",sFS[rs].c_str(),typeS[nt].c_str()));
    //
    // 		//write on tree
    //
    // 		//temp_rebin[rs][nt]->Write();
    //     //temp_rebin[rs][nt]->SetDirectory(outputFile);
    // 		std::cout<< endl << rs << "  " << nt << endl;
    //                 if (nt == 3) nt=7;//???
    //     }
    //
    //          }
    // tnew->Write();
    // outputFile->Close();
    //
    // std::cout<< endl << "Output file correctly closed" << endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    }}
