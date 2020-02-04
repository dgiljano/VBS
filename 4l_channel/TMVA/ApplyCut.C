 #include <cstdlib>
 #include <vector>
 #include <iostream>
 #include <map>
 #include <string>
 
 #include "TFile.h"
 #include "TTree.h"
 #include "TString.h"
 #include "TSystem.h"
 #include "TROOT.h"
 #include "TStopwatch.h"
 
 #include "TMVA/Tools.h"
 #include "TMVA/Reader.h"
 #include "TMVA/MethodCuts.h"
 
 using namespace TMVA;
 
 void ApplyCut( TString myMethodList = "" )
 {
 
    //---------------------------------------------------------------
    // This loads the library
    TMVA::Tools::Instance();
 
    // Default MVA methods to be trained + tested
    std::map<std::string,int> Use;
 
    // Cut optimisation
    Use["Cuts"]            = 1;
    Use["CutsD"]           = 1;
    Use["CutsPCA"]         = 0;
    Use["CutsGA"]          = 0;
    Use["CutsSA"]          = 0;
    //
    // 1-dimensional likelihood ("naive Bayes estimator")
    Use["Likelihood"]      = 1;
    Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
    Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
    Use["LikelihoodKDE"]   = 0;
    Use["LikelihoodMIX"]   = 0;
    //
    // Mutidimensional likelihood and Nearest-Neighbour methods
    Use["PDERS"]           = 1;
    Use["PDERSD"]          = 0;
    Use["PDERSPCA"]        = 0;
    Use["PDEFoam"]         = 1;
    Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
    Use["KNN"]             = 1; // k-nearest neighbour method
    //
    // Linear Discriminant Analysis
    Use["LD"]              = 1; // Linear Discriminant identical to Fisher
    Use["Fisher"]          = 0;
    Use["FisherG"]         = 0;
    Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
    Use["HMatrix"]         = 0;
    //
    // Function Discriminant analysis
    Use["FDA_GA"]          = 1; // minimisation of user-defined function using Genetics Algorithm
    Use["FDA_SA"]          = 0;
    Use["FDA_MC"]          = 0;
    Use["FDA_MT"]          = 0;
    Use["FDA_GAMT"]        = 0;
    Use["FDA_MCMT"]        = 0;
    //
    // Neural Networks (all are feed-forward Multilayer Perceptrons)
    Use["MLP"]             = 0; // Recommended ANN
    Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
    Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
    Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
    Use["TMlpANN"]         = 0; // ROOT's own ANN
    Use["DNN_CPU"] = 0;         // CUDA-accelerated DNN training.
    Use["DNN_GPU"] = 0;         // Multi-core accelerated DNN.
    //
    // Support Vector Machine
    Use["SVM"]             = 1;
    //
    // Boosted Decision Trees
    Use["BDT"]             = 1; // uses Adaptive Boost
    Use["BDTG"]            = 0; // uses Gradient Boost
    Use["BDTB"]            = 0; // uses Bagging
    Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
    Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
    //
    // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
    Use["RuleFit"]         = 1;
    // ---------------------------------------------------------------
    Use["Plugin"]          = 0;
    Use["Category"]        = 0;
    Use["SVM_Gauss"]       = 0;
    Use["SVM_Poly"]        = 0;
    Use["SVM_Lin"]         = 0;
 
    std::cout << std::endl;
    std::cout << "==> Start ApplyCut" << std::endl;
 
    // Select methods (don't look at this code - not of interest)
    if (myMethodList != "") {
       for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
 
       std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
       for (UInt_t i=0; i<mlist.size(); i++) {
          std::string regMethod(mlist[i]);
 
          if (Use.find(regMethod) == Use.end()) {
             std::cout << "Method \"" << regMethod
                       << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
             for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
                std::cout << it->first << " ";
             }
             std::cout << std::endl;
             return;
          }
          Use[regMethod] = 1;
       }
    }
 
    // --------------------------------------------------------------------------------------------------
 
    // Create the Reader object
 
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

    Float_t mjj, deta_jj, m4l, eta_star_Z1, eta_star_Z2, R_pt_hard, R_pt_jet, weight;
    //float_t abs_etajet_min, abs_etajet_max, abs_etalep_min, abs_etalep_max, delta_phi_ZZ, rapidity_Z1, rapidity_Z2, rapidity_j1, rapidity_j2, pt_Z1, pt_Z2, pt_l3, jet1_qg_tagger, jet2_qg_tagger, eta_j1, eta_j2, pt_jet1, pt_jet2, eta_j_sum, mjj_over_detajj, abs_etajet_sum;

    reader->AddVariable( "mjj", &mjj);
    reader->AddVariable( "deta_jj", &deta_jj);
    reader->AddVariable( "m4l", &m4l);
    reader->AddVariable( "eta_star_Z1", &eta_star_Z1);
    reader->AddVariable( "eta_star_Z2", &eta_star_Z2);
    reader->AddVariable( "R_pt_hard", &R_pt_hard);
    reader->AddVariable( "R_pt_jet", &R_pt_jet);

    /*reader->AddVariable( "abs_etajet_min", &abs_etajet_min);
    reader->AddVariable( "abs_etajet_max", &abs_etajet_max);
    reader->AddVariable( "abs_etalep_min", &abs_etalep_min);
    reader->AddVariable( "abs_etalep_max", &abs_etalep_max);
    reader->AddVariable( "delta_phi_ZZ", &delta_phi_ZZ);
    reader->AddVariable( "rapidity_Z1", &rapidity_Z1);
    reader->AddVariable( "rapidity_Z2", &rapidity_Z2);
    reader->AddVariable( "rapidity_j1", &rapidity_j1);
    reader->AddVariable( "rapidity_j2", &rapidity_j2);
    reader->AddVariable( "pt_Z1", &pt_Z1);
    reader->AddVariable( "pt_Z2", &pt_Z2);
    reader->AddVariable( "pt_l3", &pt_l3);
    reader->AddVariable( "jet1_qg_tagger", &jet1_qg_tagger);
    reader->AddVariable( "jet2_qg_tagger", &jet2_qg_tagger);
    reader->AddVariable( "eta_j1", &eta_j1);
    reader->AddVariable( "eta_j2", &eta_j2);
    reader->AddVariable( "pt_jet1", &pt_jet1);
    reader->AddVariable( "pt_jet2", &pt_jet2);
    reader->AddVariable( "eta_j_sum", &eta_j_sum);
    reader->AddVariable( "mjj_over_detajj", &mjj_over_detajj);
    reader->AddVariable( "abs_etajet_sum", &abs_etajet_sum);*/
 
    // Book the MVA methods
 
    TString dir    = "dataset/weights/";
    TString prefix = "TMVAClassification";
 
    // Book method(s)
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
       if (it->second) {
          TString methodName = TString(it->first) + TString(" method");
          TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
          reader->BookMVA( methodName, weightfile );
       }
    }
 
    // Prepare input tree (this must be replaced by your data source)
    // in this example, there is a toy tree with signal and one with background events
    // we'll later on use only the "signal" events for the test in this example.
    //
    TFile *input = TFile::Open( "MVA.root" );

    std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;
 
    // Event loop
 
    std::cout << "--- Selecting samples" << std::endl;
    TTree* t_signal = (TTree*)input->Get("signal");
    TTree* t_qq = (TTree*)input->Get("qq");
    TTree* t_gg = (TTree*)input->Get("gg");
    TTree* t_ttzwwz = (TTree*)input->Get("ttzwwz");
    TTree* t_zx = (TTree*)input->Get("zx");
    

    t_signal->SetBranchAddress( "mjj", &mjj);
    t_signal->SetBranchAddress( "deta_jj", &deta_jj);
    t_signal->SetBranchAddress( "m4l", &m4l);
    t_signal->SetBranchAddress( "eta_star_Z1", &eta_star_Z1);
    t_signal->SetBranchAddress( "eta_star_Z2", &eta_star_Z2);
    t_signal->SetBranchAddress( "R_pt_hard", &R_pt_hard);
    t_signal->SetBranchAddress( "R_pt_jet", &R_pt_jet);
    t_signal->SetBranchAddress( "weight", &weight);

    /*t_signal->SetBranchAddress( "abs_etajet_min", &abs_etajet_min);
    t_signal->SetBranchAddress( "abs_etajet_max", &abs_etajet_max);
    t_signal->SetBranchAddress( "abs_etalep_min", &abs_etalep_min);
    t_signal->SetBranchAddress( "abs_etalep_max", &abs_etalep_max);
    t_signal->SetBranchAddress( "delta_phi_ZZ", &delta_phi_ZZ);
    t_signal->SetBranchAddress( "rapidity_Z1", &rapidity_Z1);
    t_signal->SetBranchAddress( "rapidity_Z2", &rapidity_Z2);
    t_signal->SetBranchAddress( "rapidity_j1", &rapidity_j1);
    t_signal->SetBranchAddress( "rapidity_j2", &rapidity_j2);
    t_signal->SetBranchAddress( "pt_Z1", &pt_Z1);
    t_signal->SetBranchAddress( "pt_Z2", &pt_Z2);
    t_signal->SetBranchAddress( "pt_l3", &pt_l3);
    t_signal->SetBranchAddress( "jet1_qg_tagger", &jet1_qg_tagger);
    t_signal->SetBranchAddress( "jet2_qg_tagger", &jet2_qg_tagger);
    t_signal->SetBranchAddress( "eta_j1", &eta_j1);
    t_signal->SetBranchAddress( "eta_j2", &eta_j2);
    t_signal->SetBranchAddress( "pt_jet1", &pt_jet1);
    t_signal->SetBranchAddress( "pt_jet2", &pt_jet2);
    t_signal->SetBranchAddress( "eta_j_sum", &eta_j_sum);
    t_signal->SetBranchAddress( "mjj_over_detajj", &mjj_over_detajj);
    t_signal->SetBranchAddress( "abs_etajet_sum", &abs_etajet_sum);*/

    //--------------------------------------------------------

    t_qq->SetBranchAddress( "mjj", &mjj);
    t_qq->SetBranchAddress( "deta_jj", &deta_jj);
    t_qq->SetBranchAddress( "m4l", &m4l);
    t_qq->SetBranchAddress( "eta_star_Z1", &eta_star_Z1);
    t_qq->SetBranchAddress( "eta_star_Z2", &eta_star_Z2);
    t_qq->SetBranchAddress( "R_pt_hard", &R_pt_hard);
    t_qq->SetBranchAddress( "R_pt_jet", &R_pt_jet);
    t_qq->SetBranchAddress( "weight", &weight);

    /*t_qq->SetBranchAddress( "abs_etajet_min", &abs_etajet_min);
    t_qq->SetBranchAddress( "abs_etajet_max", &abs_etajet_max);
    t_qq->SetBranchAddress( "abs_etalep_min", &abs_etalep_min);
    t_qq->SetBranchAddress( "abs_etalep_max", &abs_etalep_max);
    t_qq->SetBranchAddress( "delta_phi_ZZ", &delta_phi_ZZ);
    t_qq->SetBranchAddress( "rapidity_Z1", &rapidity_Z1);
    t_qq->SetBranchAddress( "rapidity_Z2", &rapidity_Z2);
    t_qq->SetBranchAddress( "rapidity_j1", &rapidity_j1);
    t_qq->SetBranchAddress( "rapidity_j2", &rapidity_j2);
    t_qq->SetBranchAddress( "pt_Z1", &pt_Z1);
    t_qq->SetBranchAddress( "pt_Z2", &pt_Z2);
    t_qq->SetBranchAddress( "pt_l3", &pt_l3);
    t_qq->SetBranchAddress( "jet1_qg_tagger", &jet1_qg_tagger);
    t_qq->SetBranchAddress( "jet2_qg_tagger", &jet2_qg_tagger);
    t_qq->SetBranchAddress( "eta_j1", &eta_j1);
    t_qq->SetBranchAddress( "eta_j2", &eta_j2);
    t_qq->SetBranchAddress( "pt_jet1", &pt_jet1);
    t_qq->SetBranchAddress( "pt_jet2", &pt_jet2);
    t_qq->SetBranchAddress( "eta_j_sum", &eta_j_sum);
    t_qq->SetBranchAddress( "mjj_over_detajj", &mjj_over_detajj);
    t_qq->SetBranchAddress( "abs_etajet_sum", &abs_etajet_sum);*/

    //--------------------------------------------------------

    t_gg->SetBranchAddress( "mjj", &mjj );
    t_gg->SetBranchAddress( "deta_jj", &deta_jj );
    t_gg->SetBranchAddress( "m4l", &m4l );
    t_gg->SetBranchAddress( "eta_star_Z1", &eta_star_Z1 );
    t_gg->SetBranchAddress( "eta_star_Z2", &eta_star_Z2 );
    t_gg->SetBranchAddress( "R_pt_hard", &R_pt_hard );
    t_gg->SetBranchAddress( "R_pt_jet", &R_pt_jet );
    t_gg->SetBranchAddress( "weight", &weight);

    /*t_gg->SetBranchAddress( "abs_etajet_min", &abs_etajet_min);
    t_gg->SetBranchAddress( "abs_etajet_max", &abs_etajet_max);
    t_gg->SetBranchAddress( "abs_etalep_min", &abs_etalep_min);
    t_gg->SetBranchAddress( "abs_etalep_max", &abs_etalep_max);
    t_gg->SetBranchAddress( "delta_phi_ZZ", &delta_phi_ZZ);
    t_gg->SetBranchAddress( "rapidity_Z1", &rapidity_Z1);
    t_gg->SetBranchAddress( "rapidity_Z2", &rapidity_Z2);
    t_gg->SetBranchAddress( "rapidity_j1", &rapidity_j1);
    t_gg->SetBranchAddress( "rapidity_j2", &rapidity_j2);
    t_gg->SetBranchAddress( "pt_Z1", &pt_Z1);
    t_gg->SetBranchAddress( "pt_Z2", &pt_Z2);
    t_gg->SetBranchAddress( "pt_l3", &pt_l3);
    t_gg->SetBranchAddress( "jet1_qg_tagger", &jet1_qg_tagger);
    t_gg->SetBranchAddress( "jet2_qg_tagger", &jet2_qg_tagger);
    t_gg->SetBranchAddress( "eta_j1", &eta_j1);
    t_gg->SetBranchAddress( "eta_j2", &eta_j2);
    t_gg->SetBranchAddress( "pt_jet1", &pt_jet1);
    t_gg->SetBranchAddress( "pt_jet2", &pt_jet2);
    t_gg->SetBranchAddress( "eta_j_sum", &eta_j_sum);
    t_gg->SetBranchAddress( "mjj_over_detajj", &mjj_over_detajj);
    t_gg->SetBranchAddress( "abs_etajet_sum", &abs_etajet_sum);*/

    //-------------------------------------------------------

    t_ttzwwz->SetBranchAddress( "mjj", &mjj );
    t_ttzwwz->SetBranchAddress( "deta_jj", &deta_jj );
    t_ttzwwz->SetBranchAddress( "m4l", &m4l );
    t_ttzwwz->SetBranchAddress( "eta_star_Z1", &eta_star_Z1 );
    t_ttzwwz->SetBranchAddress( "eta_star_Z2", &eta_star_Z2 );
    t_ttzwwz->SetBranchAddress( "R_pt_hard", &R_pt_hard );
    t_ttzwwz->SetBranchAddress( "R_pt_jet", &R_pt_jet );
    t_ttzwwz->SetBranchAddress( "weight", &weight);

    /*t_ttzwwz->SetBranchAddress( "abs_etajet_min", &abs_etajet_min);
    t_ttzwwz->SetBranchAddress( "abs_etajet_max", &abs_etajet_max);
    t_ttzwwz->SetBranchAddress( "abs_etalep_min", &abs_etalep_min);
    t_ttzwwz->SetBranchAddress( "abs_etalep_max", &abs_etalep_max);
    t_ttzwwz->SetBranchAddress( "delta_phi_ZZ", &delta_phi_ZZ);
    t_ttzwwz->SetBranchAddress( "rapidity_Z1", &rapidity_Z1);
    t_ttzwwz->SetBranchAddress( "rapidity_Z2", &rapidity_Z2);
    t_ttzwwz->SetBranchAddress( "rapidity_j1", &rapidity_j1);
    t_ttzwwz->SetBranchAddress( "rapidity_j2", &rapidity_j2);
    t_ttzwwz->SetBranchAddress( "pt_Z1", &pt_Z1);
    t_ttzwwz->SetBranchAddress( "pt_Z2", &pt_Z2);
    t_ttzwwz->SetBranchAddress( "pt_l3", &pt_l3);
    t_ttzwwz->SetBranchAddress( "jet1_qg_tagger", &jet1_qg_tagger);
    t_ttzwwz->SetBranchAddress( "jet2_qg_tagger", &jet2_qg_tagger);
    t_ttzwwz->SetBranchAddress( "eta_j1", &eta_j1);
    t_ttzwwz->SetBranchAddress( "eta_j2", &eta_j2);
    t_ttzwwz->SetBranchAddress( "pt_jet1", &pt_jet1);
    t_ttzwwz->SetBranchAddress( "pt_jet2", &pt_jet2);
    t_ttzwwz->SetBranchAddress( "eta_j_sum", &eta_j_sum);
    t_ttzwwz->SetBranchAddress( "mjj_over_detajj", &mjj_over_detajj);
    t_ttzwwz->SetBranchAddress( "abs_etajet_sum", &abs_etajet_sum);*/

    //-------------------------------------------------------

    t_zx->SetBranchAddress( "mjj", &mjj );
    t_zx->SetBranchAddress( "deta_jj", &deta_jj );
    t_zx->SetBranchAddress( "m4l", &m4l );
    t_zx->SetBranchAddress( "eta_star_Z1", &eta_star_Z1 );
    t_zx->SetBranchAddress( "eta_star_Z2", &eta_star_Z2 );
    t_zx->SetBranchAddress( "R_pt_hard", &R_pt_hard );
    t_zx->SetBranchAddress( "R_pt_jet", &R_pt_jet );
    t_zx->SetBranchAddress( "weight", &weight);

    /*t_zx->SetBranchAddress( "abs_etajet_min", &abs_etajet_min);
    t_zx->SetBranchAddress( "abs_etajet_max", &abs_etajet_max);
    t_zx->SetBranchAddress( "abs_etalep_min", &abs_etalep_min);
    t_zx->SetBranchAddress( "abs_etalep_max", &abs_etalep_max);
    t_zx->SetBranchAddress( "delta_phi_ZZ", &delta_phi_ZZ);
    t_zx->SetBranchAddress( "rapidity_Z1", &rapidity_Z1);
    t_zx->SetBranchAddress( "rapidity_Z2", &rapidity_Z2);
    t_zx->SetBranchAddress( "rapidity_j1", &rapidity_j1);
    t_zx->SetBranchAddress( "rapidity_j2", &rapidity_j2);
    t_zx->SetBranchAddress( "pt_Z1", &pt_Z1);
    t_zx->SetBranchAddress( "pt_Z2", &pt_Z2);
    t_zx->SetBranchAddress( "pt_l3", &pt_l3);
    t_zx->SetBranchAddress( "jet1_qg_tagger", &jet1_qg_tagger);
    t_zx->SetBranchAddress( "jet2_qg_tagger", &jet2_qg_tagger);
    t_zx->SetBranchAddress( "eta_j1", &eta_j1);
    t_zx->SetBranchAddress( "eta_j2", &eta_j2);
    t_zx->SetBranchAddress( "pt_jet1", &pt_jet1);
    t_zx->SetBranchAddress( "pt_jet2", &pt_jet2);
    t_zx->SetBranchAddress( "eta_j_sum", &eta_j_sum);
    t_zx->SetBranchAddress( "mjj_over_detajj", &mjj_over_detajj);
    t_zx->SetBranchAddress( "abs_etajet_sum", &abs_etajet_sum);*/

    //-------------------------------------------------------

    TH1F* histo_signal = new TH1F("bkg_vbs", "histo_signal", 25, -0.3, 0.5);
    TH1F* histo_qq = new TH1F("bkg_qqzz", "histo_qq", 25, -0.3, 0.5);
    TH1F* histo_gg = new TH1F("bkg_ggzz", "histo_gg", 25, -0.3, 0.5);
    TH1F* histo_ttzwwz = new TH1F("bkg_ttzwzz", "histo_ttzwwz", 25, -0.3, 0.5);
    TH1F* histo_zx = new TH1F("bkg_zjet", "histo_zx", 25, -0.3, 0.5);
    TH1F* histo_data = new TH1F("data_obs", "data_obs", 25, -0.3, 0.5);
 
 
    std::cout << "--- Processing: " << t_signal->GetEntries() << " events" << std::endl;
    TStopwatch sw;
    sw.Start();

    TH1F* histo_qq_weight = new TH1F("histo_qq_weight", "histo_qq_weight", 25, -0.5, -0.5);

    for (Long64_t ievt=0; ievt<t_signal->GetEntries();ievt++) 
    {
        t_signal->GetEntry(ievt);
        float tBDT_signal = reader->EvaluateMVA("BDT method");
        histo_signal->Fill(tBDT_signal, weight);
    }

    for (Long64_t ievt=0; ievt<t_qq->GetEntries();ievt++) 
    {
        t_qq->GetEntry(ievt);
        float tBDT_qq = reader->EvaluateMVA("BDT method");
        histo_qq->Fill(tBDT_qq, weight);

        histo_qq_weight->Fill(weight);
    }

    for (Long64_t ievt=0; ievt<t_gg->GetEntries();ievt++) 
    {
        t_gg->GetEntry(ievt);
        float tBDT_gg = reader->EvaluateMVA("BDT method");
        histo_gg->Fill(tBDT_gg, weight);
    }

    for (Long64_t ievt=0; ievt<t_ttzwwz->GetEntries();ievt++) 
    {
        t_ttzwwz->GetEntry(ievt);
        float tBDT_ttzwwz = reader->EvaluateMVA("BDT method");
        histo_ttzwwz->Fill(tBDT_ttzwwz, weight);
    }

    for (Long64_t ievt=0; ievt<t_zx->GetEntries();ievt++) 
    {
        t_zx->GetEntry(ievt);
        float tBDT_zx = reader->EvaluateMVA("BDT method");
        histo_zx->Fill(tBDT_zx, weight);
    }
 
    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();

    /*histo_signal->Draw();
    histo_qq->Draw("same");
    histo_gg->Draw("same");
    histo_ttzwwz->Draw("same");
    histo_zx->Draw("same");*/

    //TFile *bdt_stack_histos = new TFile("bdt_stack_histos.root", "RECREATE");
    TFile *bdt_stack_histos = new TFile("BDT_stack_2018.root", "RECREATE");
    histo_signal->Write();
    histo_qq->Write();
    histo_gg->Write();
    histo_ttzwwz->Write();
    histo_zx->Write();
    histo_data->Write();

    // ------------------------------------ stacking histograms --------------------------------------------------------
    //THStack *hs = new THStack();
    char filetitle[300];
	sprintf(filetitle,"CMS Preliminary                                                                35.9 fb^{-1}");  
    THStack *hs = new THStack("",filetitle);

    histo_signal->SetFillColor(kMagenta);
    histo_qq->SetFillColor(kCyan);
    histo_gg->SetFillColor(kBlue);
    histo_ttzwwz->SetFillColor(kYellow);
    histo_zx->SetFillColor(kGreen);

    histo_ttzwwz->Add(histo_ttzwwz, histo_zx, 1, 1);
    histo_qq->Add(histo_ttzwwz, histo_qq, 1, 1);
    histo_gg->Add(histo_qq, histo_gg, 1, 1);
    histo_signal->Add(histo_gg, histo_signal);

    hs->Add(histo_signal,"hist");
	hs->Add(histo_gg,"hist");
	hs->Add(histo_qq,"hist");
    hs->Add(histo_ttzwwz,"hist");
	hs->Add(histo_zx,"hist");

    /*hs->Add(histo_zx,"hist");
    hs->Add(histo_ttzwwz,"hist");
    hs->Add(histo_gg,"hist");
    hs->Add(histo_qq,"hist");
    hs->Add(histo_signal,"hist");*/

    TCanvas *c1 = new TCanvas("c1","example",800,1000);
    //gPad->SetLogy();
    c1->cd();

    hs->Draw("nostack");

    TLegend *legend = new TLegend(0.7,0.65,0.895,0.88);
	legend->SetTextFont(132);
	legend->SetTextSize(0.02);
	legend->AddEntry(histo_zx,"Z+X","f");
    legend->AddEntry(histo_ttzwwz,"t#bar{t}Z, WWZ","f");
	legend->AddEntry(histo_qq,"q#bar{q}#rightarrowZZ","f");
	legend->AddEntry(histo_gg,"gg#rightarrowZZ","f");
	legend->AddEntry(histo_signal,"VBS","f");
	legend->SetBorderSize(0);
    legend->Draw();

    hs->GetXaxis()->SetTitle("BDT response");
    hs->GetYaxis()->SetTitle("Event/bin");
    hs->GetYaxis()->SetTitleOffset(1.4);
    hs->GetXaxis()->SetTitleOffset(1.4);

    TCanvas *c2 = new TCanvas();
    c2->cd();
    histo_qq_weight->Draw();
 
 
    std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;
 
    delete reader;
 
    std::cout << "==> ApplyCut is done!" << std::endl << std::endl;
 }
 
 int main( int argc, char** argv )
 {
    TString methodList;
    for (int i=1; i<argc; i++) {
       TString regMethod(argv[i]);
       if(regMethod=="-b" || regMethod=="--batch") continue;
       if (!methodList.IsNull()) methodList += TString(",");
       methodList += regMethod;
    }
    ApplyCut(methodList);
    return 0;
 }