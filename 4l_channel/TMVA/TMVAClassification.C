// @(#)root/tmva $Id$
/**********************************************************************************
 * Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)                     *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set of classifiers is used.                      *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 * Launch the GUI via the command:                                                *
 *                                                                                *
 *    root -l ./TMVAGui.C                                                         *
 *                                                                                *
 * You can also compile and run the example with the following commands           *
 *                                                                                *
 *    make                                                                        *
 *    ./TMVAClassification <Methods>                                              *
 *                                                                                *
 * where: <Methods> = "method1 method2"                                           *
 *        are the TMVA classifier names                                           *
 *                                                                                *
 * example:                                                                       *
 *    ./TMVAClassification Fisher LikelihoodPCA BDT                               *
 *                                                                                *
 * If no method given, a default set is of classifiers is used                    *
 **********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/DataLoader.h"

int TMVAClassification( TString myMethodList = "" )
{
   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the
   // corresponding lines from .rootrc

   // methods to be processed can be given as an argument; use format:
   //
   // mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
   //
   // if you like to use a method via the plugin mechanism, we recommend using
   //
   // mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
   // (an example is given for using the BDT as plugin (see below),
   // but of course the real application is when you write your own
   // method based)

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Cut optimisation
   Use["Cuts"]            = 1;
   Use["CutsD"]           = 1;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // 
   // --- 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 1;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 1;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 1;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 1; // k-nearest neighbour method
   //
   // --- Linear Discriminant Analysis
   Use["LD"]              = 1; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 1; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   //
   // --- Support Vector Machine 
   Use["SVM"]             = 1;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost          // TO KORISTIN
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 1;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return 1;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVA.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is 
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
                                               

      TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");

   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   dataloader->AddVariable( "mjj",                "mjj", "units", 'F' );
   dataloader->AddVariable( "deta_jj",            "deta_jj", "units", 'F' );
   dataloader->AddVariable( "m4l",                "m4l", "units", 'F' );
   dataloader->AddVariable( "eta_star_Z1",        "eta_star_Z1", "units", 'F' );
   dataloader->AddVariable( "eta_star_Z2",        "eta_star_Z2", "units", 'F' );
   dataloader->AddVariable( "R_pt_hard",          "R_pt_hard", "units", 'F' );
   dataloader->AddVariable( "R_pt_jet",           "R_pt_jet", "units", 'F' );

   /*dataloader->AddVariable( "abs_etajet_min",           "abs_etajet_min", "units", 'F' );
   dataloader->AddVariable( "abs_etajet_max",           "abs_etajet_max", "units", 'F' );
   dataloader->AddVariable( "abs_etalep_min",           "abs_etalep_min", "units", 'F' );
   dataloader->AddVariable( "abs_etalep_max",           "abs_etalep_max", "units", 'F' );
   dataloader->AddVariable( "delta_phi_ZZ",           "delta_phi_ZZ", "units", 'F' );
   dataloader->AddVariable( "rapidity_Z1",           "rapidity_Z1", "units", 'F' );
   dataloader->AddVariable( "rapidity_Z2",           "rapidity_Z2", "units", 'F' );
   dataloader->AddVariable( "rapidity_j1",           "rapidity_j1", "units", 'F' );
   dataloader->AddVariable( "rapidity_j2",           "rapidity_j2", "units", 'F' );
   dataloader->AddVariable( "pt_Z1",           "pt_Z1", "units", 'F' );
   dataloader->AddVariable( "pt_Z2",           "pt_Z2", "units", 'F' );
   dataloader->AddVariable( "pt_l3",           "pt_l3", "units", 'F' );
   dataloader->AddVariable( "jet1_qg_tagger",           "jet1_qg_tagger", "units", 'F' );
   dataloader->AddVariable( "jet2_qg_tagger",           "jet2_qg_tagger", "units", 'F' );
   //dataloader->AddVariable( "dbkg_kin",           "dbkg_kin", "units", 'F' );
   dataloader->AddVariable( "eta_j1",           "eta_j1", "units", 'F' );
   dataloader->AddVariable( "eta_j2",           "eta_j2", "units", 'F' );
   dataloader->AddVariable( "pt_jet1",           "pt_jet1", "units", 'F' );
   dataloader->AddVariable( "pt_jet2",           "pt_jet2", "units", 'F' );
   dataloader->AddVariable( "eta_j_sum",           "eta_j_sum", "units", 'F' );
   dataloader->AddVariable( "mjj_over_detajj",           "mjj_over_detajj", "units", 'F' );
   dataloader->AddVariable( "abs_etajet_sum",           "abs_etajet_sum", "units", 'F' );*/

   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   //factory->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
   //factory->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );

   // Read training and test data
   // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
   
   TFile *input = TFile::Open( "MVA.root" );
   
   std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;
   
   // --- Register the training and test trees


   TTree *signal     = (TTree*)input->Get("signal");              // MVA on gen samples
   TTree *background = (TTree*)input->Get("qq");          // MVA on gen samples

   
   // global event weights per tree (see below for setting event-wise weights)
   //Float_t signalWeight     = 1.0;
   //Float_t backgroundWeight = 1.0;
   
   // You can add an arbitrary number of signal or background trees
   //dataloader->AddSignalTree    ( signal,     signalWeight     );
   //dataloader->AddBackgroundTree( background, backgroundWeight );



   std::vector<Double_t> vars(7);
   //Signal
   float signal_mjj, signal_deta_jj, signal_m4l, signal_eta_star_Z1, signal_eta_star_Z2, signal_R_pt_hard, signal_R_pt_jet;
   float signal_abs_etajet_min, signal_abs_etajet_max, signal_abs_etalep_min, signal_abs_etalep_max, signal_delta_phi_ZZ, signal_rapidity_Z1, signal_rapidity_Z2, signal_rapidity_j1, signal_rapidity_j2, signal_pt_Z1, signal_pt_Z2, signal_pt_l3, signal_jet1_qg_tagger, signal_jet2_qg_tagger, signal_eta_j1, signal_eta_j2, signal_pt_jet1, signal_pt_jet2, signal_eta_j_sum, signal_mjj_over_detajj, signal_abs_etajet_sum;
   float signal_weight;

   signal->SetBranchAddress("mjj", &signal_mjj);
   signal->SetBranchAddress("deta_jj", &signal_deta_jj);
   signal->SetBranchAddress("m4l", &signal_m4l);
   signal->SetBranchAddress("eta_star_Z1", &signal_eta_star_Z1);
   signal->SetBranchAddress("eta_star_Z2", &signal_eta_star_Z2);
   signal->SetBranchAddress("R_pt_hard", &signal_R_pt_hard);
   signal->SetBranchAddress("R_pt_jet", &signal_R_pt_jet);
   signal->SetBranchAddress("weight", &signal_weight);

   /*signal->SetBranchAddress( "abs_etajet_min", &signal_abs_etajet_min);
   signal->SetBranchAddress( "abs_etajet_max", &signal_abs_etajet_max);
   signal->SetBranchAddress( "abs_etalep_min", &signal_abs_etalep_min);
   signal->SetBranchAddress( "abs_etalep_max", &signal_abs_etalep_max);
   signal->SetBranchAddress( "delta_phi_ZZ", &signal_delta_phi_ZZ);
   signal->SetBranchAddress( "rapidity_Z1", &signal_rapidity_Z1);
   signal->SetBranchAddress( "rapidity_Z2", &signal_rapidity_Z2);
   signal->SetBranchAddress( "rapidity_j1", &signal_rapidity_j1);
   signal->SetBranchAddress( "rapidity_j2", &signal_rapidity_j2);
   signal->SetBranchAddress( "pt_Z1", &signal_pt_Z1);
   signal->SetBranchAddress( "pt_Z2", &signal_pt_Z2);
   signal->SetBranchAddress( "pt_l3", &signal_pt_l3);
   signal->SetBranchAddress( "jet1_qg_tagger", &signal_jet1_qg_tagger);
   signal->SetBranchAddress( "jet2_qg_tagger", &signal_jet2_qg_tagger);
   signal->SetBranchAddress( "eta_j1", &signal_eta_j1);
   signal->SetBranchAddress( "eta_j2", &signal_eta_j2);
   signal->SetBranchAddress( "pt_jet1", &signal_pt_jet1);
   signal->SetBranchAddress( "pt_jet2", &signal_pt_jet2);
   signal->SetBranchAddress( "eta_j_sum", &signal_eta_j_sum);
   signal->SetBranchAddress( "mjj_over_detajj", &signal_mjj_over_detajj);
   signal->SetBranchAddress( "abs_etajet_sum", &signal_abs_etajet_sum);*/


   for (int i = 0; i < signal->GetEntries(); i++)
   {
      signal->GetEntry(i);

      vars[0] = signal_mjj;
      vars[1] = signal_deta_jj;
      vars[2] = signal_m4l;
      vars[3] = signal_eta_star_Z1;
      vars[4] = signal_eta_star_Z2;
      vars[5] = signal_R_pt_hard;
      vars[6] = signal_R_pt_jet;

      /*vars[7] = signal_abs_etajet_min;
      vars[8] = signal_abs_etajet_max;
      vars[9] = signal_abs_etalep_min;
      vars[10] = signal_abs_etalep_max;
      vars[11] = signal_delta_phi_ZZ;
      vars[12] = signal_rapidity_Z1;
      vars[13] = signal_rapidity_Z2;
      vars[14] = signal_rapidity_j1;
      vars[15] = signal_rapidity_j2;
      vars[16] = signal_pt_Z1;
      vars[17] = signal_pt_Z2;
      vars[18] = signal_pt_l3;
      vars[19] = signal_jet1_qg_tagger;
      vars[20] = signal_jet2_qg_tagger;
      vars[21] = signal_eta_j1;
      vars[22] = signal_eta_j2;
      vars[23] = signal_pt_jet1;
      vars[24] = signal_pt_jet2;
      vars[25] = signal_eta_j_sum;
      vars[26] = signal_mjj_over_detajj;
      vars[27] = signal_abs_etajet_sum;*/

      //if (signal_weight < 0) continue;
      if (i < signal->GetEntries()/2.0) dataloader->AddSignalTrainingEvent( vars, signal_weight );
      else                              dataloader->AddSignalTestEvent    ( vars, signal_weight );

   }

   //Background
   float background_mjj, background_deta_jj, background_m4l, background_eta_star_Z1, background_eta_star_Z2, background_R_pt_hard, background_R_pt_jet;
   float background_abs_etajet_min, background_abs_etajet_max, background_abs_etalep_min, background_abs_etalep_max, background_delta_phi_ZZ, background_rapidity_Z1, background_rapidity_Z2, background_rapidity_j1, background_rapidity_j2, background_pt_Z1, background_pt_Z2, background_pt_l3, background_jet1_qg_tagger, background_jet2_qg_tagger, background_eta_j1, background_eta_j2, background_pt_jet1, background_pt_jet2, background_eta_j_sum, background_mjj_over_detajj, background_abs_etajet_sum;
   float background_weight;

   background->SetBranchAddress("mjj", &background_mjj);
   background->SetBranchAddress("deta_jj", &background_deta_jj);
   background->SetBranchAddress("m4l", &background_m4l);
   background->SetBranchAddress("eta_star_Z1", &background_eta_star_Z1);
   background->SetBranchAddress("eta_star_Z2", &background_eta_star_Z2);
   background->SetBranchAddress("R_pt_hard", &background_R_pt_hard);
   background->SetBranchAddress("R_pt_jet", &background_R_pt_jet);
   background->SetBranchAddress("weight", &background_weight);

   /*background->SetBranchAddress( "abs_etajet_min", &background_abs_etajet_min);
   background->SetBranchAddress( "abs_etajet_max", &background_abs_etajet_max);
   background->SetBranchAddress( "abs_etalep_min", &background_abs_etalep_min);
   background->SetBranchAddress( "abs_etalep_max", &background_abs_etalep_max);
   background->SetBranchAddress( "delta_phi_ZZ", &background_delta_phi_ZZ);
   background->SetBranchAddress( "rapidity_Z1", &background_rapidity_Z1);
   background->SetBranchAddress( "rapidity_Z2", &background_rapidity_Z2);
   background->SetBranchAddress( "rapidity_j1", &background_rapidity_j1);
   background->SetBranchAddress( "rapidity_j2", &background_rapidity_j2);
   background->SetBranchAddress( "pt_Z1", &background_pt_Z1);
   background->SetBranchAddress( "pt_Z2", &background_pt_Z2);
   background->SetBranchAddress( "pt_l3", &background_pt_l3);
   background->SetBranchAddress( "jet1_qg_tagger", &background_jet1_qg_tagger);
   background->SetBranchAddress( "jet2_qg_tagger", &background_jet2_qg_tagger);
   background->SetBranchAddress( "eta_j1", &background_eta_j1);
   background->SetBranchAddress( "eta_j2", &background_eta_j2);
   background->SetBranchAddress( "pt_jet1", &background_pt_jet1);
   background->SetBranchAddress( "pt_jet2", &background_pt_jet2);
   background->SetBranchAddress( "eta_j_sum", &background_eta_j_sum);
   background->SetBranchAddress( "mjj_over_detajj", &background_mjj_over_detajj);
   background->SetBranchAddress( "abs_etajet_sum", &background_abs_etajet_sum);*/

   for (int i = 0; i < background->GetEntries(); i++)
   {
      background->GetEntry(i);

      vars[0] = background_mjj;
      vars[1] = background_deta_jj;
      vars[2] = background_m4l;
      vars[3] = background_eta_star_Z1;
      vars[4] = background_eta_star_Z2;
      vars[5] = background_R_pt_hard;
      vars[6] = background_R_pt_jet;

      /*vars[7] = background_abs_etajet_min;
      vars[8] = background_abs_etajet_max;
      vars[9] = background_abs_etalep_min;
      vars[10] = background_abs_etalep_max;
      vars[11] = background_delta_phi_ZZ;
      vars[12] = background_rapidity_Z1;
      vars[13] = background_rapidity_Z2;
      vars[14] = background_rapidity_j1;
      vars[15] = background_rapidity_j2;
      vars[16] = background_pt_Z1;
      vars[17] = background_pt_Z2;
      vars[18] = background_pt_l3;
      vars[19] = background_jet1_qg_tagger;
      vars[20] = background_jet2_qg_tagger;
      vars[21] = background_eta_j1;
      vars[22] = background_eta_j2;
      vars[23] = background_pt_jet1;
      vars[24] = background_pt_jet2;
      vars[25] = background_eta_j_sum;
      vars[26] = background_mjj_over_detajj;
      vars[27] = background_abs_etajet_sum;*/

      //if (background_weight < 0) continue;
      if (i < background->GetEntries()/2.0) dataloader->AddBackgroundTrainingEvent( vars, background_weight );
      else                                  dataloader->AddBackgroundTestEvent    ( vars, background_weight );
   }

   
   // To give different trees for training and testing, do as follows:
   //    factory->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
   //    factory->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );
   
   // Use the following code instead of the above two or four lines to add signal and background
   // training and test events "by hand"
   // NOTE that in this case one should not give expressions (such as "var1+var2") in the input
   //      variable definition, but simply compute the expression before adding the event
   //
   //     // --- begin ----------------------------------------------------------
   //     std::vector<Double_t> vars( 4 ); // vector has size of number of input variables
   //     Float_t  treevars[4], weight;
   //     
   //     // Signal
   //     for (UInt_t ivar=0; ivar<4; ivar++) signal->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   //     for (UInt_t i=0; i<signal->GetEntries(); i++) {
   //        signal->GetEntry(i);
   //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //        // add training and test events; here: first half is training, second is testing
   //        // note that the weight can also be event-wise
   //        if (i < signal->GetEntries()/2.0) factory->AddSignalTrainingEvent( vars, signalWeight );
   //        else                              factory->AddSignalTestEvent    ( vars, signalWeight );
   //     }
   //   
   //     // Background (has event weights)
   //     background->SetBranchAddress( "weight", &weight );
   //     for (UInt_t ivar=0; ivar<4; ivar++) background->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   //     for (UInt_t i=0; i<background->GetEntries(); i++) {
   //        background->GetEntry(i);
   //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //        // add training and test events; here: first half is training, second is testing
   //        // note that the weight can also be event-wise
   //        if (i < background->GetEntries()/2) factory->AddBackgroundTrainingEvent( vars, backgroundWeight*weight );
   //        else                                factory->AddBackgroundTestEvent    ( vars, backgroundWeight*weight );
   //     }
         // --- end ------------------------------------------------------------
   //
   // --- end of tree registration 

   // Set individual event weights (the variables must exist in the original TTree)
   //    for signal    : factory->SetSignalWeightExpression    ("weight1*weight2");
   //    for background: factory->SetBackgroundWeightExpression("weight1*weight2");
   //factory->SetBackgroundWeightExpression( "weight" );

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

   // Tell the factory how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used 
   // for training, and the other half for testing:
   //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
   // To also specify the number of testing events, use:
   //    factory->PrepareTrainingAndTestTree( mycut,
   //                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

   // ---- Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Cut optimisation
   if (Use["Cuts"])
      factory->BookMethod(dataloader, TMVA::Types::kCuts, "Cuts",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

   if (Use["CutsD"])
      factory->BookMethod(dataloader, TMVA::Types::kCuts, "CutsD",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

   if (Use["CutsPCA"])
      factory->BookMethod(dataloader, TMVA::Types::kCuts, "CutsPCA",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

   if (Use["CutsGA"])
      factory->BookMethod(dataloader, TMVA::Types::kCuts, "CutsGA",
                           "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

   if (Use["CutsSA"])
      factory->BookMethod(dataloader, TMVA::Types::kCuts, "CutsSA",
                           "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   // Likelihood ("naive Bayes estimator")
   if (Use["Likelihood"])
      factory->BookMethod(dataloader, TMVA::Types::kLikelihood, "Likelihood",
                           "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

   // Decorrelated likelihood
   if (Use["LikelihoodD"])
      factory->BookMethod(dataloader, TMVA::Types::kLikelihood, "LikelihoodD",
                           "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );

   // PCA-transformed likelihood
   if (Use["LikelihoodPCA"])
      factory->BookMethod(dataloader, TMVA::Types::kLikelihood, "LikelihoodPCA",
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 

   // Use a kernel density estimator to approximate the PDFs
   if (Use["LikelihoodKDE"])
      factory->BookMethod(dataloader, TMVA::Types::kLikelihood, "LikelihoodKDE",
                           "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

   // Use a variable-dependent mix of splines and kernel density estimator
   if (Use["LikelihoodMIX"])
      factory->BookMethod(dataloader, TMVA::Types::kLikelihood, "LikelihoodMIX",
                           "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

   // Test the multi-dimensional probability density estimator
   // here are the options strings for the MinMax and RMS methods, respectively:
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
   if (Use["PDERS"])
      factory->BookMethod(dataloader, TMVA::Types::kPDERS, "PDERS",
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

   if (Use["PDERSD"])
      factory->BookMethod(dataloader, TMVA::Types::kPDERS, "PDERSD",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

   if (Use["PDERSPCA"])
      factory->BookMethod(dataloader, TMVA::Types::kPDERS, "PDERSPCA",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

   // Multi-dimensional likelihood estimator using self-adapting phase-space binning
   if (Use["PDEFoam"])
      factory->BookMethod(dataloader, TMVA::Types::kPDEFoam, "PDEFoam",
                           "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

   if (Use["PDEFoamBoost"])
      factory->BookMethod(dataloader, TMVA::Types::kPDEFoam, "PDEFoamBoost",
                           "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"])
      factory->BookMethod(dataloader, TMVA::Types::kKNN, "KNN",
                           "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

   // H-Matrix (chi2-squared) method
   if (Use["HMatrix"])
      factory->BookMethod(dataloader, TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );

   // Linear discriminant (same as Fisher discriminant)
   if (Use["LD"])
      factory->BookMethod(dataloader, TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher discriminant (same as LD)
   if (Use["Fisher"])
      factory->BookMethod(dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher with Gauss-transformed input variables
   if (Use["FisherG"])
      factory->BookMethod(dataloader, TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

   // Composite classifier: ensemble (tree) of boosted Fisher classifiers
   if (Use["BoostedFisher"])
      factory->BookMethod(dataloader, TMVA::Types::kFisher, "BoostedFisher", 
                           "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );

   // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"])
      factory->BookMethod(dataloader, TMVA::Types::kFDA, "FDA_MC",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

   if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod(dataloader, TMVA::Types::kFDA, "FDA_GA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

   if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod(dataloader, TMVA::Types::kFDA, "FDA_SA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   if (Use["FDA_MT"])
      factory->BookMethod(dataloader, TMVA::Types::kFDA, "FDA_MT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

   if (Use["FDA_GAMT"])
      factory->BookMethod(dataloader, TMVA::Types::kFDA, "FDA_GAMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

   if (Use["FDA_MCMT"])
      factory->BookMethod(dataloader, TMVA::Types::kFDA, "FDA_MCMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   if (Use["MLP"])
      factory->BookMethod(dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

   if (Use["MLPBFGS"])
      factory->BookMethod(dataloader, TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

   if (Use["MLPBNN"])
      factory->BookMethod(dataloader, TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

   // CF(Clermont-Ferrand)ANN
   if (Use["CFMlpANN"])
      factory->BookMethod(dataloader, TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

   // Tmlp(Root)ANN
   if (Use["TMlpANN"])
      factory->BookMethod(dataloader, TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

   // Support Vector Machine
   if (Use["SVM"])
      factory->BookMethod(dataloader, TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

   if (Use["BDTB"]) // Bagging
      factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );

   if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
      factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTMitFisher",
                           "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"])
      factory->BookMethod(dataloader, TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

   // For an example of the category classifier usage, see: TMVAClassificationCategory

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

   // ---- STILL EXPERIMENTAL and only implemented for BDT's ! 
   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","FitGA");

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

   return 0;
}

int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList; 
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(","); 
      methodList += regMethod;
   }
   return TMVAClassification(methodList); 
}
