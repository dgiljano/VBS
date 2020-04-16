void prepare_one_bin()
{
    TFile *file = new TFile("./all_contributions_2016.root");

    TFile *weight_file = new TFile("/home/llr/cms/giljanovic/scratch/VBS/CMSSW_10_2_15/src/vbs_analysis/4l_channel/TMVA/MVA.root");
    TTree *tree = (TTree*)weight_file->Get("signal");
    Float_t weight;
    tree->SetBranchAddress("weight", &weight);

    TH1F *diboson_original = (TH1F*)file->Get("diboson");
    Int_t n_entries_histo = diboson_original->GetEntries();

    TH1F *diboson = new TH1F("diboson","diboson", 1, 0, 1);
    diboson->Sumw2();

    Float_t sum = 0;
    for (int i = 0; i < tree->GetEntries(); i++)
    {
        tree->GetEntry(i);

        //cout << i+1 << "\t" << weight << endl;
        sum += weight;
        diboson->Fill(0.5, weight);
    }

    cout << "suma " << sum << endl;
    diboson->Draw("E1");
}