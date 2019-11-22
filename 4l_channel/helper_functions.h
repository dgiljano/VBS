double static e = TMath::E();

void create_electron_and_muon_objects(vector<TLorentzVector> &electrons, vector<TLorentzVector> &muons, vector<short> &electrons_charge, vector<short> &muons_charge, vector<short> *lepId, vector<float> *lepPt, vector<float> *lepEta, vector<float> *lepPhi)
{    
    for (int ilep = 0; ilep < lepId->size(); ilep++)
    {
        double eta = lepEta->at(ilep);
        double theta = 2 * TMath::ATan(pow(e,-lepEta->at(ilep)));
        double phi = lepPhi->at(ilep);
        double pt = lepPt->at(ilep);

        double px = pt * TMath::Cos(phi);
        double py = pt * TMath::Sin(phi);
        double pz = pt * TMath::SinH(eta);

        double E = TMath::Sqrt(px*px + py*py + pz*pz);
            
        if (abs(lepId->at(ilep)) == 11) //e
        {
            //cout << "Elektron" << endl;
            TLorentzVector el;
            el.SetPxPyPzE(px, py, pz, E);

            electrons.push_back(el);
            electrons_charge.push_back(lepId->at(ilep));
        }
        else //u
        {
            //cout << "Mion" << endl;
            TLorentzVector mu;
            mu.SetPxPyPzE(px, py, pz, E);

            muons.push_back(mu);
            muons_charge.push_back(lepId->at(ilep));
        }
    }
}

void build_ZZ_pair(vector<TLorentzVector> &electrons, vector<TLorentzVector> &muons, vector<short> &electrons_charge, vector<short> &muons_charge, TLorentzVector &Z1, TLorentzVector &Z2)
{
    TLorentzVector Z1_temp, Z2_temp;

    if (electrons.size() == 2 && muons.size() == 2) //2e+2u
    {
        if (electrons_charge.at(0) != electrons_charge.at(1)) Z1_temp = electrons.at(0) + electrons.at(1);
        if (muons_charge.at(0) != muons_charge.at(1)) Z2_temp = muons.at(0) + muons.at(1);

        // ------------------ choosing Z1 and Z2 --------------------

        if (abs(Z1_temp.M() - 91.1876) < abs(Z2_temp.M() - 91.1876))
        {
            Z1 = Z1_temp;
            Z2 = Z2_temp;
        }
        else
        {
            Z1 = Z2_temp;
            Z2 = Z1_temp;
        }
    }
    if (electrons.size() == 4) //4e
    {
        vector<int> i_ele; // electron indicees
        vector<int> i_pos; // positron indicees
        vector<TLorentzVector> Z_temp;
        
        for (int i = 0; i < 4; i++)
        {
            if (electrons_charge.at(i) == -11)
                i_ele.push_back(i);
            else 
                i_pos.push_back(i);
        }
        Z_temp.push_back(electrons.at(i_ele.at(0)) + electrons.at(i_pos.at(0)));
        Z_temp.push_back(electrons.at(i_ele.at(0)) + electrons.at(i_pos.at(1)));
        Z_temp.push_back(electrons.at(i_ele.at(1)) + electrons.at(i_pos.at(0)));
        Z_temp.push_back(electrons.at(i_ele.at(1)) + electrons.at(i_pos.at(1)));

        // ----------------- sort Zs by the closest to nominal Z mass -----------------

        for (int i = 0; i < 4; i++)
        {
            for (int j = i; j < 4; j++)
            {
                TLorentzVector temp;
                if (abs(Z_temp.at(i).M() - 91.1876) > abs(Z_temp.at(j).M() - 91.1876)) 
                {
                    temp = Z_temp.at(i);
                    Z_temp.at(i) = Z_temp.at(j);
                    Z_temp.at(j) = temp;
                }
            }
        }        
        Z1 = Z_temp.at(0);
        Z2 = Z_temp.at(1);
    }

    if (muons.size() == 4) //4u
    {
        vector<int> i_uminus; // u- indicees
        vector<int> i_uplus; // u+ indicees
        vector<TLorentzVector> Z_temp;
        
        for (int i = 0; i < 4; i++)
        {
            if (muons_charge.at(i) == -13)
                i_uminus.push_back(i);
            else 
                i_uplus.push_back(i);
        }
        Z_temp.push_back(muons.at(i_uminus.at(0)) + muons.at(i_uplus.at(0)));
        Z_temp.push_back(muons.at(i_uminus.at(0)) + muons.at(i_uplus.at(1)));
        Z_temp.push_back(muons.at(i_uminus.at(1)) + muons.at(i_uplus.at(0)));
        Z_temp.push_back(muons.at(i_uminus.at(1)) + muons.at(i_uplus.at(1)));

        // ----------------- sort Zs by the closest to nominal Z mass -----------------

        for (int i = 0; i < 4; i++)
        {
            for (int j = i; j < 4; j++)
            {
                TLorentzVector temp;
                if (abs(Z_temp.at(i).M() - 91.1876) > abs(Z_temp.at(j).M() - 91.1876)) 
                {
                    temp = Z_temp.at(i);
                    Z_temp.at(i) = Z_temp.at(j);
                    Z_temp.at(j) = temp;
                }
            }
        }

        Z1 = Z_temp.at(0);
        Z2 = Z_temp.at(1);
    }
}

void clear_vectors(vector<TLorentzVector> &electrons, vector<TLorentzVector> &muons, vector<short> &electrons_charge, vector<short> &muons_charge)
{
    electrons.clear();
    muons.clear();
    electrons_charge.clear();
    muons_charge.clear();
}

void calculate_Zeppenfeld_Z(TLorentzVector &Z1, TLorentzVector &Z2, vector<float> *JetEta, float &eta_Z1_star, float &eta_Z2_star)
{
    eta_Z1_star = Z1.Eta() - (JetEta->at(0) + JetEta->at(1))/2;
    eta_Z2_star = Z2.Eta() - (JetEta->at(0) + JetEta->at(1))/2;
}

void calculate_R_pt_hard(TLorentzVector &Z1, TLorentzVector &Z2, vector<float> *JetEta, vector<float> *JetPhi, vector<float> *JetPt, float &R_pt_hard)
{
    float j1_px = JetPt->at(0) * TMath::Cos(JetPhi->at(0));
    float j1_py = JetPt->at(0) * TMath::Sin(JetPhi->at(0));
    float j1_pz = JetPt->at(0) * TMath::SinH(JetEta->at(0));

    float j2_px = JetPt->at(1) * TMath::Cos(JetPhi->at(1));
    float j2_py = JetPt->at(1) * TMath::Sin(JetPhi->at(1));
    float j2_pz = JetPt->at(1) * TMath::SinH(JetEta->at(1));

    float vx = j1_px + j2_px + Z1.Px() + Z2.Px();
    float vy = j1_py + j2_py + Z1.Py() + Z2.Py();
    float vz = j1_pz + j2_pz + Z1.Pz() + Z2.Pz();

    TVector3 v_hard(vx, vy, vz);
    float mag_hard = JetPt->at(0) + JetPt->at(1) + Z1.Pt() + Z2.Pt();

    R_pt_hard = v_hard.Perp()/mag_hard;
}

void calculate_R_pt_jet(vector<float> *JetEta, vector<float> *JetPhi, vector<float> *JetPt, float &R_pt_jet)
{
    float j1_px = JetPt->at(0) * TMath::Cos(JetPhi->at(0));
    float j1_py = JetPt->at(0) * TMath::Sin(JetPhi->at(0));
    float j1_pz = JetPt->at(0) * TMath::SinH(JetEta->at(0));

    float j2_px = JetPt->at(1) * TMath::Cos(JetPhi->at(1));
    float j2_py = JetPt->at(1) * TMath::Sin(JetPhi->at(1));
    float j2_pz = JetPt->at(1) * TMath::SinH(JetEta->at(1));

    float vx = j1_px + j2_px;
    float vy = j1_py + j2_py;
    float vz = j1_pz + j2_pz;

    TVector3 v_jet(vx, vy, vz);
    float mag_jet = JetPt->at(0) + JetPt->at(1);

    R_pt_jet = v_jet.Perp()/mag_jet;
}