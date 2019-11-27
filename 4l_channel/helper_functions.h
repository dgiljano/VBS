double static e = TMath::E();
const double PI  =3.141592653589793238463;

void create_electron_and_muon_objects(vector<TLorentzVector> &electrons, vector<TLorentzVector> &muons, vector<short> &electrons_charge, vector<short> &muons_charge, vector<short> *lepId, vector<float> *lepPt, vector<float> *lepEta, vector<float> *lepPhi)
{    
    for (int ilep = 0; ilep < lepId->size(); ilep++)
    {
        float eta = lepEta->at(ilep);
        float theta = 2 * TMath::ATan(pow(e,-lepEta->at(ilep)));
        float phi = lepPhi->at(ilep);
        float pt = lepPt->at(ilep);

        float px = pt * TMath::Cos(phi);
        float py = pt * TMath::Sin(phi);
        float pz = pt * TMath::SinH(eta);

        float E = TMath::Sqrt(px*px + py*py + pz*pz);
            
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

void calculate_min_max_jet_eta(vector<float> *JetEta, float &abs_etajet_min, float &abs_etajet_max)
{
    vector<float> *JetEta_ordered = JetEta;

    // --------------- sort |jet eta| --------------------
    for (int i = 0; i < JetEta_ordered->size(); i++)
    {
        for (int j = i; j < JetEta_ordered->size(); j++)
        {
            float temp;
            if (abs(JetEta_ordered->at(i)) < abs(JetEta_ordered->at(j))) 
            {
                temp = JetEta_ordered->at(i);
                JetEta_ordered->at(i) = JetEta_ordered->at(j);
                JetEta_ordered->at(j) = temp;
            }
        }
    }
    //--------------------------------------------------

    abs_etajet_min = abs(JetEta_ordered->at(JetEta_ordered->size()-1));
    abs_etajet_max = abs(JetEta_ordered->at(0));
}

void calculate_min_max_lepton_eta(vector<TLorentzVector> &electrons, vector<TLorentzVector> &muons, float &abs_etalep_min, float &abs_etalep_max)
{
    vector<TLorentzVector> leptons_ordered;

    for (int i = 0; i < electrons.size(); i++)
    {
        leptons_ordered.push_back(electrons.at(i));
    }
    for (int i = 0; i < muons.size(); i++)
    {
        leptons_ordered.push_back(muons.at(i));
    }

    // --------------- sort |lepton eta| --------------------
    for (int i = 0; i < leptons_ordered.size(); i++)
    {
        for (int j = i; j < leptons_ordered.size(); j++)
        {
            TLorentzVector temp;
            if (abs(leptons_ordered.at(i).Eta()) < abs(leptons_ordered.at(j).Eta())) 
            {
                temp = leptons_ordered.at(i);
                leptons_ordered.at(i) = leptons_ordered.at(j);
                leptons_ordered.at(j) = temp;
            }
        }
    }
    //--------------------------------------------------

    abs_etalep_min = abs(leptons_ordered.at(leptons_ordered.size()-1).Eta());
    abs_etalep_max = abs(leptons_ordered.at(0).Eta());
}

void calculate_dphi_ZZ(TLorentzVector &Z1, TLorentzVector &Z2, float &delta_phi_ZZ)
{
    if (Z1.Phi() >= 0 && Z2.Phi() >= 0)
    {
        delta_phi_ZZ = TMath::Abs(Z1.Phi() - Z2.Phi());
    }
    else if (Z1.Phi() <= 0 && Z2.Phi() <= 0)
    {
        delta_phi_ZZ = TMath::Abs(Z1.Phi() - Z2.Phi());
    }              
    else
    {
        if (Z1.Phi() < 0)
        {
            float positive_phi = Z1.Phi() + 2*PI;

            float dphi = positive_phi - Z2.Phi();

            if (dphi > PI) dphi = 2*PI - dphi;

            delta_phi_ZZ = dphi;
        }
        else
        {
            float positive_phi = Z2.Phi() + 2*PI;

            float dphi = positive_phi - Z1.Phi();

            if (dphi > PI) dphi = 2*PI - dphi;

            delta_phi_ZZ = dphi;
        }
    }
}

void calculate_rapidity_Z1_Z2(TLorentzVector &Z1, TLorentzVector &Z2, float &rapidity_Z1, float &rapidity_Z2)
{
    rapidity_Z1 = Z1.Rapidity();
    rapidity_Z2 = Z2.Rapidity();
}

void calculate_rapidity_j1_j2(vector<float> *JetEta, vector<float> *JetPhi, vector<float> *JetPt, float &rapidity_j1, float &rapidity_j2)
{
    float j1_px = JetPt->at(0) * TMath::Cos(JetPhi->at(0));
    float j1_py = JetPt->at(0) * TMath::Sin(JetPhi->at(0));
    float j1_pz = JetPt->at(0) * TMath::SinH(JetEta->at(0));
    float E_j1 = TMath::Sqrt(j1_px*j1_px + j1_py*j1_py + j1_pz*j1_pz);

    float j2_px = JetPt->at(1) * TMath::Cos(JetPhi->at(1));
    float j2_py = JetPt->at(1) * TMath::Sin(JetPhi->at(1));
    float j2_pz = JetPt->at(1) * TMath::SinH(JetEta->at(1));
    float E_j2 = TMath::Sqrt(j2_px*j2_px + j2_py*j2_py + j2_pz*j2_pz);

    rapidity_j1 = 0.5 * TMath::Log((E_j1 + j1_pz)/(E_j1 - j1_pz));
    rapidity_j2 = 0.5 * TMath::Log((E_j2 + j2_pz)/(E_j2 - j2_pz));
}

void calculate_pt_Z1_Z2_l3(TLorentzVector &Z1, TLorentzVector &Z2, vector<float> *lepPt, float &pt_Z1, float &pt_Z2, float &pt_l3)
{
    pt_Z1 = Z1.Pt();
    pt_Z2 = Z2.Pt();
    pt_l3 = lepPt->at(2);
}

// overloaded for Z+X
void build_ZZ_pair(vector<TLorentzVector> &electrons, vector<TLorentzVector> &muons, vector<short> &electrons_charge, vector<short> &muons_charge, TLorentzVector &Z1, TLorentzVector &Z2, bool &ZX)
{
    TLorentzVector Z1_temp, Z2_temp;

    if (electrons.size() == 2 && muons.size() == 2) //2e+2u
    {
        if (electrons_charge.at(0) != electrons_charge.at(1)) 
        {
            Z1 = electrons.at(0) + electrons.at(1);
            Z2 = muons.at(0) + muons.at(1);
        }
        else
        {
            Z1 = muons.at(0) + muons.at(1);
            Z2 = electrons.at(0) + electrons.at(1);
        }
    }
    if (electrons.size() == 4) //4e
    {
        vector<TLorentzVector> Z_temp;

        vector<TLorentzVector> all_ele;
        vector<TLorentzVector> all_pos;
        vector<TLorentzVector> Z1_candidates;
        vector<TLorentzVector> Z2_candidates;
        
        for (int i = 0; i < 4; i++)
        {
            if (electrons_charge.at(i) == -11)
                all_ele.push_back(electrons.at(i));
            else 
                all_pos.push_back(electrons.at(i));
        }

        if (all_ele.size() == 1) // I have 1 electron and 3 positrons
        {
            if ((all_ele.at(0) + all_pos.at(0)).M() > 12 && (all_ele.at(0) + all_pos.at(0)).M() < 120) // I have to check if I can try to combine these two into Z: I have to have 12 < mll < 120
            {
                Z1_candidates.push_back(all_ele.at(0) + all_pos.at(0));
                Z2_candidates.push_back(all_pos.at(1) + all_pos.at(2));
            }
            if ((all_ele.at(0) + all_pos.at(1)).M() > 12 && (all_ele.at(0) + all_pos.at(1)).M() < 120)
            {
                Z1_candidates.push_back(all_ele.at(0) + all_pos.at(1));
                Z2_candidates.push_back(all_pos.at(0) + all_pos.at(2));
            }
            if ((all_ele.at(0) + all_pos.at(2)).M() > 12 && (all_ele.at(0) + all_pos.at(2)).M() < 120)
            {
                Z1_candidates.push_back(all_ele.at(0) + all_pos.at(2));
                Z2_candidates.push_back(all_pos.at(0) + all_pos.at(1));
            }

            // ----------------- pairing Zs ------------------

            if (Z1_candidates.size() == 1)
            {
                Z1 = Z1_candidates.at(0);
                Z2 = Z2_candidates.at(0);
            }

            for (int i = 0; i < Z1_candidates.size(); i++)
            {
                for (int j = i; j < Z1_candidates.size(); j++)
                {
                    TLorentzVector temp_Z1, temp_Z2;
                    if (abs(Z1_candidates.at(i).M() - 91.1876) > abs(Z1_candidates.at(j).M() - 91.1876)) 
                    {
                        temp_Z1 = Z1_candidates.at(i);
                        Z1_candidates.at(i) = Z1_candidates.at(j);
                        Z1_candidates.at(j) = temp_Z1;

                        temp_Z2 = Z2_candidates.at(i);
                        Z2_candidates.at(i) = Z2_candidates.at(j);
                        Z2_candidates.at(j) = temp_Z2;
                    }
                }
            }

            Z1 = Z1_candidates.at(0);
            Z2 = Z2_candidates.at(0);
        }
        if (all_ele.size() == 3) // I have 3 electrons and 1 positron
        {
            if ((all_pos.at(0) + all_ele.at(0)).M() > 12 && (all_pos.at(0) + all_ele.at(0)).M() < 120) // I have to check if I can try to combine these two into Z: I have to have 12 < mll < 120
            {
                Z1_candidates.push_back(all_pos.at(0) + all_ele.at(0));
                Z2_candidates.push_back(all_ele.at(1) + all_ele.at(2));
            }
            if ((all_pos.at(0) + all_ele.at(1)).M() > 12 && (all_pos.at(0) + all_ele.at(1)).M() < 120)
            {
                Z1_candidates.push_back(all_pos.at(0) + all_ele.at(1));
                Z2_candidates.push_back(all_ele.at(0) + all_ele.at(2));
            }
            if ((all_pos.at(0) + all_ele.at(2)).M() > 12 && (all_pos.at(0) + all_ele.at(2)).M() < 120)
            {
                Z1_candidates.push_back(all_pos.at(0) + all_ele.at(2));
                Z2_candidates.push_back(all_ele.at(0) + all_ele.at(1));
            }

            // ----------------- pairing Zs ------------------

            if (Z1_candidates.size() == 1)
            {
                Z1 = Z1_candidates.at(0);
                Z2 = Z2_candidates.at(0);
            }

            for (int i = 0; i < Z1_candidates.size(); i++)
            {
                for (int j = i; j < Z1_candidates.size(); j++)
                {
                    TLorentzVector temp_Z1, temp_Z2;
                    if (abs(Z1_candidates.at(i).M() - 91.1876) > abs(Z1_candidates.at(j).M() - 91.1876)) 
                    {
                        temp_Z1 = Z1_candidates.at(i);
                        Z1_candidates.at(i) = Z1_candidates.at(j);
                        Z1_candidates.at(j) = temp_Z1;

                        temp_Z2 = Z2_candidates.at(i);
                        Z2_candidates.at(i) = Z2_candidates.at(j);
                        Z2_candidates.at(j) = temp_Z2;
                    }
                }
            }

            Z1 = Z1_candidates.at(0);
            Z2 = Z2_candidates.at(0);
        }
    }
    if (muons.size() == 4) //4u
    {
        vector<TLorentzVector> Z_temp;

        vector<TLorentzVector> all_uminus;
        vector<TLorentzVector> all_uplus;
        vector<TLorentzVector> Z1_candidates;
        vector<TLorentzVector> Z2_candidates;
        
        for (int i = 0; i < 4; i++)
        {
            if (muons_charge.at(i) == -13)
                all_uminus.push_back(muons.at(i));
            else 
                all_uplus.push_back(muons.at(i));
        }


        if (all_uminus.size() == 1) // I have 1 muon and 3 anti-muons
        {
            if ((all_uminus.at(0) + all_uplus.at(0)).M() > 12 && (all_uminus.at(0) + all_uplus.at(0)).M() < 120) // I have to check if I can try to combine these two into Z: I have to have 12 < mll < 120
            {
                Z1_candidates.push_back(all_uminus.at(0) + all_uplus.at(0));
                Z2_candidates.push_back(all_uplus.at(1) + all_uplus.at(2));
            }
            if ((all_uminus.at(0) + all_uplus.at(1)).M() > 12 && (all_uminus.at(0) + all_uplus.at(1)).M() < 120)
            {
                Z1_candidates.push_back(all_uminus.at(0) + all_uplus.at(1));
                Z2_candidates.push_back(all_uplus.at(0) + all_uplus.at(2));
            }
            if ((all_uminus.at(0) + all_uplus.at(2)).M() > 12 && (all_uminus.at(0) + all_uplus.at(2)).M() < 120)
            {
                Z1_candidates.push_back(all_uminus.at(0) + all_uplus.at(2));
                Z2_candidates.push_back(all_uplus.at(0) + all_uplus.at(1));
            }

            // ----------------- pairing Zs ------------------

            if (Z1_candidates.size() == 1)
            {
                Z1 = Z1_candidates.at(0);
                Z2 = Z2_candidates.at(0);
            }

            for (int i = 0; i < Z1_candidates.size(); i++)
            {
                for (int j = i; j < Z1_candidates.size(); j++)
                {
                    TLorentzVector temp_Z1, temp_Z2;
                    if (abs(Z1_candidates.at(i).M() - 91.1876) > abs(Z1_candidates.at(j).M() - 91.1876)) 
                    {
                        temp_Z1 = Z1_candidates.at(i);
                        Z1_candidates.at(i) = Z1_candidates.at(j);
                        Z1_candidates.at(j) = temp_Z1;

                        temp_Z2 = Z2_candidates.at(i);
                        Z2_candidates.at(i) = Z2_candidates.at(j);
                        Z2_candidates.at(j) = temp_Z2;
                    }
                }
            }

            Z1 = Z1_candidates.at(0);
            Z2 = Z2_candidates.at(0);
        }
        if (all_uminus.size() == 3) // I have 3 muons and 1 anti-muon
        {
            if ((all_uplus.at(0) + all_uminus.at(0)).M() > 12 && (all_uplus.at(0) + all_uminus.at(0)).M() < 120) // I have to check if I can try to combine these two into Z: I have to have 12 < mll < 120
            {
                Z1_candidates.push_back(all_uplus.at(0) + all_uminus.at(0));
                Z2_candidates.push_back(all_uminus.at(1) + all_uminus.at(2));
            }
            if ((all_uplus.at(0) + all_uminus.at(1)).M() > 12 && (all_uplus.at(0) + all_uminus.at(1)).M() < 120)
            {
                Z1_candidates.push_back(all_uplus.at(0) + all_uminus.at(1));
                Z2_candidates.push_back(all_uminus.at(0) + all_uminus.at(2));
            }
            if ((all_uplus.at(0) + all_uminus.at(2)).M() > 12 && (all_uplus.at(0) + all_uminus.at(2)).M() < 120)
            {
                Z1_candidates.push_back(all_uplus.at(0) + all_uminus.at(2));
                Z2_candidates.push_back(all_uminus.at(0) + all_uminus.at(1));
            }

            // ----------------- pairing Zs ------------------

            if (Z1_candidates.size() == 1)
            {
                Z1 = Z1_candidates.at(0);
                Z2 = Z2_candidates.at(0);
            }

            for (int i = 0; i < Z1_candidates.size(); i++)
            {
                for (int j = i; j < Z1_candidates.size(); j++)
                {
                    TLorentzVector temp_Z1, temp_Z2;
                    if (abs(Z1_candidates.at(i).M() - 91.1876) > abs(Z1_candidates.at(j).M() - 91.1876)) 
                    {
                        temp_Z1 = Z1_candidates.at(i);
                        Z1_candidates.at(i) = Z1_candidates.at(j);
                        Z1_candidates.at(j) = temp_Z1;

                        temp_Z2 = Z2_candidates.at(i);
                        Z2_candidates.at(i) = Z2_candidates.at(j);
                        Z2_candidates.at(j) = temp_Z2;
                    }
                }
            }

            Z1 = Z1_candidates.at(0);
            Z2 = Z2_candidates.at(0);
        }
    }
}