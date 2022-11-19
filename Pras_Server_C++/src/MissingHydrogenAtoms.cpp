/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
********************************************************************************************************************************/

#include "MissingHydrogenAtoms.hpp"

void AddMissingHydrogenAtons::arginine_h(Residue &res, Residue &res_p, Residue &ntermini)
    {
        V1 pos_N,pos_CA,pos_C,pos_CB,pos_CG,pos_CD,pos_NE,pos_CZ,
           pos_NH1, pos_NH2, pos_Cp, HA,H,HE,HH11,HH12,HH21,HH22,
           HB1, HB2, HG1, HG2, HD1, HD2;
        V2 HB1_HB2, HG1_HG2, HD1_HD2;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                else if(res.res_atom[i]  == "CB")  pos_CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG")  pos_CG = res.res_pos[i];
                else if (res.res_atom[i] == "CD")  pos_CD = res.res_pos[i];
                else if (res.res_atom[i] == "NE")  pos_NE = res.res_pos[i];
                else if (res.res_atom[i] == "CZ")  pos_CZ = res.res_pos[i];
                else if (res.res_atom[i] == "NH1") pos_NH1= res.res_pos[i];
                else if (res.res_atom[i] == "NH2") pos_NH2= res.res_pos[i];            
            }

        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
            }

        H         =  class5(pos_CA,pos_N,pos_Cp,1.01);
        HA        =  class3(pos_CB,pos_N,pos_CA);
        HE        =  class5(pos_CD,pos_NE,pos_CZ,1.01);
        HH11      =  class4(pos_NE,pos_CZ,pos_NH1);
        HH12      =  class4(pos_NH2,pos_CZ,pos_NH1);
        HH21      =  class4(pos_NE,pos_CZ,pos_NH2);
        HH22      =  class4(pos_NH1,pos_CZ,pos_NH2);
        HB1_HB2   =  class2(pos_CA,pos_CG,pos_CB);
        HG1_HG2   =  class2(pos_CD,pos_CB,pos_CG);
        HD1_HD2   =  class2(pos_NE,pos_CG,pos_CD);

        HB1       =  HB1_HB2[0];
        HB2       =  HB1_HB2[1];
        HG1       =  HG1_HG2[0];
        HG2       =  HG1_HG2[1];
        HD1       =  HD1_HD2[0];
        HD2       =  HD1_HD2[1];
      
        V2 arg_h = {HA, HE, HH11, HH12, HH21, HH22, HB1, HB2, HG1, HG2, HD1, HD2, H};
        RC arg_a = {"HA","HE","1HH1","2HH1","1HH2","2HH2","HB1","HB2","HG1","HG2","HD1","HD2","H"};
        for (int i = 0; i<arg_h.size(); i++)
            {
                res.res_pos.push_back(arg_h[i]);
                res.res_atom.push_back(arg_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
    }

void AddMissingHydrogenAtons::isoleucine_h(Residue &res, Residue &res_p, Residue &ntermini)
    {
        V1 pos_N,pos_CA,pos_C,pos_CB, pos_CG1, pos_CG2, pos_CD1,
           pos_Cp, HA,H,HB,HG11,HG12,HG21,HG22,HG23,HD1,HD2,HD3;
        V2 HG11_HG12;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                else if(res.res_atom[i] == "CB")   pos_CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG1") pos_CG1= res.res_pos[i];
                else if (res.res_atom[i] == "CG2") pos_CG2= res.res_pos[i];
                else if (res.res_atom[i] == "CD1") pos_CD1 = res.res_pos[i];           
            }

        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
            }

        H     =  class5(pos_CA,pos_N,pos_Cp,1.01);   
        HA             =  class3(pos_CB,pos_N,pos_CA);
        HB             =  class3(pos_CA,pos_CG1,pos_CB);
        HG11_HG12      =  class2(pos_CB,pos_CD1, pos_CG1);
        HG21           =  class1(pos_CA,pos_CB,pos_CG2, 1.09, 180.0,109.5);
        HG22           =  class1(pos_CA,pos_CB,pos_CG2, 1.09,-60.0, 109.5);
        HG23           =  class1(pos_CA,pos_CB,pos_CG2, 1.09, 60.0, 109.5);
        HD1            =  class1(pos_CB,pos_CG1,pos_CD1,1.09,-180.0,109.5);
        HD2            =  class1(pos_CB,pos_CG1,pos_CD1,1.09,-60.0, 109.5);
        HD3            =  class1(pos_CB,pos_CG1,pos_CD1,1.09, 60.0, 109.5);
        HG11           =  HG11_HG12[0];
        HG12           =  HG11_HG12[1];

        V2 ile_h = {HA, HB, HG11, HG12, HG21, HG22, HG23, HD1, HD2, HD3, H};
        RC ile_a = {"HA","HB","1HG1","2HG1","1HG2","2HG2","3HG2","HD1","HD2","HD3","H"};
        for (int i = 0; i<ile_h.size(); i++)
            {
                res.res_pos.push_back(ile_h[i]);
                res.res_atom.push_back(ile_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
    }

void AddMissingHydrogenAtons::leucine_h(Residue &res, Residue &res_p, Residue &ntermini)
    {
        V1 pos_N,pos_CA,pos_C,pos_CB, pos_CG, pos_CD1, pos_CD2,
           pos_Cp, HA,H,HB1,HB2,HG,HD11,HD12,HD13,HD21,HD22,HD23;

        V2 HB1_HB2;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                else if(res.res_atom[i] == "CB")   pos_CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG")  pos_CG= res.res_pos[i];
                else if (res.res_atom[i] == "CD1") pos_CD1= res.res_pos[i];
                else if (res.res_atom[i] == "CD2") pos_CD2 = res.res_pos[i];           
            }

        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
            }

        H              =  class5(pos_CA,pos_N,pos_Cp,1.01);
        HA             =  class3(pos_CB,pos_N,pos_CA);
        HB1_HB2        =  class2(pos_CA,pos_CG,pos_CB);
        HG             =  class3(pos_CD1,pos_CB,pos_CG);
        HD11           =  class1(pos_CB,pos_CG,pos_CD1,1.09, 180.0,109.5);
        HD12           =  class1(pos_CB,pos_CG,pos_CD1,1.09, -60.0,109.5);
        HD13           =  class1(pos_CB,pos_CG,pos_CD1,1.09,  60.0,109.5);
        HD21           =  class1(pos_CB,pos_CG,pos_CD2,1.09,-180.0,109.5);
        HD22           =  class1(pos_CB,pos_CG,pos_CD2,1.09, -60.0,109.5);
        HD23           =  class1(pos_CB,pos_CG,pos_CD2, 1.09, 60.0,109.5);
        HB1            =  HB1_HB2[0];
        HB2            =  HB1_HB2[1];

        V2 leu_h = {HA, HB1, HB2, HG, HD11, HD12, HD13, HD21, HD22, HD23, H};
        RC leu_a = {"HA","HB1","HB2","HG","1HD1","2HD1","3HD1","1HD2","2HD2","3HD2","H"};
        for (int i = 0; i<leu_h.size(); i++)
            {
                res.res_pos.push_back(leu_h[i]);
                res.res_atom.push_back(leu_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
    }

void AddMissingHydrogenAtons::tryptophan_h(Residue &res, Residue &res_p, Residue &ntermini)
    {
        V1 pos_N,pos_CA,pos_C,pos_CB, pos_CG, pos_CD1, pos_CD2,
           pos_CE2, pos_CE3, pos_NE1, pos_CZ2, pos_CZ3, pos_CH2,
           pos_Cp, HA,H, HB1, HB2, HD1, HE1, HE3, HZ2, HZ3, HH2;

        V2 HB1_HB2;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                else if(res.res_atom[i]  == "CB")  pos_CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG")  pos_CG = res.res_pos[i];
                else if (res.res_atom[i] == "CD1") pos_CD1= res.res_pos[i];
                else if (res.res_atom[i] == "CD2") pos_CD2= res.res_pos[i];  
                else if(res.res_atom[i]  == "CE2") pos_CE2= res.res_pos[i];
                else if (res.res_atom[i] == "CE3") pos_CE3= res.res_pos[i];
                else if (res.res_atom[i] == "NE1") pos_NE1= res.res_pos[i];
                else if (res.res_atom[i] == "CZ2") pos_CZ2= res.res_pos[i];
                else if (res.res_atom[i] == "CZ3") pos_CZ3= res.res_pos[i];
                else if (res.res_atom[i] == "CH2") pos_CH2= res.res_pos[i];           
            }

        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
            }

        H              =  class5(pos_CA,pos_N,pos_Cp,1.01);
        HA             =  class3(pos_CB,pos_N,pos_CA);
        HB1_HB2        =  class2(pos_CA,pos_CG,pos_CB);
        HD1            =  class5(pos_NE1,pos_CD1,pos_CG,1.08);
        HE1            =  class5(pos_CE2,pos_NE1,pos_CD1,1.01);
        HE3            =  class5(pos_CZ3,pos_CE3,pos_CD2,1.08);
        HZ2            =  class5(pos_CE2,pos_CZ2,pos_CH2,1.08);
        HZ3            =  class5(pos_CH2,pos_CZ3,pos_CE3,1.08);
        HH2            =  class5(pos_CZ2,pos_CH2,pos_CZ3,1.08);
        HB1            =  HB1_HB2[0];
        HB2            =  HB1_HB2[1];

        V2 trp_h = {HA, HB1, HB2, HD1, HE1, HE3, HZ2, HZ3, HH2, H};
        RC trp_a = {"HA","HB1","HB2","HD1","HE1","HE3","HZ2","HZ3","HH2","H"};
        for (int i = 0; i<trp_h.size(); i++)
            {
                res.res_pos.push_back(trp_h[i]);
                res.res_atom.push_back(trp_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
    }

void AddMissingHydrogenAtons::phenylalanine_h(Residue &res, Residue &res_p, Residue &ntermini)
    {

        V1 pos_N,pos_CA,pos_C,pos_CB, pos_CG, pos_CD1, pos_CD2,
           pos_CE1, pos_CE2, pos_CZ,pos_Cp, HA,H, HB1, HB2, HD1, 
           HD2, HE1, HE2, HZ;

        V2 HB1_HB2;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                else if(res.res_atom[i]  == "CB")  pos_CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG")  pos_CG = res.res_pos[i];
                else if (res.res_atom[i] == "CD1") pos_CD1= res.res_pos[i];
                else if (res.res_atom[i] == "CD2") pos_CD2= res.res_pos[i];  
                else if(res.res_atom[i]  == "CE1") pos_CE1= res.res_pos[i];
                else if (res.res_atom[i] == "CE2") pos_CE2= res.res_pos[i];
                else if (res.res_atom[i] == "CZ")  pos_CZ = res.res_pos[i];          
            }

        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
            }

        H              =  class5(pos_CA,pos_N,pos_Cp,1.01);
        HA             =  class3(pos_CB,pos_N,pos_CA);
        HB1_HB2        =  class2(pos_CA,pos_CG,pos_CB);
        HD1            =  class5(pos_CE1,pos_CD1,pos_CG,1.08);
        HD2            =  class5(pos_CE2,pos_CD2,pos_CG,1.08);
        HE1            =  class5(pos_CZ,pos_CE1,pos_CD1,1.08);
        HE2            =  class5(pos_CD2,pos_CE2,pos_CZ,1.08);
        HZ             =  class5(pos_CE2,pos_CZ,pos_CE1,1.08);
        HB1            =  HB1_HB2[0];
        HB2            =  HB1_HB2[1];

        V2 phe_h = {HA, HB1, HB2, HD1, HD2, HE1, HE2, HZ, H};
        RC phe_a = {"HA","HB1","HB2","HD1","HD2","HE1", "HE2", "HZ", "H"};
        for (int i = 0; i<phe_h.size(); i++)
            {
                res.res_pos.push_back(phe_h[i]);
                res.res_atom.push_back(phe_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
    }

void AddMissingHydrogenAtons::asparagine_h(Residue &res, Residue &res_p, Residue &ntermini)
    {
        V1 pos_N,pos_CA,pos_C,pos_CB, pos_CG, pos_ND2, 
        pos_OD1, pos_Cp, HA, H, HB1, HB2, HD12, HD22;

        V2 HB1_HB2;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                else if(res.res_atom[i] == "CB")   pos_CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG")  pos_CG= res.res_pos[i];
                else if (res.res_atom[i] == "ND2") pos_ND2= res.res_pos[i];
                else if (res.res_atom[i] == "OD1") pos_OD1 = res.res_pos[i];           
            }

        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
            }

        H         =  class5(pos_CA,pos_N,pos_Cp,1.01);
        HA        =  class3(pos_CB,pos_N,pos_CA);
        HB1_HB2   =  class2(pos_CA,pos_CG,pos_CB);
        HD12      =  class4(pos_CB,pos_CG,pos_ND2);
        HD22      =  class4(pos_OD1,pos_CG,pos_ND2);
        HB1       =  HB1_HB2[0];
        HB2       =  HB1_HB2[1];

        V2 asn_h = {HA, HB1, HB2, HD12, HD22, H};
        RC asn_a = {"HA","HB1","HB2","1HD2","2HD2","H"};
        for (int i = 0; i<asn_h.size(); i++)
            {
                res.res_pos.push_back(asn_h[i]);
                res.res_atom.push_back(asn_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
    }

void AddMissingHydrogenAtons::alanine_h(Residue &res, Residue &res_p, Residue &ntermini)
    {
        V1 pos_N, pos_CA, pos_C, pos_CB, pos_Cp, HA, H, HB1, HB2, HB3;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                else if(res.res_atom[i] == "CB")   pos_CB = res.res_pos[i];
            }

        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
            }

            H              =  class5(pos_CA,pos_N,pos_Cp,1.01);
            HA             =  class3(pos_CB,pos_N,pos_CA);
            HB1            =  class1(pos_C,pos_CA,pos_CB,1.09,55.8,109.5);
            HB2            =  class1(pos_C,pos_CA,pos_CB,1.09,175.8,109.5);
            HB3            =  class1(pos_C,pos_CA,pos_CB,1.09,-64.2,109.5);

        V2 ala_h = {HA, HB1, HB2, HB3, H};
        RC ala_a = {"HA","HB1","HB2","HB3","H"};
        for (int i = 0; i<ala_h.size(); i++)
            {
                res.res_pos.push_back(ala_h[i]);
                res.res_atom.push_back(ala_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
    }

void AddMissingHydrogenAtons::glycine_h(Residue &res, Residue &res_p, Residue &ntermini)
    {
        V1 pos_N, pos_CA, pos_C, pos_Cp, HA1, HA2, H;
        V2 HA1_HA2;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
            }

        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
            }

            H         =  class5(pos_CA,pos_N,pos_Cp,1.01);
            HA1_HA2   =  class2(pos_C,pos_N,pos_CA);
            HA1       =  HA1_HA2[0];
            HA2       =  HA1_HA2[1];

        V2 gly_h = {HA1, HA2, H};
        RC gly_a = {"HA1","HA2", "H"};
        for (int i = 0; i<gly_h.size(); i++)
            {
                res.res_pos.push_back(gly_h[i]);
                res.res_atom.push_back(gly_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_C, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_C, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_C, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
    }

void AddMissingHydrogenAtons::valine_h(Residue &res, Residue &res_p, Residue &ntermini)
{
    V1 pos_N,pos_CA,pos_C,pos_CB, pos_CG1, pos_CG2,pos_Cp, 
       HA, H, HB,  HG11,  HG12,  HG13,  HG21,  HG22, HG23;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                else if(res.res_atom[i] == "CB")   pos_CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG1") pos_CG1= res.res_pos[i];
                else if (res.res_atom[i] == "CG2") pos_CG2= res.res_pos[i];           
            }

        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
            }

        H              =  class5(pos_CA,pos_N,pos_Cp,1.01);
        HA             =  class3(pos_CB,pos_N,pos_CA);
        HB             =  class3(pos_CA,pos_CG2,pos_CB);
        HG11           =  class1(pos_CA,pos_CB,pos_CG1,1.09,180.0,109.5);
        HG12           =  class1(pos_CA,pos_CB,pos_CG1,1.09,-60.0,109.5);
        HG13           =  class1(pos_CA,pos_CB,pos_CG1,1.09,60.0,109.5);
        HG21           =  class1(pos_CG1,pos_CB,pos_CG2,1.09,-58.1,109.5);
        HG22           =  class1(pos_CG1,pos_CB,pos_CG2,1.09,61.9,109.5);
        HG23           =  class1(pos_CG1,pos_CB,pos_CG2,1.09,-178.1,109.5);

        V2 val_h = {HA, HB, HG11, HG12, HG13, HG21, HG22, HG23, H};
        RC val_a = {"HA","HB","1HG1","2HG1","3HG1","1HG2","2HG2","3HG2","H"};
        for (int i = 0; i<val_h.size(); i++)
            {
                res.res_pos.push_back(val_h[i]);
                res.res_atom.push_back(val_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
    }

void AddMissingHydrogenAtons::lysine_h(Residue &res, Residue &res_p, Residue &ntermini)
    {
        V1 pos_N,pos_CA,pos_C,pos_CB, pos_CG, pos_CD,pos_CE, pos_NZ,pos_Cp, 
           HA, H, HB1, HB2, HG1,  HG2,  HD1, HD2, HE1,  HE2, HZ1, HZ2, HZ3;

        V2 HB1_HB2, HG1_HG2, HD1_HD2, HE1_HE2;

            for (int i = 0; i<res.res_atom.size(); i++)
                {
                    if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                    else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                    else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                    else if(res.res_atom[i] == "CB")   pos_CB = res.res_pos[i];
                    else if (res.res_atom[i] == "CG")  pos_CG = res.res_pos[i];
                    else if (res.res_atom[i] == "CD")  pos_CD = res.res_pos[i];
                    else if (res.res_atom[i] == "CE")  pos_CE = res.res_pos[i];
                    else if (res.res_atom[i] == "NZ")  pos_NZ = res.res_pos[i];             
                }

            for (int i = 0; i<res_p.res_atom.size(); i++)
                {
                    if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
                }

            H              =  class5(pos_CA,pos_N,pos_Cp,1.01);
            HA             =  class3(pos_CB,pos_N,pos_CA);
            HB1_HB2        =  class2(pos_CA,pos_CG,pos_CB);
            HG1_HG2        =  class2(pos_CB,pos_CD,pos_CG);
            HD1_HD2        =  class2(pos_CE,pos_CG,pos_CD);
            HE1_HE2        =  class2(pos_NZ,pos_CD,pos_CE);
            HZ1            =  class1(pos_CD,pos_CE,pos_NZ,1.01,180.0,109.3);
            HZ2            =  class1(pos_CD,pos_CE,pos_NZ,1.01,-60.0,109.5);
            HZ3            =  class1(pos_CD,pos_CE,pos_NZ,1.01,60.0,109.5);
            HB1            = HB1_HB2[0];
            HB2            = HB1_HB2[1];
            HG1            = HG1_HG2[0];
            HG2            = HG1_HG2[1];
            HD1            = HD1_HD2[0];
            HD2            = HD1_HD2[1];
            HE1            = HE1_HE2[0];
            HE2            = HE1_HE2[1];

        V2 lys_h = {HA, HB1, HB2, HG1, HG2, HD1, HD2, HE1, HE2, HZ1, HZ2, HZ3, H};
        RC lys_a = {"HA","HB1","HB2","HG1","HG2","HD1","HD2","HE1","HE2","HZ1","HZ2","HZ3","H"};
        for (int i = 0; i<lys_h.size(); i++)
            {
                res.res_pos.push_back(lys_h[i]);
                res.res_atom.push_back(lys_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
    }

void AddMissingHydrogenAtons::aspartate_h(Residue &res, Residue &res_p, Residue &ntermini)
    {
        V1 pos_N,pos_CA,pos_C,pos_CB,pos_CG,pos_OD1, 
           pos_OD2,   pos_Cp,   HA,   H,  HB1,  HB2;

        V2 HB1_HB2;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                else if(res.res_atom[i] == "CB")   pos_CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG")  pos_CG= res.res_pos[i];          
            }

        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
            }

            H         = class5(pos_CA,pos_N,pos_Cp,1.01);
            HA        = class3(pos_CB,pos_N,pos_CA);
            HB1_HB2   = class2(pos_CA,pos_CG,pos_CB);
            HB1       = HB1_HB2[0];
            HB2       = HB1_HB2[1];

        V2 asp_h = {HA, HB1, HB2, H};
        RC asp_a = {"HA","HB1","HB2","H"};
        for (int i = 0; i<asp_h.size(); i++)
            {
                res.res_pos.push_back(asp_h[i]);
                res.res_atom.push_back(asp_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
    }

void AddMissingHydrogenAtons::glutamine_h(Residue &res, Residue &res_p, Residue &ntermini)
    {
        V1 pos_N,pos_CA,pos_C,pos_CB,pos_CG,pos_CD, pos_NE2, 
        pos_OE1, pos_Cp, HA, H, HB1,HB2,HG1,HG2, HE21, HE22;

        V2 HB1_HB2, HG1_HG2;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                else if(res.res_atom[i] == "CB")   pos_CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG")  pos_CG = res.res_pos[i];
                else if(res.res_atom[i] == "CD")   pos_CD = res.res_pos[i];
                else if (res.res_atom[i] == "NE2") pos_NE2= res.res_pos[i]; 
                else if (res.res_atom[i] == "OE1") pos_OE1= res.res_pos[i];         
            }

        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
            }

            H              =  class5(pos_CA,pos_N,pos_Cp,1.01);
            HA             =  class3(pos_CB,pos_N,pos_CA);
            HB1_HB2        =  class2(pos_CA,pos_CG,pos_CB);
            HG1_HG2        =  class2(pos_CB,pos_CD,pos_CG);
            HE21           =  class4(pos_CG,pos_CD,pos_NE2);
            HE22           =  class4(pos_OE1,pos_CD,pos_NE2);
            HB1            = HB1_HB2[0];
            HB2            = HB1_HB2[1];
            HG1            = HG1_HG2[0];
            HG2            = HG1_HG2[1];

        V2 gln_h = {HA, HB1, HB2, HG1, HG2, HE21, HE22, H};
        RC gln_a = {"HA","HB1","HB2","HG1","HG2","1HE2","2HE2","H"};
        for (int i = 0; i<gln_h.size(); i++)
            {
                res.res_pos.push_back(gln_h[i]);
                res.res_atom.push_back(gln_a[i]);
            }
        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
    }

void AddMissingHydrogenAtons::glutamate_h(Residue &res, Residue &res_p, Residue &ntermini)
    {
        V1 pos_N,pos_CA,pos_C,pos_CB,pos_CG,pos_CD, 
           pos_Cp,  HA,  H,  HB1,  HB2,  HG1,  HG2;

        V2 HB1_HB2, HG1_HG2;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                else if(res.res_atom[i] == "CB")   pos_CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG")  pos_CG = res.res_pos[i];
                else if(res.res_atom[i] == "CD")   pos_CD = res.res_pos[i];         
            }

        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
            }

            H              =  class5(pos_CA,pos_N,pos_Cp,1.01);
            HA             =  class3(pos_CB,pos_N,pos_CA);
            HB1_HB2        =  class2(pos_CA,pos_CG,pos_CB);
            HG1_HG2        =  class2(pos_CB,pos_CD,pos_CG);
            HB1            = HB1_HB2[0];
            HB2            = HB1_HB2[1];
            HG1            = HG1_HG2[0];
            HG2            = HG1_HG2[1];

        V2 glu_h = {HA, HB1, HB2, HG1, HG2, H};
        RC glu_a = {"HA","HB1","HB2","HG1","HG2","H"};
        for (int i = 0; i<glu_h.size(); i++)
            {
                res.res_pos.push_back(glu_h[i]);
                res.res_atom.push_back(glu_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
    }

void AddMissingHydrogenAtons::histidine_h(Residue &res, Residue &res_p, Residue &ntermini)
    {
        V1 pos_N,pos_CA,pos_C,pos_CB,pos_CG,pos_CD2,pos_ND1, 
           pos_CE1,pos_NE2,pos_Cp,HA,H,HB1,HB2,HD2,HE1, HE2;

        V2 HB1_HB2;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                else if(res.res_atom[i] == "CB")   pos_CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG")  pos_CG = res.res_pos[i];
                else if(res.res_atom[i] == "CD2")  pos_CD2= res.res_pos[i];  
                else if(res.res_atom[i] == "ND1")  pos_ND1= res.res_pos[i];
                else if (res.res_atom[i] == "CD2") pos_CD2= res.res_pos[i];
                else if(res.res_atom[i] == "CE1")  pos_CE1= res.res_pos[i]; 
                else if(res.res_atom[i] == "NE2")  pos_NE2= res.res_pos[i];      
            }

        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
            }

            H              =  class5(pos_CA,pos_N,pos_Cp,1.01);
            HA             =  class3(pos_CB,pos_N,pos_CA);
            HB1_HB2        =  class2(pos_CA,pos_CG,pos_CB);
            HD2            =  class5(pos_NE2,pos_CD2,pos_CG,1.08);
            HE1            =  class5(pos_NE2,pos_CE1,pos_ND1,1.08);
            HE2            =  class5(pos_CD2,pos_NE2,pos_CE1,1.01);
            HB1            = HB1_HB2[0];
            HB2            = HB1_HB2[1];

        V2 his_h = {HA, HB1, HB2, HD2, HE1, HE2, H};
        RC his_a = {"HA","HB1","HB2","HD2","HE1","HE2","H"};
        for (int i = 0; i<his_h.size(); i++)
            {
                res.res_pos.push_back(his_h[i]);
                res.res_atom.push_back(his_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
     } 

void AddMissingHydrogenAtons::proline_h(Residue &res, Residue &res_p, Residue &ntermini)
    {   
        V1 pos_N,pos_CA, pos_C, pos_CB, pos_CG, pos_CD, 
           HA, HB1, HB2, HD1, HD2, HG1, HG2;

        V2 HB1_HB2, HG1_HG2, HD1_HD2;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                else if(res.res_atom[i]  == "CB")  pos_CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG")  pos_CG = res.res_pos[i];
                else if(res.res_atom[i]  == "CD")  pos_CD = res.res_pos[i];        
            }

            HA             =  class3(pos_CB,pos_N,pos_CA);
            HB1_HB2        =  class2(pos_CA,pos_CG,pos_CB);
            HG1_HG2        =  class2(pos_CB,pos_CD,pos_CG);
            HD1_HD2        =  class2(pos_CG,pos_N,pos_CD);
            HB1            = HB1_HB2[0];
            HB2            = HB1_HB2[1];
            HG1            = HG1_HG2[0];
            HG2            = HG1_HG2[1];
            HD1            = HD1_HD2[0];
            HD2            = HD1_HD2[1];

        V2 pro_h = {HA, HB1, HB2, HG1, HG2, HD1, HD2};
        RC pro_a = {"HA","HB1","HB2", "HG1","HG2","HD1","HD2"};
        for (int i = 0; i<pro_h.size(); i++)
            {
                res.res_pos.push_back(pro_h[i]);
                res.res_atom.push_back(pro_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V2  H1_H2 = class2(pos_CA, pos_CD, pos_N, 1.01);
              V1  H1 = H1_H2[0];
              V1  H2 = H1_H2[1];
              res.res_pos.push_back(H1);   res.res_atom.push_back("H1");
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
            }
    }

void AddMissingHydrogenAtons::metheonine_h(Residue &res, Residue &res_p, Residue &ntermini)
    {
        V1 pos_N, pos_CA, pos_C, pos_CB, pos_CG, pos_SD, pos_CE, 
           pos_Cp, HA, H, HB1,  HB2,  HG1,  HG2,  HE1, HE2, HE3;

        V2 HB1_HB2, HG1_HG2;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                else if(res.res_atom[i]  == "CB")  pos_CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG")  pos_CG = res.res_pos[i];
                else if(res.res_atom[i]  == "SD")  pos_SD = res.res_pos[i]; 
                else if(res.res_atom[i]  == "CE")  pos_CE = res.res_pos[i];        
            }

        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
            }

            H              =  class5(pos_CA,pos_N,pos_Cp,1.01);
            HA             =  class3(pos_CB,pos_N,pos_CA);
            HB1_HB2        =  class2(pos_CA,pos_CG,pos_CB);
            HG1_HG2        =  class2(pos_CB,pos_SD,pos_CG);
            HE1            =  class1(pos_CG,pos_SD,pos_CE,1.09,180.0,109.4);
            HE2            =  class1(pos_CG,pos_SD,pos_CE,1.09,-60.0,109.5);
            HE3            =  class1(pos_CG,pos_SD,pos_CE,1.09,60.0,109.5);
            HB1            = HB1_HB2[0];
            HB2            = HB1_HB2[1];
            HG1            = HG1_HG2[0];
            HG2            = HG1_HG2[1];

        V2 met_h = {HA, HB1, HB2, HG1, HG2, HE1, HE2, HE3, H};
        RC met_a = {"HA","HB1","HB2","HG1","HG2","HE1","HE2","HE3","H"};
        for (int i = 0; i<met_h.size(); i++)
            {
                res.res_pos.push_back(met_h[i]);
                res.res_atom.push_back(met_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
    }

void AddMissingHydrogenAtons::threonine_h(Residue &res, Residue &res_p, Residue &ntermini, vector<Residue> &chains)
    {
        V1 pos_N,pos_CA,pos_C,pos_CB,pos_CG2,pos_OG1, 
           pos_Cp, HA, HB, H, HG1, HG12, HG22, HG23;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                else if(res.res_atom[i]  == "CB")  pos_CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG2") pos_CG2= res.res_pos[i];
                else if(res.res_atom[i]  == "OG1") pos_OG1= res.res_pos[i];        
            }

        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
            }

        H              =  class5(pos_CA,pos_N,pos_Cp,1.01);
        HA             =  class3(pos_CB,pos_N,pos_CA);
        HG1            =  class6_Thr(pos_OG1, pos_CB, pos_CG2);
        HB             =  class3(pos_CA,pos_OG1,pos_CB);
        HG12           =  class1(pos_OG1,pos_CB,pos_CG2,1.09,60.5,109.4);
        HG22           =  class1(pos_OG1,pos_CB,pos_CG2,1.09,-179.5,109.5);
        HG23           =  class1(pos_OG1,pos_CB,pos_CG2,1.09,-59.5,109.5);
        float di       =  calcTorsionAngle(pos_CA, pos_CB, pos_OG1, HG1);
        HG1            = optmizeH(pos_CA, pos_CB, pos_OG1, HG1, 0.96, di, 0.41, 0., 0., res, chains);

        V2 thr_h = {HG1, HA, HB, HG12, HG22, HG23, H};
        RC thr_a = {"HG1","HA","HB","1HG2","2HG2","3HG2","H"};
        for (int i = 0; i<thr_h.size(); i++)
            {
                res.res_pos.push_back(thr_h[i]);
                res.res_atom.push_back(thr_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
    }

void AddMissingHydrogenAtons::serine_h(Residue &res, Residue &res_p, Residue &ntermini, vector<Residue> &chains)
{
    V1 pos_N,pos_CA,pos_C,pos_CB,pos_OG, 
       pos_Cp, HA, HB1, HB2, HG, H;

    V2 HB1_HB2;

    for (int i = 0; i<res.res_atom.size(); i++)
        {
            if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
            else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
            else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
            else if(res.res_atom[i]  == "CB")  pos_CB = res.res_pos[i];
            else if(res.res_atom[i]  == "OG")  pos_OG = res.res_pos[i];        
        }

    for (int i = 0; i<res_p.res_atom.size(); i++)
        {
            if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
        }

        H              =  class5(pos_CA,pos_N,pos_Cp,1.01);
        HA             =  class3(pos_CB,pos_N,pos_CA);
        HB1_HB2        =  class2(pos_CA,pos_OG,pos_CB);
        HG             =  class6_Ser(pos_OG, pos_CB, HB1_HB2[0]);
        HB1            =  HB1_HB2[0];
        HB2            =  HB1_HB2[1];
        float di       =  calcTorsionAngle(pos_CA,pos_CB,pos_OG,HG);
        HG             =  optmizeH(pos_CA, pos_CB, pos_OG, HG, 0.96, di, 0.41, 0., 0., res, chains);

        V2 ser_h = {HA, HB1, HB2, HG, H};
        RC ser_a = {"HA","HB1","HB2","HG","H"};
        for (int i = 0; i<ser_h.size(); i++)
            {
                res.res_pos.push_back(ser_h[i]);
                res.res_atom.push_back(ser_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
}

void AddMissingHydrogenAtons::tyrosine_h(Residue &res, Residue &res_p, Residue &ntermini, vector<Residue> &chains)
{
    V1 pos_N,pos_CA,pos_C,pos_CB,pos_CG,pos_CD1,pos_CD2,pos_CE1,pos_CE2,
    pos_CZ, pos_OH, pos_Cp,  HA, HB1, HB2, HD1, HD2, HE1, HE2, HH, H;

    V2 HB1_HB2;

    for (int i = 0; i<res.res_atom.size(); i++)
        {
            if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
            else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
            else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
            else if(res.res_atom[i]  == "CB")  pos_CB = res.res_pos[i];
            else if (res.res_atom[i] == "CG")  pos_CG = res.res_pos[i];
            else if(res.res_atom[i]  == "CD1") pos_CD1= res.res_pos[i];  
            else if(res.res_atom[i]  == "CD2") pos_CD2= res.res_pos[i]; 
            else if(res.res_atom[i]  == "CE1") pos_CE1= res.res_pos[i];  
            else if(res.res_atom[i]  == "CE2") pos_CE2= res.res_pos[i];  
            else if(res.res_atom[i]  == "CZ")  pos_CZ = res.res_pos[i];  
            else if(res.res_atom[i]  == "OH")  pos_OH = res.res_pos[i];      
        }

    for (int i = 0; i<res_p.res_atom.size(); i++)
        {
            if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
        }

        H              =  class5(pos_CA,pos_N,pos_Cp,1.01);
        HA             =  class3(pos_CB,pos_N,pos_CA);
        HB1_HB2        =  class2(pos_CA,pos_CG,pos_CB);
        HD1            =  class5(pos_CG,pos_CD1,pos_CE1,1.08);
        HD2            =  class5(pos_CE2,pos_CD2,pos_CG,1.08);
        HE1            =  class5(pos_CZ,pos_CE1,pos_CD1,1.08);
        HE2            =  class5(pos_CZ,pos_CE2,pos_CD2,1.08);
        HH             =  class6_Tyr(pos_OH, pos_CZ, pos_CE2);
        HB1            =  HB1_HB2[0];
        HB2            =  HB1_HB2[1];
        float di       =  calcTorsionAngle(pos_CE2, pos_CZ, pos_OH, HH);
        HH             =  optmizeH(pos_CE2, pos_CZ, pos_OH, HH, 0.96, di, 0.41, 0., 0., res, chains);

        V2 tyr_h = {HA, HB1, HB2, HD1, HD2, HE1 ,HE2, HH, H};
        RC tyr_a = {"HA","HB1","HB2","HD1","HD2","HE1","HE2","HH","H"};
        for (int i = 0; i<tyr_h.size(); i++)
            {
                res.res_pos.push_back(tyr_h[i]);
                res.res_atom.push_back(tyr_a[i]);
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
}

void AddMissingHydrogenAtons::cystine_h(Residue &res, Residue &res_p, vector<Residue> &chains, Residue &ntermini)
    {
        V1 pos_N,pos_CA,pos_C,pos_CB,pos_SG,
           pos_Cp, HA, HB1, HB2, HG, H;

        V2 HB1_HB2;

        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")        pos_N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA")  pos_CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")   pos_C  = res.res_pos[i];
                else if(res.res_atom[i]  == "CB")  pos_CB = res.res_pos[i];
                else if (res.res_atom[i] == "SG")  pos_SG = res.res_pos[i];       
            }

        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "C")  pos_Cp   = res_p.res_pos[i];
            }

            H              =  class5(pos_CA,pos_N,pos_Cp,1.01);
            HA             =  class3(pos_CB,pos_N,pos_CA);
            HB1_HB2        =  class2(pos_CA,pos_SG,pos_CB);
            HB1            = HB1_HB2[0];
            HB2            = HB1_HB2[1];

        string ds = "";
        for (int j = 0; j<chains.size(); j++)
            {
             if (chains[j].res_type == "CYS" && chains[j].c_id == res.c_id) 
             {
                if (chains[j].res_num != res.res_num )
                 {

                    float dist = _distance(pos_SG, chains[j].res_pos[5]);
                     if ( dist <= 3.)
                        { 
                            ds = "disulfide";
                            break;  
                        } 
                 }                
             }
            }
        if (ds != "disulfide")
            {
                HG       = class6_Cys(pos_SG, pos_CB, pos_CA);
                float di = calcTorsionAngle(pos_CA,pos_CB,pos_SG,HG);
                HG       = optmizeH(pos_CA, pos_CB, pos_SG, HG, 1.34, di, 0.19, 0.11, 0.07, res, chains);
                V2 cys_h = {HA, HB1, HB2, HG, H};
                RC cys_a = {"HA","HB1","HB2","HG","H"};
                for (int i = 0; i<cys_h.size(); i++)
                    {
                        res.res_pos.push_back(cys_h[i]);
                        res.res_atom.push_back(cys_a[i]);
                    }              
            }
        else
            {
                V2 cys_h = {HA, HB1, HB2, H};
                RC cys_a = {"HA","HB1","HB2","H"};
                for (int i = 0; i<cys_h.size(); i++)
                    {
                        res.res_pos.push_back(cys_h[i]);
                        res.res_atom.push_back(cys_a[i]);
                    }
            }

        if ((res.res_num==ntermini.res_num&&res.c_id==ntermini.c_id&&res.res_insert==ntermini.res_insert)||res.c_id!=res_p.c_id)
            {
              V1  H1 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.5, 179.6);
              V1  H2 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, -60.4);
              V1  H3 = calcCoordinate(pos_CB, pos_CA, pos_N, 1.01, 109.6, 60.2);
              res.res_pos.back() = H1;   res.res_atom.back() = "H1";
              res.res_pos.push_back(H2); res.res_atom.push_back("H2");
              res.res_pos.push_back(H3); res.res_atom.push_back("H3");
            }
    }

void AddMissingHydrogenAtons::his_note(string resNo, string chainID)
    {
        resNo.erase(remove_if(resNo.begin(),resNo.end(),::isspace),resNo.end());
        ofstream file("log.txt", ios_base::app);
        file<<"HIS No."<< resNo <<" in chain " 
             <<chainID<< " is protonated (+1 charge)";
        file<<"(ref @ www.protein-science.com).\n\n";   
    }

void AddMissingHydrogenAtons::get_hd1(vector<Residue> &HIS, vector<Residue> &chains)
    {   V1 pos_CG, pos_ND1, pos_CE1;
        int prot = floor(HIS.size()/5);

        for (int i = 0; i < prot; i++)
        {
            for (int j = 0; j <HIS[i].res_atom.size(); j++)
            {
                if     (HIS[i].res_atom[j] == "CG")   pos_CG  = HIS[i].res_pos[j]; 
                else if(HIS[i].res_atom[j] == "ND1")  pos_ND1 = HIS[i].res_pos[j];
                else if(HIS[i].res_atom[j] == "CE1")  pos_CE1 = HIS[i].res_pos[j];               
            }

            V1 HD1 =  class5(pos_CE1, pos_ND1, pos_CG, 1.01);
            for (int k = 0; k < chains.size(); k++)
            {
                if (chains[k].res_num == HIS[i].res_num && chains[k].c_id == HIS[i].c_id)
                {
                    chains[k].res_atom.push_back("HD1");
                    chains[k].res_pos.push_back(HD1);
                    his_note(HIS[i].res_num, HIS[i].c_id);                    
                }
            }
        }

        HIS.clear();
    }

void AddMissingHydrogenAtons::prot_his(vector<Residue> &chains)
    {
        vector<Residue> HIS = {};
        for(int i = 0; i < chains.size(); i++)
        {
            if (chains[i].res_type == "HIS" && HIS.size() == 0)
            {
                HIS.push_back(chains[i]);
            }
            else if (chains[i].res_type == "HIS" && chains[i].c_id == HIS.back().c_id)
            {
                HIS.push_back(chains[i]);
            }
            else if (chains[i].res_type == "HIS" && chains[i].c_id != HIS.back().c_id)
            {
                get_hd1(HIS, chains);
                HIS.push_back(chains[i]); 
            }
        }

        if (HIS.size()>= 5)
            {
                get_hd1(HIS, chains);
            }
    }
