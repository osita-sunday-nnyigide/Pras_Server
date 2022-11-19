/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
********************************************************************************************************************************/

#include "MissingHeavyAtoms.hpp"

V1 FixMissingAtoms::c_termini(V1 N, V1 CA, V1 C, V1 O)
    {
        float dih = calcTorsionAngle(N, CA, C, O);
        return calcCoordinate(N, CA, C, 1.25,122.5,180.0+dih);
    }

float FixMissingAtoms::get_chi(string a1,string a2,string a3,string a4,Residue &res,vector<Residue> &faspr)
    {
        V1 pos1, pos2, pos3, pos4;
        for (int i=0; i<faspr.size(); i++)
        {
            if (faspr[i].res_num == res.res_num && faspr[i].c_id == res.c_id)
            {
                for (int j=0; j<faspr[i].res_atom.size(); j++)
                {
                    if (faspr[i].res_atom[j] == a1)       pos1 = faspr[i].res_pos[j]; 
                    else if (faspr[i].res_atom[j] == a2)  pos2 = faspr[i].res_pos[j];
                    else if (faspr[i].res_atom[j] == a3)  pos3 = faspr[i].res_pos[j];
                    else if (faspr[i].res_atom[j] == a4)  pos4 = faspr[i].res_pos[j];
                }
            try
                {
                    if (pos1.empty() || pos2.empty() || pos3.empty() || pos4.empty()) 
                        {
                             throw 1;
                        }
                    else return calcTorsionAngle(pos1, pos2, pos3, pos4);
                }
            catch(int i)
                {
                    cout<<"PRAS failed to get chi for residue No."<<res.res_num<<"in chain "
                    <<res.c_id<<".\nThis happens if there is a point mutation and" 
                    "\nyou chosed low occupancy whereas FASPR used high occupancy.\n"
                    "The solution is to delete high occupancy residue in the PDB submitted\nto" 
                    " FASPR." <<" PRAS has terminated abnormally*********************"<<endl;
                    exit(1);
                }
            }

        }

    }

BB FixMissingAtoms::addBackbone(V1 N, V1 CA, V1 C, Residue nxtres, float psi)
    {   
        V1 O = {};
        if (psi >= 100.0)
            {
                V1 di = {-30.0, 140.0, -140.0, 40.0};
                for (auto i : di)
                    {
                        O = calcCoordinate(N, CA, C, 1.23, 120.5, i);
                        if (_distance(O, nxtres.res_pos[0]) > 2.0) break;                  
                    }
            }
        else if (psi > 0.0 && psi < 100.0)
            {
                V1 di = {-140.0, 140.0, -30.0, 40.0}; 
                for (auto i : di) 
                    {
                        O = calcCoordinate(N, CA, C, 1.23, 120.5, i);
                        if (_distance(O, nxtres.res_pos[0]) > 2.0) break;                   
                    } 
            }
        else
            {
                V1 di = {140.0, -140.0, -30.0, 40.0};
                for (auto i : di)
                    {
                        O = calcCoordinate(N, CA, C, 1.23, 120.5, i);
                        if (_distance(O, nxtres.res_pos[0]) > 2.0) break;               
                    }              
            }
        return O;
    }

void FixMissingAtoms::repairSer(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {
        V1 N, CA, C, O, CB, OG, Np, CAp, Cp, OXT; 
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "OG") OG = res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            }
        RC require_chi = {"OG"};
        for (auto i: missing_atoms)
            { 
                if (i == "OG")
                    {
                        if(!faspr.empty())
                            {
                               float N_CA_CB_OG_diangle = get_chi("N", "CA", "CB", "OG", res, faspr);
                               OG = calcCoordinate(N, CA, CB, 1.417, 110.773, N_CA_CB_OG_diangle); 
                            }
                        else
                            {
                                float N_CA_CB_OG_diangle = -63.3;
                                OG = calcCoordinate(N, CA, CB, 1.417, 110.773, N_CA_CB_OG_diangle); 
                            }                   
                    }
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if (!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, OG, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "OG", "OXT"};
            }
            else
            {
                res.res_pos = {N, CA, C, O, CB, OG};
                res.res_atom= {"N", "CA", "C", "O", "CB", "OG"};
            }
    }

void FixMissingAtoms::repairAla(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres)
    {   
        V1 N, CA, C, O, CB, Np, CAp, Cp, OXT;
        
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {

                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);

                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if (!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "OXT"};            
            }
        else
            {
                res.res_pos = {N, CA, C, O, CB};
                res.res_atom= {"N", "CA", "C", "O", "CB"};            
            }    

    }

void FixMissingAtoms::repairGly(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres)
    {
        V1 N, CA, C, O, Np, CAp, Cp, OXT;
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                } 
        if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if (!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, OXT};
                res.res_atom= {"N", "CA", "C", "O", "OXT"};           
            }
        else
            {
                res.res_pos = {N, CA, C, O};
                res.res_atom= {"N", "CA", "C", "O"};           
            }
    }

void FixMissingAtoms::repairAsn(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {
        V1 N, CA, C, O, CB, CG, OD1, ND2, Np, CAp, Cp, OXT; 
        float CA_CB_CG_OD1_diangle, CA_CB_CG_ND2_diangle;
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG") CG = res.res_pos[i];
                else if (res.res_atom[i] == "ND2")ND2= res.res_pos[i];
                else if (res.res_atom[i] == "OD1")OD1= res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            } 
        RC require_chi = {"CG", "OD1"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                        if(!faspr.empty())
                            {
                                float CA_CB_CG_OD1_diangle = get_chi("CA", "CB", "CG", "OD1", res, faspr);
                                float N_CA_CB_CG_diangle   = get_chi("N", "CA", "CB", "CG", res, faspr);
                                float CA_CB_CG_ND2_diangle = 180.0 + CA_CB_CG_OD1_diangle;
                                CG  = calcCoordinate(N, CA, CB, 1.52, 112.62, N_CA_CB_CG_diangle);
                                OD1 = calcCoordinate(CA, CB, CG, 1.23, 120.85, CA_CB_CG_OD1_diangle);
                                ND2 = calcCoordinate(CA, CB, CG, 1.33,116.48, CA_CB_CG_ND2_diangle);
                                goto CTER; 
                            }
                        else
                            {
                                float CA_CB_CG_OD1_diangle = -58.3;
                                float N_CA_CB_CG_diangle   = -65.5;
                                float CA_CB_CG_ND2_diangle = 180.0 + CA_CB_CG_OD1_diangle;
                                CG  = calcCoordinate(N, CA, CB, 1.52, 112.62, N_CA_CB_CG_diangle);
                                OD1 = calcCoordinate(CA, CB, CG, 1.23, 120.85, CA_CB_CG_OD1_diangle);
                                ND2 = calcCoordinate(CA, CB, CG, 1.33,116.48, CA_CB_CG_ND2_diangle);
                                goto CTER;                             
                            }
                    }
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "ND2") != missing_atoms.end()))
            {
                CA_CB_CG_OD1_diangle = calcTorsionAngle(CA, CB, CG, OD1);
                CA_CB_CG_ND2_diangle = 180.0 + CA_CB_CG_OD1_diangle;
                ND2 = calcCoordinate(CA, CB, CG, 1.33, 116.48, CA_CB_CG_ND2_diangle);
            }
        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if(!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, CG, OD1, ND2, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2", "OXT"};         
            }
        else 
            {
                res.res_pos = {N, CA, C, O, CB, CG, OD1, ND2};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"};            
            }
    }

void FixMissingAtoms::repairGlu(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {
        V1 N, CA, C, O, CB, CG, CD, OE1, OE2, Np, CAp, Cp, OXT;
        float CB_CG_CD_OE1_diangle, CB_CG_CD_OE2_diangle; 
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG") CG = res.res_pos[i];
                else if (res.res_atom[i] == "CD") CD = res.res_pos[i];
                else if (res.res_atom[i] == "OE1")OE1= res.res_pos[i];
                else if (res.res_atom[i] == "OE2")OE2= res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            } 

        RC require_chi = {"CG", "CD", "OE1"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                        if(!faspr.empty())
                            {
                                float CB_CG_CD_OE1_diangle = get_chi("CB","CG","CD","OE1",res,faspr);
                                float CA_CB_CG_CD_diangle  = get_chi("CA","CB","CG","CD",res,faspr);
                                float N_CA_CB_CG_diangle   = get_chi("N","CA","CB","CG",res,faspr);
                                float CB_CG_CD_OE2_diangle = 180.0 + CB_CG_CD_OE1_diangle;
                                CG  = calcCoordinate(N, CA, CB, 1.52, 113.82,N_CA_CB_CG_diangle);
                                CD  = calcCoordinate(CA, CB, CG,1.52, 119.02,  CA_CB_CG_CD_diangle);
                                OE1 = calcCoordinate(CB, CG, CD, 1.25, 119.02, CB_CG_CD_OE1_diangle);
                                OE2 = calcCoordinate(CB, CG, CD, 1.25, 118.08,CB_CG_CD_OE2_diangle);
                                goto CTER;
                            }
                        else
                            {
                                float CB_CG_CD_OE1_diangle = -6.2;
                                float CA_CB_CG_CD_diangle  = -179.8;
                                float N_CA_CB_CG_diangle   = -63.8;
                                float CB_CG_CD_OE2_diangle = 180.0 + CB_CG_CD_OE1_diangle;
                                CG  = calcCoordinate(N, CA, CB, 1.52, 113.82,N_CA_CB_CG_diangle);
                                CD  = calcCoordinate(CA, CB, CG,1.52, 119.02,  CA_CB_CG_CD_diangle);
                                OE1 = calcCoordinate(CB, CG, CD, 1.25, 119.02, CB_CG_CD_OE1_diangle);
                                OE2 = calcCoordinate(CB, CG, CD, 1.25, 118.08,CB_CG_CD_OE2_diangle);  
                                goto CTER;                           
                            }
                    }
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "OE2") != missing_atoms.end()))
            {
                CB_CG_CD_OE1_diangle = calcTorsionAngle(CB, CG, CD, OE1);
                CB_CG_CD_OE2_diangle = 180.0 + CB_CG_CD_OE1_diangle;
                OE2 = calcCoordinate(CB, CG, CD, 1.25, 118.08,CB_CG_CD_OE2_diangle);
            }
        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if (!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD, OE1, OE2, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2", "OXT"}; 
            }
        else
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD, OE1, OE2};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"};             
            }
    }

void FixMissingAtoms::repairAsp(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {
        V1 N, CA, C, O, CB, CG, CD, OD1, OD2, Np, CAp, Cp, OXT; 
        float N_CA_CB_CG_diangle, CA_CB_CG_OD1_diangle, CA_CB_CG_OD2_diangle;
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG") CG = res.res_pos[i];
                else if (res.res_atom[i] == "OD1")OD1= res.res_pos[i];
                else if (res.res_atom[i] == "OD2")OD2= res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            } 
        RC require_chi = {"CG", "OD1"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                        if(!faspr.empty())
                            {
                                CA_CB_CG_OD1_diangle = get_chi("CA","CB","CG","OD1",res,faspr);
                                N_CA_CB_CG_diangle   = get_chi("N","CA","CB","CG",res,faspr);
                                CA_CB_CG_OD2_diangle = 180 + CA_CB_CG_OD1_diangle;
                                CG  = calcCoordinate(N, CA, CB, 1.52, 113.06, N_CA_CB_CG_diangle);
                                OD1 = calcCoordinate(CA, CB, CG, 1.25, 119.22, CA_CB_CG_OD1_diangle);
                                OD2 = calcCoordinate(CA, CB, CG, 1.25, 118.21, CA_CB_CG_OD2_diangle); 
                                goto CTER;                           
                            }
                        else
                            {
                                CA_CB_CG_OD1_diangle = -46.7;
                                N_CA_CB_CG_diangle   = -66.4;
                                CA_CB_CG_OD2_diangle = 180 + CA_CB_CG_OD1_diangle;
                                CG  = calcCoordinate(N, CA, CB, 1.52, 113.06, N_CA_CB_CG_diangle);
                                OD1 = calcCoordinate(CA, CB, CG, 1.25, 119.22, CA_CB_CG_OD1_diangle);
                                OD2 = calcCoordinate(CA, CB, CG, 1.25, 118.21, CA_CB_CG_OD2_diangle);
                                goto CTER;                            
                            }
                    }
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "OD2") != missing_atoms.end()))
            {
                CA_CB_CG_OD1_diangle = calcTorsionAngle(CA, CB, CG, OD1);
                CA_CB_CG_OD2_diangle = 180 + CA_CB_CG_OD1_diangle;
                OD2 = calcCoordinate(CA, CB, CG, 1.25, 118.21, CA_CB_CG_OD2_diangle);
            }
        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if(!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, CG, OD1, OD2, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2", "OXT"}; 
            }
        else
            {
                res.res_pos = {N, CA, C, O, CB, CG, OD1, OD2};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"};             
            }
    }

void FixMissingAtoms::repairGln(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {
        V1 N, CA, C, O, CB, CG, CD, OE1, NE2, Np, CAp, Cp, OXT; 
        float N_CA_CB_CG_diangle, CB_CG_CD_OE1_diangle, CA_CB_CG_CD_diangle, CB_CG_CD_NE2_diangle;
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG") CG = res.res_pos[i];
                else if (res.res_atom[i] == "CD") CD = res.res_pos[i];
                else if (res.res_atom[i] == "OE1")OE1= res.res_pos[i];
                else if (res.res_atom[i] == "NE2")NE2= res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            }
        RC require_chi = {"CG", "CD", "OE1"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                       if(!faspr.empty())
                           {
                               N_CA_CB_CG_diangle   = get_chi("N","CA","CB","CG",res,faspr);
                               CB_CG_CD_OE1_diangle = get_chi("CB","CG","CD","OE1",res,faspr);
                               CA_CB_CG_CD_diangle  = get_chi("CA","CB","CG","CD",res,faspr);
                               CB_CG_CD_NE2_diangle = 180.0 + CB_CG_CD_OE1_diangle;
                               CG  = calcCoordinate(N, CA, CB, 1.52, 113.75,N_CA_CB_CG_diangle);
                               CD  = calcCoordinate(CA, CB, CG, 1.52, 112.78, CA_CB_CG_CD_diangle);
                               OE1 = calcCoordinate(CB, CG, CD, 1.24, 120.86, CB_CG_CD_OE1_diangle);
                               NE2 = calcCoordinate(CB, CG, CD, 1.33, 116.50,CB_CG_CD_NE2_diangle);
                               goto CTER;
                           }
                       else
                           {
                               N_CA_CB_CG_diangle   = -60.2;
                               CB_CG_CD_OE1_diangle = -50.5;
                               CA_CB_CG_CD_diangle  = -69.6;
                               CB_CG_CD_NE2_diangle = 180.0 + CB_CG_CD_OE1_diangle;
                               CG  = calcCoordinate(N, CA, CB, 1.52, 113.75,N_CA_CB_CG_diangle);
                               CD  = calcCoordinate(CA, CB, CG, 1.52, 112.78, CA_CB_CG_CD_diangle);
                               OE1 = calcCoordinate(CB, CG, CD, 1.24, 120.86, CB_CG_CD_OE1_diangle);
                               NE2 = calcCoordinate(CB, CG, CD, 1.33, 116.50,CB_CG_CD_NE2_diangle);
                               goto CTER;
                           }
                    }
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "NE2") != missing_atoms.end()))
            {
                CB_CG_CD_OE1_diangle = calcTorsionAngle(CB, CG, CD, OE1);
                CB_CG_CD_NE2_diangle = 180.0 + CB_CG_CD_OE1_diangle;
                NE2 = calcCoordinate(CB, CG, CD, 1.33, 116.50, CB_CG_CD_NE2_diangle);
            }
        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if(!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD, OE1, NE2, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2", "OXT"};
            }
        else
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD, OE1, NE2};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"};          
            }
    }

void FixMissingAtoms::repairLys(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {
        V1 N, CA, C, O, CB, CG, CD, CE, NZ, Np, CAp, Cp, OXT; 
        float N_CA_CB_CG_diangle, CA_CB_CG_CD_diangle, CB_CG_CD_CE_diangle, CG_CD_CE_NZ_diangle;
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG") CG = res.res_pos[i];
                else if (res.res_atom[i] == "CD") CD = res.res_pos[i];
                else if (res.res_atom[i] == "CE") CE = res.res_pos[i];
                else if (res.res_atom[i] == "NZ") NZ = res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            }
        RC require_chi = {"CG", "CD", "CE", "NZ"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                        if(!faspr.empty())
                            {
                                N_CA_CB_CG_diangle  = get_chi("N","CA","CB","CG",res,faspr);
                                CA_CB_CG_CD_diangle = get_chi("CA","CB","CG","CD",res,faspr);
                                CB_CG_CD_CE_diangle = get_chi("CB","CG","CD","CE",res,faspr);
                                CG_CD_CE_NZ_diangle = get_chi("CG","CD","CE","NZ",res,faspr);                        
                                CG = calcCoordinate(N, CA, CB, 1.52, 113.83, N_CA_CB_CG_diangle);
                                CD = calcCoordinate(CA, CB, CG, 1.52, 111.79, CA_CB_CG_CD_diangle);
                                CE = calcCoordinate(CB, CG, CD, 1.46, 111.68, CB_CG_CD_CE_diangle);
                                NZ = calcCoordinate(CG, CD, CE, 1.33, 124.79, CG_CD_CE_NZ_diangle);
                            }
                        else
                            {
                                N_CA_CB_CG_diangle  = -64.5;
                                CA_CB_CG_CD_diangle = -178.1;
                                CB_CG_CD_CE_diangle = -179.6;
                                CG_CD_CE_NZ_diangle = 179.6;                        
                                CG = calcCoordinate(N, CA, CB, 1.52, 113.83, N_CA_CB_CG_diangle);
                                CD = calcCoordinate(CA, CB, CG, 1.52, 111.79, CA_CB_CG_CD_diangle);
                                CE = calcCoordinate(CB, CG, CD, 1.46, 111.68, CB_CG_CD_CE_diangle);
                                NZ = calcCoordinate(CG, CD, CE, 1.33, 124.79, CG_CD_CE_NZ_diangle);
                            }
                    }
            }
        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if(!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD, CE, NZ, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ","OXT"};
            }
        else
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD, CE, NZ};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"};            
            }
    }

void FixMissingAtoms::repairArg(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {
        V1 N, CA, C, O, CB, CG, CD, NE, CZ, NH1, NH2, Np, CAp, Cp, OXT; 
        float N_CA_CB_CG_diangle, CA_CB_CG_CD_diangle, CB_CG_CD_NE_diangle, CG_CD_NE_CZ_diangle;
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG") CG = res.res_pos[i];
                else if (res.res_atom[i] == "CD") CD = res.res_pos[i];
                else if (res.res_atom[i] == "NE") NE = res.res_pos[i];
                else if (res.res_atom[i] == "CZ") CZ = res.res_pos[i];
                else if (res.res_atom[i] == "NH1")NH1= res.res_pos[i];
                else if (res.res_atom[i] == "NH2")NH2= res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            }
        RC require_chi = {"CG", "CD", "NE", "CZ"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                        if(!faspr.empty())
                            {
                                N_CA_CB_CG_diangle  = get_chi("N","CA","CB","CG",res,faspr);
                                CA_CB_CG_CD_diangle = get_chi("CA","CB","CG","CD",res,faspr);
                                CB_CG_CD_NE_diangle = get_chi("CB","CG","CD","NE",res,faspr);
                                CG_CD_NE_CZ_diangle = get_chi("CG","CD","NE","CZ",res,faspr);
                                CG  = calcCoordinate(N, CA, CB, 1.52, 113.83, N_CA_CB_CG_diangle);
                                CD  = calcCoordinate(CA, CB, CG, 1.52, 111.79, CA_CB_CG_CD_diangle);
                                NE  = calcCoordinate(CB, CG, CD, 1.46, 111.68, CB_CG_CD_NE_diangle);
                                CZ  = calcCoordinate(CG, CD, NE, 1.33, 124.79, CG_CD_NE_CZ_diangle);
                                NH1 = calcCoordinate(CD, NE, CZ, 1.33, 120.64, 0.0);
                                NH2 = calcCoordinate(CD, NE, CZ, 1.33, 119.63, 180.0);
                                goto CTER;
                            }
                        else
                            {
                                N_CA_CB_CG_diangle  = -65.2;
                                CA_CB_CG_CD_diangle = -179.2;
                                CB_CG_CD_NE_diangle = -179.3;
                                CG_CD_NE_CZ_diangle = -178.7;
                                CG  = calcCoordinate(N, CA, CB, 1.52, 113.83, N_CA_CB_CG_diangle);
                                CD  = calcCoordinate(CA, CB, CG, 1.52, 111.79, CA_CB_CG_CD_diangle);
                                NE  = calcCoordinate(CB, CG, CD, 1.46, 111.68, CB_CG_CD_NE_diangle);
                                CZ  = calcCoordinate(CG, CD, NE, 1.33, 124.79, CG_CD_NE_CZ_diangle);
                                NH1 = calcCoordinate(CD, NE, CZ, 1.33, 120.64, 0.0);
                                NH2 = calcCoordinate(CD, NE, CZ, 1.33, 119.63, 180.0);
                                goto CTER;
                            }
                    }
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "NH1") != missing_atoms.end()))
            {
                NH1 = calcCoordinate(CD, NE, CZ, 1.33, 120.64, 0.0);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "NH2") != missing_atoms.end()))
            {
                NH2 = calcCoordinate(CD, NE, CZ, 1.33, 119.63, 180.0);
            }

        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if(!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD, NE, CZ, NH1, NH2, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2", "OXT"};            
            }
        else
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD, NE, CZ, NH1, NH2};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"};            
            }
    }

void FixMissingAtoms::repairTyr(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {
        V1 N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, OH, Np, CAp, Cp, OXT;
        float CA_CB_CG_CD1_diangle, N_CA_CB_CG_diangle, CA_CB_CG_CD2_diangle;
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG") CG = res.res_pos[i];
                else if (res.res_atom[i] == "CD1")CD1= res.res_pos[i];
                else if (res.res_atom[i] == "CD2")CD2= res.res_pos[i];
                else if (res.res_atom[i] == "CE1")CE1= res.res_pos[i];
                else if (res.res_atom[i] == "CE2")CE2= res.res_pos[i];
                else if (res.res_atom[i] == "CZ")  CZ= res.res_pos[i];
                else if (res.res_atom[i] == "OH")  OH= res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            }
        RC require_chi = {"CG", "CD1"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                        if(!faspr.empty())
                            {
                                CA_CB_CG_CD1_diangle = get_chi("CA","CB","CG","CD1",res,faspr);
                                N_CA_CB_CG_diangle   = get_chi("N","CA","CB","CG",res,faspr);
                                CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle + 180.0;
                                CG  = calcCoordinate(N, CA, CB, 1.51, 113.8,N_CA_CB_CG_diangle);
                                CD1 = calcCoordinate(CA, CB, CG, 1.39, 120.98, CA_CB_CG_CD1_diangle);
                                CD2 = calcCoordinate(CA, CB, CG, 1.39, 120.82,CA_CB_CG_CD2_diangle);
                                CE1 = calcCoordinate(CB, CG, CD1, 1.39, 120.0,180.0);
                                CE2 = calcCoordinate(CB, CG, CD2, 1.39, 120.0, 180.0);
                                CZ  = calcCoordinate(CG, CD1, CE1, 1.39, 120.0, 0.0);
                                OH  = calcCoordinate(CD1, CE1, CZ, 1.39, 119.78, 180.0);
                                goto CTER;
                            }
                        else
                            {
                                CA_CB_CG_CD1_diangle = 93.1;
                                N_CA_CB_CG_diangle   = -64.3;
                                CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle + 180.0;
                                CG  = calcCoordinate(N, CA, CB, 1.51, 113.8,N_CA_CB_CG_diangle);
                                CD1 = calcCoordinate(CA, CB, CG, 1.39, 120.98, CA_CB_CG_CD1_diangle);
                                CD2 = calcCoordinate(CA, CB, CG, 1.39, 120.82,CA_CB_CG_CD2_diangle);
                                CE1 = calcCoordinate(CB, CG, CD1, 1.39, 120.0,180.0);
                                CE2 = calcCoordinate(CB, CG, CD2, 1.39, 120.0, 180.0);
                                CZ  = calcCoordinate(CG, CD1, CE1, 1.39, 120.0, 0.0);
                                OH  = calcCoordinate(CD1, CE1, CZ, 1.39, 119.78, 180.0);
                                goto CTER;                            
                            }
                    }
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CD2") != missing_atoms.end()))
            {
                CA_CB_CG_CD1_diangle = calcTorsionAngle(CA, CB, CG, CD1);
                CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle + 180.0;
                CD2 = calcCoordinate(CA, CB, CG, 1.39, 120.82,CA_CB_CG_CD2_diangle);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CE1") != missing_atoms.end()))
            {
                CE1 = calcCoordinate(CB, CG, CD1, 1.39, 120.0,180.0);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CE2") != missing_atoms.end()))
            {
                CE2 = calcCoordinate(CB, CG, CD2, 1.39, 120.0, 180.0);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CZ") != missing_atoms.end()))
            {
                CZ = calcCoordinate(CG, CD1, CE1, 1.39, 120.0, 0.0);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "OH") != missing_atoms.end()))
            {
                OH = calcCoordinate(CD1, CE1, CZ, 1.39, 119.78, 180.0);
            }
        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if(!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, OH, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "OXT"};            
            }
        else
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, OH};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"};            
            }
    }

void FixMissingAtoms::repairTrp(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {
        V1  N, CA, C, O, CB, CG, CD1, CD2, NE1, OXT, CE2, CE3, CZ2, CZ3, CH2, Np, CAp, Cp; 
        float N_CA_CB_CG_diangle, CA_CB_CG_CD1_diangle, CA_CB_CG_CD2_diangle;
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG") CG = res.res_pos[i];
                else if (res.res_atom[i] == "CD1")CD1= res.res_pos[i];
                else if (res.res_atom[i] == "CD2")CD2= res.res_pos[i];
                else if (res.res_atom[i] == "NE1")NE1= res.res_pos[i];
                else if (res.res_atom[i] == "CE2")CE2= res.res_pos[i];
                else if (res.res_atom[i] == "CE3")CE3= res.res_pos[i];
                else if (res.res_atom[i] == "CZ2")CZ2= res.res_pos[i];
                else if (res.res_atom[i] == "CZ3")CZ3= res.res_pos[i];
                else if (res.res_atom[i] == "CH2")CH2= res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            }
        RC require_chi = {"CG", "CD1"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                        if(!faspr.empty())
                            {
                                N_CA_CB_CG_diangle   = get_chi("N","CA","CB","CG",res,faspr);
                                CA_CB_CG_CD1_diangle = get_chi("CA","CB","CG","CD1",res,faspr);
                                CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle - 180.0;
                                CG  = calcCoordinate(N, CA, CB, 1.50, 114.10, N_CA_CB_CG_diangle);
                                CD1 = calcCoordinate(CA, CB, CG, 1.37, 127.07, CA_CB_CG_CD1_diangle);
                                CD2 = calcCoordinate(CA, CB, CG, 1.43, 126.66, CA_CB_CG_CD2_diangle);
                                NE1 = calcCoordinate(CB, CG, CD1, 1.38, 108.5, 180.0);
                                CE2 = calcCoordinate(CB, CG, CD2, 1.40, 108.5, 180.0);
                                CE3 = calcCoordinate(CB, CG, CD2, 1.40, 133.83, 0.0);
                                CZ2 = calcCoordinate(CG, CD2, CE2, 1.40, 120.0, 180.0);
                                CZ3 = calcCoordinate(CG, CD2, CE3, 1.40, 120.0, 180.0);
                                CH2 = calcCoordinate(CD2, CE2, CZ2, 1.40, 120.0, 0.0);
                                goto CTER;
                            }
                        else
                            {
                                N_CA_CB_CG_diangle   = -66.4;
                                CA_CB_CG_CD1_diangle = 96.3;
                                CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle - 180.0;
                                CG  = calcCoordinate(N, CA, CB, 1.50, 114.10, N_CA_CB_CG_diangle);
                                CD1 = calcCoordinate(CA, CB, CG, 1.37, 127.07, CA_CB_CG_CD1_diangle);
                                CD2 = calcCoordinate(CA, CB, CG, 1.43, 126.66, CA_CB_CG_CD2_diangle);
                                NE1 = calcCoordinate(CB, CG, CD1, 1.38, 108.5, 180.0);
                                CE2 = calcCoordinate(CB, CG, CD2, 1.40, 108.5, 180.0);
                                CE3 = calcCoordinate(CB, CG, CD2, 1.40, 133.83, 0.0);
                                CZ2 = calcCoordinate(CG, CD2, CE2, 1.40, 120.0, 180.0);
                                CZ3 = calcCoordinate(CG, CD2, CE3, 1.40, 120.0, 180.0);
                                CH2 = calcCoordinate(CD2, CE2, CZ2, 1.40, 120.0, 0.0);
                                goto CTER;
                            }
                    }
            }

        if ((find(missing_atoms.begin(),missing_atoms.end(), "CD2") != missing_atoms.end()))
            {
                CA_CB_CG_CD1_diangle = calcTorsionAngle(CA, CB, CG, CD1);
                CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle - 180.0;

                CD2 = calcCoordinate(CA, CB, CG, 1.43, 126.66,CA_CB_CG_CD2_diangle);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "NE1") != missing_atoms.end()))
            {
                NE1 = calcCoordinate(CB, CG, CD1, 1.38, 108.5, 180.0);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CE2") != missing_atoms.end()))
            {
                CE2 = calcCoordinate(CB, CG, CD2, 1.40, 108.5, 180.0);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CE3") != missing_atoms.end()))
            {
                CE3 = calcCoordinate(CB, CG, CD2, 1.40, 133.83, 0.);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CZ2") != missing_atoms.end()))
            {
                CZ2 = calcCoordinate(CG, CD2, CE2, 1.40, 120.0, 180.0);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CZ3") != missing_atoms.end()))
            {
                CZ3 = calcCoordinate(CG, CD2, CE3, 1.40, 120.0, 180.0);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CH2") != missing_atoms.end()))
            {
                CH2 = calcCoordinate(CD2, CE2, CZ2, 1.40, 120.0, 0.0);
            }
        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if(!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2", "OXT"};            
            }
        else
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"};            
            }
    }

void FixMissingAtoms::repairPhe(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {
        V1  N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, Np, CAp, Cp, OXT; 
        float N_CA_CB_CG_diangle, CA_CB_CG_CD1_diangle, CA_CB_CG_CD2_diangle;
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG") CG = res.res_pos[i];
                else if (res.res_atom[i] == "CD1")CD1= res.res_pos[i];
                else if (res.res_atom[i] == "CD2")CD2= res.res_pos[i];
                else if (res.res_atom[i] == "CE1")CE1= res.res_pos[i];
                else if (res.res_atom[i] == "CE2")CE2= res.res_pos[i];
                else if (res.res_atom[i] == "CZ")  CZ= res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            }
        RC require_chi = {"CG", "CD1"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                        if(!faspr.empty())
                            {
                                N_CA_CB_CG_diangle   = get_chi("N","CA","CB","CG",res,faspr);
                                CA_CB_CG_CD1_diangle = get_chi("CA","CB","CG","CD1",res,faspr);
                                CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle - 180.0;
                                CG  = calcCoordinate(N, CA, CB, 1.50, 113.85, N_CA_CB_CG_diangle);
                                CD1 = calcCoordinate(CA, CB, CG, 1.39, 120.0, CA_CB_CG_CD1_diangle);
                                CD2 = calcCoordinate(CA, CB, CG, 1.39, 120.0, CA_CB_CG_CD2_diangle);
                                CE1 = calcCoordinate(CB, CG, CD1, 1.39, 120.0, 180.0);
                                CE2 = calcCoordinate(CB, CG, CD2, 1.39, 120.0, 180.0);
                                CZ  = calcCoordinate(CG, CD1, CE1, 1.39, 120.0, 0.0);
                                goto CTER;
                            }
                        else
                            {
                                N_CA_CB_CG_diangle   = -64.7;
                                CA_CB_CG_CD1_diangle = 93.3;
                                CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle - 180.0;
                                CG  = calcCoordinate(N, CA, CB, 1.50, 113.85, N_CA_CB_CG_diangle);
                                CD1 = calcCoordinate(CA, CB, CG, 1.39, 120.0, CA_CB_CG_CD1_diangle);
                                CD2 = calcCoordinate(CA, CB, CG, 1.39, 120.0, CA_CB_CG_CD2_diangle);
                                CE1 = calcCoordinate(CB, CG, CD1, 1.39, 120.0, 180.0);
                                CE2 = calcCoordinate(CB, CG, CD2, 1.39, 120.0, 180.0);
                                CZ  = calcCoordinate(CG, CD1, CE1, 1.39, 120.0, 0.0);
                                goto CTER;
                            }
                    }
                }

        if ((find(missing_atoms.begin(),missing_atoms.end(), "CD2") != missing_atoms.end()))
            {
                CA_CB_CG_CD1_diangle = calcTorsionAngle(CA, CB, CG, CD1);
                CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle - 180.0;
                CD2 = calcCoordinate(CA, CB, CG, 1.39, 120.0,CA_CB_CG_CD2_diangle);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CE1") != missing_atoms.end()))
            {
                CE1 = calcCoordinate(CB, CG, CD1, 1.39, 120.0, 180.0);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CE2") != missing_atoms.end()))
            {
                CE2 = calcCoordinate(CB, CG, CD2, 1.39, 120.0, 180.0);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CZ") != missing_atoms.end()))
            {
                CZ = calcCoordinate(CG, CD1, CE1, 1.39, 120.0, 0.0);
            }
        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if(!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OXT"};            
            }
        else
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"};            
            }
    }

void FixMissingAtoms::repairHis(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {
        V1  N, CA, C, O, CB, CG, ND1, CD2, CE1, NE2, Np, CAp, Cp, OXT; 
        float CA_CB_CG_ND1_diangle, N_CA_CB_CG_diangle, CA_CB_CG_CD2_diangle;
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG") CG = res.res_pos[i];
                else if (res.res_atom[i] == "ND1")ND1= res.res_pos[i];
                else if (res.res_atom[i] == "CD2")CD2= res.res_pos[i];
                else if (res.res_atom[i] == "CE1")CE1= res.res_pos[i];
                else if (res.res_atom[i] == "NE2")NE2= res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            }
        RC require_chi = {"CG", "ND1"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                        if(!faspr.empty())
                            {
                                CA_CB_CG_ND1_diangle = get_chi("CA","CB","CG","ND1",res,faspr);
                                N_CA_CB_CG_diangle   = get_chi("N","CA","CB","CG",res,faspr);
                                CA_CB_CG_CD2_diangle = 180.0 + CA_CB_CG_ND1_diangle;
                                CG  = calcCoordinate(N, CA, CB, 1.49, 113.74, N_CA_CB_CG_diangle);
                                ND1 = calcCoordinate(CA, CB, CG, 1.38, 122.85, CA_CB_CG_ND1_diangle);
                                CD2 = calcCoordinate(CA, CB, CG, 1.35, 130.61, CA_CB_CG_CD2_diangle);
                                CE1 = calcCoordinate(CB, CG, ND1, 1.32, 108.5, 180.0);
                                NE2 = calcCoordinate(CB, CG, CD2, 1.35, 108.5, 180.0);
                                goto CTER;
                            }
                        else
                            {
                                CA_CB_CG_ND1_diangle = -75.7;
                                N_CA_CB_CG_diangle   = -63.2;
                                CA_CB_CG_CD2_diangle = 180.0 + CA_CB_CG_ND1_diangle;
                                CG  = calcCoordinate(N, CA, CB, 1.49, 113.74, N_CA_CB_CG_diangle);
                                ND1 = calcCoordinate(CA, CB, CG, 1.38, 122.85, CA_CB_CG_ND1_diangle);
                                CD2 = calcCoordinate(CA, CB, CG, 1.35, 130.61, CA_CB_CG_CD2_diangle);
                                CE1 = calcCoordinate(CB, CG, ND1, 1.32, 108.5, 180.0);
                                NE2 = calcCoordinate(CB, CG, CD2, 1.35, 108.5, 180.0);
                                goto CTER;
                            }
                    }
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CD2") != missing_atoms.end()))
            {
                CA_CB_CG_ND1_diangle = calcTorsionAngle(CA, CB, CG, ND1);
                CA_CB_CG_CD2_diangle = 180.0 + CA_CB_CG_ND1_diangle;
                CD2 = calcCoordinate(CA, CB, CG, 1.35, 130.61, CA_CB_CG_CD2_diangle);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CE1") != missing_atoms.end()))
            {
                CE1 = calcCoordinate(CB, CG, ND1, 1.32, 108.5, 180.0);
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "NE2") != missing_atoms.end()))
            {
                NE2 = calcCoordinate(CB, CG, CD2, 1.35, 108.5, 180.0);
            }
        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if(!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, CG, ND1, CD2, CE1, NE2, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2", "OXT"};            
            }
        else
            {
                res.res_pos = {N, CA, C, O, CB, CG, ND1, CD2, CE1, NE2};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"};            
            }
    }

void FixMissingAtoms::repairPro(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {   
        V1  N, CA, C, O, CB, CG, CD, Np, CAp, Cp, OXT;
        float N_CA_CB_CG_diangle, CA_CB_CG_CD_diangle; 
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG") CG = res.res_pos[i];
                else if (res.res_atom[i] == "CD") CD = res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            }
        RC require_chi = {"CG", "CD"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                        if(!faspr.empty())
                            {
                                N_CA_CB_CG_diangle  = get_chi("N","CA","CB","CG",res,faspr);
                                CA_CB_CG_CD_diangle = get_chi("CA","CB","CG","CD",res,faspr);
                                CG = calcCoordinate(N, CA, CB, 1.49, 104.21, N_CA_CB_CG_diangle);
                                CD = calcCoordinate(CA, CB, CG, 1.50, 105.03, CA_CB_CG_CD_diangle);
                            }
                        else
                            {
                                N_CA_CB_CG_diangle  = 29.6;
                                CA_CB_CG_CD_diangle = -34.8;
                                CG = calcCoordinate(N, CA, CB, 1.49, 104.21, N_CA_CB_CG_diangle);
                                CD = calcCoordinate(CA, CB, CG, 1.50, 105.03, CA_CB_CG_CD_diangle);
                            }
                    }
            }
        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if(!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD", "OXT"};            
            }
        else
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD"};            
            }
    }

void FixMissingAtoms::repairMet(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {
        V1  N, CA, C, O, CB, CG, SD, CE, Np, CAp, Cp, OXT; 
        float N_CA_CB_CG_diangle, CA_CB_CG_SD_diangle, CB_CG_SD_CE_diangle;
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG") CG = res.res_pos[i];
                else if (res.res_atom[i] == "SD") SD = res.res_pos[i];
                else if (res.res_atom[i] == "CE") CE = res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            }
        RC require_chi = {"CG", "SD","CE"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                        if(!faspr.empty())
                            {
                                N_CA_CB_CG_diangle  = get_chi("N","CA","CB","CG",res,faspr);
                                CA_CB_CG_SD_diangle = get_chi("CA","CB","CG","SD",res,faspr);
                                CB_CG_SD_CE_diangle = get_chi("CB","CG","SD","CE",res,faspr);
                                CG = calcCoordinate(N, CA, CB, 1.52, 113.68, N_CA_CB_CG_diangle);
                                SD = calcCoordinate(CA, CB, CG, 1.81, 112.69, CA_CB_CG_SD_diangle);
                                CE = calcCoordinate(CB, CG, SD, 1.79, 100.61, CB_CG_SD_CE_diangle);
                            }
                        else
                            {
                                N_CA_CB_CG_diangle  = -64.4;
                                CA_CB_CG_SD_diangle = -179.6;
                                CB_CG_SD_CE_diangle = 70.1;
                                CG = calcCoordinate(N, CA, CB, 1.52, 113.68, N_CA_CB_CG_diangle);
                                SD = calcCoordinate(CA, CB, CG, 1.81, 112.69, CA_CB_CG_SD_diangle);
                                CE = calcCoordinate(CB, CG, SD, 1.79, 100.61, CB_CG_SD_CE_diangle);
                            }
                    }
            }
        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if(!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, CG, SD, CE, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "SD", "CE", "OXT"};            
            }
        else
            {
                res.res_pos = {N, CA, C, O, CB, CG, SD, CE};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "SD", "CE"};            
            }
    }

void FixMissingAtoms::repairIle(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {
        V1  N, CA, C, O, CB, CG1, CG2, CD1, Np, CAp, Cp, OXT; 
        float N_CA_CB_CG1_diangle, CA_CB_CG1_CD1_diangle, N_CA_CB_CG2_diangle;
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG1")CG1= res.res_pos[i];
                else if (res.res_atom[i] == "CD1")CD1= res.res_pos[i];
                else if (res.res_atom[i] == "CG2")CG2= res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            }
        RC require_chi = {"CG1", "CD1"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                        if(!faspr.empty())
                            {
                                N_CA_CB_CG1_diangle   = get_chi("N","CA","CB","CG1",res,faspr);
                                CA_CB_CG1_CD1_diangle = get_chi("CA","CB","CG1","CD1",res,faspr);
                                N_CA_CB_CG2_diangle   = N_CA_CB_CG1_diangle-120.0;
                                CG1 = calcCoordinate(N, CA, CB, 1.527, 110.7, N_CA_CB_CG1_diangle);
                                CG2 = calcCoordinate(N, CA, CB, 1.527, 110.4, N_CA_CB_CG2_diangle);
                                CD1 = calcCoordinate(CA, CB, CG1, 1.52, 113.97, CA_CB_CG1_CD1_diangle);
                                goto CTER;
                            }
                        else
                            {
                                N_CA_CB_CG1_diangle   = 59.7;
                                CA_CB_CG1_CD1_diangle = 169.8;
                                N_CA_CB_CG2_diangle   = N_CA_CB_CG1_diangle-120.0;
                                CG1 = calcCoordinate(N, CA, CB, 1.527, 110.7, N_CA_CB_CG1_diangle);
                                CG2 = calcCoordinate(N, CA, CB, 1.527, 110.4, N_CA_CB_CG2_diangle);
                                CD1 = calcCoordinate(CA, CB, CG1, 1.52, 113.97, CA_CB_CG1_CD1_diangle);
                                goto CTER;
                            }
                    }
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CG2") != missing_atoms.end()))
            {
                N_CA_CB_CG1_diangle = calcTorsionAngle(N, CA, CB, CG1);
                N_CA_CB_CG2_diangle = N_CA_CB_CG1_diangle-120.0;
                CG2 = calcCoordinate(N, CA, CB, 1.527, 110.4, N_CA_CB_CG2_diangle);
            }

        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if(!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, CG1, CG2, CD1, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1", "OXT"};            
            }
        else
            {
                res.res_pos = {N, CA, C, O, CB, CG1, CG2, CD1};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"};            
            }
    }

void FixMissingAtoms::repairLeu(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {
        V1  N, CA, C, O, CB, CG, CD1, CD2, Np, CAp, Cp, OXT;
        float N_CA_CB_CG_diangle, CA_CB_CG_CD1_diangle, CA_CB_CG_CD2_diangle;
        
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG") CG = res.res_pos[i];
                else if (res.res_atom[i] == "CD1")CD1= res.res_pos[i];
                else if (res.res_atom[i] == "CD2")CD2= res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            }
        RC require_chi = {"CG", "CD1"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                        if(!faspr.empty())
                            {
                                N_CA_CB_CG_diangle   = get_chi("N","CA","CB","CG",res,faspr);
                                CA_CB_CG_CD1_diangle = get_chi("CA","CB","CG","CD1",res,faspr);
                                CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle + 120.0;
                                CG  = calcCoordinate(N, CA, CB, 1.53, 116.10, N_CA_CB_CG_diangle);
                                CD1 = calcCoordinate(CA, CB, CG, 1.524, 112.50,CA_CB_CG_CD1_diangle);
                                CD2 = calcCoordinate(CA, CB, CG, 1.525, 112.50, CA_CB_CG_CD2_diangle);
                                goto CTER;
                            }
                        else
                            {
                                N_CA_CB_CG_diangle   = -60.1;
                                CA_CB_CG_CD1_diangle = 174.9;
                                CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle + 120.0;
                                CG  = calcCoordinate(N, CA, CB, 1.53, 116.10, N_CA_CB_CG_diangle);
                                CD1 = calcCoordinate(CA, CB, CG, 1.524, 112.50,CA_CB_CG_CD1_diangle);
                                CD2 = calcCoordinate(CA, CB, CG, 1.525, 112.50, CA_CB_CG_CD2_diangle);
                                goto CTER;
                            }
                    }
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CD2") != missing_atoms.end()))
            {
                CA_CB_CG_CD1_diangle = calcTorsionAngle(CA, CB, CG, CD1);
                CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle+120.0;
                CD2 = calcCoordinate(CA, CB, CG, 1.525, 112.50, CA_CB_CG_CD2_diangle);
            }
        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if(!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD1, CD2, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "OXT"};            
            }
        else
            {
                res.res_pos = {N, CA, C, O, CB, CG, CD1, CD2};
                res.res_atom= {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"};            
            }
    }

void FixMissingAtoms::repairVal(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {
        V1  N, CA, C, O, CB, CG1, CG2, Np, CAp, Cp, OXT; 
        float N_CA_CB_CG1_diangle, N_CA_CB_CG2_diangle;
        
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "CG1")CG1= res.res_pos[i];
                else if (res.res_atom[i] == "CG2")CG2= res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            }
        RC require_chi = {"CG1"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                        if(!faspr.empty())
                            {
                                N_CA_CB_CG1_diangle = get_chi("N","CA","CB","CG1",res,faspr);
                                N_CA_CB_CG2_diangle = N_CA_CB_CG1_diangle + 120;
                                CG1 = calcCoordinate(N, CA, CB, 1.527, 110.7, N_CA_CB_CG1_diangle);
                                CG2 = calcCoordinate(N, CA, CB, 1.527, 110.4, N_CA_CB_CG2_diangle);
                                goto CTER;
                            }
                        else
                            {
                                N_CA_CB_CG1_diangle = 177.2;
                                N_CA_CB_CG2_diangle = N_CA_CB_CG1_diangle + 120;
                                CG1 = calcCoordinate(N, CA, CB, 1.527, 110.7, N_CA_CB_CG1_diangle);
                                CG2 = calcCoordinate(N, CA, CB, 1.527, 110.4, N_CA_CB_CG2_diangle);
                                goto CTER;
                            }
                    }
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CG2") != missing_atoms.end()))
            {
                N_CA_CB_CG1_diangle = calcTorsionAngle(N, CA, CB, CG1);
                N_CA_CB_CG2_diangle = N_CA_CB_CG1_diangle + 120.0;
                CG2 = calcCoordinate(N, CA, CB, 1.527, 110.4, N_CA_CB_CG2_diangle);
            }
        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if (!OXT.empty())
        {
            res.res_pos = {N, CA, C, O, CB, CG1, CG2, OXT};
            res.res_atom= {"N", "CA", "C", "O", "CB", "CG1", "CG2","OXT"};
        }
        else
        {
            res.res_pos = {N, CA, C, O, CB, CG1, CG2};
            res.res_atom= {"N", "CA", "C", "O", "CB", "CG1", "CG2"};
        }
    }

void FixMissingAtoms::repairThr(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {
        V1  N, CA, C, O, CB, OG1, CG2, Np, CAp, Cp, OXT;
        float N_CA_CB_OG1_diangle, N_CA_CB_CG2_diangle;
        
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "OG1")OG1= res.res_pos[i];
                else if (res.res_atom[i] == "CG2")CG2= res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            }

        RC require_chi = {"OG1"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                        if(!faspr.empty())
                            {
                                N_CA_CB_OG1_diangle = get_chi("N","CA","CB","OG1",res,faspr);
                                N_CA_CB_CG2_diangle = N_CA_CB_OG1_diangle - 120.0;
                                OG1 = calcCoordinate(N, CA, CB, 1.43, 109.18, N_CA_CB_OG1_diangle);
                                CG2 = calcCoordinate(N, CA, CB, 1.53, 111.13, N_CA_CB_CG2_diangle);
                                goto CTER;
                            }
                        else
                            {
                                N_CA_CB_OG1_diangle = 60.0;
                                N_CA_CB_CG2_diangle = N_CA_CB_OG1_diangle - 120.0;
                                OG1 = calcCoordinate(N, CA, CB, 1.43, 109.18, N_CA_CB_OG1_diangle);
                                CG2 = calcCoordinate(N, CA, CB, 1.53, 111.13, N_CA_CB_CG2_diangle);
                                goto CTER;
                            }
                    }
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CG2") != missing_atoms.end()))
            {
                N_CA_CB_OG1_diangle = calcTorsionAngle(N, CA, CB, OG1);
                N_CA_CB_CG2_diangle = N_CA_CB_OG1_diangle - 120.0;
                CG2 = calcCoordinate(N, CA, CB, 1.53, 111.13, N_CA_CB_CG2_diangle);
            }
        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if(!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, OG1, CG2, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "OG1", "CG2", "OXT"};            
            }
        else
            {
                res.res_pos = {N, CA, C, O, CB, OG1, CG2};
                res.res_atom= {"N", "CA", "C", "O", "CB", "OG1", "CG2"};            
            }
    }

void FixMissingAtoms::repairCys(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr)
    {  
        V1  N, CA, C, O, CB, SG, Np, CAp, Cp, OXT;
        float N_CA_CB_SG_diangle; 
        for (int i = 0; i<res.res_atom.size(); i++)
            {
                if (res.res_atom[i] == "N")       N  = res.res_pos[i];
                else if (res.res_atom[i] == "CA") CA = res.res_pos[i];
                else if (res.res_atom[i] == "C")  C  = res.res_pos[i];
                else if (res.res_atom[i] == "O")  O  = res.res_pos[i];
                else if (res.res_atom[i] == "CB") CB = res.res_pos[i];
                else if (res.res_atom[i] == "SG") SG = res.res_pos[i];
                else if (res.res_atom[i] == "OXT")OXT= res.res_pos[i];
            }
        for (int i = 0; i<res_p.res_atom.size(); i++)
            {
                if (res_p.res_atom[i] == "N")       Np   = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "CA") CAp  = res_p.res_pos[i];
                else if (res_p.res_atom[i] == "C")  Cp   = res_p.res_pos[i];
            }
        if ((find(missing_atoms.begin(),missing_atoms.end(), "O") != missing_atoms.end()))
                {
                    float psi = calcTorsionAngle(Np, CAp, Cp, N);
                    O = addBackbone(N, CA, C, nxtres, psi);
                }  
        if ((find(missing_atoms.begin(),missing_atoms.end(), "CB") != missing_atoms.end()))
            {
                CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66);
            }
        RC require_chi = {"SG"};
        for (auto i: missing_atoms)
            {
                if ((find(require_chi.begin(),require_chi.end(), i) != require_chi.end()))
                    {
                        if(!faspr.empty())
                            {
                                N_CA_CB_SG_diangle = get_chi("N","CA","CB","SG",res,faspr);
                                SG = calcCoordinate(N, CA, CB, 1.81, 113.82, N_CA_CB_SG_diangle);
                            }
                        else
                            {
                                N_CA_CB_SG_diangle = -62.2;
                                SG = calcCoordinate(N, CA, CB, 1.81, 113.82, N_CA_CB_SG_diangle);
                            }
                    }
            }
        CTER:if ((find(missing_atoms.begin(),missing_atoms.end(), "OXT") != missing_atoms.end()))
            {
                OXT = c_termini(N, CA, C, O);
            }
        if(!OXT.empty())
            {
                res.res_pos = {N, CA, C, O, CB, SG, OXT};
                res.res_atom= {"N", "CA", "C", "O", "CB", "SG", "OXT"};            
            }
        else
            {
                res.res_pos = {N, CA, C, O, CB, SG};
                res.res_atom= {"N", "CA", "C", "O", "CB", "SG"};            
            }
    }