/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
********************************************************************************************************************************/

#include "PotentialEnergy.hpp"
#include "ForcefieldParam.hpp"

V2 GetPotentialEnergy::rotate(V1 pos1, V1 pos2, V1 pos3, V1 pos4, float BL, float di)
    {
        float angle = 109.5;
        V2 opt_h = {}; V1 range = {0., 60., 120., 180., 240., 300.};

        for (auto i : range)
        {
            V1 point_h1 = calcCoordinate(pos1, pos2, pos3, BL, angle, di);
            V1 point_h2 = recalcCoordinate(pos2, pos3, point_h1, BL);
            opt_h.push_back(point_h2);
            di+= i;           
        }
        return opt_h;    
    }

V1 GetPotentialEnergy::optmizeH(V1 pos1, V1 pos2, V1 pos3, V1 pos4, float BL, float di, float H_pch, float H_sig, float H_eps, Residue &res, vector<Residue> &chains)
    {
        vector<float> etot = {};
        V2 opt_h = rotate(pos1, pos2, pos3, pos4, BL, di);

        for (int i=0; i<opt_h.size(); i++)
        {
            float poten_ener = 0.;

            for (int j=0; j<chains.size(); j++)
            {
                if (chains[j].res_num != res.res_num && chains[j].c_id == res.c_id)
                {
                    for (int k=0; k<chains[j].res_atom.size(); k++)
                    {
                        if (find(heay_atoms.begin(),heay_atoms.end(), chains[j].res_atom[k]) != heay_atoms.end())
                        {
                            float dist = _distance(chains[j].res_pos[k], opt_h[i]);
                            if (dist <= 2.5)
                            {
                               float  r_die     = 138.94;
                               float A_pch      = ffparam[chains[j].res_type][chains[j].res_atom[k]][0];
                               float A_sig      = ffparam[chains[j].res_type][chains[j].res_atom[k]][1];
                               float A_eps      = ffparam[chains[j].res_type][chains[j].res_atom[k]][2];
                               float r          = 0.1*dist;
                               float coulomb    = (r_die*(H_pch*A_pch))/r;
                               float vdw_eps    = 4*(sqrt(H_eps*A_eps));
                               float vdw_sigm   = (pow(((H_sig+A_sig)*0.5/r),12))-(pow(((H_sig+A_sig)*0.5/r),6));
                               poten_ener       += coulomb+(vdw_eps*vdw_sigm);
                            }
                        }

                    }
                }   
            }

            etot.push_back(poten_ener) ;           
        }
        if (!etot.empty())
            {
                auto i = find(etot.begin(), etot.end(), *min_element(etot.begin(), etot.end()));
                int index = i - etot.begin();
                return opt_h[index];
            }
        else
        {
            return pos4;
        }
    }


