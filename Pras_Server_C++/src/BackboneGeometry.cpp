/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
********************************************************************************************************************************/

#include "LinearAlgebra.hpp"
#include "BackboneGeometry.hpp"

VVF GetBackboneGeometry::pi_helix()
  {
    VF hb_dist; VVF O_coord, N_coord, all_hb_dist; VVVF all_O, all_N;
    for (int i= 0; i<chains.size(); i++)
      {
        for (int j = 0; j<chains[i].res_atom.size(); j++)
            {
              if      (chains[i].res_atom[j] == "N")N_coord.push_back(chains[i].res_pos[j]);
              else if (chains[i].res_atom[j] == "O")O_coord.push_back(chains[i].res_pos[j]);
            }
        if (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
          {
            if (N_coord.size()>5)
              {
                VVF (N_coord.begin()+5, N_coord.end()).swap(N_coord);
                VVF (O_coord.begin(), O_coord.end()-5).swap(O_coord);
                all_O.push_back(O_coord);
                all_N.push_back(N_coord);
                O_coord.clear();N_coord.clear();              
              }
            else
            {
                all_O.push_back(O_coord);
                all_N.push_back(N_coord);
                O_coord.clear();N_coord.clear(); 
            }
          }
      }
    for (int i= 0; i<all_O.size(); i++)
      {
        if (all_O[i].size() >5)
           {
              for (int j = 0; j<all_O[i].size(); j++)
                {
                  hb_dist.push_back(_distance(all_O[i][j],all_N[i][j]));
                }
              for(auto i: {100., 100., 100., 100., 100.})hb_dist.push_back(i);
              all_hb_dist.push_back(hb_dist);hb_dist.clear();          
           }
        else
           {
              for(int k =0; k< all_O[i].size(); k++)hb_dist.push_back(100.);
              all_hb_dist.push_back(hb_dist);hb_dist.clear();  
           }
      }
    return all_hb_dist;
  }

VVF GetBackboneGeometry::alpha_helix()
  {
    VF hb_dist; VVF O_coord, N_coord, all_hb_dist; VVVF all_O, all_N;
    for (int i= 0; i<chains.size(); i++)
      {
        for (int j = 0; j<chains[i].res_atom.size(); j++)
            {
              if      (chains[i].res_atom[j] == "N")N_coord.push_back(chains[i].res_pos[j]);
              else if (chains[i].res_atom[j] == "O")O_coord.push_back(chains[i].res_pos[j]);
            }
        if (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
          {
            if (N_coord.size()>4)
              {
                VVF (N_coord.begin()+4, N_coord.end()).swap(N_coord);
                VVF (O_coord.begin(), O_coord.end()-4).swap(O_coord);
                all_O.push_back(O_coord);
                all_N.push_back(N_coord);
                O_coord.clear();N_coord.clear();
              }
            else
              {
                all_O.push_back(O_coord);
                all_N.push_back(N_coord);
                O_coord.clear();N_coord.clear();  
              }
          }
      }
    for (int i= 0; i<all_O.size(); i++)
      {
        if (all_O[i].size() >4)
          {
            for (int j = 0; j<all_O[i].size(); j++)
              {
                hb_dist.push_back(_distance(all_O[i][j],all_N[i][j]));
              }
            for(auto i: {100., 100., 100., 100.})hb_dist.push_back(i);
            all_hb_dist.push_back(hb_dist);hb_dist.clear();          
          }
        else
          {
            for(int k =0; k< all_O[i].size(); k++)hb_dist.push_back(100.);
            all_hb_dist.push_back(hb_dist);hb_dist.clear(); 
          }
      }
    return all_hb_dist;
  }

VVF GetBackboneGeometry::_310_helix()
  {
    VF hb_dist; VVF O_coord, N_coord, all_hb_dist; VVVF all_O, all_N;
    for (int i= 0; i<chains.size(); i++)
      {
        for (int j = 0; j<chains[i].res_atom.size(); j++)
            {
              if      (chains[i].res_atom[j] == "N")N_coord.push_back(chains[i].res_pos[j]);
              else if (chains[i].res_atom[j] == "O")O_coord.push_back(chains[i].res_pos[j]);
            }
        if (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
          {
            if (N_coord.size()<= 2)
              {
                cout<<"Cannot assign secondary structure to a dipeptide\n"<<endl;
                exit(1);
              }
            VVF (N_coord.begin()+3, N_coord.end()).swap(N_coord);
            VVF (O_coord.begin(), O_coord.end()-3).swap(O_coord);
            all_O.push_back(O_coord);
            all_N.push_back(N_coord);
            O_coord.clear();N_coord.clear(); 
          }
      }
    for (int i= 0; i<all_O.size(); i++)
      {
        for (int j = 0; j<all_O[i].size(); j++)
          {
            hb_dist.push_back(_distance(all_O[i][j],all_N[i][j]));
          }
        for(auto i: {100., 100., 100.})hb_dist.push_back(i);
        all_hb_dist.push_back(hb_dist);hb_dist.clear();
      }
    return all_hb_dist;
  }

VVF GetBackboneGeometry::beta_turn()
  {
    VF cn_dist; VVF C_coord, N_coord, all_cn_dist; VVVF all_C, all_N;
    for (int i= 0; i<chains.size(); i++)
      {
        for (int j = 0; j<chains[i].res_atom.size(); j++)
            {
              if      (chains[i].res_atom[j] == "N")N_coord.push_back(chains[i].res_pos[j]);
              else if (chains[i].res_atom[j] == "C")C_coord.push_back(chains[i].res_pos[j]);
            }
        if (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
          {
            if (N_coord.size()<= 2)
              {
                cout<<"Cannot assign secondary structure to a dipeptide\n"<<endl;
                exit(1);
              }
            VVF (N_coord.begin()+3, N_coord.end()).swap(N_coord);
            VVF (C_coord.begin(), C_coord.end()-3).swap(C_coord);
            all_C.push_back(C_coord);
            all_N.push_back(N_coord);
            C_coord.clear();N_coord.clear();
          }
      }
    for (int i= 0; i<all_C.size(); i++)
      {
       for (int j = 0; j<all_C[i].size(); j++)
        {
          cn_dist.push_back(_distance(all_C[i][j],all_N[i][j]));
        }
       for(auto i: {100., 100., 100.})cn_dist.push_back(i);
       all_cn_dist.push_back(cn_dist);cn_dist.clear();         
      }
    return all_cn_dist;
  }

VVF GetBackboneGeometry::amideDihedral()
{
    VVF O_coord ; VVF O_coord2 ;VVF C_coord;VVF C_coord2 ;
    VVVF all_O; VVVF all_C; VVVF all_O2; VVVF all_C2;
    for (int i= 0; i<chains.size(); i++)
      {
        for (int j = 0; j<chains[i].res_atom.size(); j++)
            {
              if      (chains[i].res_atom[j] == "C")C_coord.push_back(chains[i].res_pos[j]);
              else if (chains[i].res_atom[j] == "O")O_coord.push_back(chains[i].res_pos[j]);
            }
        if (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
          {
            O_coord2 = O_coord; C_coord2 = C_coord;
            O_coord2.erase(O_coord2.begin()+0);
            C_coord2.erase(C_coord2.begin()+0);
            O_coord.pop_back();
            C_coord.pop_back();
            all_C.push_back(C_coord);all_C2.push_back(C_coord2);
            all_O.push_back(O_coord);all_O2.push_back(O_coord2);
            O_coord.clear();C_coord.clear();O_coord2.clear();C_coord2.clear();
          }
      }
    VVF all_amide_di;
    for (int i= 0; i<all_O.size(); i++)
      {
        VF amide_di = {360.};
        for (int j = 0; j<all_O[i].size(); j++)
          {
            amide_di.push_back(calcTorsionAngle(all_O[i][j], all_C[i][j], all_C2[i][j], all_O2[i][j]));
          }
        all_amide_di.push_back(amide_di);amide_di.clear();
      }
    return all_amide_di;
}

VVF GetBackboneGeometry::psiDihedral()
{
    VVF N_coord ; VVF N_coord2 ;VVF CA_coord;VVF C_coord ;
    VVVF all_N; VVVF all_C; VVVF all_N2; VVVF all_CA;
    for (int i= 0; i<chains.size(); i++)
      {
        for (int j = 0; j<chains[i].res_atom.size(); j++)
            {
              if      (chains[i].res_atom[j] == "N")N_coord.push_back(chains[i].res_pos[j]);
              else if (chains[i].res_atom[j] == "CA")CA_coord.push_back(chains[i].res_pos[j]);
              else if (chains[i].res_atom[j] == "C")C_coord.push_back(chains[i].res_pos[j]);
            }
        if (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
          {
            N_coord2 = N_coord;
            N_coord2.erase(N_coord2.begin()+0);
            N_coord.pop_back();
            CA_coord.pop_back();
            C_coord.pop_back();
            all_C.push_back(C_coord);all_CA.push_back(CA_coord);
            all_N.push_back(N_coord);all_N2.push_back(N_coord2);
            N_coord.clear();C_coord.clear();N_coord2.clear();CA_coord.clear();
          }
      }
    VVF all_psi_di;
    for (int i= 0; i<all_N.size(); i++)
      {
        VF psi_di = {};
        for (int j = 0; j<all_N[i].size(); j++)
          {
            psi_di.push_back(calcTorsionAngle(all_N[i][j], all_CA[i][j], all_C[i][j], all_N2[i][j]));
          }
        psi_di.push_back(360.);
        all_psi_di.push_back(psi_di);psi_di.clear();
      }
    return all_psi_di;
}

VVF GetBackboneGeometry::phiDihedral()
{
    VVF N_coord ; VVF C_coord2 ;VVF CA_coord;VVF C_coord ;
    VVVF all_N; VVVF all_C; VVVF all_C2; VVVF all_CA;
    for (int i= 0; i<chains.size(); i++)
      {
        for (int j = 0; j<chains[i].res_atom.size(); j++)
            {
              if      (chains[i].res_atom[j] == "N")N_coord.push_back(chains[i].res_pos[j]);
              else if (chains[i].res_atom[j] == "CA")CA_coord.push_back(chains[i].res_pos[j]);
              else if (chains[i].res_atom[j] == "C")C_coord.push_back(chains[i].res_pos[j]);
            }
        if (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
          {
            C_coord2 = C_coord;
            C_coord2.erase(C_coord2.begin()+0);
            N_coord.erase(N_coord.begin()+0);
            CA_coord.erase(CA_coord.begin()+0);
            C_coord.pop_back();
            all_C.push_back(C_coord);all_CA.push_back(CA_coord);
            all_N.push_back(N_coord);all_C2.push_back(C_coord2);
            N_coord.clear();C_coord.clear();C_coord2.clear();CA_coord.clear();
          }
      }
    VVF all_phi_di;
    for (int i= 0; i<all_N.size(); i++)
      {
        VF phi_di = {360.};
        for (int j = 0; j<all_N[i].size(); j++)
          {
            phi_di.push_back(calcTorsionAngle(all_C[i][j], all_N[i][j], all_CA[i][j], all_C2[i][j]));
          }
        all_phi_di.push_back(phi_di);phi_di.clear();
      }
    return all_phi_di;
}

VS GetBackboneGeometry::chainID()
{
    VS chain_list;
    for (int i= 0; i<chains.size(); i++)
      {
        if (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
          {
            chain_list.push_back(chains[i].c_id);
          }
      }
    return chain_list;
}

VVS GetBackboneGeometry::resName()
{
  VS res_list; VVS all_res_list;

    for (int i= 0; i<chains.size(); i++)
      {
        res_list.push_back(chains[i].res_type);

        if (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
          {
            all_res_list.push_back(res_list);
            res_list.clear();
          }
      }
  return all_res_list;
}

VVS GetBackboneGeometry::dihedralType()
{
  VS di_type; VVS all_di_type;

  for (int i= 0; i<chains.size(); i++)
  {
      if (chains[i].res_type == "PRO") di_type.push_back("Proline");
      else if (chains[i].res_type == "GLY") di_type.push_back("Glycine");
      else if (chains[(i+1)%chains.size()].res_type == "PRO" && chains[i].c_id == chains[(i+1)%chains.size()].c_id)
        {
          if (i != chains.size()-1)
            {
              di_type.push_back("Pre-Pro");
            }
          else di_type.push_back("General");
        }
      else di_type.push_back("General");
      if (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
        {
          all_di_type.push_back(di_type);
          di_type.clear();
        }    
  }
  return all_di_type;
}

VVS GetBackboneGeometry::resNumber()
{
    VS res_No; VVS all_res_No;
    for (int i= 0; i<chains.size(); i++)
      {
        chains[i].res_num.erase(remove_if(chains[i].res_num.begin(),chains[i].res_num.end(),::isspace),chains[i].res_num.end());
        res_No.push_back(chains[i].res_num);

        if (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
          {
            all_res_No.push_back(res_No);
            res_No.clear();
          }
      }
  return all_res_No;
}

