/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
********************************************************************************************************************************/

#include "AddHydrogen.hpp"

void  AddHydrogenAtoms::addH(string hprot)
  {
    cout<<"\n**Now fixing all hydrogen atoms...\n"<<endl;
    vector<Residue> copy = chains;
    for (int i = 0; i < chains.size(); i++)
      {
        if (chains[i].res_type == "ARG")
          {
            arginine_h(chains[i], chains[fabs(i-1)], chains.front());
          }
        else if (chains[i].res_type == "ILE")
          {
            isoleucine_h(chains[i], chains[fabs(i-1)], chains.front());
          }
        else if (chains[i].res_type == "LEU")
          {
            leucine_h(chains[i], chains[fabs(i-1)], chains.front());
          }
        else if (chains[i].res_type == "TRP")
          {
            tryptophan_h(chains[i], chains[fabs(i-1)], chains.front());
          }
        else if (chains[i].res_type == "PHE")
          {
            phenylalanine_h(chains[i], chains[fabs(i-1)], chains.front());
          }
        else if (chains[i].res_type == "ASN")
          {
            asparagine_h(chains[i], chains[fabs(i-1)], chains.front());
          }
        else if (chains[i].res_type == "ALA")
          {
            alanine_h(chains[i], chains[fabs(i-1)], chains.front());
          }
        else if (chains[i].res_type == "GLY")
          {
            glycine_h(chains[i], chains[fabs(i-1)], chains.front());
          }
        else if (chains[i].res_type == "VAL")
          {
            valine_h(chains[i], chains[fabs(i-1)], chains.front());
          }
        else if (chains[i].res_type == "LYS")
          {
            lysine_h(chains[i], chains[fabs(i-1)], chains.front());
          }
        else if (chains[i].res_type == "ASP")
          {
            aspartate_h(chains[i], chains[fabs(i-1)], chains.front());
          }
        else if (chains[i].res_type == "GLN")
          {
            glutamine_h(chains[i], chains[fabs(i-1)], chains.front());
          }
        else if (chains[i].res_type == "HIS")
          {
            histidine_h(chains[i], chains[fabs(i-1)], chains.front());
          }
        else if (chains[i].res_type == "PRO")
          {
            proline_h(chains[i], chains[fabs(i-1)], chains.front());
          }
        else if (chains[i].res_type == "MET")
          {
            metheonine_h(chains[i], chains[fabs(i-1)], chains.front());
          }
        else if (chains[i].res_type == "GLU")
          {
            glutamate_h(chains[i], chains[fabs(i-1)], chains.front());
          }
        else if (chains[i].res_type == "THR")
          {
            threonine_h(chains[i], chains[fabs(i-1)], chains.front(), copy);
          }
        else if (chains[i].res_type == "SER")
          {
            serine_h(chains[i], chains[fabs(i-1)], chains.front(), copy);
          }
        else if (chains[i].res_type == "TYR")
          {
            tyrosine_h(chains[i], chains[fabs(i-1)], chains.front(), copy);
          }
        else if (chains[i].res_type == "CYS")
          {
            cystine_h(chains[i], chains[fabs(i-1)], copy, chains.front());
          }
      }
      if (hprot == "yes")
      {
        prot_his(chains);
      }
  }