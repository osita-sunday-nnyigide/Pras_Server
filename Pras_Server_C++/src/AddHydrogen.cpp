/*******************************************************************************************************************************
Copyright (c) 2022 Osita Sunday Nnyigide (osita@protein-science.com)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy,
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
        else if (chains[i].res_type == "HIS" || chains[i].res_type == "HSD" ||chains[i].res_type == "HSP" || chains[i].res_type == "HSE" || chains[i].res_type == "HIE")
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