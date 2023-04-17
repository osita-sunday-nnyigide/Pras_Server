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

#include "ReadMaster.hpp"

Atom::Atom(string atom_num,
      string alt_conform,
      string atom_type,
      string res_type,
      string chain_id,
      string res_num,
      string model,
      string res_insert,
      float occupancy,
      float x,
      float y,
      float z)
  {
  a_num         = atom_num;
  a_conform     = alt_conform;
  a_type        = atom_type;
  r_type        = res_type;
  c_id          = chain_id;
  r_num         = res_num;
  m_num         = model;
  r_insert      = res_insert;
  a_occupancy   = occupancy;
  x_coord       = x;
  y_coord       = y;
  z_coord       = z;
  pos           = {x,y,z};
  }
Residue::Residue(
    string residue_type,
    string insert,
    string residue_num,
    string chain_id,
    float occupancy,
    vector<vector<float>> residue_pos,
    vector<string> atom_type,
    vector<float> r_ocpancy)
{
res_type      = residue_type;
res_pos       = residue_pos;
res_atom      = atom_type;
res_ocpancy   = occupancy;
res_num       = residue_num;
c_id          = chain_id;
res_insert    = insert;
_ocpancy      = r_ocpancy;
}
void Residue::appendAtom(Atom &atom, string rotamer)
{
  if (find(res_atom.begin(),res_atom.end(), atom.a_type) == res_atom.end())
    {
    res_atom.push_back(atom.a_type);
    res_pos.push_back(atom.pos);
    _ocpancy.push_back(atom.a_occupancy);
    }
  else if ((find(res_atom.begin(),res_atom.end(), atom.a_type) != res_atom.end()) &&  rotamer != "yes")
    {
      if (_ocpancy.back() < atom.a_occupancy)
        {
          res_atom.back() = atom.a_type;
          res_pos.back()  = atom.pos;
        }
    }
  else if ((find(res_atom.begin(),res_atom.end(), atom.a_type) != res_atom.end()) && rotamer == "yes")
    {
      if (_ocpancy.back() < atom.a_occupancy)
        {
          res_atom.back() = atom.a_type;
          res_pos.back()  = atom.pos;
        }
    }
}

void Chain::appendResidue(Residue &res, float ocpancy, string mutation)
  {
    if (_chain[0].empty())
      {
        _chain[0].push_back(res);

      }

    if (!_chain[0].empty())
      {

        if ((_chain[0].back().res_num != res.res_num) || ((_chain[0].back().res_insert != res.res_insert)))
          {
          _chain[0].push_back(res);
          }

          else if ((_chain[0].back().res_num == res.res_num) && mutation != "yes")
            {
            if ((float) _chain[0].back().res_ocpancy > (float) res.res_ocpancy)
              {
               _chain[0].back() = res;
              }
            }

          else if ((_chain[0].back().res_num == res.res_num) && mutation == "yes")
            {
            if ((float) _chain[0].back().res_ocpancy < (float) res.res_ocpancy)
              {
               _chain[0].back() = res;
              }
            }
      }
  }
void Chain::insertAtom(Atom &atom, string rotamer)
  {
  _chain[0].back().appendAtom(atom, rotamer);
  }
Atom read_atom(string &PDBline)
{
using namespace std;
Atom atom;
atom.a_num         = PDBline.substr(6,5);
atom.a_conform     = PDBline.substr(16,1);
atom.a_type        = PDBline.substr(12,4);
atom.r_type        = PDBline.substr(17,3);
atom.c_id          = PDBline.substr(21,1);
atom.r_num         = PDBline.substr(22,4);
atom.r_insert      = PDBline.substr(26,1);
atom.x_coord       = stof(PDBline.substr(30,8));
atom.y_coord       = stof(PDBline.substr(38,8));
atom.z_coord       = stof(PDBline.substr(46,8));
atom.a_occupancy   = stof(PDBline.substr(55,5));
atom.pos           = {atom.x_coord,atom.y_coord,atom.z_coord};
atom.a_type.erase(remove_if(atom.a_type.begin(), atom.a_type.end(), ::isspace), atom.a_type.end());

if (atom.r_insert == " ")
  {
    atom.r_insert = "";
  }
return atom;
}

void Chains::getChains(string PDB, string PDB2, string rotamer, string mutation, string ligand)
{

  string PDBline;

  if (PDB == "noPDB") PDB = PDB2;
  ifstream MyReadFile(PDB);

  string res_num = "-1"; string res_insert = " "; string res_type = " ";
  Chain chain;

  while (getline(MyReadFile, PDBline)) {

    if (PDBline.substr(0,4) == "ATOM")
    {
        string atm = PDBline.substr(12,4);

        atm.erase(remove_if(atm.begin(), atm.end(), ::isspace), atm.end());
        if (find(heay_atoms.begin(), heay_atoms.end(),atm) != heay_atoms.end())
        {
          Atom atoms;
          atoms     = read_atom(PDBline);
          atoms.rot = rotamer;
          if ((res_num!=atoms.r_num) || (res_insert!=atoms.r_insert) || (res_num==atoms.r_num && res_type!=atoms.r_type))
            {
              Residue residue = Residue(atoms.r_type,atoms.r_insert,atoms.r_num,atoms.c_id,atoms.a_occupancy);
              residue.mutation = mutation;
              chain.appendResidue(residue,atoms.a_occupancy,residue.mutation);
              res_num = atoms.r_num;
              res_insert = atoms.r_insert;
              res_type = atoms.r_type;
            }
          if (res_type==atoms.r_type && res_type == chain._chain[0].back().res_type)
            {
              chain.insertAtom(atoms, rotamer);
            }
        }
      }
    if (ligand == "yes" && PDBline.substr(0,6) == "HETATM" && PDBline.substr(17,3) != "HOH")
      {
        _hetatm.push_back(PDBline);
      }
    if (PDBline.substr(0,3) == "END" || MyReadFile.peek() == EOF)
      {
        for (int i = 0; i< chain._chain[0].size(); i++)
          {
            if (PDB == PDB2)
              {
                faspr.push_back(chain._chain[0][i]);
              }
            else
              {
                  chains.push_back(chain._chain[0][i]);
              }
          }
    chain._chain[0].clear();
    MyReadFile.close();
      }
  }
}