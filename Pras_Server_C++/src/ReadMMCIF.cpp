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

#include "ReadMMCIF.hpp"

Atom read_atom(VS &PDBline)
  {
  Atom atom;
  atom.a_num         = PDBline[1];
  atom.a_type        = PDBline[3];
  atom.r_type        = PDBline[5];
  atom.c_id          = PDBline[6];
  atom.r_num         = PDBline[8];
  atom.r_insert      = PDBline[9];
  atom.m_num         = PDBline[20];
  atom.x_coord       = stof(PDBline[10]);
  atom.y_coord       = stof(PDBline[11]);
  atom.z_coord       = stof(PDBline[12]);
  atom.a_occupancy   = stof(PDBline[13]);
  atom.pos           = {atom.x_coord,atom.y_coord,atom.z_coord};

  if (atom.r_insert == " ")
    {
      atom.r_insert = "";
    }
  return atom;
  }

void ReadPDBMMCIF::getMMCIFChains(string PDB, string rotamer, string mutation, string ligand)
	{
		chains.clear();
	    string PDBline;
	    ifstream MyReadFile(PDB);
	    string res_num = "-1"; string res_insert = " "; string res_type = " ";
	    string model_no = " "; int i; VS model; Chain chain;
	    while (getline(MyReadFile, PDBline))
	    {
	      if (PDBline.substr(0,4) == "ATOM")
		      {
		      	size_t idx = 0;
		      	string sep = " ";
		      	VS atomline = {};
		      	while (( idx = PDBline.find (sep)) != string::npos)
			      	{
						string field = PDBline.substr(0, PDBline.find(sep));
						if(field.find_first_not_of(" ") != string::npos)
							{
								atomline.push_back(field);
							}
				    	PDBline.erase(0, idx + sep.length());
					}
					if (find(heay_atoms.begin(),heay_atoms.end(),atomline[3]) != heay_atoms.end())
						{
	            			Atom atoms;
	            			atoms     = read_atom(atomline);
	            			atoms.rot = rotamer;
	            			model.push_back(atoms.m_num);
							atomline.clear();
				            if ((res_num!=atoms.r_num) || (res_insert!=atoms.r_insert) || (res_num==atoms.r_num && res_type!=atoms.r_type))
				              {
				                Residue residue = Residue(atoms.r_type,atoms.r_insert,atoms.r_num,atoms.c_id,atoms.a_occupancy);
				                residue.mutation = mutation;
				                if (atoms.m_num == model[0]) chain.appendResidue(residue,atoms.a_occupancy,residue.mutation);
				                else goto LIG;
				                res_num    = atoms.r_num;
				                res_insert = atoms.r_insert;
				                res_type   = atoms.r_type;
				                model_no   = atoms.m_num;
				              }
				            if (res_type==atoms.r_type && res_type==chain._chain[0].back().res_type)
				              {
				                chain.insertAtom(atoms, rotamer);
				              }
					    }
					atomline.clear();
		      }
    LIG:  	if (ligand == "yes" && PDBline.substr(0,6) == "HETATM" && PDBline.find("HOH") == string::npos)
	        {
	          _hetatm.push_back(PDBline);
	        }
      	if (MyReadFile.peek() == EOF)
	        {
	          for (i = 0; i< chain._chain[0].size(); i++)
	            {
	              chains.push_back(chain._chain[0][i]);
	            }
	      	chain._chain[0].clear();
	      	model.clear();
	      	MyReadFile.close();
	        }
		  }
	}