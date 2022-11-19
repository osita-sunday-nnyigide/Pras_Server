/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
********************************************************************************************************************************/

#ifndef RCIF_H
#define RCIF_H

#include "ReadMaster.hpp" 

using namespace std;

typedef vector<string> VS;

class ReadPDBMMCIF : public Chains, public Chain, public Residue, public Atom
	{
		public:
		void getMMCIFChains(string PDB, string rotamer, string mutation, string ligand);
	};

Atom read_atom(VS &PDBline);

#endif