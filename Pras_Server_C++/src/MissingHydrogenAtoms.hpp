/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
********************************************************************************************************************************/

#ifndef MHDA_H
#define MHDA_H

#include "CheckPDBatoms.hpp"
#include "LinearAlgebra.hpp"

using namespace std;

class AddMissingHydrogenAtons: public CheckMissingAtoms
{
	public:
	void arginine_h(Residue &res, Residue &res_p, Residue &ntermini);
	void isoleucine_h(Residue &res, Residue &res_p, Residue &ntermini);
	void leucine_h(Residue &res, Residue &res_p, Residue &ntermini);
	void tryptophan_h(Residue &res, Residue &res_p, Residue &ntermini);
	void phenylalanine_h(Residue &res, Residue &res_p, Residue &ntermini);
	void asparagine_h(Residue &res, Residue &res_p, Residue &ntermini);
	void alanine_h(Residue &res, Residue &res_p, Residue &ntermini);
	void glycine_h(Residue &res, Residue &res_p, Residue &ntermini);
	void valine_h(Residue &res, Residue &res_p, Residue &ntermini);
	void lysine_h(Residue &res, Residue &res_p, Residue &ntermini);
	void aspartate_h(Residue &res, Residue &res_p, Residue &ntermini);
	void glutamine_h(Residue &res, Residue &res_p, Residue &ntermini);
	void glutamate_h(Residue &res, Residue &res_p, Residue &ntermini);
	void histidine_h(Residue &res, Residue &res_p, Residue &ntermini);
	void proline_h(Residue &res, Residue &res_p, Residue &ntermini);
	void metheonine_h(Residue &res, Residue &res_p, Residue &ntermini);
	void threonine_h(Residue &res, Residue &res_p, Residue &ntermini, vector<Residue> &chains);
	void serine_h(Residue &res, Residue &res_p, Residue &ntermini, vector<Residue> &chains);
	void tyrosine_h(Residue &res, Residue &res_p, Residue &ntermini, vector<Residue> &chains);
	void cystine_h(Residue &res, Residue &res_p, vector<Residue> &chains, Residue &ntermini);
	void his_note(string resNo, string chainID);
	void get_hd1(vector<Residue> &HIS, vector<Residue> &chains);
	void prot_his(vector<Residue> &chains);	
};

#endif