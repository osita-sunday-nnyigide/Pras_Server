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