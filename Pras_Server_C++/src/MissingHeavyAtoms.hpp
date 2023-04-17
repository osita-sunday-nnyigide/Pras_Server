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

#ifndef MHA_H
#define MHA_H

#include "LinearAlgebra.hpp"
#include "PotentialEnergy.hpp"

using namespace std;

typedef vector<string> RC;
typedef vector<float>  BB;

class FixMissingAtoms : public GetPotentialEnergy
{
	public:
	BB c_termini(BB N, BB CA, BB C, BB O);
	float get_chi(string a1,string a2,string a3,string a4,Residue &res,vector<Residue> &faspr);
	BB addBackbone(BB N, BB CA, BB C, Residue nxtres, float psi);
	void repairSer(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairAla(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres);
	void repairGly(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres);
	void repairAsn(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairGlu(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairAsp(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairGln(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairLys(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairArg(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairTyr(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairTrp(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairPhe(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairHis(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairPro(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairMet(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairIle(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairLeu(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairVal(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairThr(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);
	void repairCys(Residue &res, Residue &res_p, RC missing_atoms, Residue &nxtres, vector<Residue> &faspr);

};

#endif