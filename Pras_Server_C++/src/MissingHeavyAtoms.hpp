/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
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