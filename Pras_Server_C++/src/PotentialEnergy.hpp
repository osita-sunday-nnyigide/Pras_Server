/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
********************************************************************************************************************************/

#ifndef PE_H
#define PE_H

#include "ReadMMCIF.hpp"
#include "LinearAlgebra.hpp"

using namespace std;

typedef vector<float> V1;
typedef vector<vector<float> > V2;

class GetPotentialEnergy:public ReadPDBMMCIF
{
	public:
	V2 rotate(V1 pos1, V1 pos2, V1 pos3, V1 pos4, float BL, float di);
	V1 optmizeH(V1 pos1, V1 pos2, V1 pos3, V1 pos4, float BL, float di, float H_pch,float H_sig, float H_eps,Residue &res,vector<Residue> &chains);
};


#endif
