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
