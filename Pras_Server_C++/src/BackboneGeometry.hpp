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

#ifndef BBG_H
#define BBG_H

#include "ReadMMCIF.hpp"

using namespace std;

typedef vector<int> VI;
typedef vector<float> VF;
typedef vector<string> VS;
typedef vector<vector<float>> VVF;
typedef vector<vector<string>> VVS;
typedef vector<vector<vector<float>>> VVVF;

class GetBackboneGeometry: public ReadPDBMMCIF
  {
    public:
    VVF pi_helix();
    VVF alpha_helix();
    VVF _310_helix();
    VVF beta_turn();
    VVF amideDihedral();
    VVF psiDihedral();
    VVF phiDihedral();
    VS  chainID();
    VVS resName();
    VVS dihedralType();
    VVS resNumber();
  };

#endif

