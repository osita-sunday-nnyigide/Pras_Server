/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
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

