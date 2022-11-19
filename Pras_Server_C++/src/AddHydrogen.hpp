/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
********************************************************************************************************************************/

#ifndef ADDH_H
#define ADDH_H

#include "MissingHydrogenAtoms.hpp"

class AddHydrogenAtoms: public AddMissingHydrogenAtons
{
  public: 
  void  addH(string hprot);
};

#endif