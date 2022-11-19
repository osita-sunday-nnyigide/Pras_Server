/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
********************************************************************************************************************************/

#ifndef CPDBA_H
#define CPDBA_H

#include <iterator>
#include "MissingHeavyAtoms.hpp"

using namespace std;

class CheckMissingAtoms : public FixMissingAtoms
{
	public:
	void repairNotes(string num, string id, vector<string> v);
	void missing_backbone(string num, string id);
	void cite_pras();
	int FixHeavyAtoms();
}; 
#endif