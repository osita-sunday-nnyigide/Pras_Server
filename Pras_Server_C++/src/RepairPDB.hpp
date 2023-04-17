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

#include "AddHydrogen.hpp"

using namespace std;

class OutputStructure : public AddHydrogenAtoms
  {
    public:
    vector<string> atom_site =
                  {"group_PDB",
                 "id",
                 "type_symbol",
                 "label_atom_id",
                 "label_alt_id",
                 "label_comp_id",
                 "label_asym_id",
                 "label_entity_id",
                 "label_seq_id",
                 "pdbx_PDB_ins_code",
                     "Cartn_x",
                     "Cartn_y",
                     "Cartn_z",
                     "occupancy",
                     "B_iso_or_equiv",
                     "pdbx_formal_charge",
                     "auth_seq_id",
                     "auth_comp_id",
                     "auth_asym_id",
                     "auth_atom_id",
                     "pdbx_PDB_model_num"};

    vector<string> symbol = {"C",
                             "H",
                             "N",
                             "O",
                             "S",
                             "#",
                             "loop_"};

    char _element(string atom_name);
    string atomType(string name);
    void writeOutput(char format, string output);
  };