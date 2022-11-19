/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
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
    void writeOutput(char format);
  };