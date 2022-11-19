/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
********************************************************************************************************************************/

#ifndef RMPDB_H
#define RMPDB_H

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

using namespace std;

class Atom
{
 
 public:
  string a_num;
  string a_conform;
  string a_type;
  string r_type;
  string c_id;
  string r_num;  
  string r_insert;
  string m_num;
  float a_occupancy;
  string rot;
  float x_coord;
  float y_coord;
  float z_coord; 
  vector<float> pos;

  Atom(string atom_num = " ", 
      string alt_conform = " ", 
      string atom_type = " ",
      string res_type = " ",
      string chain_id = " ",
      string res_num = " ",
      string model   = " ",
      string res_insert = " ",
      float occupancy = 0.0,
      float x = 0.,
      float y = 0.,
      float z = 0.);
};

class Residue
{
 
 public:
  float res_ocpancy;
  string res_type;
  string c_id;
  string res_num;  
  string res_insert;
  string mutation;
  vector<string> res_atom;
  vector<vector<float>> res_pos;
  vector<string> _atom;
  vector<float> _ocpancy = {};

  Residue(
      string residue_type = " ", 
      string insert = " ",
      string residue_num = " ",
      string chain_id = " ",
      float occupancy = 0.0,
      vector<vector<float>> residue_pos = {},
      vector<string> atom_type = {},
      vector<float> r_ocpancy = {});

  void appendAtom(Atom &atom, string rotamer);

};

class Chain
{
  public:
  vector<vector<Residue>> _chain = {{}};
  void appendResidue(Residue &res, float ocpancy, string mutation);

  void insertAtom(Atom &atom, string rotamer);

};

Atom read_atom(string &PDBline);

class Chains
{
  public:
  vector<string> heay_atoms  = {

    "N",  "CA",  "C",   "O",   "CB",  "CG", 
    "OD1","OD2", "OXT", "CD",  "CE",  "NZ", 
    "OG", "CG1", "CG2", "CD1", "CD2", "OG1", 
    "CE2","CE3", "NE1", "CZ2", "CZ3", "CH2", 
    "NE", "CZ",  "NH1", "NH2", "ND2", "NE2", 
    "OE1","SG",  "OE2", "ND1", "CE1", "SD", 
    "OH"
                                          };
  vector<string> _hetatm = {};                                          
  vector<Residue> chains = {};
  vector<Residue> faspr  = {};
  void getChains(string PDB, string PDB2, string rotamer, string mutation, string ligand);

};
#endif