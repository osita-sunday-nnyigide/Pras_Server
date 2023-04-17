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

#include "CheckPDBatoms.hpp"

void CheckMissingAtoms::repairNotes(string num, string id, vector<string> v)
    {
    	num.erase(remove_if(num.begin(),num.end(),::isspace),num.end());
        ofstream file("log.txt", ios_base::app);
        ostream_iterator<string> cout_it(file, " ");
        file<<"Residue No."<< num <<" in chain "
             <<id<< " has missing atom(s) [";
             copy(v.begin(), v.end(), cout_it);
        file<<"]. We have fixed it.\n\n";
    }

void CheckMissingAtoms::missing_backbone(string num, string id)
	{
		cout<<"Residue No."<< num <<" in chain "
			 <<id<< " has missing backbone atoms";
		cout<<"\nprogram has terminated abnormally."
			  " For a .cif file,\n"
			  "label_seq_id is reported\n\n";
	}

void CheckMissingAtoms::cite_pras(string fname)
	{
		ofstream file("log.txt", ios_base::app);
		file<<"When using files generated by this program in a publication,"
		        " please cite this program as O.S. Nnyigide, T.O. Nnyigide,"
		        " S.G. Lee, K. Hyun. Protein Repair and Analysis Server: A Web\n"
		        "Server to Repair PDB Structures, Add Missing Heavy"
		        " Atoms and Hydrogen Atoms, and Assign Secondary Structures by"
		        " Amide Interactions. J. Chem. Inf. Model., 2022, 62, 4232–4246.\n";
		file<<string(183,'*')<<"\n";
		file<<"PRAS 1.0.11. This is a PRAS-generated log file."
		        " For your information, all missing or fixed atoms"
		        " and other relevant information concerning the repair"
		        " are appended below\n";
		file<<string(183,'*')<<"\n";
		file<<"This is the log file for "<<fname<<"\n\n";
	}

int CheckMissingAtoms::FixHeavyAtoms(string fname)
	{
	CheckMissingAtoms::cite_pras(fname);
	if (faspr.empty())
		{
			ofstream file("log.txt", ios_base::app);
			file<<"No FASPR output PDB supplied. Default chi will be used\n\n";
		}
	cout<<"\n**Now checking and fixing missing heavy atoms..."<<endl;
    for (int i= 0; i<chains.size(); i++)
            {
			    if (chains[i].res_type == "ASP")
	              {
	              	RC all_atoms = { }; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
			             all_atoms = {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2","OXT"};
			        else all_atoms = {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"};

		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(), chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
				        }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairAsp(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
			    else if (chains[i].res_type == "ALA")
	              {
	              	RC all_atoms = { }; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	 	 all_atoms = {"N", "CA", "C", "O", "CB","OXT"};
	              	else all_atoms = {"N", "CA", "C", "O", "CB"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairAla(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()]);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
			    else if (chains[i].res_type == "ARG")
	              {
	              	RC all_atoms = { }; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
			             all_atoms = {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2","OXT"};
			        else all_atoms = {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairArg(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
                else if (chains[i].res_type == "ASN")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	     all_atoms = {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2", "OXT"};
	              	else all_atoms = {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairAsn(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
                else if (chains[i].res_type == "CYS")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	     all_atoms = {"N", "CA", "C", "O", "CB", "SG", "OXT"};
	              	else all_atoms = {"N", "CA", "C", "O", "CB", "SG"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairCys(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
                else if (chains[i].res_type == "GLU")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	     all_atoms = {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2", "OXT"};
	              	else all_atoms = {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairGlu(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
                else if (chains[i].res_type == "GLN")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	     all_atoms = {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2", "OXT"};
	              	else all_atoms = {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairGln(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
                else if (chains[i].res_type == "GLY")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	     all_atoms = {"N", "CA", "C", "O", "OXT"};
	              	else all_atoms = {"N", "CA", "C", "O"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairGly(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()]);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
                else if (chains[i].res_type == "HIS" || chains[i].res_type == "HSD" ||chains[i].res_type == "HSP" || chains[i].res_type == "HSE" || chains[i].res_type == "HIE")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              		 all_atoms  = {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2", "OXT"};
	              	else all_atoms  = {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairHis(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
                else if (chains[i].res_type == "ILE")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	     all_atoms = {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1", "OXT"};
	              	else all_atoms = {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairIle(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
                else if (chains[i].res_type == "LYS")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	     all_atoms  = {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ", "OXT"};
	              	else all_atoms  = {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"};

		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairLys(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
                else if (chains[i].res_type == "LEU")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	     all_atoms  = {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "OXT"};
	              	else all_atoms  = {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num,chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairLeu(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
                else if (chains[i].res_type == "MET")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	     all_atoms  = {"N", "CA", "C", "O", "CB", "CG", "SD", "CE", "OXT"};
	              	else all_atoms  = {"N", "CA", "C", "O", "CB", "CG", "SD", "CE"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairMet(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
                else if (chains[i].res_type == "PRO")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	     all_atoms  = {"N", "CA", "C", "O", "CB", "CG", "CD", "OXT"};
	              	else all_atoms  = {"N", "CA", "C", "O", "CB", "CG", "CD"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairPro(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
                else if (chains[i].res_type == "PHE")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	     all_atoms  = {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OXT"};
	              	else all_atoms  = {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairPhe(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
                else if (chains[i].res_type == "SER")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	     all_atoms  = {"N", "CA", "C", "O", "CB", "OG", "OXT"};
	              	else all_atoms  = {"N", "CA", "C", "O", "CB", "OG"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairSer(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
				else if (chains[i].res_type == "THR")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	     all_atoms  = {"N", "CA", "C", "O", "CB", "OG1", "CG2", "OXT"};
	              	else all_atoms  = {"N", "CA", "C", "O", "CB", "OG1", "CG2"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairThr(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
				else if (chains[i].res_type == "TYR")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	     all_atoms  = {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "OXT"};
	              	else all_atoms  = {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairTyr(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
				else if (chains[i].res_type == "TRP")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	     all_atoms  = {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2", "OXT"};
	              	else all_atoms  = {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairTrp(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
				else if (chains[i].res_type == "VAL")
	              {
	              	RC all_atoms = {}; RC missing = {};
	              	if   (chains[i].c_id != chains[(i+1)%chains.size()].c_id || (chains[i].res_num == chains.back().res_num && chains[i].c_id == chains.back().c_id))
	              	     all_atoms  = {"N", "CA", "C", "O", "CB", "CG1", "CG2", "OXT"};
	              	else all_atoms  = {"N", "CA", "C", "O", "CB", "CG1", "CG2"};
		            for (string k: all_atoms)
			            {
			            	if ((find(chains[i].res_atom.begin(),chains[i].res_atom.end(), k) == chains[i].res_atom.end()))
				            	{
				            		missing.push_back(k);
				            	}
			            }
			        if (!missing.empty())
		            	{
		            		for (auto l: missing)
			            		{
			            			if ( l == "N" || l == "CA" || l == "C")
				            			{
				            				missing_backbone(chains[i].res_num, chains[i].c_id);
				            				exit(1);
				            			}
			            		}
			            	repairVal(chains[i], chains[fabs(i-1)], missing, chains[(i+1)%chains.size()], faspr);
			            	repairNotes(chains[i].res_num, chains[i].c_id, missing);
		            	}
	              }
            }
	}
