/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
********************************************************************************************************************************/

#include "RepairPDB.hpp"

char OutputStructure::_element(string atom_name)
  {

    for (auto i :atom_name)
      {
        if (!isdigit(i))
          {
            return i; 
            break;
          } 
      }
  }

string OutputStructure::atomType(string name)
  {
      if (name.size() == 1)      return " "+name+"  ";
      else if (name.size() == 2) return " "+name+" ";
      else if (name.size() == 3) return " "+name;
      else if (name.size() == 4) return  name;
  }

void OutputStructure:: writeOutput(char format)
  {
    int atom_num; int cnum;
    if (format == 'b'|| format =='t')
      {
        ofstream ofile;
        int atom_num = 1; 
        ofile.open("out.pdb");
        ofile<<setiosflags(ios::fixed)<<setprecision(3);
        ofile<<"# structure genared by PRAS"<<endl;
        for (int i=0;i<chains.size();i++)
          {
            for(int j=0;j<chains[i].res_atom.size();j++)
              {

                ofile<<"ATOM  "<<setw(5)<<atom_num<<" ";
                ofile<<atomType(chains[i].res_atom[j])
                <<" "<<chains[i].res_type<<" "
                <<chains[i].c_id<<setw(4)<<chains[i].res_num
                <<chains[i].res_insert<<"   ";
                for(int k=0;k<3;k++)
                  {
                    ofile<<setw(8)<<chains[i].res_pos[j][k];
                  }
                ofile<<"  1.00  0.00           "
                <<_element(chains[i].res_atom[j])<<"  "<<endl;
                atom_num++;
              }
              if(i==chains.size()-1)
                  {
                    goto TER;
                  }
              else
                  {
                    if(chains[i].c_id!=chains[i+1].c_id)
                            {
                    TER:          ofile<<"TER   "<<setw(5)<<atom_num<<"      "
                                  <<chains[i].res_type<<" "<<chains[i].c_id<<setw(4)
                                  <<chains[i].res_num<<endl;
                            }
                  } 
          }
        if (_hetatm.size()>0)
          {
            for (auto i:_hetatm)ofile<<i<<endl;
          }
        ofile<<"END";
        ofile.close();
      }
    else if (format == 'f')
      {
        cnum = 1; atom_num = 1;
        ofstream ofile;
        ofile.open("out.cif");
        ofile<<setiosflags(ios::fixed)<<setprecision(3);
        ofile<<     "# structure generated by PRAS\n"<<"#"
                    <<"\n"<<"data_xxxx\n" <<"#"<<"\n"
                    <<"_entry.id xxxx\n"<<"#"<<"\n"<<"loop_\n";
        ofile<<"_atom_type.symbol\n";

        for (auto i:symbol)
            {
                ofile<<i<<"\n";
            }
        for (auto i:atom_site) 
          {
            ofile<<"_atom_site."<<i<<"\n";
          }
        for (int i=0;i<chains.size();i++)
          {
            if (chains[i].c_id != chains[fabs(i-1)].c_id)
              {
                cnum+=1;
              }
            for(int j=0;j<chains[i].res_atom.size();j++)
              {

                ofile<<"ATOM  "<<atom_num<<" "<<_element(chains[i].res_atom[j])<<" ";
                ofile<<setw(4)<<atomType(chains[i].res_atom[j])<<" "<<". "<<chains[i].res_type
                <<" "<<chains[i].c_id<<" "<<cnum 
                <<" "<<chains[i].res_num<<" "<<chains[i].res_insert<<" ";
                for(int k=0;k<3;k++)
                  {
                    ofile<<setw(8)<<chains[i].res_pos[j][k];
                  }
                ofile<<" 1.00 0.00 0"<<" "<<chains[i].res_num<<" "<<chains[i].res_type
                <<" "<<chains[i].c_id<<" "<<setw(4)<<atomType(chains[i].res_atom[j])<<" "<<"1 "
                <<endl;
                atom_num++;
              } 
          }
        if (_hetatm.size()>0)
          {
            for (auto i:_hetatm)ofile<<i<<endl;
          }
        ofile<<"#";
        ofile.close();
      } 
  }