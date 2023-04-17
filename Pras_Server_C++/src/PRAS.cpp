

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
#include "RepairPDB.hpp"
#include "SecondaryStructure.hpp"

using namespace std;

int main(int argc, char** argv)

{
    cout<<string(60,'*')<<"\n";
    cout<<"            PRAS 1.0.11"<<endl;
    cout<<"Protein Repair and Analysis Server"<<endl;
    cout<<endl;
    cout<<"Copyright (c) 2022 Osita Sunday Nnyigide"<<endl;
    cout<<"Complex Fluid Laboratory"<<endl;
    cout<<"School of Chemical and Biomolecular Engineering"<<endl;
    cout<<"Pusan National University"<<endl;
    cout<<"Email:osita@protein-science.com"<<endl;
    cout<<string(60,'*')<<"\n";
	clock_t start_time,end_time;
    float duration;
    start_time = clock();
	string PDB, faspr, rotamer, mutation, ligand, addh, prot, ss, output, log; int i;
	OutputStructure PRAS;AnalyzeProtein SSA;
	for(i=1;i<argc-1;i++)
		{
	   if(argv[i][0]=='-')
		   {
		     if(argv[i][1]=='i')
			     {
			       i++;
			       PDB=argv[i];
			     }
		     else if(argv[i][1]=='f')
			     {
			       i++;
			       faspr=argv[i];
			     }
		     else if(argv[i][1]=='r')
			     {
			       i++;
			       rotamer=argv[i];
			     }
		     else if(argv[i][1]=='m')
			     {
			       i++;
			       mutation=argv[i];
			     }
		     else if(argv[i][1]=='l')
			     {
			       i++;
			       ligand=argv[i];
			     }
		     else if(argv[i][1]=='h')
			     {
			       i++;
			       addh=argv[i];
			     }
		     else if(argv[i][1]=='p')
			     {
			       i++;
			       prot=argv[i];
			     }
		     else if(argv[i][1]=='s')
			     {
			       i++;
			       ss=argv[i];
			     }
		     else if(argv[i][1]=='o')
			     {
			       i++;
			       output=argv[i];
			     }
		     else if(argv[i][1]=='k')
			     {
			       i++;
			       log=argv[i];
			     }
		   }
		}
	if (log.empty()) remove("log.txt");
	char ftype = PDB.back();
	if ((PDB,faspr,rotamer,mutation,ligand,addh,prot,ss).empty()) goto ERR;
	if (ftype=='b'||ftype=='t')PRAS.getChains(PDB,"nofaspr",rotamer,mutation,ligand);
	else if (ftype=='f')         PRAS.getMMCIFChains(PDB, rotamer, mutation, ligand);
	PRAS.getChains("noPDB",faspr,rotamer,mutation,ligand);
	PRAS.FixHeavyAtoms(PDB);
	if (addh=="yes")PRAS.addH(prot);
	if (output.empty()) output = "no";
	PRAS.writeOutput(ftype,output);
	if (ss=="yes") SSA.SecondaryStructure(ftype,output);
    end_time = clock();
    duration = (float)(end_time-start_time)/CLOCKS_PER_SEC;
    cout<<"\nThis program took "<<duration<<" seconds\n"<<endl;
    return 0;

	ERR:  cout<<"Error! PRAS requires 8 arguments"<<endl;
		  cout<<"Usage: ./PRAS -i in.pdb -f faspr.pdb -r yes -m yes -l no -h no -p no -s no"<<endl;
		  cout<<"-i (input)     enter input PDB to repair and/or analyze"<<endl;
		  cout<<"-f (faspr)     enter input PDB as obtained from FASPR or a dummy value to use default chi"<<endl;
		  cout<<"-r (rotamer)   enter 'yes' to use high occupancy atoms or a dummy value for low occupancy"<<endl;
		  cout<<"-m (mutation)  enter 'yes' to use high occupancy atoms or a dummy value for low occupancy"<<endl;
		  cout<<"-l (ligand)    enter 'yes' to keep non-water ligands or a dummy value to ignore them"<<endl;
		  cout<<"-h (hydrogen)  enter 'yes' to add hydrogen atoms or a dummy value to do otherwise"<<endl;
		  cout<<"-p (protonate) enter 'yes' to protonate 20%" " of histidines or a dummy value to do otherwise"<<endl;
		  cout<<"-s (structure) enter 'yes' to assign secondary structure elements or a dummy value to do otherwise"<<endl;
		  exit(1);
}
