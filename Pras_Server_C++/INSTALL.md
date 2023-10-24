***************************************************************************************************************************
This is a brief description of how to run a C++ version of PRAS server. 
***************************************************************************************************************************
INSTALLATION

On a linux machine,

g++ -std=c++17 -ffast-math -O3 src/*.cpp -o PRAS

Where src is the folder containing all the source code

-std=c++17 flag tells your compiler to use at least a C++ 17 version, without which this program will not work

-ffast-math -O3 flag tells your compiler to break some rules of floating point operations specified in the IEEE754 standard. 

***************************************************** ABOUT FAST-MATH ******************************************************
With the -ffast-math -O3 flag, this program will be roughly 12 times faster than the python version for any given task,
and without it, the program will be only 3 times faster than the python version. The whole point of fast-math is trading 
off speed with correctness. However, for what this program does, the effect is negligible. E.g., the total number of 
class 6 H-bonds analyzed for the selected proteins in the PRAS paper was reproduced using this C++ version of PRAS server. 
The minor difference is that the orientation differs slightly in some cases owing to the reason stated. For all other 
calculations, the results are the same because it is only in optimization of class 6 H-atoms with a potential energy 
function that this program has to deal with huge floating point arithmetic.
****************************************************************************************************************************

USAGE

./PRAS -i in.pdb -f faspr.pdb -r yes -m yes -l no -h no -p no -s no

-i (input)     enter input PDB to repair and/or analyze

-f (faspr)     enter input PDB as obtained from FASPR or a dummy value to use default chi

-r (rotamer)   enter yes to use high occupancy atoms or a dummy value for low occupancy

-m (mutation)  enter yes to use high occupancy atoms or a dummy value for low occupancy

-l (ligand)    enter yes to keep non-water ligands or a dummy value to ignore them

-h (hydrogen)  enter yes to add hydrogen atoms or a dummy value to do otherwise

-p (protonate) enter yes to protonate 20% of histidines or a dummy value to do otherwise

-s (structure) enter yes to assign secondary structure elements or a dummy value to do otherwise

OPTIONAL

-o (output) enter your preferred output name to keep processed file when processing multiple files, otherwise it will 
            be named out.pdb and removed if you process another PDB
            
-k (keep)   enter yes to keep the log file. All fixed PDBs will be mentioned by name in the log file so that you 
            can identify the ones whose atoms were fixed.  

To understand these OPTIONAL flags, execute the shell.sh to process 713 PDBs (see COMPREHENSIVE TEST WITH 713 PDB IDs). It has been programmed to name all the output files

by PDB 4 letter code

PLOTTING GRAPH

In C++, plotting graph is not as easy as in python. So we use the output of the phi-psi data and the percentage of secondary 
structure written by the C++ version of PRAS and plot the values with a python code named PlotStructure.py that is included 
with this distribution. So, if you want to process a PDB structure and plot the data, the correct usage will be:

./PRAS -i in.pdb -f faspr.pdb -r yes -m yes -l no -h no -p no -s yes;python3 PlotStructure.py jet
 

COMPREHENSIVE TEST WITH 713 PDB IDs

A python document named PDBID.py that contains a list of 713 PDB IDs is included with this distribution. The reason for this
test is to ensure the program works correctly so that if you encounter an error when using a PDB structure from a third-party
then the problem may likely be with the PDB. Of course, there may still be PDB from protein data bank that is damaged beyond
repair that this program cannot fix or with errors whose fixes have not been covered programmatically in PRAS server. 

To do this test, a folder named pdb has been supplied with this distribution and contains the following python documents; 
test.py, Donwnload.py and PDBID.py. cd to this directory and execute the test.py to download 713 PDBs. The speed will depend on
your internet connection speed (implementing ftp protocol in C++ is not as easy as in python so we use python here again).   

When the download is complete, remove the python files and _pycache_ created in pdb folder so that only PDB structures are left.
Now make the shell.sh script supplied with this distribution executable and execute it as follows:

chmod +x shell.sh

time ./shell.sh

This will call PRAS and pass all default arguments required to repair and analyze the 713 proteins. PlotStructure.py is 
commented out in shell.sh. Do not forget to set the -s flag to yes if you intend to plot the secondary structure.
In windows subsystem for linux, one must have X-server installed and running before using graphics. 
Detailed information on X-server is given in the README supplied with the python distribution of PRAS server.

