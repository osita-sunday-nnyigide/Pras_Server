Protein Repair and Analysis (PRAS) Server is a tool to repair PDB or mmCIF structures, add missing heavy 

atoms and hydrogen atoms and assign secondary structures by amide-amide interactions of the backbone

You can use the server online at https://www.protein-science.com/ or on your local machine following the 

instructions given below

The program consists of:
```python
Tools to automatically download PDB or mmCIF structures

Tools to generate different rotamers if present in a PDB or mmCIF structure file

Tools to replace missing heavy atoms or entire side-chain

Tools to add hydrogen atoms, including optimization of those with rotational freedom

Tools to assign secondary structure elements by amide-amide interactions of the backbone

Tools to draw 4 Ramachandran types (i.e., general, glycine, proline and pre-proline)
```

## PRAS installation

`pip install Pras-Server==1.2.1`

## PRAS usage

A comprehensive test with 82 or 494 PDB or mmCIF structures may be performed to ensure the program works correctly.

Since PRAS can download these structures automatically, 2 lists are provided which contain all the PDB IDs.

To do this test, make a folder in any directory of your choosing where these structures would be downloaded or

stored. Whilst in the folder, create a python document (e.g. pras_test.py). Enter the following code and execute:

```python
#computational time depends on
#user's internet download speed
import time
from Pras_Server.RunType import InitRunType
from Pras_Server.PDBID import  _82_pdbs, _494_pdbs

startTime = time.time()

fixing = InitRunType( rotamer="", 
					mutation="", 
					pdb_faspr="", 
					keep_ligand="", 
					chain_no="",
					addh = False,
					ss=False,
					raman=False,
					ofname=False,
					pdbid=_82_pdbs,
		     			his_p=False
					)

# fixed PDB is saved with name=PDB ID+_out.pdb
# addh, use default or set to True to add hydrogen
# pdbid, must be a list of PDB IDs or set to False
# ss=secondary structure, set to True if you need it
# raman=ramachandran plot, set to True if you need it
# rotamer, use default or "yes" to use low occupancy conformers
# mutation, use default or "yes" to use low occupancy conformers
# keep_ligand, use default or "yes" to keep non-water ligands
# chain_no, use default or string number (e.g., "1") to process a specific chain
# ofname, output file name, use default. To use your own name read documentaion below
# his_p, protonate 20% of HIS. A fraction of HIS may be charged at pH of 7. Use default
# Execute print(InitRunType.ProcessOther().__doc__) or see RunType.py doc instructions
# pdb_faspr, use default or enter a pdb output from faspr, reason is stated in pras paper

fixing.ProcessWithoutDefaultUsingPDBID()
#You can print(InitRunType.ProcessWithoutDefaultUsingPDBID.__doc__) 
#to learn more about the member function
print ('The program took {0} second !'.format(time.time() - startTime))
```

The above will download and analyze 82 or 494 protein structures. One can modify the code to process either of the lists.

If you set ss=True and rama=True, you may encounter failed to allocate bitmap error (depending on your computer's memory).

Note that the above test will work in both windows and linux systems. For linux you just need to do python3 pras_test.py

For users that do not want to do the above comprehensive test, instructions are given below to process a

single protein structure in a windows or linux environment.

For winddows, download a protein (e.g toxin II, 1aho.pdb, or let PRAS download it automatically)

Open a Python document (e.g., example.py) and copy 1aho.pdb to the directory where example.py is located

Copy and paste the following code in example.py and execute (you must have permission to write in this directory)

PRAS will replace all missing atoms of toxin II, assign secondary structure and draw 4 types of ramachandran plots


```python
import time
from Pras_Server.RunType import InitRunType

startTime = time.time()

fixing = InitRunType(rotamer="", 
					mutation="", 
					pdb_faspr="", 
					keep_ligand="", 
					chain_no="",
					addh = False,
					ss=False,
					raman=False,
					ofname=False,
					pdbid=False,
		     			his_p=False
					)

fixing.fname = '1aho.pdb'

fixing.ProcessWithDefault()
#You can print(InitRunType.ProcessWithDefault.__doc__) 
#to learn more about the member function
print ('The program took {0} second !'.format(time.time() - startTime))
```

For linux, do the same as above except for the example code where you should copy the code below (the annotation given above applies below):

```python
import sys

if len(sys.argv) == 9:
	import time
	from Pras_Server.RunType import InitRunType
	startTime = time.time()

	rotamer=sys.argv[2] 
	mutation=sys.argv[3]
	pdb_faspr=sys.argv[4]
	keep_ligand=sys.argv[5]
	chain_no=sys.argv[6] 
	ofname=sys.argv[7]
	his_p=sys.argv[8]
	fixing=InitRunType(
					rotamer, 
					mutation, 
					pdb_faspr, 
					keep_ligand, 
					chain_no,
					addh=False,
					ss=False,
					raman=False,
					ofname=False,
					pdbid=False,
					his_p=False
					)
	fixing.fname=sys.argv[1]
	fixing.ProcessWithDefault()
	print ('The program took {0} second !'.format(time.time() - startTime))

else:
	print("PRAS takes 9 compulsory arguments." 
		  " Execute the code below on your shell prompt to learn more\n")
	print("printf \"from Pras_Server.PRAS import" 
		  " repairPDB\\nprint(repairPDB.__doc__)\\n\\n()\" | python3\n")
```

Then, cd to the directory where example.py is located and enter the following argument on your shell prompt:

`python3 example.py 1aho.pdb "" "" "" "" "" "" ""`


## WINDOWS SUBSYSTEM FOR LINUX (WSL)

In order to run Linux GUI applications e.g., the plots for

secondary structure assignment and 4 Ramachandran types

using Windows Subsystem for Linux (WSL), you must install X server for Windows.

Thus, you need to:

Install X server for Windows

Configure bash to tell GUIs to use the local X server

For X server, install VcXsrv which is open source by downloading from https://sourceforge.net/projects/vcxsrv/

Configure bash to use the local X server. In bash run:

`echo "export DISPLAY=localhost:0.0" >> ~/.bashrc`

To have the configuration changes take effect, restart bash, or run:

`. ~/.bashrc`

Then open VcXsrv from your taskbar (you should send the icon to taskbar for easy access).
Note that VcXsrv must be open/running each time you use plotting tools in WSL

PRAS has been peer reviewed and published. Please cite PRAS as:

O.S. Nnyigide, T.O. Nnyigide, S.G. Lee, K. Hyun. Protein Repair and Analysis Server: A Web Server to Repair PDB Structures, Add Missing Heavy Atoms and Hydrogen Atoms, and Assign Secondary Structures by Amide Interactions.
J. Chem. Inf. Model., 2022, 62, 4232–4246.

## License
[MIT](https://choosealicense.com/licenses/mit/)

