Pras Server is a library to repair PDB or mmCIF structures, add missing heavy atoms and hydrogen atoms and assign secondary structure by amide-amide interactions of the backbone

You can use the server online at https://www.protein-science.com/ or on your local machine following the instructions given below

The library consists of:
```python
Tools to automatically download PDB or mmCIF structures

Tools to generate different rotamers if present in a PDB or mmCIF structure file

Tools to replace missing heavy atoms or entire side-chain

Tools to add hydrogen atoms, including optimization of those with rotational freedom

Tools to assign secondary structure elements by amide-amide interactions of the backbone

Tools to draw 4 Ramachandran types (i.e., general, glycine, proline and pre-proline)
```

## PRAS installation

`pip install Pras-Server==1.0.7`

## PRAS usage

A comprehensive test with 82 or 494 PDB or mmCIF structures may be performed to ensure the program works correctly.

Since PRAS can download these structures automatically, 2 lists are provided which contain all the PDB IDs.

To do this test, make a folder in any directory of your choosing where these structures would be downloaded or 

stored. Whilst in the folder, create a python document (e.g. testpras.py). Enter the following code and execute:

```python
#computational time depends on
#user's internet download speed
from Pras_Server import test
```

The above will download and analyze 82 or 494 protein structures. One can modify test.py to process either of the lists

or comment the SS assignment and Ramachandran plots to save memory (if you  encounter failed to allocate bitmap)

For users that do not want to do the above comprehensive test, instructions are given below to process a 

single protein structure in a windows or linux environment.

For winddows, download a protein (e.g toxin II, 1aho.pdb, or let PRAS download it automatically)

Open a Python document (e.g., example.py) and copy 1aho.pdb to the directory where example.py is located

Copy and paste the following code in example.py and execute (you must have permission to write in this directory)

PRAS will replace all missing atoms of toxin II, assign secondary structure and draw 4 types of ramachandran plots


```python
#computational time depends on what you do
import time
startTime = time.time()

#this function replaces missing heavy atoms and adds H-atoms. 
#note that PRAS will always replace all H-atoms.
from Pras_Server.PRAS import repairPDB 

#this function replaces only missing heavy atoms.
#always use this function if you only need to replace heavy atoms
from Pras_Server.FixHeavyAtoms import fixheavyAtoms
 
#this function draws the 4 ramachandran types
from Pras_Server.RamaChandra import ramachandranTypes 

#this function assigns the secondary structure elements
from Pras_Server.SecondaryStructure import assignStructure 

#print(fixheavyAtoms.__doc__) to understand the arguments
fixheavyAtoms(pdb_pras='1aho.pdb', rotamer="", mutation="",
	      pdb_faspr="", keep_ligand="", chain_no="")

#out_no_h.pdb is the repaired PDB file written by PRAS. 
#you must have write permission in this directory
assignStructure('out_no_h.pdb')

#out_no_h.pdb is the same as above
ramachandranTypes('out_no_h.pdb')

print ('The program took {0} second !'.format(time.time() - startTime))
```

For linux, do the same as above except for the example code where you should copy the code below (the annotation given above applies below):

```python
import sys
if len(sys.argv) == 7:
	import time
	from Pras_Server.PRAS import repairPDB
	from Pras_Server.FixHeavyAtoms import fixheavyAtoms
	from Pras_Server.RamaChandra import ramachandranTypes
	from Pras_Server.SecondaryStructure import assignStructure
	startTime = time.time()

	fixheavyAtoms(sys.argv[1], sys.argv[2], sys.argv[3],
		      sys.argv[4], sys.argv[5], sys.argv[6])

	assignStructure('out_no_h.pdb')

	ramachandranTypes('out_no_h.pdb')
	
	print ('The program took {0} second !'.format(time.time() - startTime))

else:
	print('PRAS takes 6 compulsory arguments. Execute the code below on your'
	     ' shell prompt to learn more\n')

	print('printf \'from Pras_Server.PRAS import' 
	      ' repairPDB\\nprint(repairPDB.__doc__)\\n\\n()\' | python3\n')
```

Then, cd to the directory where example.py is located and enter the following argument on your shell prompt:

`python3 example.py 1aho.pdb "" "" "" "" ""`

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

Bear in mind that PRAS uses a particular name for the output file and will automatically remove the file during another run. To process multiple files at the same time, concatenate the file name with PRAS dafault name in the appropriate .py file in your installation folder. Alternatively you can download the whole PRAS source code at https://www.protein-science.com/ and modify the software as deemed fit before installation.

PRAS is being considered for publication and the link to the paper will be provided as a reference when published. In the meantime you can cite PRAS as:

O.S. Nnyigide, T.O. Nnyigide, S.G. Lee, K. Hyun. PRAS: A Web Server to Repair PDB Structures, Add Missing Heavy Atoms and Hydrogen Atoms and Assign Secondary Structure by Amide Interactions. Submitted

## License
[MIT](https://choosealicense.com/licenses/mit/)

