Protein Repair and Analysis (PRAS) Server is a library to repair PDB or mmCIF structures, add missing heavy atoms and hydrogen atoms and assign secondary structures by amide-amide interactions of the backbone

The program was originally written in python programming language. Comments are not provided in the C++ source code because it follows the logic of the python source code where extensive comments are provided

The difference between this and the python counterpart is the processing speed. The C++ version is written for those processing hundreds of PDB structures or a very large PDB/mmCIF structure.  For this reason, the option to select a single chain is limited to the python version which is reasonably fast for that kind of job. In summary, the C++ version is 9 to 12 times faster than the python version.

You can use the PRAS server online at https://www.protein-science.com/ or on your local machine following the instructions in the INSTALL file included with this distribution.

PRAS has been peer reviewed and published. Please cite PRAS as:

O.S. Nnyigide, T.O. Nnyigide, S.G. Lee, K. Hyun. Protein Repair and Analysis Server: A Web Server to Repair PDB Structures, Add Missing Heavy Atoms and Hydrogen Atoms, and Assign Secondary Structures by Amide Interactions.
J. Chem. Inf. Model., 2022, 62, 4232–4246.

## License
[MIT](https://choosealicense.com/licenses/mit/)

