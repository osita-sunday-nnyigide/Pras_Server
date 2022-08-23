# computational time varies
import time

# this function adds both heavy and H-atoms
from Pras_Server.PRAS import repairPDB

# this function replaces only missing heavy atoms
from Pras_Server.FixHeavyAtoms import fixheavyAtoms

# this function draws the 4 ramachandran types
from Pras_Server.RamaChandra import ramachandranTypes

# this function assigns the secondary structure elements
from Pras_Server.SecondaryStructure import assignStructure

startTime = time.time()


#################################################################################
# PRAS automatically downloads the structure as specified if it does not exists #
# in the directory where this script is located.                                #
# User must have permission to write in this directory.                         #
# Processing speed will include the internet connection and download speed      #
#################################################################################

#  .cif must be passed explicitly as below because the default behaviour is .pdb format,
#  if you need to process structure in .pdb format you do not need to specify the format.
#  rotamer="", always use this default option except you need low occupancy conformers
#  mutation="", same as above if there is point mutation, use residue with highest occupancy
#  pdb_faspr="", use this default option except you have provided such pdb obtained from faspr
#  keep_ligand="", use this default option except you want to keep ligands in the structure
#  chain_no="", use this default to process all chains or provide an integer or string number
fixheavyAtoms( pdb_pras='1aho.cif', rotamer="", mutation="", pdb_faspr="", keep_ligand="", chain_no="")

# out_no_h.cif is the repaired PDB file written by PRAS, if you use repairPDB() it becomes out_with_h.cif
assignStructure('out_no_h.cif')

# out_no_h.cif or out_with_h.cif is same as above, cmap options are viridis, magma, inferno, jet, plasma.
ramachandranTypes('out_no_h.cif', cmap = "")

print ('The program took {0} second !'.format(time.time() - startTime))
