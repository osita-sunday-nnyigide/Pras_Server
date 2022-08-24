#!/usr/bin/env python

__doc__ = """

This program requires python 3.6 or higher.

This module, test.py can be used to test PRAS

with a total of 576 PDB IDs. There are  two lists

that have been imported,  _82_pdbs, and _494_pdbs

"""

__author__     = "Osita Sunday Nnyigide"

__copyright__  = "Copyright 2022, Osita Sunday Nnyigide"

__credits__    = ["Tochukwu Olunna Nnyigide", "Lee Sun-Gu", "Hyun Kyu"]

__license__    = "MIT"

__version__    = "1.0.8"

__maintainer__ = "Osita Sunday Nnyigide"

__email__      = "osita@protein-science.com"

__status__     = "Production"

__date__       = "May 11, 2022"


import time
from Pras_Server.PRAS import repairPDB
from Pras_Server.PDBID import  _82_pdbs, _494_pdbs
from Pras_Server.FixHeavyAtoms import fixheavyAtoms
from Pras_Server.RamaChandra import ramachandranTypes
from Pras_Server.SecondaryStructure import assignStructure

			# INSTRUCTION !!!!!!!!!!!!!!

##################################################################################
# Note that some of the PDB files have up to 24 chains and each chain will       #
# be analyzed (SS assignment and ramachandran plots).PRAS overwrites each        #
# ouput file, except plots with numbers that the subsequent ones do not get      #
# up to (plot with number ending 7 or 8 cannot be overwritten except another     #
# structure has up to such number of chains, then it overwrites the previous).   #
# The output .pdb or .cif is always overwritten b/c each file is saved with the  #
# same name.                                                                     #
#                                                                                #
#                                                                                #
# ALWAYS REMEMBER THAT PRAS WRITTES IN APPEND MODE, YOU SHOULD ONLY KEEP         #
# PREVIOUS OUTPUT PDB OR mmCIF FILE IF THE NAME IS NOT THE SAME AS ANOTHER!      #
#                                                                                #
# Depending on your computer's memory, processing the entire 494 PDB IDs may     #
# cause your computer to fail to allocate bitmap. So you can either comment      #
# assignStructure() and ramachandranTypes() or continue each time from where     #
# bitmap allocation failed. During our test we processed the whole 494 PDB IDs,  #
# assigned the SS elements and plotted the ramachandran types                    #
#                                                                                #
# Perhaps you should start the test with _82_pdbs which are distinct from        #
# _494_pdbs. When combined together, there are 576 PDB IDs!                      #
#                                                                                #
# Acess to internet is required for this test and the processing speed is subject#
# to your internet download speed.                                               #
##################################################################################



startTime = time.time()


for i in _82_pdbs:

		#  .cif must be passed explicitly as below because the default behaviour is .pdb format,
		#  if you need to process all in .pdb format you do not need to specify the format
		#  rotamer="", always use this default option except you need low occupancy conformers
		#  mutation="", same as above if there is point mutation, use residue with highest occupancy
		#  pdb_faspr="", use this default option except you have provided such pdb obtained from faspr
		#  keep_ligand="", use this default option except you want to keep ligands in the structure
        #  chain_no="", use this default to process all chains or provide an integer or string number
		fixheavyAtoms( i+'.cif', rotamer="", mutation="", pdb_faspr="", keep_ligand="", chain_no="")

		# out_no_h.pdb is the output
		# PDB file written by PRAS
		assignStructure('out_no_h.cif')

		# out_no_h.pdb same as above
		ramachandranTypes('out_no_h.cif')

		print("fixed {}".format(i))


print ('The program took {0} second !'.format(time.time() - startTime))
