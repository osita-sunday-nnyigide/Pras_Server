#!/usr/bin/env python

__doc__ = """

This program requires python 3.6 or higher.

This module, RunType.py can be used to run PRAS

with different options. 

"""

__author__     = "Osita Sunday Nnyigide"

__copyright__  = "Copyright 2022, Osita Sunday Nnyigide"

__credits__    = ["Tochukwu Olunna Nnyigide", "Lee Sun-Gu", "Hyun Kyu"]

__license__    = "MIT"

__version__    = "1.2.1"

__maintainer__ = "Osita Sunday Nnyigide"

__email__      = "osita@protein-science.com"

__status__     = "Production"

__date__       = "May 11, 2022"

import os
import sys
from Pras_Server.PRAS import repairPDB
from Pras_Server.FixHeavyAtoms import fixheavyAtoms
from Pras_Server.RamaChandra import ramachandranTypes
from Pras_Server.SecondaryStructure import assignStructure

class InitRunType:
	"""

	"""
	fname=None
	def __init__(self, rotamer="", mutation="", pdb_faspr="", keep_ligand="", chain_no="",addh=False, ss=False, raman=False, ofname=False, pdbid=False, his_p=False):
		self.ss=ss
		self.addh=addh 
		self.raman=raman         
		self.rotamer=rotamer
		self.mutation=mutation
		self.chain_no=chain_no
		self.pdb_faspr=pdb_faspr
		self.keep_ligand=keep_ligand
		self.ofname=ofname
		self.pdbid=pdbid
		self.his_p=his_p

	def ProcessWithoutDefault(self):
		"""
		Use this member function if you want to process all the PDBs in the directory where you called this program.
		Here, the name of the output file will be the name of your input file minus ending .pdb plus _out.
		In the case of a rerun PRAS will remove previous output files because PRAS writes in append mode

		The log file has the name of each PDB repaired so that you will know which PDB structure has missing atoms
		In the case of a rerun the log file will be removed for the reasons stated above.

		As an example, copy and run the below code in a python document in your directory with a PDB file.

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

		fixing.ProcessWithoutDefault()

		print ('The program took {0} second !'.format(time.time() - startTime))

		"""
		for i in os.listdir(os.getcwd()):
			if i[-8:-4] == '_out' or i == 'log.txt': os.remove(i)

		for i in os.listdir(os.getcwd()):
			if i.endswith(".pdb") or i.endswith(".ent"):
				if self.addh:repairPDB( i, self.rotamer, self.mutation, self.pdb_faspr, self.keep_ligand, self.chain_no, i[:-4]+'_out.pdb', self.his_p)
				else:fixheavyAtoms( i, self.rotamer, self.mutation, self.pdb_faspr, self.keep_ligand, self.chain_no, i[:-4]+'_out.pdb')
				if self.ss:assignStructure(i[:-4]+'_out.pdb')
				if self.raman:ramachandranTypes(i[:-4]+'_out.pdb')
				print("fixed {}".format(i))

			elif  i.endswith(".cif"):
				if self.addh:repairPDB( i, self.rotamer, self.mutation, self.pdb_faspr, self.keep_ligand, self.chain_no, i[:-4]+'_out.cif', self.his_p)
				else:fixheavyAtoms( i, self.rotamer, self.mutation, self.pdb_faspr, self.keep_ligand, self.chain_no, i[:-4]+'_out.cif')
				if self.ss:assignStructure(i[:-4]+'_out.cif')
				if self.raman:ramachandranTypes(i[:-4]+'_out.cif')
				print("fixed {}".format(i))

	def ProcessWithoutDefaultUsingPDBID(self):
		"""
		Use this member function if you want to process all the PDBs in a PDB ID list supplied to PRAS.
		Here, the name of the output file will be the name of your input PDB ID plus _out.
		In the case of a rerun PRAS will remove previous output files because PRAS writes in append mode

		The log file has the name of each PDB repaired so that you will know which PDB structure has missing atoms
		In the case of a rerun the log file will be removed for the reasons stated above.

		As an example, copy and run the below code in a python document in your current working directory.

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
							pdbid=['1crn','1hen'],
							his_p=False
							)

		fixing.ProcessWithoutDefaultUsingPDBID()

		print ('The program took {0} second !'.format(time.time() - startTime))

		"""
		if not isinstance(self.pdbid, list):
			print("You need to supply a list of PDB IDs. PRAS has terminated abnormally")
			sys.exit(1)
		for i in os.listdir(os.getcwd()):
			if i[-8:-4] == '_out' or i == 'log.txt': os.remove(i)

		for i in self.pdbid:
			if not i.endswith(".cif"):
				if self.addh:repairPDB( i, self.rotamer, self.mutation, self.pdb_faspr, self.keep_ligand, self.chain_no, i+'_out.pdb', self.his_p)
				else:fixheavyAtoms( i, self.rotamer, self.mutation, self.pdb_faspr, self.keep_ligand, self.chain_no, i+'_out.pdb')
				if self.ss:assignStructure(i+'_out.pdb')
				if self.raman:ramachandranTypes(i+'_out.pdb')
				print("fixed {}".format(i+'.pdb'))

			elif  i.endswith(".cif"):
				if not self.ofname:self.ofname=i[:-4]+'_out.cif'
				if self.addh:repairPDB( i, self.rotamer, self.mutation, self.pdb_faspr, self.keep_ligand, self.chain_no, i[:-4]+'_out.cif', self.his_p)
				else:fixheavyAtoms( i, self.rotamer, self.mutation, self.pdb_faspr, self.keep_ligand, self.chain_no, i[:-4]+'_out.cif')
				if self.ss:assignStructure(i[:-4]+'_out.cif')
				if self.raman:ramachandranTypes(i[:-4]+'_out.cif')
				print("fixed {}".format(i))

	def ProcessWithDefault(self):
		"""
		Use this member function if you want to process a single PDB file using PRAS default output name.
		In the case of a rerun PRAS will remove previous output files because PRAS writes in append mode

		As an example, copy and run the below code in a python document in your directory with a PDB file 

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

		fixing.fname = '1crn.pdb' # change this to whatever your file is named

		fixing.ProcessWithDefault()

		print ('The program took {0} second !'.format(time.time() - startTime))

		"""
		if not self.fname:
			print("You need to supply a file name. PRAS has terminated abnormally")
			sys.exit(1)
		if self.ofname:
			print("Output file name must be set to False. PRAS has terminated abnormally")
			sys.exit(1)
		if self.fname.endswith(".pdb") or self.fname.endswith(".ent"):
			if self.addh:repairPDB(self.fname, self.rotamer, self.mutation, self.pdb_faspr, self.keep_ligand, self.chain_no, self.ofname, self.his_p)
			else:fixheavyAtoms(self.fname, self.rotamer, self.mutation, self.pdb_faspr, self.keep_ligand, self.chain_no)
			if self.ss:assignStructure('out.pdb')
			if self.raman:ramachandranTypes('out.pdb')
			print("fixed {}".format(self.fname))

		if self.fname.endswith(".cif"):
			if self.addh:repairPDB(self.fname, self.rotamer, self.mutation, self.pdb_faspr, self.keep_ligand, self.chain_no, self.ofname, self.his_p)
			else:fixheavyAtoms(self.fname, self.rotamer, self.mutation, self.pdb_faspr, self.keep_ligand, self.chain_no)
			if self.ss:assignStructure('out.cif')
			if self.raman:ramachandranTypes('out.cif')
			print("fixed {}".format(self.fname))

	def ProcessOther(self):
		"""
		Use this member function if you want to give a specific name to your output file.
		Here, you will loop over all your PDB structures in the directory where you called 
		this program; that is to say that you are using a loop by yourself and for each PDB, 
		you instantiate InitRunType class and supply the ofname, say xxxx.pdb.

		Warning: PRAS writes in append mode and you're responsibe for deleting the output files in case of a rerun!!

		As an example, copy and run the below code in a python document in your directory with many PDBs

		import os
		import sys
		import time
		from Pras_Server.RunType import InitRunType

		startTime = time.time()
		out = ['xxxx.pdb','zzzz.pdb','ssss.pdb'] # size of this list MUST be equal to the total number of PDB structures in your working directory
		k=0
		for i in os.listdir(os.getcwd()):
			if i.endswith(".pdb") or i.endswith(".ent"):
				try:
					out[k]
				except:
					print('Index error. Size of the name list is small')
					sys.exit()
				fixing = InitRunType(rotamer="", 
									mutation="", 
									pdb_faspr="", 
									keep_ligand="", 
									chain_no="",
									addh = False,
									ss=False,
									raman=False,
									ofname=out[k],
									pdbid=False,
									his_p=False
									)
				fixing.fname = i
				fixing.ProcessOther()
				k+=1
		print ('The program took {0} second !'.format(time.time() - startTime))

		"""
		if not self.fname:
			print("No PDB structure file found in your working direcotry. PRAS has terminated abnormally")
			sys.exit()
		if self.fname.endswith(".pdb") or self.fname.endswith(".ent"):
			if self.addh:repairPDB(self.fname, self.rotamer, self.mutation, self.pdb_faspr, self.keep_ligand, self.chain_no, self.ofname, self.his_p)
			else:fixheavyAtoms(self.fname, self.rotamer, self.mutation, self.pdb_faspr, self.keep_ligand, self.chain_no,self.ofname)
			if self.ss:assignStructure(self.ofname)
			if self.raman:ramachandranTypes(self.ofname)
			print("fixed {}".format(self.ofname))

		if self.fname.endswith(".cif"):
			if self.addh:repairPDB(self.fname, self.rotamer, self.mutation, self.pdb_faspr, self.keep_ligand, self.chain_no, self.ofname, self.his_p)
			else:fixheavyAtoms(self.fname, self.rotamer, self.mutation, self.pdb_faspr, self.keep_ligand, self.chain_no, self.ofname)
			if self.ss:assignStructure(self.ofname)
			if self.raman:ramachandranTypes(self.ofname)
			print("fixed {}".format(self.ofname))
