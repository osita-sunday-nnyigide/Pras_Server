#!/usr/bin/env python

__doc__ = """

This program requires python 3.6 or higher.

This module, SecondaryStructure.py has the function

that is used to assign secondary structure elements.

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

import sys
import itertools
from . import ReadStructure
import matplotlib.pyplot as plt

def assignStructure(filename=False):
	"""
	This program assigns secondary structure elements and
	writes all data of structure assignment in a text file
	and plots the value in TIF format

	This funtion & the dependency (ReadStructure.py) read PRAS
	output PDB file in a quick and complex manner. Since it is a
	PRAS output there will not be error (missing atoms, bad format, etc).

	In general, this program can be used for any standard PDB file.
	However, if there are missing carbonyl O atom or a bad format,
	this program will run into trouble. Thus, it is better to run the PDB
	through PRAS first and then use the PRAS output here.


	Arguments
	----------
	filename: PRAS generated PDB structure file (.pdb or .cif)

	Returns
	-------
	None:     writes a text file containing detailed secondary
	          structure and plots a figure  in TIF format
	"""
	if not filename:
		filename="out.pdb" 
	_format = filename[-4:]

	if _format == '.cif':
		try:
			lines = open(filename, 'r').readlines()
			rf              = [i.strip("\n") for i in lines if i[:4] == 'ATOM']
			chains          = list(dict.fromkeys([i.split()[6] for i in rf]))
			atoms_per_chain = [[j for j in rf if j.split()[6] == chains[i]]\
							  for i in range(len(chains))]
		except:
			print('{} could not be read, check the file and format.'
				 '\n' 'Perhaps it should be a pdb but you specified cif.'
				  ' Secondary''\n''structure could not be assigned'.format(filename))
			sys.exit(1)

	elif _format == '.pdb':
		try:
			lines   = open(filename, 'r').readlines()
			rf      = [i.strip("\n") for i in lines if i[:4] == 'ATOM' or i[0:3] == "TER"]
			chains  = [i for i,j in enumerate(rf) if j[:3].strip(' ')  == "TER"]
		except:
			print('{} could not be read, check the file and format.'
				 '\n' 'Perhaps it should be a cif but you specified pdb.'
				  ' Secondary''\n''structure could not be assigned'.format(filename))
			sys.exit(1)

	else:
		print('Fatal error! Please check file type')
		sys.exit(1)

	count = 0
	for n in range(len(chains)):
		chain_i = [i for i in list(itertools.chain(*[rf[count:i]for i,j in enumerate(rf) if i == chains[n]]))]\
			   if filename[-4:] != '.cif' else atoms_per_chain[n]

		ahelix_hb    = ReadStructure.alpha_helix(chain_i, _format)
		amide_dih    = ReadStructure.amideDihedral(chain_i, _format)
		_310helix_hb = ReadStructure._310helix(chain_i, _format)
		bturn_nhb    = ReadStructure.beta_turn(chain_i,  _format)
		phi_dih      = ReadStructure.phiDihedral(chain_i, _format)
		psi_dih      = ReadStructure.psiDihedral(chain_i, _format)
		pihelix_hb   = ReadStructure.pi_helix(chain_i, _format)

		nextres      = ReadStructure.nextRes(chain_i, _format)
		chain        = ReadStructure.resChain(chain_i, _format)
		resn         = ReadStructure.resName(chain_i, _format)
		resNo        = ReadStructure.resNumber(chain_i, _format)

		"""
		initially assign pi helix, alpha helix, 310 helix and the rest Coil
		using the recognition criterion as described in the PRAS paper
		"""
		v = list(zip(ahelix_hb,_310helix_hb,amide_dih,bturn_nhb,phi_dih,psi_dih,pihelix_hb))
		SSA = []
		for i in range(len(v)):

			if  (v[i][2] > -10 and v[i][2] < 45 and v[i-1][6] < 3.5 and v[(i+1)%len(v)][2] > -10\
				and v[(i+1)%len(v)][2] < 45 and v[(i+2)%len(v)][2] > -10 and v[(i+2)%len(v)][2] < 45)\
				and ((v[i-1][6] < v[i-1][0] + 0.1) or (v[i][6] < v[i][0] + 0.1 and v[(i+1)%len(v)][6] <\
				v[(i+1)%len(v)][0] + 0.1)):
				SSA.append('Pi_helix')

			elif  (v[i][2] > -10 and v[i][2] < 45 and v[i-1][0] < 3.5 and v[(i+1)%len(v)][2] > -10\
				  and v[(i+1)%len(v)][2] < 45 and v[(i+2)%len(v)][2] > -10 and v[(i+2)%len(v)][2] < 45)\
			      and ((v[i-1][0] < v[i-1][1] + 0.1) or (v[i][0] < v[i][1] + 0.1 and\
				  v[(i+1)%len(v)][0] < v[(i+1)%len(v)][1] + 0.1)):
				SSA.append('Alpha_helix')

			elif  (v[i][2] > -10 and v[i][2] < 45 and v[i-1][1] < 3.5 and v[(i+1)%len(v)][2] > -10 and\
				  v[(i+1)%len(v)][2] < 45 and v[(i+2)%len(v)][2] > -10 and v[(i+2)%len(v)][2] < 45) and\
			      ((v[i-1][1] < v[i-1][0] + 0.1) or (v[i][1] < v[i][0] + 0.1 and v[(i+1)%len(v)][1] <\
			      v[(i+1)%len(v)][0] + 0.1)):
				SSA.append('310_helix')

			else:
				SSA.append('Coil')

		"""
		extend each pi helix region by 5,
		the 5th (first resi) is assigned already
		"""
		extend_pihelix = []
		for i, extend_helix in enumerate(SSA):
			_helix_dict ={}
			if (SSA[i] == 'Pi_helix') and (SSA[(i+1)% len(SSA)] != 'Pi_helix') and len(SSA[i:])>4:
				_helix_dict[extend_helix] = i
				extend_pihelix.append(_helix_dict)
		for helix in extend_pihelix:
			for key, value in helix.items() :
				SSA[value+1]   = 'Pi_helix'
				SSA[value+2]   = 'Pi_helix'
				SSA[value+3]   = 'Pi_helix'
				SSA[value+4]   = 'Pi_helix'

		"""
		extend each alpha helix region by 4,
		the 4th (first resi) is assigned already
		"""
		extend_ahelix = []
		for i, extend_helix in enumerate(SSA):
			_helix_dict ={}
			if (SSA[i] == 'Alpha_helix') and (SSA[(i+1)% len(SSA)] != 'Alpha_helix'\
				and SSA[(i+1)% len(SSA)] != 'Pi_helix') and len(SSA[i:])>3:
				_helix_dict[extend_helix] = i
				extend_ahelix.append(_helix_dict)

		for helix in extend_ahelix:
			for key, value in helix.items() :
				SSA[value+1]   = 'Alpha_helix'
				SSA[value+2]   = 'Alpha_helix'
				SSA[value+3]   = 'Alpha_helix'

		"""
		extend each 310helix region by 3,
		the 3rd (first resi) is assigned already
		"""
		extend_3helix = []
		for i, extend_310helix in enumerate(SSA):
			_helix_dict ={}
			if (SSA[i] == '310_helix') and (SSA[(i+1)% len(SSA)] != '310_helix') and len(SSA[i:])>2:
				_helix_dict[extend_310helix] = i
				extend_3helix.append(_helix_dict)
		for three10helix in extend_3helix:
			for key, value in three10helix.items() :
				SSA[value+1]   = '310_helix'
				SSA[value+2]   = '310_helix'

		# assign beta-strand. Check incursion into helices
		_strand_list = []
		for i in range(len(v)):
			_strand_dict ={}
			if  (((v[i][2] > 173 and v[i][2] < 181) or (v[i][2] < -140)) and ((v[(i+1)%len(v)][2] >=\
			 	173 and v[(i+1)%len(v)][2] < 181) or (v[(i+1)%len(v)][2] <= -140))) or (v[i][4] < -120\
			 	and v[i][5] > 120 and SSA[i] != 'Alpha_helix'and SSA[i] != '310_helix' and SSA[i]\
			 	!= 'Pi_helix') and (SSA[(i+1)% len(SSA)] != 'Alpha_helix' and SSA[(i+1)% len(SSA)]\
			 	!= 'Pi_helix' and SSA[(i+1)% len(SSA)] != '310_helix') and len(SSA[i:])>2:
				_strand_dict[SSA[i]] = i
				_strand_list.append(_strand_dict)
		for beta_strand in _strand_list:
			for key, value in beta_strand.items():
				SSA[value]     = 'Strand'
				SSA[value+1]   = 'Strand'

		# assign beta-turn. Check incursion into helices but not beta-strand
		_turn_list = []
		for i, beta_turn in enumerate(SSA):
			_turn_dict ={}
			if  (v[i][1] <= 3.5) and (SSA[i]!= 'Alpha_helix' and SSA[i] != '310_helix' and SSA[i] != 'Pi_helix')\
			 	and (SSA[(i+1)% len(SSA)] != 'Alpha_helix' and SSA[(i+1)% len(SSA)]!='Pi_helix'and SSA[(i+1)% len(SSA)]\
			 	!= '310_helix') and len(SSA[i:])>2:
				_turn_dict[beta_turn] = i
				_turn_list.append(_turn_dict)
		for beta_turn in _turn_list:
			for key, value in beta_turn.items() :
				SSA[value]   = 'Turn_hb'
				SSA[value+1] = 'Turn_hb'

		# assign nhb beta-turn. Check incursion into helices, hb turn and beta-strand
		_nhbturn_list = []
		for i, nhb_beta_turn in enumerate(SSA):
			_nhbturn_dict ={}
			if  (v[i][3] <= 4.7) and (SSA[i]!= 'Alpha_helix' and SSA[i] != '310_helix' and SSA[i] != 'Pi_helix' and SSA[i]\
				!= 'Turn_hb'and SSA[i] != 'Strand') and (SSA[(i+1)% len(SSA)] != 'Alpha_helix' and SSA[(i+1)% len(SSA)] \
				!= 'Pi_helix' and SSA[(i+1)% len(SSA)]!= '310_helix') and len(SSA[i:])>2:
				_nhbturn_dict[nhb_beta_turn] = i
				_nhbturn_list.append(_nhbturn_dict)
		for nhb_beta_turn in _nhbturn_list:
			for key, value in nhb_beta_turn.items() :
				SSA[value]   = 'Turn_nhb'
				SSA[value+1] = 'Turn_nhb'

		# assign pp helix. Check incursion into the others
		_pphelix_list = []
		for i, pp_helix in enumerate(SSA):
			_pphelix_dict ={}
			if (v[i][2] < -150 and v[i][2] < -80) and (v[i][4] < -65 and v[i][4] > -85 and v[i][5] < 175 and v[i][5] > 140)\
			   and (v[i-1][4] < -65 and v[i-1][4] > -85 and v[i-1][5] < 175 and v[i-1][5] > 140) and (SSA[i] != 'Alpha_helix')\
			   and (SSA[i] != '310_helix') and (SSA[i] != 'Pi_helix') and (SSA[i] != 'Turn_hb') and (SSA[i] != 'Turn_nhb') \
			   and (SSA[i] != 'Strand') and len(v[i:])>1:
				_pphelix_dict[pp_helix] = i
				_pphelix_list.append(_pphelix_dict)
		for pp_helix in _pphelix_list:
			for key, value in pp_helix.items():
				SSA[value]     = 'Polyp_helix'

		# count all assigned SS elements and write to output file
		total_coil       =  SSA.count('Coil')
		total_nhb_turn   =  SSA.count('Turn_nhb')
		total_hb_turn    =  SSA.count('Turn_hb')
		total_strand     =  SSA.count('Strand')
		total_alphahelix =  SSA.count('Alpha_helix')
		total_310helix   =  SSA.count('310_helix')
		total_pihelix    =  SSA.count('Pi_helix')
		total_pphelix    =  SSA.count('Polyp_helix')
		pdb_code         =  filename


		with open ('sec_strc_table_'+'chain'+str(n+1)+'_'+filename[:-4]+'.txt', 'w') as f:
			f.write('When using files generated by this program in a publication, please cite this program as "O.S. Nnyigide, T.O. Nnyigide,'
				' S.G. Lee, K. Hyun. Protein Repair and Analysis Server: A Web' +'\n'+'Server to Repair PDB Structures, Add Missing Heavy'
				' Atoms and Hydrogen Atoms, and Assign Secondary Structures by Amide Interactions. J. Chem. Inf. Model., 2022, 62, 4232–4246.'+'\n')
			f.write('#'*183 +'\n');f.write('Secondary Structure Assignment by Amide-Amide Interactions of the Backbone.'
			' This is the summary of the measurements used to assign the secondary structure elements' +'\n');f.write('#'*183+'\n')
			f.write('Residues in PDB \t\t' 'Amide_dihedral\t' 'Phi_dihedral\t' 'Psi_dihedral\t' 'O(i)->N(i+4)\t' 'O(i)->N(i+3)\t'\
				    'O(i)->N(i+5)\t' 'C(i)->N(i+3)\t' 'Dihedral_type\t''Secondary structure\n')

			for (chain, resn, resNo, amide_dih, phi, psi, ahelix, _310helix,pihelix, beta_turn, res, ss) in zip(chain,resn,resNo,\
				amide_dih,phi_dih,psi_dih,ahelix_hb,_310helix_hb,pihelix_hb,bturn_nhb,nextres,SSA):
				f.write("%s:Chain%s:%s%s\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8s\t%8s\n" %(pdb_code, chain, resn, resNo,\
				 		amide_dih, phi, psi, ahelix, _310helix, pihelix, beta_turn, res, ss))

			data=[(total_alphahelix/len(SSA))*100,(total_310helix/len(SSA))*100,(total_pihelix/len(SSA))*100,(total_pphelix/len(SSA))*100,\
			 	 (total_strand/len(SSA))*100,((total_hb_turn+total_nhb_turn)/len(SSA))*100,(total_coil/len(SSA))*100]

			SS=['α-helix = {0:.1f}%'.format(data[0]), '3T-helix = {0:.1f}%'.format(data[1]), 'π-helix = {0:.1f}%'.format(data[2]),\
				  'Pp-helix = {0:.1f}%'.format(data[3]), 'β-strand = {0:.1f}%'.format(data[4]), 'Turn = {0:.1f}%'.format(data[5]),\
			 	  'Coil = {0:.1f}%'.format(data[6])]

			x, y = [], []
			for i,j in enumerate(SS):
			  if data[i] != 0.000:
			    x.append(data[i])
			    y.append(SS[i])

			# Create plot
			fig = plt.figure(figsize =(10, 7))
			plt.pie(x, colors = list('rgbkymc'))
			plt.legend(labels = y, loc='upper center',bbox_to_anchor = (1.1, 1.0),prop={'size': 15})

			# save plot
			plt.savefig('sec_strc_plot_'+'chain'+str(n+1)+'_'+filename[:-4]+'.tif')
			plt.close(fig)

			if _format == '.pdb':
				if n == 0:
					count+= chains[n]+1
				else:
					count+= (chains[n]-chains[n-1])
