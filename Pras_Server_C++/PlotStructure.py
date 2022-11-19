#!/usr/bin/env python

__author__     = "Osita Sunday Nnyigide"

__copyright__  = "Copyright 2022, Osita Sunday Nnyigide"

__credits__    = ["Tochukwu Olunna Nnyigide", "Lee Sun-Gu", "Hyun Kyu"]

__license__    = "MIT"

__version__    = "1.0.10"

__maintainer__ = "Osita Sunday Nnyigide"

__email__      = "osita@protein-science.com"

__status__     = "Production"

__date__       = "November 13, 2022"

import os
import sys
import time
import itertools
import numpy as np
import matplotlib.pyplot as plt
from pkg_resources import resource_stream
startTime = time.time()
def ramachandranAndSecondaryStructure(cmap = ""):
	for fname in os.listdir(os.getcwd()):
		if fname.startswith("rama_ssa"):
			data = open(fname, 'r').readlines()
			gen = list(itertools.chain(*[line[7:].strip().split('\t') for line in data if line.startswith("General")]))
			gly = list(itertools.chain(*[line[7:].strip().split('\t') for line in data if line.startswith("Glycine")]))
			pro = list(itertools.chain(*[line[7:].strip().split('\t') for line in data if line.startswith("Proline")]))
			pre = list(itertools.chain(*[line[7:].strip().split('\t') for line in data if line.startswith("Pre-Pro")]))

			ah        = [float(line[11:].strip()) for line in data if line.startswith("Alpha_helix")]
			ph        = [float(line[8:].strip())  for line in data if line.startswith("Pi_helix")]
			tth       = [float(line[9:].strip())  for line in data if line.startswith("310_helix")]
			pph       = [float(line[11:].strip()) for line in data if line.startswith("Polyp_helix")]
			strand    = [float(line[6:].strip())  for line in data if line.startswith("Strand")]
			turn      = [float(line[4:].strip())  for line in data if line.startswith("Turn")]
			coil      = [float(line[4:].strip())  for line in data if line.startswith("Coil")]

			general_phi, general_psi = [float(i.split()[0]) for i in gen], [float(i.split()[1]) for i in gen]
			glycine_phi, glycine_psi = [float(i.split()[0]) for i in gly], [float(i.split()[1]) for i in gly]
			proline_phi, proline_psi = [float(i.split()[0]) for i in pro], [float(i.split()[1]) for i in pro]
			pre_pro_phi, pre_pro_psi = [float(i.split()[0]) for i in pre], [float(i.split()[1]) for i in pre]


			cmap = sys.argv[1] if len(sys.argv) > 0 else 'viridis'
			alpha=0.75; dpi=100
			fig, ax = plt.subplots(2,2,figsize=(9.5, 9), dpi=dpi)
			Z = np.fromfile(resource_stream(__name__, 'data/KD.dat'))
			Z = np.reshape(Z, (100, 100))

			ax[0, 0].set_aspect('equal')
			ax[0, 0].set_xlabel('φ')
			ax[0, 0].set_ylabel('ψ')
			ax[0, 0].set_xlim(-180, 180)
			ax[0, 0].set_ylim(-180, 180)
			ax[0, 0].set_xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180], minor=False)
			ax[0, 0].set_yticks([-180, -135, -90, -45, 0, 45, 90, 135, 180], minor=False)
			ax[0, 0].set_title("General")

			ax[1, 0].set_aspect('equal')
			ax[1, 0].set_xlabel('φ')
			ax[1, 0].set_ylabel('ψ')
			ax[1, 0].set_xlim(-180, 180)
			ax[1, 0].set_ylim(-180, 180)
			ax[1, 0].set_xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180], minor=False)
			ax[1, 0].set_yticks([-180, -135, -90, -45, 0, 45, 90, 135, 180], minor=False)
			ax[1, 0].set_title("Proline")

			ax[0, 1].set_aspect('equal')
			ax[0, 1].set_xlabel('φ')
			ax[0, 1].set_ylabel('ψ')
			ax[0, 1].set_xlim(-180, 180)
			ax[0, 1].set_ylim(-180, 180)
			ax[0, 1].set_xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180], minor=False)
			ax[0, 1].set_yticks([-180, -135, -90, -45, 0, 45, 90, 135, 180], minor=False)
			ax[0, 1].set_title("Glycine")

			ax[1, 1].set_aspect('equal')
			ax[1, 1].set_xlabel('φ')
			ax[1, 1].set_ylabel('ψ')
			ax[1, 1].set_xlim(-180, 180)
			ax[1, 1].set_ylim(-180, 180)
			ax[1, 1].set_xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180], minor=False)
			ax[1, 1].set_yticks([-180, -135, -90, -45, 0, 45, 90, 135, 180], minor=False)
			ax[1, 1].set_title("Pre-proline")

			data = np.log(np.rot90(Z))

			ax[0, 0].imshow(data, cmap=plt.get_cmap(cmap), extent=[-180, 180, -180, 180], alpha=alpha)
			ax[1, 0].imshow(data, cmap=plt.get_cmap(cmap), extent=[-180, 180, -180, 180], alpha=alpha)
			ax[0, 1].imshow(data, cmap=plt.get_cmap(cmap), extent=[-180, 180, -180, 180], alpha=alpha)
			ax[1, 1].imshow(data, cmap=plt.get_cmap(cmap), extent=[-180, 180, -180, 180], alpha=alpha)

			data = np.rot90(np.fliplr(Z))
			ax[0, 0].contour(data, colors='k', linewidths=0.5,
			           levels=[10 ** i for i in range(-7, 0)],
			           antialiased=True, extent=[-180, 180, -180, 180], alpha=0.65)
			ax[1, 0].contour(data, colors='k', linewidths=0.5,
			           levels=[10 ** i for i in range(-7, 0)],
			           antialiased=True, extent=[-180, 180, -180, 180], alpha=0.65)
			ax[0, 1].contour(data, colors='k', linewidths=0.5,
			           levels=[10 ** i for i in range(-7, 0)],
			           antialiased=True, extent=[-180, 180, -180, 180], alpha=0.65)
			ax[1, 1].contour(data, colors='k', linewidths=0.5,
			           levels=[10 ** i for i in range(-7, 0)],
			           antialiased=True, extent=[-180, 180, -180, 180], alpha=0.65)

			ax[0,0].scatter(general_phi, general_psi, marker='.', s=3, c="k")
			ax[1,0].scatter(proline_phi, proline_psi, marker='.', s=3, c="k")
			ax[0,1].scatter(glycine_phi, glycine_psi, marker='.', s=3, c="k")
			ax[1,1].scatter(pre_pro_phi, pre_pro_psi, marker='.', s=3, c="k")

			plt.savefig("ramachandra_chain"+fname[-5:-4]+'.jpg')
			plt.close(fig)

			data=[ah[0]if ah else 0.,tth[0]if tth else 0.,ph[0]if ph else 0.,pph[0]if pph else 0.,\
			     strand[0] if strand else 0. ,turn[0] if turn else 0., coil[0] if coil else 0. ]

			SS=['α-helix = {0:.1f}%'.format(data[0]), '3T-helix = {0:.1f}%'.format(data[1]), 
			   'π-helix = {0:.1f}%'.format(data[2]),'Pp-helix = {0:.1f}%'.format(data[3]), \
			   'β-strand = {0:.1f}%'.format(data[4]), 'Turn = {0:.1f}%'.format(data[5]),\
			   'Coil = {0:.1f}%'.format(data[6])]

			x, y = [], []
			for i,j in enumerate(SS):
			  if data[i] != 0.:
			    x.append(data[i])
			    y.append(SS[i])

			fig = plt.figure(figsize =(10, 7))
			plt.pie(x, colors = list('rgbkymc'))
			plt.legend(labels = y, loc='upper center',bbox_to_anchor = (1.1, 1.0),prop={'size': 15})
			plt.savefig('sec_strc_plot_'+'chain'+fname[-5:-4]+'.jpg')
			plt.close(fig)

			os.remove(fname)

print("**Analysing secondary strcutre...\n")

ramachandranAndSecondaryStructure()

print ('The program took {0} second !'.format(time.time() - startTime))