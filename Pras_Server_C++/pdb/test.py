from PDBID import _713_PDBs
from DownloadPDB import downloadFile

for i in _713_PDBs:
	downloadFile(i+".cif")