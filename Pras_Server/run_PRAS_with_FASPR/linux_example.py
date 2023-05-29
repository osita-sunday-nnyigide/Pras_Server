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
	print("PRAS takes 9 compulsory arguments. Execute the code below on your shell prompt to learn more\n")
	print("printf \"from Pras_Server.PRAS import repairPDB\\nprint(repairPDB.__doc__)\\n\\n()\" | python3\n")
