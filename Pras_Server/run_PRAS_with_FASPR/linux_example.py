import sys
if len(sys.argv) == 7:
	import time
	from Pras_Server.PRAS import repairPDB
	from Pras_Server.FixHeavyAtoms import fixheavyAtoms
	from Pras_Server.RamaChandra import ramachandranTypes
	from Pras_Server.SecondaryStructure import assignStructure
	startTime = time.time()
	
	fixheavyAtoms(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4], sys.argv[5], sys.argv[6])

	#assignStructure('out_no_h.pdb')

	#ramachandranTypes('out_no_h.pdb')
	
	print ('The program took {0} second !'.format(time.time() - startTime))

else:
	print("PRAS takes 6 compulsory arguments. Execute the code below on your shell prompt to learn more\n")

	print("printf \"from Pras_Server.PRAS import repairPDB\\nprint(repairPDB.__doc__)\\n\\n()\" | python3\n")