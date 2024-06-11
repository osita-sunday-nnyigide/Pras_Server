#!/usr/bin/env python

__doc__ = """

This program requires python 3.6 or higher.

This module, CheckPDBatoms.py has the function

that is used to check for missing heavy atoms.

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
from .MissingHeavyAtoms import *
from .ReadMaster import getChains
from .ForcefieldParam import param
from .ReadmmCIF import get_mmcif_chains


def keepLigands(chains, chain_no, _format):

    """
    The default behaviour of PRAS is to generate
    clean protein structures. This is because ligands
    are specially treated in docking or MD simulation.

    However, if for any reason one intends to keep
    the ligands (non-water ligands), this function will
    extract the ligand coordinates to be appended to
    the particular chain it belongs to in the PDB
    structure file

    Arguments
    ----------
    chains:          all the chains in the PDB structure

    chain_no:        chain number. If user selects a  particular
                     chain then the ligand in that particular
                     chain is exracted, otherwise all ligands
                     will be extracted

    _format:         the file format. If cif, the HETATM
                     format is different compared with pdb

    Returns
    -------
    A nested list:   index 0 = chains, index 1 = ligand
    """

    ligand = chains[1]

    if chain_no:
        try:
            chain_no = int(chain_no) if type(chain_no) is str else chain_no
        except:
            invalidChain(chain_no)

        all_chains = [i[-1][0] for i in chains[0]]
        ligand_per_chain = [[] for _ in range(len(chains[0]))]

        for i in range(len(all_chains)):
            for j in ligand:

                if _format != '.cif':
                    if j[21] == all_chains[i]:
                        ligand_per_chain[i].append(j)
                else:
                    # check if _atom_site.auth_asym_id is provided
                    # which seems to have priority based on some
                    # .cif files where _atom_site.label_asym_id is
                    # different from _atom_site.auth_asym_id for
                    # the same entry.
                    try:
                        if j.split()[18] == all_chains[i]:
                            ligand_per_chain[i].append(j)

                    # if _atom_site.auth_asym_id is not provided,
                    # use _atom_site.label_asym_id which is always
                    # provided (_atom_site.label_asym_id is compulsory!).
                    except:
                        if j.split()[6] == all_chains[i]:
                            ligand_per_chain[i].append(j)

        ligand = ligand_per_chain[chain_no-1]

        return [chains[0], ligand]

    return [chains[0], ligand]

def repairNotes(resNo,chain,missing):

    """
    PRAS reminds users of residues whose atoms are fixed.
    This will help user check the quality of the fix.

    Arguments
    ----------
    resNo:   the amino acid residue number
    chain:   the chain or chain number
    missing: a list of the missing atoms

    Returns
    -------
    None:    it writes to a text file
    """
    with open('log.txt', 'a') as log:
        log.write('The residue No {} in chain {} has '
                  'missing atom(s) {}. We have fixed it'
                  .format(resNo, chain,missing)+'\n\n')

def non_standard_residue(resNo,resn,chain):

    """
    The expected standard amino-acid residue names
    are given in a dictionary at ForcefieldParam.py.

    If user is using PDB file obtained from
    GROMACS, CHARMM, AMBER or other MD programs,
    histidine may be named HSD, HSE, HIE or HSP.
    PRAS will use it as normal HIS because
    HIS protonation is JUST as described in the PRAS paper

    If PRAS encounters a non-standard residue, the program
    will be terminated and user will be informed

    Arguments
    ----------
    resNo:   the amino acid residue number
    resn:    the amino acid residue name
    chain:   the chain or chain number

    Returns
    -------
    None:    a premature termination of the program
    """
    with open('non_standard_residue.txt', 'w') as f:
        f.write('The residue {}  named {} in chain {} is a'
                ' non-standard residue. Program has terminated'
                ' abnormally'.format(resNo,resn,chain)+'\n\n')

        sys.exit(1)

def invalidChain(chain_no):

    """
    This function checks if an invalid chain is selected
    and then terminate the program.

    Returns
    -------
    None:   a premature termination of the program
    """
    with open('invalid_chain.txt', 'w') as f:
        f.write('You selected chain No. {} which does not exist.'
            ' Chain No. must be a positive number that is less'
            ' than or equal the total number of chain in the PDB'
            ' structure file'.format(chain_no)+'\n')

        sys.exit(1)

def missing_backbone_atom(resNo,chain):

    """
    PRAS can fix missing backbone O but demands that N, CA, & C be present.
    This function will terminate the program if a main-chain atom is missing.

    Arguments
    ----------
    resNo:   the amino acid residue number
    chain:   the chain or chain number

    Returns
    -------
    None:    a premature termination of the program
    """
    with open('missing_backbone_atom.txt', 'w') as f:
        f.write('The residue No {} in chain {} is missing a'
            ' backbone heavy atom. Program has terminated'
            ' abnormally. Note that for .cif the reported'
            ' sequence number is the label_seq_id which may'
            ' or may not be the same as auth_seq_id'.
            format(resNo,chain)+'\n\n')

        sys.exit(1)

def invalid_pdb_file():

    """
    PRAS will read PDBs with as much flaws as covered in the program.
    Unfortunately, if the file is severely damaged or illegal, PRAS will fail.

    Returns
    -------
    None:   a premature termination of the program
    """
    with open('invalid_pdb_file.txt', 'w') as f:
        f.write('Error occured while reading the PDB file'+'\n')
        f.write('For a .pdb format, ensure that each chain ends'
               ' with TER and the last chain TER+END.'+'\n')
        f.write('Also, visit www.protein-science.com and'
               ' click "How-tos" for more information')

        sys.exit(1)

def checkpdbAtoms (pdb_pras, rotamer, mutation, pdb_faspr, keep_ligand, chain_no, ofname=False):

    """
    This function obtains all chains in the PDB file as read by ReadMAster.py or ReadmmCIF.py
    It then loops through each chain and checks if a residue has missing heavy atom(s).
    If it has, this function will replace it by calling MissingHeavyAtoms.py
    See PRAS.py for the detailed explanation of the arguments passed.

    Arguments
    ----------
    pdb_pras:    the PDB file you intend to repair/analyze

    rotamer:     supply "no" if you need to generate lower occupancy conformers,
                 if not supply "". PRAS uses atoms with highest occupancy by default.

    mutation:    supply "no" if you need to generate lower occupancy conformers,
                 if not supply "". PRAS uses the residue with highest occupancy by default
                 if two residues are given same residue number

    pdb_faspr:   the output PDB obtained by running your PDB with FASPR

    keep_ligand: by default, ligands and entries with HETATM are ignored because in most
                 MD simulation and docking studies the ligand is specially treated. Moreover,
                 PRAS does not repair ligands (ligand structure is diverse compared with protein
                 with 20 repeating common amino acids!). However, user can supply any string as
                 the argument if user intends to keep ligands.

    chain_no   : by default, all the chains in the PDB will be processed. However, if user intends to
                 use a specific chain, user should provide an integer or a string number as the argument.
                 PRAS will map the integer to alphabets. Chain can start with any letter. Thus,
                 1 = the first chain, 2 = the second chain, etc. If user supplies a string number PRAS
                 will convert it to integer but if alphabet error will be generated.

    ofname     : the output file name, default name is out.pdb or out.cif

    Returns
    -------
    A list:    a list that contains all chains with complete heavy atom of each residue
    """


    """
    Remove previous log file in
    the case of another run
    """
    if not ofname:
        try:
            os.remove('log.txt')
        except:
            pass

    """
    Now initialize the list to grab all
    chains with all missing atoms fixed.
    """
    ter_all = []

    with open('log.txt', 'a', encoding="utf8") as log:

        log.write('When using files generated by this program in a publication,'
                ' please cite this program as "O.S. Nnyigide, T.O. Nnyigide,'
                ' S.G. Lee, K. Hyun. Protein Repair and Analysis Server: A Web' +'\n'+
                'Server to Repair PDB Structures, Add Missing Heavy'
                ' Atoms and Hydrogen Atoms, and Assign Secondary Structures by'
                ' Amide Interactions. J. Chem. Inf. Model., 2022, 62, 4232–4246.'+'\n')

        log.write('#'*183+'\n')

        log.write('PRAS 1.2.1. This is a PRAS-generated log file.'
                ' For your information, all missing or fixed atoms'
                ' and other relevant information concerning the repair'
                ' are appended below '+'\n')

        log.write('#'*183+'\n\n')
        log.write('This is the log for {}'.format(pdb_pras)+'\n\n')

    """
    Read user's input file obtained from faspr but if
    it is not readable inform user and use default chi
    from Dunbrack 2011 rotamer library
    """
    try:
        faspr = getChains("", "", "", pdb_faspr, keep_ligand)
        faspr = faspr[0]
    except:
        faspr = []

    if not faspr:
        with open('log.txt', 'a') as log:
            log.write('No FASPR PDB output supplied.'
                      ' Default chi will be used'+'\n\n'
                      )

    """
    Read user's input file to be repaired or analyzed
    but if it is not readable inform user and terminate
    the program
    """
    try:
        # the default behaviour is pdb file. To use cif, indicate it in the file name.
        # the program will automatically download the file if it is not in the user's
        # directory where PRAS is run.
        if pdb_pras[-4:] != '.cif':
            chains = getChains(
                                pdb_pras, rotamer, mutation,
                                pdb_faspr, keep_ligand

                                )
        else:
            chains = get_mmcif_chains(

                                    pdb_pras, rotamer,
                                    mutation, keep_ligand
                                 )

    # most errors come from using .pdb or .cif
    # files that do not follow the specified standard.
    # In most cases this is a file obtained from a
    # third-party tool.
    except:
        chains = []

    if not chains or not chains[0]:
        invalid_pdb_file()

    """
    If user spicified to keep ligands, keepLigands()
    will be called to add ligands (HETATM) to the
    particular chain it belongs. If there are no ligands
    but user specified one, PRAS will ignore user rather
    than raise exception
    """
    if keep_ligand:
        try:
            chains, ligand = keepLigands(chains, chain_no, pdb_pras[-4:])
        except:
           chains, ligand = chains[0], []
    else:
        chains = chains[0]

    """
    Specify the particular chain if you need to, otherwise all
    chains will be processed.The chain is specified as an integer
    or a string number. PRAS will map the integer to alphabets.
    Chain  can start with any letter. Thus, 1 = the first chain,
    2 = the second chain, etc.
    """
    if chain_no:
        try:
            chain_no = int(chain_no) if type(chain_no) is str else chain_no
            chains = [chains[chain_no-1]]
        except:
            invalidChain(chain_no)

    # now loop through the whole chain and fix the missing
    # atoms in each chain by calling MissingHeavyAtoms.py

    for n in range(len(chains)):
        resn,atoms,atmpos,resNo,chain_id = chains[n][0],\
        chains[n][1],chains[n][2],chains[n][3],chains[n][4]
        resseq = [i[:3] for i in resn]

        # first, check if cterminal is missing a backbone atom and terminate the program
        if not  ('N' in atoms[-1] and 'CA' in atoms[-1] and 'C' in atoms[-1]):
            missing_backbone_atom(resNo[-1],n+1)

        # now loop through all the chains, checking and fixing missing heavy  atoms
        for k in range(len(atoms)):

            try:
                x = list(param[resseq[k]].keys())[:-1]
            except:
                if any(resseq[k].startswith(i) for i in ['HSD','HSP','HSE','HIE']):
                    resseq[k] = 'HIS' # replace with HIS
                    x = list(param[resseq[k]].keys())[:-1]
                else:
                    non_standard_residue(resNo[k],resseq[k],n+1) # non-standard residues not allowed


            if resseq[k] == 'ASP':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0], atmpos[k-1][1],atmpos[k-1][2], atmpos[k][0]]
                    atoms[k],atmpos[k]=repairAsp(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'ASN':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0], atmpos[k-1][1],atmpos[k-1][2], atmpos[k][0]]
                    atoms[k],atmpos[k]=repairAsn(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'ALA':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0],atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                    atoms[k],atmpos[k]=repairAla(atoms[k],atmpos[k],atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'ARG':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0],atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                    atoms[k],atmpos[k]=repairArg(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'LEU':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0],atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                    atoms[k],atmpos[k]=repairLeu(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'ILE':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):

                    if 'CD1' in missing and 'CD' in atoms[k]: # CD1 is named CD by GROMACS, CHARMM, etc.
                        del missing[missing.index('CD1')]
                        atoms[k][atoms[k].index('CD')] = 'CD1'

                    if len(missing) != 0:
                        psi = [atmpos[k-1][0],atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                        atoms[k],atmpos[k]=repairIle(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                        repairNotes(resNo[k], n+1,missing)
                    else:
                        pass

                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'PRO':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0],atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                    atoms[k],atmpos[k]=repairPro(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'SER':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0],atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                    atoms[k],atmpos[k]=repairSer(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'VAL':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0],atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                    atoms[k],atmpos[k]=repairVal(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'LYS':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0],atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                    atoms[k],atmpos[k]=repairLys(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'TYR':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0],atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                    atoms[k],atmpos[k]=repairTyr(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'TRP':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0],atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                    atoms[k],atmpos[k]=repairTrp(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'MET':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0], atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                    atoms[k],atmpos[k]=repairMet(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'HIS':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0],atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                    atoms[k],atmpos[k]=repairHis(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'GLY':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0], atmpos[k-1][1],atmpos[k-1][2], atmpos[k][0]]
                    atoms[k],atmpos[k]=repairGly(atoms[k],atmpos[k],atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'CYS':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0],atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                    atoms[k],atmpos[k]=repairCys(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'GLU':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0],atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                    atoms[k],atmpos[k]=repairGlu(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'PHE':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0],atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                    atoms[k],atmpos[k]=repairPhe(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'GLN':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0],atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                    atoms[k],atmpos[k]=repairGln(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)


            elif resseq[k] == 'THR':
                missing=[i for i in x if not i in atoms[k]]
                if not missing:
                    pass
                elif missing and not ('N' in missing or 'CA' in missing or 'C' in missing):
                    psi = [atmpos[k-1][0],atmpos[k-1][1],atmpos[k-1][2],atmpos[k][0]]
                    atoms[k],atmpos[k]=repairThr(atoms[k],atmpos[k],resNo[k],n,missing,faspr,atmpos[(k+1)%len(atmpos)],psi)
                    repairNotes(resNo[k], n+1,missing)
                elif ('N' in missing or 'CA' in missing or 'C' in missing):
                    missing_backbone_atom(resNo[k],n+1)

        # add C-terminal O [OXT] if it is missing
        if not 'OXT' in atoms[-1]:
            atmpos[-1].extend(c_termini(atoms[-1],atmpos[-1]))
            atoms[-1].extend(['OXT'])
            repairNotes(resNo[-1], n+1,['OXT'])

        # get all CYS SG coordinates for disulfide bond checks.
        # notice that this is done when all missing atoms are fixed.
        atom_pos =  list(zip(atoms,atmpos));sg_coord = []
        for i in range(len(atom_pos)):
            for j,k in enumerate(atom_pos[i][0]):
                if k == 'SG':
                    sg_coord.append(atom_pos[i][1][j])

        ter_all.append([resn,atoms,atmpos,resNo,resseq,sg_coord,chain_id])

    return  (ter_all if not keep_ligand else  [ter_all,ligand])
