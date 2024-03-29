#!/usr/bin/env python3

__doc__ = """

This program requires python 3.6 or higher.

This module, PRAS.py is the main module that

has the function that is used to fix both

missing hydrogen and heavy atoms.

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
import itertools
from .MissingHydrogenAtoms import *
from .CheckPDBatoms import checkpdbAtoms

# The dictionary below points to the data for ATOM line.
# A given category and attribute is also known as an
# mmCIF token. Values or attributes with comments are
# additional information from the author. They are not
# always provided.

atom_site = {
    "_atom_site": [
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_alt_id",
        "label_comp_id",
        "label_asym_id",
        "label_entity_id",
        "label_seq_id",
        "pdbx_PDB_ins_code",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "B_iso_or_equiv",
        "pdbx_formal_charge",
        "auth_seq_id",  # an alternative identifier for _atom_site.label_seq_id
        "auth_comp_id", # an alternative identifier for _atom_site.label_comp_id
        "auth_asym_id", # an alternative identifier for _atom_site.label_asym_id
        "auth_atom_id", # an alternative identifier for _atom_site.label_atom_id
        "pdbx_PDB_model_num",
    ]
}


# Below is the atom_type.symbol category
# '#'and 'loop_' show end and beginning
# of another category, respectively.

symbol = ['C','H','N','O','S','#','loop_']

def writeHeader(_format, ofname):
    """
    This function writes the appropriate header
    for a given category in tabular style. For
    a .pdb file only a comment is written

    Arguments
    ----------
    _format: the file format (.pdb or .cif)

    ofname     : the output file name, default name is out.pdb or out.cif

    Returns
    -------
    None:    it writes to a .cif or .pdb file
    """

    if _format == 'cif':

        with open(ofname, 'a')  as f:

            f.write('# structure generated by PRAS'+'\n'+'#'
                    +'\n'+'data_xxxx'+'\n'+'#'+'\n'+
                    '_entry.id xxxx'+'\n'+'#'+'\n'+'loop_'+'\n')

            f.write('_atom_type.symbol'+'\n')

            for i in symbol:
                f.write(i+'\n')

            for i in atom_site["_atom_site"]:
                f.write('_atom_site.'+i+'\n')
    else:
        with open(ofname, 'a')  as f:
            f.write('# structure generated by PRAS'+'\n')

def atomName(element):
    """
    This function returns the
    element symbol corresponding
    to the atom type

    Arguments
    ----------
    name:   the standard PDB atom name

    Returns
    -------
    string: the element name of the atom
            type for writing to a pdb or
            cif file
    """

    element = [i for i in element if i.isalpha()]

    return element[0]

def atomType(name):
    """
    This function returns the
    format for writing the atom
    name to a pdb file or cif file.

    For cif, the format is not essential
    because data values have no specific
    column width as in pdb file

    Arguments
    ----------
    name:   the standard PDB atom name

    Returns
    -------
    string: the atom name and format
            for writing to a PDB file
    """

    if len(name) == 1:
        name = " %s  " % name
    elif len(name) == 2:
        name = " %s " % name
    elif len(name) == 3:
        name = " %s" % name
    elif len(name) == 4:
        name = "%s" % name
    return name

def insertRes(res, _format = None):
    """
    This function returns the letter for
    residue insertion to be written to a
    PDB file (.pdb of .cif). In the case
    of .cif the code is usually "?" except
    it exists.

    Arguments
    ----------
    res:    the res name + res No + insertion
            as concatenated by PRAS when reading
            input PDB file

    Returns
    -------
    string: the insertion code/letter if it
            exists else a blank space (for .pdb)
            or "?" for .cif
    """

    if not _format:
        if res[-1].isalpha():
            return res[-1]
        else:
            return " "
    else:
        if res[-1].isalpha():
            return res[-1]
        else:
            return "?"


def repairPDB(pdb_pras,rotamer,mutation,pdb_faspr,keep_ligand,chain_no, ofname=False,his_p=False):
    """
    This function steps over each residue of a chain and obtains from
    MissingHydrogenAtoms.py all atoms of the residue to be wrttien to a new PDB file.

    If the residue is not the last and the next residue is not PRO the returned list
    will contain a backbone H otherwise it will have no backbone H.

    For Cys residue, SG_coords are used to check for disulfide bonds.
    If disulfide bond exists it obtains no HG, otherwise it obtains HG

    Arguments (this function takes 6 compulsory arguments)
    ----------
    pdb_pras:    the syntactically correct PDB file to repair/analyze (.pdb or .cif format)

    rotamer:     by default, atoms with the highest occupancy are written to PDB file in
                 the case of rotamers. If you want low occupancy instead, supply "no" or
                 any string as the argument, otherwise supply "" as argument.

    mutation:    when there is point mutation, two residues with unequal
                 occupany (i.e. the residue atoms) are given same residue number.
                 By default the residue with the highest occupancy is
                 written to PDB file. If you want low occupancy instead supply "no"
                 or any string as the argument, otherwise supply "" as argument.

    pdb_faspr:   the PDB file obtained from FASPR (by running the same PDB supplied to PRAS).
                 FASPR is a free and open source side-chain packer, thanks to Huang et al.
                 You can obtain it at https://github.com/tommyhuangthu/FASPR

                 Note that if there are no missing atoms that require flexible chi to fix then
                 FASPR is not required. If there are such atoms and you do not supply FASPR PDB file,
                 PRAS will use the default or most probable chi from Dunbrack 2011 rotamer library
                 to fix the atoms. So in any case PRAS will run but will notify you if FASPR PDB is
                 not supplied.

                 Although FASPR is more accurate than most state-of-the-art side-chain packers,
                 it is not infinitely accurate and sometimes the default chi from PRAS is the
                 right conformation for the amino acid residue side-chain. It is advised that
                 you compare both methods when necessary.

                 Be mindful of the fact that FASPR may be less flexible with reading PDB files and
                 has no mechanism to process .cif files or use lower conformers in the case of
                 rotamers/mutation. To avoid manual editing, pass the PDB through PRAS first.
                 Also, for repairing .cif, you can use the equivalent of that in .pdb obtained from
                 FASPR since it cannot process a .cif file. The chain PRAS obtains from .pdb is the
                 same as .cif if both are the same structure, even though it reads .cif differently!

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

    his_p     : If True, 20% of histidines will be protonated (+1 charge). It only matters if there are more than 4 HIS

    Returns
    -------
    None:        repaired and hydrogenated PDB structure will be written.
                 Ensure you have write permission in the directory where you
                 run the program
    """

    # remove previous files written by PRAS b/c
    # it writes in the append mode
    for i in os.listdir(os.getcwd()):
        if i == 'out.pdb' or i == 'out.cif':
            os.remove(i)

    nchains = checkpdbAtoms(pdb_pras,rotamer,mutation,pdb_faspr, keep_ligand, chain_no,ofname)

    if keep_ligand:
        nchains, ligand = nchains

    chains = list(itertools.chain(*[i[-1] for i in nchains]))

    # atom numbering starts from the first chain
    number=1
    for n in range(len(nchains)):
        resn,atom,atmpos,resNo,resseq,sg_coord,chain = nchains[n]
        x = list(zip(resn,atom,atmpos))
        res_pos,h_pos,atom_name = [],[],[]

        # steps over the residues and if the next residue is not PRO
        # it obtains backbone H but if PRO or the last residue in the
        # chain it does not obtain backbone H. The H is later added
        # to residue i+1
        for i,j in  enumerate(x):
            if j[0].startswith('ARG'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name = arginine_h(x[i],i,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name = arginine_h(x[i])
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('ALA'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=alanine_h(x[i],i,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=alanine_h(x[i])
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('ASP'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=aspartate_h(x[i],i,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=aspartate_h(x[i])
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('ASN'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=asparagine_h(x[i],i,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=asparagine_h(x[i])
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('GLU'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=glutamate_h(x[i],i,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=glutamate_h(x[i])
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('GLN'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=glutamine_h(x[i],i,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=glutamine_h(x[i])
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('GLY'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=glycine_h(x[i],i,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=glycine_h(x[i])
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('HIS') or any(j[0].startswith(hs) for hs in ['HSD','HSP','HSE','HIE']):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=histidine_h(x[i],i,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=histidine_h(x[i])
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('ILE'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=isoleucine_h(x[i],i,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=isoleucine_h(x[i])
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('LEU'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=leucine_h(x[i],i,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=leucine_h(x[i])
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('LYS'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=lysine_h(x[i],i,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=lysine_h(x[i])
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('MET'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=metheonine_h(x[i],i,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=metheonine_h(x[i])
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('PRO'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=proline_h(x[i],i,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=proline_h(x[i])
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('TRP'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=tryptophan_h(x[i],i,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=tryptophan_h(x[i])
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('VAL'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=valine_h(x[i],i,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=valine_h(x[i])
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('PHE'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=phenylalanine_h(x[i],i,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=phenylalanine_h(x[i])
                    res_pos.extend([apos]);atom_name.extend([name])

                """
                The remaining residues (SER, THR, TYR and free CYS) have what we regard as
                class6 H-atoms with rotational freedom. These H-atoms are optimized with a
                potential energy function that may affect the computational time. You can go
                to PotentialEnergy.py and look for "optimize()" function and decrease/increase
                the rotation interval. Currently, it is set to 60 degrees but in the online PRAS
                server it is set to a lower value.

                If you are not interested in the optimization, go to MissingHydrogenAtoms.py
                locate SER, THR, TYR and CYS and comment the lines as appropriate
                """
            elif j[0].startswith('SER'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=serine_h(x[i],i,resNo[i],atom,atmpos,resNo,resseq,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=serine_h(x[i],i,resNo[i],atom,atmpos,resNo,resseq)
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('THR'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=threonine_h(x[i],i,resNo[i],atom,atmpos,resNo,resseq,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=threonine_h(x[i],i,resNo[i],atom,atmpos,resNo,resseq)
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('TYR'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    apos,hpos,name=tyrosine_h(x[i],i,resNo[i],atom,atmpos,resNo,resseq,x[(i+1)%len(x)])
                    res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    apos,name=tyrosine_h(x[i],i,resNo[i],atom,atmpos,resNo,resseq)
                    res_pos.extend([apos]);atom_name.extend([name])

            elif j[0].startswith('CYS'):
                if x[(i+1)% len(x)][0][:3] != 'PRO' and x[i] != x[-1]:
                    if checkDisulfide(x[i][2][5],sg_coord) == 'is_bond':
                        apos,hpos,name=isDisulfide(x[i],i,x[(i+1)%len(x)])
                        res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                    else:
                        apos,hpos,name=notDisulfide(x[i],i,resNo[i],atom,atmpos,resNo,resseq,x[(i+1)%len(x)])
                        res_pos.extend([apos]);h_pos.extend([hpos]);atom_name.extend([name])
                else:
                    if checkDisulfide(x[i][2][5],sg_coord) == 'is_bond':
                        apos,name=isDisulfide(x[i],i)
                        res_pos.extend([apos]);atom_name.extend([name])
                    else:
                        apos,name=notDisulfide(x[i],i,resNo[i],atom,atmpos,resNo,resseq)
                        res_pos.extend([apos])
                        atom_name.extend([name])

        """
        Next code line randomly protonates 20% of all HIS residues
        See explanation online at www.protein-science.com
        """
        if his_p: res_pos, atom_name = prot_his(x,res_pos,atom_name,n)

        # Add backbone hydrogen.
        # It belongs to resi i+1
        for i,j in enumerate(h_pos):
            res_pos[h_pos[i][1]+1].extend([h_pos[i][0]])
            atom_name[h_pos[i][1]+1].extend('H')

        # Add n-terminal hydrogen, if proline only two H
        # If not PRO, three H
        if x[0][0][:3] == 'PRO':
            res_pos[0].extend(ntermini_pro(x[0],resn[0][:3]))
            atom_name[0].extend(['H1','H2'])
        else:
            res_pos[0].extend(ntermini_notpro(x[0],resn[0][:3]))
            atom_name[0].extend(['H1','H2','H3'])

        # chek if user specified cif format
        if pdb_pras[-4:] == '.cif':
            if not ofname:
                ofname = "out.cif"

            # steps over each residue of .cif file,
            # then steps over each residue atom,
            # and writes to a PDB file in .cif format
            with open(ofname, 'a') as f:
                if n == 0:
                    # write a header for the .cif output file
                    writeHeader('cif',ofname)

                for k,l in enumerate (res_pos):
                    atm, charge, model_No = atom_name[k],'0','1' # last two values are cosmetics.
                    for i,j in enumerate(l):
                        f.write("%6s %-4s %s %s %s %3s %1s %s %-3s %s %8.3f %8.3f %8.3f %6.2f %6.2f %s %s %s %s %s %s"\
                        %('ATOM  ',str(number), atomName(atm[i]), atomType(atm[i]), '.', resseq[k],\
                        chains[n], str(n+1), str(resNo[k]),insertRes(resn[k], 'cif'),\
                        l[i][0], l[i][1], l[i][2], 1.00, 0.00, charge, str(resNo[k]), resseq[k], chains[n],\
                        atomType(atm[i]), model_No)+'\n')
                        number+=1

                if n == len(nchains)-1:
                    # ligand is written at the end which is
                    # the case with most .cif files, except
                    # it is covalently bonded to the protein.
                    # If the HETATM is bonded covalently to the
                    # protein, user should not opt to keep
                    # the ligand as PRAS is not meant to
                    # repair ligands.
                    if (keep_ligand and ligand):
                        for line in ligand:
                            f.write(line+'\n')
                    f.write("#")
                    f.close()

        else:
            # steps over each residue of .pdb file,
            # then steps over each residue atom
            # and writes to a PDB file in .pdb format
            if not ofname:
                ofname = "out.pdb"
            with open(ofname, 'a')  as f:
                if n == 0:
                    # write a comment for the .pdb output file
                    writeHeader('pdb',ofname)

                for k,l in enumerate (res_pos):
                    atm = atom_name[k]
                    for i,j in enumerate(l):
                        f.write("%6s%5s %4s %-4s%1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f           %1s" \
                        %('ATOM  ',str(number)[-5:],atomType(atm[i]),resseq[k],\
                        chains[n],str(resNo[k])[-4:],insertRes(resn[k]),\
                        l[i][0],l[i][1],l[i][2],1.00,0.00, atomName(atm[i]))+'\n')
                        number+=1

                # if there are more than 1 chain write,
                # TER at the end and if the last chain
                # write TER and then END
                if n != len(nchains)-1:
                    f.write('{:4s}' '{:2s}' '{:5d}''{:s}' '{:4s}''{:s}''{:3s}' '{:2s}''{:4d}\n'.format\
                        ('TER',' ',number,' ',' ', ' ',resseq[len(x)-1],' '+chains[n],resNo[len(x)-1]))
                    number+=1
                else:
                    f.write('{:4s}' '{:2s}' '{:5d}''{:s}' '{:4s}''{:s}''{:3s}' '{:2s}''{:4d}\n'.format\
                        ('TER',' ',number,' ',' ', ' ',resseq[len(x)-1],' '+chains[n],resNo[len(x)-1]))

                    # for .pdb, the ligand is always written after TER except for a
                    # HETATM that is covalently bonded to the protein. Obviously, anything
                    # other than the 20 common amino acids that is bonded to the protein
                    # covalently (e.g., O-SULFO-L-TYROSINE) is ignored by PRAS if user
                    # did not opt to keep ligands, otherwise the chemistry of the added
                    # H-atoms will be flawed. That is, such HETATM would clash with the
                    # added H-atoms b/c they are not considered when adding H-atoms
                    if (keep_ligand and ligand):
                        for line in ligand:
                            f.write(line+'\n')

                    f.write('END')
                    f.close()
