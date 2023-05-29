#!/usr/bin/env python

__doc__ = """

This program requires python 3.6 or higher.

This module, ReadStructure.py is used to

read .pdb or .cif files written by Pras_Server.

The purpose is to return data used to assign

secondary structure or draw Ramanchandran plots

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

import numpy as np
from numpy import cross, dot
from math import sqrt, acos,  degrees

def calcTorsionAngle(coord1, coord2, coord3, coord4):
    """
    Given 4 known atom coordinates this function
    calculates the dihedral angle defining the atoms

    Arguments
    ----------
    coord1: the first 3D point
    coord2: the second 3D point
    coord3: the third 3D point
    coord4: the fourth 3D point

    Returns
    -------
    Float: the calculated angle in degrees
    """

    coord1,coord2,coord3,coord4 =\
    np.array(coord1), np.array(coord2),\
    np.array(coord3),np.array(coord4)

    bvec12 = coord1 - coord2
    bvec32 = coord3 - coord2
    bvec43 = coord4 - coord3

    pvec13 = cross(bvec12, bvec32)
    pvec24 = cross(bvec43, bvec32)

    projec = dot(pvec13, pvec24)
    sqrd13 = dot(pvec13, pvec13)
    sqrd24 = dot(pvec24, pvec24)

    cosine = projec / sqrt(sqrd13*sqrd24)
    cosine = min(1.0, max(-1.0, cosine))
    radian = acos(cosine)

    if dot(pvec13, cross(pvec24, bvec32)) < 0:
        radian = -radian

    return degrees(radian)

def calcDistance (m, n):
    """
    This function calculates the distance between two points

    Arguments
    ----------
    m:   the first 3D point
    n:   the second 3D point

    Returns
    -------
    Float: the calculated distance
    """

    return (sqrt(sum([(m[z] - n[z])**2 for z in range(3)])))

def alpha_helix(chain_i, _format):
    """
    This function generates measurements required to assign alpha-helix

    Arguments
    ----------
    chain_i: a list that represents the entire single chain of a PDB file

    _format: the file format

    Returns
    -------
    A list: the calculated measurements
    """

    O_list = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "O"] if _format ==\
             '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "O"]

    N_list = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "N"] if _format ==\
             '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "N"]

    # get atom coordinates depending on the input tile format
    O_coord = [[float(O[30+8*i:38+8*i]) for i in range(3)] for O in O_list] if _format ==\
             '.pdb' else [[float(O[10+i]) for i in range(3)] for O in O_list]

    N_coord = [[float(N[30+8*i:38+8*i]) for i in range(3)] for N in N_list] if _format ==\
              '.pdb' else [[float(N[10+i]) for i in range(3)] for N in N_list]

    O_coord = [np.array(i) for i in O_coord]
    N_coord = [np.array(i) for i in N_coord]

    del O_coord[(len(N_list)-4):]
    del N_coord[:4]

    hbond_distance = []
    for k in range(len(N_coord)):
    	distance = calcDistance(O_coord[k],N_coord[k])
    	hbond_distance.append(distance)

    last_four_res = [100, 100, 100, 100]
    hbond_distance.extend(last_four_res)

    return hbond_distance

def _310helix(chain_i, _format):
    """
    This function generates measurements required to assign 310-helix

    Arguments
    ----------
    chain_i: a list that represents the entire single chain of a PDB file

    _format: the file format

    Returns
    -------
    A list: the calculated measurements
    """

    O_list = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "O"] if _format ==\
             '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "O"]

    N_list = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "N"] if _format ==\
             '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "N"]

    # get atom coordinates depending on the input tile format
    O_coord = [[float(O[30+8*i:38+8*i]) for i in range(3)] for O in O_list] if _format ==\
             '.pdb' else [[float(O[10+i]) for i in range(3)] for O in O_list]

    N_coord = [[float(N[30+8*i:38+8*i]) for i in range(3)] for N in N_list] if _format ==\
              '.pdb' else [[float(N[10+i]) for i in range(3)] for N in N_list]

    O_coord = [np.array(i) for i in O_coord]
    N_coord = [np.array(i) for i in N_coord]

    del O_coord[(len(N_list)-3):]
    del N_coord[:3]

    hbond_distance = []
    for k in range(len(N_coord)):
    	distance = calcDistance(O_coord[k],N_coord[k])
    	hbond_distance.append(distance)

    last_three_res = [100, 100, 100]
    hbond_distance.extend(last_three_res)

    return hbond_distance

def pi_helix(chain_i, _format):
    """
    This function generates measurements required to assign pi-helix

    Arguments
    ----------
    chain_i: a list that represents the entire single chain of a PDB file

    _format: the file format

    Returns
    -------
    A list: the calculated measurements
    """

    O_list = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "O"] if _format ==\
             '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "O"]

    N_list = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "N"] if _format ==\
             '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "N"]

    # get atom coordinates depending on the input tile format
    O_coord = [[float(O[30+8*i:38+8*i]) for i in range(3)] for O in O_list] if _format ==\
             '.pdb' else [[float(O[10+i]) for i in range(3)] for O in O_list]

    N_coord = [[float(N[30+8*i:38+8*i]) for i in range(3)] for N in N_list] if _format ==\
              '.pdb' else [[float(N[10+i]) for i in range(3)] for N in N_list]

    O_coord = [np.array(i) for i in O_coord]
    N_coord = [np.array(i) for i in N_coord]

    del O_coord[(len(N_list)-5):]
    del N_coord[:5]

    hbond_distance = []

    for k in range(len(N_coord)):
    	distance = calcDistance(O_coord[k],N_coord[k])
    	hbond_distance.append(distance)

    last_five_res = [100, 100, 100, 100, 100]
    hbond_distance.extend(last_five_res)

    return hbond_distance

def beta_turn(chain_i, _format):
    """
    This function generates measurements required to assign beta-turn

    Arguments
    ----------
    chain_i: a list that represents the entire single chain of a PDB file

    _format: the file format

    Returns
    -------
    A list: the calculated measurements
    """

    C_list = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "C"] if _format ==\
             '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "C"]

    N_list = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "N"] if _format ==\
             '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "N"]

    # get atom coordinates depending on the input tile format
    C_coord = [[float(C[30+8*i:38+8*i]) for i in range(3)] for C in C_list] if _format ==\
              '.pdb' else [[float(C[10+i]) for i in range(3)] for C in C_list]

    N_coord = [[float(N[30+8*i:38+8*i]) for i in range(3)] for N in N_list] if _format ==\
              '.pdb' else [[float(N[10+i]) for i in range(3)] for N in N_list]

    C_coord = [np.array(i) for i in C_coord]
    N_coord = [np.array(i) for i in N_coord]

    del C_coord[(len(N_list)-3):]
    del N_coord[:3]

    hbond_distance = []

    for k in range(len(N_coord)):
    	distance = calcDistance(C_coord[k],N_coord[k])
    	hbond_distance.append(distance)

    last_three_res = [100, 100, 100]
    hbond_distance.extend(last_three_res)

    return hbond_distance

def amideDihedral(chain_i, _format):
    """
    This function generates amide dihedral angle required to assign Sec. Struc.

    Arguments
    ----------
    chain_i: a list that represents the entire single chain of a PDB file

    _format: the file format

    Returns
    -------
    A list: the calculated dihedral angle
    """

    O_list  = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "O"] if _format ==\
              '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "O"]

    O_list2 = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "O"] if _format ==\
              '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "O"]

    C_list  = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "C"] if _format ==\
              '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "C"]

    C_list2 = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "C"] if _format ==\
              '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "C"]

    # get atom coordinates depending on the input tile format
    O_coord  = [[float(O[30+8*i:38+8*i]) for i in range(3)] for O in O_list]  if _format ==\
               '.pdb' else [[float(O[10+i]) for i in range(3)] for O in O_list]

    O_coord2 = [[float(O[30+8*i:38+8*i]) for i in range(3)] for O in O_list2] if _format ==\
               '.pdb' else [[float(O[10+i]) for i in range(3)] for O in O_list2]

    C_coord  = [[float(C[30+8*i:38+8*i]) for i in range(3)] for C in C_list]  if _format ==\
               '.pdb' else [[float(C[10+i]) for i in range(3)] for C in C_list]

    C_coord2 = [[float(C[30+8*i:38+8*i]) for i in range(3)] for C in C_list2] if _format ==\
               '.pdb' else [[float(C[10+i]) for i in range(3)] for C in C_list2]

    del O_coord2[0]
    del O_coord[-1]
    del C_coord2[0]
    del C_coord[-1]

    amide_angles = [360]

    for k in range(len(O_coord)):
    	angle = calcTorsionAngle(O_coord[k],C_coord[k], C_coord2[k],O_coord2[k])
    	amide_angles.append(angle)

    return amide_angles

def psiDihedral(chain_i, _format):
    """
    This function generates psi dihedral angle required to assign Sec. Struc.

    Arguments
    ----------
    chain_i: a list that represents the entire single chain of a PDB file

    _format: the file format

    Returns
    -------
    A list: the calculated dihedral angle
    """

    N_list   = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "N"] if _format ==\
              '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "N"]

    N_list2  = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "N"] if _format ==\
              '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "N"]

    Ca_list  = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "CA"] if _format ==\
              '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "CA"]

    C_list   = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "C"] if _format ==\
              '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "C"]

    # get atom coordinates depending on the input tile format
    N_coord  = [[float(N[30+8*i:38+8*i]) for i in range(3)] for N in N_list] if _format ==\
               '.pdb' else [[float(N[10+i]) for i in range(3)] for N in N_list]

    N_coord2 = [[float(N[30+8*i:38+8*i]) for i in range(3)] for N in N_list2] if _format ==\
               '.pdb' else [[float(N[10+i]) for i in range(3)] for N in N_list2]

    Ca_coord = [[float(Ca[30+8*i:38+8*i]) for i in range(3)] for Ca in Ca_list] if _format ==\
               '.pdb' else [[float(Ca[10+i]) for i in range(3)] for Ca in Ca_list]

    C_coord  = [[float(C[30+8*i:38+8*i]) for i in range(3)] for C in C_list] if _format ==\
               '.pdb' else [[float(C[10+i]) for i in range(3)] for C in C_list]

    del N_coord2[0]
    del N_coord[-1]
    del Ca_coord[-1]
    del C_coord[-1]

    psi_angles = []
    for k in range(len(N_coord)):
    	angle = calcTorsionAngle(N_coord[k],Ca_coord[k],C_coord[k], N_coord2[k])
    	psi_angles.append(angle)

    psi_angles.append(360)

    return psi_angles

def phiDihedral(chain_i, _format):
    """
    This function generates phi dihedral angle required to assign Sec. Struc.

    Arguments
    ----------
    chain_i: a list that represents the entire single chain of a PDB file

    _format: the file format

    Returns
    -------
    A list: the calculated dihedral angle

    """

    chain_i = chain_i[2:]

    C_list   = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "C"] if _format ==\
              '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "C"]

    N_list   = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "N"] if _format ==\
              '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "N"]

    Ca_list  = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "CA"] if _format ==\
              '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "CA"]

    C_list2  = [l for l in chain_i if l[0:6] == "ATOM  " and l[13:16].strip(' ') == "C"] if _format ==\
              '.pdb' else [l.split() for l in chain_i if l.split()[0] == "ATOM" and l.split()[3] == "C"]

    # get atom coordinates depending on the input tile format
    C_coord  = [[float(C[30+8*i:38+8*i]) for i in range(3)] for C in C_list] if _format ==\
               '.pdb' else [[float(C[10+i]) for i in range(3)] for C in C_list]

    N_coord  = [[float(N[30+8*i:38+8*i]) for i in range(3)] for N in N_list] if _format ==\
               '.pdb' else [[float(N[10+i]) for i in range(3)] for N in N_list]

    Ca_coord = [[float(Ca[30+8*i:38+8*i]) for i in range(3)] for Ca in Ca_list] if _format ==\
               '.pdb' else [[float(Ca[10+i]) for i in range(3)] for Ca in Ca_list]

    C_coord2 = [[float(C[30+8*i:38+8*i]) for i in range(3)] for C in C_list2] if _format ==\
               '.pdb' else [[float(C[10+i]) for i in range(3)] for C in C_list2]

    del C_coord[-1]
    del C_coord2[0]

    phi_angles = [360]

    for k in range(len(C_coord)):
    	angle = calcTorsionAngle(C_coord[k], N_coord[k], Ca_coord[k], C_coord2[k])
    	phi_angles.append(angle)

    return phi_angles

def resChain(chain_i, _format):
    """
    This function generates the chain ID

    Arguments
    ----------
    chain_i: a list that represents the entire single chain of a PDB file

    _format: the file format

    Returns
    -------
    A list: the chain ID
    """

    chain_list   =  []

    if _format == '.pdb':
        for lines in chain_i:
            if lines[13:16].rstrip() == 'C':
                chain   = lines[21:22]
                chain_list.append(chain)
    else:
        for lines in chain_i:
            if lines.split()[3] == 'C':
                chain   = lines.split()[6]
                chain_list.append(chain)

    return chain_list

def resName(chain_i, _format):
    """
    This function generates the three letter residue code

    Arguments
    ----------
    chain_i: a list that represents the entire single chain of a PDB file

    _format: the file format

    Returns
    -------
    A list: the three letter code
    """

    resname_list =  []

    if _format == '.pdb':
        for lines in chain_i:
            if lines[13:16].rstrip() == 'C':
                resname = lines[17:20]
                resname_list.append(resname)
    else:
        for lines in chain_i:
            if lines.split()[3] == 'C':
                resname = lines.split()[5]
                resname_list.append(resname)

    return resname_list

def nextRes(chain_i, _format):
    """
    This function generates the Ramachandran type

    Arguments
    ----------
    chain_i: a list that represents the entire single chain of a PDB file

    _format: the file format

    Returns
    -------
    A list: the Ramachandran type
    """

    resn =  []

    if _format == '.pdb':
        for lines in chain_i:
            if lines[13:16].rstrip() == 'C':
                restype = lines[17:20]
                resn.append(restype)
    else:
        for lines in chain_i:
            if lines.split()[3] == 'C':
                restype = lines.split()[5]
                resn.append(restype)

    restype_list = []

    for i, restype in enumerate(resn):
        if resn[i] == 'PRO':
            restype = 'Proline'
            restype_list.append(restype)
        elif resn[i] == 'GLY':
            restype = 'Glycine'
            restype_list.append(restype)
        elif resn[(i+1)% len(resn)] =='PRO' and i != (len(resn)-1):
            restype = 'Pre-Pro'
            restype_list.append(restype)
        else:
            restype = 'General'
            restype_list.append(restype)

    return restype_list

def resNumber(chain_i, _format):
    """
    This function generates the residue number

    Arguments
    ----------
    chain_i: a list that represents the entire single chain of a PDB file

    _format: the file format

    Returns
    -------
    A list: the residue number
    """

    resnum_list  =  []

    if _format == '.pdb':
        for lines in chain_i:
            if lines[13:16].rstrip() == 'C':
                resnum  = lines[23:26]
                resnum_list.append(resnum)

    else:
        for lines in chain_i:
            if lines.split()[3] == 'C':
                resnum  = lines.split()[8]
                resnum_list.append(resnum)

    return resnum_list
