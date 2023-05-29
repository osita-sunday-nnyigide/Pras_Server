#!/usr/bin/env python

__doc__ = """

This program requires python 3.6 or higher.

This module, MissingHeavyAtoms.py has all

the functions required to fix missing heavy

atoms. These functions represent the 20 common

amino acid residues and other extra functions.

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
from math import sqrt, acos,degrees, pi

def calcDistance(m, n):

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

def calcCoordinate(a, b, c, bond_len, theta, di_angle):

    """
    Given 3 known atom coordinates this function
    obtains the 4th atom coordinate that satisfies
    the input geometry constraints.

    Arguments
    ----------
    a:        the first 3D point
    b:        the second 3D point
    c:        the third 3D point
    bond_len: the bond length of c and the 4th atom
    theta:    the angle formed by b, c & the 4th atom
    di_angle: the dihedral angle defining the 4 atoms

    Returns
    -------
    A list: the calculated 4th atom coordinate
    """

    a,b,c    = np.array(a),np.array(b),np.array(c)
    di_angle = di_angle*pi/180.0

    u = c - b
    x = a - b

    v = ((x) - (np.dot((x), u)/np.dot(u,u))*u)

    w = np.cross(u, x)

    q = (v/np.linalg.norm(v))*np.cos(di_angle)
    e = (w/np.linalg.norm(w))*np.sin(di_angle)

    pos_temp2 = np.array((b + (q+e)))

    u1 = b - c
    y1 = pos_temp2 - c

    mag_y1  = sqrt((y1[0]**2)+(y1[1]**2)+(y1[2]**2))
    mag_u1  = sqrt((u1[0]**2)+(u1[1]**2)+(u1[2]**2))

    theta_bcd = acos(np.dot(u1,y1)/(mag_u1*mag_y1))
    rotate    = theta - degrees(theta_bcd)

    z  = np.cross(u1, y1)
    n  = z/sqrt((z[0]**2)+(z[1]**2)+(z[2]**2))

    pos_ini = c + y1*np.cos(rotate*pi/180) +\
              (np.cross(n,y1)*np.sin(rotate*pi/180))\
              + n*(np.dot(n,y1))*(1-np.cos(rotate*pi/180))

    pos_BL = ((pos_ini-c)*(bond_len/np.linalg.norm(pos_ini-c)))+c

    return pos_BL.tolist()

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

    coord1,coord2,coord3,coord4 = np.array(coord1), \
    np.array(coord2),np.array(coord3),np.array(coord4)

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

def get_chi(atom_1, atom_2, atom_3, atom_4, resNum, chainNum,faspr):

    """
    This function calculates flexible dihedral angle (chi)
    if a PDB output obtained from FASPR is supplied, if not
    this function is never used!

    Arguments
    ----------
    atom_1:   the first atom name
    atom_2:   the second atom name
    atom_3:   the third atom name
    atom_4:   the fourth atom name
    resNum:   the residue number
    chainNum: the chain number
    faspr   : the chain as obtained from ReadMaster.py

    Returns
    -------
    Float: the calculated angle in degrees
    """

    out_list = faspr
    for i in range(len(out_list)):
        if i == chainNum:
            resn,atom,atmpos,resNo,chain = out_list[i]
            for i in range(len(resNo)):
                if resNo[i] == resNum:
                    pos1 = atmpos[i][atom[i].index(atom_1)]
                    pos2 = atmpos[i][atom[i].index(atom_2)]
                    pos3 = atmpos[i][atom[i].index(atom_3)]
                    pos4 = atmpos[i][atom[i].index(atom_4)]

                    return calcTorsionAngle(pos1,pos2,pos3,pos4)

def addBackbone(N, CA, C, atmpos, psi):

    """
    This function calculates the backbone O coordinate
    and then checks for van der Waals conflict. If there is conflict,
    it uses next favourable di angle until VDW conflict is resolved

    Arguments
    ----------
    N:      the N atom coordinate
    CA:     the CA atom coordinate
    C:      the C atom coordinate
    atmpos: the amino acid that may clash with O
    psi   : the coordinates of atoms required to calculate psi

    Returns
    -------
    A list: the calculated 3D coordinate of atom O
    """

    psi=calcTorsionAngle(psi[0],psi[1],psi[2],psi[3])

    if psi >= 100:
        di = [-30,140,-140,40]
        for i in di:
            O = calcCoordinate(N,CA,C,1.23,120.5,i)
            if calcDistance(O,atmpos[0]) > 2.0:
                break

    elif psi > 0 and psi < 100:
        di = [-140,140,-30,40]
        for i in di:
            O = calcCoordinate(N,CA,C,1.23,120.5,i)
            if calcDistance(O,atmpos[0]) > 2.0:
                break
    else:
        di = [140,-140,-30,40]
        for i in di:
            O = calcCoordinate(N,CA,C,1.23,120.5,i)
            if calcDistance(O,atmpos[0]) > 2.0:
                break

    return O

def c_termini(atom,pos):

    """
    This function will generate C-ter oxygen if it is missing.

    Arguments
    ----------
    atom: a list containing the atom names of the residue
    pos:  a list containing the atom coordinates

    Returns
    -------
    A list: a list containing the coordinates of the c-ter O
    """

    atoms      = atom
    coord      = pos

    for j,k in enumerate(atoms):

        if atoms[j] == 'C':

             pos_C = coord[j]

        elif atoms[j] == 'CA':

            pos_CA = coord[j]

        elif atoms[j] == 'N':

            pos_N = coord[j]

        elif atoms[j] == 'O':
            pos_O = coord[j]

    dih = calcTorsionAngle(pos_N,pos_CA,pos_C,pos_O)

    OXT   = calcCoordinate(pos_N,pos_CA,pos_C,1.25,122.5,180+dih)

    return [OXT]

def repairSer(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function plus the rest represent 20 amino acid residues.
    Annotation is provided in detail here only  b/c all the
    functions are essentially the same, that is, they take
    similar arguments except Ala and Gly that have no side-chain.

    The first priority of PRAS is to keep resolved atoms and replace
    only the missing ones. Each atom that requires flexible chi to fix
    is indicated. If such atom is missing, all other atoms are fixed
    using chi supplied by FASPR or default PRAS chi taken from literature
    if FASPR PDB is not provided. If atoms that are missing do not require
    flexible chi, FASPR is not needed. All constraints (bond/di angle, bond
    length are taken from Dunbrack 2011 rotamer library.

    Arguments
    ----------
    atoms:   a list containing atom name for the residue
    coord:   a list of coordinates of all atoms for the residue
    resNo:   the residue number of the residue
    chain:   the chain number the residue belongs to
    missing: the list of all missing atoms
    faspr:   the chain of FASPR output PDB as obtained from ReadMaster.py
    nextres: the coordinates of the next residue in the chain
    psi   :  the coordinates of atoms required to calculate psi

    Returns
    -------
    A list of lists: index 0 is the atom name and index 1 the atom coordinates
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.66)

    if faspr:

        try:
            N_CA_CB_OG_diangle = get_chi('N', 'CA', 'CB', 'OG',resNo, chain,faspr)

        except:

            N_CA_CB_OG_diangle = []

    """
    We think that PRAS is more flexible than faspr
    when it comes to reading PDB files (especially
    in regards to generating different conformers
    when there are point mutations and/or rotamers.

    So, in rare cases where the PDB file contains
    point mutation and user opts for residue with
    lower occupancy, faspr-named atoms  may not
    correspond to PRAS. Therefore, this maneuver
    ensures that if faspr is supplied and PRAS fail to
    get required chi using the regular atom names,
    a default chi will be used.

    """
    if faspr and not N_CA_CB_OG_diangle:

        N_CA_CB_OG_diangle = -63.3

    if not faspr:

        N_CA_CB_OG_diangle = -63.3

    # for SER, fixing OG, requires flexible chi
    require_chi = ['OG']

    for i in missing:

        if i in require_chi:

            OG = calcCoordinate(N, CA, CB, 1.417, 110.773, N_CA_CB_OG_diangle)

            atompos = [N, CA, C, O, CB, OG]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'OG']

            return [atoms, atompos]


    OG = coord[atoms.index('OG')]

    atompos = [N, CA, C, O, CB, OG]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'OG']

    return [atoms, atompos]

def repairAla(atoms,coord,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the ALA residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O   = coord[atoms.index('O')]
    except:
        O   = addBackbone(N,CA,C,nextres,psi)

    try:
        CB  = coord[atoms.index('CB')]
    except:
        CB  = calcCoordinate(N,C,CA, 1.52,109.5,122.69)

    atompos = [N, CA, C, O, CB]

    atoms   = ['N', 'CA', 'C', 'O', 'CB']

    return [atoms, atompos]

def repairGly(atoms,coord,nextres='',psi=''):

    """
    This function adds missing [O]
    atom of the GLY residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O   = coord[atoms.index('O')]
    except:
        O   = addBackbone(N,CA,C,nextres,psi)

    atompos = [N, CA, C, O]

    atoms   = ['N', 'CA', 'C', 'O']

    return [atoms, atompos]

def repairAsn(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the ASN residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O  = coord[atoms.index('O')]
    except:
        O  = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52,109.50, 123.23)

    if faspr:

        try:
            CA_CB_CG_OD1_diangle = get_chi('CA', 'CB', 'CG','OD1', resNo, chain,faspr)

            N_CA_CB_CG_diangle   = get_chi('N', 'CA', 'CB','CG', resNo, chain,faspr)

        except:

            CA_CB_CG_OD1_diangle, N_CA_CB_CG_diangle = [],[]

    if  faspr and not all([CA_CB_CG_OD1_diangle, N_CA_CB_CG_diangle]):

        CA_CB_CG_OD1_diangle = -58.3

        N_CA_CB_CG_diangle   = -65.5

    if not faspr:

        CA_CB_CG_OD1_diangle = -58.3

        N_CA_CB_CG_diangle   = -65.5

    CA_CB_CG_ND2_diangle = 180.0 + CA_CB_CG_OD1_diangle

    # for ARG, fixing CG,OD1 require flexible chi
    require_chi = ['CG','OD1']

    for i in missing:

        if i in require_chi:

            CG  = calcCoordinate(N, CA, CB, 1.52,112.62, N_CA_CB_CG_diangle)

            OD1 = calcCoordinate(CA, CB, CG, 1.23,120.85, CA_CB_CG_OD1_diangle)

            ND2 = calcCoordinate(CA, CB, CG, 1.33,116.48, CA_CB_CG_ND2_diangle)

            atompos = [N, CA, C, O, CB, CG, OD1, ND2]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2']

            return [atoms, atompos]

    # this residue has only one atom that does not require flexible chi
    # for execution to arrive here it means only that atom is missing
    # so we just calculate it
    CG = coord[atoms.index('CG')]

    OD1 = coord[atoms.index('OD1')]

    CA_CB_CG_OD1_diangle = calcTorsionAngle(CA,CB,CG,OD1)

    CA_CB_CG_ND2_diangle = 180.0 + CA_CB_CG_OD1_diangle

    ND2 = calcCoordinate(CA, CB, CG, 1.33,116.48, CA_CB_CG_ND2_diangle)

    atompos = [N, CA, C, O, CB, CG, OD1, ND2]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2']

    return [atoms, atompos]

def repairGlu(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the GLU residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.87)

    if faspr:

        try:

            CB_CG_CD_OE1_diangle = get_chi('CB', 'CG', 'CD', 'OE1', resNo, chain,faspr)

            CA_CB_CG_CD_diangle = get_chi('CA', 'CB', 'CG', 'CD', resNo, chain,faspr)

            N_CA_CB_CG_diangle = get_chi('N', 'CA', 'CB', 'CG', resNo, chain,faspr)

        except:

            CB_CG_CD_OE1_diangle, CA_CB_CG_CD_diangle, N_CA_CB_CG_diangle = [],[],[]


    if  faspr and not all([CB_CG_CD_OE1_diangle, CA_CB_CG_CD_diangle, N_CA_CB_CG_diangle]):

        CB_CG_CD_OE1_diangle = -6.2

        CA_CB_CG_CD_diangle = -179.8

        N_CA_CB_CG_diangle = -63.8

    if  not faspr:

        CB_CG_CD_OE1_diangle = -6.2

        CA_CB_CG_CD_diangle = -179.8

        N_CA_CB_CG_diangle = -63.8

    CB_CG_CD_OE2_diangle = 180.0 + CB_CG_CD_OE1_diangle

    # for GLU, fixing CG,CD,OE1 require flexible chi
    require_chi = ['CG','CD','OE1']

    for i in missing:

        if i in require_chi:

            CG = calcCoordinate(N, CA, CB, 1.52, 113.82,N_CA_CB_CG_diangle)

            CD = calcCoordinate(CA, CB, CG,1.52, 119.02,  CA_CB_CG_CD_diangle)

            OE1 = calcCoordinate(CB, CG, CD, 1.25, 119.02, CB_CG_CD_OE1_diangle)

            OE2 = calcCoordinate(CB, CG, CD, 1.25, 118.08,CB_CG_CD_OE2_diangle)

            atompos = [N, CA, C, O, CB, CG, CD, OE1, OE2]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2']

            return [atoms, atompos]


    CG = coord[atoms.index('CG')]

    CD = coord[atoms.index('CD')]

    OE1 = coord[atoms.index('OE1')]

    CB_CG_CD_OE1_diangle = calcTorsionAngle(CB,CG,CD,OE1)

    CB_CG_CD_OE2_diangle = 180.0 + CB_CG_CD_OE1_diangle

    OE2 = calcCoordinate(CB, CG, CD, 1.25, 118.08,CB_CG_CD_OE2_diangle)

    atompos = [N, CA, C, O, CB, CG, CD, OE1, OE2]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2']

    return [atoms, atompos]

def repairAsp(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the ASP residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.82)

    if faspr:

        try:
            CA_CB_CG_OD1_diangle = get_chi('CA', 'CB', 'CG','OD1',resNo, chain,faspr)

            N_CA_CB_CG_diangle = get_chi('N', 'CA', 'CB', 'CG',resNo, chain,faspr)

        except:

            CA_CB_CG_OD1_diangle, N_CA_CB_CG_diangle = [],[]


    if faspr and not all([CA_CB_CG_OD1_diangle, N_CA_CB_CG_diangle]):

        CA_CB_CG_OD1_diangle = -46.7

        N_CA_CB_CG_diangle = -66.4

    if not faspr:

        CA_CB_CG_OD1_diangle = -46.7

        N_CA_CB_CG_diangle = -66.4

    CA_CB_CG_OD2_diangle = 180 + CA_CB_CG_OD1_diangle

    # for ASP, fixing CG,OD1 require flexible chi
    require_chi = ['CG','OD1']

    for i in missing:

        if i in require_chi:

            CG = calcCoordinate(N, CA, CB, 1.52, 113.06, N_CA_CB_CG_diangle)

            OD1 = calcCoordinate(CA, CB, CG, 1.25, 119.22, CA_CB_CG_OD1_diangle)

            OD2 = calcCoordinate(CA, CB, CG, 1.25, 118.21, CA_CB_CG_OD2_diangle)

            atompos = [N, CA, C, O, CB, CG, OD1, OD2]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2']

            return [atoms, atompos]


    CG = coord[atoms.index('CG')]

    OD1 = coord[atoms.index('OD1')]

    CA_CB_CG_OD1_diangle = calcTorsionAngle(CA,CB,CG,OD1)

    CA_CB_CG_OD2_diangle = 180 + CA_CB_CG_OD1_diangle

    OD2 = calcCoordinate(CA, CB, CG, 1.25, 118.21, CA_CB_CG_OD2_diangle)

    atompos = [N, CA, C, O, CB, CG, OD1, OD2]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2']

    return [atoms, atompos]

def repairGln(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the GLN residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.81)

    if faspr:

        try:
            CB_CG_CD_OE1_diangle = get_chi('CB', 'CG', 'CD', 'OE1', resNo, chain,faspr)

            CA_CB_CG_CD_diangle = get_chi('CA', 'CB', 'CG', 'CD', resNo, chain,faspr)

            N_CA_CB_CG_diangle = get_chi('N', 'CA', 'CB', 'CG', resNo, chain,faspr)

        except:

            CB_CG_CD_OE1_diangle, CA_CB_CG_CD_diangle, N_CA_CB_CG_diangle = [],[],[]

    if faspr and not all([CB_CG_CD_OE1_diangle, CA_CB_CG_CD_diangle, N_CA_CB_CG_diangle]):

        N_CA_CB_CG_diangle = -60.2

        CB_CG_CD_OE1_diangle = -50.5

        CA_CB_CG_CD_diangle = -69.6

    if not faspr:

        N_CA_CB_CG_diangle = -60.2

        CB_CG_CD_OE1_diangle = -50.5

        CA_CB_CG_CD_diangle = -69.6

    CB_CG_CD_NE2_diangle = 180.0 + CB_CG_CD_OE1_diangle

    # for GLN, fixing CG,CD,OE1 require flexible chi
    require_chi = ['CG','CD','OE1']

    for i in missing:

        if i in require_chi:

            CG = calcCoordinate(N, CA, CB, 1.52, 113.75,N_CA_CB_CG_diangle)

            CD = calcCoordinate(CA, CB, CG, 1.52, 112.78, CA_CB_CG_CD_diangle)

            OE1 = calcCoordinate(CB, CG, CD, 1.24, 120.86, CB_CG_CD_OE1_diangle)

            NE2 = calcCoordinate(CB, CG, CD, 1.33, 116.50,CB_CG_CD_NE2_diangle)

            atompos = [N, CA, C, O, CB, CG, CD, OE1, NE2]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2']

            return [atoms, atompos]


    CG  = coord[atoms.index('CG')]

    CD  = coord[atoms.index('CD')]

    OE1 = coord[atoms.index('OE1')]

    CB_CG_CD_OE1_diangle = calcTorsionAngle(CB,CG,CD,OE1)

    CB_CG_CD_NE2_diangle = 180.0 + CB_CG_CD_OE1_diangle

    NE2 = calcCoordinate(CB, CG, CD, 1.33, 116.50,CB_CG_CD_NE2_diangle)

    atompos = [N, CA, C, O, CB, CG, CD, OE1, NE2]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2']

    return [atoms, atompos]

def repairLys(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the LYS residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.76)

    if faspr:

        try:

            N_CA_CB_CG_diangle  = get_chi('N', 'CA', 'CB', 'CG', resNo, chain,faspr)

            CA_CB_CG_CD_diangle = get_chi('CA', 'CB', 'CG', 'CD', resNo, chain,faspr)

            CB_CG_CD_CE_diangle = get_chi('CB', 'CG', 'CD', 'CE', resNo, chain,faspr)

            CG_CD_CE_NZ_diangle = get_chi('CG', 'CD', 'CE', 'NZ', resNo, chain,faspr)

        except:

            N_CA_CB_CG_diangle,CA_CB_CG_CD_diangle, CB_CG_CD_CE_diangle, CG_CD_CE_NZ_diangle = [],[],[],[]

    if faspr and not all([N_CA_CB_CG_diangle,CA_CB_CG_CD_diangle, CB_CG_CD_CE_diangle, CG_CD_CE_NZ_diangle]):

        N_CA_CB_CG_diangle = -64.5

        CA_CB_CG_CD_diangle = -178.1

        CB_CG_CD_CE_diangle = -179.6

        CG_CD_CE_NZ_diangle = 179.6

    if not faspr:

        N_CA_CB_CG_diangle = -64.5

        CA_CB_CG_CD_diangle = -178.1

        CB_CG_CD_CE_diangle = -179.6

        CG_CD_CE_NZ_diangle = 179.6

    # for LYS, fixing CG,CD,CE,NZ require flexible chi
    require_chi = ['CG', 'CD', 'CE', 'NZ']

    for i in missing:

        if i in require_chi:

            CG = calcCoordinate(N, CA, CB, 1.52, 113.83, N_CA_CB_CG_diangle)

            CD = calcCoordinate(CA, CB, CG, 1.52, 111.79, CA_CB_CG_CD_diangle)

            CE = calcCoordinate(CB, CG, CD, 1.46, 111.68, CB_CG_CD_CE_diangle)

            NZ = calcCoordinate(CG, CD, CE, 1.33, 124.79, CG_CD_CE_NZ_diangle)

            atompos = [N, CA, C, O, CB, CG, CD, CE, NZ]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ']

            return [atoms, atompos]


    CG = coord[atoms.index('CG')]

    CD = coord[atoms.index('CD')]

    CE = coord[atoms.index('CE')]

    NZ = coord[atoms.index('NZ')]

    atompos = [N, CA, C, O, CB, CG, CD, CE, NZ]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ']

    return [atoms, atompos]

def repairArg(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the ARG residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.76)

    if faspr:

        try:
            N_CA_CB_CG_diangle = get_chi('N', 'CA', 'CB', 'CG', resNo, chain,faspr)

            CA_CB_CG_CD_diangle = get_chi('CA', 'CB', 'CG', 'CD', resNo, chain,faspr)

            CB_CG_CD_NE_diangle = get_chi('CB', 'CG','CD', 'NE', resNo, chain,faspr)

            CG_CD_NE_CZ_diangle = get_chi( 'CG','CD', 'NE', 'CZ', resNo, chain,faspr)

        except:

            N_CA_CB_CG_diangle, CA_CB_CG_CD_diangle, CB_CG_CD_NE_diangle, CG_CD_NE_CZ_diangle = [],[],[],[]

    if faspr and not all([N_CA_CB_CG_diangle, CA_CB_CG_CD_diangle, CB_CG_CD_NE_diangle, CG_CD_NE_CZ_diangle]):

        N_CA_CB_CG_diangle = -65.2

        CA_CB_CG_CD_diangle = -179.2

        CB_CG_CD_NE_diangle = -179.3

        CG_CD_NE_CZ_diangle = -178.7

    if not faspr:

        N_CA_CB_CG_diangle = -65.2

        CA_CB_CG_CD_diangle = -179.2

        CB_CG_CD_NE_diangle = -179.3

        CG_CD_NE_CZ_diangle = -178.7

    # for ARG, fixing CG,CD,NE,CZ require flexible chi
    require_chi = ['CG', 'CD', 'NE', 'CZ']

    for i in missing:

        if i in require_chi:

            CG = calcCoordinate(N, CA, CB, 1.52, 113.83, N_CA_CB_CG_diangle)

            CD = calcCoordinate(CA, CB, CG, 1.52, 111.79, CA_CB_CG_CD_diangle)

            NE = calcCoordinate(CB, CG, CD, 1.46, 111.68, CB_CG_CD_NE_diangle)

            CZ = calcCoordinate(CG, CD, NE, 1.33, 124.79, CG_CD_NE_CZ_diangle)

            NH1 = calcCoordinate(CD, NE, CZ, 1.33, 120.64, 0.0)

            NH2 = calcCoordinate(CD, NE, CZ, 1.33, 119.63, 180.0)

            atompos = [N, CA, C, O, CB, CG, CD, NE, CZ, NH1, NH2]

            atoms   = ['N','CA','C','O','CB','CG','CD','NE','CZ','NH1','NH2']

            return [atoms, atompos]


    # this residue has more than one atom that do not require flexible chi
    # for execution to arrive here it means some or all are missing
    # so we try to check which is missing first
    CG = coord[atoms.index('CG')]

    CD = coord[atoms.index('CD')]

    NE = coord[atoms.index('NE')]

    CZ = coord[atoms.index('CZ')]

    try:
        NH1 = coord[atoms.index('NH1')]
    except:
        NH1 = calcCoordinate(CD, NE, CZ, 1.33, 120.64, 0.0)

    try:
        NH2 = coord[atoms.index('NH2')]
    except:
        NH2 = calcCoordinate(CD, NE, CZ, 1.33, 119.63, 180.0)

    atompos = [N, CA, C, O, CB, CG, CD, NE, CZ, NH1, NH2]

    atoms   = ['N','CA','C','O','CB','CG','CD','NE','CZ','NH1','NH2']

    return [atoms, atompos]

def repairTyr(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the TYR residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.60)

    if faspr:

        try:

            CA_CB_CG_CD1_diangle = get_chi('CA', 'CB', 'CG', 'CD1', resNo, chain,faspr)

            N_CA_CB_CG_diangle = get_chi('N', 'CA', 'CB', 'CG', resNo, chain,faspr)

        except:

            CA_CB_CG_CD1_diangle, N_CA_CB_CG_diangle = [],[]

    if faspr and not all([CA_CB_CG_CD1_diangle, N_CA_CB_CG_diangle]):

        CA_CB_CG_CD1_diangle = 93.1

        N_CA_CB_CG_diangle = -64.3

    if not faspr:

        CA_CB_CG_CD1_diangle = 93.1

        N_CA_CB_CG_diangle = -64.3

    CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle + 180.0

    # for TYR, fixing CG,CD1 require flexible chi
    require_chi = ['CG', 'CD1']

    for i in missing:

        if i in require_chi:

            CG = calcCoordinate(N, CA, CB, 1.51, 113.8,N_CA_CB_CG_diangle)

            CD1 = calcCoordinate(CA, CB, CG, 1.39, 120.98, CA_CB_CG_CD1_diangle)

            CD2 = calcCoordinate(CA, CB, CG, 1.39, 120.82,CA_CB_CG_CD2_diangle)

            CE1 = calcCoordinate(CB, CG, CD1, 1.39, 120.0,180.0)

            CE2 = calcCoordinate(CB, CG, CD2, 1.39, 120.0, 180.0)

            CZ = calcCoordinate(CG, CD1, CE1, 1.39, 120.0, 0.0)

            OH = calcCoordinate(CD1, CE1, CZ, 1.39, 119.78, 180.0)

            atompos = [N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, OH]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2','CZ', 'OH']

            return [atoms, atompos]


    CG = coord[atoms.index('CG')]

    CD1 = coord[atoms.index('CD1')]

    # first test which atom is present and then decide
    try:
        CD2 = coord[atoms.index('CD2')]

    except:

        CA_CB_CG_CD1_diangle = calcTorsionAngle(CA,CB,CG,CD1)

        CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle + 180.0

        CD2 = calcCoordinate(CA, CB, CG, 1.39, 120.82,CA_CB_CG_CD2_diangle)

    try:
        CE1 = coord[atoms.index('CE1')]
    except:
        CE1 = calcCoordinate(CB, CG, CD1, 1.39, 120.0,180.0)

    try:
        CE2 = coord[atoms.index('CE2')]
    except:
        CE2 = calcCoordinate(CB, CG, CD2, 1.39, 120.0, 180.0)

    try:
        CZ = coord[atoms.index('CZ')]
    except:
        CZ = calcCoordinate(CG, CD1, CE1, 1.39, 120.0, 0.0)

    try:
        OH = coord[atoms.index('OH')]
    except:
            OH = calcCoordinate(CD1, CE1, CZ, 1.39, 119.78, 180.0)

    atompos = [N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, OH]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH']

    return [atoms, atompos]

def repairTrp(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the TRP residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.61)

    if faspr:

        try:
            CA_CB_CG_CD1_diangle = get_chi('CA', 'CB', 'CG', 'CD1', resNo, chain,faspr)

            N_CA_CB_CG_diangle = get_chi('N', 'CA', 'CB', 'CG', resNo, chain,faspr)

        except:

            CA_CB_CG_CD1_diangle, N_CA_CB_CG_diangle = [],[]

    if faspr and not all([CA_CB_CG_CD1_diangle, N_CA_CB_CG_diangle]):

        N_CA_CB_CG_diangle = -66.4

        CA_CB_CG_CD1_diangle = 96.3

    if not faspr:

        N_CA_CB_CG_diangle = -66.4

        CA_CB_CG_CD1_diangle = 96.3

    CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle - 180.0

    # for TRP, fixing CG,CD1 require flexible chi
    require_chi = ['CG', 'CD1']

    for i in missing:

        if i in require_chi:

            CG = calcCoordinate(N, CA, CB, 1.50, 114.10, N_CA_CB_CG_diangle)

            CD1 = calcCoordinate(CA, CB, CG, 1.37, 127.07, CA_CB_CG_CD1_diangle)

            CD2 = calcCoordinate(CA, CB, CG, 1.43, 126.66, CA_CB_CG_CD2_diangle)

            NE1 = calcCoordinate(CB, CG, CD1, 1.38, 108.5, 180.0)

            CE2 = calcCoordinate(CB, CG, CD2, 1.40, 108.5, 180.0)

            CE3 = calcCoordinate(CB, CG, CD2, 1.40, 133.83, 0.0)

            CZ2 = calcCoordinate(CG, CD2, CE2, 1.40, 120.0, 180.0)

            CZ3 = calcCoordinate(CG, CD2, CE3, 1.40, 120.0, 180.0)

            CH2 = calcCoordinate(CD2, CE2, CZ2, 1.40, 120.0, 0.0)

            atompos = [N, CA, C, O, CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2','CE3', 'CZ2','CZ3', 'CH2']

            return [atoms, atompos]


    CG = coord[atoms.index('CG')]

    CD1 = coord[atoms.index('CD1')]

    # first test which atom is present and then decide
    try:

        CD2 = coord[atoms.index('CD2')]

    except:

        CA_CB_CG_CD1_diangle = calcTorsionAngle(CA,CB,CG,CD1)

        CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle - 180.0

        CD2 = calcCoordinate(CA, CB, CG, 1.43, 126.66,CA_CB_CG_CD2_diangle)

    try:
        NE1 = coord[atoms.index('NE1')]
    except:
        NE1 = calcCoordinate(CB, CG, CD1, 1.38, 108.5, 180.0)

    try:
        CE2 = coord[atoms.index('CE2')]
    except:
        CE2 = calcCoordinate(CB, CG, CD2, 1.40, 108.5, 180.0)

    try:
        CE3 = coord[atoms.index('CE3')]
    except:
        CE3 = calcCoordinate(CB, CG, CD2, 1.40, 133.83, 0.0)

    try:
        CZ2 = coord[atoms.index('CZ2')]
    except:
        CZ2 = calcCoordinate(CG, CD2, CE2, 1.40, 120.0, 180.0)

    try:
        CZ3 = coord[atoms.index('CZ3')]
    except:
        CZ3 = calcCoordinate(CG, CD2, CE3, 1.40, 120.0, 180.0)

    try:
        CH2 = coord[atoms.index('CH2')]
    except:
        CH2 = calcCoordinate(CD2, CE2, CZ2, 1.40, 120.0, 0.0)

    atompos = [N, CA, C, O, CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2','CZ3', 'CH2']

    return [atoms, atompos]

def repairPhe(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the PHE residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.61)

    if faspr:

        try:

            CA_CB_CG_CD1_diangle = get_chi('CA', 'CB', 'CG','CD1', resNo, chain,faspr)

            N_CA_CB_CG_diangle = get_chi('N', 'CA', 'CB', 'CG', resNo, chain,faspr)

        except:

            CA_CB_CG_CD1_diangle, N_CA_CB_CG_diangle = [],[]

    if faspr and not all([CA_CB_CG_CD1_diangle, N_CA_CB_CG_diangle]):

        N_CA_CB_CG_diangle = -64.7

        CA_CB_CG_CD1_diangle = 93.3

    if not faspr:

        N_CA_CB_CG_diangle = -64.7

        CA_CB_CG_CD1_diangle = 93.3

    CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle - 180.0

    # for PHE, fixing CG,CD1 require flexible chi
    require_chi = ['CG', 'CD1']

    for i in missing:

        if i in require_chi:

            CG = calcCoordinate(N, CA, CB, 1.50, 113.85, N_CA_CB_CG_diangle)

            CD1 = calcCoordinate(CA, CB, CG, 1.39, 120.0, CA_CB_CG_CD1_diangle)

            CD2 = calcCoordinate(CA, CB, CG, 1.39, 120.0, CA_CB_CG_CD2_diangle)

            CE1 = calcCoordinate(CB, CG, CD1, 1.39, 120.0, 180.0)

            CE2 = calcCoordinate(CB, CG, CD2, 1.39, 120.0, 180.0)

            CZ = calcCoordinate(CG, CD1, CE1, 1.39, 120.0, 0.0)

            atompos = [N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']

            return [atoms, atompos]


    CG = coord[atoms.index('CG')]

    CD1 = coord[atoms.index('CD1')]

    # first test which atom is present and then decide
    try:

        CD2 = coord[atoms.index('CD2')]

    except:

        CA_CB_CG_CD1_diangle = calcTorsionAngle(CA,CB,CG,CD1)

        CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle - 180.0

        CD2 = calcCoordinate(CA, CB, CG, 1.39, 120.0,CA_CB_CG_CD2_diangle)

    try:
        CE1 = coord[atoms.index('CE1')]
    except:
        CE1 = calcCoordinate(CB, CG, CD1, 1.39, 120.0, 180.0)

    try:
        CE2 = coord[atoms.index('CE2')]
    except:
        CE2 = calcCoordinate(CB, CG, CD2, 1.39, 120.0, 180.0)

    try:
        CZ = coord[atoms.index('CZ')]
    except:
        CZ = calcCoordinate(CG, CD1, CE1, 1.39, 120.0, 0.0)

    atompos = [N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']

    return [atoms, atompos]

def repairHis(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the HIS residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.67)

    if faspr:

        try:

            CA_CB_CG_ND1_diangle = get_chi('CA', 'CB', 'CG','ND1', resNo, chain,faspr)

            N_CA_CB_CG_diangle = get_chi('N', 'CA', 'CB', 'CG', resNo, chain,faspr)

        except:

            N_CA_CB_CG_diangle, CA_CB_CG_ND1_diangle = [],[]

    if faspr and not all([N_CA_CB_CG_diangle, CA_CB_CG_ND1_diangle]):

        CA_CB_CG_ND1_diangle = -75.7

        N_CA_CB_CG_diangle = -63.2

    if not faspr:

        CA_CB_CG_ND1_diangle = -75.7

        N_CA_CB_CG_diangle = -63.2

    CA_CB_CG_CD2_diangle = 180.0 + CA_CB_CG_ND1_diangle

    # for HIS, fixing CG,ND1 require flexible chi
    require_chi = ['CG', 'ND1']

    for i in missing:

        if i in require_chi:

            CG = calcCoordinate(N, CA, CB, 1.49, 113.74, N_CA_CB_CG_diangle)

            ND1 = calcCoordinate(CA, CB, CG, 1.38, 122.85, CA_CB_CG_ND1_diangle)

            CD2 = calcCoordinate(CA, CB, CG, 1.35, 130.61, CA_CB_CG_CD2_diangle)

            CE1 = calcCoordinate(CB, CG, ND1, 1.32, 108.5, 180.0)

            NE2 = calcCoordinate(CB, CG, CD2, 1.35, 108.5, 180.0)

            atompos = [N, CA, C, O, CB, CG, ND1, CD2, CE1, NE2]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2']

            return [atoms, atompos]


    CG = coord[atoms.index('CG')]

    ND1 = coord[atoms.index('ND1')]

    # first test which atom is present and then decide
    try:
        CD2 = coord[atoms.index('CD2')]
    except:
        CA_CB_CG_ND1_diangle = calcTorsionAngle(CA, CB, CG, ND1)

        CA_CB_CG_CD2_diangle = 180.0 + CA_CB_CG_ND1_diangle

        CD2 = calcCoordinate(CA, CB, CG, 1.35, 130.61, CA_CB_CG_CD2_diangle)

    try:
        CE1 = coord[atoms.index('CE1')]
    except:
        CE1 = calcCoordinate(CB, CG, ND1, 1.32, 108.5, 180.0)

    try:
        NE2 = coord[atoms.index('NE2')]
    except:
        NE2 = calcCoordinate(CB, CG, CD2, 1.35, 108.5, 180.0)

    atompos = [N, CA, C, O, CB, CG, ND1, CD2, CE1, NE2]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2']

    return [atoms, atompos]

def repairPro(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the PRO residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 115.30)

    if faspr:

        try:
            N_CA_CB_CG_diangle = get_chi('N', 'CA', 'CB', 'CG', resNo, chain,faspr)

            CA_CB_CG_CD_diangle = get_chi('CA', 'CB', 'CG','CD', resNo, chain,faspr)

        except:

            N_CA_CB_CG_diangle, CA_CB_CG_CD_diangle = [],[]

    if faspr and not all([N_CA_CB_CG_diangle, CA_CB_CG_CD_diangle]) :

        N_CA_CB_CG_diangle = 29.6

        CA_CB_CG_CD_diangle = -34.8

    if not faspr :

        N_CA_CB_CG_diangle = 29.6

        CA_CB_CG_CD_diangle = -34.8

    # for PRO, fixing CG,CD require flexible chi
    require_chi = ['CG', 'CD']

    for i in missing:

        if i in require_chi:

            CG = calcCoordinate(N, CA, CB, 1.49, 104.21, N_CA_CB_CG_diangle)

            CD = calcCoordinate(CA, CB, CG, 1.50, 105.03, CA_CB_CG_CD_diangle)

            atompos = [N, CA, C, O, CB, CG, CD]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD']

            return [atoms, atompos]


    CG = coord[atoms.index('CG')]

    CD = coord[atoms.index('CD')]

    atompos = [N, CA, C, O, CB, CG, CD]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD']

    return [atoms, atompos]

def repairMet(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the MET residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.67)

    if faspr:

        try:

            N_CA_CB_CG_diangle = get_chi('N', 'CA', 'CB', 'CG', resNo, chain,faspr)

            CA_CB_CG_SD_diangle = get_chi('CA', 'CB', 'CG','SD', resNo, chain,faspr)

            CB_CG_SD_CE_diangle = get_chi('CB', 'CG','SD', 'CE', resNo, chain,faspr)

        except:

            N_CA_CB_CG_diangle, CA_CB_CG_SD_diangle, CB_CG_SD_CE_diangle = [],[],[]

    if faspr and not all([N_CA_CB_CG_diangle, CA_CB_CG_SD_diangle, CB_CG_SD_CE_diangle]):

        N_CA_CB_CG_diangle = -64.4

        CA_CB_CG_SD_diangle = -179.6

        CB_CG_SD_CE_diangle = 70.1

    if not faspr:

        N_CA_CB_CG_diangle = -64.4

        CA_CB_CG_SD_diangle = -179.6

        CB_CG_SD_CE_diangle = 70.1

    # for MET, fixing CG,SD, CE require flexible chi
    require_chi = ['CG', 'SD','CE']

    for i in missing:

        if i in require_chi:

            CG = calcCoordinate(N, CA, CB, 1.52, 113.68, N_CA_CB_CG_diangle)

            SD = calcCoordinate(CA, CB, CG, 1.81, 112.69, CA_CB_CG_SD_diangle)

            CE = calcCoordinate(CB, CG, SD, 1.79, 100.61, CB_CG_SD_CE_diangle)

            atompos = [N, CA, C, O, CB, CG, SD, CE]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE']

            return [atoms, atompos]


    CG = coord[atoms.index('CG')]

    SD = coord[atoms.index('SD')]

    CE = coord[atoms.index('CE')]

    atompos = [N, CA, C, O, CB, CG, SD, CE]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE']

    return [atoms, atompos]

def repairIle(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the ILE residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 123.23)

    if faspr:

        try:
            N_CA_CB_CG1_diangle = get_chi('N', 'CA', 'CB', 'CG1', resNo, chain,faspr)

            CA_CB_CG1_CD1_diangle = get_chi('CA', 'CB', 'CG1','CD1',  resNo, chain,faspr)

        except:

            N_CA_CB_CG1_diangle, CA_CB_CG1_CD1_diangle = [],[]

    if faspr and not all([N_CA_CB_CG1_diangle, CA_CB_CG1_CD1_diangle]):

        N_CA_CB_CG1_diangle = 59.7

        CA_CB_CG1_CD1_diangle = 169.8

    if not faspr:

        N_CA_CB_CG1_diangle = 59.7

        CA_CB_CG1_CD1_diangle = 169.8

    N_CA_CB_CG2_diangle = N_CA_CB_CG1_diangle-120

    # for ILE, fixing CG1,CD1 require flexible chi
    require_chi = ['CG1', 'CD1']

    for i in missing:

        if i in require_chi:

            CG1 = calcCoordinate(N, CA, CB, 1.527, 110.7, N_CA_CB_CG1_diangle)

            CG2 = calcCoordinate(N, CA, CB, 1.527, 110.4, N_CA_CB_CG2_diangle)

            CD1 = calcCoordinate(CA, CB, CG1, 1.52, 113.97, CA_CB_CG1_CD1_diangle)

            atompos = [N, CA, C, O, CB, CG1, CG2, CD1]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1']

            return [atoms, atompos]


    CG1 = coord[atoms.index('CG1')]

    CD1 = coord[atoms.index('CD1')]

    N_CA_CB_CG1_diangle = calcTorsionAngle(N,CA,CB,CG1)

    N_CA_CB_CG2_diangle = N_CA_CB_CG1_diangle-120

    CG2 = calcCoordinate(N, CA, CB, 1.527, 110.4, N_CA_CB_CG2_diangle)

    atompos = [N, CA, C, O, CB, CG1, CG2, CD1]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1']

    return [atoms, atompos]

def repairLeu(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the LEU residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.49)

    if faspr:

        try:

            N_CA_CB_CG_diangle = get_chi('N', 'CA', 'CB', 'CG', resNo, chain,faspr)

            CA_CB_CG_CD1_diangle = get_chi( 'CA', 'CB', 'CG','CD1', resNo, chain,faspr)

        except:

            N_CA_CB_CG_diangle, CA_CB_CG_CD1_diangle = [],[]

    if faspr and not all([N_CA_CB_CG_diangle, CA_CB_CG_CD1_diangle]):

        N_CA_CB_CG_diangle = -60.1

        CA_CB_CG_CD1_diangle = 174.9

    if not faspr:

        N_CA_CB_CG_diangle = -60.1

        CA_CB_CG_CD1_diangle = 174.9

    CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle + 120

    # for LEU, fixing CG,CD1 require flexible chi
    require_chi = ['CG', 'CD1']

    for i in missing:

        if i in require_chi:

            CG = calcCoordinate(N, CA, CB, 1.53, 116.10, N_CA_CB_CG_diangle)

            CD1 = calcCoordinate(CA, CB, CG, 1.524, 112.50,CA_CB_CG_CD1_diangle)

            CD2 = calcCoordinate(CA, CB, CG, 1.525, 112.50, CA_CB_CG_CD2_diangle)

            atompos = [N, CA, C, O, CB, CG, CD1, CD2]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2']

            return [atoms, atompos]


    CG = coord[atoms.index('CG')]

    CD1 = coord[atoms.index('CD1')]

    CA_CB_CG_CD1_diangle = calcTorsionAngle(CA,CB,CG,CD1)

    CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle+120

    CD2 = calcCoordinate(CA, CB, CG, 1.525, 112.50, CA_CB_CG_CD2_diangle)

    atompos = [N, CA, C, O, CB, CG, CD1, CD2]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2']

    return [atoms, atompos]

def repairVal(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the VAL residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 123.23)

    if faspr:

        try:
            N_CA_CB_CG1_diangle = get_chi('N', 'CA', 'CB', 'CG1',resNo, chain,faspr)

        except:

            N_CA_CB_CG1_diangle = []

    if faspr and not N_CA_CB_CG1_diangle:

        N_CA_CB_CG1_diangle = 177.2

    if not faspr:

        N_CA_CB_CG1_diangle = 177.2

    N_CA_CB_CG2_diangle = N_CA_CB_CG1_diangle + 120

    # for VAL, fixing CG1 requires flexible chi
    require_chi = ['CG1']

    for i in missing:

        if i in require_chi:

            CG1 = calcCoordinate(N, CA, CB, 1.527, 110.7, N_CA_CB_CG1_diangle)

            CG2 = calcCoordinate(N, CA, CB, 1.527, 110.4, N_CA_CB_CG2_diangle)

            atompos = [N, CA, C, O, CB, CG1, CG2]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2']

            return [atoms, atompos]


    CG1 = coord[atoms.index('CG1')]

    N_CA_CB_CG1_diangle = calcTorsionAngle(N,CA,CB,CG1)

    N_CA_CB_CG2_diangle = N_CA_CB_CG1_diangle + 120

    CG2 = calcCoordinate(N, CA, CB, 1.527, 110.4, N_CA_CB_CG2_diangle)

    atompos = [N, CA, C, O, CB, CG1, CG2]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2']

    return [atoms, atompos]

def repairThr(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the THR residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 123.10)

    if faspr:

        try:

            N_CA_CB_OG1_diangle = get_chi('N', 'CA', 'CB', 'OG1',resNo, chain,faspr)

        except:

            N_CA_CB_OG1_diangle = []

    if faspr and not N_CA_CB_OG1_diangle:

        N_CA_CB_OG1_diangle = 60.0

    if not faspr:

        N_CA_CB_OG1_diangle = 60.0

    N_CA_CB_CG2_diangle = N_CA_CB_OG1_diangle - 120

    # for THR, fixing OG1 requires flexible chi
    require_chi = ['OG1']

    for i in missing:

        if i in require_chi:

            OG1 = calcCoordinate(N, CA, CB, 1.43, 109.18, N_CA_CB_OG1_diangle)

            CG2 = calcCoordinate(N, CA, CB, 1.53, 111.13, N_CA_CB_CG2_diangle)

            atompos = [N, CA, C, O, CB, OG1, CG2]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2']

            return [atoms, atompos]

    OG1 = coord[atoms.index('OG1')]

    N_CA_CB_OG1_diangle = calcTorsionAngle(N,CA,CB,OG1)

    N_CA_CB_CG2_diangle = N_CA_CB_OG1_diangle - 120

    CG2 = calcCoordinate(N, CA, CB, 1.53, 111.13, N_CA_CB_CG2_diangle)

    atompos = [N, CA, C, O, CB, OG1, CG2]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2']

    return [atoms, atompos]

def repairCys(atoms,coord,resNo,chain,missing,faspr,nextres='',psi=''):

    """
    This function adds missing heavy
    atoms of the CYS residue.
    See repairSer() for full annotation
    """

    N  = coord[atoms.index('N')]
    CA = coord[atoms.index('CA')]
    C  = coord[atoms.index('C')]

    try:
        O = coord[atoms.index('O')]
    except:
        O = addBackbone(N,CA,C,nextres,psi)

    try:
        CB = coord[atoms.index('CB')]
    except:
        CB = calcCoordinate(N, C, CA, 1.52, 109.5, 122.50)

    if faspr:

        try:
            N_CA_CB_SG_diangle = get_chi('N', 'CA', 'CB', 'SG', resNo, chain,faspr)

        except:

            N_CA_CB_SG_diangle = []

    if faspr and not N_CA_CB_SG_diangle:

        N_CA_CB_SG_diangle = -62.2

    if not faspr:

        N_CA_CB_SG_diangle = -62.2

    # for CYS, fixing SG requires flexible chi
    require_chi = ['SG']

    for i in missing:

        if i in require_chi:

            SG = calcCoordinate(N, CA, CB, 1.81, 113.82, N_CA_CB_SG_diangle)

            atompos = [N, CA, C, O, CB, SG]

            atoms   = ['N', 'CA', 'C', 'O', 'CB', 'SG']

            return [atoms, atompos]


    SG = coord[atoms.index('SG')]

    atompos = [N, CA, C, O, CB, SG]

    atoms   = ['N', 'CA', 'C', 'O', 'CB', 'SG']

    return [atoms, atompos]
