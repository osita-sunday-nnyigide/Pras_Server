#!/usr/bin/env python

__doc__ = """

This program requires python 3.6 or higher.

This module, MissingHydrogenAtoms.py has all

the functions required to fix missing hydrogen

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

import copy
from math import sqrt
from . import LinearAlgebra
from . import PotentialEnergy

def arginine_h(ires,i="",nextres=None):
    """
    This function plus the next 19 represent 20 amino acid residues.
    Annotation is provided in detail here only  b/c all the
    functions are essentially the same, that is, they take
    similar arguments except SER, THR, TYR and free CYS that have class6 H-atoms
    to be optimized with a potential energy function. The last 6 functions
    are different from these 20 and annotations are provided

    This function adds missing hydrogen atoms of the ARG residue.

    Arguments
    ----------
    ires:    a list containing the residue name,atoms,atom coordinates
    i:       integer, the index that identifies this residue
    nextres: same as ires but for the next residue in the chain

    Extra Arguments for SER, THR, TYR and free CYS
    ----------
    resNo:   the residue number of the residue (not the index)
    atoms2:  a list of all atoms in the chain
    atmpos2: a list of all coordinates of all atoms in the chain
    resNo2:  a list of all residue numbers in the chain
    resseq2: a list of all residue names in the chain

    Returns
    -------
    A list: index 0 is the atom coordinates,
            index 1 is backbone H if the next residue is not PRO
            index 2 is the atom name, however,
            if the next residue is PRO, index 1 is the atom name
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_CG     = coord[atoms.index('CG')]
    pos_CD     = coord[atoms.index('CD')]
    pos_NE     = coord[atoms.index('NE')]
    pos_CZ     = coord[atoms.index('CZ')]
    pos_NH1    = coord[atoms.index('NH1')]
    pos_NH2    = coord[atoms.index('NH2')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H         =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA        =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HE        =  LinearAlgebra.class5(pos_CD,pos_NE,pos_CZ,1.01)
        HH11      =  LinearAlgebra.class4(pos_NE,pos_CZ,pos_NH1)
        HH12      =  LinearAlgebra.class4(pos_NH2,pos_CZ,pos_NH1)
        HH21      =  LinearAlgebra.class4(pos_NE,pos_CZ,pos_NH2)
        HH22      =  LinearAlgebra.class4(pos_NH1,pos_CZ,pos_NH2)
        HB1, HB2  =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
        HG1, HG2  =  LinearAlgebra.class2(pos_CD,pos_CB,pos_CG)
        HD1, HD2  =  LinearAlgebra.class2(pos_NE,pos_CG,pos_CD)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HE','1HH1','2HH1','1HH2','2HH2','HB1','HB2','HG1','HG2','HD1','HD2'])
        coord_copy.extend([HA,HE,HH11,HH12,HH21,HH22,HB1,HB2,HG1,HG2,HD1,HD2])
        return coord_copy, [H,i], atoms_copy

    HA        =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HE        =  LinearAlgebra.class5(pos_CD,pos_NE,pos_CZ,1.01)
    HH11      =  LinearAlgebra.class4(pos_NE,pos_CZ,pos_NH1)
    HH12      =  LinearAlgebra.class4(pos_NH2,pos_CZ,pos_NH1)
    HH21      =  LinearAlgebra.class4(pos_NE,pos_CZ,pos_NH2)
    HH22      =  LinearAlgebra.class4(pos_NH1,pos_CZ,pos_NH2)
    HB1, HB2  =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
    HG1, HG2  =  LinearAlgebra.class2(pos_CD,pos_CB,pos_CG)
    HD1, HD2  =  LinearAlgebra.class2(pos_NE,pos_CG,pos_CD)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA','HE','1HH1','2HH1','1HH2','2HH2','HB1','HB2','HG1','HG2','HD1','HD2'])
    coord_copy.extend([HA,HE,HH11,HH12,HH21,HH22,HB1,HB2,HG1,HG2,HD1,HD2])
    return coord_copy, atoms_copy


def isoleucine_h(ires,i="",nextres=None):
    """
    This function adds missing hydrogen
    atoms of the ILE residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_CG1    = coord[atoms.index('CG1')]
    pos_CG2    = coord[atoms.index('CG2')]
    pos_CD1    = coord[atoms.index('CD1')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB             =  LinearAlgebra.class3(pos_CA,pos_CG1,pos_CB)
        HG11,HG12      =  LinearAlgebra.class2(pos_CB,pos_CD1, pos_CG1)
        HG21           =  LinearAlgebra.calcCoordinate(pos_CA,pos_CB,pos_CG2, 1.09, 180.0,109.5)
        HG22           =  LinearAlgebra.calcCoordinate(pos_CA,pos_CB,pos_CG2, 1.09,-60.0, 109.5)
        HG23           =  LinearAlgebra.calcCoordinate(pos_CA,pos_CB,pos_CG2, 1.09, 60.0, 109.5)
        HD1            =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG1,pos_CD1,1.09,-180.0,109.5)
        HD2            =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG1,pos_CD1,1.09,-60.0, 109.5)
        HD3            =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG1,pos_CD1,1.09, 60.0, 109.5)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB','1HG1','2HG1','1HG2','2HG2','3HG2','HD1','HD2','HD3'])
        coord_copy.extend([HA,HB,HG11,HG12,HG21,HG22,HG23,HD1,HD2,HD3])
        return coord_copy, [H,i], atoms_copy

    HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HB             =  LinearAlgebra.class3(pos_CA,pos_CG1,pos_CB)
    HG11,HG12      =  LinearAlgebra.class2(pos_CB,pos_CD1, pos_CG1)
    HG21           =  LinearAlgebra.calcCoordinate(pos_CA,pos_CB,pos_CG2, 1.09, 180.0,109.5)
    HG22           =  LinearAlgebra.calcCoordinate(pos_CA,pos_CB,pos_CG2, 1.09,-60.0, 109.5)
    HG23           =  LinearAlgebra.calcCoordinate(pos_CA,pos_CB,pos_CG2, 1.09, 60.0, 109.5)
    HD1            =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG1,pos_CD1,1.09,-180.0,109.5)
    HD2            =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG1,pos_CD1,1.09,-60.0, 109.5)
    HD3            =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG1,pos_CD1,1.09, 60.0, 109.5)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA','HB','1HG1','2HG1','1HG2','2HG2','3HG2','HD1','HD2','HD3'])
    coord_copy.extend([HA,HB,HG11,HG12,HG21,HG22,HG23,HD1,HD2,HD3])
    return coord_copy, atoms_copy

def leucine_h(ires,i="",nextres=None):
    """
    This function adds missing hydrogen
    atoms of the LEU residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_CG     = coord[atoms.index('CG')]
    pos_CD1    = coord[atoms.index('CD1')]
    pos_CD2    = coord[atoms.index('CD2')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
        HG             =  LinearAlgebra.class3(pos_CD1,pos_CB,pos_CG)
        HD11           =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG,pos_CD1,1.09, 180.0,109.5)
        HD12           =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG,pos_CD1,1.09, -60.0,109.5)
        HD13           =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG,pos_CD1,1.09,  60.0,109.5)
        HD21           =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG,pos_CD2,1.09,-180.0,109.5)
        HD22           =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG,pos_CD2,1.09, -60.0,109.5)
        HD23           =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG,pos_CD2, 1.09, 60.0,109.5)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB1','HB2','HG','1HD1','2HD1','3HD1','1HD2','2HD2','3HD2'])
        coord_copy.extend([HA,HB1,HB2,HG,HD11,HD12,HD13,HD21,HD22,HD23])
        return coord_copy, [H,i], atoms_copy

    HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
    HG             =  LinearAlgebra.class3(pos_CD1,pos_CB,pos_CG)
    HD11           =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG,pos_CD1, 1.09,180.0,109.5)
    HD12           =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG,pos_CD1, 1.09,-60.0,109.5)
    HD13           =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG,pos_CD1, 1.09, 60.0,109.5)
    HD21           =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG,pos_CD2,1.09,-180.0,109.5)
    HD22           =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG,pos_CD2, 1.09,-60.0,109.5)
    HD23           =  LinearAlgebra.calcCoordinate(pos_CB,pos_CG,pos_CD2,  1.09,60.0,109.5)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA','HB1','HB2','HG','1HD1','2HD1','3HD1','1HD2','2HD2','3HD2'])
    coord_copy.extend([HA,HB1,HB2,HG,HD11,HD12,HD13,HD21,HD22,HD23])
    return coord_copy, atoms_copy

def tryptophan_h(ires,i="",nextres=None):
    """
    This function adds missing hydrogen
    atoms of the TRP residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_CG     = coord[atoms.index('CG')]
    pos_CD1    = coord[atoms.index('CD1')]
    pos_CD2    = coord[atoms.index('CD2')]
    pos_CE2    = coord[atoms.index('CE2')]
    pos_CE3    = coord[atoms.index('CE3')]
    pos_NE1    = coord[atoms.index('NE1')]
    pos_CZ2    = coord[atoms.index('CZ2')]
    pos_CZ3    = coord[atoms.index('CZ3')]
    pos_CH2    = coord[atoms.index('CH2')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
        HD1            =  LinearAlgebra.class5(pos_NE1,pos_CD1,pos_CG,1.08)
        HE1            =  LinearAlgebra.class5(pos_CE2,pos_NE1,pos_CD1,1.01)
        HE3            =  LinearAlgebra.class5(pos_CZ3,pos_CE3,pos_CD2,1.08)
        HZ2            =  LinearAlgebra.class5(pos_CE2,pos_CZ2,pos_CH2,1.08)
        HZ3            =  LinearAlgebra.class5(pos_CH2,pos_CZ3,pos_CE3,1.08)
        HH2            =  LinearAlgebra.class5(pos_CZ2,pos_CH2,pos_CZ3,1.08)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB1','HB2','HD1','HE1','HE3','HZ2','HZ3','HH2'])
        coord_copy.extend([HA,HB1,HB2,HD1,HE1,HE3,HZ2,HZ3,HH2])
        return coord_copy, [H,i], atoms_copy

    HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
    HD1            =  LinearAlgebra.class5(pos_NE1,pos_CD1,pos_CG,1.01)
    HE1            =  LinearAlgebra.class5(pos_CE2,pos_NE1,pos_CD1,1.08)
    HE3            =  LinearAlgebra.class5(pos_CZ3,pos_CE3,pos_CD2,1.08)
    HZ2            =  LinearAlgebra.class5(pos_CE2,pos_CZ2,pos_CH2,1.08)
    HZ3            =  LinearAlgebra.class5(pos_CH2,pos_CZ3,pos_CE3,1.08)
    HH2            =  LinearAlgebra.class5(pos_CZ2,pos_CH2,pos_CZ3,1.08)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA','HB1','HB2','HD1','HE1','HE3','HZ2','HZ3','HH2'])
    coord_copy.extend([HA,HB1,HB2,HD1,HE1,HE3,HZ2,HZ3,HH2])
    return coord_copy, atoms_copy

def phenylalanine_h(ires,i="",nextres=None):
    """
    This function adds missing hydrogen
    atoms of the PHE residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_CG     = coord[atoms.index('CG')]
    pos_CD1    = coord[atoms.index('CD1')]
    pos_CD2    = coord[atoms.index('CD2')]
    pos_CE1    = coord[atoms.index('CE1')]
    pos_CE2    = coord[atoms.index('CE2')]
    pos_CZ     = coord[atoms.index('CZ')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
        HD1            =  LinearAlgebra.class5(pos_CE1,pos_CD1,pos_CG,1.08)
        HD2            =  LinearAlgebra.class5(pos_CE2,pos_CD2,pos_CG,1.08)
        HE1            =  LinearAlgebra.class5(pos_CZ,pos_CE1,pos_CD1,1.08)
        HE2            =  LinearAlgebra.class5(pos_CD2,pos_CE2,pos_CZ,1.08)
        HZ             =  LinearAlgebra.class5(pos_CE2,pos_CZ,pos_CE1,1.08)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB1','HB2','HD1','HD2','HE1','HE2','HZ'])
        coord_copy.extend([HA,HB1,HB2,HD1,HD2,HE1,HE2,HZ])
        return coord_copy, [H,i], atoms_copy

    HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
    HD1            =  LinearAlgebra.class5(pos_CE1,pos_CD1,pos_CG,1.08)
    HD2            =  LinearAlgebra.class5(pos_CE2,pos_CD2,pos_CG,1.08)
    HE1            =  LinearAlgebra.class5(pos_CZ,pos_CE1,pos_CD1,1.08)
    HE2            =  LinearAlgebra.class5(pos_CD2,pos_CE2,pos_CZ,1.08)
    HZ             =  LinearAlgebra.class5(pos_CE2,pos_CZ,pos_CE1,1.08)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA','HB1','HB2','HD1','HD2','HE1','HE2','HZ'])
    coord_copy.extend([HA,HB1,HB2,HD1,HD2,HE1,HE2,HZ])
    return coord_copy, atoms_copy

def asparagine_h(ires,i="",nextres=None):
    """
    This function adds missing hydrogen
    atoms of the ASN residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_CG     = coord[atoms.index('CG')]
    pos_ND2    = coord[atoms.index('ND2')]
    pos_OD1    = coord[atoms.index('OD1')]

    if nextres !=None:
        nextres_CA= nextres[2][1]
        nextres_N = nextres[2][0]
        H         =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA        =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB1,HB2   =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
        HD12      =  LinearAlgebra.class4(pos_CB,pos_CG,pos_ND2)
        HD22      =  LinearAlgebra.class4(pos_OD1,pos_CG,pos_ND2)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB1','HB2','1HD2','2HD2'])
        coord_copy.extend([HA,HB1,HB2,HD12,HD22])
        return coord_copy, [H,i], atoms_copy

    HA        =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HB1,HB2   =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
    HD12      =  LinearAlgebra.class4(pos_CB,pos_CG,pos_ND2)
    HD22      =  LinearAlgebra.class4(pos_OD1,pos_CG,pos_ND2)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA','HB1','HB2','1HD2','2HD2'])
    coord_copy.extend([HA,HB1,HB2,HD12,HD22])
    return coord_copy, atoms_copy

def alanine_h(ires,i="",nextres=None):
    """
    This function adds missing hydrogen
    atoms of the ALA residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB1            =  LinearAlgebra.calcCoordinate(pos_C,pos_CA,pos_CB,1.09,55.8,109.5)
        HB2            =  LinearAlgebra.calcCoordinate(pos_C,pos_CA,pos_CB,1.09,175.8,109.5)
        HB3            =  LinearAlgebra.calcCoordinate(pos_C,pos_CA,pos_CB,1.09,-64.2,109.5)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA', 'HB1', 'HB2', 'HB3'])
        coord_copy.extend([HA,HB1,HB2,HB3])
        return coord_copy,[H,i], atoms_copy

    HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HB1            =  LinearAlgebra.calcCoordinate(pos_C,pos_CA,pos_CB,1.09,55.8,109.5)
    HB2            =  LinearAlgebra.calcCoordinate(pos_C,pos_CA,pos_CB,1.09,175.8,109.5)
    HB3            =  LinearAlgebra.calcCoordinate(pos_C,pos_CA,pos_CB,1.09,-64.2,109.5)
    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA', 'HB1', 'HB2', 'HB3'])
    coord_copy.extend([HA,HB1,HB2,HB3])
    return coord_copy, atoms_copy

def glycine_h(ires,i="",nextres=None):
    """
    This function adds missing hydrogen
    atoms of the GLY residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H         =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA1,HA2   =  LinearAlgebra.class2(pos_C,pos_N,pos_CA)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA1', 'HA2'])
        coord_copy.extend([HA1,HA2])
        return coord_copy, [H,i], atoms_copy

    HA1,HA2   =  LinearAlgebra.class2(pos_C,pos_N,pos_CA)
    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA1', 'HA2'])
    coord_copy.extend([HA1,HA2])
    return coord_copy, atoms_copy

def valine_h(ires,i="",nextres=None):
    """
    This function adds missing hydrogen
    atoms of the VAL residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_CG1    = coord[atoms.index('CG1')]
    pos_CG2    = coord[atoms.index('CG2')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB             =  LinearAlgebra.class3(pos_CA,pos_CG2,pos_CB)
        HG11           =  LinearAlgebra.calcCoordinate(pos_CA,pos_CB,pos_CG1,1.09,180.0,109.5)
        HG12           =  LinearAlgebra.calcCoordinate(pos_CA,pos_CB,pos_CG1,1.09,-60.0,109.5)
        HG13           =  LinearAlgebra.calcCoordinate(pos_CA,pos_CB,pos_CG1,1.09,60.0,109.5)
        HG21           =  LinearAlgebra.calcCoordinate(pos_CG1,pos_CB,pos_CG2,1.09,-58.1,109.5)
        HG22           =  LinearAlgebra.calcCoordinate(pos_CG1,pos_CB,pos_CG2,1.09,61.9,109.5)
        HG23           =  LinearAlgebra.calcCoordinate(pos_CG1,pos_CB,pos_CG2,1.09,-178.1,109.5)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB','1HG1','2HG1','3HG1','1HG2','2HG2','3HG2'])
        coord_copy.extend([HA,HB,HG11,HG12,HG13,HG21,HG22,HG23])
        return coord_copy, [H,i], atoms_copy

    HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HB             =  LinearAlgebra.class3(pos_CA,pos_CG2,pos_CB)
    HG11           =  LinearAlgebra.calcCoordinate(pos_CA,pos_CB,pos_CG1,1.09,180.0,109.5)
    HG12           =  LinearAlgebra.calcCoordinate(pos_CA,pos_CB,pos_CG1,1.09,-60.0,109.5)
    HG13           =  LinearAlgebra.calcCoordinate(pos_CA,pos_CB,pos_CG1,1.09,60.0,109.5)
    HG21           =  LinearAlgebra.calcCoordinate(pos_CG1,pos_CB,pos_CG2,1.09,-58.1,109.5)
    HG22           =  LinearAlgebra.calcCoordinate(pos_CG1,pos_CB,pos_CG2,1.09,61.9,109.5)
    HG23           =  LinearAlgebra.calcCoordinate(pos_CG1,pos_CB,pos_CG2,1.09,-178.1,109.5)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA','HB','1HG1','2HG1','3HG1','1HG2','2HG2','3HG2'])
    coord_copy.extend([HA,HB,HG11,HG12,HG13,HG21,HG22,HG23])
    return coord_copy, atoms_copy

def lysine_h(ires,i="",nextres=None):
    """
    This function adds missing hydrogen
    atoms of the LYS residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_CG     = coord[atoms.index('CG')]
    pos_CD     = coord[atoms.index('CD')]
    pos_CE     = coord[atoms.index('CE')]
    pos_NZ     = coord[atoms.index('NZ')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
        HG1,HG2        =  LinearAlgebra.class2(pos_CB,pos_CD,pos_CG)
        HD1,HD2        =  LinearAlgebra.class2(pos_CE,pos_CG,pos_CD)
        HE1,HE2        =  LinearAlgebra.class2(pos_NZ,pos_CD,pos_CE)
        HZ1            =  LinearAlgebra.calcCoordinate(pos_CD,pos_CE,pos_NZ,1.01,180.0,109.3)
        HZ2            =  LinearAlgebra.calcCoordinate(pos_CD,pos_CE,pos_NZ,1.01,-60.0,109.5)
        HZ3            =  LinearAlgebra.calcCoordinate(pos_CD,pos_CE,pos_NZ,1.01,60.0,109.5)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB1','HB2','HG1','HG2','HD1','HD2','HE1','HE2','HZ1','HZ2','HZ3'])
        coord_copy.extend([HA,HB1,HB2,HG1,HG2,HD1,HD2,HE1,HE2,HZ1,HZ2,HZ3])
        return coord_copy, [H,i], atoms_copy

    HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
    HG1,HG2        =  LinearAlgebra.class2(pos_CB,pos_CD,pos_CG)
    HD1,HD2        =  LinearAlgebra.class2(pos_CE,pos_CG,pos_CD)
    HE1,HE2        =  LinearAlgebra.class2(pos_NZ,pos_CD,pos_CE)
    HZ1            =  LinearAlgebra.calcCoordinate(pos_CD,pos_CE,pos_NZ,1.01,180.0,109.3)
    HZ2            =  LinearAlgebra.calcCoordinate(pos_CD,pos_CE,pos_NZ,1.01,-60.0,109.5)
    HZ3            =  LinearAlgebra.calcCoordinate(pos_CD,pos_CE,pos_NZ,1.01,60.0,109.5)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA','HB1','HB2','HG1','HG2','HD1','HD2','HE1','HE2','HZ1','HZ2','HZ3'])
    coord_copy.extend([HA,HB1,HB2,HG1,HG2,HD1,HD2,HE1,HE2,HZ1,HZ2,HZ3])
    return coord_copy, atoms_copy

def aspartate_h(ires,i="",nextres=None):
    """
    This function adds missing hydrogen
    atoms of the ASP residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_CG     = coord[atoms.index('CG')]
    pos_OD1    = coord[atoms.index('OD1')]
    pos_OD2    = coord[atoms.index('OD2')]

    if nextres !=None:
        nextres_CA= nextres[2][1]
        nextres_N = nextres[2][0]
        H         = LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA        = LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB1,HB2   = LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA', 'HB1', 'HB2'])
        coord_copy.extend([HA,HB1,HB2])
        return coord_copy, [H,i], atoms_copy

    HA      = LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HB1,HB2 = LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA', 'HB1', 'HB2'])
    coord_copy.extend([HA,HB1,HB2])
    return coord_copy, atoms_copy

def glutamine_h(ires,i="",nextres=None):
    """
    This function adds missing hydrogen
    atoms of the GLN residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_CG     = coord[atoms.index('CG')]
    pos_CD     = coord[atoms.index('CD')]
    pos_NE2    = coord[atoms.index('NE2')]
    pos_OE1    = coord[atoms.index('OE1')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
        HG1,HG2        =  LinearAlgebra.class2(pos_CB,pos_CD,pos_CG)
        HE21           =  LinearAlgebra.class4(pos_CG,pos_CD,pos_NE2)
        HE22           =  LinearAlgebra.class4(pos_OE1,pos_CD,pos_NE2)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB1','HB2','HG1','HG2','1HE2','2HE2'])
        coord_copy.extend([HA,HB1,HB2,HG1,HG2,HE21,HE22])
        return coord_copy, [H,i], atoms_copy

    HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
    HG1,HG2        =  LinearAlgebra.class2(pos_CB,pos_CD,pos_CG)
    HE21           =  LinearAlgebra.class4(pos_CG,pos_CD,pos_NE2)
    HE22           =  LinearAlgebra.class4(pos_OE1,pos_CD,pos_NE2)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA','HB1','HB2','HG1','HG2','1HE2','2HE2'])
    coord_copy.extend([HA,HB1,HB2,HG1,HG2,HE21,HE22])
    return coord_copy, atoms_copy

def glutamate_h(ires,i="",nextres=None):
    """
    This function adds missing hydrogen
    atoms of the GLU residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_CG     = coord[atoms.index('CG')]
    pos_CD     = coord[atoms.index('CD')]
    pos_OE1    = coord[atoms.index('OE1')]
    pos_OE2    = coord[atoms.index('OE2')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
        HG1,HG2        =  LinearAlgebra.class2(pos_CB,pos_CD,pos_CG)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB1','HB2','HG1','HG2'])
        coord_copy.extend([HA,HB1,HB2,HG1,HG2])
        return coord_copy, [H,i], atoms_copy

    HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
    HG1,HG2        =  LinearAlgebra.class2(pos_CB,pos_CD,pos_CG)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA','HB1','HB2','HG1','HG2'])
    coord_copy.extend([HA,HB1,HB2,HG1,HG2])
    return coord_copy, atoms_copy

def histidine_h(ires,i="",nextres=None):
    """
    This function adds missing hydrogen
    atoms of the HIS residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_CG     = coord[atoms.index('CG')]
    pos_CD2    = coord[atoms.index('CD2')]
    pos_ND1    = coord[atoms.index('ND1')]
    pos_CE1    = coord[atoms.index('CE1')]
    pos_NE2    = coord[atoms.index('NE2')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
        HD2            =  LinearAlgebra.class5(pos_NE2,pos_CD2,pos_CG,1.08)
        HE1            =  LinearAlgebra.class5(pos_NE2,pos_CE1,pos_ND1,1.08)
        HE2            =  LinearAlgebra.class5(pos_CD2,pos_NE2,pos_CE1,1.01)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB1','HB2','HD2','HE1','HE2'])
        coord_copy.extend([HA,HB1,HB2,HD2,HE1,HE2])
        return coord_copy, [H,i], atoms_copy

    HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
    HD2            =  LinearAlgebra.class5(pos_NE2,pos_CD2,pos_CG,1.08)
    HE1            =  LinearAlgebra.class5(pos_NE2,pos_CE1,pos_ND1,1.08)
    HE2            =  LinearAlgebra.class5(pos_CD2,pos_NE2,pos_CE1,1.01)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA','HB1','HB2','HD2','HE1','HE2'])
    coord_copy.extend([HA,HB1,HB2,HD2,HE1,HE2])
    return coord_copy, atoms_copy

def proline_h(ires,i="",nextres=None):
    """
    This function adds missing hydrogen
    atoms of the PRO residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_CG     = coord[atoms.index('CG')]
    pos_CD     = coord[atoms.index('CD')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
        HG1,HG2        =  LinearAlgebra.class2(pos_CB,pos_CD,pos_CG)
        HD1,HD2        =  LinearAlgebra.class2(pos_CG,pos_N,pos_CD)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA', 'HB1', 'HB2', 'HG1', 'HG2', 'HD1', 'HD2'])
        coord_copy.extend([HA,HB1,HB2,HG1,HG2,HD1,HD2])
        return coord_copy,[H,i], atoms_copy


    HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
    HG1,HG2        =  LinearAlgebra.class2(pos_CB,pos_CD,pos_CG)
    HD1,HD2        =  LinearAlgebra.class2(pos_CG,pos_N,pos_CD)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA', 'HB1', 'HB2', 'HG1', 'HG2', 'HD1', 'HD2'])
    coord_copy.extend([HA,HB1,HB2,HG1,HG2,HD1,HD2])
    return coord_copy, atoms_copy


def metheonine_h(ires,i="",nextres=None):
    """
    This function adds missing hydrogen
    atoms of the MET residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_CG     = coord[atoms.index('CG')]
    pos_SD     = coord[atoms.index('SD')]
    pos_CE     = coord[atoms.index('CE')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
        HG1,HG2        =  LinearAlgebra.class2(pos_CB,pos_SD,pos_CG)
        HE1            =  LinearAlgebra.calcCoordinate(pos_CG,pos_SD,pos_CE,1.09,180.0,109.4)
        HE2            =  LinearAlgebra.calcCoordinate(pos_CG,pos_SD,pos_CE,1.09,-60.0,109.5)
        HE3            =  LinearAlgebra.calcCoordinate(pos_CG,pos_SD,pos_CE,1.09,60.0,109.5)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB1','HB2','HG1','HG2','HE1','HE2','HE3'])
        coord_copy.extend([HA,HB1,HB2,HG1,HG2,HE1,HE2,HE3])
        return coord_copy, [H,i], atoms_copy

    HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
    HG1,HG2        =  LinearAlgebra.class2(pos_CB,pos_SD,pos_CG)
    HE1            =  LinearAlgebra.calcCoordinate(pos_CG,pos_SD,pos_CE,1.09,180.0,109.4)
    HE2            =  LinearAlgebra.calcCoordinate(pos_CG,pos_SD,pos_CE,1.09,-60.0,109.5)
    HE3            =  LinearAlgebra.calcCoordinate(pos_CG,pos_SD,pos_CE,1.09,60.0,109.5)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA','HB1','HB2','HG1','HG2','HE1','HE2','HE3'])
    coord_copy.extend([HA,HB1,HB2,HG1,HG2,HE1,HE2,HE3])
    return coord_copy, atoms_copy

def isDisulfide(ires,i,nextres=None):
    """
    This function adds missing hydrogen
    atoms of the CYS (Cystine) residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_SG     = coord[atoms.index('SG')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_SG,pos_CB)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB1','HB2'])
        coord_copy.extend([HA,HB1,HB2])
        return coord_copy, [H,i], atoms_copy

    else:
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_SG,pos_CB)
        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB1','HB2'])
        coord_copy.extend([HA,HB1,HB2])
        return coord_copy, atoms_copy

def notDisulfide(ires,i,resNo,atoms2,atmpos2,resNo2,resseq2,nextres=None):
    """
    This function adds missing hydrogen
    atoms of the CYS (Cysteine) residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_SG     = coord[atoms.index('SG')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HG             =  LinearAlgebra.cys_HG(pos_SG, pos_CB, pos_CA)
        HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_SG,pos_CB)
        di             =  LinearAlgebra.calcTorsionAngle(pos_CA,pos_CB,pos_SG,HG)
        xyz            =  [pos_CA,pos_CB,pos_SG,HG,1.34,di,0.19,0.11,0.07]
        HG            =  PotentialEnergy.optmizeH(resNo,atoms2,atmpos2,resNo2,resseq2,xyz)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB1','HB2','HG'])
        coord_copy.extend([HA,HB1,HB2,HG])
        return coord_copy, [H,i], atoms_copy

    else:
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HG             =  LinearAlgebra.cys_HG(pos_SG, pos_CB, pos_CA)
        HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_SG,pos_CB)
        di             =  LinearAlgebra.calcTorsionAngle(pos_CA,pos_CB,pos_SG,HG)
        xyz            =  [pos_CA,pos_CB,pos_SG,HG,1.34,di,0.19,0.11,0.07]
        HG            =  PotentialEnergy.optmizeH(resNo,atoms2,atmpos2,resNo2,resseq2,xyz)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB1','HB2','HG'])
        coord_copy.extend([HA,HB1,HB2,HG])
        return coord_copy, atoms_copy

def threonine_h(ires,i,resNo,atoms2,atmpos2,resNo2,resseq2,nextres=None):
    """
    This function adds missing hydrogen
    atoms of the THR residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_CG2    = coord[atoms.index('CG2')]
    pos_OG1    = coord[atoms.index('OG1')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HG1            =  LinearAlgebra.thr_HG1(pos_OG1,pos_CB,pos_CG2)
        HB             =  LinearAlgebra.class3(pos_CA,pos_OG1,pos_CB)
        HG12           =  LinearAlgebra.calcCoordinate(pos_OG1,pos_CB,pos_CG2,1.09,60.5,109.4)
        HG22           =  LinearAlgebra.calcCoordinate(pos_OG1,pos_CB,pos_CG2,1.09,-179.5,109.5)
        HG23           =  LinearAlgebra.calcCoordinate(pos_OG1,pos_CB,pos_CG2,1.09,-59.5,109.5)
        di             =  LinearAlgebra.calcTorsionAngle(pos_CA,pos_CB,pos_OG1,HG1)
        xyz            =  [pos_CA,pos_CB,pos_OG1,HG1,0.96,di,0.41,0.0,0.0]
        HG1            = PotentialEnergy.optmizeH(resNo,atoms2,atmpos2,resNo2,resseq2,xyz)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HG1','HA','HB','1HG2','2HG2','3HG2'])
        coord_copy.extend([HG1,HA,HB,HG12,HG22,HG23])
        return coord_copy,[H,i], atoms_copy


    HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HG1            =  LinearAlgebra.thr_HG1(pos_OG1,pos_CB,pos_CG2)
    HB             =  LinearAlgebra.class3(pos_CA,pos_OG1,pos_CB)
    HG12           =  LinearAlgebra.calcCoordinate(pos_OG1,pos_CB,pos_CG2,1.09,60.5,109.4)
    HG22           =  LinearAlgebra.calcCoordinate(pos_OG1,pos_CB,pos_CG2,1.09,-179.5,109.5)
    HG23           =  LinearAlgebra.calcCoordinate(pos_OG1,pos_CB,pos_CG2,1.09,-59.5,109.5)
    di             =  LinearAlgebra.calcTorsionAngle(pos_CA,pos_CB,pos_OG1,HG1)
    xyz            =  [pos_CA,pos_CB,pos_OG1,HG1,0.96,di,0.41,0.0,0.0]
    HG1            =  PotentialEnergy.optmizeH(resNo,atoms2,atmpos2,resNo2,resseq2,xyz)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HG1','HA','HB','1HG2','2HG2','3HG2'])
    coord_copy.extend([HG1,HA,HB,HG12,HG22,HG23])
    return coord_copy, atoms_copy


def serine_h(ires,i,resNo,atoms2,atmpos2,resNo2,resseq2,nextres=None):
    """
    This function adds missing hydrogen
    atoms of the SER residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_OG     = coord[atoms.index('OG')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_OG,pos_CB)
        HG             =  LinearAlgebra.ser_HG(pos_OG,pos_CB,pos_CA)
        di             =  LinearAlgebra.calcTorsionAngle(pos_CA,pos_CB,pos_OG,HG)
        xyz            =  [pos_CA,pos_CB,pos_OG,HG,0.96,di,0.41,0.0,0.0]
        HG             =  PotentialEnergy.optmizeH(resNo,atoms2,atmpos2,resNo2,resseq2,xyz)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB1','HB2','HG'])
        coord_copy.extend([HA,HB1,HB2,HG])
        return coord_copy,[H,i], atoms_copy

    HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_OG,pos_CB)
    HG             =  LinearAlgebra.ser_HG(pos_OG,pos_CB,pos_CA)
    di             =  LinearAlgebra.calcTorsionAngle(pos_CA,pos_CB,pos_OG,HG)
    xyz            =  [pos_CA,pos_CB,pos_OG,HG,0.96,di,0.41,0.0,0.0]
    HG             =  PotentialEnergy.optmizeH(resNo,atoms2,atmpos2,resNo2,resseq2,xyz)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA','HB1','HB2','HG'])
    coord_copy.extend([HA,HB1,HB2,HG])
    return coord_copy, atoms_copy


def tyrosine_h(ires,i,resNo,atoms2,atmpos2,resNo2,resseq2,nextres=None):
    """
    This function adds missing hydrogen
    atoms of the TYR residue.
    See arginine() for full annotation
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    pos_O      = coord[atoms.index('O')]
    pos_CB     = coord[atoms.index('CB')]
    pos_CG     = coord[atoms.index('CG')]
    pos_CD1    = coord[atoms.index('CD1')]
    pos_CD2    = coord[atoms.index('CD2')]
    pos_CE1    = coord[atoms.index('CE1')]
    pos_CE2    = coord[atoms.index('CE2')]
    pos_CZ     = coord[atoms.index('CZ')]
    pos_OH     = coord[atoms.index('OH')]

    if nextres !=None:
        nextres_CA = nextres[2][1]
        nextres_N  = nextres[2][0]
        H              =  LinearAlgebra.class5(nextres_CA,nextres_N,pos_C,1.01)
        HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
        HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
        HD1            =  LinearAlgebra.class5(pos_CG,pos_CD1,pos_CE1,1.08)
        HD2            =  LinearAlgebra.class5(pos_CE2,pos_CD2,pos_CG,1.08)
        HE1            =  LinearAlgebra.class5(pos_CZ,pos_CE1,pos_CD1,1.08)
        HE2            =  LinearAlgebra.class5(pos_CZ,pos_CE2,pos_CD2,1.08)
        HH             =  LinearAlgebra.tyr_HH(pos_OH,pos_CZ,pos_CE2)
        di             =  LinearAlgebra.calcTorsionAngle(pos_CE2,pos_CZ,pos_OH,HH)
        xyz            =  [pos_CE2,pos_CZ,pos_OH,HH,0.96,di,0.37,0.0,0.0]
        HH             =  PotentialEnergy.optmizeH(resNo,atoms2,atmpos2,resNo2,resseq2,xyz)

        atoms_copy = copy.deepcopy(atoms)
        coord_copy = copy.deepcopy(coord)
        atoms_copy.extend(['HA','HB1','HB2','HD1','HD2','HE1','HE2','HH'])
        coord_copy.extend([HA,HB1,HB2,HD1,HD2,HE1,HE2,HH])
        return coord_copy,[H,i], atoms_copy

    HA             =  LinearAlgebra.class3(pos_CB,pos_N,pos_CA)
    HB1,HB2        =  LinearAlgebra.class2(pos_CA,pos_CG,pos_CB)
    HD1            =  LinearAlgebra.class5(pos_CG,pos_CD1,pos_CE1,1.08)
    HD2            =  LinearAlgebra.class5(pos_CE2,pos_CD2,pos_CG,1.08)
    HE1            =  LinearAlgebra.class5(pos_CZ,pos_CE1,pos_CD1,1.08)
    HE2            =  LinearAlgebra.class5(pos_CZ,pos_CE2,pos_CD2,1.08)
    HH             =  LinearAlgebra.tyr_HH(pos_OH,pos_CZ,pos_CE2)
    di             =  LinearAlgebra.calcTorsionAngle(pos_CE2,pos_CZ,pos_OH,HH)
    xyz            =  [pos_CE2,pos_CZ,pos_OH,HH,0.96,di,0.37,0.0,0.0]
    HH             =  PotentialEnergy.optmizeH(resNo,atoms2,atmpos2,resNo2,resseq2,xyz)

    atoms_copy = copy.deepcopy(atoms)
    coord_copy = copy.deepcopy(coord)
    atoms_copy.extend(['HA','HB1','HB2','HD1','HD2','HE1','HE2','HH'])
    coord_copy.extend([HA,HB1,HB2,HD1,HD2,HE1,HE2,HH])
    return coord_copy, atoms_copy

def checkDisulfide(sg_i,coord):
    """
    This function checks whether CYS is cystine or cysteine.

    Arguments
    ----------
    sg_i:    the 3D coordinate of the CYS SG being investigated
    coord:   the coordinate of all CYSs

    Returns
    -------
    A string: if cystine, return it is a disulfide bond else return it is not
    """
    coord = [i for i in coord if i != sg_i]
    for i in coord:
        r = sqrt(sum([(i[k]- sg_i[k])**2 for k in range(3)]))
        if r <= 3.0:
            return 'is_bond'

    return 'not_bond'

def ntermini_notpro(ires,resn2):
    """
    This function will generate N-ter hydrogens if it is not a PRO.

    Arguments
    ----------
    ires:    a list of list containing the atom names and coordinates
    resn2:   a list containing the residue names

    Returns
    -------
    A list: a list containing the coordinates of the H-atoms
    """
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_C      = coord[atoms.index('C')]
    try:
        pos_CB = coord[atoms.index('CB')]
    except:
        pass

    if resn2 != 'GLY':
        H1 = LinearAlgebra.calcCoordinate(pos_CB,pos_CA,pos_N,1.01,179.6,109.5)
        H2 = LinearAlgebra.calcCoordinate(pos_CB,pos_CA,pos_N,1.01,-60.4,109.6)
        H3 = LinearAlgebra.calcCoordinate(pos_CB,pos_CA,pos_N,1.01,60.2,109.6)

    elif resn2 == 'GLY':
        H1 = LinearAlgebra.calcCoordinate(pos_C,pos_CA,pos_N,1.01,179.6,109.5)
        H2 = LinearAlgebra.calcCoordinate(pos_C,pos_CA,pos_N,1.01,-60.4,109.6)
        H3 = LinearAlgebra.calcCoordinate(pos_C,pos_CA,pos_N,1.01,60.2,109.6)

    return [H1,H2,H3]

def ntermini_pro(ires,resn2):
    """
    This function will generate N-ter hydrogens if it is a PRO.

    Arguments
    ----------
    ires:    a list of list containing the atom names and coordinates
    resn2:   a list containing the residue names

    Returns
    -------
    A list: a list containing the coordinates of the H-atoms
    """
    atoms      = ires[1]
    coord      = ires[2]
    atoms,coord= ires[1], ires[2]
    pos_N      = coord[atoms.index('N')]
    pos_CA     = coord[atoms.index('CA')]
    pos_CD     = coord[atoms.index('CD')]

    H1,H2 = LinearAlgebra.class2(pos_CA,pos_CD,pos_N,1.01)

    return [H1,H2]

def his_note(resNo, chainNo):
    """
    This function informs user
    whenever HIS is protonated
    on both delta and epsilon N

    Arguments
    ----------
    resNo:   the residue number

    chainNo: the chain number

    Returns
    -------
    None:   writes to a log.txt file

    """
    if resNo[-1].isalpha():
        resNo = resNo[:-1]
    with open('log.txt', 'a') as f:
        f.write('HIS {} in chain {} is protonated (+1 charge)'\
        ' (ref @ www.protein-science.com)'
        .format(resNo,chainNo)+'\n\n')

def get_hd1(ires):
    """
    This function adds hydrogen atom to ND1 of
    histidine making it to have +1 charge.
    input argument has the same meaning as previous
    """
    atoms,coord= ires[1], ires[2]
    pos_CG     = coord[atoms.index('CG')]
    pos_ND1    = coord[atoms.index('ND1')]
    pos_CE1    = coord[atoms.index('CE1')]

    HD1        =  LinearAlgebra.class5(pos_CE1,pos_ND1,pos_CG,1.01)

    return [HD1]

def prot_his(x,res_pos,atom_name,chainNo):
    """
    This function checks if there are more than 4 HIS resi in order to add
    hydrogen atom to ND1 of histidine making it to have +1 charge.

    ----------
    x:         a list containing all residue name,atoms,atom coordinates

    res_pos:   a list containing all residue atom positions

    atom_name: a list containing all residue atom names

    chainNo:   the chain number of the chain

    Returns
    -------
    A list:     index 0 = updated res_pos and index 1 = updated atom_name
    """
    his  = [j for i,j in enumerate(x) if x[i][0][:3] == 'HIS']
    num  = [res_pos[i] for i,j in enumerate(x) if x[i][0][:3] == 'HIS']
    name = [atom_name[i] for i,j in enumerate(x) if x[i][0][:3] == 'HIS']

    if len(his) <= 4:
        pass
    else:
        prot = len(his)//5

        for i in range(prot):
            res_pos[res_pos.index(num[i])].extend(get_hd1(his[i]))
            atom_name[atom_name.index(name[i])].extend(['HD1'])
            his_note(his[i][0][3:],chainNo+1)

    return [res_pos, atom_name]
