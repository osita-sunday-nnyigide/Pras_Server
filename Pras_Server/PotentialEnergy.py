#!/usr/bin/env python

__doc__ = """

This program requires python 3.6 or higher.

This module, PotentialEnergy.py has the function(s)

that is used to calculate the potential energy between

class6 hydrogen atoms with rotational freedom and the

rest of protein heavy atoms.

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
from .ForcefieldParam import param
from math import sqrt, acos,pi, degrees

def calcDistance(m, n):
    """
    This function calculates the distance in nm between two points

    Arguments
    ----------
    m:   the first 3D point
    n:   the second 3D point

    Returns
    -------
    Float: the calculated distance
    """

    return (sqrt(sum([(m[z] - n[z])**2 for z in range(3)])))/10

def recalcCoordinate(b,c,d,bond_len):
    """
    This function is used to re-calculate the class6
    H-atom after rotation of the dihedral angle

    Arguments
    ----------
    b:        the second 3D point
    c:        the third  3D point
    d:        the 4th 3D point
    bond_len: the bond length of c and d

    Returns
    -------
    A list: the calculated coordinate
    """

    b,c,d = np.array(b),np.array(c),np.array(d)
    theta = 107.0

    u1 = b - c
    y1 = d - c

    mag_y1  = sqrt((y1[0]**2)+(y1[1]**2)+(y1[2]**2))
    mag_u1  = sqrt((u1[0]**2)+(u1[1]**2)+(u1[2]**2))

    theta_bcd = acos(np.dot(u1,y1)/(mag_u1*mag_y1))
    rotate    = theta - degrees(theta_bcd)

    z  = np.cross(u1, y1)
    n  = z/sqrt((z[0]**2)+(z[1]**2)+(z[2]**2))

    pos_ini = c + y1*np.cos(rotate*pi/180) \
              +(np.cross(n,y1)*np.sin(rotate*pi/180))\
              + n*(np.dot(n,y1))*(1-np.cos(rotate*pi/180))

    pos_BL = ((pos_ini-c)*(0.96/np.linalg.norm(pos_ini-c)))+c

    return pos_BL.tolist()

def calcCoordinate(a,b,c,bond_len,di_angle,theta):
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

    d = (b + (q+e))

    u1 = b - c
    y1 = d - c

    mag_y1  = sqrt((y1[0]**2)+(y1[1]**2)+(y1[2]**2))
    mag_u1  = sqrt((u1[0]**2)+(u1[1]**2)+(u1[2]**2))

    theta_bcd = acos(np.dot(u1,y1)/(mag_u1*mag_y1))
    rotate    = theta - degrees(theta_bcd)

    z  = np.cross(u1, y1)
    n  = z/sqrt((z[0]**2)+(z[1]**2)+(z[2]**2))

    pos_ini = c + y1*np.cos(rotate*pi/180) \
              +(np.cross(n,y1)*np.sin(rotate*pi/180))\
              + n*(np.dot(n,y1))*(1-np.cos(rotate*pi/180))

    pos_BL = ((pos_ini-c)*(bond_len/np.linalg.norm(pos_ini-c)))+c

    return pos_BL.tolist()

def _optmize(pos_CA,pos_CB,pos_OG1,pos_HG1,BL,di):
    """
    This function rotates class6 H-atoms through 360 deg

    Arguments (variables do not always correspond to atom names)
    ----------
    pos_OG1: the coordinate of the first atom
    pos_CB:  the coordinate of the second atom
    pos_CA:  the coordinate of the third atom
    pos_HG1: the coordinate of the H-atom
    BL     : the bond length of the H-atom
    di     : the dihedral angle that defines the H-atom

    Returns
    -------
    A list: the coordinate corresponding to H after each rotation
    """

    CB_OG1_HG1_angle = 109.5
    opt_h = []

    for i in range(0, 360, 60):
        point_h = calcCoordinate(pos_CA,pos_CB,pos_OG1,BL,CB_OG1_HG1_angle,di)
        point_h = recalcCoordinate(pos_CB,pos_OG1,point_h,BL)
        opt_h.extend([point_h])
        di+= i

    return opt_h

def optmizeH(num,atom,atmpos,resNo,resn,data):
    """
    For each rotation above, this funtion calculates
    the potential energy b/w the H-atom in question
    and other heavy atoms within the interaction distance.

    Arguments
    ----------
    num:    the residue number of the residue
    atom:   a list of the atoms of each residue in the chain
    atmpos: a list of the coordinates of each atom of each residue
    resNo:  a list of all residue numbers in the chain
    resn:   a list of all residue names in the chain
    data:   a list of coordinates of the H-atom and other 3 atoms,
            H-atom bond length and dihedral angle that defines it,
            partial charge, sigma and epsilon

    Returns
    -------
    A list: the coordinate corresponding to the lowest potential energy
    """

    total_ener = []
    opt_h = _optmize(data[0],data[1],data[2],data[3],data[4],data[5])
    for z in range(len(resNo)):
        if resNo[z] ==num:
            _pos  = [pos for i, pos in enumerate(atmpos) if i != z]
            _name = [atm for i, atm in enumerate(atom) if i   != z]
            _res  = [res for i, res in enumerate(resn) if i   != z]
            for h in opt_h:
                poten_ener = 0
                for k,l in enumerate(_pos):
                    for n,m in enumerate(_pos[k]):
                        if sqrt(sum([(h[q]-m[q])**2 for q in range(3)])) <= 2.5:
                            r_die     = 138.94
                            H_pcharge = data[6]
                            A_pcharge = param[_res[k]][_name[k][n]][0]
                            """
                            for the AMBER forcefield the sigma and epsilon
                            of class6 H of SER, TYR and THR = 0 and that of CYS is
                            small compared to heavy atoms so the VDW portion of this
                            equation can be commented out. However, it has no observable
                            effect on computaional speed. So it is allowed to stay incase
                            one decides to use a different forcefield (e.g., CHARMM)
                            that has values for them. In that case, user should change the
                            data in ForcefieldParam.py
                            """
                            H_sigma,H_epsilon = data[7],data[8]
                            A_sigma   = param[_res[k]][_name[k][n]][1]
                            A_epsilon = param[_res[k]][_name[k][n]][2]
                            r         = calcDistance(h,_pos[k][n])
                            coulomb   = (r_die*(H_pcharge*A_pcharge))/r
                            vdw_eps   = 4*(sqrt(H_epsilon*A_epsilon))
                            vdw_sigm  = (((H_sigma+A_sigma)*0.5/r)**12)\
                            -(((H_sigma+A_sigma)*0.5/r)**6)
                            poten_ener+= coulomb+(vdw_eps*vdw_sigm)
                total_ener.extend([poten_ener])

    # if no acceptor is interacting,
    # return the initial H-coord
    try:
        min_energy = total_ener.index(min(total_ener))
        min_H_pos = opt_h[min_energy]
        return  min_H_pos
    except:
        return data[3]
