#!/usr/bin/env python

__doc__ = """

This program requires python 3.6 or higher.

This module, LinearAlgebra.py has all the

algebraic functions used to fix missing

hydrogen atoms.

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
from math import sqrt, cos, sin, acos,pi, degrees

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

    u1      = b - c
    y1      = d - c
    mag_y1  = sqrt((y1[0]**2)+(y1[1]**2)+(y1[2]**2))
    mag_u1  = sqrt((u1[0]**2)+(u1[1]**2)+(u1[2]**2))

    theta_bcd = acos(np.dot(u1,y1)/(mag_u1*mag_y1))
    rotate    = theta - degrees(theta_bcd)

    z  = np.cross(u1, y1)
    n  = z/sqrt((z[0]**2)+(z[1]**2)+(z[2]**2))

    pos_ini = c + y1*np.cos(rotate*pi/180) \
              +(np.cross(n,y1)*np.sin(rotate*pi/180))\
              + n*(np.dot(n,y1))*(1-np.cos(rotate*pi/180))

    pos_BL  = ((pos_ini-c)*(bond_len/np.linalg.norm(pos_ini-c)))+c

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

def rotMatrix(axis,theta):
    """
    This function generates 3D rotation matrix
    for a rotation about an axis.

    Arguments
    ----------
    axis:  the vector about which rotation is applied
    theta: the amount in degrees to rotate the vector

    Returns
    -------
    Nested list: the matrix
    """
    matrix = [[0. for j in range(3)] for i in range(3)]

    axis_length = sqrt((axis[0]**2 + axis[1]**2 + axis[2]**2))
    xNorm = axis[0]/axis_length
    yNorm = axis[1]/axis_length
    zNorm = axis[2]/axis_length

    sin_theta = sin(theta)
    cos_theta = cos(theta)
    one_costheta = 1.0 - cos_theta

    matrix[0][0] = cos_theta                + xNorm*xNorm*one_costheta
    matrix[0][1] = xNorm*yNorm*one_costheta - zNorm*sin_theta
    matrix[0][2] = xNorm*zNorm*one_costheta + yNorm*sin_theta
    matrix[1][0] = xNorm*yNorm*one_costheta + zNorm*sin_theta
    matrix[1][1] = cos_theta                + yNorm*yNorm*one_costheta
    matrix[1][2] = yNorm*zNorm*one_costheta - xNorm*sin_theta
    matrix[2][0] = xNorm*zNorm*one_costheta - yNorm*sin_theta
    matrix[2][1] = yNorm*zNorm*one_costheta + xNorm*sin_theta
    matrix[2][2] = cos_theta                + zNorm*zNorm*one_costheta

    return matrix

def arbRot(vector,axis,theta):
    """
    This function generates a new position vector after
    rotating a vector by an angle through an axis.

    Arguments
    ----------
    vector: the vector to rotate
    axis:   the vector about which rotation is applied
    theta:  the amount in degrees to rotate the vector

    Returns
    -------
    A list: the coordinates of the new position vector
    """
    matrix = rotMatrix(axis,theta)

    return [sum([vector[j]*matrix[j][i] for j in range(3)]) for i in range(3)]

def class2(pos_CA,pos_CG,pos_CB,bond_len=""):
    """
    Given a tetrahedron with two known vertices this function
    obtains the other vertices corresponding to class2 H-atoms.

    Arguments (variables do not always correspond to atom names)
    ----------
    pos_CA:   one of the vertices of the tetrahedron
    pos_CG:   one of the vertices of the tetrahedron
    pos_CB:   the center of the tetrahedron
    bond_len: bond length specified explicitly for H-atoms of PRO at N-ter

    Returns
    -------
    A list : index 0 the first H-atom coordinate
             index 1 the second H-atom coordinate
    """
    pos_CA,pos_CG,pos_CB = np.array(pos_CA),np.array(pos_CG),np.array(pos_CB)

    m = 0.5*(pos_CA+pos_CG)
    n = m + 2*(pos_CB-m)

    hb1  = n - (np.cross(pos_CA-m, pos_CB-m))/np.linalg.norm(pos_CB-m)
    hb2  = n + (np.cross(pos_CA-m, pos_CB-m))/np.linalg.norm(pos_CB-m)

    if bond_len:
        point_hb1 = ((hb1 - pos_CB)*(1.01/np.linalg.norm(hb1-pos_CB)))+pos_CB
        point_hb2 = ((hb2 - pos_CB)*(1.01/np.linalg.norm(hb2-pos_CB)))+pos_CB
        return point_hb1.tolist(), point_hb2.tolist()

    point_hb1 = ((hb1 - pos_CB)*(1.09/np.linalg.norm(hb1-pos_CB)))+pos_CB
    point_hb2 = ((hb2 - pos_CB)*(1.09/np.linalg.norm(hb2-pos_CB)))+pos_CB

    return point_hb1.tolist(), point_hb2.tolist()

def class3(pos_CB,pos_N,pos_CA):
    """
    Given a tetrahedron with three known vertices this function
    obtains the other vertex corresponding to a class3 H-atom.

    Arguments (variables do not always correspond to atom names)
    ----------
    pos_CB: one of the vertices of the tetrahedron
    pos_N:  one of the vertices of the tetrahedron
    pos_CA  the center of the tetrahedron

    Returns
    -------
    A list: the H-atom coordinate
    """
    pos_CB,pos_N,pos_CA = np.array(pos_CB),np.array(pos_N),np.array(pos_CA)

    m = 0.5*(pos_CB+pos_N)
    n = m + 2*(pos_CA-m)

    ha       = n + (np.cross(pos_CB-m, pos_CA-m))/np.linalg.norm(pos_CA-m)
    point_ha = ((ha - pos_CA)*(1.09/np.linalg.norm(ha-pos_CA)))+pos_CA

    return point_ha.tolist()

def class4(pos_NE,pos_CZ,pos_NH1):
    """
    This function translates vectors defining heavy atoms
    that form a trangle with class4 H-atoms to obtain coordinates
    of the H-atoms

    Arguments (variables do not always correspond to atom names)
    ----------
    pos_NE:  e.g. Arg NE attached to CZ that forms one vertex of the triangle
    pos_CZ:  the other vertex of the triangle
    pos_NH1: the center of the triangle

    Returns
    -------
    A list: a H-atom coordinate
    """
    pos_NE,pos_CZ,pos_NH1 = np.array(pos_NE),np.array(pos_CZ),np.array(pos_NH1)

    pos_hh11   = ((pos_NE- pos_CZ) + (pos_NH1 - pos_CZ))+pos_NH1
    point_1hh1 = ((pos_hh11 - pos_NH1)*(1.01/np.linalg.norm(pos_hh11 - pos_NH1)))+pos_NH1

    return point_1hh1.tolist()

def class5(pos_CD,pos_NE,pos_CZ,bond_len):
    """
    This function translates vectors defining heavy atoms
    that form a trangle with class5 H-atoms to obtain coordinates
    of the H-atom. Further check is performed (eqn. 3-5 in PRAS paper)
    to ensure the bond angle is accurate.

    Arguments (variables do not always correspond to atom names)
    ----------
    pos_CD:   one vertex of the triangle
    pos_NE:   the center of the triangle
    pos_CZ    the other vertex of the triangle
    bond_len: bond length as specified explicitly

    Returns
    -------
    A list: a H-atom coordinate
    """
    pos_CD,pos_NE,pos_CZ = np.array(pos_CD),np.array(pos_NE),np.array(pos_CZ)

    pos_he = ((pos_NE-pos_CD) + (pos_NE-pos_CZ))+pos_NE

    u1 = pos_CZ - pos_NE
    y1 = pos_he - pos_NE
    u2 = pos_CD - pos_NE

    mag_y1  = sqrt((y1[0]**2)+(y1[1]**2)+(y1[2]**2))
    mag_u1  = sqrt((u1[0]**2)+(u1[1]**2)+(u1[2]**2))
    mag_u2  = sqrt((u2[0]**2)+(u2[1]**2)+(u2[2]**2))

    theta_bcd = acos(np.dot(u1,y1)/(mag_u1*mag_y1))

    theta     = (degrees(acos(np.dot(u1,y1)/(mag_u1*mag_y1))\
                + acos(np.dot(u2,y1)/(mag_u2*mag_y1))))/2

    rotate    = theta - degrees(theta_bcd)

    z  = np.cross(u1, y1)
    n  = z/sqrt((z[0]**2)+(z[1]**2)+(z[2]**2))

    pos_ini = pos_NE + y1*np.cos(rotate*pi/180)\
              +(np.cross(n,y1)*np.sin(rotate*pi/180))\
              + n*(np.dot(n,y1))*(1-np.cos(rotate*pi/180))

    point_he = ((pos_ini-pos_NE)*(bond_len/np.linalg.norm(pos_ini-pos_NE)))+pos_NE

    return point_he.tolist()

def ser_HG(pos_OG,pos_CB,pos_CA):
    """
    This function  applies rotation to generate new coordinate for H-atom

    Arguments
    ----------
    pos_OG: the coordinate for the atom OG
    pos_CB: the coordinate for the atom CB
    pos_CA  the coordinate for the atom CA

    Returns
    -------
    A list: a H-atom coordinate
    """

    HB1,HB2 = class2(pos_CA,pos_OG,pos_CB) # will get HB1

    og_cb_vector  = [pos_CB[i] - pos_OG[i] for i in range(3)]
    cb_hb1_vector = [HB1[i] - pos_CB[i] for i in range(3)]

    rotation_amount = -240.2*(pi/180.)

    rotated = arbRot(og_cb_vector, cb_hb1_vector, rotation_amount)
    pos_h   = [rotated[i] + pos_OG[i] for i in range(3)]

    vec1  = [pos_h[i] - pos_OG[i] for i in range(3)]
    const = 0.96/sqrt(sum([(pos_h[i] - pos_OG[i])**2 for i in range(3)]))
    vec3  = [vec1[i]*const for i in range(3)]

    pos_hg = [vec3[i]+pos_OG[i] for i in range(3)]

    return pos_hg

def thr_HG1(pos_OG1,pos_CB,pos_CG2):
    """
    This function  applies rotation to generate new coordinate for H-atom

    Arguments
    ----------
    pos_OG1: the coordinate for the atom OG1
    pos_CB:  the coordinate for the atom CB
    pos_CG2  the coordinate for the atom CG2

    Returns
    -------
    A list: a H-atom coordinate
    """
    og1_cb_vector = [pos_CB[i] - pos_OG1[i] for i in range(3)]
    cb_cg2_vector = [pos_CG2[i] - pos_CB[i] for i in range(3)]

    rotation_amount = -243.2*(pi/180.)

    rotated = arbRot(og1_cb_vector, cb_cg2_vector, rotation_amount)
    pos_h   = [rotated[i] + pos_OG1[i] for i in range(3)]

    vec1  = [pos_h[i] - pos_OG1[i] for i in range(3)]
    const = 0.96/sqrt(sum([(pos_h[i] - pos_OG1[i])**2 for i in range(3)]))
    vec3  = [vec1[i]*const for i in range(3)]

    pos_hg1 = [vec3[i]+pos_OG1[i] for i in range(3)]

    return pos_hg1

def cys_HG(pos_SG,pos_CB,pos_CA):
    """
    This function  applies rotation to generate new coordinate for H-atom

    Arguments
    ----------
    pos_SG: the coordinate for the atom SG
    pos_CB: the coordinate for the atom CB
    pos_CA  the coordinate for the atom CA

    Returns
    -------
    A list: a H-atom coordinate
    """
    og1_cb_vector = [pos_CB[i] - pos_SG[i]for i in range(3)]
    cb_cg2_vector = [pos_CA[i] - pos_CB[i]for i in range(3)]

    rotation_amount = -243.2*(pi/180.)

    rotated = arbRot(og1_cb_vector, cb_cg2_vector, rotation_amount)
    pos_h   = [rotated[i] + pos_SG[i] for i in range(3)]

    vec1  = [pos_h[i] - pos_SG[i] for i in range(3)]
    const = 1.3/sqrt(sum([(pos_h[i] - pos_SG[i])**2 for i in range(3)]))
    vec3  = [vec1[i]*const for i in range(3)]

    pos_hg1 = [vec3[i]+pos_SG[i] for i in range(3)]

    return pos_hg1

def tyr_HH(pos_OH,pos_CZ,pos_CE2):
    """
    This function  applies rotation to generate new coordinate for H-atom

    Arguments
    ----------
    pos_OH:  the coordinate for the atom OH
    pos_CZ:  the coordinate for the atom CZ
    pos_CE2  the coordinate for the atom CE2

    Returns
    -------
    A list: a H-atom coordinate
    """
    oh_cz_vector  = [pos_CZ[i] - pos_OH[i] for i in range(3)]
    cz_ce2_vector = [pos_CE2[i] - pos_CZ[i] for i in range(3)]

    rotation_amount = -220.2*(pi/180.)

    rotated = arbRot(oh_cz_vector, cz_ce2_vector, rotation_amount)
    pos_h   = [rotated[i] + pos_OH[i]for i in range(3)]

    vec1  = [pos_h[i] - pos_OH[i] for i in range(3)]
    const = 0.96/sqrt(sum([(pos_h[i] - pos_OH[i])**2 for i in range(3)]))
    vec3  = [vec1[i]*const for i in range(3)]

    pos_hh = [vec3[i]+pos_OH[i] for i in range(3)]

    return pos_hh
