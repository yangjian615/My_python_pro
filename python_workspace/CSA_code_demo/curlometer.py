#!/usr/bin/env python
# coding=utf-8

import numpy as np
import math

def delta(ref, i):
    delrefi = i - ref
    return delrefi

def curlometer(d1, d2, d3, d4):
    
    km2m = 1e3
    nT2T = 1e-9
    mu0 = (4*math.pi)*1e-7
    
    C1R = np.array([d1['x'], d1['y'], d1['z']])*km2m
    C1B = np.array([d1['Bx'], d1['By'], d1['Bz']])*nT2T
    C2R = np.array([d2['x'], d2['y'], d2['z']])*km2m
    C2B = np.array([d2['Bx'], d2['By'], d2['Bz']])*nT2T
    C3R = np.array([d3['x'], d3['y'], d3['z']])*km2m
    C3B = np.array([d3['Bx'], d3['By'], d3['Bz']])*nT2T
    C4R = np.array([d4['x'], d4['y'], d4['z']])*km2m
    C4B = np.array([d4['Bx'], d4['By'], d4['Bz']])*nT2T
    
    delB14 = delta(C4B, C1B)
    delB24 = delta(C4B, C2B)
    delB34 = delta(C4B, C3B)
    delR14 = delta(C4R, C1R)
    delR24 = delta(C4R, C2R)
    delR34 = delta(C4R, C3R)

# J

    # Have to 'convert' this to a matrix to be able to get the inverse.
    R = np.matrix(([np.cross(delR14, delR24), np.cross(delR24, delR34),\
         np.cross(delR14, delR34)]))
    Rinv = R.I

    # I(average) matrix:
    Iave = ([np.dot(delB14, delR24) - np.dot(delB24, delR14)],\
        [np.dot(delB24, delR34) - np.dot(delB34, delR24)],\
        [np.dot(delB14, delR34) - np.dot(delB34, delR14)])

    JJ = (Rinv*Iave)/mu0
                  
# div B
    lhs = np.dot(delR14, np.cross(delR24, delR34))

    rhs = np.dot(delB14, np.cross(delR24, delR34)) + \
        np.dot(delB24, np.cross(delR34, delR14)) + \
        np.dot(delB34, np.cross(delR14, delR24))

    divB = rhs/lhs

# div B / curl B
    curlB = JJ*mu0
    magcurlB = math.sqrt(curlB[0]**2 + curlB[1]**2 + curlB[2]**2)
    divBbycurlB = divB/magcurlB

    return [JJ, divB, divBbycurlB]