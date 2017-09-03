# region Description
# test curlometer
#

import numpy as np
import math


def delta(ref, i):
    delrefi = i - ref
    return delrefi


# Constants and conversions
mu0 = (4 * math.pi) * 1e-7
km2m = 1e3
nT2T = 1e-9
# Data example in GSE (2001-06-11T20:05:00.100Z)
C1B = np.array([-0.074, -4.924, -1.178]) * nT2T
C1R = np.array([-17048.9, -121511.3, 12814.3]) * km2m
C2B = np.array([0.575, -0.864, 0.440]) * nT2T
C2R = np.array([-18907.9, -122063.9, 13353.4]) * km2m
C3B = np.array([-0.239, -2.992, 0.959]) * nT2T
C3R = np.array([-18363.2, -121961.3, 11480.5]) * km2m
C4B = np.array([0.670, -3.871, -0.616]) * nT2T
C4R = np.array([-17470.0, -123201.5, 12519.8]) * km2m
delR14 = delta(C4R, C1R)
delR24 = delta(C4R, C2R)
delR34 = delta(C4R, C3R)
delB14 = delta(C4B, C1B)
delB24 = delta(C4B, C2B)
delB34 = delta(C4B, C3B)

# Calculate J

# Have to °Øconvert°Ø this to a matrix to be able to get the inverse.
R = np.matrix(([np.cross(delR14, delR24), np.cross(delR24, delR34),
                np.cross(delR14, delR34)]))
# The inverse:
Rinv = R.I

# I(average) matrix (note the shape):
Iave = ([np.dot(delB14, delR24) - np.dot(delB24, delR14)],
        [np.dot(delB24, delR34) - np.dot(delB34, delR24)],
        [np.dot(delB14, delR34) - np.dot(delB34, delR14)])
# The dot is equivalent to * here
JJ = np.dot(Rinv, Iave) / mu0
print(JJ)
# Calculate div B
lhs = np.dot(delR14, np.cross(delR24, delR34))
rhs = np.dot(delB14, np.cross(delR24, delR34)) + \
      np.dot(delB24, np.cross(delR34, delR14)) + \
      np.dot(delB34, np.cross(delR14, delR24))
divB = rhs / lhs
print(divB)
curlB = JJ * mu0
magcurlB = math.sqrt(curlB[0] ** 2 + curlB[1] ** 2 + curlB[2] ** 2)
divBbycurlB = abs(divB) / magcurlB
print(divBbycurlB)
# endregion





# region Description

import numpy as np  # for matrix calculations
import math  # for pi and sqrt
import glob  # for sensible listdir()
import spacepy.pycdf as pycdf  # for reading CDF files
from copy import deepcopy  # for obtaining variables in CEF files import matplotlib.pyplot as plt # for plotting
import datetime as dt  # for dates
from matplotlib import dates  # for formatting axes
import sys
sys.path.append('/Users/yangjian/CEFLIB/PYTHON')  #‘ÿ»Îceflib.so
from ceflib import *

import spacepy.time as spt
import datetime as dt



# User-defined variables:

# Path to data:
path = '/Users/yangjian/Desktop/case_20010731_201417_201617_overview/CAA/'
# Filename for output current density:
outfile = 'test_J.txt'

# Plot filenames:
BJQFileName = 'test_BJQ.png'
GeomFileName = 'test_Geom.png'
# X-axis labels:
XAxisLabel = 'Time on 11th June 2001'
# Desired resolution of data for the curlometer calculation
window = 0.2  # in seconds, 0.2 = minimum


def delta(ref, i):
    delrefi = i - ref
    return delrefi


def curlometer(d1, d2, d3, d4):


    km2m = 1e3
    nT2T = 1e-9
    mu0 = (4 * math.pi) * 1e-7
    C1R = np.array([d1[3], d1[4], d1[5]]) * km2m
    C1B = np.array([d1[0], d1[1], d1[2]]) * nT2T
    C2R = np.array([d2[3], d2[4], d2[5]]) * km2m
    C2B = np.array([d2[0], d2[1], d2[2]]) * nT2T
    C3R = np.array([d3[3], d3[4], d3[5]]) * km2m
    C3B = np.array([d3[0], d3[1], d3[2]]) * nT2T
    C4R = np.array([d4[3], d4[4], d4[5]]) * km2m
    C4B = np.array([d4[0], d4[1], d4[2]]) * nT2T

    delB14 = delta(C4B, C1B)
    delB24 = delta(C4B, C2B)
    delB34 = delta(C4B, C3B)
    delR14 = delta(C4R, C1R)
    delR24 = delta(C4R, C2R)
    delR34 = delta(C4R, C3R)

    # J

    # Have to 'convert' this to a matrix to be able to get the inverse.
    R = np.matrix(([np.cross(delR14, delR24), np.cross(delR24, delR34),
                    np.cross(delR14, delR34)]))
    Rinv = R.I
    # I(average) matrix:
    Iave = ([np.dot(delB14, delR24) - np.dot(delB24, delR14)],
            [np.dot(delB24, delR34) - np.dot(delB34, delR24)],
            [np.dot(delB14, delR34) - np.dot(delB34, delR14)])

    JJ = (Rinv * Iave) / mu0

    # div B
    lhs = np.dot(delR14, np.cross(delR24, delR34))
    rhs = np.dot(delB14, np.cross(delR24, delR34)) + \
          np.dot(delB24, np.cross(delR34, delR14)) + \
          np.dot(delB34, np.cross(delR14, delR24))
    divB = abs(rhs) / abs(lhs)
    # div B / curl B
    curlB = JJ * mu0
    magcurlB = math.sqrt(curlB[0] ** 2 + curlB[1] ** 2 + curlB[2] ** 2)
    divBbycurlB = divB / magcurlB

    return [JJ, divB, divBbycurlB]
# End of curlometer function

'''Read in all the data using pycdf.read'''

cluster = ['C' + str(x) for x in range(1, 5)]

time = {}
B = {}
pos = {}

for c in cluster:
    folder = c + '_CP_FGM_5VPS/*.cdf'
    filename = glob.glob(path + folder)
    print(filename)

    cdf = pycdf.CDF(filename[0])
    time[c] = cdf['time_tags__'+c+'_CP_FGM_5VPS'][...] # in milli-seconds
    B[c] = cdf['B_vec_xyz_gse__'+c+'_CP_FGM_5VPS'][...]
    pos[c] = cdf['sc_pos_xyz_gse__'+c+'_CP_FGM_5VPS'][...]



'''Align all the data with the time by using a dictionary with
the time in milliseconds as the key'''

clean = {}
for c in cluster:

    for i, p in enumerate(time[c]):

        if p not in clean.keys():
            clean[int(p)] = {}

        clean[p][c] = [B[c][i][0],\
                       B[c][i][1],\
                       B[c][i][2],\
                       pos[c][i][0],\
                       pos[c][i][1],\
                       pos[c][i][2]]

mintime,maxtime = min(clean.keys()), max(clean.keys())
# Time array (min, max, step)
tarr = range(mintime, maxtime, int(window * 1000))
nwin = len(tarr)

Jave = np.zeros(nwin, dtype=[('time', float), ('Jx', float),\
                                     ('Jy', float), ('Jz', float),\
                                     ('divB', float),\
                                     ('divBcurlB', float)])

for i, t in enumerate(tarr):

            if len(clean[t]) == 4:
                onej = curlometer(clean[t]['C1'], clean[t]['C2'],\
                                     clean[t]['C3'], clean[t]['C4'])

                Jave['time'][i] = t / 1000
                Jave['Jx'][i] = onej[0][0]

                Jave['Jy'][i] = onej[0][1]
                Jave['Jz'][i] = onej[0][2]

                Jave['divB'][i] = onej[1]
                Jave['divBcurlB'][i] = onej[2]
            else:
                Jave['time'][i] = t / 1000

                Jave['Jx'][i] = np.nan
                Jave['Jy'][i] = np.nan

                Jave['Jz'][i] = np.nan
                Jave['divB'][i] = np.nan

                Jave['divBcurlB'][i] = np.nan
'''Write all results out to file, tarr is already sorted'''

with open(outfile, 'w') as f:
    for j in Jave:
        outstring = str(dt.datetime.utcfromtimestamp(j['time'])) + \
                            ',' + str(j['Jx']) + ',' + str(j['Jy']) + \
                            ',' + str(j['Jz']) + ',' + str(j['divBcurlB']) + '\n'
        f.write(outstring)

'''Pull out the mag field used for the calculation'''

Magnpt = {}
for c in cluster:
    Bx, By, Bz, Bmag = [], [], [], []
    for p in tarr:
        if c in clean[p].keys():
            Bx.append(clean[p][c][0])
            By.append(clean[p][c][1])
            Bz.append(clean[p][c][2])
            Bmag.append(math.sqrt(clean[p][c][0] ** 2 + clean[p][c][1] ** 2 + clean[p][c][2] ** 2))
        else:
                Bx.append(np.nan)
                By.append(np.nan)
                Bz.append(np.nan)
                Bmag.append(np.nan)
        Magnpt[c] = [Bx, By, Bz, Bmag]

'''Take times and put as date into list'''




# endregion
