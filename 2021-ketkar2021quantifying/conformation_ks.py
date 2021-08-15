# This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#!/usr/bin/python

# Script:  conformation_ks.py
# Purpose: analyze conformation of polymers
# Syntax:  conformation_ks.py
# Author:  Y. Seo from Lisa's msd code & Youngmi's circling.py, modified by Kevin Shen, Jan, 2020

# derived from fortran code
# -------------------------------------------------------

import sys
import string
import numpy as np
from numpy import *
from math import *
import fileinput

# INPUT PARAMETERs
infiles = ['equil2.lammpstrj']
fT = 0.53
fA = 0.39
nbeads = 180
npoly = 330  # need this for com calc
nskip = 0
nconf = 201
nlayers = 3

# hntaper = int(fT*nbeads*0.5)
pureA = int(nbeads*(fA-fT/2))
pureB = int(nbeads*(1-fA-fT/2))
hntaper = (nbeads-pureA-pureB)/2
print('taper, hntaper, pureA, pureB:', fT, hntaper, pureA, pureB)
conf = -1
fold1 = [0]*nconf
fold2 = [0]*nconf
stretch1 = [0]*nconf
stretch2 = [0]*nconf
config = [[0 for x in range(nconf)] for x in range(npoly)]


def e2e():
    fb = 0
    lb = 0
    for jj in range(0, npoly):
        for k in range(jj*nbeads+1, jj*nbeads+nbeads+1):
            if k % nbeads == 0:  # last bead of the chain
                fb = zc[k]
            elif k == 1 or (k-1) % nbeads == 0:  # first bead of the chain
                lb = zc[k]
        square_e2ezdistance[jj] = ((fb-lb)**2)


def segcom():
    for jj in range(0, npoly):
        p1z = 0
        t1z = 0
        p2z = 0
        t2z = 0
        for k in range(jj*nbeads+1, jj*nbeads+1+pureA):
            p1z += zc[k]/pureA
        for k in range(jj*nbeads+1+pureA, jj*nbeads+1+pureA+hntaper):
            t1z += zc[k]/hntaper
        for k in range(jj*nbeads+1+pureA+hntaper, jj*nbeads+1+pureA+hntaper*2):
            t2z += zc[k]/hntaper
        for k in range(jj*nbeads+1+pureA+hntaper*2, jj*nbeads+1+pureA+pureB+hntaper*2):
            p2z += zc[k]/pureB
        t1zp1z[jj] = t1z - p1z
        p1zt1z[jj] = p1z - t1z
        t1zt2z[jj] = t1z - t2z
        t2zp2z[jj] = t2z - p2z


IN = fileinput.input(infiles)
natoms = 0

for loopnum in range(0, nskip):
    IN.readline()
    IN.readline()
    IN.readline()
    line = IN.readline()      # number of atoms
    fields = string.split(line)
    natoms = int(fields[0])
    IN.readline()
    IN.readline()
    IN.readline()
    IN.readline()
    IN.readline()
    for j in range(natoms):
        IN.readline()

for loopnum in range(0, nconf):  # read starting timestep

    IN.readline()
    line = IN.readline()      # time step
    fields = string.split(line)
    step = int(fields[0])
    IN.readline()
    line = IN.readline()      # number of atoms
    fields = string.split(line)

    # reads number of atoms from first configuration; don't support changing num atoms
    natoms = int(fields[0])
    dim = natoms+1
    if loopnum == 0:
        zc = np.zeros(dim, np.float32)
        nbead = np.zeros(npoly+1)
        zcm = np.zeros(npoly+1)
        zi = np.zeros(dim)
        zr = np.zeros(dim)
        typea = [0]*dim
        mol = [0]*dim
        square_e2ezdistance = [0]*npoly
        t1zp1z = [0]*npoly
        p1zt1z = [0]*npoly
        t1zt2z = [0]*npoly
        t2zp2z = [0]*npoly
        #
    IN.readline()
    line = IN.readline()
    [xm, xp] = map(float, line.split())
    line = IN.readline()
    [ym, yp] = map(float, line.split())
    line = IN.readline()
    [zm, zp] = map(float, line.split())
    line = IN.readline()
    zbox = zp - zm
    Lz = zbox/nlayers
    conf += 1
    for j in range(1, dim):
        line = IN.readline()
        [ii, molj, typej, q, x1, x2, x3, n1, n2, n3] = string.split(line)
        k = int(ii)
        typea[k] = int(typej)
        mol[k] = int(molj)
        zc[k] = zbox*float(x3) + int(n3)*zbox

    e2e()
    segcom()

    for jj in range(0, npoly):  # polymer categorization step
        if sqrt(square_e2ezdistance[jj]) > Lz:  # stretching2
            config[jj][conf] = 4
        elif t1zp1z[jj]*t2zp2z[jj] > 0:  # bridging
            config[jj][conf] = 3
        elif p1zt1z[jj]*t1zt2z[jj] > 0:  # stretching1
            config[jj][conf] = 1
        else:  # folding
            config[jj][conf] = 2

config2 = [[0 for x in range(6)] for x in range(npoly)]
for p in range(0, npoly):
    for t in range(0, nconf):
        config2[p][config[p][t]-1] += 1

for p in range(0, npoly):
    for t in range(0, nconf):
        config2[p][5] = max(config2[p])
        config2[p][4] = config2[p].index(max(config2[p])) + 1
    # print(" ".join(str(e) for e in config2[p]))

stretch1id = []
foldingid = []
bridgingid = []
stretch2id = []
for p in range(0, npoly):
    if config2[p][4] == 1:
        stretch1id.append(p+1)
    elif config2[p][4] == 2:
        foldingid.append(p+1)
    elif config2[p][4] == 3:
        bridgingid.append(p+1)
    else:
        stretch2id.append(p+1)

print('=================Conformation Summary=================')
print('stretch1, folding, bridging, stretch2')
print(len(stretch1id), len(foldingid), len(bridgingid), len(stretch2id))
print('stretch1:', stretch1id)
print('folding:', foldingid)
print('bridging:', bridgingid)
print('stretch2', stretch2id)
print('====================End of Summary====================')
