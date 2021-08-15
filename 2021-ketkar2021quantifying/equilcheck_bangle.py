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

# Script:   equilcheck_bangle.py
# Author:   Y. Seo from Lisa's msd code & Youngmi's circling.py
# Edited:   M. Fan Jan.2019, K.-H. Shen 2020
# python equilcheck_bangle.py

# derived from fortran code
# -------------------------------------------------------

import sys
import string
import fileinput
import numpy as np
import pandas as pd
from math import *

# INPUT PARAMETERs
NBEADS = 180
NLAYERS = 3
NSKIP = 1
NCONF = 200
NPOLY = 330
paras = NBEADS, NLAYERS, NSKIP, NCONF, NPOLY


def bangle(paras, xc, yc, zc, bond_angle):
    ''' Returns the unit vector of the vector.'''
    def unit_vector(vector):
        return vector / np.linalg.norm(vector)

    nbeads, nlayers, nskip, nconf, npoly = paras
    for jj in range(0, npoly):
        avg_angle = 0
        for k in range(jj*nbeads+1, jj*nbeads+nbeads+1):
            # not the first or last bead of the chain
            if k % nbeads != 0 or (k-1) % nbeads != 0:
                v1 = (xc[k]-xc[k-1], yc[k]-yc[k-1], zc[k]-zc[k-1])
                v2 = (xc[k+1]-xc[k], yc[k+1]-yc[k], zc[k+1]-zc[k])
                v1_u = unit_vector(v1)
                v2_u = unit_vector(v2)
                avg_angle += np.arccos(np.clip(np.dot(v1_u,
                                       v2_u), -1.0, 1.0))/(nbeads-2)
        bond_angle[jj] = avg_angle


def e2e(paras, xc, yc, zc, square_e2edistance):
    nbeads, nlayers, nskip, nconf, npoly = paras
    for jj in range(0, npoly):
        for k in range(jj*nbeads+1, jj*nbeads+nbeads+1):
            if k % nbeads == 0:  # last bead of the chain
                # print k, xc[k], yc[k], zc[k]
                x_1st = xc[k]
                y_1st = yc[k]
                z_1st = zc[k]
            elif k == 1 or (k-1) % nbeads == 0:  # first bead of the chain
                # print k, xc[k], yc[k], zc[k]
                x_last = xc[k]
                y_last = yc[k]
                z_last = zc[k]
        square_e2edistance[jj] = (
            ((x_1st-x_last)**2)+((y_1st-y_last)**2)+((z_1st-z_last)**2))


def radgy(paras, xc, yc, zc, xcm, ycm, zcm, rg):
    nbeads, nlayers, nskip, nconf, npoly = paras
    for jj in range(0, npoly):
        for kk in range(1, nbeads+1):
            k = nbeads*jj + kk
            zcm[jj] += zc[k]/nbeads
            ycm[jj] += yc[k]/nbeads
            xcm[jj] += xc[k]/nbeads
    for jj in range(0, npoly):
        for kk in range(1, nbeads+1):
            k = nbeads*jj + kk
            rg[jj] += (((xc[k]-xcm[jj])**2)+((yc[k]-ycm[jj])**2) +
                       ((zc[k]-zcm[jj])**2))/nbeads


def main(infile, outfile, paras):
    nbeads, nlayers, nskip, nconf, npoly = paras
    IN = fileinput.input(infile)
    natoms = 0
    with open(outfile, 'w') as f:
        f.write("step, lz, ndens, sqRee, sqRg\n")
        for xx in range(nskip):
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
            IN.readline()
            IN.readline()
            IN.readline()
            IN.readline()
            IN.readline()
            for j in range(1, dim):
                IN.readline()

        for xx in range(nconf):
            IN.readline()
            line = IN.readline()      # time step
            fields = string.split(line)
            step = int(fields[0])
            IN.readline()
            line = IN.readline()      # number of atoms
            fields = string.split(line)

            if xx == 0:
                # reads number of atoms from first configuration; don't support changing num atoms
                natoms = int(fields[0])
                dim = natoms+1
                xc = np.zeros(dim, np.float32)
                yc = np.zeros(dim, np.float32)
                zc = np.zeros(dim, np.float32)
                typea = [0]*dim
                mol = [0]*dim

            xcm = np.zeros(npoly+1)
            ycm = np.zeros(npoly+1)
            zcm = np.zeros(npoly+1)
            square_e2edistance = [0]*npoly
            rg = [0]*npoly
            bond_angle = [0]*npoly

            IN.readline()
            line = IN.readline()
            [xm, xp] = map(np.float, line.split())
            line = IN.readline()
            [ym, yp] = map(np.float, line.split())
            line = IN.readline()
            [zm, zp] = map(np.float, line.split())
            line = IN.readline()
            index = string.split(line)[2:]
            xbox = xp - xm
            ybox = yp - ym
            zbox = zp - zm
            vol = xbox*ybox*zbox
            ndens = natoms/vol
            lz = zbox/nlayers
            for j in range(1, dim):
                line = IN.readline()
                if "q" in set(index):
                    [ii, molj, typej, qj, x1, x2, x3,
                        n1, n2, n3] = string.split(line)
                else:
                    [ii, molj, typej, x1, x2, x3, n1,
                        n2, n3] = string.split(line)
                k = int(ii)
                typea[k] = int(typej)
                mol[k] = int(molj)

                xc[k] = xbox*np.float(x1) + int(n1)*xbox
                yc[k] = ybox*np.float(x2) + int(n2)*ybox
                zc[k] = zbox*np.float(x3) + int(n3)*zbox

            sqe2edist = 0
            sqrg = 0
            avgbangle = 0

            bangle(paras, xc, yc, zc, bond_angle)
            e2e(paras, xc, yc, zc, square_e2edistance)
            radgy(paras, xc, yc, zc, xcm, ycm, zcm, rg)
            for jj in range(0, npoly):
                avgbangle += bond_angle[jj]/npoly
                sqe2edist += square_e2edistance[jj]/npoly
                sqrg += rg[jj]/npoly
            f.write("{}, {}, {}, {}, {}, {}\n".format(
                step, ndens, lz, sqe2edist, sqrg, avgbangle))
    fileinput.close()


if __name__ == "__main__":
    infile = 'equil.lammpstrj'
    outfile = 'equilcheck.csv'
    main(infile, outfile, paras)
