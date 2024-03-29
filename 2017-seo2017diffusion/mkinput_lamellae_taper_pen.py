#!/usr/bin/python

# This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

import itertools
import math
import string
import sys
from random import *

from numpy import *

nbeads = 40
hnbeads = nbeads / 2
nmonomersperpoly = 1
npoly = 1920
nmon = npoly * nbeads
hnpoly = npoly / 2
nlayers = 4

percentTaper = 30.0 / 100.0  # set taper length
npurebeads = int(nbeads * (1 - percentTaper))
ntapers = nbeads - npurebeads

saltratio = 0.06
npen = int(npoly * hnbeads * saltratio)
while npen % nlayers != 0:
    npen = npen + 1
dim = nmon + npen
minsep = 1.0
dens = 0.88
cdens = 0.1
bond = 0.97

atype = ("00", "A", "B", "C")

# define monomer sequence
pureA = [1] * (npurebeads / 2)
pureB = [2] * (npurebeads / 2)
taperAB = [0] * ntapers
j = 0
for ij in range(ntapers):
    j = j + 1
    taperAB[j - 1] = [2] * int(round(npoly * (j - 0.5) / ntapers)) + [1] * (
        npoly - int(round(npoly * (j - 0.5) / ntapers))
    )
    # taperAB[j-1] = [1]*int(round(npoly*(j-0.5)/ntapers)) + [2]*(npoly-int(round(npoly*(j-0.5)/ntapers))) #for inverse taper
    # randAB[j-1] = [1]*int(round(npoly*0.5)) + [2]*(npoly-int(round(npoly*0.5))) #for random midblock
    shuffle(taperAB[j - 1])
monomerTypeTaper = [[row[i] for row in taperAB] for i in range(npoly)]

# files
INPUT_LAMMPS = open("input.lammps", "w")

# parameters for grafting polymers
gnpoly = npoly / nlayers
gntot = gnpoly * hnbeads
surf = gnpoly / 2 / cdens
gvol = (dim / dens - surf * bond * (nlayers * 2 - 1)) / nlayers / 2
ghx = surf ** 0.5
ghy = surf ** 0.5
ghz = gvol / surf
ghx2 = ghx / 2.0
ghy2 = ghy / 2.0
ghz2 = ghz / 2.0
gvol = ghx * ghy * ghz
gnbonds = gntot - gnpoly

# simulation cell parameters
ntypes = 3
vol = dim / dens
surf = npoly / nlayers / cdens
hx = ghx
hy = ghy
hz = ghz * nlayers * 2 + (nlayers * 2 - 1) * bond
hx2 = hx / 2.0
hy2 = hy / 2.0
hz2 = hz / 2.0
vol = hx * hy * hz
nbonds = nmon - npoly

print
print "Total number of beads:", dim
print "Number of chains =", npoly
print "Number of Beads in monomer =", nbeads
print "Taper percent =", percentTaper * 100, "%"
print "Number of penetrants =", npen
print "Number of atoms types = ", ntypes
print " "
print "dens = ", dens
print "vol = ", vol
print " "
print "metric: %10.4f %10.4f %10.4f\n\n" % (hx, hy, hz)

# init position variables
xc = zeros(dim + 1, float32)
yc = zeros(dim + 1, float32)
zc = zeros(dim + 1, float32)
cx = zeros(dim + 1)
cy = zeros(dim + 1)
cz = zeros(dim + 1)

# Build polymers
rg2ave = 0.0
rgave = 0.0
rend2ave = 0.0
typeb = [0] * (dim + 1)
molnum = [0] * (dim + 1)
k = 0

for nl in range(0, nlayers):
    for ix in range((hnpoly / nlayers) * nl * 2, (hnpoly / nlayers) * (nl * 2 + 1)):
        lengthcurrentpoly = 0
        sequence = pureA + monomerTypeTaper[ix] + pureB
        for iy in range(nmonomersperpoly):
            currentmonomer = ix * nmonomersperpoly + iy
            seqnum = 0
            seq = sequence
            for iz in seq:
                seqnum = seqnum + 1
                k = k + 1
                kk = k + hnbeads - 1 - ((k - 1) % hnbeads) * 2
                lengthcurrentpoly = lengthcurrentpoly + 1
                if iy == 0 and seqnum == 1:
                    typeb[k] = iz
                    molnum[kk] = ix + 1
                    k1 = k
                    xc[kk] = random() * ghx
                    yc[kk] = random() * ghy
                    zc[kk] = 0.0000 + (ghz * nl * 2 + nl * bond * 2)
                elif seqnum <= hnbeads:
                    typeb[k] = iz
                    molnum[kk] = ix + 1
                    theta = random() * 2 * pi
                    dz = random() * 2 - 1
                    dx = sqrt(1 - dz ** 2) * cos(theta)
                    dy = sqrt(1 - dz ** 2) * sin(theta)
                    r = sqrt(dx * dx + dy * dy + dz * dz)
                    scale = bond / r
                    dx = scale * dx
                    dy = scale * dy
                    dz = scale * dz
                    xc[kk] = xc[kk + 1] + dx
                    yc[kk] = yc[kk + 1] + dy
                    zc[kk] = zc[kk + 1] + dz
                    if zc[kk] >= (ghz + (ghz * nl * 2 + nl * bond * 2)) or zc[kk] <= (
                        ghz * nl * 2 + nl * bond * 2
                    ):
                        zc[kk] = zc[kk + 1] - dz
                else:
                    typeb[k] = iz
                    molnum[k] = ix + 1
                    xc[k] = xc[kk - hnbeads]
                    yc[k] = yc[kk - hnbeads]
                    zc[k] = (
                        zc[kk - hnbeads]
                        - (
                            zc[kk - hnbeads]
                            - zc[kk - hnbeads + ((k - hnbeads - 1) % hnbeads)]
                            + bond / 2
                        )
                        * 2
                    )

    for ix in range((hnpoly / nlayers) * (nl * 2 + 1), (hnpoly / nlayers) * (nl * 2 + 2)):
        sequence = pureA + monomerTypeTaper[ix] + pureB
        for iy in range(nmonomersperpoly):
            currentmonomer = ix * nmonomersperpoly + iy
            seqnum = 0
            seq = sequence
            for iz in seq:
                seqnum = seqnum + 1
                k = k + 1
                kk = k + hnbeads - 1 - ((k - 1) % hnbeads) * 2
                lengthcurrentpoly = lengthcurrentpoly + 1
                if iy == 0 and seqnum == 1:
                    typeb[k] = iz
                    molnum[kk] = ix + 1
                    k1 = k
                    xc[kk] = random() * ghx
                    yc[kk] = random() * ghy
                    zc[kk] = ghz + (ghz * nl * 2 + nl * bond * 2)
                elif seqnum <= hnbeads:
                    typeb[k] = iz
                    molnum[kk] = ix + 1
                    theta = random() * 2 * pi
                    dz = random() * 2 - 1
                    dx = sqrt(1 - dz ** 2) * cos(theta)
                    dy = sqrt(1 - dz ** 2) * sin(theta)
                    r = sqrt(dx * dx + dy * dy + dz * dz)
                    scale = bond / r
                    dx = scale * dx
                    dy = scale * dy
                    dz = scale * dz
                    xc[kk] = xc[kk + 1] + dx
                    yc[kk] = yc[kk + 1] + dy
                    zc[kk] = zc[kk + 1] + dz
                    if zc[kk] >= (ghz + (ghz * nl * 2 + nl * bond * 2)) or zc[kk] <= (
                        ghz * nl * 2 + nl * bond * 2
                    ):
                        zc[kk] = zc[kk + 1] - dz
                else:
                    typeb[k] = iz
                    molnum[k] = ix + 1
                    xc[k] = xc[kk - hnbeads]
                    yc[k] = yc[kk - hnbeads]
                    zc[k] = (
                        zc[kk - hnbeads]
                        - (
                            zc[kk - hnbeads]
                            - zc[kk - hnbeads + ((k - hnbeads - 1) % hnbeads)]
                            - bond / 2
                        )
                        * 2
                    )

for k in xrange(1, nmon + 1):
    if xc[k] > hx:
        cx[k] = int(xc[k] / hx)
        xc[k] = xc[k] - cx[k] * hx - hx2
    elif xc[k] < 0.0:
        cx[k] = -int((-xc[k] + hx) / hx)
        xc[k] = xc[k] - cx[k] * hx - hx2
    else:
        cx[k] = 0
        xc[k] = xc[k] - hx2

    if yc[k] > hy:
        cy[k] = int(yc[k] / hy)
        yc[k] = yc[k] - cy[k] * hy - hy2
    elif yc[k] < 0.0:
        cy[k] = -int((-yc[k] + hy) / hy)
        yc[k] = yc[k] - cy[k] * hy - hy2
    else:
        cy[k] = 0
        yc[k] = yc[k] - hy2

    if zc[k] > hz:
        cz[k] = int(zc[k] / hz)
        zc[k] = zc[k] - cz[k] * hz - hz2
    elif zc[k] < 0.0:
        cz[k] = -int((-zc[k] + hz) / hz)
        zc[k] = zc[k] - cz[k] * hz - hz2
    else:
        cz[k] = 0
        zc[k] = zc[k] - hz2

print "Polymers built."

for nl in range(0, nlayers):
    for ii in range(1 + (npen / nlayers) * nl, 1 + (npen / nlayers) * (nl + 1)):
        k = ii + nmon
        xc[k] = random() * ghx
        yc[k] = random() * ghy
        zc[k] = random() * ghz + (ghz * nl * 2 + nl * bond * 2)
        typeb[k] = 3
        if xc[k] > hx:
            cx[k] = int(xc[k] / hx)
            xc[k] = xc[k] - cx[k] * hx - hx2
        elif xc[k] < 0.0:
            cx[k] = -int((-xc[k] + hx) / hx)
            xc[k] = xc[k] - cx[k] * hx - hx2
        else:
            cx[k] = 0
            xc[k] = xc[k] - hx2
        if yc[k] > hy:
            cy[k] = int(yc[k] / hy)
            yc[k] = yc[k] - cy[k] * hy - hy2
        elif yc[k] < 0.0:
            cy[k] = -int((-yc[k] + hy) / hy)
            yc[k] = yc[k] - cy[k] * hy - hy2
        else:
            cy[k] = 0
            yc[k] = yc[k] - hy2
        if zc[k] > hz:
            cz[k] = int(zc[k] / hz)
            zc[k] = zc[k] - cz[k] * hz - hz2
        elif zc[k] < 0.0:
            cz[k] = -int((-zc[k] + hz) / hz)
            zc[k] = zc[k] - cz[k] * hz - hz2
        else:
            cz[k] = 0
            zc[k] = zc[k] - hz2

print "Penetrants complete."


# OUTPUT headers ---------------------------------------------------------------
INPUT_LAMMPS.write("# Kremer-Grest type Lamellae with Penetrants\n")
INPUT_LAMMPS.write("# 4 Layers of lamellae\n")
INPUT_LAMMPS.write("# Number of polymers: %1i\n" % npoly)
INPUT_LAMMPS.write("# Number of beads per polymer: %1i\n" % nbeads)
INPUT_LAMMPS.write("# Percent taper: %1.1f\n" % percentTaper)
INPUT_LAMMPS.write("# Salt ratio: %1.2f\n" % saltratio)
INPUT_LAMMPS.write("# Number of penetrants: %1i\n" % npen)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atoms\n" % dim)
INPUT_LAMMPS.write("%10i    bonds\n" % nbonds)
INPUT_LAMMPS.write("%10i    angles\n" % 0)
INPUT_LAMMPS.write("%10i    dihedrals\n" % 0)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atom types\n" % ntypes)
INPUT_LAMMPS.write("%10i    bond types\n" % 2)
INPUT_LAMMPS.write("%10i    angle types\n" % 0)
INPUT_LAMMPS.write("%10i    dihedral types\n" % 0)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write(" %16.8f %16.8f   xlo xhi\n" % (-hx2, hx2))
INPUT_LAMMPS.write(" %16.8f %16.8f   ylo yhi\n" % (-hy2, hy2))
INPUT_LAMMPS.write(" %16.8f %16.8f   zlo zhi\n" % (-hz2, hz2))
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Atoms\n")
INPUT_LAMMPS.write("\n")

# END OUTPUT headers -----------------------------------------------------------
mass = 1.0

i = 0
imol = 0
ibond = 0

# positions
for i in xrange(1, dim + 1):
    if itype != 3:
        imol = molnum[i]
    elif itype == 3:
        imol = i - nmon + npoly
    INPUT_LAMMPS.write(
        "%6i %6i %2i %9.4f %9.4f %9.4f %6i %6i %6i\n"
        % (i, imol, typeb[i], xc[i], yc[i], zc[i], cx[i], cy[i], cz[i])
    )

# Bonds
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Bonds\n")
INPUT_LAMMPS.write("\n")
for i in xrange(1, nmon):
    if molnum[i + 1] == molnum[i]:
        ibond = ibond + 1
        j = i + 1
        if typeb[i] != typeb[j]:
            INPUT_LAMMPS.write("%8i  2 %8i %8i\n" % (ibond, i, j))
        else:
            INPUT_LAMMPS.write("%8i  1 %8i %8i\n" % (ibond, i, j))

# Masses
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Masses\n")
INPUT_LAMMPS.write("\n")

for ii in xrange(1, ntypes + 1):
    INPUT_LAMMPS.write("%3i  1.0\n" % ii)

# Close files
INPUT_LAMMPS.close()

print "LAMMPS output complete."
