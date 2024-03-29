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


import math
import string
import sys
from random import *

from numpy import *

iseed = 9328
nbeads = 100
percentTaper = 50 / 100.0
fractionA = 0.33
npurebeads1 = round(nbeads * (fractionA - percentTaper * 0.5))
npurebeads2 = round(nbeads * (1 - percentTaper) - npurebeads1)
npurebeadsA = int(npurebeads1)
npurebeadsB = int(npurebeads2)
ntapers = int(nbeads * percentTaper)
nmonomersperpoly = 1
npoly = 800
hnpoly = int(npoly / 2)
dens = 0.85
nlayers = 2
bond = 0.97
seed(iseed)

neutralizedfraction = 1
chargeonunneutralized = 0
atype = ("00", "CB", "NB")

# define monomer sequence
pureA = [1] * npurebeadsA
pureB = [2] * npurebeadsB

j = 0
taperAB = [0] * (ntapers)
for ij in range(ntapers):
    j = j + 1
    taperAB[j - 1] = [1] * int(round(npoly * (j - 0.5) / ntapers)) + [2] * (
        npoly - int(round(npoly * (j - 0.5) / ntapers))
    )
    shuffle(taperAB[j - 1])
monomerTypeTaper = [[row[i] for row in taperAB] for i in range(npoly)]

# files
INPUT_LAMMPS = open("input.lammps", "w")

nmonomers = nmonomersperpoly * npoly
ntypes = 2
ntot = nmonomers * nbeads
vol = (ntot) / dens
hz = bond * (nbeads * nlayers - 1)
side = (vol / hz) ** (1.0 / 2.0)
hx = side
hy = side
dim = ntot + 1

# simulation cell parameters
hx2 = hx / 2.0
hy2 = hy / 2.0
hz2 = hz / 2.0
vol = hx * hy * hz
nbonds = ntot - npoly

print "Total number of beads:", ntot
print "Number of chains =", npoly
print "beads in monomer =", nbeads
print "monomers total =", nmonomers
print "Number of atoms types = ", ntypes
print "seed = ", iseed
print " "
print "dens = ", dens
print "vol = ", vol
print " "
print "metric: %10.4f %10.4f %10.4f\n\n" % (hx, hy, hz)

xc = zeros(dim, float32)
yc = zeros(dim, float32)
zc = zeros(dim, float32)
cx = zeros(dim)
cy = zeros(dim)
cz = zeros(dim)

rg2ave = 0.0
rgave = 0.0
rend2ave = 0.0
typeb = [0] * dim
seq = [0] * nbeads
molnum = [0] * dim
n = 0
k = 0

for ix in xrange(hnpoly):
    lengthcurrentpoly = 0
    n = n + 1
    sequence = pureA + monomerTypeTaper[n - 1] + pureB
    for iy in range(nmonomersperpoly):
        currentmonomer = ix * nmonomersperpoly + iy
        seqnum = 0
        seq = sequence
        for iz in seq:
            seqnum = seqnum + 1
            k = k + 1
            lengthcurrentpoly = lengthcurrentpoly + 1
            typeb[k] = iz
            molnum[k] = ix + 1
            if iy == 0 and seqnum == 1:
                k1 = k
                xc[k] = random() * hx
                yc[k] = random() * hy
                zc[k] = 0
            else:
                xc[k] = xc[k - 1]
                yc[k] = yc[k - 1]
                zc[k] = zc[k - 1] + bond
for ix in xrange(hnpoly, npoly):
    lengthcurrentpoly = 0
    n = n + 1
    sequence = pureA + monomerTypeTaper[n - 1] + pureB
    for iy in range(nmonomersperpoly):
        currentmonomer = ix * nmonomersperpoly + iy
        seqnum = 0
        seq = sequence
        for iz in seq:
            seqnum = seqnum + 1
            k = k + 1
            lengthcurrentpoly = lengthcurrentpoly + 1
            typeb[k] = iz
            molnum[k] = ix + 1
            if iy == 0 and seqnum == 1:
                k1 = k
                xc[k] = random() * hx
                yc[k] = random() * hy
                zc[k] = hz
            else:
                xc[k] = xc[k - 1]
                yc[k] = yc[k - 1]
                zc[k] = zc[k - 1] - bond


for k in xrange(1, ntot + 1):
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


# input.lammps header
INPUT_LAMMPS.write("#50% Inverse Tapers - L* Structure\n")
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atoms\n" % ntot)
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
for i in xrange(1, dim):
    imol = molnum[i]
    INPUT_LAMMPS.write(
        "%6i %6i %2i %9.4f %9.4f %9.4f %6i %6i %6i\n"
        % (i, imol, typeb[i], xc[i], yc[i], zc[i], cx[i], cy[i], cz[i])
    )

# bonds
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Bonds\n")
INPUT_LAMMPS.write("\n")

for i in xrange(1, ntot):
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
