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
# Script:  mkinput_lamellae_ions.py
# Purpose: Make input file of block copolymers with ions
# Syntax:  python mkinput_lamellae_ions.py
# Author:  Mark Stevens
# Modified: Youngmi Seo 09/2016, Kevin Shen 12/2018

#  LAMMPS types (all neutral)
#  1  soft polymer beads (A)
#  2  hard polymer beads (B)
#  3  (+) charged beads (C)
#  4  (-) charged beads (D)


import itertools
import math
import string
import sys
from random import *

from numpy import *

# -------------------------------------------------------------------------

## INPUT PARAMETERs
nbeads = 40  # number of beads in monomer
fA = 0.5
hnbeads = nbeads / 2
npoly = 1280  # number of polymers
hnpoly = npoly / 2  # half number of polymers
nmon = npoly * nbeads
nlayers = 4  # number of lamellar layers
ion_conc = 0.208  # [+]/[A]
nions = int(nbeads * npoly * fA * ion_conc) * 2  # both + and - (assuming the system is neutral)
while (nions / 2.0) % 4 != 0:
    nions = nions + 1
ncations = nanions = nions / 2
dim = nmon + nions
dens = 0.87  # bead number density
cdens = 0.18  # coverage density on a planar surface
bond = 0.97  # bond length. depends on bond potential, but close to 1 is good enough

# define monomer
sequence = [1] * int(nbeads * fA) + [2] * int(nbeads * (1 - fA))

# dictionary between type and charge
charge = {1: 0, 2: 0, 3: +1, 4: -1}
# END INPUT Parameters ------------------------------------------------------


# files
INPUT_LAMMPS = open("input.lammps", "w")

# constants for grafting polymers
gnpoly = npoly / nlayers  # Number of grafting polymers in one layer
gntot = gnpoly * hnbeads  # Number of monomers in one layer
surf = gnpoly / 2 / cdens  # xy surface area
gvol = (dim / dens - surf * bond * (nlayers * 2 - 1)) / nlayers / 2
ghx = surf ** 0.5
ghy = surf ** 0.5
ghz = gvol / surf

# simulation cell parameters
ntypes = 4
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
print "[+]/[A] =", ion_conc
print "Number of +/- ions =", nions
print "Number of atoms types = ", ntypes
print " "
print "Geometry:"
print "dens = ", dens
print "vol = ", vol
print " "
print "metric: %10.4f %10.4f %10.4f\n\n" % (hx, hy, hz)

# Init position variables
xc = zeros(dim + 1, float32)
yc = zeros(dim + 1, float32)
zc = zeros(dim + 1, float32)
cx = zeros(dim + 1)
cy = zeros(dim + 1)
cz = zeros(dim + 1)
typeb = [0] * (dim + 1)
molnum = [0] * (dim + 1)
k = 0


# Build polymers
for nl in range(0, nlayers):
    for ix in range((hnpoly / nlayers) * nl * 2, (hnpoly / nlayers) * (nl * 2 + 1)):
        seqnum = 0
        seq = sequence
        for iz in seq:
            seqnum += 1
            k += 1
            kk = k + hnbeads - 1 - ((k - 1) % hnbeads) * 2
            if seqnum == 1:
                typeb[k] = iz
                molnum[k] = ix + 1
                k1 = k
                xc[kk] = random() * ghx
                yc[kk] = random() * ghy
                zc[kk] = 0.0000 + (ghz * nl * 2 + nl * bond * 2)
            elif seqnum <= hnbeads:  # pick random direction; scale to be bond length
                typeb[k] = iz
                molnum[k] = ix + 1
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
                zc[kk] = zc[kk + 1] + dz  # don't need to shift, already done for the first bead
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
        seqnum = 0
        seq = sequence
        for iz in seq:
            seqnum += 1
            k += 1
            kk = k + hnbeads - 1 - ((k - 1) % hnbeads) * 2
            if seqnum == 1:
                typeb[k] = iz
                molnum[k] = ix + 1
                k1 = k
                xc[kk] = random() * ghx
                yc[kk] = random() * ghy
                zc[kk] = ghz + (ghz * nl * 2 + nl * bond * 2)
            elif seqnum <= hnbeads:  # pick random direction; scale to be bond length
                typeb[k] = iz
                molnum[k] = ix + 1
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

# Build ions
for nl in range(0, nlayers):
    for ii in range(1 + (ncations / nlayers) * nl, 1 + (ncations / nlayers) * (nl + 1)):
        k = ii + nmon
        xc[k] = random() * ghx
        yc[k] = random() * ghy
        zc[k] = random() * ghz + (ghz * nl * 2 + nl * bond * 2)
        typeb[k] = 3

for nl in range(0, nlayers):
    for jj in range(1 + (nanions / nlayers) * nl, 1 + (nanions / nlayers) * (nl + 1)):
        k = ii + nmon + jj
        xc[k] = random() * ghx
        yc[k] = random() * ghy
        zc[k] = random() * ghz + (ghz * nl * 2 + nl * bond * 2)
        typeb[k] = 4

# PBC #added -hx2 to everything to center box
for k in xrange(1, dim + 1):
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

print "Polymers and ions built."

# OUTPUT headers ---------------------------------------------------------------
INPUT_LAMMPS.write("# Kremer-Grest type Lamellae with ions YS 09/2016\n")
INPUT_LAMMPS.write("# 4 Layers of lamellae\n")
INPUT_LAMMPS.write("# Number of polymers: %1i\n" % npoly)
INPUT_LAMMPS.write("# Number of beads per polymer: %1i\n" % nbeads)
INPUT_LAMMPS.write("# ion fraction, [+]/[A]: %1.3f\n" % ion_conc)
INPUT_LAMMPS.write("# Number of ions: %1i\n" % nions)
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

# Positions
for i in xrange(1, dim + 1):
    itype = typeb[i]
    if itype != 3 and itype != 4:
        imol = molnum[i]
    else:
        imol = i - nmon + npoly
    INPUT_LAMMPS.write(
        "%6i %6i %2i %6.2f %9.4f %9.4f %9.4f %6i %6i %6i\n"
        % (i, imol, typeb[i], charge[typeb[i]], xc[i], yc[i], zc[i], cx[i], cy[i], cz[i])
    )

# Bonds
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Bonds\n")
INPUT_LAMMPS.write("\n")
ibond = 0
for i in xrange(1, nmon):
    # if not at the end of the polymer
    if molnum[i + 1] == molnum[i]:
        ibond = ibond + 1  # the bond number
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
