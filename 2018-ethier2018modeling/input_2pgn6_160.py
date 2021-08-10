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


# Script to generate polymer grafted nanoparticle on a surface
#
# input paramters are
#  npoly - number of chains
#  N - number of nanoparticles
#  L - boxsize (set to make density of polymer 0.85)
#  constraint -
#    wall: assumes that the z=0 plane is a wall

# edited by Jeff Ethier 06/30/2016


#  lammps types
#  1         rigid monomer on nanoparticle
#  2         all other monomers
#  3         nanoparticles

import itertools
import math
import string
import sys
from random import *

from numpy import *

# -------------------------------------------------------------------------
# -----101 beads, 1 monomers, 1 poly = 101 beads, 101 bead/polymer (You can change this below)------#

# INPUT PARAMETERs
iseed = 9328  # random number seed
nmonomersperpoly = 1
bondp = 1.0  # bond length. depends on bond potential, but close to 1 is good enough
angletypes = 0

# Nanoparticle parameters
N = 2  # number of nanoparticles
R = 5.5  # radius of nanoparticles
Vnp = N * ((4.0 / 3.0) * pi * ((R - 0.5) ** 3))  # total volume of nanoparticles
gdens = 0.60  # grafting density of nanoparticle
npoly = (
    int(gdens * pi * (2 * (R - 0.5)) ** 2) * N
)  # number of polymers (calculated from grafting density)
Vmon = (4.0 / 3.0) * pi * (0.5 ** 3)
corevp = 0.032
nbeads = int((1 - corevp) * Vnp / corevp / Vmon / npoly)
while nbeads % 10 > 0:
    if nbeads % 10 < 2.5 or (nbeads % 10 < 7.5 and nbeads % 10 > 5):
        nbeads -= 1
    else:
        nbeads += 1
    if nbeads % 5 == 0:
        break

# define monomer
if nbeads == nbeads:
    sequence = [1] * nbeads
else:
    sys.exit("define sequence of monomers.")

# END INPUT Parameters ------------------------------------------------------

# files
INPUT_LAMMPS = open("input.lammps", "w")
nmonomers = nmonomersperpoly * npoly
nmon = nmonomers * nbeads
# constants
ntypes = 3  # graft bead, monomer on chain, and nanoparticles
ntot = nmonomers * nbeads + N
dim = ntot + 1
nanoi = dim
nanomolnum = npoly + 1
massNP = (0.85) * (4.0 / 3.0) * pi * ((R - 0.5) ** 3)

masses = [1.0, 1.0, massNP]

# simulation cell parameters
Vmon = (4.0 / 3.0) * pi * (0.5 ** 3) * nmon
densTarget = 0.85  # bead density target
V = nmon / densTarget + Vnp
L = (V) ** (1.0 / 3.0)
nanoSize = (2 * (R - 0.5)) / L
dens = 0.85 * (
    1 - (N * (((4.0 / 3.0) * pi * (nanoSize ** 3) / (2 ** 3))))
)  # average bead density across entire box, accounting for nanoparticle
Lx = int(L) + 10
Ly = int(L) + 10
Lz = int(L)
Lx2 = Lx / 2.0
Ly2 = Ly / 2.0
Lz2 = Lz / 2.0
vol = Lx * Ly * Lz
nbonds = ntot - npoly - N
k = 0
dim = ntot + 1
molnum = [0] * dim
volmon = (4.0 / 3.0) * pi * nmon * (0.5) ** 3
corevol = Vnp / (Vnp + volmon)

for ix in xrange(npoly):
    lengthcurrentpoly = 0
    for iy in range(nmonomersperpoly):
        currentmonomer = ix * nmonomersperpoly + iy
        seq = sequence
        seqnum = 0
        for iz in seq:
            seqnum = seqnum + 1
            k = k + 1
            lengthcurrentpoly = lengthcurrentpoly + 1
            molnum[k] = ix + 1

# OUTPUT headers ---------------------------------------------------------------
# input.lammps header
INPUT_LAMMPS.write("#Polymer Grafted Nanoparticles JE 03/2016\n")
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atoms\n" % ntot)
INPUT_LAMMPS.write("%10i    bonds\n" % nbonds)
INPUT_LAMMPS.write("%10i    angles\n" % 0)
INPUT_LAMMPS.write("%10i    dihedrals\n" % 0)
# INPUT_LAMMPS.write("%10i    impropers\n" % 0)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atom types\n" % ntypes)
INPUT_LAMMPS.write("%10i    bond types\n" % 1)
INPUT_LAMMPS.write("%10i    angle types\n" % angletypes)
INPUT_LAMMPS.write("%10i    dihedral types\n" % 0)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write(" %16.8f %16.8f   xlo xhi\n" % (0, Lx * 3))
INPUT_LAMMPS.write(" %16.8f %16.8f   ylo yhi\n" % (0, Ly * 3))
INPUT_LAMMPS.write(" %16.8f %16.8f   zlo zhi\n" % (0, Lz))
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Atoms\n")
INPUT_LAMMPS.write("\n")

# END OUTPUT headers -----------------------------------------------------------


print(" ")
print("Total number of beads:", ntot)
print("Number of chains =", npoly)
print("beads in polymer =", nbeads)
print("monomers total =", nmonomers)
print("nanoparticles =", N)
print("nanoparticle core volume fraction = ", corevol)
print("seed = ", iseed)
print(" ")
print("Geometry:")
print("dens = ", dens)
print("vol = ", V)
print("Radius of nanoparticle: ", R - 0.5)
print(" ")
# print("metric: %10.4f %10.4f %10.4f\n\n" % (hx, hy, hz))

# init position variables
xc = zeros(dim, float32)
yc = zeros(dim, float32)
zc = zeros(dim, float32)
cx = zeros(dim)
cy = zeros(dim)
cz = zeros(dim)

# Nanoparticle Positions---------------------------------------------------------------------------------------------------------------
RX_NEW = [Lx + Lx2, Lx + Lx2 + 2 * R + 10]
RY_NEW = [Ly + Ly2, Ly + Ly2 + 2 * R + 10]
RZ_NEW = [Lz2, Lz2]
N = 2
for I in xrange(0, N):
    RX_NEW[I] = RX_NEW[I]
    RY_NEW[I] = RY_NEW[I]
    RZ_NEW[I] = RZ_NEW[I]
NPs = []
Mnp = 2
for I in xrange(0, N):
    NPs.append([RX_NEW[I], RY_NEW[I], RZ_NEW[I]])

# ------------------------------------------------------------------------------------------------------------------------------------------
# Polymer Initial Positions------------------------------------------------------------------------------------------------------------------------------
# build a list of NPs to check against for each octant
# note that this includes 2 extra periodic copies of each NP per dimension, e.g. x-L, x, and x+L, so that the periodic box need not be considered later
NP_octant = []
sphere_dist = 0
r = [0.0, 0.0, 0.0]
for octant in xrange(8):
    NP_octant.append([])
    xlo = -R - 1 + (octant % 2) * Lx / 2
    xhi = Lx / 2 + R + 1 + (octant % 2) * Lx / 2
    ylo = -R - 1 + ((octant / 2) % 2) * Ly / 2
    yhi = Ly / 2 + R + 1 + ((octant / 2) % 2) * Ly / 2
    zlo = -R - 1 + (octant / 4) * Lz / 2
    zhi = Lz / 2 + R + 1 + (octant / 4) * Lz / 2

    for np in xrange(Mnp):
        for x in [NPs[np][0], NPs[np][0] + Lx, NPs[np][0] - Lx]:
            for y in [NPs[np][1], NPs[np][1] + Ly, NPs[np][1] - Ly]:
                for z in [NPs[np][2], NPs[np][2] + Lz, NPs[np][2] - Lz]:
                    if xlo < x and x < xhi and ylo < y and y < yhi and zlo < z and z < zhi:
                        NP_octant[octant].append([x, y, z])

# loop over all the chains
for numNP in range(0, N):
    grafts = [[0, 0, 0]]
    for m in xrange(npoly / 2):
        # place the starting location, as long as it's not inside a NP
        ir = [0, 0, 0]
        dz = (2 * R) * random() - R
        theta = 2 * pi * random()
        ds = sqrt((R ** 2) - (dz * dz))
        dx = ds * cos(theta)
        dy = ds * sin(theta)
        if numNP > 0:
            r[0] = dx + Lx + Lx2 + R * 2 + 10
            r[1] = dy + Ly + Ly2 + R * 2 + 10
            r[2] = dz + Lz2
        else:
            r[0] = dx + Lx + Lx2
            r[1] = dy + Ly + Ly2
            r[2] = dz + Lz2
        while (
            min(
                (r[0] - grafts[s][0]) ** 2 + (r[1] - grafts[s][1]) ** 2 + (r[2] - grafts[s][2]) ** 2
                for s in range(len(grafts))
            )
            < 0.9
        ):
            dz = (2 * R) * random() - R
            theta = 2 * pi * random()
            ds = sqrt((R ** 2) - (dz * dz))
            dx = ds * cos(theta)
            dy = ds * sin(theta)
            if numNP > 0:
                r[0] = dx + Lx + Lx2 + R * 2 + 10
                r[1] = dy + Ly + Ly2 + R * 2 + 10
                r[2] = dz + Lz2
            else:
                r[0] = dx + Lx + Lx2
                r[1] = dy + Ly + Ly2
                r[2] = dz + Lz2
        if numNP > 0:
            m = m + (npoly / 2)
        k = m * nbeads + 1
        xc[k] = r[0]
        yc[k] = r[1]
        zc[k] = r[2]
        grafts.append([r[0], r[1], r[2]])

        INPUT_LAMMPS.write(
            "%6i %6i %2i %9.4f %9.4f %9.4f %6i %6i %6i\n"
            % (m * (nbeads) + 1, m + 1, 1, xc[k], yc[k], zc[k], cx[k], cy[k], cz[k])
        )

        # initialize the previous step (this makes the first step free to go in any direction, so long as the required angle is < 90)
        dx_prev = 0
        dy_prev = 0
        dz_prev = 0
        monomers = []
        monomers.append([r[0], r[1], r[2]])
        n = 0
        # loop over all the links in the chain
        conflict = False
        while len(monomers) < (nbeads):
            # choose a random point on a unit sphere
            dz = 2 * random() - 1.0
            theta = 2 * pi * random()
            ds = sqrt(1 - dz * dz)
            dx = ds * cos(theta)
            dy = ds * sin(theta)
            # while the angle/dot product is too low, retry
            # while (dx_prev or dy_prev or dz_prev) and dx*dx_prev + dy*dy_prev + dz*dz_prev < 0.9:
            while dx * dx_prev + dy * dy_prev + dz * dz_prev < -0.5:
                dz = 2 * random() - 1.0
                theta = 2 * pi * random()
                ds = sqrt(1 - dz * dz)
                dx = ds * cos(theta)
                dy = ds * sin(theta)

            # move bead by dx, dy, dz
            r[0] += dx
            r[1] += dy
            r[2] += dz

            # check for conflicts with existing NPs
            conflict = False
            for np in xrange(len(NPs)):
                dx2 = min(
                    (r[0] - NPs[np][0]) ** 2,
                    (r[0] - NPs[np][0] + (Lx)) ** 2,
                    (r[0] - NPs[np][0] - (Lx)) ** 2,
                )  # , (r[0]-NPs[np][0]+(2*Lx))**2, (r[0]-NPs[np][0]-(2*Lx))**2, (r[0]-NPs[np][0]+(3*Lx))**2, (r[0]-NPs[np][0]-(3*Lx))**2 )
                dy2 = min(
                    (r[1] - NPs[np][1]) ** 2,
                    (r[1] - NPs[np][1] + (Ly)) ** 2,
                    (r[1] - NPs[np][1] - (Ly)) ** 2,
                )  # , (r[1]-NPs[np][1]+(2*Ly))**2, (r[1]-NPs[np][1]-(2*Ly))**2, (r[1]-NPs[np][1]+(3*Ly))**2, (r[1]-NPs[np][1]-(3*Ly))**2 )
                dz2 = min(
                    (r[2] - NPs[np][2]) ** 2,
                    (r[2] - NPs[np][2] + (Lz)) ** 2,
                    (r[2] - NPs[np][2] - (Lz)) ** 2,
                )  # , (r[2]-NPs[np][2]+(2*Lz))**2, (r[2]-NPs[np][2]-(2*Lz))**2, (r[2]-NPs[np][2]+(3*Lz))**2, (r[2]-NPs[np][2]-(3*Lz))**2 )
                dist = dx2 + dy2 + dz2
                if dist < (R) ** 2:
                    conflict = True
                    break

            # check for conflicts with overlapping beads in current chain
            for mon in xrange(len(monomers)):
                dist = (
                    (r[0] - monomers[mon][0]) ** 2
                    + (r[1] - monomers[mon][1]) ** 2
                    + (r[2] - monomers[mon][2]) ** 2
                )
                if dist < 1.0:
                    conflict = True
                    break
                if r[2] < 2.0 or r[2] > Lz - 2.0:
                    conflict = True

            # if no overlaps, write data to file
            if not conflict:
                monomers.append([r[0], r[1], r[2]])
                n = n + 1
                k = m * (nbeads) + 1 + n
                dx_prev = dx
                dy_prev = dy
                dz_prev = dz
                xc[k] = r[0]
                yc[k] = r[1]
                zc[k] = r[2]
                # enforce periodic BCs
                ## if (xc[k] > Lx):
                ##     cx[k] = int(xc[k]/Lx)
                ##     xc[k] = xc[k] - cx[k]*Lx - Lx2
                ## elif (xc[k] < 0.0):
                ##     cx[k] = -int((-xc[k]+Lx)/Lx)
                ##     xc[k] = xc[k] - cx[k]*Lx - Lx2
                ## else:
                ##     cx[k] = 0
                ##     xc[k] = xc[k] - Lx2
                ## if (yc[k] > Ly):
                ##     cy[k] = int(yc[k]/Ly)
                ##     yc[k] = yc[k] - cy[k]*Ly - Ly2
                ## elif (yc[k] < 0.0):
                ##     cy[k] = -int((-yc[k]+Ly)/Ly)
                ##     yc[k] = yc[k] - cy[k]*Ly - Ly2
                ## else:
                ##     cy[k] = 0
                ##     yc[k] = yc[k] - Ly2
                ## if (zc[k] > Lz):
                ##     cz[k] = int(zc[k]/Lz)
                ##     zc[k] = zc[k] - cz[k]*Lz - Lz2
                ## elif (zc[k] < 0.0):
                ##     cz[k] = -int((-zc[k]+Lz)/Lz)
                ##     zc[k] = zc[k] - cz[k]*Lz - Lz2
                ## else:
                ##     cz[k] = 0
                ##     zc[k] = zc[k] - Lz2

                INPUT_LAMMPS.write(
                    "%6i %6i %2i %9.4f %9.4f %9.4f %6i %6i %6i\n"
                    % (m * (nbeads) + 1 + n, m + 1, 2, xc[k], yc[k], zc[k], cx[k], cy[k], cz[k])
                )

            # if overlap, reverse step and try again
            if conflict:
                r[0] -= dx
                r[1] -= dy
                r[2] -= dz

# output the NP locations
for np in xrange(Mnp):
    INPUT_LAMMPS.write(
        "%6i %6i %2i %9.4f %9.4f %9.4f %6i %6i %6i\n"
        % (
            npoly * (nbeads) + np + 1,
            npoly + np + 1,
            3,
            NPs[np][0],
            NPs[np][1],
            NPs[np][2],
            0,
            0,
            0,
        )
    )


print("Polymers built.")

# Bonds--------------------------------------------------------------------
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Bonds\n")
INPUT_LAMMPS.write("\n")
jbond1 = zeros(nbonds + 1)
jbond2 = zeros(nbonds + 1)
pbond1 = zeros(nbonds + 1)
pbond2 = zeros(nbonds + 1)
ibond = 0
pbond = 0
i0 = 0

for i in xrange(1, ntot - N):
    # if not at the end of the polymer
    if molnum[i + 1] == molnum[i]:
        ibond = ibond + 1  # the bond number
        j = i + 1
        INPUT_LAMMPS.write("%8i 1 %8i %8i\n" % (ibond, i, j))
# ---------------------------------------------------------------------------

# Masses
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Masses\n")
INPUT_LAMMPS.write("\n")

for ii in range(1, ntypes + 1):
    INPUT_LAMMPS.write("%3i  %.1f \n" % (ii, masses[ii - 1]))
    print(ii)


INPUT_LAMMPS.close()
print("LAMMPS output complete.")
