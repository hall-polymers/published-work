#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#!/usr/bin/python
# Script:  mkinput.py
# Purpose: Make input file
# Syntax:  mkinput.py
# Example: mkinput.py
# Author:  Mark Stevens
# Modified: Lisa Hall 1/10, Kevin Shen 1/2016
# Change to box center of 0 (-L/2 to L/2)
# Chains are random walk.

import sys, string, math
from numpy import *
from random import *
#-------------------------------------------------------------------------

##  40 beads, 1 monomers, 800 poly = 32000 beads, 40 bead/monomer


# INPUT PARAMETERs
nbeads = 20                 # number of beads in monomer (3, 5, 7, or 9 for ionomer project)
npoly =  800                # number of polymers      
dens = 1                # bead density
bond = 0.97                 # bond length. depends on bond potential, but close to 1 is good enough
ion_conc = 0.052
nanion = int(npoly*nbeads*ion_conc)
ncation = int(npoly*nbeads*ion_conc)

# dictionary between type and charge (if neutralizedfraction<1, some - charges will be overwritten as 0)
charge  = {1:0, 2:0, 3:+1, 4:-1}

# define monomer
sequence = [1]*int(nbeads)

# files
INPUT_LAMMPS = open('input_homoions.lammps', 'w')

# constants
ntypes = 4 # PEO + anions + cations (no type 2)
ntot = npoly*nbeads
vol = (ntot + nanion + ncation)/dens
dim = ntot + nanion + ncation + 1

# simulation cell parameters
hx = vol**(1./3.)
hy = vol**(1./3.)
hz = vol**(1./3.)
hx2 = hx/2.
hy2 = hy/2.
hz2 = hz/2.
vol = hx*hy*hz
nbonds = ntot-npoly

print 
print "Total number of beads: ", ntot + nanion + ncation
print "Number of chains = ", npoly
print "beads in monomer = ", nbeads
print "Number of ions (anions and cations) = ", ncation + nanion
print "Number of atoms types = ",ntypes
print " "
print "Geometry:"
print "dens = ", dens
print "vol = ", vol
print " "
print "metric: %10.4f %10.4f %10.4f\n\n" % (hx, hy, hz)

# init position variables
xc=zeros(dim,float32)
yc=zeros(dim,float32)
zc=zeros(dim,float32)
cx=zeros(dim)
cy=zeros(dim)
cz=zeros(dim)

# Build polymers
typeb=[0]*dim
molnum=[0]*dim
k=0

for ix in xrange(npoly):
    seqnum = 0
    seq = sequence
    for iz in seq:
            seqnum += 1
            k = k + 1
            typeb[k] = iz
            molnum[k] = ix + 1
            if seqnum == 1:
                xc[k] = random()*hx
                yc[k] = random()*hy
                zc[k] = random()*hz
            else:
            # pick random direction; scale to be bond length 
                theta = random()*2*pi
                dz = random()*2 - 1
                dx = sqrt(1-dz**2)*cos(theta)
                dy = sqrt(1-dz**2)*sin(theta)
                r = sqrt(dx*dx+dy*dy+dz*dz)
                scale = bond/r
                dx = scale*dx
                dy = scale*dy
                dz = scale*dz
                xc[k] = xc[k-1] + dx
                yc[k] = yc[k-1] + dy
                zc[k] = zc[k-1] + dz

# BUILD PENETRANTS
for ncat in range(ncation):
    mol = ntot + ncat + 1
    typeb[mol] = 3
    xc[mol] = random()*hx
    yc[mol] = random()*hy
    zc[mol] = random()*hz

for nan in range(nanion):
    mol = ntot + ncation + nan + 1
    typeb[mol] = 4
    xc[mol] = random()*hx
    yc[mol] = random()*hy
    zc[mol] = random()*hz

# PBC #added -hx2 to everything to center box LMH; correct previous version error where it's different (not just the negative) depending on if/elif loop below for y and z

for k in xrange(1,dim):
    if (xc[k] > hx):
        cx[k] = int(xc[k]/hx)
        xc[k] = xc[k] - cx[k]*hx - hx2
    elif (xc[k] < 0.0):
        cx[k] = -int((-xc[k]+hx)/hx)
        xc[k] = xc[k] - cx[k]*hx - hx2
    else:
        cx[k] = 0
        xc[k] = xc[k] - hx2
    if (yc[k] > hy):
        cy[k] = int(yc[k]/hy)
        yc[k] = yc[k] - cy[k]*hy - hy2
    elif (yc[k] < 0.0):
        cy[k] = -int((-yc[k]+hy)/hy)
        yc[k] = yc[k] - cy[k]*hy - hy2
    else:
        cy[k] = 0
        yc[k] = yc[k] - hy2
    if (zc[k] > hz):
        cz[k] = int(zc[k]/hz)
        zc[k] = zc[k] - cz[k]*hz - hz2
    elif (zc[k] < 0.0):
        cz[k] = -int((-zc[k]+hz)/hz)
        zc[k] = zc[k] - cz[k]*hz - hz2
    else:
        cz[k] = 0
        zc[k] = zc[k] - hz2

# OUTPUT headers ---------------------------------------------------------------

# input.lammps header 
INPUT_LAMMPS.write("#Homopolymers with ions KS 6/2018\n")
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atoms\n" %     int(ntot+ncation+nanion))
INPUT_LAMMPS.write("%10i    bonds\n" %     nbonds)
INPUT_LAMMPS.write("%10i    angles\n" %     0)
INPUT_LAMMPS.write("%10i    dihedrals\n" % 0)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atom types\n" % ntypes)
INPUT_LAMMPS.write("%10i    bond types\n" % 1)
INPUT_LAMMPS.write("%10i    angle types\n" % 0)
INPUT_LAMMPS.write("%10i    dihedral types\n" % 0)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write(" %16.8f %16.8f   xlo xhi\n" % (-hx2,hx2))
INPUT_LAMMPS.write(" %16.8f %16.8f   ylo yhi\n" % (-hy2,hy2))
INPUT_LAMMPS.write(" %16.8f %16.8f   zlo zhi\n" % (-hz2,hz2))
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Atoms\n")
INPUT_LAMMPS.write("\n")

# END OUTPUT headers -----------------------------------------------------------
# Atoms output
mass = 1.0

# Polymers
i = 0
imol = 0

# positions 
for i in xrange(1,dim):
    itype = typeb[i]
    # could use a dictionary here between type and segname
    if itype == 1:
        imol = molnum[i]
    else:
        #imol = npoly+1 #the molecule number for all counterions is the same; it's more like a group number
        imol = i-ntot+npoly 
    INPUT_LAMMPS.write("%6i %6i %2i %6.2f %9.4f %9.4f %9.4f %6i %6i %6i\n" % (i, imol, typeb[i], charge[typeb[i]], xc[i], yc[i], zc[i], cx[i], cy[i], cz[i]))

# Bonds
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Bonds\n")
INPUT_LAMMPS.write("\n")
ibond=0
for i in xrange(1,ntot):
    #if not at the end of the polymer
    if molnum[i+1] == molnum[i]:
        ibond = ibond+1 #the bond number
        j=i+1
        INPUT_LAMMPS.write("%8i  1 %8i %8i\n" % (ibond,i,j))

# Masses
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Masses\n")
INPUT_LAMMPS.write("\n")
for ii in xrange(1,ntypes+1):
    INPUT_LAMMPS.write("%3i  1.0\n" % ii)

#Close files
INPUT_LAMMPS.close()

print "LAMMPS output complete."



