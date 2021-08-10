#!/usr/bin/env python

#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

'''
Script:  mkinput_gy_penetrants.py
Purpose: script to generate gyroid phase artificially
Syntax: "python mkinput_gy_penetrants.py"
Author: Kevin Shen(2016/3/17, modified: 2016/7/18)
'''

import sys, string, math, itertools
from numpy import *
from random import *

###INPUT PARAMETERS
Nbeads = 40                     # Number of beads in monomer
Npoly = 1920	                # Number of polymers
fA = input("Please enter fA: ")	# Volume fraction of the minority component
Asegs = int(Nbeads*fA)
Npenetrants = int(Npoly*Nbeads*fA/16)

###CELL PARAMETERS
dens = 0.87                    # 0.865 for chi = 0.43 (attractive LJ)
bond = 0.965
accRange = 0.005     	        ### FIXME: What should be the acceptable range of g?
g0 = (1-fA)/0.067               # From Jon's Matlab result
Ntot = Nbeads*Npoly+Npenetrants # Total number of beads 
vol = Ntot/dens
lbox = vol**(1./3.)
lbox2 = lbox/2.
Nbonds = Ntot-Npoly-Npenetrants
ntypes = 3                      # Number of bead types
nbtypes = 2                     # Number of bond types
a = lbox/(2*pi)                 # Parameter deciding the gyroid lattice spacing 


# Define monomer sequence: diblock copolymers
sequence = [1]*Asegs + [2]*(Nbeads-Asegs)

# Files
INPUT_LAMMPS = open('input_diblock.lammps', 'w')

# Initial variables
xc=zeros(Ntot+1,float32)
yc=zeros(Ntot+1,float32)
zc=zeros(Ntot+1,float32)
cx=zeros(Ntot+1)
cy=zeros(Ntot+1)
cz=zeros(Ntot+1)
typeb=[0]*(Ntot+1)
molnum=[0]*(Ntot+1)
k=0
gc =[0]*2

def randomChainWalk():
    global dx, dy, dz
    theta = random()*2*pi         # pick random direction; scale to be bond length
    dz = random()*2 - 1 
    dx = sqrt(1-dz**2)*cos(theta)
    dy = sqrt(1-dz**2)*sin(theta)                    
    r = sqrt(dx*dx+dy*dy+dz*dz)
    scale = bond/r
    dx = scale*dx
    dy = scale*dy
    dz = scale*dz

def gyroid(x, y, z):
    x = x/a
    y = y/a
    z = z/a
    cx = cos(2*x)
    cy = cos(2*y)
    cz = cos(2*z)
    g1 = 10.0*(cos(x) * sin(y) + cos(y) * sin(z) + cos(z) * sin(x)) - 0.5*(cx*cy + cy*cz + cz*cx)

    x = -x
    y = -y
    z = -z
    cx = cos(2*x)
    cy = cos(2*y)
    cz = cos(2*z)
    g2 = 10.0*(cos(x) * sin(y) + cos(y) * sin(z) + cos(z) * sin(x)) - 0.5*(cx*cy + cy*cz + cz*cx)

    maxg = max(g1, g2)
    if maxg == g1:
        g = 0
    else:
        g = 1

    return maxg, g


###BUILD POLYMERS
for mol in range(Npoly):
    print mol
    seqnum = 0
    for iz in sequence:
        seqnum = seqnum + 1
        k = k + 1
        kk = k + (Asegs-1) - (((k%Nbeads)-1) % Asegs)*2
        if seqnum == 1:
            typeb[kk] = iz
            molnum[kk] = mol + 1
            xc[kk] = random()*lbox
            yc[kk] = random()*lbox
            zc[kk] = random()*lbox
            while abs(gyroid(xc[kk], yc[kk], zc[kk])[0] - g0) > accRange or gc[gyroid(xc[kk], yc[kk], zc[kk])[1]] == Npoly/2:
                xc[kk] = random()*lbox
                yc[kk] = random()*lbox
                zc[kk] = random()*lbox
            if gyroid(xc[kk], yc[kk], zc[kk])[1] == 0:
                gc[0] += 1
            else:
                gc[1] += 1
            # print gc

        elif seqnum <= Asegs:             # forward walk  
            typeb[kk] = iz
            molnum[kk] = mol + 1
            randomChainWalk()
            xc[kk] = xc[kk+1] + dx
            yc[kk] = yc[kk+1] + dy
            zc[kk] = zc[kk+1] + dz        # don't need to shift, already done for the first bead
            while gyroid(xc[kk], yc[kk], zc[kk])[0] < g0:     ### FIXME: Is this OK? 
                randomChainWalk()
                xc[kk] = xc[kk+1] + dx
                yc[kk] = yc[kk+1] + dy
                zc[kk] = zc[kk+1] + dz
        else:                             # backward walk
            typeb[k] = iz
            molnum[k] = mol + 1
            randomChainWalk()
            xc[k] = xc[k-1] + dx
            yc[k] = yc[k-1] + dy
            zc[k] = zc[k-1] + dz        
            while gyroid(xc[k], yc[k], zc[k])[0] > g0:     ### FIXME: Is this OK?
                randomChainWalk()
                xc[k] = xc[k-1] + dx
                yc[k] = yc[k-1] + dy
                zc[k] = zc[k-1] + dz 

###BUILD PENETRANTS
for npen in range(Npenetrants):
    mol = npen + Ntot - Npenetrants + 1
    # print mol
    typeb[mol] = 3
    xc[mol] = random()*lbox
    yc[mol] = random()*lbox
    zc[mol] = random()*lbox
    while gyroid(xc[mol], yc[mol], zc[mol])[0] < g0:     ### FIXME: Is this OK? 
        xc[mol] = random()*lbox
        yc[mol] = random()*lbox
        zc[mol] = random()*lbox

###IMAGE FLAGS
for k in xrange(1,Ntot+1):
    if (xc[k] > lbox):
        cx[k] = int(xc[k]/lbox)
        xc[k] = xc[k] - cx[k]*lbox - lbox2
    elif (xc[k] < 0.0):
        cx[k] = -int((-xc[k]+lbox)/lbox)
        xc[k] = xc[k] - cx[k]*lbox - lbox2
    else:
        cx[k] = 0
        xc[k] = xc[k] - lbox2
        
    if (yc[k] > lbox):
        cy[k] = int(yc[k]/lbox)
        yc[k] = yc[k] - cy[k]*lbox - lbox2
    elif (yc[k] < 0.0):
        cy[k] = -int((-yc[k]+lbox)/lbox)
        yc[k] = yc[k] - cy[k]*lbox - lbox2
    else:
        cy[k] = 0
        yc[k] = yc[k] - lbox2
        
    if (zc[k] > lbox):
        cz[k] = int(zc[k]/lbox)
        zc[k] = zc[k] - cz[k]*lbox - lbox2
    elif (zc[k] < 0.0):
        cz[k] = -int((-zc[k]+lbox)/lbox)
        zc[k] = zc[k] - cz[k]*lbox - lbox2
    else:
        cz[k] = 0
        zc[k] = zc[k] - lbox2

print "Polymers built."
print "Npoly, fA:", Npoly, ",", fA
print "Total number of beads:", Ntot
print "Number of penetrants:", Npenetrants
print "Box size:", lbox

###OUTPUT Headers ---------------------------------------------------------------

INPUT_LAMMPS.write("#Diblock gyroid phase generated by Kevin Shen (2016/03)\n")
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atoms\n" % Ntot)
INPUT_LAMMPS.write("%10i    bonds\n" % Nbonds)
INPUT_LAMMPS.write("%10i    angles\n" % 0)
INPUT_LAMMPS.write("%10i    dihedrals\n" % 0)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atom types\n" % ntypes)
INPUT_LAMMPS.write("%10i    bond types\n" % nbtypes)
INPUT_LAMMPS.write("%10i    angle types\n" % 0)
INPUT_LAMMPS.write("%10i    dihedral types\n" % 0)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write(" %16.8f %16.8f   xlo xhi\n" % (-lbox2,lbox2))
INPUT_LAMMPS.write(" %16.8f %16.8f   ylo yhi\n" % (-lbox2,lbox2))
INPUT_LAMMPS.write(" %16.8f %16.8f   zlo zhi\n" % (-lbox2,lbox2))
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Atoms\n")
INPUT_LAMMPS.write("\n")

###END OUTPUT Headers -----------------------------------------------------------

mass = 1.0
i = 0
imol = 0

# Positions 
for i in xrange(1,Ntot+1):
    itype = typeb[i]
    if itype != 3:
        imol = molnum[i]
    elif itype == 3:       #for the case there are penetrants or ions added
        imol = i-Ntot+Npoly+Npenetrants  
    INPUT_LAMMPS.write("%6i %6i %2i %9.4f %9.4f %9.4f %6i %6i %6i\n" % (i, imol, typeb[i], xc[i], yc[i], zc[i], cx[i], cy[i], cz[i]))

# Bonds
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Bonds\n")
INPUT_LAMMPS.write("\n")
jbond1 = zeros(Nbonds+1)
jbond2 = zeros(Nbonds+1)
pbond1 = zeros(Nbonds+1)
pbond2 = zeros(Nbonds+1)
ibond=0
for i in xrange(1,Ntot-Npenetrants):       
        if molnum[i+1] == molnum[i]:           #if not at the end of the polymer
            ibond = ibond+1 #bond number
            j=i+1
            if typeb[i] != typeb[j]:
                INPUT_LAMMPS.write("%8i  2 %8i %8i\n" % (ibond,i,j))
            else: 
                INPUT_LAMMPS.write("%8i  1 %8i %8i\n" % (ibond,i,j))


# Masses
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Masses\n")
INPUT_LAMMPS.write("\n")

for ii in xrange(1,ntypes+1):
    INPUT_LAMMPS.write("%3i  1.0\n" % ii)

# Close files
INPUT_LAMMPS.close()

print "LAMMPS output complete."
