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
Script: mkinput_cyl_penetrants.py
Purpose: script to generate cylinder phase artificially
Syntax: "python mkinput_cyl_penetrants.py"
Author: Kevin Shen(2016/7/20)
'''

import sys, string, math, itertools
from numpy import *
from random import *

###INPUT PARAMETERS
Nbeads = 100                 # Number of beads in monomer
Npoly = 1920				# Number of polymers
fA = input("please enter the fA you want: ")					# Volume fraction of the minority component
NlayersX = 3								#layers in x direction
NlayersY = 2 								#layers in y direction
Ncyl = NlayersX*NlayersY*2 					#no of cylinders (FCC type)
Nunit = Npoly/NlayersX/NlayersY				#polymer units in each layer
Asegs = int(Nbeads*fA)
Npenetrants = int(Npoly*Nbeads*fA/16)		#[A]:[C] = 16:1
Nunit_pen = Npenetrants/NlayersX/NlayersY
Nunit_pen_left = Npenetrants%(NlayersX*NlayersY)
print "Nunit_pen_left =", Nunit_pen_left

###CELL PARAMETERS
dens = 0.87                 # Density of beads  ###FIXME: What's the density?
bond = 0.965
Ntot = Nbeads*Npoly+Npenetrants         # Total number of beads 
vol = Ntot/dens
zbox = 40.                  # FIXME: Need to change to coverage density someday.
xbox = ((vol*NlayersX)/((3**(0.5))*zbox*NlayersY))**(0.5)  # vol = zbox*xbox*ybox
ybox = (xbox*NlayersY*(3**(0.5)))/NlayersX
xbox2 = xbox/2.
ybox2 = ybox/2.
zbox2 = zbox/2.
Nbonds = Ntot-Npoly-Npenetrants
ntypes = 3                  # Number of monomer types
nbtypes = 2                 # Number of bond types
radius = (xbox*ybox*fA/Ncyl/pi)**(1./2.)     # (xbox * ybox) : 12*pi*r^2 = 1 : fA
rsquare = radius**2

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
nc=nc2=0

def progressBar(i, end_val, bar_length = 50):
    percent = float(i) / end_val
    hashes = '=' * int(round(percent * bar_length))
    spaces = ' ' * (bar_length - len(hashes))
    sys.stdout.write("\rPercent: [%s] %.1f%%, %i/%i" % (hashes + spaces, percent * 100, i+1, end_val))
    sys.stdout.flush()

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

def inCylinder(xcenter, ycenter, x, y):
    if (x - xcenter)**2 + (y - ycenter)**2 <= rsquare:
        return True
    return False

def inOtherCylinders(xcenter, ycenter, x, y):
    if (x - xcenter - xbox2/NlayersX)**2 + (y - ycenter - ybox2/NlayersY)**2 <= rsquare:
        return True
    elif (x - xcenter + xbox2/NlayersX)**2 + (y - ycenter - ybox2/NlayersY)**2 <= rsquare:
        return True
    elif (x - xcenter - xbox2/NlayersX)**2 + (y - ycenter + ybox2/NlayersY)**2 <= rsquare:
        return True
    elif (x - xcenter + xbox2/NlayersX)**2 + (y - ycenter + ybox2/NlayersY)**2 <= rsquare:
        return True
    elif (x - xcenter - xbox/NlayersX)**2 + (y - ycenter)**2 <= rsquare:
        return True
    elif (x - xcenter + xbox/NlayersX)**2 + (y - ycenter)**2 <= rsquare:
        return True
    return False

###BUILD POLYMERS
for layerX in range(NlayersX):
    for layerY in range(NlayersY):
        for mol in range(Nunit):
            seqnum = 0
            if mol < Nunit/2:
                xcenter = xbox2/NlayersX + layerX*xbox/NlayersX
                ycenter = ybox2/NlayersY + layerY*ybox/NlayersY
            else:
                xcenter = xbox/NlayersX + layerX*xbox/NlayersX
                ycenter = ybox/NlayersY + layerY*ybox/NlayersY
            mol += nc * Nunit
            # print mol
            progressBar(mol, Npoly)
            for iz in sequence:
                seqnum = seqnum + 1
                k = k + 1
                kk = k + (Asegs-1) - (((k%Nbeads)-1) % Asegs)*2
                if seqnum == 1:
                    typeb[kk] = iz
                    molnum[kk] = mol + 1
                    theta = random()*2*pi
                    xc[kk] = xcenter + radius*cos(theta)
                    yc[kk] = ycenter + radius*sin(theta)
                    zc[kk] = random()*zbox
                elif seqnum <= Asegs:             # forward walk
                    typeb[kk] = iz
                    molnum[kk] = mol + 1
                    randomChainWalk()
                    xc[kk] = xc[kk+1] + dx
                    yc[kk] = yc[kk+1] + dy
                    zc[kk] = zc[kk+1] + dz
                    while inCylinder(xcenter, ycenter, xc[kk], yc[kk]) == False:
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
                    while inCylinder(xcenter, ycenter, xc[k], yc[k]) == True or inOtherCylinders(xcenter, ycenter, xc[k], yc[k]) == True:
                        randomChainWalk()
                        xc[k] = xc[k-1] + dx
                        yc[k] = yc[k-1] + dy
                        zc[k] = zc[k-1] + dz
        nc += 1

###BUILD PENETRANTS
for layerX in range(NlayersX):
    for layerY in range(NlayersY):
        for mol in range(Nunit_pen):
            if mol < Nunit_pen/2:
                xcenter = xbox2/NlayersX + layerX*xbox/NlayersX
                ycenter = ybox2/NlayersY + layerY*ybox/NlayersY
            else:
                xcenter = xbox/NlayersX + layerX*xbox/NlayersX
                ycenter = ybox/NlayersY + layerY*ybox/NlayersY
            mol += nc2 * Nunit_pen + Ntot - Npenetrants + 1
            #print mol
            typeb[mol] = 3
            theta = random()*2*pi
            xc[mol] = xcenter + random()*radius*cos(theta)
            yc[mol] = ycenter + random()*radius*sin(theta)
            zc[mol] = random()*zbox
            last_mol = mol
        nc2 += 1

for pen_left in range(Nunit_pen_left):
    mol = last_mol + pen_left + 1
    typeb[mol] = 3
    xc[mol] = xc[pen_left*Nunit_pen+1]
    yc[mol] = yc[pen_left*Nunit_pen+1]
    zc[mol] = zc[pen_left*Nunit_pen+1]

###IMAGE FLAGS
for k in xrange(1,Ntot+1):
    if (xc[k] > xbox):
        cx[k] = int(xc[k]/xbox)
        xc[k] = xc[k] - cx[k]*xbox - xbox2
    elif (xc[k] < 0.0):
        cx[k] = -int((-xc[k]+xbox)/xbox)
        xc[k] = xc[k] - cx[k]*xbox - xbox2
    else:
        cx[k] = 0
        xc[k] = xc[k] - xbox2
    if (yc[k] > ybox):
        cy[k] = int(yc[k]/ybox)
        yc[k] = yc[k] - cy[k]*ybox - ybox2
    elif (yc[k] < 0.0):
        cy[k] = -int((-yc[k]+ybox)/ybox)
        yc[k] = yc[k] - cy[k]*ybox - ybox2
    else:
        cy[k] = 0
        yc[k] = yc[k] - ybox2
        
    if (zc[k] > zbox):
        cz[k] = int(zc[k]/zbox)
        zc[k] = zc[k] - cz[k]*zbox - zbox2
    elif (zc[k] < 0.0):
        cz[k] = -int((-zc[k]+zbox)/zbox)
        zc[k] = zc[k] - cz[k]*zbox - zbox2
    else:
        cz[k] = 0
        zc[k] = zc[k] - zbox2

print "\nPolymers built."
print "Npoly, fA:", Npoly, fA
print "Total number of beads:", Ntot
print "Number of penetrants:", Npenetrants
print "radius:", radius
print "xbox:", xbox
print "ybox:", ybox
print "zbox:", zbox

###OUTPUT Headers ---------------------------------------------------------------

INPUT_LAMMPS.write("#Diblock cylinder phase generated by Kevin Shen (2016/03)\n")
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
INPUT_LAMMPS.write(" %16.8f %16.8f   xlo xhi\n" % (-xbox2,xbox2))
INPUT_LAMMPS.write(" %16.8f %16.8f   ylo yhi\n" % (-ybox2,ybox2))
INPUT_LAMMPS.write(" %16.8f %16.8f   zlo zhi\n" % (-zbox2,zbox2))
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
