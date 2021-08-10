#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#Read Z1 coordinates from Z1 output files and output double kink locations as lammpstrj file

#Author: Jeff Ethier 8/25/16

from numpy import *
from math import *
import fileinput
import sys,string
from string import join

#Parameters
delr = 0.05
maxbin = 3000
npoly = 2256
nbeads = 160
NP = 12
nmon = npoly*nbeads + NP
filenum=0

infiles = ['coords_sp_full_plus_beadinfo_and_entanglements']
file = "coords_entangle_check.lammpstrj" #%(filenum)
OUT = open(file, 'a+')
IN = fileinput.input(infiles)
print "Reading coordinate file..."
line = IN.readline() #number of polymers
totpoly = int(line)
line = IN.readline() #box dimensions
fields = string.split(line)
Lx = float(fields[0])
Ly = float(fields[1])
Lz = float(fields[2])
Lx2 = Lx/2
Ly2 = Ly/2
Lz2 = Lz/2
#print Lx,Ly,Lz
count = 0
xc=[]
yc=[]
zc=[]
en=[]
cx=[]
cy=[]
cz=[]
typep = []
typemon = []
molnum = []
typeb = 0
molnum = 0
k = 0
ent = 0
nps = 13
pairs = 0
polynum = []
for k in range(0,totpoly):
    beads = IN.readline()
    numbeads = int(beads)
    if numbeads == 2 and k>=2256:
        pairs += 1
        if pairs==73:
            nps += 1
            pairs = 1
        typeb += 1
        typenp = nps
    elif k < (188):
        typeb += 1
        typenp = 1
    elif k >= (188) and k < (376):
        typeb += 1
        typenp = 2
    elif k >= (376) and k < (564):
        typeb += 1
        typenp = 3
    elif k >= (564) and k < (752):
        typeb += 1
        typenp = 4
    elif k >= (752) and k < (940):
        typeb += 1
        typenp = 5
    elif k >= (940) and k < (1128):
        typeb += 1
        typenp = 6
    elif k >= (1128) and k < (1316):
        typeb += 1
        typenp = 7
    elif k >= (1316) and k < (1504):
        typeb += 1
        typenp = 8
    elif k >= (1504) and k < (1692):
        typeb += 1
        typenp = 9
    elif k >= (1692) and k < (1880):
        typeb += 1
        typenp = 10
    elif k >= (1880) and k < (2068):
        typeb += 1
        typenp = 11
    elif k >= (2068) and k < (2256):
        typeb += 1
        typenp = 12

    for j in range(0,numbeads):
        count += 1
        m = count
        typep.append(typeb)
        typemon.append(typenp)
        polynum.append(k+1)
        line = IN.readline()
        fields = string.split(line)
        en.append(int(fields[4]))
        x = float(fields[0])
        y = float(fields[1])
        z = float(fields[2])
        xc.append(float(x))
        yc.append(float(y))
        zc.append(float(z))

entx = []
enty = []
entz = []
typepart = []
typepoly = []
numkinks = 0
count2 = 0
for k in range(0,count):
    if en[k] == 1 and typemon[k] < 13:
        numkinks = numkinks + 1
        entx.append(xc[k])
        enty.append(yc[k])
        entz.append(zc[k])
        typepart.append(typemon[k])
        typepoly.append(polynum[k])


countnps = 0
for k in range(0,count):
    if typemon[k] >= 13:
        countnps += 1
#print count
numatoms=0
totatoms = numatoms*2 + countnps
step=0
OUT.write("ITEM: TIMESTEP\n")
OUT.write("%i\n" %step)
OUT.write("ITEM: NUMBER OF ATOMS\n")
OUT.write("%s\n" %totatoms)
OUT.write("ITEM: BOX BOUNDS pp pp pp\n")
OUT.write("%s %s \n" % (-Lx2,Lx2))
OUT.write("%s %s \n" % (-Ly2,Ly2))
OUT.write("%s %s \n" % (-Lz2,Lz2))
OUT.write("ITEM: ATOMS id mol type xu yu zu\n")

hist = [0]*maxbin
numtwokinks = 0
numselftwokinks = 0
entangle = 0
atomid = 0
for k in range(0,numkinks):
    for l in range(k+1,numkinks):
        dist1=[]
        dx1 = min((entx[k]-entx[l])**2, (entx[k]-entx[l]-Lx)**2, (entx[k]-entx[l]+Lx)**2)
        dy1 = min((enty[k]-enty[l])**2, (enty[k]-enty[l]-Ly)**2, (enty[k]-enty[l]+Ly)**2)
        dz1 = min((entz[k]-entz[l])**2, (entz[k]-entz[l]-Lz)**2, (entz[k]-entz[l]+Lz)**2)
        dist1.append(dx1+dy1+dz1)
        #if sqrt(min(dist1)) < Lx and sqrt(min(dist1)) > -Lx:
        rij = sqrt(min(dist1))
        #print sqrt(min(dist1))
        #rij = sqrt(min(dist1))
        BIN = int(rij/delr) + 1
        if BIN == 1 and typepart[k] != typepart[l]:
            numtwokinks += 1
            atomid += 1
            OUT.write("%i %i %i %f %f %f\n" % (atomid, typepoly[k], 1, entx[k], enty[k], entz[k]))
            atomid += 1
            OUT.write("%i %i %i %f %f %f\n" % (atomid, typepoly[l], 1, entx[l], enty[l], entz[l]))

print numtwokinks
for k in range(0,count):
    if typemon[k] >= 13:
        atomid += 1
        OUT.write("%i %i %i %f %f %f\n" % (atomid, max(typepoly)+typemon[k]-NP, 2, xc[k], yc[k], zc[k]))

#print "Done."
OUT.close


