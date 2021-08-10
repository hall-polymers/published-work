#Read Z1 coordinates from Z1 output files and calculate interparticle entanglements

#Author: Jeff Ethier 8/25/16

#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


from numpy import *
from math import *
import fileinput
import sys,string
from string import join

#Parameters
delr = 0.05
maxbin = 3000
npoly = 376
nbeads = 188
NP = 2
nmon = npoly*nbeads + NP

infiles = ['coords_sp_full_plus_beadinfo_and_entanglements'] #This file is a result from the Z1 algorithm
file2 = 'numtwokinks_check.txt'
OUT2 = open(file2, 'a+')
IN = fileinput.input(infiles)
print "Reading coordinate file..."
line = IN.readline() #number of polymers
totpoly = int(line)
print totpoly
line = IN.readline() #box dimensions
fields = string.split(line)
Lx = float(fields[0])
Ly = float(fields[1])
Lz = float(fields[2])
Lx2 = Lx/2
Ly2 = Ly/2
Lz2 = Lz/2
print Lx,Ly,Lz
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
typeb = 0
k = 0
ent = 0
for k in range(0,totpoly):
    beads = IN.readline()
    numbeads = int(beads)
    if numbeads == 2:
        typeb += 1
        typenp = 3
    elif k < (npoly/2) and numbeads != 2:
        typeb += 1
        typenp = 1
    else:
        typeb += 1
        typenp = 2
    for j in range(0,numbeads):
        count += 1
        m = count
        typep.append(typeb)
        typemon.append(typenp)
        line = IN.readline()
        fields = string.split(line)
        en.append(int(fields[4]))
        x = float(fields[0])
        y = float(fields[1])
        z = float(fields[2])
        xc.append(float(x+Lx2))
        yc.append(float(y+Ly2))
        zc.append(float(z+Lz2))

entx = []
enty = []
entz = []
typepart = []
numkinks = 0
count2 = 0
for k in range(0,count):
    if en[k] == 1 and typemon[k] < 3:
        numkinks = numkinks + 1
        entx.append(xc[k])
        enty.append(yc[k])
        entz.append(zc[k])
        typepart.append(typemon[k])

print numkinks
hist = [0]*maxbin
numtwokinks = 0
numselftwokinks = 0
entangle = 0
for k in range(0,numkinks):
    for l in range(k+1,numkinks):
        dist1=[]
        dx1 = min((entx[k]-entx[l])**2, (entx[k]-entx[l]-Lx)**2, (entx[k]-entx[l]+Lx)**2)
        dy1 = min((enty[k]-enty[l])**2, (enty[k]-enty[l]-Ly)**2, (enty[k]-enty[l]+Ly)**2)
        dz1 = min((entz[k]-entz[l])**2, (entz[k]-entz[l]-Lz)**2, (entz[k]-entz[l]+Lz)**2)
        dist1.append(dx1+dy1+dz1)
        rij = sqrt(min(dist1))
        BIN = int(rij/delr) + 1
        if BIN < maxbin:
            hist[BIN] = hist[BIN] + 1
        if BIN == 1 and typepart[k] != typepart[l]:
            numtwokinks = numtwokinks + 1
        elif BIN == 1 and typepart[k] == typepart[l]:
            numselftwokinks = numselftwokinks + 1


stresult = ' '.join(str(hist[x]) for x in range(0,len(hist)))
#print stresult
print "Number of two kinks:%s %s %s" %(numtwokinks,numselftwokinks,hist[1])
OUT2.write("%s %s \n" %(numtwokinks,numselftwokinks))
#OUT.write(stresult)
#OUT.close
OUT2.close

