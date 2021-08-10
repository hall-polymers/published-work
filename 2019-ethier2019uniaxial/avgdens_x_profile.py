#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#Compute mean squared end to end distance
#Run: python meansend.py < equilnvt.lammpstrj

import sys,string
from numpy import *
from math import *

nconf = 1 #number of configurations to average over
iskip = 329
npoly_part = 2256
nbeads = 160

#Initial conditions
conf = 0
natoms = npoly_part*nbeads
delta = 0.5
maxbin = 4000
avgdens = [0]*(maxbin+1)
file = 'avgdens_profile_x_frame_%d.txt' %(nconf+iskip) 
def avgdens_brush():
    global avgdens
    count = 0
    for k in range(1,dim):
        if typea[k] == 1 or typea[k] == 2:
            BIN = int(xc[k]/delta) + 1
            if BIN <= maxbin:
                avgdens[BIN] = avgdens[BIN] + 1
    return avgdens         

for i in range(0,iskip+1):
    print "skip",i
    sys.stdin.readline()
    sys.stdin.readline() #timestep
    sys.stdin.readline()
    line = sys.stdin.readline() # number of atoms
    fields = string.split(line) #reads the next line
    natoms = int(fields[0])  #reads number of atoms from first configuration; don't support changing num atoms
    dim=natoms+1
    sys.stdin.readline()
    sys.stdin.readline() #xbox bounds  
    sys.stdin.readline() #ybox bounds
    sys.stdin.readline() #zbox bounds
    sys.stdin.readline()
    for j in range(1,dim):
        sys.stdin.readline()

    conf = conf + 1
    #print conf

OUT = open(file, 'w')

print "reading config file..."
istart = iskip+1

OUT = open(file, 'w')
OUT.write("configuration #    <Density(z)>\n")

xbox1 = 0
ybox1 = 0
zbox1 = 0
for i in range(istart,nconf+iskip+1):
    print i
    sys.stdin.readline()
    sys.stdin.readline() #timestep
    sys.stdin.readline()
    line = sys.stdin.readline() # number of atoms
    fields = string.split(line) #reads the next line
    natoms = int(fields[0])  #reads number of atoms from first configuration; don't support changing num atoms
    dim=natoms+1
    xc=zeros(dim,float32)
    yc=zeros(dim,float32)
    zc=zeros(dim,float32)
    cx=zeros(dim)
    cy=zeros(dim)
    cz=zeros(dim)
    typea=[0]*dim
    mol=[0]*dim
    sys.stdin.readline()
    line = sys.stdin.readline() #xbox bounds  
    [xm,xp] = map(float,line.split())
    line = sys.stdin.readline() #ybox bounds
    [ym,yp] = map(float,line.split())
    line = sys.stdin.readline() #zbox bounds
    [zm,zp] = map(float,line.split())
    line = sys.stdin.readline()
    xbox = (xp-xm)
    ybox = (yp-ym)
    zbox = (zp-zm)
    xbox1 += (xp - xm)/nconf
    ybox1 += (yp - ym)/nconf
    zbox1 += (zp - zm)/nconf
    vol = xbox*ybox*zbox
    xbox2 = xbox/2.0
    ybox2 = ybox/2.0
    zbox2 = zbox/2.0
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3] = string.split(line) #atom number, molecule number, bead type, x,y,z, image flags
        k=int(ii) #sets the index k equal to the atom number
        typea[k]=int(typej)
        mol[k] = int(molj) #sets the chain number for each bead equal to molecule number
        if mol[k] >= 380:
            mol[k] -= 1
        xc[k] = (float(x1))
        yc[k] = (float(x2)) 
        zc[k] = (float(x3))
        #should have assigned each bead a chain number

    avgdens_brush()
#    avgzNP()
    conf = conf + 1
    #print conf
for jj in range(1,len(avgdens)):
    avgdens[jj] = float(avgdens[jj])/nconf/(ybox1*37)
    OUT.write("%8.4f %8.4f\n" % (((jj*delta)/xbox), avgdens[jj]))
    #print totavgsquaredrg
#print avgz/nconf    
OUT.close()
