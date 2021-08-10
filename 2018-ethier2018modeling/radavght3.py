#Radially averaged height profile
#J. Ethier
#July 20, 2016
#python radavght.py < equilnvt.lammpstrj


#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


import sys,string
from numpy import *
from math import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

#fig = plt.figure()
#ax = fig.gca(projection='2d')
nconf = 5 #number of configurations to average over
iskip = 5995
npoly = 376
nbeads = 160
NP = 2    #number of nanoparticles
R = 5    #radius of NP

#Initial conditions
conf = 0
natoms = npoly*nbeads+NP

file = 'radavght_Eps4.txt'
thresh = 0.25
countconf = 0
xlength = 61
ylength = 61
xgrid = [[0 for x in range(xlength)] for y in range(ylength)]
ygrid = [[0 for x in range(xlength)] for y in range(ylength)]
zgrid = [[0 for x in range(xlength)] for y in range(ylength)]

def radavght():
    x = xc[natoms]
    y = yc[natoms]
    xstep = 2
    ystep = 2
    zstep = 0.1
    for ydir in range(0,ylength):
        x = xc[natoms]
        y = yc[natoms]
        if ydir == 0:
            y += 0
        elif ydir > 0.5*(ylength):
            y = yc[natoms]
            y -= ystep*(ydir-(ylength+1)*0.5+1)
        else:
            y = yc[natoms]
            y += ystep*(ydir)
        for xdir in range(0,xlength):
 #       print rstart,rend
            if xdir > 0.5*(xlength):
                x = xc[natoms]
                x -= xstep*(xdir-0.5*(xlength+1)+1)
            elif xdir == 0:
                x += 0
            else:
                x = xc[natoms]
                x += xstep*xdir
            z = zp + 2.0
            density = 0
            conflict = False
            while density < thresh:
                z -= zstep
                count = 0
                for jj in range(0,npoly):
                    for kk in range(1,nbeads+1):
                        k = nbeads*jj + kk
                        #check for monomers in grid
                        if xc[k] < (x+1) and xc[k] > (x-1) and yc[k] < (y+1) and yc[k] > (y-1) and zc[k] < (z + 1) and zc[k] > (z - 1):
                            count += 1
                        else:
                            count += 0

                volume = 8.0
                density = float(count/volume)
                if z <= 0 and density == 0:
                    conflict = True
                    break
    
            zgrid[xdir][ydir] += (z)/nconf
            ygrid[xdir][ydir] = (y-yc[natoms])
            xgrid[xdir][ydir] = (x-xc[natoms])
                
            print xgrid[xdir][ydir],ygrid[xdir][ydir],zgrid[xdir][ydir]

for i in range(0,iskip+1):
    print "skip",i
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
    zcm=[0]*npoly
    xcm=[0]*npoly
    ycm=[0]*npoly
    sys.stdin.readline()
    line = sys.stdin.readline() #xbox bounds  
    [xm,xp] = map(float,line.split())
    line = sys.stdin.readline() #ybox bounds
    [ym,yp] = map(float,line.split())
    line = sys.stdin.readline() #zbox bounds
    [zm,zp] = map(float,line.split())
    line = sys.stdin.readline()
    for j in range(1,dim):
        line = sys.stdin.readline()
        #[ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line) #atom number, molecule number, bead type, x,y,z, image flags

    conf = conf + 1
    #print conf

print "reading config file..."
istart = iskip+1

OUT = open(file, 'w')
OUT.write("Lateral X   Lateral Y     Height\n")

for i in range(istart,nconf+iskip+1):
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
    zcm=[0]*npoly
    xcm=[0]*npoly
    ycm=[0]*npoly
    rg=[0]*npoly
    totavgsquaredrg = 0
    sys.stdin.readline()
    line = sys.stdin.readline() #xbox bounds  
    [xm,xp] = map(float,line.split())
    line = sys.stdin.readline() #ybox bounds
    [ym,yp] = map(float,line.split())
    line = sys.stdin.readline() #zbox bounds
    [zm,zp] = map(float,line.split())
    line = sys.stdin.readline()
    xbox = xp - xm
    ybox = yp - ym
    zbox = zp - zm
    vol = xbox*ybox*zbox
    xbox2 = xbox/2.0
    ybox2 = ybox/2.0
    zbox2 = zbox/2.0
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line) #atom number, molecule number, bead type, x,y,z, image flags
        k=int(ii) #sets the index k equal to the atom number
        typea[k]=int(typej)
        mol[k] = int(molj) #sets the chain number for each bead equal to molecule number
        xc[k] = xbox*float(x1) + int(n1)*xbox
        yc[k] = ybox*float(x2) + int(n2)*ybox
        zc[k] = zbox*float(x3) + int(n3)*zbox
            #should have assigned each bead a chain number
    
    radavght()
    print "Configuration complete"
    conf = conf + 1
    
for j in range(len(xgrid)):
    for k in range(len(ygrid)):
        OUT.writelines("%8.4f %8.4f %8.4f\n" % (xgrid[j][k], ygrid[j][k], zgrid[j][k]))

OUT.close()


