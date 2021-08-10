#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

# output_Z1: outputs the read in chains in a format appropriate for the Z1 algorithm
# see: http://polyphys-s01.ethz.ch/cgi-bin/Z1

import sys,string
from numpy import *
from math import *
import numpy as np

nconf = 500 #number of configurations to average over
iskip = 0
npoly = 2256
nbeads = 160
NP = 12    #number of nanoparticles
R = 5    #radius of NP
fname = 'Eps3_%d.Z1'
#Initial conditions
conf = 0
natoms = npoly*nbeads+NP

def output_Z1(filename):
	# if there's a filename redirect standard out there	
    f = open(filename,'w')
    sys.stdout = f
	
	# may have to shift so that the origin is at the center of the box
    x_shift = 0#-(xp-xm)/2
    y_shift = 0#-(yp-ym)/2 
    z_shift = 0#-(zp-zm)/2 

	# M
    print npoly+1

	# box_x box_y box_z
    print Lx, Ly, Lz

	# N1 N2 N3 .. NM
    print mols

    for j in range(1,dim):
        if typea[j] < 3:
            print xc[j]+x_shift, yc[j]+y_shift, zc[j]+z_shift
    for k in range(1,dim):
        if typea[k] == 3:
            print xc[k]+x_shift, yc[k]+y_shift, zc[k]+z_shift 
	# close the file when done
    if filename:
        f.close()
        sys.stdout = sys.__stdout__
        
for i in range(0,iskip+1):
    #print "skip",i
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

conf = 0
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
    Lx = xbox
    Ly = ybox
    Lz = zbox
    Lx2 = Lx/2
    Ly2 = Ly/2
    Lz2 = Lz/2
    polymer =[]
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3] = string.split(line) #atom number, molecule number, bead type, x,y,z, image flags
        k=int(ii) #sets the index k equal to the atom number
        typea[k]=int(typej)
        mol[k] = int(molj) #sets the chain number for each bead equal to molecule number
        xc[k] = float(x1)+xm
        yc[k] = float(x2)+ym
        zc[k] = float(x3)+zm 
    for m in range(1,dim):
        if typea[m] < 3 and mol[m] != mol[m-1]:
            polymer.append(nbeads)
    polymer.append(NP)
    #for k in range(1,dim):
    #    if (xc[k] > Lx):
    #        cx[k] = int(xc[k]/xbox)
    #        xc[k] = xc[k] - cx[k]*xbox
    #    elif (xc[k] < 0.0):
    #        cx[k] = -int((-xc[k]+xbox)/xbox)
    #        xc[k] = xc[k] - cx[k]*xbox
    #    else:
    #        cx[k] = 0
    #        xc[k] = xc[k]
    #    if (yc[k] > Ly):
    #        cy[k] = int(yc[k]/ybox)
    #        yc[k] = yc[k] - cy[k]*ybox
    #    elif (yc[k] < 0.0):
    #        cy[k] = -int((-yc[k]+ybox)/ybox)
    #        yc[k] = yc[k] - cy[k]*ybox 
    #    else:
    #        cy[k] = 0
    #        yc[k] = yc[k]
    #    if (zc[k] > Lz):
    #        cz[k] = int(zc[k]/zbox)
    #        zc[k] = zc[k] - cz[k]*zbox 
    #    elif (zc[k] < 0.0):
    #        cz[k] = -int((-zc[k]+zbox)/zbox)
    #        zc[k] = zc[k] - cz[k]*zbox 
    #    else:
    #        cz[k] = 0
    #        zc[k] = zc[k]

    mols = ' '.join(str(e) for e in polymer)
    filename = fname %(i-iskip)
    output_Z1(filename)
    #print vol



