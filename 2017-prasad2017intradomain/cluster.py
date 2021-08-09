#!/usr/bin/python

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
from random import *
from itertools import chain
from heapq import *

nconf = 1
iskip = 0
deltr = 2.0   
ncylinder = 12
nbeads = 100
npoly = 2400
fA = 0.25
numA = int(fA*nbeads)

for i in range(0,iskip):
    print i
    sys.stdin.readline()
    sys.stdin.readline()
    sys.stdin.readline()
    line = sys.stdin.readline()
    fields = string.split(line)

    if i==0:
        natoms = int(fields[0])
        dim=natoms+1
        xc=zeros(dim,float32)
        yc=zeros(dim,float32)
        zc=zeros(dim,float32)
        xu=zeros(dim,float32)
        yu=zeros(dim,float32)
        zu=zeros(dim,float32)
        ix=zeros(dim)
        iy=zeros(dim)
        iz=zeros(dim)
        typea=[0]*dim
        mol=[0]*dim
        flag=[0]*dim
        conf = 0
        
    sys.stdin.readline()
    line = sys.stdin.readline()      
    [xm,xp] = map(float,line.split())
    line = sys.stdin.readline()      
    [ym,yp] = map(float,line.split())
    line = sys.stdin.readline()      
    [zm,zp] = map(float,line.split())
    line = sys.stdin.readline()
    xbox = xp - xm
    ybox = yp - ym
    zbox = zp - zm
    vol = xbox*ybox*zbox
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line)

for i in range(iskip,nconf+iskip):
    #print i
    sys.stdin.readline()
    line = sys.stdin.readline()
    fields = string.split(line)
    timestep = int(fields[0])
    sys.stdin.readline()
    line = sys.stdin.readline()
    fields = string.split(line)
    if i==0:
        natoms = int(fields[0])
        dim=natoms+1
        xc=zeros(dim,float32)
        yc=zeros(dim,float32)
        zc=zeros(dim,float32)
        xu=zeros(dim,float32)
        yu=zeros(dim,float32)
        zu=zeros(dim,float32)
        ix=zeros(dim)
        iy=zeros(dim)
        iz=zeros(dim)
        typea=[0]*dim
        mol=[0]*dim
        flag=[0]*dim
        conf = 0
    sys.stdin.readline()
    line = sys.stdin.readline()      
    [xm,xp] = map(float,line.split())
    line = sys.stdin.readline()      
    [ym,yp] = map(float,line.split())
    line = sys.stdin.readline()      
    [zm,zp] = map(float,line.split())
    line = sys.stdin.readline()
    xbox = xp - xm
    ybox = yp - ym
    zbox = zp - zm
    xbox2 = xbox/2
    ybox2 = ybox/2
    zbox2 = zbox/2   
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line)
        k = int(ii)
        flag[k] = k
        mol[k] = int(molj)
        typea[k] = int(typej)
        ix[k] = int(n1)
        iy[k] = int(n2)
        iz[k] = int(n3)
        xu[k] = float(x1)
        yu[k] = float(x2)
        zu[k] = float(x3)
        xc[k] = (float(x1)-0.5)*xbox
        yc[k] = (float(x2)-0.5)*ybox
        zc[k] = (float(x3)-0.5)*zbox
    for l in range(0,npoly):
        for q in range(1,numA+1):
            p=l*nbeads+q
            print p
            for m in range(p+1,dim):
                if typea[m] == 1:
                     rsq = sqrt(((xc[p]-xc[m])-xbox*round((xc[p]-xc[m])/xbox))**2+((yc[p]-yc[m])-ybox*round((yc[p]-yc[m])/ybox))**2+((zc[p]-zc[m])-zbox*round((zc[p]-zc[m])/zbox))**2)
                     if rsq < deltr:
                         if flag[m] != flag[p]:
                             it = flag[m]
                             for ic in xrange(1,dim):
                                 if typea[ic]==1 and flag[ic]==it:
                                     flag[ic]=flag[p]                       
    print nlargest(ncylinder+1,sorted(set(flag)))
    
    for n in range(1,dim):
        if typea[n]==1:
            for kk in range(0,ncylinder):
                if flag[n]==nlargest(ncylinder+1,sorted(set(flag)))[kk]:
                    flag[n]=kk+3                    
    print nlargest(ncylinder+1,sorted(set(flag)))


INPUT_LAMMPS = open('equil_clustered_skip2998.lammpstrj', 'w')
INPUT_LAMMPS.write("ITEM: TIMESTEP\n")
INPUT_LAMMPS.write("%i\n" % (timestep))
INPUT_LAMMPS.write("ITEM: NUMBER OF ATOMS\n")
INPUT_LAMMPS.write("%i\n" % (natoms))
INPUT_LAMMPS.write("ITEM: BOX BOUNDS pp pp pp\n")
INPUT_LAMMPS.write("%.4f %.4f\n" % (-xbox2,xbox2))
INPUT_LAMMPS.write("%.4f %.4f\n" % (-ybox2,ybox2))
INPUT_LAMMPS.write("%.4f %.4f\n" % (-zbox2,zbox2))
INPUT_LAMMPS.write("ITEM: ATOMS id mol type xs ys zs ix iy iz \n")
for i in range(1,dim):
    INPUT_LAMMPS.write("%i %i %i %.8f %.8f %.8f %i %i %i\n" % (i, mol[i], flag[i], xu[i], yu[i], zu[i], ix[i], iy[i], iz[i]))
