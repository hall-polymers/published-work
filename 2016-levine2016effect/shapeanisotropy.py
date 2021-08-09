#!/opt/local/bin/python

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

nconf = 30 
iskip = 1969
npoly = 480
nbeads = 40
conf = 0
totavgsa = 0

file = 'ShapeAnisotropy.txt'

istart = 0
avgsa = [0]*npoly
sa = [0]*npoly

def shapeanisotropy():
    for jj in range(0,npoly):
        for kk in range(1,nbeads+1):
           k = nbeads*jj + kk
           zcm[jj] += zc[k]/nbeads
           ycm[jj] += yc[k]/nbeads
           xcm[jj] += xc[k]/nbeads
    for jj in range(0,npoly):
        for kk in range(1,nbeads+1):
            k = nbeads*jj + kk
            ixx[jj] += (yc[k]-ycm[jj])**2+(zc[k]-zcm[jj])**2
            iyy[jj] += (xc[k]-xcm[jj])**2+(zc[k]-zcm[jj])**2
            izz[jj] += (xc[k]-xcm[jj])**2+(yc[k]-ycm[jj])**2
            ixy[jj] += (xc[k]-xcm[jj])*(yc[k]-ycm[jj])
            ixz[jj] += (xc[k]-xcm[jj])*(zc[k]-zcm[jj])
            iyz[jj] += (yc[k]-ycm[jj])*(zc[k]-zcm[jj])  
    for jj in range(0,npoly):
        a[jj] = [[ixx[jj],ixy[jj],ixz[jj]],[ixy[jj],iyy[jj],iyz[jj]],[ixz[jj],iyz[jj],izz[jj]]]
        eigenvalues[jj] = linalg.eigvals(a[jj])
        maxeigvalue[jj] = max(eigenvalues[jj])
        mineigvalue[jj] = min(eigenvalues[jj])
        mideigvalue[jj] = median(eigenvalues[jj])
    for jj in range(0,npoly):
        sa[jj] = 4-12*(((maxeigvalue[jj]*mideigvalue[jj])+(maxeigvalue[jj]*mineigvalue[jj])+(mideigvalue[jj]*mineigvalue[jj]))/((maxeigvalue[jj]+mideigvalue[jj]+mineigvalue[jj])**2))
     
for i in range(0,iskip+1):
    sys.stdin.readline()
    sys.stdin.readline() 
    sys.stdin.readline()
    line = sys.stdin.readline()
    fields = string.split(line) 
    natoms = int(fields[0])
    dim=natoms+1
    xc=zeros(dim,float32)
    yc=zeros(dim,float32)
    zc=zeros(dim,float32)
    zci=zeros(dim,float32)
    cx=zeros(dim)
    cy=zeros(dim)
    cz=zeros(dim)
    typea=[0]*dim
    mol=[0]*dim
    chain=[0]*dim
    imf=[0]*dim
    zcm=[0]*npoly
    xcm=[0]*npoly
    ycm=[0]*npoly
    rgxx=[0]*npoly
    ixx=[0]*npoly
    rgyy=[0]*npoly
    iyy=[0]*npoly
    rgzz=[0]*npoly
    izz=[0]*npoly
    rgxy=[0]*npoly
    ixy=[0]*npoly
    rgxz=[0]*npoly
    ixz=[0]*npoly
    rgyz=[0]*npoly
    iyz=[0]*npoly
    a=[[0,0,0],[0,0,0],[0,0,0]]*npoly
    eigenvalues=[0]*npoly
    maxeigvalue = [0]*npoly
    mineigvalue = [0]*npoly
    mideigvalue = [0]*npoly
    sys.stdin.readline()
    line = sys.stdin.readline() 
    [xm,xp] = map(float,line.split())
    line = sys.stdin.readline() 
    [ym,yp] = map(float,line.split())
    line = sys.stdin.readline() 
    [zm,zp] = map(float,line.split())
    line = sys.stdin.readline()
    for j in range(1,dim):
        line = sys.stdin.readline()

    conf = conf + 1
    print conf

print "reading config file..."
istart = iskip+1

for i in range(istart,nconf+iskip+1):
    sys.stdin.readline()
    sys.stdin.readline() 
    sys.stdin.readline()
    line = sys.stdin.readline() 
    fields = string.split(line) 
    natoms = int(fields[0])
    dim=natoms+1
    xc=zeros(dim,float32)
    yc=zeros(dim,float32)
    zc=zeros(dim,float32)
    zci=zeros(dim,float32)
    cx=zeros(dim)
    cy=zeros(dim)
    cz=zeros(dim)
    typea=[0]*dim
    mol=[0]*dim
    chain=[0]*dim
    imf=[0]*dim
    zcm=[0]*npoly
    xcm=[0]*npoly
    ycm=[0]*npoly
    rgxx=[0]*npoly
    rgyy=[0]*npoly
    rgzz=[0]*npoly
    rgxy=[0]*npoly
    rgxz=[0]*npoly
    rgyz=[0]*npoly
    a=[[0,0,0],[0,0,0],[0,0,0]]*npoly
    eigenvalues=[0]*npoly
    maxeigvalue = [0]*npoly
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
    xbox2 = xbox/2.0
    ybox2 = ybox/2.0
    zbox2 = zbox/2.0
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line)
        k=int(ii)
        chain[k] = int(molj)
        xc[k] = xbox*float(x1) 
        yc[k] = ybox*float(x2)
        zc[k] = zbox*float(x3)
        xc[k] = xbox*float(x1) + int(n1)*xbox
        yc[k] = ybox*float(x2) + int(n2)*ybox
        zc[k] = zbox*float(x3) + int(n3)*zbox
            
    shapeanisotropy()
    conf = conf + 1
    print conf
    for jj in range(0,npoly):
        avgsa[jj] += sa[jj]/nconf


OUT = open(file, 'w')
OUT.write("#%7i\n" % (conf))
OUT.write("polymer# SA\n")
for jj in range(0, npoly):
    totavgsa += avgsa[jj]/npoly
    OUT.write("%8.4f %8.4f\n" % (jj, avgsa[jj]))
OUT.write("Total average")
OUT.write("%8.4f\n" % (totavgsa))
OUT.close()
