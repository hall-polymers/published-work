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

nconf =   30 
iskip =   1969 
mg    =   80
numg1 = mg+1
d1 = zeros(numg1,float32) 
d2 = zeros(numg1,float32)

file  = 'denprof.txt'

def den():
    for ii in range(1,natoms+1):
        zr = zc[ii]
        if zr >= zcm:
            ig = int((zr - zcm) / drg) + 1
            if ig < len(d1):
                if typea[ii] == 1:
                    d1[ig] = d1[ig] + 1.0
                elif typea[ii] == 2:
                    d2[ig] = d2[ig] + 1.0
        elif zr < zcm:
            ig = int((zr+zbox - zcm) / drg) + 1
            if ig < len(d1):
                if typea[ii] == 1:
                    d1[ig] = d1[ig] + 1.0
                elif typea[ii] == 2:
                    d2[ig] = d2[ig] + 1.0

natoms = 0
num1 = 0
num2 = 0
cirxcm = 0
cirycm = 0
thetacm = 0
zcm = 0
conf = 0

for i in range(0,iskip+1): 
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
        zc1=zeros(dim,float32)
        theta=zeros(dim,float32)
        cirx=zeros(dim,float32)
        ciry=zeros(dim,float32)
        cx=zeros(dim)
        cy=zeros(dim)
        cz=zeros(dim)
        typea=[0]*dim
        mol=[0]*dim
        
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
    zbox4 = zbox/4.0
    cirxcm = 0  
    cirycm = 0
    for j in range(1,dim):
        line = sys.stdin.readline()


print "Reading config file...."
istart = iskip+1

for kconf in range(istart,nconf+iskip+1): 
    sys.stdin.readline()
    sys.stdin.readline()     
    sys.stdin.readline()
    line = sys.stdin.readline()    
    fields = string.split(line)
    num = fields[0]
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
    zbox4 = zbox/4.0
    cirxcm = 0   
    cirycm = 0
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line)
        k=int(ii)
        typea[k] = int(typej)
        mol[k] = int(molj)
        xc[k] = xbox*float(x1) 
        yc[k] = ybox*float(x2)
        zc[k] = zbox*float(x3)
        if zc[k] <= zbox4:           
            zc1[k] = zc[k]
        elif zc[k] > zbox4 and zc[k] <= 2*zbox4:
            zc1[k] = zc[k] - zbox4
        elif zc[k] > 2*zbox4 and zc[k] <= 3*zbox4:
            zc1[k] = zc[k] - 2*zbox4
        elif zc[k] > 3*zbox4 and zc[k] <= 4*zbox4:
            zc1[k] = zc[k] - 3*zbox4
        theta[k] = (2*pi/zbox4)*zc1[k] 
        cirx[k] = (zbox4/(2*pi))*sin(theta[k]) 
        ciry[k] = (zbox4/(2*pi))*cos(theta[k])
        
        if typea[k] == 1:
            cirxcm += cirx[k]/(natoms*0.5)
            cirycm += ciry[k]/(natoms*0.5)
        thetacm = atan2(-cirxcm,-cirycm) + pi 
        zcm = (zbox4/(2*pi))*thetacm


    drg = zbox/mg

    den()
    
    conf = conf + 1
    dd1 = d1 / nconf
    dd2 = d2 / nconf

    OUT = open(file, 'w')
    OUT.write("#%7i\n" % (conf))
    OUT.write("z num(A) num(B)\n")
    for ig in range(1,numg1):
        OUT.write("%8.4f %8.4f %8.4f\n" % (ig,dd1[ig],dd2[ig]))
    OUT.close()



