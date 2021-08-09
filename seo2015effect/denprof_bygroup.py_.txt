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

nconf =    30 
npoly =    480
nbeads =   40
hnbeads =  int(nbeads/2.)
iskip =    1969 
mg    =    80 
numb1 = [0]*npoly
numb2 = [0]*npoly
binvol = 0
conf = 0

d11 = zeros(mg+1,float32)
d12 = zeros(mg+1,float32)
d21 = zeros(mg+1,float32)
d22 = zeros(mg+1,float32)
d31 = zeros(mg+1,float32)
d32 = zeros(mg+1,float32)

file  = 'denprof_bygroup.txt'

def classifytaper():
    for ii in range(0,npoly):
        for k in range(0,nbeads):
            sequence1[k] = typea[k+1+nbeads*ii]
            if sequence1[k] == 1:
                numb1[ii] = numb1[ii] + 1
            elif sequence1[k] == 2:
                numb2[ii] = numb2[ii] + 1
        sequence2[ii] = sequence1
        if numb1[ii] > (hnbeads+10):
            numb1[ii] = numb1[ii]/2

def circling():
    cirxcm = 0
    cirycm = 0
    thetacm = 0
    zcm = 0
    for ii in range(0,natoms):
        cx = cirx[ii+1]
        cy = ciry[ii+1]
        if typea[ii+1] == 1:
            cirxcm += cx/(natoms*0.5)
            cirycm += cy/(natoms*0.5)            
        thetacm = atan2(-cirxcm,-cirycm) + pi
        zcm = (zbox4/(2*pi))*thetacm
        
    for ii in range(0,natoms):
        jj = int(ii/nbeads)
        binvol = xbox*ybox*drg
        if numb1[jj] > hnbeads:
            if zc[ii+1] >= zcm:
                ig = int((zc[ii+1] - zcm) / drg) + 1
                if typea[ii+1] == 1:
                    d11[ig] = d11[ig] + 1.0
                elif typea[ii+1] == 2:
                    d12[ig] = d12[ig] + 1.0
            elif zc[ii+1] < zcm:
                ig = int((zc[ii+1]+zbox - zcm) / drg) + 1
                if typea[ii+1] == 1:
                    d11[ig] = d11[ig] + 1.0
                elif typea[ii+1] == 2:
                    d12[ig] = d12[ig] + 1.0
        elif numb1[jj] == hnbeads:
            if zc[ii+1] >= zcm:
                ig = int((zc[ii+1] - zcm) / drg) + 1
                if typea[ii+1] == 1:
                    d21[ig] = d21[ig] + 1.0
                elif typea[ii+1] == 2:
                    d22[ig] = d22[ig] + 1.0
            elif zc[ii+1] < zcm:
                ig = int((zc[ii+1]+zbox - zcm) / drg) + 1
                if typea[ii+1] == 1:
                    d21[ig] = d21[ig] + 1.0
                elif typea[ii+1] == 2:
                    d22[ig] = d22[ig] + 1.0
        elif numb1[jj] < hnbeads:
            if zc[ii+1] >= zcm:
                ig = int((zc[ii+1] - zcm) / drg) + 1
                if typea[ii+1] == 1:
                    d31[ig] = d31[ig] + 1.0
                elif typea[ii+1] == 2:
                    d32[ig] = d32[ig] + 1.0
            elif zc[ii+1] < zcm:
                ig = int((zc[ii+1]+zbox - zcm) / drg) + 1
                if typea[ii+1] == 1:
                    d31[ig] = d31[ig] + 1.0
                elif typea[ii+1] == 2:
                    d32[ig] = d32[ig] + 1.0
            
natoms = 0
num1 = 0
num2 = 0

for i in range(0,iskip+1):
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
        zc1=zeros(dim,float32)
        cx=zeros(dim)
        cy=zeros(dim)
        cz=zeros(dim)
        typea=[0]*dim
        mol=[0]*dim
        theta=[0]*dim
        cirx=[0]*dim
        ciry=[0]*dim
        sequence1=[0]*nbeads
        sequence2=[0]*npoly
        jj=0
        
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
    zbox2 = zbox/2
    zbox4 = zbox/4   
    for j in range(1,dim):
        line = sys.stdin.readline()

print "Reading config file...."

istart = iskip + 1
for kconf in range(istart,nconf+iskip+1):
    print kconf
    sys.stdin.readline()
    sys.stdin.readline()     
    sys.stdin.readline()
    line = sys.stdin.readline()    
    fields = string.split(line)
    num = int(fields[0])
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
    zbox4 = zbox/4
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

    drg = zbox/mg
    numg1 = mg+1
    
    classifytaper()
    circling()
    conf = conf + 1
    dd11 = d11 / nconf
    dd12 = d12 / nconf
    dd21 = d21 / nconf
    dd22 = d22 / nconf
    dd31 = d31 / nconf
    dd32 = d32 / nconf

    OUT = open(file, 'w')
    OUT.write("configuration number\n")
    OUT.write("#%7i\n" % (conf))
    OUT.write("z den(A1) den(B1) den(A2) den(B2) den(A3) den(B3)\n")
    for ig in range(1,numg1):
        OUT.write("%8.0f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n" % (ig,dd11[ig],dd12[ig],dd21[ig],dd22[ig],dd31[ig],dd32[ig]))
    OUT.close()


