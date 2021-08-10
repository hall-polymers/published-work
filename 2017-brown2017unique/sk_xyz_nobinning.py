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

nconf =   47
iskip =    0 
dr    =    0.2  

file  = 'Sk_xyz_nobinning'

def Sk():
    jj=0
    for ll in range(0,numk): 
        for mm in range(0,numk):
            for nn in range(0,numk):
                print ll,mm,nn
                s1c=s1s=s2c=s2s=s2pc=s2ps=0      
                kx=ll*dkx
                ky=mm*dky
                kz=nn*dkz
                kk=sqrt(kx*kx+ky*ky+kz*kz)
                if kk > 0:
                    jj=jj+1
                    for ii in range(1,dim): 
                        kdotr=kx*xc[ii]+ky*yc[ii]+kz*zc[ii]
                        if (typea[ii] == 1):
                            s1c = s1c + cos(kdotr)
                            s1s = s1s + sin(kdotr)
                        elif (typea[ii] == 2):
                            s2c = s2c + cos(kdotr)
                            s2s = s2s + sin(kdotr)
                            if (ii%40)>=25 or (ii%40)==0:  #pure B
                                s2pc = s2pc + cos(kdotr)
                                s2ps = s2ps + sin(kdotr)
                    kcc[jj]=kk
                    Sk11[jj]+=(s1c*s1c+s1s*s1s)/num1
                    Sk12[jj]+=(s1c*s2c+s1s*s2s)/((num1+num2)/2)
                    Sk22[jj]+=(s2c*s2c+s2s*s2s)/num2
                    Sk22p[jj]+=(s2pc*s2pc+s2ps*s2ps)/num2p

natoms = 0
num1 = 0
num2 = 0
num2p = 0
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
    xbox2 = xbox/2.0
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line)
        k=int(ii)
        typea[k] = int(typej)
        if i==0:
            if typea[k] == 1:
                num1=num1+1
            elif typea[k] == 2:
                num2=num2+1
                if (k%40)>=25 or (k%40)==0:  #pure B
                    num2p=num2p+1
        mol[k] = int(molj)
        xc[k] = xbox*(float(x1)-0.5) 
        yc[k] = ybox*(float(x2)-0.5)
        zc[k] = zbox*(float(x3)-0.5)

dkx = 2*pi/xbox
dky = 2*pi/ybox
dkz = 2*pi/zbox
numk = 10
numg1=numk+1
point = numk**3

kountSk = 0
kcc = zeros(point,float32)
Sk11 = zeros(point,float32)  
Sk12 = zeros(point,float32)
Sk22 = zeros(point,float32)
Sk22p = zeros(point,float32)
SS11 = zeros(point,float32)  
SS12 = zeros(point,float32)
SS22 = zeros(point,float32)
SS22p = zeros(point,float32)
SStot = zeros(point,float32)

Sk()
kountSk = 1

OUT = open(file, 'w')
OUT.write("#%7i\n" % (kountSk))
OUT.write("k S11 S12 S22 S22p Stot\n")
for ig in range(1,point):
    SStot[ig] = Sk11[ig] + 2*Sk12[ig] + Sk22[ig]
    OUT.write("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n" % (kcc[ig],Sk11[ig],Sk12[ig],Sk22[ig],Sk22p[ig],SStot[ig]))
OUT.close()

print "Reading config file...."

istart = iskip + 1
for kconf in range(istart,nconf+iskip):
    print kconf
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
    xbox2 = xbox/2.0
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line)
        k=int(ii)
        typea[k] = int(typej)
        mol[k] = int(molj)
        xc[k] = xbox*(float(x1)-0.5)
        yc[k] = ybox*(float(x2)-0.5)
        zc[k] = zbox*(float(x3)-0.5)

    Sk()
    kountSk = kountSk + 1
    for i in range(0,point):
        SS11[i] = Sk11[i] / nconf
        SS12[i] = Sk12[i] / nconf
        SS22[i] = Sk22[i] / nconf
        SS22p[i] = Sk22p[i] / nconf

    OUT = open(file, 'w')
    OUT.write("#%7i\n" % (kountSk))
    OUT.write("k S11 S12 S22 S22p Stot\n")
    for ig in range(1,point):
        SStot[ig] = SS11[ig] + 2*SS12[ig] + SS22[ig]
        OUT.write("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n" % (kcc[ig],SS11[ig],SS12[ig],SS22[ig],SS22p[ig],SStot[ig]))
    OUT.close()


