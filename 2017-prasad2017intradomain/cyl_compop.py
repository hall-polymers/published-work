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

nconf = 100
iskip = 399
mg = 20
nbeads = 100
npoly = 2400
fA = 0.25
ncyl = 12

file  = 'cyl_compOP.txt'

def bond():
    for ii in range(0,npoly):
        for jj in range(2,nbeads):
            kk = ii*nbeads+jj
            if typea[kk] == 2:
                pp = typea[ii*nbeads+1]-3
                for m in range(0,ncyl):
                    ig = int(sqrt(((xcmid[pp][kk] - xcom[m])-xbox*round((xcmid[pp][kk] - xcom[m])/xbox))**2+((ycmid[pp][kk] - ycom[m])-ybox*round((ycmid[pp][kk] - ycom[m])/ybox))**2)/drg) + 1
                    if ig < len(d1):
                        d2[ig] = d2[ig] + 1
            else:
                pp = typea[kk]-3
                for m in range(0,ncyl):
                    ig = int(sqrt(((xcmid[pp][kk] - xcom[m])-xbox*round((xcmid[pp][kk] - xcom[m])/xbox))**2+((ycmid[pp][kk] - ycom[m])-ybox*round((ycmid[pp][kk] - ycom[m])/ybox))**2)/drg) + 1
                    if ig < len(d1):
                        d1[ig] = d1[ig] + 1

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
    conf=conf+1
    print conf
    #
    if i==0:
        natoms = int(fields[0]) 
        dim=natoms+1
        xc=zeros(dim,float32)
        yc=zeros(dim,float32)
        zc=zeros(dim,float32)
        xu=zeros(dim,float32)
        yu=zeros(dim,float32)
        zu=zeros(dim,float32)
        xcmid=[[0 for x in range(dim)] for x in range(12)]
        ycmid=[[0 for x in range(dim)] for x in range(12)]
        thetax=[[0 for x in range(dim)] for x in range(12)]
        thetay=[[0 for x in range(dim)] for x in range(12)]
        cirxx=[[0 for x in range(dim)] for x in range(12)]
        cirxy=[[0 for x in range(dim)] for x in range(12)]
        ciryx=[[0 for x in range(dim)] for x in range(12)]
        ciryy=[[0 for x in range(dim)] for x in range(12)]
        cx=zeros(dim)
        cy=zeros(dim)
        cz=zeros(dim)
        typea=[0]*dim
        mol=[0]*dim
        xcom=[0]*12
        ycom=[0]*12
        #
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
            
drg = (xbox/6)/mg
numg = mg
numg1 = numg+1

dd1 = zeros(numg1,float32)
dd2 = zeros(numg1,float32)

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
    d1 = zeros(numg1,float32)
    d2 = zeros(numg1,float32)
    nm1 = zeros(numg1,float32)
    numvec=[0]*12
    cirxcmx=[0]*12
    cirycmx=[0]*12
    cirxcmy=[0]*12
    cirycmy=[0]*12
    thetacmx=[0]*12
    thetacmy=[0]*12
    xcom=[0]*12
    ycom=[0]*12
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line)
        k=int(ii)
        typea[k] = int(typej)
        mol[k] = int(molj)
        xc[k] = xbox*float(x1)
        yc[k] = ybox*float(x2)
        zc[k] = zbox*float(x3)
        xu[k] = xbox*(float(x1)+int(n1))
        yu[k] = ybox*(float(x2)+int(n2))
        zu[k] = zbox*(float(x3)+int(n3))

    for n in range(0,npoly):
        for m in range(2,int(nbeads*fA+1)):
            kk = n*nbeads+m
            pp = typea[kk]-3
            numvec[pp] += 1
            xcmid[pp][kk] = (xu[kk]+xu[kk-1])/2
            ycmid[pp][kk] = (yu[kk]+yu[kk-1])/2
            if xcmid[pp][kk] > xbox:              
                xcmid[pp][kk] = xcmid[pp][kk] - xbox*int(xcmid[pp][kk]/xbox)
            elif xcmid[pp][kk] < 0:
                xcmid[pp][kk] = xcmid[pp][kk] + xbox*int((-xcmid[pp][kk]+xbox)/xbox)
            if ycmid[pp][kk] > ybox:                 
                ycmid[pp][kk] = ycmid[pp][kk] - ybox*int(ycmid[pp][kk]/ybox)
            elif ycmid[pp][kk] < 0:
                 ycmid[pp][kk] = ycmid[pp][kk] + ybox*int((-ycmid[pp][kk]+ybox)/ybox)
            thetax[pp][kk] = (2*pi/xbox)*xcmid[pp][kk]
            thetay[pp][kk] = (2*pi/ybox)*ycmid[pp][kk]
            cirxx[pp][kk] = (xbox/(2*pi))*sin(thetax[pp][kk])
            cirxy[pp][kk] = (xbox/(2*pi))*cos(thetax[pp][kk])
            ciryx[pp][kk] = (ybox/(2*pi))*sin(thetay[pp][kk])
            ciryy[pp][kk] = (ybox/(2*pi))*cos(thetay[pp][kk])
        for l in range(int(nbeads*fA+1),nbeads):   
            kk = n*nbeads+l
            pp = typea[n*nbeads+1]-3 
            xcmid[pp][kk] = (xu[kk]+xu[kk+1])/2
            ycmid[pp][kk] = (yu[kk]+yu[kk+1])/2
            if xcmid[pp][kk] > xbox:                 
                xcmid[pp][kk] = xcmid[pp][kk] - xbox*int(xcmid[pp][kk]/xbox)
            elif xcmid[pp][kk] < 0:
                xcmid[pp][kk] = xcmid[pp][kk] + xbox*int((-xcmid[pp][kk]+xbox)/xbox)
            if ycmid[pp][kk] > ybox:                 
                ycmid[pp][kk] = ycmid[pp][kk] - ybox*int(ycmid[pp][kk]/ybox)
            elif ycmid[pp][kk] < 0:
                 ycmid[pp][kk] = ycmid[pp][kk] + ybox*int((-ycmid[pp][kk]+ybox)/ybox)

    for p in range(0,12):
        cirxcmx[p] = sum(cirxx[p])/numvec[p]
        cirycmx[p] = sum(cirxy[p])/numvec[p]
        cirxcmy[p] = sum(ciryx[p])/numvec[p]
        cirycmy[p] = sum(ciryy[p])/numvec[p]
        thetacmx[p] = atan2(-cirxcmx[p],-cirycmx[p]) + pi
        thetacmy[p] = atan2(-cirxcmy[p],-cirycmy[p]) + pi
        xcom[p] = (xbox/(2*pi))*thetacmx[p]
        ycom[p] = (ybox/(2*pi))*thetacmy[p]

        
    drg = (xbox/6)/mg
    numg = mg
    numg1 = numg+1
   
    conf = conf + 1
    print conf
    bond()
    dd1 += d1 / nconf
    dd2 += d2 / nconf

    OUT = open(file, 'w')
    OUT.write("#%7i\n" % (conf))
    OUT.write("z comp(A) comp(B)\n")
    for ig in range(1,numg1):
        OUT.write("%8.4f %8.4f %8.4f\n" % (ig,dd1[ig],dd2[ig]))
    OUT.close()


