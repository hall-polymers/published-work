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

nconf = 99
iskip = 0
mg = 160
nlayers = 4
nbeads = 100
npoly = 1600

file  = 'vectorOP.txt'

def bond():
    for ii in range(0,npoly):
        for jj in range(2,nbeads):
            xr = yr = zr = zv = 0
            zrc = 0
            if jj < (nbeads/2+1):   #A monomer
                xr = (xu[ii*nbeads+jj] - xu[ii*nbeads+jj-1])**2
                yr = (yu[ii*nbeads+jj] - yu[ii*nbeads+jj-1])**2
                zr = (zu[ii*nbeads+jj] - zu[ii*nbeads+jj-1])**2
                zv = (zu[ii*nbeads+jj] - zu[ii*nbeads+jj-1])/sqrt(xr+yr+zr)
                xv = (xu[ii*nbeads+jj] - xu[ii*nbeads+jj-1])/sqrt(xr+yr+zr)
                yv = (yu[ii*nbeads+jj] - yu[ii*nbeads+jj-1])/sqrt(xr+yr+zr)
                zrc = (zu[ii*nbeads+jj] + zu[ii*nbeads+jj-1])/2
                if zrc > zbox:                   
                    zrc = zrc - zbox*int(zrc/zbox)
                elif zrc < 0:
                    zrc = zrc + zbox*int((-zrc+zbox)/zbox)
                if zrc >= zcm:
                    ig = int((zrc - zcm) / drg) + 1
                    if ig < len(d1par):
                        d1par[ig] = d1par[ig] + xv + yv  #par
                        d1perp[ig] = d1perp[ig] + zv     #perp
                        nm1[ig] = nm1[ig] + 1.0
                elif zrc < zcm:
                    ig = int((zrc+zbox - zcm) / drg) + 1
                    if ig < len(d1par):
                        d1par[ig] = d1par[ig] + xv + yv  #par
                        d1perp[ig] = d1perp[ig] + zv     #perp
                        nm1[ig] = nm1[ig] + 1.0
            else:   #B monomer
                xr = (xu[ii*nbeads+jj] - xu[ii*nbeads+jj+1])**2
                yr = (yu[ii*nbeads+jj] - yu[ii*nbeads+jj+1])**2
                zr = (zu[ii*nbeads+jj] - zu[ii*nbeads+jj+1])**2
                zv = (zu[ii*nbeads+jj] - zu[ii*nbeads+jj+1])/sqrt(xr+yr+zr)
                xv = (xu[ii*nbeads+jj] - xu[ii*nbeads+jj-1])/sqrt(xr+yr+zr)
                yv = (yu[ii*nbeads+jj] - yu[ii*nbeads+jj-1])/sqrt(xr+yr+zr)
                zrc = (zc[ii*nbeads+jj] + zc[ii*nbeads+jj+1])/2
                if zrc > zbox:                   
                    zrc = zrc - zbox*int(zrc/zbox)
                elif zrc < 0:
                    zrc = zrc + zbox*int((-zrc+zbox)/zbox)
                if zrc >= zcm:
                    ig = int((zrc - zcm) / drg) + 1
                    if ig < len(d1par):
                        d2par[ig] = d2par[ig] + xv + yv  #par
                        d2perp[ig] = d2perp[ig] + zv     #perp
                        nm1[ig] = nm1[ig] + 1.0
                elif zrc < zcm:
                    ig = int((zrc+zbox - zcm) / drg) + 1
                    if ig < len(d1par):
                        d2par[ig] = d2par[ig] + xv + yv  #par
                        d2perp[ig] = d2perp[ig] + zv     #perp
                        nm1[ig] = nm1[ig] + 1.0


    for ig in range(1,numg1):
        d1perp[ig] = d1perp[ig]/nm1[ig]
        d1par[ig] = d1par[ig]/nm1[ig]/2
        d2perp[ig] = d2perp[ig]/nm1[ig]
        d2par[ig] = d2par[ig]/nm1[ig]/2

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
    
    if i==0:
        natoms = int(fields[0])
        dim=natoms+1
        xc=zeros(dim,float32)
        yc=zeros(dim,float32)
        zc=zeros(dim,float32)
        xu=zeros(dim,float32)
        yu=zeros(dim,float32)
        zu=zeros(dim,float32)
        zuc=zeros(dim,float32)
        zuc1=zeros(dim,float32)
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
    zboxi = zbox/nlayers
    cirxcm = 0
    cirycm = 0
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line)
            
drg = zbox/mg
numg = mg
numg1 = numg+1

dd1perp = zeros(numg1,float32)
dd1par = zeros(numg1,float32)
dd2perp = zeros(numg1,float32)
dd2par = zeros(numg1,float32)

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
    zboxi = zbox/nlayers
    cirxcm = 0
    cirycm = 0
    d1perp = zeros(numg1,float32)
    d1par = zeros(numg1,float32)
    d2perp = zeros(numg1,float32)
    d2par = zeros(numg1,float32)
    nm1 = zeros(numg1,float32)
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
        for m in range(2,nbeads/2+1): 
            kk = n*nbeads+m
            zuc[kk] = (zu[kk]+zu[kk-1])/2
            if zuc[kk] > zbox:                   
                zuc[kk] = zuc[kk] - zbox*int(zuc[kk]/zbox)
            elif zuc[kk] < 0:
                zuc[kk] = zuc[kk] + zbox*int((-zuc[kk]+zbox)/zbox)  
            for ii in range(0,nlayers):
                if zboxi*ii < zuc[kk] <= zboxi*(ii+1):
                    zuc1[kk] = zuc[kk] - zboxi*ii
                
            theta[kk] = (2*pi/zboxi)*zuc1[kk]
            cirx[kk] = (zboxi/(2*pi))*sin(theta[kk])
            ciry[kk] = (zboxi/(2*pi))*cos(theta[kk])
            cirxcm += cirx[kk]/(natoms*0.5)
            cirycm += ciry[kk]/(natoms*0.5)
    thetacm = atan2(-cirxcm,-cirycm) + pi
    zcm = (zboxi/(2*pi))*thetacm

    drg = zbox/mg
    numg = mg
    numg1 = numg+1
   
    conf = conf + 1
    print conf
    bond()
    dd1perp += d1perp / nconf
    dd1par += d1par / nconf
    dd2perp += d2perp / nconf
    dd2par += d2par / nconf

    #    Output
    OUT = open(file, 'w')
    OUT.write("#%7i\n" % (conf))
    OUT.write("z vecOP(A)_perp vecOP(B)_perp vecOP(A)_par vecOP(B)_par\n")
    for ig in range(1,numg1):
        OUT.write("%8.4f %8.4f %8.4f %8.4f %8.4f\n" % (ig,dd1perp[ig],dd2perp[ig],dd1par[ig],dd2par[ig]))
    OUT.close()

# end of loop


