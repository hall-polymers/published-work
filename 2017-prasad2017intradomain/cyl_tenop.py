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

nconf = 399 
iskip = 100
mg = 20 
nbeads = 100
npoly = 2400
fA = 0.25
ncyl = 12

file  = 'cyl_tenOP.txt'

def bond():
    for ii in range(0,npoly):
        for jj in range(2,nbeads):
            kk = ii*nbeads+jj
            inter = ii*nbeads+int(nbeads*fA)
            if typea[kk] == 2:
                xr = (xu[kk]-xu[kk+1])
                yr = (yu[kk]-yu[kk+1])
                zr = (zu[kk]-zu[kk+1])
                xmid = (xu[kk]+xu[kk+1])/2
                ymid = (yu[kk]+yu[kk+1])/2
                if xmid > xbox:               
                    xmid = xmid - xbox*int(xmid/xbox)
                elif xmid < 0:
                    xmid = xmid + xbox*int((-xmid+xbox)/xbox)
                if ymid > ybox:
                    ymid = ymid - ybox*int(ymid/ybox)
                elif ymid < 0:
                    ymid = ymid + ybox*int((-ymid+ybox)/ybox)
                pp = typea[ii*nbeads+1]-3
                xradial = (xmid-xcom[pp])-xbox*round((xmid-xcom[pp])/xbox)
                yradial = (ymid-ycom[pp])-ybox*round((ymid-ycom[pp])/ybox)
                zradial = 0
                uxradial = xradial/sqrt(xradial**2+yradial**2+zradial**2)
                uyradial = yradial/sqrt(xradial**2+yradial**2+zradial**2)
                thetatan = arctan(-uxradial/uyradial)
                uxtan = cos(thetatan)
                uytan = sin(thetatan)
                xv = (xr/sqrt(xr**2+yr**2+zr**2))
                yv = (yr/sqrt(xr**2+yr**2+zr**2))
                zv = (zr/sqrt(xr**2+yr**2+zr**2))
                trr = xv*uxradial + yv*uyradial
                ttheta = xv*uxtan + yv*uytan
                for m in range(0,ncyl):
                    ig = int(sqrt(((xcmid[pp][kk] - xcom[m])-xbox*round((xcmid[pp][kk] - xcom[m])/xbox))**2+((ycmid[pp][kk] - ycom[m])-ybox*round((ycmid[pp][kk] - ycom[m])/ybox))**2)/drg) + 1
                    if ig < len(d1r):
                        d2r[ig] = d2r[ig] + trr**2 - float(1./3.)          #q_rr
                        d2t[ig] = d2t[ig] + ttheta**2 - float(1./3.)       #q_tt
                        d2z[ig] = d2z[ig] + zv**2 - float(1./3.)           #q_zz
                        nm1[ig] = nm1[ig] + 1.0
            else:
                xr = (xu[kk]-xu[kk-1])
                yr = (yu[kk]-yu[kk-1])
                zr = (zu[kk]-zu[kk-1])
                xmid = (xu[kk]+xu[kk-1])/2
                ymid = (yu[kk]+yu[kk-1])/2
                if xmid > xbox:                
                    xmid = xmid - xbox*int(xmid/xbox)
                elif xmid < 0:
                    xmid = xmid + xbox*int((-xmid+xbox)/xbox)
                if ymid > ybox:               
                    ymid = ymid - ybox*int(ymid/ybox)
                elif ymid < 0:
                    ymid = ymid + ybox*int((-ymid+ybox)/ybox)
                pp = typea[kk]-3
                xradial = (xmid-xcom[pp])-xbox*round((xmid-xcom[pp])/xbox)
                yradial = (ymid-ycom[pp])-ybox*round((ymid-ycom[pp])/ybox)
                zradial = 0
                uxradial = xradial/sqrt(xradial**2+yradial**2+zradial**2)
                uyradial = yradial/sqrt(xradial**2+yradial**2+zradial**2)
                thetatan = arctan(-uxradial/uyradial)
                uxtan = cos(thetatan)
                uytan = sin(thetatan)
                xv = (xr/sqrt(xr**2+yr**2+zr**2))
                yv = (yr/sqrt(xr**2+yr**2+zr**2))
                zv = (zr/sqrt(xr**2+yr**2+zr**2))
                trr = xv*uxradial + yv*uyradial
                ttheta = xv*uxtan + yv*uytan
                for m in range(0,ncyl):
                    ig = int(sqrt(((xcmid[pp][kk] - xcom[m])-xbox*round((xcmid[pp][kk] - xcom[m])/xbox))**2+((ycmid[pp][kk] - ycom[m])-ybox*round((ycmid[pp][kk] - ycom[m])/ybox))**2)/drg) + 1
                    if ig < len(d1r):
                        d1r[ig] = d1r[ig] + trr**2 - float(1./3.)         #q_rr
                        d1t[ig] = d1t[ig] + ttheta**2 - float(1./3.)      #q_tt
                        d1z[ig] = d1z[ig] + zv**2 - float(1./3.)          #q_zz
                        nm1[ig] = nm1[ig] + 1.0

    for ig in range(1,numg1):
        if nm1[ig] != 0:
            d1r[ig] = d1r[ig]/nm1[ig]
            d1t[ig] = d1t[ig]/nm1[ig]
            d1z[ig] = d1z[ig]/nm1[ig]
            d2r[ig] = d2r[ig]/nm1[ig]
            d2t[ig] = d2t[ig]/nm1[ig]
            d2z[ig] = d2z[ig]/nm1[ig]

#initialize
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
            
drg = zbox/mg
numg = mg
numg1 = numg+1

dd1r = zeros(numg1,float32)
dd1t = zeros(numg1,float32)
dd1z = zeros(numg1,float32)
dd2r = zeros(numg1,float32)
dd2t = zeros(numg1,float32)
dd2z = zeros(numg1,float32)

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
    d1r = zeros(numg1,float32)
    d1t = zeros(numg1,float32)
    d1z = zeros(numg1,float32)
    d2r = zeros(numg1,float32)
    d2t = zeros(numg1,float32)
    d2z = zeros(numg1,float32)
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
    dd1r += d1r / nconf
    dd1t += d1t / nconf
    dd1z += d1z / nconf
    dd2r += d2r / nconf
    dd2t += d2t / nconf
    dd2z += d2z / nconf

    OUT = open(file, 'w')
    OUT.write("#%7i\n" % (conf))
    OUT.write("r Qrr(A) Qtt(A) Qzz(A) Qrr(B) Qtt(B) Qzz(B)\n")
    for ig in range(1,numg1):
        OUT.write("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n" % (ig,dd1r[ig],dd1t[ig],dd1z[ig],dd2r[ig],dd2t[ig],dd2z[ig]))
    OUT.close()


