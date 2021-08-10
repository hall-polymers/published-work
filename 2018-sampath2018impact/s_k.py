#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#!/usr/bin/python
 
# Script:  gr_corners.py
# Purpose: calculate average g(r)s between types 1,2,and 3, averaging over configs in dump file after iskip
# Syntax:  gr_corners.py < filename 
# Example: gr_corners.py < test.dump (dump file with scaled coordinates)
# Author:  Lisa Hall
# Now calculates into the corners of the box to L/2*sqrt(2)
#derived from Mark Stevens' code
# -------------------------------------------------------

import sys,string
from numpy import *
from math import * 

# INPUT PARAMETERs
nconf =   2 # number of configurations to average over
#polylength= 30  #now based on type and not polymer number
#npoly =    100
iskip =      0 #number of configurations to skip before taking data
dr    =    0.2  #sets number of bins in S(k)


file  = 'S_k_1'


def Sk():
    OUT = open('allS_1yz', 'a')
    OUT.write( "ll, mm, nn, kk, s1c*s1c+s1s*s1s,s1c*s2c+s1s*s2s,s1c*s3c+s1s*s3s,s1c*s4c+s1s*s4s,\
            s2c*s2c+s2s*s2s,s2c*s3c+s2s*s3s,s2c*s4c+s2s*s4s,s3c*s3c+s3s*s3s,s3c*s4c+s3s*s4s,s4c*s4c+s4s*s4s\n")   
    for ll in range(0,1): #jj indexes in q/k
      for mm in range(0,numg1):
       for nn in range(0,numg1):
        s1c=s1s=s2c=s2s=s3c=s3s=s4c=s4s=0
        kx=ll*dk
        ky=mm*dk
        kz=nn*dk
        kk=sqrt(kx*kx+ky*ky+kz*kz)
        jj=int(kk/dk)
        if jj<numg1 and jj>0:
         nk[jj]=nk[jj]+1
         kcc[jj]=kcc[jj]+kk
         for ii in range(1,dim): 
            kdotr=kx*xc[ii]+ky*yc[ii]+kz*zc[ii]          
            if typea[ii] == 1:
                s1c = s1c + cos(kdotr)
                s1s = s1s + sin(kdotr)
            elif (typea[ii] == 2):
                s2c = s2c + cos(kdotr)
                s2s = s2s + sin(kdotr)
            elif typea[ii] == 3:
                s3c = s3c + cos(kdotr)
                s3s = s3s + sin(kdotr)
            elif (typea[ii] == 4):
                s4c = s4c + cos(kdotr)
                s4s = s4s + sin(kdotr)

         S11[jj]+=s1c*s1c+s1s*s1s
         S12[jj]+=s1c*s2c+s1s*s2s
         S13[jj]+=s1c*s3c+s1s*s3s
         S14[jj]+=s1c*s4c+s1s*s4s
         S22[jj]+=s2c*s2c+s2s*s2s
         S23[jj]+=s2c*s3c+s2s*s3s
         S24[jj]+=s2c*s4c+s2s*s4s
         S33[jj]+=s3c*s3c+s3s*s3s
         S34[jj]+=s3c*s4c+s3s*s4s
         S44[jj]+=s4c*s4c+s4s*s4s
         if nn<20 and mm <20 and ll < 1:
            
             OUT.write( "%i, %i, %i, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f\n" % (ll, mm, nn, kk, s1c*s1c+s1s*s1s,s1c*s2c+s1s*s2s,s1c*s3c+s1s*s3s,s1c*s4c+s1s*s4s,\
                                                                                             s2c*s2c+s2s*s2s,s2c*s3c+s2s*s3s,s2c*s4c+s2s*s4s,s3c*s3c+s3s*s3s,s3c*s4c+s3s*s4s,s4c*s4c+s4s*s4s))

    OUT.close()


#initialize
natoms = 0
num1 = 0
num2 = 0
num3 = 0
num4 = 0
for i in range(0,iskip+1):  #reads in first configuration then reads in more over top before doing the first g(r)
                           # (to skip pushoff dump configurations, for example)
    sys.stdin.readline()
    sys.stdin.readline()      # time step
    sys.stdin.readline()
    line = sys.stdin.readline()      # number of atoms
    fields = string.split(line)
    #
    if i==0:
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
    xbox2 = xbox/2.0
    for j in range(1,dim):
        line = sys.stdin.readline()
        #todo: make a list called "dumpstyle" so this can be changed easily in one place?
        [ii,molj,typej,q,x1,x2,x3,n1,n2,n3] = string.split(line)
        k=int(ii)
        typea[k] = int(typej)
        if i==0:
            if typea[k] == 1:
                num1=num1+1
            elif typea[k] == 2:
                num2=num2+1
            elif typea[k] == 3:
                num3=num3+1
            elif typea[k] == 4:
                num4=num4+1
        mol[k] = int(molj)
        xc[k] = xbox*(float(x1)-0.5) #scaled coords go from 0 to 1; want to go from -xbox/2 to xbox/2
        yc[k] = ybox*(float(x2)-0.5)
        zc[k] = zbox*(float(x3)-0.5)

# rdf stuff
#assuming a cubic box
dk = 2*pi/xbox
numk = 30 
numg1=numk+1

# Zero averages: 
kountSk = 0
nk = zeros(numg1,int)
kcc = zeros(numg1,float)
S11 = zeros(numg1,float32)  
S12 = zeros(numg1,float32)
S22 = zeros(numg1,float32)
S13 = zeros(numg1,float32)
S23 = zeros(numg1,float32)
S33 = zeros(numg1,float32)
S14 = zeros(numg1,float32)
S24 = zeros(numg1,float32)
S34 = zeros(numg1,float32)
S44 = zeros(numg1,float32)
SS11 = zeros(numg1,float32)  
SS12 = zeros(numg1,float32)
SS22 = zeros(numg1,float32)
SS13 = zeros(numg1,float32)
SS23 = zeros(numg1,float32)
SS33 = zeros(numg1,float32)
SS14 = zeros(numg1,float32)
SS24 = zeros(numg1,float32)
SS34 = zeros(numg1,float32)
SS44 = zeros(numg1,float32)

#call g(r) here for last read config (=initial config if iskip=0)
Sk()
kountSk = 1
for i in range(1,numg1):
    kcc[i] = kcc[i]/nk[i]
    SS11[i] = S11[i] / (nk[i] * (num1))
    SS12[i] = S12[i] / (nk[i] * (num1 + num2)/2)
    SS13[i] = S13[i] / (nk[i] * (num1 + num3)/2)
    SS14[i] = S14[i] / (nk[i] * (num1 + num4)/2)
    SS22[i] = S22[i] / (nk[i] * (num2))
    SS23[i] = S23[i] / (nk[i] * (num2 + num3)/2)
    SS24[i] = S24[i] / (nk[i] * (num2 + num4)/2)
    SS33[i] = S33[i] / (nk[i] * (num3))
    SS34[i] = S34[i] / (nk[i] * (num3 + num4)/2)
    SS44[i] = S44[i] / (nk[i] * (num4 + num4)/2)


#    Output
OUT = open(file, 'w')
OUT.write("#%7i\n" % (kountSk))
OUT.write("k S11 S12 S22 S13 S23 S33 S14 S24 S34 S44\n")
for ig in range(1,numg1):
    OUT.write("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n" % (kcc[ig],SS11[ig],SS12[ig],SS22[ig],SS13[ig],SS23[ig],SS33[ig],SS14[ig],SS24[ig],SS34[ig],SS44[ig]))
OUT.close()

print "Reading config file...."

istart = iskip + 1
# Read configuration from zconfig 
for kconf in range(istart,nconf+iskip): 
    sys.stdin.readline()
    sys.stdin.readline()      # time step
    sys.stdin.readline()
    line = sys.stdin.readline()      # number of atoms
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
        [ii,molj,typej,q,x1,x2,x3,n1,n2,n3] = string.split(line)
        k=int(ii)
        typea[k] = int(typej)
        mol[k] = int(molj)
        xc[k] = xbox*(float(x1)-0.5)
        yc[k] = ybox*(float(x2)-0.5)
        zc[k] = zbox*(float(x3)-0.5)

    # g(r)
    Sk()

    kountSk = kountSk + 1
    

    for i in range(1,numg1):
        kcc[i] = kcc[i]/nk[i]
        SS11[i] = S11[i] / (nk[i] * (num1))
        SS12[i] = S12[i] / (nk[i] * (num1 + num2)/2)
        SS13[i] = S13[i] / (nk[i] * (num1 + num3)/2)
        SS14[i] = S14[i] / (nk[i] * (num1 + num4)/2)
        SS22[i] = S22[i] / (nk[i] * (num2))
        SS23[i] = S23[i] / (nk[i] * (num2 + num3)/2)
        SS24[i] = S24[i] / (nk[i] * (num2 + num4)/2)
        SS33[i] = S33[i] / (nk[i] * (num3))
        SS34[i] = S34[i] / (nk[i] * (num3 + num4)/2)
        SS44[i] = S44[i] / (nk[i] * (num4 + num4)/2)


    #    Output
    OUT = open(file, 'w')
    OUT.write("#%7i\n" % (kountSk))
    OUT.write("k S11 S12 S22 S13 S23 S33 S14 S24 S34 S44\n")
    for ig in range(1,numg1):
        OUT.write("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n" % (kcc[ig],SS11[ig],SS12[ig],SS22[ig],SS13[ig],SS23[ig],SS33[ig],SS14[ig],SS24[ig],SS34[ig],SS44[ig]))
    OUT.close()

# end of loop

