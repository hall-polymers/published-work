#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

### Purpose: Calculates isotropic structure factor ###
### Syntax: python test_sk.py < equinvt.lammpstrj ###
### Author: Janani Sampath ###
### Date:  July 2015 ###
### derived from Lisa Hall's code, derived Mark Stevens' code ###

import sys,string
from numpy import *
import numpy as np
from math import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize_scalar 

# INPUT PARAMETERs
nconf =   1 # number of configurations to average over
iskip =     0 #number of configurations to skip before taking data
numk  =    25
numg1 =    numk+1


file  = 'S_k_1.csv'


# Zero averages: 

nk = zeros(numg1,int)
kcc = zeros(numg1,float)
S11 = zeros(numg1,float32)  
S44 = zeros(numg1,float32)
SS11 = zeros(numg1,float32)  
SS44 = zeros(numg1,float32)


#initialize
natoms = 0
num1 = 0
num4 = 0
sofk = []*numg1



for i in range(0,nconf): 
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
    dkx = 2*pi/xbox
    Lx=numk*dkx
    ybox = yp - ym
    dky = 2*pi/ybox
    Ly=numk*dky
    zbox = zp - zm
    dkz = 2*pi/zbox
    Lz=numk*dkz
    vol = xbox*ybox*zbox
    xbox2 = xbox/2.0
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,q,x1,x2,x3,n1,n2,n3] = string.split(line)
        k=int(ii)
        typea[k] = int(typej)
        if i==0:
            if typea[k] == 1:
                num1=num1+1
            elif typea[k] == 4:
                num4=num4+1
        mol[k] = int(molj)
        xc[k] = xbox*(float(x1)-0.5) 
        yc[k] = ybox*(float(x2)-0.5)
        zc[k] = zbox*(float(x3)-0.5)
        

    for ll in range( 0,numg1):
        for mm in range(0,numg1):
           for nn in range(0,numg1):
               s1c=s1s=s2c=s2s=s3c=s3s=s4c=s4s=0
               kx=ll*dkx
               ky=mm*dky
               kz=nn*dkz
               kk=sqrt(kx*kx+ky*ky+kz*kz)
               jj = int(kk/dkx)
               if jj < numg1 and jj > 0:
                   kcc[jj] = kcc[jj]+kk
                   nk[jj] = nk[jj]+1
                   for ii in range(1,dim):
                       kdotr=kx*xc[ii]+ky*yc[ii]+kz*zc[ii]
                       if typea[ii] == 4:
                           s4c = s4c + cos(kdotr)
                           s4s = s4s + sin(kdotr)
                       elif typea[ii] == 1:
                           s1c = s1c + cos(kdotr)
                           s1s = s1s + sin(kdotr)                                

                   S11[jj]+=(s1c*s1c+s1s*s1s)
                   S44[jj]+=(s4c*s4c+s4s*s4s)



    for k in range(1,numg1):
        kcc[k] = kcc[k]/nk[k]
        SS11[k] = S11[k] / (nk[k] * (num1))
        SS44[k] = S44[k] / (nk[k] * (num4))


    ## for ig in range(1,numg1):
    ##     line = [kcc[ig],SS44[ig]]
    ##     sofk.append(line)
        
    ## for i in range(len(sofk)):
    ##     if sofk[i] == max(sofk, key = lambda x:x[1]): # finds the maximum S(k), and then selects the points just before and after this maximum
    ##         x = [sofk[i-1][0],sofk[i][0],sofk[i+1][0]]
    ##         y = [sofk[i-1][1],sofk[i][1],sofk[i+1][1]]

    ##         def f(x, p1, p2, p3): # 2 deg polynomial fit
    ##             return p1*x*x + p2*x + p3

    ##         popt, pcov = curve_fit(f,x,y) 

    ##         fm  = lambda x: -f(x, *popt) # finds the peak for the fitted curve
    ##         r = minimize_scalar(fm)
    ##         maximum = r["x"], f(r["x"], *popt)

    ##         print maximum


            
    OUT = open(file, 'w')
    OUT.write("k, S11, S44\n")
    for ig in range(1,numg1):
        OUT.write("%8.4f, %8.4f, %8.4f\n" % (kcc[ig],SS11[ig],SS44[ig]))
    OUT.close()
