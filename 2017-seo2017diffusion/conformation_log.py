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
import fileinput

taper = 0.3
nbeads = 40
hntaper = int(taper*nbeads*0.5)
hnpure = int((1-taper)*nbeads*0.5)
print "taper, hntaper, hnpure:", taper,hntaper,hnpure
nconf =    77
npoly =    1920
nlayers =   4
velcall = 0

fold1 = [0]*nconf
fold2 = [0]*nconf
stretch1 = [0]*nconf
stretch2 = [0]*nconf

def e2e():
    num1 = num2 = num3 = num4 = num5 = num6 = [0]*npoly
    n1 = n2 = n3 = n4 = n5 = n6 = float(0)
    for jj in range(0,npoly):
        for k in range(jj*nbeads+1, jj*nbeads+nbeads+1):
            if k%nbeads == 0:
                n1 = xc[k] 
                n2 = yc[k]
                n3 = zc[k]
            elif k == 1 or (k-1)%nbeads == 0:
                n4 = xc[k]
                n5 = yc[k]
                n6 = zc[k]
        square_e2edistance[jj] = (((n1-n4)**2)+((n2-n5)**2)+((n3-n6)**2))
        square_e2ezdistance[jj] = ((n3-n6)**2)

def com():
    for jj in range(0,npoly):
        p1x = p1y = p1z = t1x = t1y = t1z = p2x = p2y = p2z = t2x = t2y = t2z = float(0)
        for k in range(jj*nbeads+1,jj*nbeads+1+hnpure):
            p1x += xc[k]/hnpure
            p1y += yc[k]/hnpure
            p1z += zc[k]/hnpure
        for k in range(jj*nbeads+1+hnpure,jj*nbeads+1+hnpure+hntaper):
            t1x += xc[k]/hntaper
            t1y += yc[k]/hntaper
            t1z += zc[k]/hntaper
        for k in range(jj*nbeads+1+hnpure+hntaper,jj*nbeads+1+hnpure+hntaper*2):
            t2x += xc[k]/hntaper
            t2y += yc[k]/hntaper
            t2z += zc[k]/hntaper
        for k in range(jj*nbeads+1+hnpure+hntaper*2,jj*nbeads+1+hnpure*2+hntaper*2):
            p2x += xc[k]/hnpure
            p2y += yc[k]/hnpure
            p2z += zc[k]/hnpure
        t1zp1z[jj] = t1z - p1z
        p1zt1z[jj] = p1z - t1z
        t1zt2z[jj] = t1z - t2z
        t2zp2z[jj] = t2z - p2z

def wholecalculation():
    global timestep, step, num1, num2, num3, num4, npoly, zc, drg, d1, d2, d3, v1t, v2t, v3t, zbox, nbeads, xcm1, ycm1, zcm1, xcm2, ycm2, zcm2, xc, yc, zc, k, ixx, iyy, izz, ixy, ixz, iyz, a, square_e2edistance, square_e2ezdistance, Lz, skip, p1x, p1y, p1z, t1x, t1y, t1z, t2x, t2y, t2z, p2x, p2y, p2z, t1zp1z, p1zt1z, t1zt2z, t2zp2z, hnpure, hntaper
    global skip
    global natoms, typea, velcall
    
    infiles = ['loglog.lammpstrj']
    file = 'composition_log.csv' #% skip
    dummy=0
    timestep=zeros(nconf,int)

    IN = fileinput.input(infiles)
    natoms = 0

    for loopnum in range(0,skip): 
                               
        IN.readline()
        line = IN.readline()     
        IN.readline()
        line = IN.readline()    
        fields = string.split(line)
        natoms = int(fields[0])
        dim=natoms+1
        IN.readline()
        line = IN.readline()      
        line = IN.readline()      
        line = IN.readline()      
        line = IN.readline()

        for j in range(1,dim):
            IN.readline()

    for loopnum in range(0,1): 
                               
        IN.readline()
        line = IN.readline() 
        fields = string.split(line)
        step = int(fields[0])
        IN.readline()
        line = IN.readline()
        fields = string.split(line)

        natoms = int(fields[0])
        dim=natoms+1
        xc=zeros(dim,float32)
        yc=zeros(dim,float32)
        zc=zeros(dim,float32)
        nbead=zeros(npoly+1)
        xcm=zeros(npoly+1)
        ycm=zeros(npoly+1)
        zcm=zeros(npoly+1)
        xi=zeros(dim)
        yi=zeros(dim)
        zi=zeros(dim)
        xr=zeros(dim)
        yr=zeros(dim)
        zr=zeros(dim)
        r2=zeros(dim)
        typea=[0]*dim
        mol=[0]*dim
        square_e2edistance = [0]*npoly
        square_e2ezdistance = [0]*npoly
        num1 = num2 = num3 = num4 = num5 = num6 = [0]*npoly
        zcm1=[0]*npoly
        xcm1=[0]*npoly
        ycm1=[0]*npoly
        zcm2=[0]*npoly
        xcm2=[0]*npoly
        ycm2=[0]*npoly
        t1zp1z=[0]*npoly
        p1zt1z=[0]*npoly
        t1zt2z=[0]*npoly
        t2zp2z=[0]*npoly
        
        IN.readline()
        line = IN.readline()      
        [xm,xp] = map(float,line.split())
        line = IN.readline()      
        [ym,yp] = map(float,line.split())
        line = IN.readline()      
        [zm,zp] = map(float,line.split())
        line = IN.readline()
        xbox = xp - xm
        ybox = yp - ym
        zbox = zp - zm
        vol = xbox*ybox*zbox
        xbox2 = xbox/2.0
        num1 = 0
        num2 = 0
        num3 = 0
        for j in range(1,dim):
            line = IN.readline()
            [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line)
            k=int(ii)
            typea[k] = int(typej)
            if typea[k] == 1:
                num1=num1+1
            elif typea[k] == 2:
                num2=num2+1
            elif typea[k] == 3:
                num3=num3+1
            mol[k] = int(molj)
            xc[k] = xbox*float(x1) + int(n1)*xbox
            yc[k] = ybox*float(x2) + int(n2)*ybox
            zc[k] = zbox*float(x3) + int(n3)*zbox

            istep=step

        istart = loopnum + 1
        for kconf in range(istart,2):
            IN.readline()
            line = IN.readline()
            fields = string.split(line)
            step = int(fields[0])
            IN.readline()
            line = IN.readline()
            fields = string.split(line)
            num = fields[0]
            IN.readline()
            line = IN.readline()      
            [xm,xp] = map(float,line.split())
            line = IN.readline()      
            [ym,yp] = map(float,line.split())
            line = IN.readline()      
            [zm,zp] = map(float,line.split())
            line = IN.readline()
            xbox = xp - xm
            ybox = yp - ym
            zbox = zp - zm
            vol = xbox*ybox*zbox
            Lz = zbox/nlayers
            nbead=zeros(npoly+1)
            for j in range(1,dim):
                line = IN.readline()
                [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line)
                k=int(ii)
                typea[k] = int(typej)
                mol[k] = int(molj)
                xc[k] = xbox*float(x1) + int(n1)*xbox
                yc[k] = ybox*float(x2) + int(n2)*ybox
                zc[k] = zbox*float(x3) + int(n3)*zbox

            velcall+=1
            print velcall
            if (step - istep) == 2:
                e2e()
                com()
                for jj in range(0,npoly):
                    if sqrt(square_e2ezdistance[jj]) > Lz:
                        stretch2[xx] += 1
                    elif t1zp1z[jj]*t2zp2z[jj] > 0:
                        fold1[xx] += 1
                    elif p1zt1z[jj]*t1zt2z[jj] > 0:
                        stretch1[xx] += 1
                    else:
                        fold2[xx] += 1
                #print xx, stretch1[xx], fold2[xx], fold1[xx], stretch2[xx]
                    
        fileinput.close()

    avgstretch1 = 0
    avgfold2 = 0
    avgfold1 = 0
    avgstretch2 = 0
    stdstretch1 = 0
    stdfold2 = 0
    stdfold1 = 0
    stdstretch2 = 0
    if xx == (nconf-1):
        OUT = open(file, 'w')
        OUT.write("nconf, stretch1, fold2, fold1, stretch2\n")
        for j in range(0,nconf):
            avgstretch1 += float(stretch1[j])/nconf
            avgfold2 += float(fold2[j])/nconf
            avgfold1 += float(fold1[j])/nconf
            avgstretch2 += float(stretch2[j])/nconf
            OUT.write("%7i, %7i, %7i, %7i, %7i\n" % (j,stretch1[j],fold2[j],fold1[j],stretch2[j]))
        stdstretch1 = std(stretch1)
        stdfold2 = std(fold2)
        stdfold1 = std(fold1)
        stdstretch2 = std(stretch2)
        OUT.write("Average and Stdev\n")
        OUT.write("%7i, %8.2f, %8.2f, %8.2f, %8.2f, %8.4f, %8.4f, %8.4f, %8.4f\n" % (j,avgstretch1,avgfold2,avgfold1,avgstretch2,stdstretch1,stdfold2,stdfold1,stdstretch2))
        OUT.close()

for xx in range(0,nconf):
    skip=xx*20
    wholecalculation()

    
