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
 
# Script:  denprof_log_mod.py
# Purpose: create a density profile of block copolymers along z direction
# Syntax:  denprof_log_mod.py  
# Author:  Y. Seo
# -------------------------------------------------------

import sys,string
from numpy import *
from math import * 
import fileinput

# INPUT PARAMETERs
nconf =    77
npoly =    768
mg    =    160
nlayers =   4
nmon = 100
velcall = 0

d1 = zeros(mg+1,float32)
d2 = zeros(mg+1,float32)
d3 = zeros(mg+1,float32)
ddnum1 = zeros(mg+1,float32)
ddnum2 = zeros(mg+1,float32)
ddnum3 = zeros(mg+1,float32)

def dennum():
    for ii in range(1,natoms+1):
        zr = zc[ii]
        if zr >= zcm1:
            ig = int((zr - zcm1) / drg) + 1
            if ig < len(d1):
                if typea[ii] == 1:
                    dnum1[ig] += 1
                elif typea[ii] == 2:
                    dnum2[ig] += 1
                elif typea[ii] == 3 or 4:
                    dnum3[ig] += 1
        elif zr < zcm1:
            ig = int((zr+zbox - zcm1) / drg) + 1
            if ig < len(d1):
                if typea[ii] == 1:
                    dnum1[ig] += 1
                elif typea[ii] == 2:
                    dnum2[ig] += 1
                elif typea[ii] == 3 or 4:
                    dnum3[ig] += 1

def wholecalculation():
    global timestep, vel, xx, yy, zz, x0, y0, z0, x0com, y0com, z0com, xcom, ycom, zcom, zcm1, step, num1, num2, num3, num4, npoly, zc, drg, d1, d2, d3, zbox, dnum1, dnum2, dnum3
    global skip
    global natoms, typea, velcall
    
    infiles = ['loglog.lammpstrj']
    file = 'denprof_log.csv' #% skip
    dummy=0

    timestep=zeros(nconf,int)

    IN = fileinput.input(infiles)
    natoms = 0

    for loopnum in range(0,skip): 
                               
        IN.readline()
        line = IN.readline()      # time step
        #print "skip timestep %s" % line
        IN.readline()
        line = IN.readline()      # number of atoms
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

    for loopnum in range(0,1): #read starting timestep
                               
        IN.readline()
        line = IN.readline()      # time step
        fields = string.split(line)
        step = int(fields[0])
        #print "read step %i" % step
        IN.readline()
        line = IN.readline()      # number of atoms
        fields = string.split(line)

        natoms = int(fields[0])  #reads number of atoms from first configuration; don't support changing num atoms
        dim=natoms+1
        xc=zeros(dim,float32)
        yc=zeros(dim,float32)
        zc=zeros(dim,float32)
        zc1=zeros(dim,float32)
        theta=zeros(dim,float32)
        cirx=zeros(dim,float32)
        ciry=zeros(dim,float32)
        xx=zeros(dim,float32)
        yy=zeros(dim,float32)
        zz=zeros(dim,float32)
        x0=zeros(dim,float32)
        y0=zeros(dim,float32)
        z0=zeros(dim,float32)
        x0com=zeros(npoly+1,float32)
        y0com=zeros(npoly+1,float32)
        z0com=zeros(npoly+1,float32)
        nbead=zeros(npoly+1)
        xi=zeros(dim)
        yi=zeros(dim)
        zi=zeros(dim)
        xr=zeros(dim)
        yr=zeros(dim)
        zr=zeros(dim)
        r2=zeros(dim)
        zscaled=zeros(dim)
        typea=[0]*dim
        mol=[0]*dim
            #
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
        zboxi = zbox/mg
        drg = zbox/mg
        vol = xbox*ybox*zbox
        xbox2 = xbox/2.0
        xcom=zeros(npoly+1,float32) #initialize these before every config read
        ycom=zeros(npoly+1,float32)
        zcom=zeros(npoly+1,float32)
        xcm=ycm=zcm=0.0
        num1 = 0
        num2 = 0
        num3 = 0
        cirxcm = 0
        cirycm = 0
        for j in range(1,dim):
            line = IN.readline()
            #todo: make a list called "dumpstyle" so this can be changed easily in one place?

            istep=step

        #print "Reading config file...."

        istart = loopnum + 1
        # Read configuration from zconfig     
        for kconf in range(istart,2):
          #try:
            IN.readline()
            line = IN.readline()      # time step
            fields = string.split(line)
            step = int(fields[0])
            IN.readline()
            line = IN.readline()      # number of atoms
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
            drg = zbox/mg
            vol = xbox*ybox*zbox
            xbox2 = xbox/2.0
            zboxi = zbox/nlayers
            xcom=zeros(npoly+1,float32) #initialize these before every config read
            ycom=zeros(npoly+1,float32)
            zcom=zeros(npoly+1,float32)
            xcm=ycm=zcm=0.0
            nbead=zeros(npoly+1)
            for j in range(1,dim):
                line = IN.readline()
                #todo: make a list called "dumpstyle" so this can be changed easily in one place?
                [ii,molj,typej,q,x1,x2,x3,n1,n2,n3] = string.split(line)
                k=int(ii)
                typea[k] = int(typej)
                mol[k] = int(molj)
                xc[k] = xbox*float(x1) #scaled coords go from 0 to 1; want to go from 0 to xbox
                yc[k] = ybox*float(x2)
                zc[k] = zbox*float(x3)
                for ii in range(0,nlayers):      #combine the four layers of lamellae
                    if zboxi*ii < zc[k] <= zboxi*(ii+1):
                        zc1[k] = zc[k] - zboxi*ii
                if typea[k] == 1 or typea[k] == 2:
                    theta[k] = (2*pi/zboxi)*zc1[k] #circling z-axis for one layer (combined with others)
                    cirx[k] = (zboxi/(2*pi))*sin(theta[k]) #x coordinate in cylindrical coordinates(xytheta), R=zbox/(2*pi)
                    ciry[k] = (zboxi/(2*pi))*cos(theta[k]) #y coordinate in cylindrical coordinates(xytheta), R=zbox/(2*pi)
        
                if typea[k] == 1:
                    cirxcm += cirx[k]/(npoly*nmon*0.5)
                    cirycm += ciry[k]/(npoly*nmon*0.5)
                    thetacm = atan2(-cirxcm,-cirycm) + pi #atan2(X,Y); return theta in positive way from (0,R), arctan(X/Y), range from -pi to pi, atan2(1,1)=pi/4 and atan2(-1,-1)=-3*pi/4; +pi to make the range from 0 to 2pi
                    zcm1 = (zboxi/(2*pi))*thetacm
                xi[k] = int(round(float(n1)))
                yi[k] = int(round(float(n2)))
                zi[k] = int(round(float(n3)))

            #print "%7i %7i %8.4f %8.4f %8.4f %8.4f\n" % (istep,step,vel[1],vel[2],vel[3],vel[4])
            #print step,istep,step-istep,velcall
            dnum1 = [0]*(mg+1)
            dnum2 = [0]*(mg+1)
            dnum3 = [0]*(mg+1)
            dennum()
            for ig in range(1,mg+1):
                ddnum1[ig] += float(dnum1[ig])/(vol/mg)/nconf
                ddnum2[ig] += float(dnum2[ig])/(vol/mg)/nconf
                ddnum3[ig] += float(dnum3[ig])/(vol/mg)/nconf
                    
          #except: break
        fileinput.close()

    # end of loop
        OUT = open(file, 'w')
        OUT.write("z, den(A), den(B), den(ion)\n")
        for i in range(1,mg+1):
            OUT.write("%7i, %8.9f, %8.9f, %8.9f\n" % (i,ddnum1[i],ddnum2[i],ddnum3[i]))
        OUT.close()

#main:
#skip=0
#wholecalculation()
for xx in range(0,nconf):
#run with skip = 0, 11, 22
    print xx
    skip=xx*20
    #skipnum=skip/200
    wholecalculation()

    
