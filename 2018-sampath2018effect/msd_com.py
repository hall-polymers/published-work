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
 
# Script:  msdions.py
# Purpose: calculate mean sq displacement of various types
# Syntax:  msdions.py < filename 
# Example: msdions.py < test.dump (dump file with scaled coordinates)
# Author:  Lisa Hall from Mark's g(r) code
 
# derived from fortran code
# -------------------------------------------------------

import sys,string
from numpy import *
from math import * 
import fileinput

# INPUT PARAMETERs
nconf =  62 #20030 # number of configurations to take data for
#stepspacing = 10000000
#polylength= 30  #now based on type and not polymer number
npoly =    800  #need this for com calc

def msdcalc():
    global msdvec1, msdvec2, msdvec3, msdvec4, msdvec5, msdvec6, msdvec7, timestep, msd, xx, yy, zz, x0, y0, z0, x0com, y0com, z0com, xcom, ycom, zcom, step, num1, num2, num3, num4, num5, num6, npoly
    global natoms, typea,msdcall
    # needs xc,yc,zc,xbox,ybox,zbox,drg
    # outputs g[]
    #from numpy import * 
    # debug
##    r2mol=zeros(numatoms,float32) #todo:change to num molecules
##    ninmol=zeros(numatoms, int)
    msd[1]=0
    msd[2]=0
    msd[3]=0
    msd[4]=0
    msd[5]=0
    msd[6]=0
    msd[7]=0
    for ii in range(1,natoms+1): 
        xr = xx[ii] - x0[ii]
        yr = yy[ii] - y0[ii]
        zr = zz[ii] - z0[ii]
        r2 = xr*xr + yr*yr + zr*zr
        msd[typea[ii]] = msd[typea[ii]]+r2
    for ii in range(1,npoly+1):
        xr = xcom[ii] - x0com[ii]
        yr = ycom[ii] - y0com[ii]
        zr = zcom[ii] - z0com[ii]
        r2 = xr*xr + yr*yr + zr*zr
        msd[7] = msd[7]+r2
    msd[1]=msd[1]/num1
    msd[2]=msd[2]/num2
    msd[3]=msd[3]/num3
    msd[4]=msd[4]/num4
    msd[5]=msd[5]/num5
    msd[6]=msd[6]/num6
    msd[7]=msd[7]/npoly
    msdvec1[msdcall]=msdvec1[msdcall]+msd[1]
    msdvec2[msdcall]=msdvec2[msdcall]+msd[2]
    msdvec3[msdcall]=msdvec3[msdcall]+msd[3]
    msdvec4[msdcall]=msdvec4[msdcall]+msd[4]
    msdvec5[msdcall]=msdvec5[msdcall]+msd[5]
    msdvec6[msdcall]=msdvec6[msdcall]+msd[6]
    msdvec7[msdcall]=msdvec7[msdcall]+msd[7]
    timestep[msdcall]=step

    #OUT = open(file, 'a')
   # OUT.write("%7i %8.4f %8.4f %8.4f\n" % (step,msd[1],msd[2],msd[3]))
   # for ii in range(1,natoms+1):
    #    OUT.write("%7i %7i %7i %7i %7i %8.4f %8.4f %8.4f\n" % (mol[ii],typea[ii],xi[ii],yi[ii],zi[ii],xc[ii],yc[ii],zc[ii]))
   # OUT.close()
# end def msdcalc

def wholecalculation():
    global msdvec1, msdvec2, msdvec3, msdvec4, msdvec5, msdvec6,msdvec7, timestep, msd, xx, yy, zz, x0, y0, z0, x0com, y0com, z0com, xcom, ycom, zcom, step, num1, num2, num3, num4, num5, num6, npoly
    global skip
    global natoms, typea,msdcall
    
    infiles = ['equilnpt1.dump']
    file = 'msdcm_R1skip%d.csv'% skip
    dummy=0

    msdcall=0
    msdvec1=zeros(nconf,float)
    msdvec2=zeros(nconf,float)
    msdvec3=zeros(nconf,float)
    msdvec4=zeros(nconf,float)
    msdvec5=zeros(nconf,float)
    msdvec6=zeros(nconf,float)
    msdvec7=zeros(nconf,float)
    timestep=zeros(nconf,int)

    IN = fileinput.input(infiles)
    natoms = 0

    for loopnum in range(0,skip):          
        IN.readline()
        line = IN.readline()      # time step
        print "skip timestep %s" % line
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
        print "read step %i" % step
        IN.readline()
        line = IN.readline()      # number of atoms
        fields = string.split(line)
        natoms = int(fields[0])  #reads number of atoms from first configuration; don't support changing num atoms
        dim=natoms+1
        xc=zeros(dim,float32)
        yc=zeros(dim,float32)
        zc=zeros(dim,float32)
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
        typea=[0]*dim
        mol=[0]*dim
        msd=zeros(8,float32)
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
        xcom=zeros(npoly+1,float32) #initialize these before every config read
        ycom=zeros(npoly+1,float32)
        zcom=zeros(npoly+1,float32)
        xcm=ycm=zcm=0.0
        num1 = 0
        num2 = 0
        num3 = 0
        num4 = 0
        num5 = 0
        num6 = 0
        for j in range(1,dim):
            line = IN.readline()
            #todo: make a list called "dumpstyle" so this can be changed easily in one place?
            [ii,molj,typej,q,x1,x2,x3,n1,n2,n3] = string.split(line)
            k=int(ii)
            typea[k] = int(typej)
            if typea[k] == 1:
                num1=num1+1
            elif typea[k] == 2:
                num2=num2+1
            elif typea[k] == 3:
                num3=num3+1
            elif typea[k] == 4:
                num4=num4+1
            elif typea[k] == 5:
                num5=num5+1
            elif typea[k] == 6:
                num6=num6+1
            mol[k] = int(molj)
            xc[k] = xbox*(float(x1)-0.5) #scaled coords go from 0 to 1; want to go from -xbox/2 to xbox/2
            yc[k] = ybox*(float(x2)-0.5)
            zc[k] = zbox*(float(x3)-0.5)
            xi[k] = int(round(float(n1)))
            yi[k] = int(round(float(n2)))
            zi[k] = int(round(float(n3)))
            xx[k] = xc[k] + xbox*xi[k]
            yy[k] = yc[k] + ybox*yi[k]
            zz[k] = zc[k] + zbox*zi[k]
        for j in range(1,dim):
                xcm = xcm + xx[j]
                ycm = ycm + yy[j]
                zcm = zcm + zz[j]
        xcm=xcm/natoms
        ycm=ycm/natoms
        zcm=zcm/natoms
        for j in range(1,dim):
                xx[j] = xx[j] - xcm
                yy[j] = yy[j] - ycm
                zz[j] = zz[j] - zcm
        for j in range(1,dim):
            #save initial config
            x0[j] = xx[j]
            y0[j] = yy[j]
            z0[j] = zz[j]
            if mol[j] < npoly+1: #implies that polymers come first in molecule number
                xcom[mol[j]] = xcom[mol[j]]+xx[j]
                ycom[mol[j]] = ycom[mol[j]]+yy[j]
                zcom[mol[j]] = zcom[mol[j]]+zz[j]
                nbead[mol[j]] = nbead[mol[j]]+1
        for m in range(1,npoly+1):
            xcom[m] = xcom[m]/float(nbead[m])
            ycom[m] = ycom[m]/float(nbead[m])
            zcom[m] = zcom[m]/float(nbead[m])
            #save initial config
            x0com[m] = xcom[m]
            y0com[m] = ycom[m]
            z0com[m] = zcom[m]

            istep=step

        msdcalc()
        msdcall+=1

        print "Reading config file...."

        istart = loopnum + 1
        # Read configuration from zconfig     
        for kconf in range(istart,nconf):
          try:
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
            vol = xbox*ybox*zbox
            xbox2 = xbox/2.0
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
                xc[k] = xbox*(float(x1)-0.5) #scaled coords go from 0 to 1; want to go from -xbox/2 to xbox/2
                yc[k] = ybox*(float(x2)-0.5)
                zc[k] = zbox*(float(x3)-0.5)
                xi[k] = int(round(float(n1)))
                yi[k] = int(round(float(n2)))
                zi[k] = int(round(float(n3)))
                xx[k] = xc[k] + xbox*xi[k]
                yy[k] = yc[k] + ybox*yi[k]
                zz[k] = zc[k] + zbox*zi[k]
            for j in range(1,dim):
                xcm = xcm + xx[j]
                ycm = ycm + yy[j]
                zcm = zcm + zz[j]
            xcm=xcm/natoms
            ycm=ycm/natoms
            zcm=zcm/natoms
            for j in range(1,dim):
                xx[j] = xx[j] - xcm
                yy[j] = yy[j] - ycm
                zz[j] = zz[j] - zcm
            for j in range(1,dim):
                if mol[j] < npoly+1: #implies that polymers come first in molecule number
                    xcom[mol[j]] = xcom[mol[j]]+xx[j]
                    ycom[mol[j]] = ycom[mol[j]]+yy[j]
                    zcom[mol[j]] = zcom[mol[j]]+zz[j]
                    nbead[mol[j]] = nbead[mol[j]]+1
            for m in range(1,npoly+1):
                xcom[m] = xcom[m]/float(nbead[m])
                ycom[m] = ycom[m]/float(nbead[m])
                zcom[m] = zcom[m]/float(nbead[m])


            # calculate msd for this config
            msdcalc()
            msdcall+=1
            #    Output
            #OUT = open(file, 'a')
            #OUT.write("%7i %8.4f %8.4f %8.4f %8.4f %8.4f\n" % (step,msd[1],msd[2],msd[3],msd[4],msd[5]))
            #OUT.close()
            #except: break
            print "%7i %7i %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n" % (istep,step,msd[1],msd[2],msd[3],msd[4],msd[5],msd[6],msd[7])
          except: break
        fileinput.close()

    # end of loop
    OUT = open(file, 'w')
    #OUT.write("msd\n")
    OUT.write("step, msd_type1, msd_type2, msd_type3, msd_type4, msd_type5, msd_type6, msd_com\n")
    for i in range(0,msdcall):
    ##    msdvec1[i]=msdvec1[i]/(nconf-i)
    ##    msdvec2[i]=msdvec2[i]/(nconf-i)
    ##    msdvec3[i]=msdvec3[i]/(nconf-i)
    ##    msdvec4[i]=msdvec4[i]/(nconf-i)
    ##    msdvec5[i]=msdvec5[i]/(nconf-i)
        OUT.write("%7i, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f\n" % (timestep[i],msdvec1[i],msdvec2[i],msdvec3[i],msdvec4[i],msdvec5[i],msdvec6[i],msdvec7[i]))
    OUT.close()

#main:
skip=0
wholecalculation()

    
