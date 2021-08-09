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
import fileinput

nconf = 1 
npoly = 480
nbeads = 100
nbonds = nbeads - 1

def wholecalculation():
    global iskip, natoms, typea, nbonds, nbeads, npoly, nconf, istart
    
    infiles = ['equil.lammpstrj']
    file = 'cmsid_%d.txt'% iskip

    IN = fileinput.input(infiles)
    natoms = 0  
      
    for i in range(0,iskip):
        IN.readline()
        IN.readline()
        IN.readline()
        line = IN.readline() 
        fields = string.split(line) 
        natoms = int(fields[0]) 
        dim=natoms+1
        xc=zeros(dim,float32)
        yc=zeros(dim,float32)
        zc=zeros(dim,float32)
        zci=zeros(dim,float32)
        cx=zeros(dim)
        cy=zeros(dim)
        cz=zeros(dim)
        typea=[0]*dim
        mol=[0]*dim
        IN.readline()
        line = IN.readline()  
        [xm,xp] = map(float,line.split())
        line = IN.readline() 
        [ym,yp] = map(float,line.split())
        line = IN.readline() 
        [zm,zp] = map(float,line.split())
        line = IN.readline()
        for j in range(1,dim):
            line = IN.readline()

    print "reading config file..."
    istart = iskip+1

    for i in range(istart,nconf+iskip+1):
        IN.readline()
        IN.readline()
        IN.readline()
        line = IN.readline() 
        fields = string.split(line) 
        natoms = int(fields[0]) 
        dim=natoms+1
        xc=zeros(dim,float32)
        yc=zeros(dim,float32)
        zc=zeros(dim,float32)
        zci=zeros(dim,float32)
        cx=zeros(dim)
        cy=zeros(dim)
        cz=zeros(dim)
        typea=[0]*dim
        mol=[0]*dim
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
        ybox2 = ybox/2.0
        zbox2 = zbox/2.0
        msidistp = [0]*(nbonds+1)   
        for j in range(1,dim):
            line = IN.readline()
            [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line) 
            k=int(ii) 
            mol[k] = molj
            typea[k] = typej
            xc[k] = xbox*float(x1) + int(n1)*xbox
            yc[k] = ybox*float(x2) + int(n2)*ybox
            zc[k] = zbox*float(x3) + int(n3)*zbox

        for n in range(0,npoly):
            for m in range(1,nbeads+1):
                if m < nbeads/2:
                    msidistp[nbeads/2-m] += (xc[n*nbeads+m]-xc[n*nbeads+nbeads/2])**2 + (yc[n*nbeads+m]-yc[n*nbeads+nbeads/2])**2 + (zc[n*nbeads+m]-zc[n*nbeads+nbeads/2])**2
                elif m > nbeads/2+1:
                    msidistp[m-nbeads/2-1] += (xc[n*nbeads+m]-xc[n*nbeads+nbeads/2+1])**2 + (yc[n*nbeads+m]-yc[n*nbeads+nbeads/2+1])**2 + (zc[n*nbeads+m]-zc[n*nbeads+nbeads/2+1])**2
    
        fileinput.close()


    OUT = open(file, 'w')
    OUT.write("#%7i\n" % (iskip))
    OUT.write("nbonds    <Rstoc^2>/n\n")
    msidist = [0]*(nbonds+1)
    for i in range(1, nbonds+1):
        msidist[i] += msidistp[i]/2/npoly
    for k in range(1,nbonds+1):
        OUT.write("%8.4f  %8.4f\n" % (k, msidist[k]/k))
    print "completed!_%d"% iskip
    OUT.close()

for xx in range(1,5):
    iskip = xx*500
    if xx == 4:
        iskip = 1990
    wholecalculation()
