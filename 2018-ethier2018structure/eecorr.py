#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

""" End to end vector correlation 
    Jeff Ethier
    November 10, 2015
    python eecorr.py < (filename).dump
    This program will compute the end to end correlation function."""

import sys,string
from numpy import *
from math import *


tconf = 504          #number of configurations in dump file (including 0)
npoly = 188
nbeads = 160   
ndata = 293         #configurations in each block
NP = 12
#Initial conditions
natoms = npoly*NP*nbeads
mdata = [[0]*1 for i in range(1,ndata+1)]


file = 'endtoendACF_10_check.txt'
OUT = open(file, 'w')

print "Reading config file..."
iconf = 0
iskip = 0
col = 0
nave = 10
nstart = 21        #reset starting configuration every nstart configurations
nstop = 211        #stop block averaging at nstop configuration
for s in range(0,tconf):
    iconf = s + nstart + 1
    sys.stdin.seek(0)
    iskip = iconf - 1
    if iconf % nstart == 1 and iconf < nstop+2:
        for i in range(0,iskip):
            print i
            sys.stdin.readline()
            sys.stdin.readline() #timestep
            sys.stdin.readline()
            line = sys.stdin.readline() # number of atoms
            fields = string.split(line) #reads the next line
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
            zcm=[0]*npoly
            xcm=[0]*npoly
            ycm=[0]*npoly
            rg=[0]*npoly
            sys.stdin.readline()
            line = sys.stdin.readline() #xbox bounds  
            [xm,xp] = map(float,line.split())
            line = sys.stdin.readline() #ybox bounds
            [ym,yp] = map(float,line.split())
            line = sys.stdin.readline() #zbox bounds
            [zm,zp] = map(float,line.split())
            line = sys.stdin.readline()
            for j in range(1,dim):
                line = sys.stdin.readline()
               # [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line) #atom number, molecule number, bead type, x,y,z, image flags
        for i in range(iskip,iconf):
            #print i
            sys.stdin.readline()
            line = sys.stdin.readline() #timestep
            fields = string.split(line)
            timestep = int(fields[0])
            sys.stdin.readline()
            line = sys.stdin.readline() # number of atoms
            fields = string.split(line) #reads the next line
            natoms = int(fields[0])  #reads number of atoms from first configuration; don't support changing num atoms
            dim=natoms+1
            xc=zeros(dim,float32)
            yc=zeros(dim,float32)
            zc=zeros(dim,float32)
            cx=zeros(dim,float32)
            cy=zeros(dim,float32)
            cz=zeros(dim,float32)
            typea=[0]*dim
            mol=[0]*dim
            sys.stdin.readline()
            line = sys.stdin.readline() #xbox bounds  
            [xm,xp] = map(float,line.split())
            line = sys.stdin.readline() #ybox bounds
            [ym,yp] = map(float,line.split())
            line = sys.stdin.readline() #zbox bounds
            [zm,zp] = map(float,line.split())
            line = sys.stdin.readline()
            xbox = xp - xm
            ybox = yp - ym
            zbox = zp - zm
            vol = xbox*ybox*zbox
            xbox2 = xbox/2.0
            ybox2 = ybox/2.0
            zbox2 = zbox/2.0
            Re0 = [0]*(npoly*NP)
            Re1 = [0]*(npoly*NP)
            R_EE = [[0 for m in range(0,3)] for i in range(0,(npoly*NP))]
            for j in range(1,dim):
                line = sys.stdin.readline()
                [ii,molj,typej,x1,x2,x3] = string.split(line) #atom number, molecule number, bead type, x,y,z, image flags
                k = int(ii)
                mol[k] = int(molj)
                if mol[k] >= 380:
                    mol[k] -= 1
                xc[k] = float(x1) - xm
                yc[k] = float(x2) - ym
                zc[k] = float(x3) - zm
                typea[k] = int(typej)
            count =0
            for i in range(0,NP):
                for j in range(0,npoly):
                    for k in range(1,nbeads+1):
                        l = (npoly*i+j)*nbeads + k + (i+1)
                        pol = npoly*i+j 
                        if (l-(i+1)) % nbeads == 1:
                            Re0[pol] = [xc[l],yc[l],zc[l]]
                        if (l-(i+1)) % nbeads == 0:
                            count += 1
                            Re1[pol] = [xc[l],yc[l],zc[l]]
        print count
        for l in range(0,(npoly*NP)):
            for n in range(0,3):
                    R_EE[l][n] = Re1[l][n]-Re0[l][n]
        RE0 = 0
        for l in range(0,npoly*NP):
            for n in range(0,3):
                RE0 += ((R_EE[l][n])**2)/(npoly*12)
        print ("R^2(%s) = %s" % (timestep,RE0))

        mylist = []
        nconf = ((nave-1)*nstart+1)
        if nconf == 1:
            nconf == 0
        print nconf
        for i in range(iconf,(tconf+1)-nconf):
            sys.stdin.readline()
            line = sys.stdin.readline() #timestep
            fields = string.split(line)
            timestep = int(fields[0])
            print timestep
            ## if timestep % nevery == 0:
            sys.stdin.readline()
            line = sys.stdin.readline() # number of atoms
            fields = string.split(line) #reads the next line
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
            zcm=[0]*npoly
            xcm=[0]*npoly
            ycm=[0]*npoly
            rg=[0]*npoly
            sys.stdin.readline()
            line = sys.stdin.readline() #xbox bounds  
            [xm,xp] = map(float,line.split())
            line = sys.stdin.readline() #ybox bounds
            [ym,yp] = map(float,line.split())
            line = sys.stdin.readline() #zbox bounds
            [zm,zp] = map(float,line.split())
            line = sys.stdin.readline()
            Rcorr = 0
            Re0 = [0]*(npoly*NP)
            Re1 = [0]*(npoly*NP)
            R_E = [[0 for m in range(0,3)] for i in range(0,(npoly*NP))]
            for j in range(1,dim):
                line = sys.stdin.readline()
                [ii,molj,typej,x1,x2,x3] = string.split(line) #atom number, molecule number, bead type, x,y,z, image flags
                k = int(ii)
                mol[k] = int(molj)
                if mol[k] >= 380:
                    mol[k] -= 1
                xc[k] = float(x1) - xm
                yc[k] = float(x2) - ym
                zc[k] = float(x3) - zm
                typea[k] = int(typej)
            xc_new=[0]
            yc_new=[0]
            zc_new=[0]
            typea_new=[0]
            mol_new=[0]
            for i in range(0,NP):
                for j in range(0,npoly):
                    for k in range(1,nbeads+1):
                        l = (npoly*i+j)*nbeads + k + (i+1)
                        pol = npoly*i+j
                        if (l-(i+1)) % nbeads == 1:
                            Re0[pol] = [xc[l],yc[l],zc[l]]
                        if (l-(i+1)) % nbeads == 0:
                            Re1[pol] = [xc[l],yc[l],zc[l]]
                        
            for l in range(0,(npoly*NP)):
                for n in range(0,3):
                    R_E[l][n] = ((Re1[l][n]-Re0[l][n]))
            for l in range(0,(npoly*NP)):
                for n in range(0,3):
                    Rcorr += (R_E[l][n]*R_EE[l][n])/(npoly*12)/RE0
            mylist.append(Rcorr)
            print Rcorr
        nave = nave-1
        print nave
        length = len(mylist)
        if length < ndata:
            for j in range(length,ndata):
                mylist.append(" ")
        for j in range(1,ndata+1):
            mdata[j-1][0] = j
        for i in range(1,ndata+1):
            mdata[i-1].append(mylist[i-1])
        for k in range(1,ndata+1):
            OUT.writelines("%s " % item for item in mdata[k-1])
            OUT.writelines('\n')
        
            
    else:
        sys.stdin.readline()
        sys.stdin.readline() #timestep
        sys.stdin.readline()
        line = sys.stdin.readline() # number of atoms
        fields = string.split(line) #reads the next line
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
        zcm=[0]*npoly
        xcm=[0]*npoly
        ycm=[0]*npoly
        rg=[0]*npoly
        sys.stdin.readline()
        line = sys.stdin.readline() #xbox bounds  
        [xm,xp] = map(float,line.split())
        line = sys.stdin.readline() #ybox bounds
        [ym,yp] = map(float,line.split())
        line = sys.stdin.readline() #zbox bounds
        [zm,zp] = map(float,line.split())
        line = sys.stdin.readline()
        for j in range(1,dim):
            line = sys.stdin.readline()
            #[ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line) #atom number, molecule number, bead type, x,y,z, image flags


OUT.close()

