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

nconf = 500
iskip = -1
step=1

xlength = arange(0,301,step)
ylength = arange(0,301,step)
Y, Z = meshgrid(xlength, ylength)
X = zeros_like(Y)

file = 'entanglement_grid_yz.txt'
def gridfunc():
	#global xlength, ylength, dim, Z
	#find location of tag particle
	for k in range(1,dim):
		if typea[k] == 1 and k%2 != 0:
			yc[k] = yc[k] + 25
			zc[k] = zc[k] + 25
			ygrid = int(yc[k]/step)
			zgrid = int(zc[k]/step)
			X[ygrid][zgrid] += 1

for i in range(0,iskip+1):
	sys.stdin.readline()
	sys.stdin.readline() #timestep
	sys.stdin.readline()
	line = sys.stdin.readline() # number of atoms
	fields = string.split(line) #reads the next line
	natoms = int(fields[0])  #reads number of atoms from first configuration; don't support changing num atoms
	dim=natoms+1
	sys.stdin.readline()
	sys.stdin.readline() #xbox bounds
	sys.stdin.readline() #ybox bounds
	sys.stdin.readline() #zbox bounds
	sys.stdin.readline()
   	for j in range(1,dim):
           line = sys.stdin.readline()

	#conf = conf + 1
	#print conf

print "reading config file..."
istart = iskip+1

OUT = open(file, 'w')
OUT.write("Lateral X   Lateral Y     Height\n")

for i in range(istart,nconf+iskip+1):
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
	sys.stdin.readline()
	line = sys.stdin.readline() #xbox bounds  
	[xm,xp] = map(float,line.split())
	line = sys.stdin.readline() #ybox bounds
	[ym,yp] = map(float,line.split())
	line = sys.stdin.readline() #zbox bounds
	[zm,zp] = map(float,line.split())
	line = sys.stdin.readline()
	Lx = xp - xm
	Ly = yp - ym
	Lz = zp - zm
	vol = Lx*Ly*Lz
	Lx2 = Lx/2.0
	Ly2 = Ly/2.0
	Lz2 = Lz/2.0
	for j in range(1,dim):
	    line = sys.stdin.readline()
	    [ii,molj,typej,x1,x2,x3] = string.split(line) #atom number, molecule number, bead type, x,y,z, image flags
	    k=int(ii) #sets the index k equal to the atom number
	    typea[k]=int(typej)
	    mol[k] = int(molj) #sets the chain number for each bead equal to molecule number
	    xc[k] = float(x1)+Lx2
	    yc[k] = float(x2)+Ly2
	    zc[k] = float(x3)+Lz2
	        #should have assigned each bead a chain number
	        
	xcom = 0
	ycom = 0
	count = 0
	for k in range(1,dim):
		if mol[k] >= 2256 and typea[k] == 2: #used molecule number 2257 (particle 1 of 12) arbitrarily
			xcom += xc[k]
			ycom += yc[k]
			count += 1
	
	xcom = xcom/count
	ycom = ycom/count
	print xcom,ycom
	for k in range(1,dim):
		xc[k]=xc[k]-xcom
		yc[k]=yc[k]-ycom
		
	for k in range(1,dim):
		if (xc[k] > Lx):
			cx[k] = int(xc[k]/Lx)
			xc[k] = xc[k] - cx[k]*Lx
		elif (xc[k] < 0.0):
			cx[k] = -int((-xc[k]+Lx)/Lx)
			xc[k] = xc[k] - cx[k]*Lx
		else:
			cx[k] = 0
			xc[k] = xc[k]
		if (yc[k] > Ly):
			cy[k] = int(yc[k]/Ly)
			yc[k] = yc[k] - cy[k]*Ly
		elif (yc[k] < 0.0):
			cy[k] = -int((-yc[k]+Ly)/Ly)
			yc[k] = yc[k] - cy[k]*Ly
		else:
			cy[k] = 0
			yc[k] = yc[k]
		if (zc[k] > Lz):
			cz[k] = int(zc[k]/Lz)
			zc[k] = zc[k] - cz[k]*Lz 
		elif (zc[k] < 0.0):
			cz[k] = -int((-zc[k]+Lz)/Lz)
			zc[k] = zc[k] - cz[k]*Lz
		else:
			cz[k] = 0
			zc[k] = zc[k]

	gridfunc()
	print "Configuration complete"
	#conf = conf+1

for i,k in enumerate(xlength):
	for j,l in enumerate(ylength):
		#Z[j][k] = Z[j][k]/float(nconf)
		#print Z[j][k]
		OUT.writelines("%i %i %i\n" % (Y[i][j],Z[i][j],X[i][j]))


