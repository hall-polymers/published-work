### Purpose: This Script will compute Hermans orientation factor (between e2e vector and direction of pull) for systems undergoing deformation ####
### Syntax: python orientation.py < (filename).dump > orientation.txt ###
### Author: Janani Sampath ###
### Date: Feb 2016 ###

#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

import sys, string
from numpy import *
from math import *
import numpy as np
import math as math


tframes = 384
linelist = []
coslist = []
natoms = 28800
nbeads = 36
n = nbeads-1
strain = [0]*tframes
f = [0]*tframes
        
  
for i in range(tframes):
    sys.stdin.readline() 
    sys.stdin.readline()
    sys.stdin.readline()
    sys.stdin.readline()
    sys.stdin.readline()
    line = sys.stdin.readline()
    [xlow,xhi] = map(float,line.split())
    if i == 0:
        xbox0 = xhi-xlow
        xbox = xbox0
    else:
        xbox = xhi-xlow
    strain[i] = (xbox-xbox0)/xbox0
    line = sys.stdin.readline()
    [ylow,yhi] = map(float,line.split())
    ybox = yhi-ylow
    line = sys.stdin.readline()
    [zlow,zhi] = map(float,line.split())
    zbox = zhi-zlow
    sys.stdin.readline()

    for j in range(natoms):
        line = sys.stdin.readline()
        [ii,molj,typej,q,xs,ys,zs,n1,n2,n3] = string.split(line)
        w=int(ii)
        x = float(xs)*xbox + (int(n1))*xbox ##unscale and unwrap coordinates
        y = float(ys)*ybox + (int(n2))*ybox
        z = float(zs)*zbox + (int(n3))*zbox
        line=[w,x,y,z,xbox]
        linelist.append(line)
           
for i in range(0,tframes*natoms,nbeads):
        x1 = linelist[i][1]
        x2 = linelist[i+n][1]
        y1 = linelist[i][2]
        y2 = linelist[i+n][2]
        z1 = linelist[i][3]
        z2 = linelist[i+n][3]
        ABx = (x2-x1) 
        ABy = (y2-y1) #E2E vecotor coordinates#      
        ABz = (z2-z1)
        AB = math.sqrt((ABx**2)+(ABy**2)+(ABz**2))
        ABdotBC = ABx*xbox
        cosT = (ABdotBC)/((AB)*xbox)
        cosT2 = cosT*cosT
        coslist.append(cosT2)

for i in range(0, len(coslist), 800):
    j=i/800
    mean = np.mean(coslist[i:i+800])
    f[j] = ((3*mean)-1)/2


print "strain orientation"
for k in range(tframes):
    print strain[k],f[k]


            
            
                
            
        
