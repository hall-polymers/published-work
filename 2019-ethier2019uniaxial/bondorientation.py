#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

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


tframes = 1
skiptoframe = 0
linelist = []
coslist = [0]*300
natoms = 360972
nbeads = 160
npoly = natoms/nbeads
strain = [0]*(300)
f = [0]*(300)
count = [0]*300

for i in range(0, skiptoframe):
    sys.stdin.readline() 
    sys.stdin.readline()
    sys.stdin.readline()
    sys.stdin.readline()
    sys.stdin.readline()
    sys.stdin.readline()
    sys.stdin.readline()
    sys.stdin.readline()
    sys.stdin.readline()

    for j in range(natoms):
        line = sys.stdin.readline()      
  
for i in range(0, tframes):
    sys.stdin.readline() 
    sys.stdin.readline()
    sys.stdin.readline()
    sys.stdin.readline()
    sys.stdin.readline()
    line = sys.stdin.readline()
    [xlow,xhi] = map(float,line.split())
    xbox = xhi-xlow
    line = sys.stdin.readline()
    [ylow,yhi] = map(float,line.split())
    ybox = yhi-ylow
    line = sys.stdin.readline()
    [zlow,zhi] = map(float,line.split())
    zbox = zhi-zlow
    sys.stdin.readline()

    for j in range(natoms):
        line = sys.stdin.readline()
        [ii,molj,typej,xu,yu,zu] = string.split(line)
        typemon = int(typej)
        w=int(ii)
        x = float(xu)  ##unscale and unwrap coordinates
        y = float(yu)
        z = float(zu) 
        line=[w,x,y,z,xbox,typemon]
        linelist.append(line)
           
for i in range(0,natoms-1):
    if linelist[i][5] == 3 or linelist[i+1][5] == 1:
        continue
    else:
        x1 = linelist[i][1]
        x2 = linelist[i+1][1]
        y1 = linelist[i][2]
        y2 = linelist[i+1][2]
        z1 = linelist[i][3]
        z2 = linelist[i+1][3]
        ABy = (y2-y1)           #bond vector coordinates 
        if abs(ABy) > 5:        #check if two beads are straddling the periodic boundary (check distance must be some number larger than 1.5 which is maximum extent of bond)
            continue            #do not include these bond vectors in the calculation
        else:
            avgloc = int((y2+y1)/2)
            count[avgloc] += 1
            ABx = (x2-x1)               #bond vector coordinates     
            if abs(ABx) > 5:            #check if two beads are straddling the periodic boundary
                if ABx < 0:
                    ABx = ABx + xbox 
                elif ABx > 0:
                    ABx = ABx - xbox
            ABz = (z2-z1)                              #bond vector coordinates
            AB = math.sqrt((ABx**2)+(ABy**2)+(ABz**2)) #bond vector between two beads
            ABdotBC = ABz*zbox                         #dot product of bond vector with (x, y, or z box coordinate)
            cosT = (ABdotBC)/((AB)*zbox)               #cosine theta (angle between bond vector and box coordinate)
            cosT2 = cosT*cosT                          #squared
            coslist[avgloc] += cosT2                   #sum of all values in (avgloc) bin

for i in range(0, len(coslist)):
    if count[i] == 0:
        continue
    else:
        mean = coslist[i]/count[i]
        f[i] = ((3*mean)-1)/2


print "strain orientation"
for k in range(0,len(f)):
    print k,f[k]


            
            
                
            
        
