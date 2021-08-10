#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


### Purpose: This script can customize and sort large dump files ###
## specifically for orientation.py ##
### Syntax: python sort.py > equilsort.txt ###
### Author: Janani Sampath ###
### Date:  February 2015 ###
 

import sys,string
from numpy import *
from math import *
import fileinput

linelist=[]
infiles=['equilnvt.lammpstrj']
IN=fileinput.input(infiles)
tframes = 2
tatoms = 39361

for step in range(tframes):
    linelist=[]
    l1=IN.readline()
    l2=IN.readline() 
    l3=IN.readline()
    l4=IN.readline()
    l5=IN.readline()
    l6=IN.readline()
    l7=IN.readline()
    l8=IN.readline()
    l9=IN.readline()

    
    for j in range(1,tatoms+1):
        line = IN.readline()
        [ii,molj,typej,qi,x1,x2,x3,n1,n2,n3] = string.split(line)
        k=int(ii)
        mol=int(molj)
        typek=int(typej)
        q = int(qi)
        x=float(x1)
        y=float(x2)
        z=float(x3)
        m1=int(n1)
        m2=int(n2)
        m3=int(n3)
        line=[k,mol,typek,q,x,y,z,m1,m2,m3]
        if (typek == 2 or typek == 1 or typek == 5 or typek == 3 or typek == 6): 
            linelist.append(line)
            
            
    linelist.sort()
            
    ## for i in range(len(linelist)-2):
    ##     if ((linelist[i][2]==4) and (linelist[i+1][2]==5) and (linelist[i+2][2]==5)):
    ##         linelist[i][2] = 1
    ##         linelist[i+1][2] = 2
    ##         linelist[i+2][2] = 3
            
            
      #print("new time step")


    print l1,
    print l2,
    print l3,
    print 28800
    print l5,
    print l6,
    print l7,
    print l8,
    print l9, 
    for line in linelist:
        print ' '.join(map(str, line))    


        
    
    
