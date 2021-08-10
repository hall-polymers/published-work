#Plot 2D image of nanoparticle from excel data file
#J.Ethier
#August 3rd, 2016
#python plotcontour.py


#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


import sys,os
import numpy as np
import pandas as pd
import pylab as py
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from math import *

yrange = 600
ymax = 600
volume = 0
vmax = 0.20
#Initialize figure
fig, ax = plt.subplots(figsize=(6,6))

#Call excel sheet 
excelsheet = './2PGN_avgdens_rad.xlsx'
sheet = 'Sheet1'
poldata=pd.read_excel(excelsheet,sheet)

#Printing title of each column
for k in poldata.keys():
    print k

#Printing the size of one of the columns
print poldata['xgrid'].size

#Printing the first row of data
print poldata['xgrid'][0],poldata['ygrid'][0],poldata['zgrid'][0]

#Calculate volume of PGN canopy
for z in range(poldata['zgrid'].size):
    if poldata['zgrid'][z] < 0:
        volume += 0 
    else:
        volume += poldata['zgrid'][z]
print volume

#Setting x and y into 1D arrays from 2D array
ydata = np.array([poldata['ygrid'][i] for i in range (yrange)])
xdata = []
poldata['xgrid'][-1] = -1
for x in range(poldata['xgrid'].size-1):
    if poldata['xgrid'][x] != poldata['xgrid'][x-1]: xdata.append(poldata['xgrid'][x])
    else: continue
xdata = np.array(xdata)
print 'xdata = ', xdata
print 'ydata = ', ydata

addx = (ymax - yrange) 
addy = (ymax - yrange)

#Initialize z values in 2D array
zplot = np.zeros(((len(xdata) + addx),(len(ydata) + addy)))
cnt=0

#Set z values into 2D array from excel sheet
for xdat in range(0,addx/2):
    for ydat in range(0,addy/2):
        zplot[xdat][ydat] = 0

for xdat in range(addx/2,len(xdata)+addx/2):
    for ydat in range(addy/2,len(ydata)+addy/2):
        if poldata['zgrid'][cnt] < 0:
            zplot[xdat][ydat] = 0
            cnt+=1
        else:
            zplot[xdat][ydat] = poldata['zgrid'][cnt]
            cnt += 1

for xdat in range(len(xdata)+addx/2,len(xdata) + addx):
    for ydat in range(len(ydata)+addy/2,len(ydata) + addy):
        zplot[xdat][ydat] = 0
        
        
#Plot z values in color grid
surf = plt.imshow(zplot,interpolation='nearest',origin='lower',extent=[-(ymax-1),(ymax-1),-(ymax-1),(ymax-1)],aspect='equal',cmap=cm.coolwarm,vmin=0,vmax=vmax)
#ax.set_aspect(1)
ax.set_xticks([])
ax.set_yticks([])
#plt.xlabel('Lateral X')
#plt.ylabel('Lateral Y')
#plt.title(sheet)
v = [0,vmax]

#ax.set_zlim(-1,25)
#fig.colorbar(set_major_locator(LinearLocator(10))
#ax.colorbar((FormatStrFormatter('%.02f')))
#surf.cmap.set_under('blue')
#surf.cmap.set_over('red')
#boundaries=[0,2,4,6,8,10,12,14,16,18]
cbar = fig.colorbar(surf,format='%s',ticks=v,extend='both',shrink=1.0,aspect=25)
cbar.ax.tick_params(labelsize=18)
plt.show()
