from numpy import *
import sys
import math
import fileinput
#from scipy.optimize import curve_fit
#import matplotlib.pyplot as plt

nc = 500
clength = 100

nconf = 50

blockint = 0
overallf = 200
filename = "4_nvt.lammpstrj"

def MSD(n_chain, chainlength, nconf):

    rs = open('tmp.data', 'w')
    tyid = [0] * 300001
    position_x = [0] * 300001
    position_y = [0] * 300001
    position_z = [0] * 300001
    position_zr = [0] * 300001
    MS = zeros((clength+1,overallf-nconf+1))
    step = 0
    n_poly = 0
    mod=1
    k=0
    t =-1
    xmin=xmax=ymax=ymin=zmax=zmin=0
    fp=fileinput.input(filename)
    while True:
        Line=fp.readline()
        if Line=='':
            break
        k=k+1       
        if mod==1:              
            if 'ITEM: NUMBER OF ATOMS' in Line:
                Line=fp.readline()
                Natom=int(Line)                         
            if 'ITEM: BOX BOUNDS' in Line:
                Line=fp.readline()
                [xmin, xmax] = map(float,Line.split())
                Line=fp.readline()                   
                [ymin, ymax] = map(float,Line.split())
                Line=fp.readline()
                [zmin, zmax] = map(float,Line.split())                    
            if "ITEM: ATOMS" in Line:                   
                step = step + 1
                mod=2   
                Line=fp.readline()  
        if "ITEM: TIMESTEP" in Line:
            mod=3                           
        if mod==2:                           
            parts = Line.split()
            #print(Line)
            atom_id, mol, atom_type, charge, xs, ys, zs, x1, y1, z1 = map(float, parts)
            x = xs * (xmax-xmin) + x1 * (xmax-xmin)
            y = ys * (ymax-ymin) + y1 * (ymax-ymin)
            zr = zs * (zmax-zmin) + z1 * (zmax-zmin)
            z = zs
            tyid[int(atom_id)] = int(atom_type)

            if (step >= nconf) and (atom_type ==3 or 2):
                position_x[int(atom_id)] = x
                position_y[int(atom_id)] = y
                position_z[int(atom_id)] = z
                position_zr[int(atom_id)] = zr
                n_poly = n_poly + 1
       
        if mod==3:
            mod=1
            if step >=overallf:
                break
               
            print("doing step: %3i" %step)
            if step > nconf:
                t+=1
                for cid in range(1, n_chain + 1):
                    begin = int((cid - 1) * chainlength +1)
                    for cl in range(1,chainlength+1):
                        end   = int((cid - 1) * chainlength + cl)                    
                        delt_x = position_x[begin] - position_x[end]
                        delt_y = position_y[begin] - position_y[end]
                        delt_z = position_zr[begin] - position_zr[end]
                        MS[cl,t] += (delt_x*delt_x + delt_y*delt_y + delt_z*delt_z)/(cl*n_chain) 
      
                
    for i in range(1,chainlength+1):
        tmp = MS[i,:]
        rs.write("%5.4f %5.4f\n"% (mean(tmp),std(tmp)))

MSD(nc, clength, nconf)
  


    
