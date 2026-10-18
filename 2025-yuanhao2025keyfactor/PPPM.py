from numpy import *
import math
import fileinput
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

nc = 150
clength = 50
nconf = 200
filename = "4_nvt.lammpstrj"

def func(x,a):
    return exp(-a*x)


def MSD(n_chain, chainlength, nconf):
    rs = open('tmp.data', 'w')

    tyid = [0] * 300001
    position_x = [0] * 300001
    position_y = [0] * 300001
    position_z = [0] * 300001
    position_zr = [0] * 300001

    Natom=0
    count = 0
    ind = 0
    step = 0
    rho = [0] * 5000
    bl = [0] * 5000

    delt_x=[0]*200
    delt_y=[0]*200
    delt_z=[0]*200
   
    
    angle=zeros((60000))
    angle_count=zeros((100))
    tmp=zeros((60000))
    ee=zeros((1000))
    c_v=0
    cc=0

    n_poly = 0
    mod=1
    k=0
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
            ind = ind + 1
            parts = Line.split()
            #print(Line)
            atom_id, mol, atom_type, charge, xs, ys, zs, x1, y1, z1 = map(float, parts)
            x = xs * (xmax-xmin) + x1 * (xmax-xmin)
            y = ys * (ymax-ymin) + y1 * (ymax-ymin)
            zr = zs * (zmax-zmin) + z1 * (zmax-zmin)
            z = zs
            tyid[int(atom_id)] = int(atom_type)

            if step >= nconf:
                position_x[int(atom_id)] = x
                position_y[int(atom_id)] = y
                position_z[int(atom_id)] = z
                position_zr[int(atom_id)] = zr
                n_poly = n_poly + 1

        ind = 0
        if mod==3:
            mod=1
            if step >= nconf:
                for cid in range(1, n_chain + 1):

                    for cl in range(1, chainlength):
                        a1 = (cid - 1) * chainlength + cl
                        a2 = (cid - 1) * chainlength + cl + 1                           
                        count = count + 1
                        delt_x[cl] = position_x[a2] - position_x[a1]
                        delt_y[cl] = position_y[a2] - position_y[a1]
                        delt_z[cl] = position_zr[a2] - position_zr[a1]


                        bl[step] = bl[step] + math.sqrt(delt_x[cl] * delt_x[cl] + delt_y[cl] * delt_y[cl] + delt_z[cl] * delt_z[cl])

                        rho[step] = Natom / ((zmax-zmin)*(ymax-ymin)*(xmax-xmin))
                    #for i in range(0,chainlength-1-1 +1):#interval between bonds
                    
                    for j  in range(1,chainlength-1-1 + 1):#starting point of bonds
                        ac1=[delt_x[j],   delt_y[j],   delt_z[j]   ]
                        ac2=[delt_x[j+1], delt_y[j+1], delt_z[j+1] ]
                        c_v+=1
                        tmp[c_v]=dot(ac1,ac2)/sqrt(dot(ac1,ac1) * dot(ac2,ac2)  ) 
                        angle[c_v]=math.acos(tmp[c_v])*180/pi
                        for i in range(0,59):
                            downlim = i*3
                            uplim = (i+1)*3
                            if (angle[c_v] <=uplim) and (angle[c_v] >downlim):
                                index=i
                        angle_count[index]+=1
                    begin = (cid - 1) * chainlength + 1
                    end = (cid) * chainlength
                    ee[step] = ee[step] + ((position_x[end] - position_x[begin]) ** 2 + \
                            (position_y[end] - position_y[begin]) ** 2 + \
                            (position_zr[end] - position_zr[begin]) ** 2)/n_chain                                              
                nstd=std(angle[angle!=0])  
                aveangle=mean(angle[angle!=0])                          

                bl[step] = bl[step] / count
                acf=mean(tmp[tmp!=0] )
                cc=c_v
                c_v=0
                             
                #rs.write('%4i %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f\n' % (step, rho[step], bl[step],ee[step] ,acf, nstd,aveangle) )
                
                count = 0
                print("doing step: %3i" %step)
    print(mean(bl[nconf:step]))
    for i in range(0,59):
        rs.write('%4i  %7.5f\n' % ( i, float(angle_count[i]/((step-nconf)*cc))) )

MSD(nc, clength, nconf)
  


    

