#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

# Script:  make_PEGDD_ionene_initial_structure.py
# Purpose: Make .data file containing initial structure for PEG/DD ionenes
# Usage: python3 make_PEGDD_ionene_initial_structure.py
# Author:  Nicholas Liesen - Based on script from Prasant Vijayaraghavan
# Created: 09/2018


import math
import os
import string
import sys
from random import random
from random import *
import matplotlib.pyplot as plt
import numpy as np  # So that I can reference the data types.
# For visualization
from mpl_toolkits import mplot3d as mpl
from numpy import *

wdir = '/home/nlognick/Dropbox/research/Share_Research_Progress/Hall_Project/nanocomp_sys_creation/PEO_system/commented_5_7_21'
os.chdir(wdir)


def quickViz(x_list, y_list, z_list):
    """ This function visualizes the cartesian coordinates of the beads in 3D space

    Parameters:
    x_list: x-coordinates for all position vectors
    y_list: y-coordinates for all position vectors
    z_list: z-coordinates for all position vectors

    Returns:
    None
    """

    figure = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(x_list, y_list, z_list, marker='.')
    plt.show()
    

def imageFlags(x, hi):
    """
    Generates an image flag for a position x passed into the function. Takes into account PBC.
    Adjusts particle positions so that they are within a box of size (-hi/2, hi/2)

    Usage: (position variable, relevant dimension length)

    Parameters:
    x: component of unwrapped position vector
    hi: box length, Lx

    Returns:
    x: component of wrapped position vector
    cx: component of image flag
    """

    hi2 = hi/2    # The factor of hi/2 will just change the coordinates from [0,hi) to [-hi/2, hi/2)
    if x > hi:
        cx = int(x/hi)    #for (hi, 2hi) cx = 1, for [2hi, 3hi) cx=2 etc.
        x = x - cx*hi - hi2
    elif x < 0.0:
        cx = -int((-x+hi)/hi)    #(-hi, 0) cx = -1, for (-2hi, -hi] cx = -2, etc.
        x = x - cx*hi - hi2
    else:
        cx = 0
        x = x - hi2    
    return x, cx


def posReject(x_raw, y_raw, z_raw, hx, hy, hz, np_r, tol):
    """ 
    Checks the PBC corrected position to make sure it doesn't overlap with the nanoparticle  
    ___which is always at the center of the box in this code___ (so really we are just measuring
    distance from the wrapped coords to the box's center). This is an inefficient way to check
    positions since we will be calculating the image flags and PBC corrected positions twice.

    In this code the NP is always located at center of the box (0, 0, 0).
    As a consequence measured distance can never exeed half the box length.
    If the particle were in any other position, calculating the distance between
    the NP and the placed bead would be more complicated due to the periodic
    boundary conditions, and distance would need to be measured differently.
    
    Usage: Input uncorrected/raw (x,y,z) positions, box dimensions, nanoparticle radius,
    a tolerance (what length the tip of the position vector must exceed the nanoparticle
    radius by).

    Parameters:
    x_raw, y_raw, z_raw: position vector in unwrapped coordinates (float, float, float)
    hx, hy, hz: box dimensions (float, float, float)
    np_r: radius of nanoparticle (float)

    Returns:
    reject_position: True/False for whether position should be rejected (bool)
    r_center: distance from box center/center of NP
    Note: This function relies on imageFlags() function to correct positions
    """   

    x_pbc, _ = imageFlags(x_raw, hx)
    y_pbc, _ = imageFlags(y_raw, hy)
    z_pbc, _ = imageFlags(z_raw, hz)
    
    r_center = sqrt(x_pbc**2 + y_pbc**2 + z_pbc**2)
    if r_center < np_r + tol:
        reject_position = True
    else:
        reject_position = False
    return reject_position, r_center


def randomWalk(first_bead, bond = None):
    """
    Takes in information on whether this is the first bead in a molecule (or not) & the bond length
    and either starts the chain, or grows the chain via a random walk. The following references describe
    the method used for the random walk.
    refs: https://mathworld.wolfram.com/SpherePointPicking.html
    https://math.stackexchange.com/questions/44689/how-to-find-a-random-axis-or-unit-vector-in-3d

    Parameters:
    first_bead: True/False depending on whether bead is the first in the molecule (bool)
    bond: bond length (float)
    Returns:

    first_bead (see above)
    dx, dy, dz: (x, y, z) position vector if first bead; otherwise it corresponds to the
    bond vector
    """   

    if first_bead:  # 1st bead placed at random location in box
        dx, dy, dz = random.random()*hx, random.random()*hy, random.random()*hz    #Return random value between (0,hi)
        first_bead = False
    else:
        # Generates vector on unit sphere surrounding previous interaction site and than
        # rescales to bond length.
        theta = random.random()*2*pi  # theta dist uniform over [0, 2pi]
        dz = random.random()*2 - 1  # uniform dist over cos(phi) from [-1,1]
        dx = sqrt(1-dz**2)*cos(theta)
        dy = sqrt(1-dz**2)*sin(theta)
        r = sqrt(dx**2+dy**2+dz**2) # Calculate length of vector
        scale = bond/r
        dx, dy, dz = scale*dx, scale*dy, scale*dz    #Rescale vector to bond length
    return first_bead, dx, dy, dz


def randomPoly(npoly, fracPEG):
    """
    Pass in the number of polymer chains you want to generate, and the fraction
    of PEG beads relative to the total beads in the chain. This function will generate
    npoly polymer chains through a random polymerization process

    Parameters:
    npoly: number of polymer chains to generate (int)
    fracPEG: fraction of PEG segments/all segments (float)

    Returns:
    seq: contains sequence of beads for all npoly chains (list of lists)
    Totallength: total number of beads contained within all polymer chains (int)
    """

    PEGsp = 23  # PEG beads in full segment; varies with molecular weight
    hPEG = 11  # PEG beads in half segment
    nPEG = fracPEG*npoly  # Total number of PEG segments / 3
    nPE = npoly - nPEG    # npoly(1-fracPEG) = Total number of PE segments / 3

    a=[2,2]
    b=[1]*PEGsp
    c=[1]*hPEG

    # Segment Definitions
    PE = [2,2,2,2]   # Polyethylene   (PEbead)4
    hPE = [2,2]       # half-polyethylene   (PEbead)2
    PEG = a+b+a       # poly(ethylene glycol)     (PE)2-(PEG)23-(PE)2
    fhPEG = c+a       # first half poly(ethylene glycol)     (PEG)23-(PE)2
    shPEG = a+c       # second half poly(ethylene glycol)    (PE)2-(PEG)12
    IG = [4,2,2,4]    # ionic group    AMbead-(PE)2-AMbead

    seq, tfhPE, tshPE, tfhPEG, tPEG, tPE, tshPEG = ([] for i in range(7))

    # Create pool of # segments per chain * number of chains (3*npoly) segments.
    # The total PEG segments are fracPEG * 3npoly. We run a loop (1, fracPEG*npoly)
    # and consider the factor of 3 by adding PEG segments into all 3 arrays (1st half,
    # middle, 2nd half). Use the same trick to add PE segments, using frac PE = 1-fracPEG.
    # Here, we count only the 3 segments which can be PE or PEG. The ammonium groups aren't considered.

    # Create a list of half and full PEG units
    for ix in range(int(nPEG)):
        tfhPEG.append(fhPEG)
        tPEG.append(PEG)
        tshPEG.append(shPEG)

    # Create list of half & full PE units
    for ix in range(int(nPE)):
        tfhPE.append(hPE)
        tshPE.append(hPE)
        tPE.append(PE)

    # Combine list of 1st, 2nd, 3rd polymer segments
    fgroup = tfhPEG + tfhPE
    sgroup = tPEG + tPE
    tgroup = tshPEG + tshPE

    # Now that we have mixed the appropriate ratios of PEG:PE mixed into each
    # of the three segments of the the polymer chain, we want to randomize our
    # arrays for each segment (fgroup, sgroup, tgroup) so that any entry in the
    # array that is selected will randomly be PE or PEG w/ the appropriate probability.
    shuffle(fgroup)
    shuffle(sgroup)
    shuffle(tgroup)

    Totallength=0
    for ix in range(npoly):
        seq1 = []
        seq1 = fgroup[ix]+IG+sgroup[ix]+IG+tgroup[ix]
        # Create a single polymer chain where 1st, 2nd, & 3rd segments are randomly selected.
        Totallength = Totallength+len(seq1) # records total # beads from all PE, PEG, and Ammonium segments
        seq.append(seq1)
    return seq, Totallength


# ---------------Input & system parameters--------------------

npoly = 1600  # num polymer chains
fracPEG = 0.3  # fraction PEG segments/all segments
ncounterions = 4*npoly  # 4 counterions per chain

# -----------------Physical properties-----------------

cisize = 1.0  # counter-ion diameter/normal-bead diameter
z_c = 1  # counterion charge

np_radius = 0
np_tol = 0.0
#np_mass = 17307.0  # Generating system w/out NP
# Moves are rejected w/in sphere of r = np_radius + np_tol

dens = 0.85  # bead density
bond = 0.97  # bond length (depends on bond potential ~ close to 1 good enough?)

# type names:    Poly-ethylene oxide, Polyethylene ,counterions, ammonium group
#            0 		  1	            2            3              4       
#ntype = 5
ntype = 4
atype = ('00', 'PEG', 'PE', 'IO', 'CG')


# dictionary between type and charge

charge = {
        1:0,
        2:0,
        3:-1,
        4:+1,
        }


"""charge = {
        1:0,
        2:0,
        3:-1,
        4:+1,
        5:0
        }
"""

mass = {
        1:1.0,
        2:1.0,
        3:1.82,
        4:1.0,
        }

"""mass = {
        1:1.0,
        2:1.0,
        3:1.82,
        4:1.0,
        5:np_mass
        }
"""

seq, Totallength = randomPoly(npoly, fracPEG)  # returns total length and list of polymer chain configurations

# -------------- Simulation Cell Parameters -------------------

ntot = Totallength + ncounterions  # Total number of beads in system excluding AU NP
#ntot_plus_np = ntot + 1

dim = ntot + 1    

# vol = (ntot-ncounterions+ncounterions*cisize**3)/dens 
vol_for_poly = ntot/dens
vol_for_np = (4./3.)*pi*np_radius**3
vol = vol_for_poly + vol_for_np

side = vol**(1./3.)
hx, hy, hz = side, side, side
hx2, hy2, hz2 = hx/2, hy/2, hz/2

# ----------------------- Random Walk for polymers & polymer bead parameter assignment -------------------------

xc, yc, zc = zeros(dim, np.float32), zeros(dim, np.float32), zeros(dim, np.float32)    # Bead Positions
cx, cy, cz = zeros(dim), zeros(dim), zeros(dim)    # Bead image flags
rc = zeros(dim, np.float32)                          # Radius of position vector for beads

# bead type, which polymer chain owns bead, bead charge
typeb, molnum, q = [0]*dim, [0]*dim, [0.0]*dim

# assign positions disregarding box size when growing the chain after the 1st bead
k = 0  # Counter for bead ID
lengthcurrentpoly = 0  # lengthcurrentpoly = length of beads on polymer chains
for ix in range(npoly):
    newPoly = True
    for iz in seq[ix]:
        k = k + 1
        lengthcurrentpoly = lengthcurrentpoly + 1
        typeb[k] = iz  # record bead type for each bead
        molnum[k] = ix + 1  # record which polymer chain each bead belongs to
        if typeb[k] == 4:  # record the charge of each bead
            q[k] = 1
        else: q[k] = 0
        if newPoly:
            reject_position = True
            while reject_position == True:
                result = randomWalk(newPoly, bond)  # Returns (newPoly, x position, y position, z position)
                x_test, y_test, z_test = result[1], result[2], result[3]
                # test the position to see if it overlaps with the nanoparticle
                reject_position, _ = posReject(x_test, y_test, z_test, hx, hy, hz, np_radius, np_tol)

            xc[k], yc[k], zc[k] = x_test, y_test, z_test  # Assign test position as final if no overlap
            newPoly = result[0]  # Returns newPoly as false
            
        else:
            reject_position = True
            while reject_position == True:  # Need to write in code to exit if the loop repeats too many times
                #result = randomWalk(newPoly, bond)
                #dx, dy, dz = result[1], result[2], result[3]
                _, dx, dy, dz = randomWalk(newPoly, bond)
                # test the position to see if it overlaps with the nanoparticle
                x_test, y_test, z_test = xc[k-1] + dx, yc[k-1] + dy, zc[k-1] + dz
                reject_position, _ = posReject(x_test, y_test, z_test, hx, hy, hz, np_radius, np_tol)

            xc[k], yc[k], zc[k] = x_test, y_test, z_test  # Assign test position as final when no overlap is present

print("Polymers Built \n\n")

# ---------------------------------------- Random Placement of Counterion Beads ---------------------------------------

no_chains = True  # flag used with randomWalk function
for ii in range(1, ncounterions+1):
    k = ii + ntot - ncounterions
    typeb[k] = 3
    q[k] = z_c*(-1.0)
    reject_position = True
    while reject_position == True:
        result = randomWalk(no_chains)
        x_test, y_test, z_test = result[1], result[2], result[3]
        reject_position, _ = posReject(x_test, y_test, z_test, hx, hy, hz, np_radius, np_tol) 
        
    xc[k], yc[k], zc[k] = x_test, y_test, z_test  # Assign test position when there is no overlap w/ NP
    rc[k] = sqrt(xc[k]**2+yc[k]**2+zc[k]**2)
    
print("Counterions Placed")

# ------------------------------------- Correct all positions using image flags ---------------------------------------

for k in range(1,ntot+1):
    xc[k], cx[k] = imageFlags(xc[k], hx)
    yc[k], cy[k] = imageFlags(yc[k],hy)
    zc[k], cz[k] = imageFlags(zc[k],hz)
    rc[k] = sqrt(xc[k]**2+yc[k]**2+zc[k]**2)  

print("Positions corrected for PBC & image flags generated")

# --------------------------------------- Visualize & Messages -----------------------------------------------

quickViz(xc[:], yc[:], zc[:])

print("Total number of particles:", ntot,"\n\n", "Number of Chains = ", npoly, "monomers total =", Totallength, "\n\n")
print("Number of counter ions = ", ncounterions)

print("dens = ", dens, ", vol = ", vol, "\n\n")
            
# --------------------------------------- Lammps Final Output ---------------------------------------------------------

nbonds = ntot - ncounterions - npoly

# input.lammps header 

with open('1kPEO_1600chains_fracp3_homo_small.data', 'w') as INPUT_LAMMPS:
    INPUT_LAMMPS.write("\n")
    #INPUT_LAMMPS.write("%10i    atoms\n" %     ntot_plus_np)
    INPUT_LAMMPS.write("%10i    atoms\n" %     ntot)
    INPUT_LAMMPS.write("%10i    bonds\n" %     nbonds)
    INPUT_LAMMPS.write("%10i    angles\n" %     0)  
    INPUT_LAMMPS.write("%10i    dihedrals\n" % 0)
    INPUT_LAMMPS.write("%10i    impropers\n" % 0)
    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write("%10i    atom types\n" %  ntype)
    INPUT_LAMMPS.write("%10i    bond types\n" % 1)
    
    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write(" %16.8f %16.8f   xlo xhi\n" % (-hx2, hx2))
    INPUT_LAMMPS.write(" %16.8f %16.8f   ylo yhi\n" % (-hy2, hy2))
    INPUT_LAMMPS.write(" %16.8f %16.8f   zlo zhi\n" % (-hz2, hz2))
    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write("Atoms\n")
    INPUT_LAMMPS.write("\n")

    # --------------------- End Output Headers --------------------

    # Polymers
    i = 0
    imol = 0
    
    # Positions
    for i in range(1, ntot+1):
        itype = typeb[i]
        aname = atype[itype]    # Recalls atom type using bead type as key
    
        if itype != 3:
            imol = molnum[i]    # Recall which chain bead belongs to
            segname = "POLY"
        elif itype == 3:
            imol = i-ntot+ncounterions+npoly  # LMH now each ion has its own molecule number
            segname = "CION"
    
        INPUT_LAMMPS.write("%6i %6i %2i %6.2f %9.4f %9.4f %9.4f %6i %6i %6i\n" % (i, imol, typeb[i], q[i], xc[i], yc[i], zc[i], cx[i], cy[i], cz[i]))
    

    #typeb[0] = 5; itype = typeb[0]
    #INPUT_LAMMPS.write("%6i %6i %2i %6.2f %9.4f %9.4f %9.4f %6i %6i %6i\n" % (ntot+1, imol+1, typeb[0], charge[itype], xc[0], yc[0], zc[0], cx[0], cy[0], cz[0]))
    # imol from previous loop will be left @ highest molecule-ID generated
    # To keep things 1-indexed they left index 0 open in the array. I have borrowed this position for my NP coordinates. I want my NP at the center
    # of the box, so the coordinates and image flags are kept as (0,0,0)

    # Bonds
    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write("Bonds\n")
    INPUT_LAMMPS.write("\n")

    ibond=0 
    for i in range(1,ntot-ncounterions+1):
            # if not at the end of the polymer
            if molnum[i+1] == molnum[i]:
                ibond = ibond+1 # the bond number
                j=i+1
                INPUT_LAMMPS.write("%8i  1 %8i %8i\n" % (ibond,i,j))
                # Bond ID, Bond type, 1st atom in bond, 2nd atom in bond
                # Note: There is only a single bond type (FENE)
    
    # Masses
    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write("Masses\n")
    INPUT_LAMMPS.write("\n")
    
    for key in mass.keys():  # Loop over all atom-type masses
        INPUT_LAMMPS.write("%3i  %10.4f\n" % (key, mass[key]))

    # For homo simulation can remove NP mass from input script

print("LAMMPS output complete.")
