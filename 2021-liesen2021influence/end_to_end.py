#!/bin/env python2

#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

# Script:  end_to_end.py
# Purpose: Calculates end-to-end autocorrelation function from a .lammpstrj file
# Usage: python3 scriptname.py $number_of_polymer_chains $fraction_of_PTMO
# Created: read_lammpstrj and end2end_autocorr function from Jonathan Brown's pppmd package (2016-10-20)
# Modified: Nicholas Liesen (Spring 2018)


import os
import sys
import numpy as np


# read_lammpstrj: read in a lammps trajectory
#
# Input: fname, num_frames
#  fname: filename string or 'stdin' (or a value that evaluates to false) for reading from standard in 
#  num_frames: optional number of frames to read before stopping, defaults to reading in all frames
#  skip_beginning: skip this many frames at the beginning of the dump file
#  skip_between: skip this many frames between saved frames
#
# Output: r, ir, timestep, box_bounds, id2type, id2mol, mol2ids
#  r: num_frames by num_atoms+1 by 3 array of wrapped and unscaled coordinates (indexed by frame number then atom id)
#  ir: num_frames by num_atoms+1 by 3 array of image flags
#  timestep: num_frames length array of timesteps
#  box_bounds: 3D array to store boundaries of the box, indexed by frame, x/y/z, then lower/upper
#  id2type, id2mol: num_atoms+1 length arrays to map atom id to type and molecule id (if available, may be None)
#  mol2ids: num_mols+1 length list of atom id arrays corresponding to the molecules (if available, may be None)
#
# NOTE: assumes that the number of atoms in the simulation is fixed
# NOTE: also assumes that the coordinates are wrapped, so x or xs type coordinates are allowed but not xu or xsu
#
def read_lammpstrj(fname, num_frames=float('inf'), skip_beginning=0, skip_between=0):
    # helper function to read in the header and return the timestep, number of atoms, and box boundaries
    def read_header(f):
        f.readline() # ITEM: TIMESTEP
        timestep = int(f.readline())

        f.readline() # ITEM: NUMBER OF ATOMS
        num_atoms = int(f.readline())

        f.readline() # ITEM: BOX BOUNDS xx yy zz
        line = f.readline()
        line = line.split()
        xlo = float(line[0])
        xhi = float(line[1])
        line = f.readline()
        line = line.split()
        ylo = float(line[0])
        yhi = float(line[1])
        line = f.readline()
        line = line.split()
        zlo = float(line[0])
        zhi = float(line[1])

        return timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi


    #allow reading from standard input
    if not fname or fname == 'stdin':
        f = sys.stdin
    else:
        f = open(fname, 'r')

    # read in the initial header
    frame = 0
    init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = read_header(f)

    # skip the beginning frames, if requested
    for skippedframe in range(skip_beginning):
        f.readline() # ITEM: ATOMS
        # loop over the atoms lines
        for atom in range(num_atoms):
            f.readline()
        init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = read_header(f)

    # preallocate arrays, if possible
    if num_frames < float('inf'):
        alloc = num_frames
        inf_frames = False
    else:
        alloc = 1
        inf_frames = True
    timestep = np.zeros(alloc, np.int) # 1D array of timesteps
    box_bounds = np.zeros([alloc,3,2], np.float) # 3D array to store boundaries of the box, indexed by frame, x/y/z, then lower/upper

    timestep[frame] = init_timestep
    box_bounds[frame][0][0] = xlo
    box_bounds[frame][0][1] = xhi
    box_bounds[frame][1][0] = ylo
    box_bounds[frame][1][1] = yhi
    box_bounds[frame][2][0] = zlo
    box_bounds[frame][2][1] = zhi

    # NOTE: using num_atoms+1 here so that the arrays are indexed by their LAMMPS atom id
    r = np.zeros([alloc, num_atoms+1, 3], np.float) # 3D array of x, y, z coordinates, r[frame][id][coordinate]
    ir = np.zeros([alloc, num_atoms+1, 3], np.int) # 3D array of x, y, z image flags, r[frame][id][coordinate]

    id2mol = np.zeros(num_atoms+1, np.int) # array to map from atom id to molecule id, builds this from the first frame, if available
    id2type = np.zeros(num_atoms+1, np.int) # array to map from atom id to type, builds this from the first frame, if available


    # separately do the first ATOMS section so that we can initialize things, build the id2mol and id2type arrays, and so that the main loop starts with reading in the header
    line = f.readline()
    line = line.split()
    id_index = line.index("id") - 2
    if "mol" in line:
        mol_index = line.index("mol") - 2
    else:
        mol_index = None
    if "type" in line:
        type_index = line.index("type") - 2
    else:
        type_index = None

    if "x" in line:
        scaled = False
        x_index = line.index("x") - 2
        y_index = line.index("y") - 2
        z_index = line.index("z") - 2
    elif "xs" in line:
        scaled = True
        x_index = line.index("xs") - 2
        y_index = line.index("ys") - 2
        z_index = line.index("zs") - 2
    else:
        print >> sys.stderr, "ERROR: x coordinate not found in lammps trajectory"
        return

    if "ix" in line:
        ix_index = line.index("ix") - 2
        iy_index = line.index("iy") - 2
        iz_index = line.index("iz") - 2
    else:
        print >> sys.stderr, "ERROR: x image flag not found in lammps trajectory"
        return

    # loop over the atoms lines for the first frame separately, the rest of the frames will be read in below
    for atom in range(num_atoms):
        line = f.readline()
        line = line.split()

        # get the atom id
        my_id = int(line[id_index])

        # x, y, z coordinates
        r[frame][my_id][0] = float(line[x_index])
        r[frame][my_id][1] = float(line[y_index])
        r[frame][my_id][2] = float(line[z_index])

        # unscale, if necessary
        if scaled:
            r[frame][my_id][0] = r[frame][my_id][0]*(box_bounds[frame][0][1]-box_bounds[frame][0][0]) + box_bounds[frame][0][0]
            r[frame][my_id][1] = r[frame][my_id][1]*(box_bounds[frame][1][1]-box_bounds[frame][1][0]) + box_bounds[frame][1][0]
            r[frame][my_id][2] = r[frame][my_id][2]*(box_bounds[frame][2][1]-box_bounds[frame][2][0]) + box_bounds[frame][2][0]

        # x, y, z image flags
        ir[frame][my_id][0] = int(line[ix_index])
        ir[frame][my_id][1] = int(line[iy_index])
        ir[frame][my_id][2] = int(line[iz_index])

        # if available, buidl the i2mol and id2type arrays
        if mol_index is not None:
            id2mol[my_id] = int(line[mol_index])
        if type_index is not None:
            id2type[my_id] = int(line[type_index])

    # build the reverse of the id2mol array
    # this is a 2D array with rows of (potentially) varying length, so nest a numpy array into a python list
    if mol_index is not None:
        num_mols = id2mol.max()	 # NTL: Finds total number of molecules
        mol2ids = [[]]  # NTL: This is a list of numpy arrays. Each entry in the list will contain an array corresponding
        # to the atom IDs which are a part of the moleculeID corresponding to the list's index where
        # the array is contained.
        for molid in range(1, num_mols+1):
            mol2ids.append(np.where(id2mol==molid)[0])
    else:
        num_mols = None
        mol2ids = None

    # loop over number of num_frames frames, if num_frames is infinite, will loop over all the frames in the file
    frame = 1 # this is the frame counter for frames actually read in
    frame_attempt = 0 # this is the actual frame count in the file (not counting the ones skipped in the beginning
    while frame < num_frames:

        frame_attempt += 1

        # try to read in a new header
        try:
            my_timestep, my_num_atoms, my_xlo, my_xhi, my_ylo, my_yhi, my_zlo, my_zhi = read_header(f)
        except:
            print >> sys.stderr, "WARNING: hit end of file when reading in", fname, "at frame", skip_beginning + frame_attempt
            break

        # skip the frame if between frames to be read in and restart the loop
        if frame_attempt%(skip_between+1) > 0:  # NTL: Only read in frames separated skip_between values from the prev read in frame.
            f.readline() # ITEM: ATOMS
            # loop over the atoms lines
            for atom in range(num_atoms):
                f.readline()
            continue

        if inf_frames:
            timestep = np.append(timestep, 0)

            box_bounds = np.concatenate( ( box_bounds, np.zeros([1,3,2],np.float) ) )

            r = np.concatenate( ( r, np.zeros([1, num_atoms+1, 3], np.float) ) )
            ir = np.concatenate( ( ir, np.zeros([1, num_atoms+1, 3], np.float) ) )

        # update the timestep and box size arrays
        timestep[frame] = my_timestep
        box_bounds[frame][0][0] = my_xlo
        box_bounds[frame][0][1] = my_xhi
        box_bounds[frame][1][0] = my_ylo
        box_bounds[frame][1][1] = my_yhi
        box_bounds[frame][2][0] = my_zlo
        box_bounds[frame][2][1] = my_zhi

        f.readline() # ITEM: ATOMS
        # loop over the atoms lines
        for atom in range(num_atoms):
            line = f.readline()
            line = line.split()

            # get the atom id
            my_id = int(line[id_index])

            # x, y, z coordinates
            r[frame][my_id][0] = float(line[x_index])
            r[frame][my_id][1] = float(line[y_index])
            r[frame][my_id][2] = float(line[z_index])

            # unscale, if necessary
            if scaled:
                r[frame][my_id][0] = r[frame][my_id][0]*(box_bounds[frame][0][1]-box_bounds[frame][0][0]) + box_bounds[frame][0][0]
                r[frame][my_id][1] = r[frame][my_id][1]*(box_bounds[frame][1][1]-box_bounds[frame][1][0]) + box_bounds[frame][1][0]
                r[frame][my_id][2] = r[frame][my_id][2]*(box_bounds[frame][2][1]-box_bounds[frame][2][0]) + box_bounds[frame][2][0]

            # x, y, z image flags
            ir[frame][my_id][0] = int(line[ix_index])
            ir[frame][my_id][1] = int(line[iy_index])
            ir[frame][my_id][2] = int(line[iz_index])

        frame += 1

    return r, ir, timestep, box_bounds, id2type, id2mol, mol2ids


# end2end_autocorr: end to end autocorrelation function
#
# Input: r, ir, box_bounds, mol2ids
#  r: unscaled (but wrapped) coordinates
#  ir: image flags
#  box_bounds: boundaries of the box
#  mol2ids: list of atom id arrays of the molecules to be used
#  (format as read in from read_lammpstrj)
#
# Output: end to end vector autocorrelation function (1D array indexed by frame count) constructued by dotting the original end to end vector to the 
#
# NOTE: all the listings in mol2ids will be used and averaged together
# NOTE: it is assumed that the end-to-end vector is the one between the lowest and highest id in each molecule (if this is not the case, you'd have to mess with mol2ids, e.g. make it only contain the ids of the two end beads)
# NOTE: scaled by the average end-to-end vector at frame 0, so that e2e_autocorr[0]=1.0
#
def end2end_autocorr(r, ir, box_bounds, mol2ids):
    frames = len(r)  # NTL: number of timesteps/frames/snapshots read in
    mols = len(mol2ids)  # NTL: number of molecules

    # preallocate e2e vector arrays and autocorr array
    e2e_t = np.zeros([mols, 3], np.float)  # NTL: Will store the 3 distance components of the end to end vector for each molecule (1st bead -> last bead) @ some time t
    e2e_0 = np.zeros([mols, 3], np.float)
    e2e_autocorr = np.zeros(frames, np.float)  #NTL: For each frame/timestep, store the autocorrelation function's value @ that time.

    # loop over time
    for t in range(frames):  # NTL: For the 2nd index 0=x, 1=y, 2=z --- for the 3rd index 0=lower & 1=upper.
        #				   x, upper  	minus	x, lower	       y, upper	    minus    y, lower
        box_size = np.array([ box_bounds[t][0][1] - box_bounds[t][0][0], box_bounds[t][1][1] - box_bounds[t][1][0], box_bounds[t][2][1] - box_bounds[t][2][0] ])

        # loop over molecules
        for molid in range(1, mols):
            # Assume that the ends of the chain have the maximum and minimum id numbers
            id1 = mol2ids[molid].min()
            id2 = mol2ids[molid].max()

            # calculate the end-to-end vector
            r1 = r[t][id1] + ir[t][id1]*box_size  # NTL: Convert position from wrapped value using image flags (unwrap coords)
            r2 = r[t][id2] + ir[t][id2]*box_size

            e2e_t[molid] = r2 - r1
            if t == 0:
                e2e_0[molid] = e2e_t[molid]

            # take dot products
            e2e_autocorr[t] += np.dot(e2e_0[molid], e2e_t[molid])  # Sum up all the dot products for each molecule in order to average them

    # scaling
    e2e_autocorr = e2e_autocorr/(mols-1)
    e2e_autocorr = e2e_autocorr/e2e_autocorr[0]  # NTL: First entry in e2e_autocorr will be <Ree(0) dot Ree(0)>. Use to calculate <Ree(t) dot Ree(0)>/<Ree(0) dot Ree(0)>

    return e2e_autocorr

in_fname = "PEO2k_75pcent_log_dump.lammpstrj"
out_fname = "PEO2k_75pcent_e2e.acf"
r, ir, timestep, boxbds, id2type, id2mol, mol2ids = read_lammpstrj(in_fname,250)

window_size = 45  # Frames
stop_after = 43  # Frames
low = 0
high = low + stop_after
frame_zero = low
blzero_steps = timestep[low:high]
eeacf = end2end_autocorr(r[low:high], ir[low:high], boxbds[low:high], mol2ids)

num_blocks = 3
for n in range(1,num_blocks):
    low = (window_size - 1)*n
    high = low + stop_after
    eeacf += end2end_autocorr(r[low:high], ir[low:high], boxbds[low:high], mol2ids)

eeacf /= float(num_blocks)

delta_t = np.zeros(len(eeacf), np.float)
for n in range(len(eeacf)):
    delta_t[n] = blzero_steps[n+frame_zero]-blzero_steps[frame_zero]  # only works up to point where log spacing restarts
    print delta_t[n], eeacf[n]

with open(out_fname, "a") as f:
    for n in range(len(eeacf)):
        f.write(str(delta_t[n])+'    '+str(eeacf[n])+"\n")
