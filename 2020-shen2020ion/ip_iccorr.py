#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#!/usr/bin/python
 
# Script:  residencetime.py
# Purpose: calculate ion pair/cage correlation function
# Syntax:  python residencetime.py 
# Author:  Kevin Shen, May 2019

import sys
import numpy as np 
import pandas as pd

# special_read (written based on Jon's pppmd script): read in a lammps trajectory, but only extract data of certain bead types
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
#
# NOTE: assumes that the number of atoms in the simulation is fixed
# NOTE: also assumes that the coordinates are wrapped, so x or xs type coordinates are allowed but not xu or xsu
#
def special_read(fname, types, num_frames=float('inf'), skip_beginning=0, skip_between=0):
    print "Reading configurations..."
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
        print "Skipping " + str(skippedframe) 
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

    # NOTE: using num_atoms+1 here so that the arrays are indexed by their LAMMPS atom id
    r_temp = np.zeros([num_atoms+1, 3], np.float) # 3D array of x, y, z coordinates, r[frame][id][coordinate]
    ir_temp = np.zeros([num_atoms+1, 3], np.int) # 3D array of x, y, z image flags, r[frame][id][coordinate]

    # id2mol = np.zeros(num_atoms+1, np.int) # array to map from atom id to molecule id, builds this from the first frame, if available
    index2type = np.zeros(num_atoms+1, np.int) # array to map from atom id to type, builds this from the first frame, if available
    index2id = np.zeros(num_atoms+1, np.int) # array to map from atom id to type, builds this from the first frame, if available
    
    num_types = 0
    # loop over the atoms lines for the first frame separately, the rest of the frames will be read in below
    for atom in range(num_atoms):
        line = f.readline()
        line = line.split()

        # get the atom id
        my_id = int(line[id_index])

        # build the index2type array
        my_type = int(line[type_index])
        index2type[my_id] = my_type
        if my_type in types:
            num_types += 1

        # x, y, z coordinates
        r_temp[my_id][0] = float(line[x_index])
        r_temp[my_id][1] = float(line[y_index])
        r_temp[my_id][2] = float(line[z_index])

        # unscale, if necessary
        if scaled:
            r_temp[my_id][0] = r_temp[my_id][0]*(box_bounds[frame][0][1]-box_bounds[frame][0][0]) + box_bounds[frame][0][0]
            r_temp[my_id][1] = r_temp[my_id][1]*(box_bounds[frame][1][1]-box_bounds[frame][1][0]) + box_bounds[frame][1][0]
            r_temp[my_id][2] = r_temp[my_id][2]*(box_bounds[frame][2][1]-box_bounds[frame][2][0]) + box_bounds[frame][2][0]

        # x, y, z image flags
        ir_temp[my_id][0] = int(line[ix_index])
        ir_temp[my_id][1] = int(line[iy_index])
        ir_temp[my_id][2] = int(line[iz_index])

    # NOTE: using num_types+1 here so that the arrays are indexed by their LAMMPS atom id
    num_not_types = num_atoms - num_types
    r = np.zeros([alloc, num_types+1, 3], np.float) # 3D array of x, y, z coordinates, r[frame][id][coordinate]
    ir = np.zeros([alloc, num_types+1, 3], np.int) # 3D array of x, y, z image flags, r[frame][id][coordinate]

    id2type = np.zeros(num_types+1, np.int) # array to map from atom id to type, builds this from the first frame, if available
    id2index =  np.zeros(num_types+1, np.int) # array to map from atom id to index, builds this from the first frame, if available

    # store the temporary data into real arrays
    my_id = 0
    for atom in range(num_atoms):
        index = atom+1
        if index2type[index] in types:
            my_id += 1
            # x, y, z coordinates
            r[frame][my_id][0] = r_temp[index][0]
            r[frame][my_id][1] = r_temp[index][1]
            r[frame][my_id][2] = r_temp[index][2]

            # x, y, z image flags
            ir[frame][my_id][0] = ir_temp[index][0]
            ir[frame][my_id][1] = ir_temp[index][1]
            ir[frame][my_id][2] = ir_temp[index][2]

            id2type[my_id] = index2type[index]
            id2index[my_id] = index
            index2id[index] = my_id

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
        if frame_attempt%(skip_between+1) > 0:
            f.readline() # ITEM: ATOMS
            # loop over the atoms lines
            for atom in range(num_atoms):
                f.readline()
            continue

        # if we don't know how many frames to read in, have to allocate more memeory for the arrays
        if inf_frames:
            timestep = np.append(timestep, 0)
            
            box_bounds = np.concatenate( ( box_bounds, np.zeros([1,3,2],np.float) ) )

            r = np.concatenate( ( r, np.zeros([1, num_types+1, 3], np.float) ) )
            ir = np.concatenate( ( ir, np.zeros([1, num_types+1, 3], np.float) ) )
        
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
            index = int(line[id_index])
            if index2type[index] in types:
                my_id = index2id[index]

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
    
        print "Reading frame {}".format(frame)
        frame += 1
    
    print 'Summary:'
    print 'Total number of atoms =', num_atoms
    print 'Total number of ions =', num_types

    return r, ir, timestep, box_bounds, id2type, id2index


def buildnlist(r, bin2id, id2bin, bins, boxsize, id2type, dist_range, nearest):
    print 'Building neighbor lists...'
    #loop over the local area of each bead, and shift at the periodic BCs
    id2neighbors = [[[] for n in xrange(len(r[0]))]for n in xrange(len(r))]
    for t in range(len(r)):
        # loop over atoms
        for atom in range(1, len(r[0])):
            shift = np.zeros(3, np.float)
            binloc = id2bin[t][atom]
            if nearest:
                dmin = dist_range**2
                dmin_id = 0
            for i in range(binloc[0]-1,binloc[0]+2):
                if i == -1:
                    i = bins[0]-1
                    shift[0] = -boxsize[0]
                elif i == bins[0]:
                    i = 0
                    shift[0] = boxsize[0]
                else:
                	shift[0] = 0
            
                for j in range(binloc[1]-1,binloc[1]+2):
                    if j == -1:
                        j = bins[1] - 1
                        shift[1] = -boxsize[1]
                    elif j == bins[1]:
                        j = 0
                        shift[1] = boxsize[1]
                    else:
                		shift[1] = 0

                    for k in range(binloc[2]-1,binloc[2]+2):
                        if k == -1:
                            k = bins[2] - 1
                            shift[2] = -boxsize[2]
                        elif k == bins[2]:
                            k = 0
                            shift[2] = boxsize[2]
                        else:
                			shift[2] = 0
                            
                        #loop over the beads in this box and calculate distance
                        for test_id in bin2id[t][i][j][k]:
                            if (test_id in id2neighbors[t][atom]) or (test_id == atom) or (id2type[test_id] == id2type[atom]):
                                continue       
                            dr = r[t][test_id] + shift - r[t][atom]
                            dr2 = dr.dot(dr)
                            if dr2 < dist_range**2:
                                if nearest:
                                    if dr2 < dmin:
                                        dmin = dr2
                                        dmin_id = test_id
                                else:
                                    id2neighbors[t][atom].append(test_id)
                                    id2neighbors[t][test_id].append(atom)

            if nearest and dmin_id != 0:
                id2neighbors[t][atom].append(dmin_id)

    return id2neighbors

# Input: r, box_bounds
#  r: unscaled (but wrapped) coordinates 
#  boxsize: box dimensions (x/y/z)
# Output: bin2id, id2bin
# NOTE: does not account for changes in box size
#
def binbox(r, boxsize, bound_lo, dist_range):
    print 'Binning box...'
    bins = np.floor(boxsize/dist_range).astype(int)
    id2bin = np.zeros([len(r), len(r[0]), 3], np.int)
    bin2id = [[[[[] for k in range(bins[2])] for j in range(bins[1])] for i in range(bins[0])] for t in range(len(r))]

    # loop over frames
    for t in range(len(r)):
        # loop over atoms
        for atom in range(1, len(r[0])):
            xbin, ybin, zbin = np.floor((r[t][atom] - bound_lo)/boxsize*bins).astype(int)
            try:
                id2bin[t][atom] = (xbin, ybin, zbin)
                bin2id[t][xbin][ybin][zbin].append(atom)
            except:
                print >> sys.stderr, 'Warning: atom {} at timestep {} slightly outside of box: '.format(atom, t)
                print >> sys.stderr, '{} > {}'.format(r[t][atom], bound_hi)

    return bin2id, id2bin, bins

# NOTE: Based on https://lammps.sandia.gov/threads/msg32219.html, atoms can be slightly outside the periodic boundary
# if such thing happens, we wrap it back
def wrap(r, box_bounds):
    print 'Wrapping...'
    bound_lo = np.array([ box_bounds[0][0][0], box_bounds[0][1][0], box_bounds[0][2][0] ])
    bound_hi = np.array([ box_bounds[0][0][1], box_bounds[0][1][1], box_bounds[0][2][1] ])
    boxsize = bound_hi-bound_lo
    for t in range(len(r)):
        for atom in range(1, len(r[0])):
            shift = np.zeros(3, np.float)
            for axis, coord in enumerate(r[t][atom]):
                if coord >= bound_hi[axis]:
                    shift[axis] = -boxsize[axis]
                elif coord < bound_lo[axis]:
                    shift[axis] = boxsize[axis]
            if np.sum(np.abs(shift)) != 0:
                r[t][atom] += shift
    return r, boxsize, bound_lo

def ipcorr(r, id2neighbors):
    # print 'Calculating ion pair correlation function...'
    h = np.zeros([len(r), len(r[0])], np.float)
    corr = np.zeros([len(r), len(r[0])], np.float)
    for t in range(len(r)):
        for atom in range(1, len(r[0])):
            # check if the list is empty (see if there is ion pair)
            if id2neighbors[t][atom]:
                h[t][atom] = 1
                if t == 0:
                    corr[t][atom] = 1
                elif id2neighbors[t][atom] == id2neighbors[0][atom]:
                    corr[t][atom] = 1
                elif not id2neighbors[0][atom]: #FIXME: a convenient but dangerous way to do this..
                    for atom_id in id2neighbors[t][atom]:
                        id2neighbors[0][atom].append(atom_id)
                    corr[t][atom] = 1
                    print 'new neighbors:' t, atom, id2neighbors[t][atom]
            # else:
            #     print 'no neighbor:', t, atom, id2neighbors[t][atom]
    return corr.mean(axis=1)/h.mean()

if __name__ == '__main__':    

    r_1stshell = 1.45
    intrvl = 38
    nBlock = 76
    blockSize = 39
    timescale = 0.005

    r, ir, timestep, box_bounds, id2type, id2index = special_read('loglog.lammpstrj', types=[3,4], num_frames=((nBlock-1)*intrvl + blockSize))
    r, boxsize, bound_lo = wrap(r, box_bounds)
    bin2id, id2bin, bins = binbox(r, boxsize, bound_lo, dist_range=r_1stshell)
    
    if not (nBlock-1)*intrvl + blockSize == len(r):
        print >> sys.stderr, 'Warning: not using all data in the dump file.'

    #### IC autocorrelation calculation with 1st shell distance
    id2neighbors = buildnlist(r, bin2id, id2bin, bins, boxsize, id2type, dist_range=r_1stshell, nearest=False)
    iccorrarr = np.zeros([nBlock, blockSize], np.float)
    for n in range(nBlock):
        lo = n*intrvl
        hi = n*intrvl + blockSize
        print 'Calculating ion cage correlation in block' , n, 'with timesteps between:', timestep[lo], timestep[hi-1]
        iccorrarr[n] = ipcorr(r[lo:hi], id2neighbors[lo:hi])
    iccorrarr = iccorrarr[~np.isnan(iccorrarr).all(axis=1)]
    
    pd.DataFrame(zip(timestep[:blockSize]*timescale, iccorrarr.mean(axis=0), iccorrarr.std(axis=0))).to_excel('iccorr_{}blockavg_1stshell.xlsx'.format(nBlock), header=['Time', 'Avg', 'Std'], index=False)

    # #### IP autocorrelation calculation with 5sigma
    id2neighbors = buildnlist(r, bin2id, id2bin, bins, boxsize, id2type, dist_range=r_1stshell, nearest=True)
    ipcorrarr = np.zeros([nBlock, blockSize], np.float)
    for n in range(nBlock):
        lo = n*intrvl
        hi = n*intrvl + blockSize
        print 'Calculating ion pair correlation in block' , n, 'with timesteps between:', timestep[lo], timestep[hi-1]
        ipcorrarr[n] = ipcorr(r[lo:hi], id2neighbors[lo:hi])
    ipcorrarr = ipcorrarr[~np.isnan(ipcorrarr).all(axis=1)]
    
    pd.DataFrame(zip(timestep[:blockSize]*timescale, ipcorrarr.mean(axis=0), ipcorrarr.std(axis=0))).to_excel('ipcorr_{}blockavg_r_1stshell.xlsx'.format(nBlock), header=['Time', 'Avg', 'Std'], index=False)

    