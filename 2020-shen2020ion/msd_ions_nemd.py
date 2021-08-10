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
 
# Script:  msd_ions_NEMD.py
# Purpose: calculate mean squared displacement of various types in various directions (parallel or perpendicular to the electric field)
# Syntax:  msd_ions_NEMD.py 
# Author:  Kevin Shen Oct. 2018
 
# derived from Jon's pppmd.py
# -------------------------------------------------------

import sys
import numpy as np

# special_read: read in a lammps trajectory, but only extract data of certain bead types
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
def special_read(fname, types, com_ref='whole', num_frames=float('inf'), skip_beginning=0, skip_between=0):
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
    box_bounds_lo = np.zeros([alloc,3], np.float) # 3D array to store lower boundaries of the box, indexed by frame, x/y/z
    box_bounds_hi = np.zeros([alloc,3], np.float) # 3D array to store upper boundaries of the box, indexed by frame, x/y/z
    box_size = np.zeros([alloc,3], np.float) # 3D array to store sizes of the box, indexed by frame, x/y/z

    timestep[frame] = init_timestep
    box_bounds_lo[frame][0] = xlo
    box_bounds_lo[frame][1] = ylo
    box_bounds_lo[frame][2] = zlo
    box_bounds_hi[frame][0] = xhi
    box_bounds_hi[frame][1] = yhi
    box_bounds_hi[frame][2] = zhi

    box_size[frame] = box_bounds_hi[frame] - box_bounds_lo[frame]
    
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
            r_temp[my_id] = r_temp[my_id]*box_size[frame] + box_bounds_lo[frame]

        # x, y, z image flags
        ir_temp[my_id][0] = int(line[ix_index])
        ir_temp[my_id][1] = int(line[iy_index])
        ir_temp[my_id][2] = int(line[iz_index])

    # NOTE: using num_types+1 here so that the arrays are indexed by their LAMMPS atom id
    num_not_types = num_atoms - num_types
    r = np.zeros([alloc, num_types+1, 3], np.float) # 3D array of x, y, z coordinates, r[frame][id][coordinate]
    ir = np.zeros([alloc, num_types+1, 3], np.int) # 3D array of x, y, z image flags, ir[frame][id][coordinate]
    box_com = np.zeros([alloc, 3], np.float) # 3D array of x, y, z box center of mass, box_com[frame][coordinate]

    id2type = np.zeros(num_types+1, np.int) # array to map from atom id to type, builds this from the first frame, if available
    id2index =  np.zeros(num_types+1, np.int) # array to map from atom id to index, builds this from the first frame, if available

    # store the temporary data into real arrays
    my_id = 0
    for atom in range(num_atoms):
        index = atom+1
        if com_ref == 'whole':
            box_com[frame] += (r_temp[index] + ir_temp[index]*box_size[frame])/num_atoms

        if index2type[index] in types:
            my_id += 1
            # x, y, z coordinates
            r[frame][my_id] = r_temp[index]

            # x, y, z image flags
            ir[frame][my_id] = ir_temp[index]

            id2type[my_id] = index2type[index]
            id2index[my_id] = index
            index2id[index] = my_id

        elif com_ref == 'polymer':
            box_com[frame] += (r_temp[index] + ir_temp[index]*box_size[frame])/num_not_types

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
            
            box_bounds_lo = np.concatenate( ( box_bounds_lo, np.zeros([1,3],np.float) ) )
            box_bounds_hi = np.concatenate( ( box_bounds_hi, np.zeros([1,3],np.float) ) )
            box_size = np.concatenate( ( box_size, np.zeros([1,3],np.float) ) )

            r = np.concatenate( ( r, np.zeros([1, num_types+1, 3], np.float) ) )
            ir = np.concatenate( ( ir, np.zeros([1, num_types+1, 3], np.float) ) )
            box_com = np.concatenate( ( box_com, np.zeros([1, 3], np.float) ) )
        
        # update the timestep and box size arrays
        timestep[frame] = my_timestep
        box_bounds_lo[frame][0] = xlo
        box_bounds_lo[frame][1] = ylo
        box_bounds_lo[frame][2] = zlo
        box_bounds_hi[frame][0] = xhi
        box_bounds_hi[frame][1] = yhi
        box_bounds_hi[frame][2] = zhi
        box_size[frame] = box_bounds_hi[frame] - box_bounds_lo[frame]

        f.readline() # ITEM: ATOMS
        # loop over the atoms lines
        r_temp = np.zeros([num_atoms+1, 3], np.float) # 3D array of x, y, z coordinates, r[frame][id][coordinate]
        ir_temp = np.zeros([num_atoms+1, 3], np.int) # 3D array of x, y, z image flags, r[frame][id][coordinate]
        for atom in range(num_atoms):
            line = f.readline()
            line = line.split()
    
            # get the atom id
            index = int(line[id_index])

            # x, y, z coordinates
            r_temp[index][0] = float(line[x_index])
            r_temp[index][1] = float(line[y_index])
            r_temp[index][2] = float(line[z_index])

            # unscale, if necessary
            if scaled:
                r_temp[index] = r_temp[index]*box_size[frame] + box_bounds_lo[frame]

            # x, y, z image flags
            ir_temp[index][0] = int(line[ix_index])
            ir_temp[index][1] = int(line[iy_index])
            ir_temp[index][2] = int(line[iz_index])

            if com_ref == 'whole':
                box_com[frame] += (r_temp[index] + ir_temp[index]*box_size[frame])/num_atoms

            if index2type[index] in types:
                my_id = index2id[index]

                # x, y, z coordinates
                r[frame][my_id] = r_temp[index]

                # x, y, z image flags
                ir[frame][my_id] = ir_temp[index]

            elif com_ref == 'polymer':
                box_com[frame] += (r_temp[index] + ir_temp[index]*box_size[frame])/num_not_types

        print "Reading frame {}".format(frame)
        frame += 1
    
    print 'Summary:'
    print 'Total number of atoms =', num_atoms
    print 'Total number of ions =', num_types

    return r, ir, timestep, box_size, id2type, id2index, box_com

# MSD: mean squared displacement
#
# Input: r, ir, box_bounds, id2type
# computes separate MSDs for each type as given in id2type, msd_dict["type"] is a dict of these MSDs
#  r: unscaled (but wrapped) coordinates 
#  ir: image flags
#  box_bounds: boundaries of the box
#  id2type: array that maps atom id to type 
#  (format as read in from read_lammpstrj)
#
# Output: msd_dict filled with the calculated MSDs 
#  each entry is an array indexed by frame gives average MSD of all beads of the different types
#
# NOTE: assumes mass = 1 for all beads
# NOTE: does not account for changes in box size
#
def MSD(r, ir, box_size, box_com, id2type, nofield_dir):
    # set up some constants
    frames = len(r)

    # preallocate msd vectors
    dr_dict_efield = {}
    msd_dict_efield = {}
    msd_dict_noefield ={}
    for type_id in set(id2type):
        dr_dict_efield[type_id] = np.zeros(frames, np.float)
        msd_dict_efield[type_id] = np.zeros(frames, np.float)
        msd_dict_noefield[type_id] = np.zeros(frames, np.float)

    # loop over frames
    for t in xrange(frames):
        # loop over atoms
        for atom in range(1, len(id2type)):
            # calculate how much the bead has moved reletive to the center of mass (note that this is a vector equation)
            r_t = r[t][atom] + ir[t][atom]*box_size[t] - box_com[t]
            r_0 = r[0][atom] + ir[0][atom]*box_size[0] - box_com[0]
            diff = r_t - r_0

            dr_dict_efield[id2type[atom]][t] += diff[0]
            # the mean squared displacement is this difference dotted with itself
            msd_dict_efield[id2type[atom]][t] += diff[0]*diff[0]
            if nofield_dir == 'yz':
                msd_dict_noefield[id2type[atom]][t] += (diff[1]*diff[1]+diff[2]*diff[2])
            elif nofield_dir == 'y':
                msd_dict_noefield[id2type[atom]][t] += diff[1]*diff[1]

    # scale MSD by the number of beads of each type, to get the average MSD
    for type_id in set(id2type):
        dr_dict_efield[type_id] = dr_dict_efield[type_id]/sum(id2type == type_id)
        msd_dict_efield[type_id] = msd_dict_efield[type_id]/sum(id2type == type_id)
        msd_dict_noefield[type_id] = msd_dict_noefield[type_id]/sum(id2type == type_id)

    del dr_dict_efield[0] # this is needed since id2type has a dummy entry of 0 at index 0 so that it is indexed by LAMMPS atom_id
    del msd_dict_efield[0]
    del msd_dict_noefield[0]

    return dr_dict_efield, msd_dict_efield, msd_dict_noefield

'''
                                simulation
!---------------------------------------------------------------------!
|________________blockSize___________________|
    intrvl    |_____________________________________________|
                intrvl     |_____________________________________________|

'''

# NOTE: (nBlock-1)*intrvl + blockSize = len(r) should always be true.
# NOTE: if intrvl = blockSize, blocks are independent of each other, elif intrvl < blockSize, blocks are overlapping.
def msd_block_avg(data, intrvl, nBlock, blockSize, timescale, comRef, nofield_dir):

    r, ir, timestep, box_size, id2type, id2index, box_com = data
    if not (nBlock-1)*intrvl + blockSize == len(r):
        print >> sys.stderr, "WARNING: missing some timesteps in the dump file."

    if nofield_dir != 'y' and nofield_dir != 'yz':
        print >> sys.stderr, "WARNING: Please sepcify the no field direction: y or yz."

    dr_efield = {}
    msd_efield = {}
    msd_noefield = {}
    dr_efield_avg = {}
    msd_efield_avg = {}
    msd_noefield_avg = {}

    for n in range(nBlock):
        low = n*intrvl
        high = n*intrvl + blockSize
        print "calculating MSD in block ", n, "with timesteps between:", timestep[low], timestep[high-1]
        dr_efield[n], msd_efield[n], msd_noefield[n] = MSD(r[low:high], ir[low:high], box_size[low:high], box_com[low:high], id2type, nofield_dir)

    for n in range(nBlock):
        for type_id in sorted(msd_efield[n]):
            if n ==0:
                dr_efield_avg[type_id] = np.zeros(blockSize, np.float)
                msd_efield_avg[type_id] = np.zeros(blockSize, np.float)
                msd_noefield_avg[type_id] = np.zeros(blockSize, np.float)
            dr_efield_avg[type_id] += dr_efield[n][type_id]/float(nBlock)
            msd_efield_avg[type_id] += msd_efield[n][type_id]/float(nBlock)
            msd_noefield_avg[type_id] += msd_noefield[n][type_id]/float(nBlock)

    totFrame = (nBlock-1)*intrvl + blockSize
    file = 'ion_%s_msdavged%d_tot%d_x.txt'% (comRef, nBlock, totFrame)
    OUT = open(file, 'w')
    OUT.write("time")
    for n in range (nBlock):
        OUT.write(" blk%i_cat" % n)
        OUT.write(" blk%i_an" % n)
    OUT.write(" avg_cat avg_an dr_cat dr_an\n")

    for t in range(blockSize):
        OUT.write("%8.5f" % (float(timestep[t+low]-timestep[low])*timescale))
        for n in range(nBlock):
            OUT.write(" %8.5f %8.5f" % (msd_efield[n][3][t], msd_efield[n][4][t]))
        OUT.write(" %8.5f %8.5f" % (msd_efield_avg[3][t], msd_efield_avg[4][t]))
        OUT.write(" %8.5f %8.5f\n" % (dr_efield_avg[3][t], dr_efield_avg[4][t]))
    OUT.close()

    file = 'ion_%s_msdavged%d_tot%d_%s.txt'% (comRef, nBlock, totFrame, nofield_dir)
        
    OUT = open(file, 'w') 
    OUT.write("time")
    for n in range (nBlock):
        OUT.write(" blk%i_cat" % n)
        OUT.write(" blk%i_an" % n)
    OUT.write(" avg_cat avg_an\n")

    for t in range(blockSize):
        OUT.write("%8.5f" % (float(timestep[t+low]-timestep[low])*timescale))
        for n in range(nBlock):
            OUT.write(" %8.5f %8.5f" % (msd_noefield[n][3][t], msd_noefield[n][4][t]))
        OUT.write(" %8.5f %8.5f\n" % (msd_noefield_avg[3][t], msd_noefield_avg[4][t]))
    OUT.close()

if __name__ == '__main__':
    
    init_skip = 0 
    com_ref = ['whole']   #option: use 'polymer' to return only the center of mass of polymers instead that of the whole system
    msd_nofield_dir = 'yz'   # options: use 'y' or 'yz' to account for msd in the direction(s) with no electric field; msd in the x direction (efield direction) will be calculated regardless of what this option is
    totalframe = 6001
    nblock_list = [1, 2, 3, 4, 5, 6, 8, 10 ,12, 15, 20, 24, 25, 30]
    for comRef in com_ref:
        data = special_read('nvt.lammpstrj', [3, 4], com_ref=comRef, skip_beginning=init_skip) 
        for nb in nblock_list:
            if (totalframe-1)%nb == 0:
                msd_block_avg(data, (totalframe-1)/nb, nb, (totalframe-1)/nb+1, 0.005, comRef, msd_nofield_dir)

