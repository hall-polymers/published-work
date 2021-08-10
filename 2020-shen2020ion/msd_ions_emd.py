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
 
# Script:  msd_ions_EMD.py
# Purpose: calculate mean squared displacement of various types
# Syntax:  msd_ions_EMD.py 
# Author:  Kevin Shen Oct. 2017
 
# derived from Jon's pppmd.py
# -------------------------------------------------------

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
		num_mols = id2mol.max()	
		mol2ids = [[]]
		for molid in range(1, num_mols+1):
			mol2ids.append(np.where(id2mol==molid)[0])
	else:
		num_mols = None
		mol2ids = None

	# loop over number of num_frames frames, if num_frames is infinite, will loop over all the frames in the file
	frame = 1 # this is the frame counter for frames actually read in
	frame_attempt = 0 # this is the actual frame count in the file (not counting the ones skipped in the beginning
	while frame < num_frames:
		print "reading frame", frame
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
def MSD(r, ir, box_bounds, id2type=[]):
	# set up some constants
	frames = len(r)
	box_size = np.array([ box_bounds[0][0][1] - box_bounds[0][0][0], box_bounds[0][1][1] - box_bounds[0][1][0], box_bounds[0][2][1] - box_bounds[0][2][0] ])

	#  allocate an array for the box center of mass which needs to be subtracted off
	box_com = np.zeros([frames,3], np.float)

	# preallocate msd vectors
	msd_dict = {}
	for type_id in set(id2type):
		msd_dict[type_id] = np.zeros(frames, np.float)

	cation_com = np.zeros([frames,3], np.float)
	anion_com = np.zeros([frames,3], np.float)
	cond = np.zeros(frames, np.float)

	# loop over frames
	for t in xrange(frames):
		# calculate the center of mass of the entire box
		for atom in range(1, len(r[0])):
			box_com[t] += r[t][atom] + ir[t][atom]*box_size
		box_com[t] = box_com[t]/(len(r[0])-1)

		# if t == 0:
		# 	for atom in [i for i, x in enumerate(id2type) if x == 3]:
		# 		cation_com_0 += r[0][atom] + ir[0][atom]*box_size - box_com[0]

		# 	for atom in [i for i, x in enumerate(id2type) if x == 4]:
		# 		anion_com_0 += r[0][atom] + ir[0][atom]*box_size - box_com[0]

		# loop over atoms
		for atom in range(1, len(id2type)):
			# calculate how much the bead has moved reletive to the center of mass (note that this is a vector equation)
			r_t = r[t][atom] + ir[t][atom]*box_size - box_com[t]
			r_0 = r[0][atom] + ir[0][atom]*box_size - box_com[0]
			diff = r_t - r_0
			# the mean squared displacement is this difference dotted with itself
			msd_dict[id2type[atom]][t] += diff.dot(diff)

			if id2type[atom] == 3:
				cation_com[t] += r_t
			elif id2type[atom] == 4:
				anion_com[t] += r_t

		cond_x = ((cation_com[t][0] - cation_com[0][0]) - (anion_com[t][0] - anion_com[0][0]))**2
		cond_y = ((cation_com[t][1] - cation_com[0][1]) - (anion_com[t][1] - anion_com[0][1]))**2
		cond_z = ((cation_com[t][2] - cation_com[0][2]) - (anion_com[t][2] - anion_com[0][2]))**2
		cond[t] += (cond_x + cond_y + cond_z)

	# scale MSD by the number of beads of each type, to get the average MSD
	for type_id in set(id2type):
		msd_dict[type_id] = msd_dict[type_id]/sum(id2type == type_id)
	del msd_dict[0] # this is needed since id2type has a dummy entry of 0 at index 0 so that it is indexed by LAMMPS atom_id

	cond = cond/(sum(id2type == 3)+sum(id2type == 4))
	
	return msd_dict, cond

'''
								simulation
!---------------------------------------------------------------------!
|________________blockSize___________________|
    intrvl	|_____________________________________________|
				intrvl 	|_____________________________________________|

'''

# NOTE: (nBlock-1)*intrvl + blockSize = len(r) should always be true.
# NOTE: if intrvl = blockSize, blocks are independent of each other, elif intrvl < blockSize, blocks are overlapping.
def cond_block_avg(intrvl, nBlock, blockSize, timescale):
	if not (nBlock-1)*intrvl + blockSize == len(r):
		print >> sys.stderr, "WARNING: missing some timesteps in the dump file."

	msd = {}
	msd_avg = {}
	cond = {}
	cond_avg = np.zeros(blockSize, np.float)

	for n in range(nBlock):
		low = n*intrvl
		high = n*intrvl + blockSize
		print "calculating MSD in block ", n, "with timesteps between:", timestep[low], timestep[high-1]
		[msd[n], cond[n]] = MSD(r[low:high], ir[low:high], boxbds[low:high], id2type)
		# print msd[n]


	for n in xrange(nBlock):
		for type_id in sorted(msd[n]):
			if n ==0:
				msd_avg[type_id] = np.zeros(blockSize, np.float)
			# print "msd_avg[type_id], msd[n][type_id]", msd_avg[type_id], msd[n][type_id]
			msd_avg[type_id] += msd[n][type_id]/float(nBlock)
		cond_avg += cond[n]/float(nBlock)

	totFrame = (nBlock-1)*intrvl + blockSize
	file = 'ion_msdavged%d_tot%d_xyz.txt'% (nBlock, totFrame)
	OUT = open(file, 'w')
	OUT.write("time")
	# for n in xrange (nBlock):
	# 	OUT.write("	block%i" % n)
	OUT.write("	avg stdev\n")

	for t in xrange(blockSize):
		OUT.write("%8.5f" % (float(timestep[t+low]-timestep[low])*timescale))
		msd_std = []
		for n in range(nBlock):
			msd_std.append((msd[n][3][t]+msd[n][4][t])/2)
			# OUT.write("	%8.5f" % (msd[n][3][t]+msd[n][4][t]))
		OUT.write("	%8.5f %8.5f\n" % ((msd_avg[3][t]+msd_avg[4][t])/2, np.std(msd_std)))
	OUT.close()

	file = 'ion_condavged%d_tot%d_xyz.txt'% (nBlock, totFrame)
	OUT = open(file, 'w')
	OUT.write("time")
	# for n in xrange (nBlock):
	# 	OUT.write("	block%i" % n)
	OUT.write("	avg\n")

	for t in xrange(blockSize):
		OUT.write("%8.5f" % (float(timestep[t+low]-timestep[low])*timescale))
		cond_std = []
		for n in range(nBlock):
			cond_std.append(cond[n][t])
			# OUT.write("	%8.5f" % cond[n][t])
		OUT.write("	%8.5f %8.5f\n" % (cond_avg[t], np.std(cond_std)))
	OUT.close()


init_skip = 0 
r, ir, timestep, boxbds, id2type, id2mol, mol2ids = read_lammpstrj('nvt.lammpstrj', skip_beginning=init_skip)

# conduction, 600000 tau
cond_block_avg(10, 600, 11, 0.005)
cond_block_avg(20, 300, 21, 0.005)
cond_block_avg(30, 200, 31, 0.005)
cond_block_avg(40, 150, 41, 0.005)
cond_block_avg(50, 120, 51, 0.005)
cond_block_avg(60, 100, 61, 0.005)
cond_block_avg(100, 60, 101, 0.005)
cond_block_avg(120, 50, 121, 0.005)

# conduction, 500000 tau
cond_block_avg(10, 500, 11, 0.005)
cond_block_avg(20, 250, 21, 0.005)
cond_block_avg(40, 125, 41, 0.005)
cond_block_avg(50, 100, 51, 0.005)
cond_block_avg(100, 50, 101, 0.005)

# conduction, 400000 tau
cond_block_avg(10, 400, 11, 0.005)
cond_block_avg(20, 200, 21, 0.005)
cond_block_avg(40, 100, 41, 0.005)
cond_block_avg(50, 80, 51, 0.005)
cond_block_avg(100, 40, 101, 0.005)

# conduction, 300000 tau
cond_block_avg(10, 300, 11, 0.005)
cond_block_avg(15, 200, 16, 0.005)
cond_block_avg(20, 150, 21, 0.005)
cond_block_avg(30, 100, 31, 0.005)
cond_block_avg(40, 75, 41, 0.005)
cond_block_avg(50, 60, 51, 0.005)
cond_block_avg(60, 50, 61, 0.005)
cond_block_avg(100, 30, 101, 0.005)

# conduction, 200000 tau
cond_block_avg(10, 200, 11, 0.005)
cond_block_avg(20, 100, 21, 0.005)
cond_block_avg(40, 50, 41, 0.005)
cond_block_avg(50, 40, 51, 0.005)
cond_block_avg(100, 20, 101, 0.005)

# conduction, 100000 tau
cond_block_avg(10, 100, 11, 0.005)
cond_block_avg(20, 50, 21, 0.005)
cond_block_avg(40, 25, 41, 0.005)
cond_block_avg(50, 20, 51, 0.005)
cond_block_avg(100, 10, 101, 0.005)

# diffusion, 600000 tau
cond_block_avg(500, 12, 501, 0.005)
cond_block_avg(1000, 6, 1001, 0.005)
cond_block_avg(1500, 4, 1501, 0.005)
cond_block_avg(2000, 3, 2001, 0.005)

# diffusion, 500000 tau
cond_block_avg(500, 10, 501, 0.005)
cond_block_avg(1000, 5, 1001, 0.005)

# diffusion, 400000 tau
cond_block_avg(500, 8, 501, 0.005)
cond_block_avg(1000, 4, 1001, 0.005)

# diffusion, 300000 tau
cond_block_avg(500, 6, 501, 0.005)
cond_block_avg(750, 4, 751, 0.005)
cond_block_avg(1000, 3, 1001, 0.005)
cond_block_avg(1500, 2, 1501, 0.005)

# diffusion, 200000 tau
cond_block_avg(500, 4, 501, 0.005)
cond_block_avg(1000, 2, 1001, 0.005)
cond_block_avg(2000, 1, 2001, 0.005)

# diffusion, 100000 tau
cond_block_avg(500, 2, 501, 0.005)
cond_block_avg(1000, 1, 1001, 0.005)

