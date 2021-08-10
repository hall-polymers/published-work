#!/usr/bin/env python

#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


# PPPMD: post-pocessing polymer molecular dynamics
# Suite of tools to do post processing calculations relevant to polymers 
#
# Modified J. Brown 2014-04-02

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
# FIXME: also assumes that the coordinates are wrapped, so x or xs type coordinates are allowed but not xu or xsu
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

	# loop over number of num_frames frames, if num_frames is infinite, will look over all the frames in the file
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
# Output: end to end vector autocorrelation function (1D array indexed by frame count)
#
# NOTE: all the listings in mol2ids will be used and averaged together
# NOTE: it is assumed that the end-to-end vector is the one between the lowest and highest id in each molecule (if this is not the case, you'd have to mess with mol2ids, e.g. make it only contain the ids of the two end beads)
# NOTE: scaled by the average end-to-end vector, so that e2e_autocorr[0]=1.0
#
def end2end_autocorr(r, ir, box_bounds, mol2ids):
	frames = len(r)

	# by default, to make the mol2ids array be indexed by molid, there is an extra "[]" element inserted at the beginning, this removes it if it's present
	if len(mol2ids[0]) == 0:
		del mol2ids[0]
	
	mols = len(mol2ids)

	# preallocate e2e vector array and autocorr array
	e2e = np.zeros([frames, mols, 3], np.float)
	e2e_autocorr = np.zeros(frames, np.float)

	# loop over time
	for t in range(frames):
		box_size = np.array([ box_bounds[t][0][1] - box_bounds[t][0][0], box_bounds[t][1][1] - box_bounds[t][1][0], box_bounds[t][2][1] - box_bounds[t][2][0] ])

		# loop over molecules
		for mol in range(mols):
			id1 = mol2ids[mol].min()
			id2 = mol2ids[mol].max()

			r1 = r[t][id1] + ir[t][id1]*box_size
			r2 = r[t][id2] + ir[t][id2]*box_size

			e2e[t][mol] = r2 - r1

	# loop over time and take dot products
	# FIXME: set t1 to zero for now
	for t1 in range(1):
		for t2 in range(t1,frames):
			dt = t2 - t1
			for mol in range(mols):
				e2e_autocorr[dt] += np.dot(e2e[t1][mol], e2e[t2][mol])

	# scaling
	e2e_autocorr = e2e_autocorr/mols
	e2e_autocorr = e2e_autocorr/e2e_autocorr[0]

	return e2e_autocorr

# MSD_chain: mean squared displacement on chains
#
# Input: r, ir, box_bounds, mol2ids
#  r: unscaled (but wrapped) coordinates 
#  ir: image flags
#  box_bounds: boundaries of the box
#  mol2ids: list of atom id arrays of the molecules to be used
#  (format as read in from read_lammpstrj)
#
# Output: msd_ccom, msd_beads, msd_ends, msd_mid
#  each is an array indexed by frame gives average MSD of all chains in mol2ids
#  msd_ccom: chain center of mass
#  msd_beads: all the beads separately
#  msd_ends: end of the chains (assumes the ends are the max/min atom ids of the molecule)
#  msd_mid: middle of the chain (assumes the middle is the middle atom id in the molecule)
#
# NOTE: all the listings in mol2ids will be used and averaged together
# FIXME: assuming equally spaced timesteps for the moment
# FIXME: does not correctly accound to changing box size
# FIXME: brute force n^2 algorithim
#
def MSD_chain(r, ir, box_bounds, mol2ids):
	frames = len(r)

	# by default, to make the mol2ids array be indexed by molid, there is an extra "[]" element inserted at the beginning, this removes it if it's present
	if len(mol2ids[0]) == 0:
		del mol2ids[0]
	
	mols = len(mol2ids)

	# preallocate the chain center of mass array
	ccom = np.zeros([frames,mols,3], np.float)
	# preallocate msd vectors
	msd_ccom = np.zeros(frames, np.float)
	msd_beads = np.zeros(frames, np.float)
	msd_ends = np.zeros(frames, np.float)
	msd_mid = np.zeros(frames, np.float)

	# first build the chain center of mass array, and unwrap the coordinates
	# loop over frames
	for t in range(frames):
		box_size = np.array([ box_bounds[t][0][1] - box_bounds[t][0][0], box_bounds[t][1][1] - box_bounds[t][1][0], box_bounds[t][2][1] - box_bounds[t][2][0] ])
		# loop over molecules
		for mol in range(mols):
			# loop over atoms in each molecule, unscale and calc CoM
			for atom in mol2ids[mol]:
				r[t][atom] = r[t][atom] + ir[t][atom]*box_size
				ir[t][atom] = np.zeros(3,np.float)
				ccom[t][mol] += r[t][atom]
			ccom[t][mol] = ccom[t][mol]/len(mol2ids[mol])


	# loop over time twice to take advantage of all data and count up differences
	for t1 in range(frames):
		for t2 in range(t1+1,frames):
			dt = t2 - t1
			for mol in range(mols):
				# chain center of mass
				ccom_diff = ccom[t2][mol] - ccom[t1][mol]
				msd_ccom[dt] += ccom_diff.dot(ccom_diff)
				
				# end beads
				id1 = mol2ids[mol].min()
				id2 = mol2ids[mol].max()
				diff1 = r[t2][id1] - r[t1][id1]
				diff2 = r[t2][id2] - r[t1][id2]
				msd_ends[dt] += diff1.dot(diff1) + diff2.dot(diff2)

				# mid bead
				id3 = mol2ids[mol][len(mol2ids[mol])/2]
				diff3 = r[t2][id3] - r[t1][id3]
				msd_mid[dt] += diff3.dot(diff3)

				# bead by bead msds
				for atom in mol2ids[mol]:
					diff = r[t2][atom] - r[t1][atom]
					msd_beads[dt] += diff.dot(diff)
	# scaling
	msd_ccom = msd_ccom/np.arange(float(frames),0,-1)/mols
	msd_ends = msd_ends/np.arange(float(frames),0,-1)/mols/2
	msd_mid = msd_mid/np.arange(float(frames),0,-1)/mols
	msd_beads = msd_beads/np.arange(float(frames),0,-1)/sum( [len(mol) for mol in mol2ids] )

	return msd_ccom, msd_beads, msd_ends, msd_mid

# MSD_chain_xyz: mean squared displacement on chains in the x y and z components
#
# Input: r, ir, box_bounds, mol2ids
#  r: unscaled (but wrapped) coordinates 
#  ir: image flags
#  box_bounds: boundaries of the box
#  mol2ids: list of atom id arrays of the molecules to be used
#  (format as read in from read_lammpstrj)
#
# Output: msd_ccom, msd_beads, msd_ends, msd_mid
#  each is an array indexed by frame gives average MSD of all chains in mol2ids
#  msd_ccom: chain center of mass
#  msd_beads: all the beads separately
#  msd_ends: end of the chains (assumes the ends are the max/min atom ids of the molecule)
#  msd_mid: middle of the chain (assumes the middle is the middle atom id in the molecule)
#
# NOTE: all the listings in mol2ids will be used and averaged together
# FIXME: assuming equally spaced timesteps for the moment
# FIXME: does not correctly accound to changing box size
# FIXME: brute force n^2 algorithim
#
def MSD_chain_xyz(r, ir, box_bounds, mol2ids):
	frames = len(r)

	# by default, to make the mol2ids array be indexed by molid, there is an extra "[]" element inserted at the beginning, this removes it if it's present
	if len(mol2ids[0]) == 0:
		del mol2ids[0]
	
	mols = len(mol2ids)

	# preallocate the chain center of mass array
	ccom = np.zeros([frames,mols,3], np.float)
	# preallocate msd vectors
	msd_ccom = np.zeros([frames, 3], np.float)
	msd_beads = np.zeros([frames, 3], np.float)
	msd_ends = np.zeros([frames, 3], np.float)
	msd_mid = np.zeros([frames, 3], np.float)

	# first build the chain center of mass array, and unwrap the coordinates
	# loop over frames
	for t in range(frames):
		box_size = np.array([ box_bounds[t][0][1] - box_bounds[t][0][0], box_bounds[t][1][1] - box_bounds[t][1][0], box_bounds[t][2][1] - box_bounds[t][2][0] ])
		# loop over molecules
		for mol in range(mols):
			# loop over atoms in each molecule, unscale and calc CoM
			for atom in mol2ids[mol]:
				r[t][atom] = r[t][atom] + ir[t][atom]*box_size
				ir[t][atom] = np.zeros(3,np.float)
				ccom[t][mol] += r[t][atom]
			ccom[t][mol] = ccom[t][mol]/len(mol2ids[mol])


	# loop over time twice to take advantage of all data and count up differences
	for t1 in range(frames):
		for t2 in range(t1+1,frames):
			dt = t2 - t1
			for mol in range(mols):
				# chain center of mass
				ccom_diff = ccom[t2][mol] - ccom[t1][mol]
				msd_ccom[dt] += ccom_diff*ccom_diff
				
				# end beads
				id1 = mol2ids[mol].min()
				id2 = mol2ids[mol].max()
				diff1 = r[t2][id1] - r[t1][id1]
				diff2 = r[t2][id2] - r[t1][id2]
				msd_ends[dt] += diff1*diff1 + diff2*diff2

				# mid bead
				id3 = mol2ids[mol][len(mol2ids[mol])/2]
				diff3 = r[t2][id3] - r[t1][id3]
				msd_mid[dt] += diff3*diff3

				# bead by bead msds
				for atom in mol2ids[mol]:
					diff = r[t2][atom] - r[t1][atom]
					msd_beads[dt] += diff*diff
	# scaling
	msd_ccom = np.transpose( np.transpose(msd_ccom)/np.arange(float(frames),0,-1) )/mols 
	msd_ends = np.transpose( np.transpose(msd_ends)/np.arange(float(frames),0,-1) )/mols/2
	msd_mid = np.transpose( np.transpose(msd_mid)/np.arange(float(frames),0,-1) )/mols 
	msd_beads = np.transpose( np.transpose(msd_beads)/np.arange(float(frames),0,-1) )/sum( [len(mol) for mol in mol2ids] )

	return msd_ccom, msd_beads, msd_ends, msd_mid

# MSD_type: mean squared displacement by type
#
# Input: r, ir, box_bounds, id2type
#  r: unscaled (but wrapped) coordinates 
#  ir: image flags
#  box_bounds: boundaries of the box
#  id2type: array to map atom id to type
#  (format as read in from read_lammpstrj)
#
# Output: 
#  msd: types by frames array of MSDs
#
# FIXME: assuming equally spaced timesteps for the moment
# FIXME: does not correctly accound to changing box size
# FIXME: brute force n^2 algorithim
#
def MSD_type(r, ir, box_bounds, id2type):
	frames = len(r)
	types = id2type.max()
	atoms = len(id2type)

	# preallocate msd array, types+1 to start counting at 1
	msd = np.zeros([types+1, frames], np.float)

	# first build the chain center of mass array, and unwrap the coordinates
	# loop over frames
	for t in range(frames):
		box_size = np.array([ box_bounds[t][0][1] - box_bounds[t][0][0], box_bounds[t][1][1] - box_bounds[t][1][0], box_bounds[t][2][1] - box_bounds[t][2][0] ])
		# loop over atoms
		for atom in range(1,atoms):
			r[t][atom] = r[t][atom] + ir[t][atom]*box_size
			ir[t][atom] = np.zeros(3,np.float)


	# loop over time twice to take advantage of all data and count up differences
	for t1 in range(frames):
		for t2 in range(t1+1,frames):
			dt = t2 - t1
			for atom in range(1,atoms):
				diff = r[t2][atom] - r[t1][atom]
				msd[id2type[atom]][dt] += diff.dot(diff)

	count_types = np.bincount(id2type)
	for atom_type in range(1,types+1):
		msd[atom_type] = msd[atom_type]/np.arange(float(frames),0,-1)/count_types[atom_type]

	return msd

# MSD_type_xyz: mean squared displacement by type in the x, y, and z directions
#
# Input: r, ir, box_bounds, id2type
#  r: unscaled (but wrapped) coordinates 
#  ir: image flags
#  box_bounds: boundaries of the box
#  id2type: array to map atom id to type
#  (format as read in from read_lammpstrj)
#
# Output: 
#  msd: types by frames array of MSDs
#
# FIXME: assuming equally spaced timesteps for the moment
# FIXME: does not correctly accound to changing box size
# FIXME: brute force n^2 algorithim
#
def MSD_type_xyz(r, ir, box_bounds, id2type):
	frames = len(r)
	types = id2type.max()
	atoms = len(id2type)

	# preallocate msd array, types+1 to start counting at 1
	msd = np.zeros([types+1, frames, 3], np.float)

	# first build the chain center of mass array, and unwrap the coordinates
	# loop over frames
	for t in range(frames):
		box_size = np.array([ box_bounds[t][0][1] - box_bounds[t][0][0], box_bounds[t][1][1] - box_bounds[t][1][0], box_bounds[t][2][1] - box_bounds[t][2][0] ])
		# loop over atoms
		for atom in range(1,atoms):
			r[t][atom] = r[t][atom] + ir[t][atom]*box_size
			ir[t][atom] = np.zeros(3,np.float)


	# loop over time twice to take advantage of all data and count up differences
	for t1 in range(frames):
		for t2 in range(t1+1,frames):
			dt = t2 - t1
			for atom in range(1,atoms):
				diff = r[t2][atom] - r[t1][atom]
				msd[id2type[atom]][dt] += diff*diff

	count_types = np.bincount(id2type)
	for atom_type in range(1,types+1):
		msd[atom_type] = np.transpose( np.transpose(msd[atom_type])/np.arange(float(frames),0,-1) )/count_types[atom_type]

	return msd

# output_Z1: outputs the read in chains in a format appropriate for the Z1 algorithm
# see: http://polyphys-s01.ethz.ch/cgi-bin/Z1
#
# Input: r, box_bounds, mol2ids
#  r: unscaled (but wrapped) coordinates 
#  box_bounds: boundaries of the box
#  mol2ids: num_mols+1 length list of atom id arrays corresponding to the molecules (if available, may be None)
#  (format as read in from read_lammpstrj)
#
# Optional inputs:
#  filename: file to output to (defaults to std out)
#  frame: which frame to use (defaults to first)
#
# Output: 
#  none, writes the polymer locations in r from mol2ids into filename
def output_Z1(r, box_bounds, mol2ids, filename=None, frame=0):
	# if there's a filename redirect standard out there	
	if filename:
		f = open(filename,'w')
		sys.stdout = f
	
	# may have to shift so that the origin is at the center of the box
	x_shift = -(box_bounds[frame][0][1] + box_bounds[frame][0][0]) / 2
	y_shift = -(box_bounds[frame][1][1] + box_bounds[frame][1][0]) / 2
	z_shift = -(box_bounds[frame][2][1] + box_bounds[frame][2][0]) / 2

	# M
	print sum( len(mollist) > 1 for mollist in mol2ids) 

	# box_x box_y box_z
	print box_bounds[frame][0][1] - box_bounds[frame][0][0], box_bounds[frame][1][1] - box_bounds[frame][1][0], box_bounds[frame][2][1] - box_bounds[frame][2][0]

	# N1 N2 N3 .. NM
	for mol in range(1, len(mol2ids)):
		if len(mol2ids[mol]) > 1: # don't include and single site molecules in the calculation
			print len(mol2ids[mol]),
	print

	# x1,1 y1,1 z1,1
	# x1,2 y1,2 z1,2
	# ...
	# x1,N1 y1,N1 z1,N1
	# x2,1 y2,1 z2,1
	# ..
	# x2,N2 y2,N2 z2,N2
	# x3,1 y3,1 z3,1
	# ..
	# xM,NM yM,NM zM,NM
	for mol in range(1, len(mol2ids)):
		if len(mol2ids[mol]) > 1:
			for my_id in mol2ids[mol]:
				print r[frame][my_id][0] + x_shift, r[frame][my_id][1] + y_shift, r[frame][my_id][2] + z_shift
	
	# close the file when done
	if filename:
		f.close()
		sys.stdout = sys.__stdout__
