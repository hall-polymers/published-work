#!/opt/local/bin/python

#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


import sys
import numpy as np

def read_lammpstrj(fname, num_frames=float('inf')):
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

	if not fname or fname == 'stdin':
		f = sys.stdin
	else:
		f = open(fname, 'r')
	
	# read in the initial header
	frame = 0
	init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = read_header(f)

	# preallocate arrays, if possible
	if num_frames < float('inf'):
		alloc = num_frames
		inf_frames = False
	else:
		alloc = 1
		inf_frames = True
	timestep = np.zeros(alloc, np.int)
	box_bounds = np.zeros([alloc,3,2], np.float) 

	timestep[frame] = init_timestep
	box_bounds[frame][0][0] = xlo
	box_bounds[frame][0][1] = xhi
	box_bounds[frame][1][0] = ylo
	box_bounds[frame][1][1] = yhi
	box_bounds[frame][2][0] = zlo
	box_bounds[frame][2][1] = zhi
	
	r = np.zeros([alloc, num_atoms+1, 3], np.float) 
	ir = np.zeros([alloc, num_atoms+1, 3], np.int) 
	id2mol = np.zeros(num_atoms+1, np.int) 
	id2type = np.zeros(num_atoms+1, np.int) 

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

		# unscal, if necessary
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
			
	if mol_index is not None:
		num_mols = id2mol.max()	
		mol2ids = [[]]
		for molid in range(1, num_mols+1):
			mol2ids.append(np.where(id2mol==molid)[0])
	else:
		num_mols = None
		mol2ids = None

	frame = 1
	while frame < num_frames:
		try:
			my_timestep, my_num_atoms, my_xlo, my_xhi, my_ylo, my_yhi, my_zlo, my_zhi = read_header(f)
		except:
			print >> sys.stderr, "WARNING: hit end of file when reading in", fname, "at frame", frame
			break

		if inf_frames:
			timestep = np.append(timestep, 0)
			
			box_bounds = np.concatenate( ( box_bounds, np.zeros([1,3,2],np.float) ) )

			r = np.concatenate( ( r, np.zeros([1, num_atoms+1, 3], np.float) ) )
			ir = np.concatenate( ( ir, np.zeros([1, num_atoms+1, 3], np.float) ) )
		
		timestep[frame] = my_timestep
		box_bounds[frame][0][0] = my_xlo
		box_bounds[frame][0][1] = my_xhi
		box_bounds[frame][1][0] = my_ylo
		box_bounds[frame][1][1] = my_yhi
		box_bounds[frame][2][0] = my_zlo
		box_bounds[frame][2][1] = my_zhi

		f.readline() # ITEM: ATOMS
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




def end2end_autocorr(r, ir, box_bounds, mol2ids):
	frames = len(r)

	if len(mol2ids[0]) == 0:
		del mol2ids[0]
	
	mols = len(mol2ids)

	e2e = np.zeros([frames, mols, 3], np.float)
	e2e_autocorr = np.zeros(frames, np.float)

	for t in range(frames):
		box_size = np.array([ box_bounds[t][0][1] - box_bounds[t][0][0], box_bounds[t][1][1] - box_bounds[t][1][0], box_bounds[t][2][1] - box_bounds[t][2][0] ])

		for mol in range(mols):
			id1 = mol2ids[mol].min()
			id2 = mol2ids[mol].max()

			r1 = r[t][id1] + ir[t][id1]*box_size
			r2 = r[t][id2] + ir[t][id2]*box_size

			e2e[t][mol] = r2 - r1

	t1=0
	for t2 in range(t1,frames):
		dt = t2 - t1
		for mol in range(mols):
			e2e_autocorr[dt] += np.dot(e2e[t1][mol], e2e[t2][mol])

	# scaling
	e2e_autocorr = e2e_autocorr/mols
	e2e_autocorr = e2e_autocorr/e2e_autocorr[0]

	return e2e_autocorr

