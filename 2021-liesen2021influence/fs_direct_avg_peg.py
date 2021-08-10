#!/usr/bin/python2

#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
 
# Script:  Fs_direct_avg.py
# Purpose: calculate self-intermediate scattering function of polymers or ions
# Syntax:  fs_direct_avg_peg.py 
# Author:  Kevin Shen Feb. 2018
# Modified: Nicholas Liesen, Spring 2021 (minor modifications to calculate
# hard and soft segment mobility in PEG/DD-ionene system)


import sys
import numpy as np


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


"""
The first part of the script deals solely with reading in the coordinates from the trajectory file.
This part of the script can probably be replaced with Jon's PPPMD function read_lammpstrj
""" 

fname = 'PEO2k_75pcent_log_dump.lammpstrj'
f = open(fname, 'r')
navg = 3  # averaging for this many different starting times
frameperblock = 43  # block size
nskip = 1  # skip this many frames between windows/blocks, added by NTL  
#num_frames = navg*frameperblock+1
num_frames = navg*frameperblock + (navg-1)*nskip  # changed by NTL, skipping frames between blocks
length_scale = [0.9, 1]
timescale = 0.005

# read in the initial header
frame = 0
init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = read_header(f)

# it is not possible to preallocate arrays (or to know number of frames in advance)

alloc = num_frames
inf_frames = False

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
if "ix" in line:
	ix_index = line.index("ix") - 2
	iy_index = line.index("iy") - 2
	iz_index = line.index("iz") - 2
else:
	print >> sys.stderr, "ERROR: x image flag not found in lammps trajectory"

num_type1 = 0
num_type2 = 0
num_ions = 0

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

	# if available, build the i2mol and id2type arrays
	if mol_index is not None:
		id2mol[my_id] = int(line[mol_index])
	if type_index is not None:
		id2type[my_id] = int(line[type_index])
		if id2type[my_id] == 1:
			num_type1 += 1
		elif id2type[my_id] == 2:
			num_type2 += 1
		else:
			num_ions += 1

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
frame = 1
print "reading input file..."
while frame < num_frames:
	print frame
	# try to read in a new header
	try:
		my_timestep, my_num_atoms, my_xlo, my_xhi, my_ylo, my_yhi, my_zlo, my_zhi = read_header(f)
	except:
		print >> sys.stderr, "WARNING: hit end of file when reading in", fname, "at frame", frame
		break

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


""" 
Finished reading .lammpstrj trajectory file. At this point the script switches over to the actual calculation of the self-intermediate
scattering function -- NTL
"""

frames = len(r)
atoms = len(id2type)

print "Reading input file done: now calculate self-intermediate scattering function!"

# loop over time
for t in range(frames):
	print "getting needed coordinates for all polymers for frame %i" % t
	box_size = np.array([ box_bounds[t][0][1] - box_bounds[t][0][0], box_bounds[t][1][1] - box_bounds[t][1][0], box_bounds[t][2][1] - box_bounds[t][2][0] ])
	for k in xrange(atoms):  # unwrapping coordinates
		r[t][k] += ir[t][k]*box_size


print "averaging for different starting times..."
# loop over time and take dot products
print "calculating for the first k vector"

num_peg_beads = num_type1
num_hard_beads = num_type2 + num_ions
print "Number PEG beads "+str(num_peg_beads)
print "Number Hard beads"+str(num_hard_beads)
for length in length_scale:
	# preallocate bond vector array and autocorr array
	fs_peg = np.zeros((navg, frameperblock))
        fs_hard = np.zeros((navg, frameperblock))
	#fs_PE = np.zeros((navg, frameperblock))
	#fs_ions = np.zeros((navg, frameperblock))

	for t in range(0, navg):  # different time origins
		print t
		#t1 = t*frameperblock
                t1 = t*frameperblock + t*nskip  # because I'm skipping frames between blocks -- NTL
		for t2 in range(t1, t1+frameperblock):
			dt = t2 - t1
			dr_peg = np.array(r[t2][id2type==1])-np.array(r[t1][id2type==1])  # type 1 is PEG -- NTL
                        dr_hard = np.array(r[t2][(id2type==2)|(id2type==3)|(id2type==4)])-np.array(r[t1][(id2type==2)|(id2type==3)|(id2type==4)])  # type 2 is PE -- NTL
			#dr_ions = np.array(r[t2][(id2type==3)|(id2type==4)])-np.array(r[t1][(id2type==3)|(id2type==4)])  # type 3 and 4 are Br- and C2N+ -- NTL
			#dr_PE = np.array(r[t2][id2type==2])-np.array(r[t1][id2type==2])

			# Calculation of Fs summation -- NTL
			fs_peg[t][dt] += sum (np.exp( -1.j * np.inner( [0,0,2*np.pi/length], dr_peg))).real/num_peg_beads/3  # calculate sum( e^[ -i [ kz( riz(t2) - riz(t1) ) ] ] ) / 3N  :: NTL 
			fs_peg[t][dt] += sum (np.exp( -1.j * np.inner( [0,2*np.pi/length,0], dr_peg))).real/num_peg_beads/3  # calculate sum( e^[ -i [ ky( riy(t2) - riy(t1) ) ] ] ) / 3N  :: NTL
			fs_peg[t][dt] += sum (np.exp( -1.j * np.inner( [2*np.pi/length,0,0], dr_peg))).real/num_peg_beads/3  # calculate sum( e^[ -i [ kx( rix(t2) - rix(t1) ) ] ] ) / 3N  :: NTL

			fs_hard[t][dt] += sum (np.exp( -1.j * np.inner( [0,0,2*np.pi/length], dr_hard))).real/num_hard_beads/3
			fs_hard[t][dt] += sum (np.exp( -1.j * np.inner( [0,2*np.pi/length,0], dr_hard))).real/num_hard_beads/3
			fs_hard[t][dt] += sum (np.exp( -1.j * np.inner( [2*np.pi/length,0,0], dr_hard))).real/num_hard_beads/3

			#fs_ions[t][dt] += sum (np.exp( -1.j * np.inner( [0,0,2*np.pi/length], dr_ions))).real/num_ions/3
			#fs_ions[t][dt] += sum (np.exp( -1.j * np.inner( [0,2*np.pi/length,0], dr_ions))).real/num_ions/3
			#fs_ions[t][dt] += sum (np.exp( -1.j * np.inner( [2*np.pi/length,0,0], dr_ions))).real/num_ions/3

			#fs_PE[t][dt] += sum (np.exp( -1.j * np.inner( [0,0,2*np.pi/length], dr_PE))).real/num_type2/3
			#fs_PE[t][dt] += sum (np.exp( -1.j * np.inner( [0,2*np.pi/length,0], dr_PE))).real/num_type2/3
			#fs_PE[t][dt] += sum (np.exp( -1.j * np.inner( [2*np.pi/length,0,0], dr_PE))).real/num_type2/3
	
	print fs_peg
	print fs_hard
	file = 'fs_of_kt_{}sigma_peg_beads.dat'.format(length)
	OUT = open(file, 'w')
	OUT.write("timesteps avg std\n")
	for dt in range(frameperblock):  # iterate through all dt values -- NTL
		OUT.write("%7.3f " % (float(timestep[dt]-timestep[0])*timescale))  # output dt value in units of simulation -- NTL
		OUT.write("%7f " % fs_peg.mean(axis=0)[dt])  # average over all blocks/origins -- NTL
		OUT.write("%7f\n" % fs_peg.std(axis=0)[dt])
		#OUT.write("%7f\n" % fs_final2[t])
	OUT.close()

	file = 'fs_of_kt_{}sigma_hard_beads.dat'.format(length)
	OUT = open(file, 'w')
	OUT.write("timesteps avg std\n")
	for dt in range(frameperblock):
		OUT.write("%7.3f " % (float(timestep[dt]-timestep[0])*timescale))
		OUT.write("%7f " % fs_hard.mean(axis=0)[dt])
		OUT.write("%7f\n" % fs_hard.std(axis=0)[dt])
	OUT.close()

	#file = 'fs_of_kt_{}sigma_ions.dat'.format(length)
	#OUT = open(file, 'w')
	#OUT.write("timesteps avg std\n")
	#for dt in range(frameperblock):
	#	OUT.write("%7.3f " % (float(timestep[dt]-timestep[0])*timescale))
	#	OUT.write("%7f " % fs_ions.mean(axis=0)[dt])
	#	OUT.write("%7f\n" % fs_ions.std(axis=0)[dt])
	#OUT.close()

	#file = 'fs_of_kt_{}sigma_PE.dat'.format(length)
	#OUT = open(file, 'w')
	#OUT.write("timesteps avg std\n")
	#for dt in range(frameperblock):
	#	OUT.write("%7.3f " % (float(timestep[dt]-timestep[0])*timescale))
	#	OUT.write("%7f " % fs_PE.mean(axis=0)[dt])
	#	OUT.write("%7f\n" % fs_PE.std(axis=0)[dt])
	#OUT.close()
