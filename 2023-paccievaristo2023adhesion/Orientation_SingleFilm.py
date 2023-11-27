# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

###Purpose: This script calculates the bond orientation function as a function of interfacial distance###
### Syntax: python orientationTF.py ###
### Author: Janani Sampath, from Jon's PPPMD script###
### Date: March 2017 ###


import sys
import numpy as np
import math
npoly = 1600
#file = bond_acf.csv
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

fname = '750.lammpstrj'
f = open(fname, 'r')
file = '750_orientation_type1.csv'

# read in the initial header
frame = 0
init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = read_header(f)

# it is not possible to preallocate arrays (or to know number of frames in advance)
num_frames = float('inf')
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
if "ix" in line:
	ix_index = line.index("ix") - 2
	iy_index = line.index("iy") - 2
	iz_index = line.index("iz") - 2

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

# build the reverse of the id2mol array
# this is a 2D array with rows of (potentially) varying length, so nest a numpy array into a python list
if mol_index is not None:
	num_mols = id2mol.max()	
	mol2ids = [[]]
	for molid in range(1, num_mols+1):
		mol2ids.append(np.where(id2mol==molid)[0]) ## mol2ids :[[1,2,3,4...*nmonomers{atom id}],...*npoly {mol id}]
else:
	num_mols = None
	mol2ids = None

# loop over number of num_frames frames, if num_frames is infinite, will look over all the frames in the file
frame = 1
#print "reading input file..."
while frame < num_frames:
	#print frame
	# try to read in a new header
	try:
		my_timestep, my_num_atoms, my_xlo, my_xhi, my_ylo, my_yhi, my_zlo, my_zhi = read_header(f)
	except:
		#print >> sys.stderr, "WARNING: hit end of file when reading in", fname, "at frame", frame
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

#print "reading input file done, now calculating bonds and assigning bonds to shells around NP"

frames = len(r)
if len(mol2ids[0]) == 0:
	del mol2ids[0]


#### this sets the number of shells ####

shell_min = [] 
Count = zlo
shell_width = 1 ## sets the size of shell ##
while Count < zhi+shell_width:
        shell_min.append(Count)
        Count = Count+shell_width
## array initialization ##
n = np.zeros([len(shell_min),frames])
cosT2 = np.zeros([len(shell_min),frames,npoly,len(max(mol2ids,key=len))]) 
f2 =  np.zeros([len(shell_min),frames])
orientation =  np.zeros([len(shell_min)])

## Calculate angle between bond vector and z box, but only for backbone beads not connected to a pendant bead (type 1)##
for i in shell_min:
    for t in range(frames):
        for mol in range(npoly):
            for atomid in range(len(mol2ids[mol])-1):
                if ((id2type[mol2ids[mol][atomid]] == 1)):
                    id1 = mol2ids[mol][atomid] 
                    id3 = mol2ids[mol][atomid+1]
                    mid_point = (r[t][id1] + r[t][id3])/2
                    if (i < mid_point[2] < i+shell_width):
                        shell = shell_min.index(i)
                        n[shell][t] = n[shell][t]+1
                        vector = r[t][id1] - r[t][id3]
                        mag = math.sqrt(sum(j**2 for j in vector))       
                        BdotZ = vector[2]*(zhi-zlo)
                        cosT = BdotZ/(mag*(zhi-zlo))
                        cosT2[shell][t][mol][atomid+1] = cosT*cosT
                
## Average over shells and frames ##

for i in range(len(shell_min)):
    for t in range(frames):
        for mol in range(npoly):
            for atomid in range(len(mol2ids[mol])-1):
                f2[i][t] += cosT2[i][t][mol][atomid+1]
        f2[i][t] = f2[i][t]/n[i][t]
        orientation[i] += f2[i][t]
    orientation[i] = orientation[i]/frames
        
    
## Write output to file
        
OUT = open(file, 'w')
OUT.write("shell, orientation\n")
for s in range(1,len(shell_min)):
    OUT.write("%7i, %7f\n" % (s, ((3*orientation[s]-1)/2)))
OUT.close()


 





























            
















        

            


    

             
            
