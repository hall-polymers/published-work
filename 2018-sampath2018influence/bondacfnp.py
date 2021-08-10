#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


###Purpose: This script calculates the bond autocorrelation function as a function of distance from a fixed nanoparticle###
###From Jon's PPPMD script###
###Dump file was first sorted using sort.py####
### Syntax: python bondacfNP.py ###
### Author: Janani Sampath ###
### Date: March 2017 ###


import sys
import numpy as np
import math
npoly = 800
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

fname = 'equilnvt.lammpstrj'
f = open(fname, 'r')
file = 'bacf_sticker.csv'

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
else:
	print >> sys.stderr, "ERROR: x coordinate not found in lammps trajectory"
if "ix" in line:
	ix_index = line.index("ix") - 2
	iy_index = line.index("iy") - 2
	iz_index = line.index("iz") - 2
else:
	print >> sys.stderr, "ERROR: x image flag not found in lammps trajectory"

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
print "reading input file..."
while frame < num_frames:
	#print frame
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

print "reading input file done, now calculating bonds and assigning bonds to shells around NP"

frames = 810
#frames = len(r)
  
if len(mol2ids[0]) == 0:
	del mol2ids[0]
    
## Radius of NP = 5; each shell is 2sigma wide
shell_min = [5,8,11,14,17]

n = np.zeros(len(shell_min))
navg = 1
## Array inititalization

bond = np.zeros([len(shell_min),frames, npoly, len(max(mol2ids,key=len)), 3], np.float) ## the arrays in mol2ids could be of different lengths as the ionomers are partially neutralized and random
bond_vector = np.zeros([len(shell_min),frames, npoly, len(max(mol2ids,key=len)), 3], np.float)
bacf_autocorr = np.zeros([len(shell_min),navg,frames], np.float) ##
bacf_autocorr_final = np.zeros([len(shell_min),frames], np.float)
 
## Fix each bond to a shell at timestep 0, this is its permanent position
 
for i in shell_min:
    for t in range(0,1):
        for mol in range(npoly):
            for molid in range(len(mol2ids[mol])-1):
                if ((id2type[mol2ids[mol][molid]] == 3) or (id2type[mol2ids[mol][molid]] == 6)):
                    continue
                if ((id2type[mol2ids[mol][molid]] == 1) or (id2type[mol2ids[mol][molid]] == 5)):                   
                    id1 = mol2ids[mol][molid] ## type 1 or 5
                    id2 = mol2ids[mol][molid+1] ## type 3 or 6
                    if (molid <= (len(mol2ids[mol]))-3): ## this is if type one is the last bead on the chain 
                        id3 = mol2ids[mol][molid+2] ## type 2
                    else:
                        id3 = 0
                    ## first bond type1 - type3
                    mid_point = (r[t][id1] + r[t][id2])/2
                    magnitude = math.sqrt(sum(j**2 for j in mid_point))
                    if (i < magnitude < i+2):
                        shell = shell_min.index(i)
                        n[shell]=n[shell]+1        
                        bond[shell][t][mol][molid+1] = r[t][id1] - r[t][id2]
                        bond_vector[shell][t][mol][molid+1] = bond[shell][t][mol][molid+1]/np.sqrt(np.dot(bond[shell][t][mol][molid+1],bond[shell][t][mol][molid+1])) ## unit vector
                    ## second bond type1 - type 2
                    if id3 != 0:
                        mid_point = (r[t][id1] + r[t][id3])/2
                        magnitude = math.sqrt(sum(j**2 for j in mid_point))
                        if (i < magnitude < i+2):
                            shell = shell_min.index(i)
                            n[shell]=n[shell]+1        
                            bond[shell][t][mol][molid] = r[t][id1] - r[t][id3]
                            bond_vector[shell][t][mol][molid] = bond[shell][t][mol][molid]/np.sqrt(np.dot(bond[shell][t][mol][molid],bond[shell][t][mol][molid])) ## unit vector
                elif (id2type[mol2ids[mol][molid]] == 2):
                    id1 = mol2ids[mol][molid]
                    id2 = mol2ids[mol][molid+1]
                    mid_point = (r[t][id1] + r[t][id2])/2
                    magnitude = math.sqrt(sum(j**2 for j in mid_point))
                    if (i < magnitude < i+2):
                        shell = shell_min.index(i)
                        n[shell]=n[shell]+1        
                        bond[shell][t][mol][molid] = r[t][id1] - r[t][id2]
                        bond_vector[shell][t][mol][molid] = bond[shell][t][mol][molid]/np.sqrt(np.dot(bond[shell][t][mol][molid],bond[shell][t][mol][molid])) ##unit vector
                     
## Now calculate bonds at different times, asigned to the shell they were in at time 0
                    
for i in range(len(shell_min)):                    
    for t in range(1,frames):
        for mol in range(npoly):
            for molid in range(len(mol2ids[mol])-1):
                if (id2type[mol2ids[mol][molid]] == 3):
                    continue
                if ((id2type[mol2ids[mol][molid]] == 1) or (id2type[mol2ids[mol][molid]] == 5)):                   
                    id1 = mol2ids[mol][molid]## type 1 or 5
                    id2 = mol2ids[mol][molid+1]## type 3 or 6
                    if (molid <= (len(mol2ids[mol]))-3):## this is if type one or five is the last bead on the chain
                        id3 = mol2ids[mol][molid+2]## type 2
                    else:
                        id3 = 0
                    bond[i][t][mol][molid+1] = r[t][id1] - r[t][id2]
                    bond_vector[i][t][mol][molid+1] = bond[i][t][mol][molid+1]/np.sqrt(np.dot(bond[i][t][mol][molid+1],bond[i][t][mol][molid+1])) ## unit vector
                    if id3 != 0:
                        bond[i][t][mol][molid] = r[t][id1] - r[t][id3]
                        bond_vector[i][t][mol][molid] = bond[i][t][mol][molid]/np.sqrt(np.dot(bond[i][t][mol][molid],bond[i][t][mol][molid])) ## unit vector
                elif (id2type[mol2ids[mol][molid]] == 2):
                    id1 = mol2ids[mol][molid]
                    id2 = mol2ids[mol][molid+1]
                    bond[i][t][mol][molid] = r[t][id1] - r[t][id2]
                    bond_vector[i][t][mol][molid] = bond[i][t][mol][molid]/np.sqrt(np.dot(bond[i][t][mol][molid],bond[i][t][mol][molid])) ## unit vector
                   
## Compute the bond autocorrelation function


print "averaging for different starting times..."

for dist in range(len(shell_min)):                
    for t in range(0,navg):
        t1 = 21*t
        for t2 in range(t1,(frames/navg)+t1):
            dt = t2 - t1
            for mol in range(npoly):
                for molid in range(len(mol2ids[mol])):
                    bacf_autocorr[dist][t][dt] += np.dot(bond_vector[dist][t1][mol][molid], bond_vector[dist][t2][mol][molid])/n[dist]



for dist in range(len(shell_min)):
    for dt in range(frames/navg):
        for avg in range(0,navg):
            bacf_autocorr_final[dist][dt] += bacf_autocorr[dist][avg][dt]/navg


## Write output to file
        
OUT = open(file, 'w')
OUT.write("shell,timesteps,BACF_type1\n")
for j in range(len(shell_min)):
    for t in range(frames/navg):
        OUT.write("%7i, %7i, %7f\n" % (j+1, (timestep[t]), (bacf_autocorr_final[j][t]))) ###
OUT.close()






























            
















        

            


    

             
            
