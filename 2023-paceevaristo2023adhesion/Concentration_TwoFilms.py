# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

### Calculate concentration profile similar to Figure 6 in https://dx.doi.org/10.1021/acs.macromol.0c01508. Profiles are calculated along the z direction with z positions ranging from zlo < 0 to zhi > 0 and |zlo| = |zhi|. Healing interface is shifted to z = 0. ###
### From Jon's PPPMD script ###
### Syntax: python2 Concentration_TwoFilms.py ###
### Author: Janani Sampath ###
### Date: Oct 2020 ###
### Modifier: Felipe F. Pacci Evaristo ###
### Date: Nov 2022 ###

import sys
import numpy as np
import math as m

np.set_printoptions(threshold = np.inf)

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

def wrap(frame, r, box_bounds):
    t = frame
    print 'Wrapping frame %d...' % (t)

    bound_lo = np.array([box_bounds[t][0][0], box_bounds[t][1][0], box_bounds[t][2][0]])
    bound_hi = np.array([box_bounds[t][0][1], box_bounds[t][1][1], box_bounds[t][2][1]])

    boxsize = bound_hi - bound_lo

    for atom in range(1, len(r[t])):
        shift = np.zeros(3, np.float)
        for axis, coord in enumerate(r[t][atom]):
            if (coord > bound_hi[axis]):
                shift[axis] = - boxsize[axis]
            elif (coord < bound_lo[axis]):
                shift[axis] = boxsize[axis]
        if (np.sum(np.abs(shift)) != 0):
            r[t][atom] += shift

    return r

print "Reading input file..."

fname = '3_50.lammpstrj'
f = open(fname, 'r')
file = '3_50_interface_total.csv'

# Set the width of the bins to the diameter of the monomers in the system, in LJ units.
bin_width = 1.

# Read in the initial header.
frame = 0
init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = read_header(f)

# It is not possible to preallocate arrays (or to know the number of frames in advance).
num_frames = float('inf')
alloc = 1
inf_frames = True

timestep = np.zeros(alloc, np.int) # 1D array of timesteps.
box_bounds = np.zeros([alloc, 3, 2], np.float) # 3D array to store boundaries of the box, indexed by frame, x/y/z, then lower/upper.
box_dimensions = np.zeros([alloc, 3], np.float) # 2D array to store dimensions of the box, indexed by frame, then x/y/z.
bin_volume = np.zeros(alloc, np.float) # 1D array to store volumes of the bins, indexed by frame.

timestep[frame] = init_timestep
box_bounds[frame][0][0] = xlo
box_bounds[frame][0][1] = xhi
box_bounds[frame][1][0] = ylo
box_bounds[frame][1][1] = yhi
box_bounds[frame][2][0] = zlo
box_bounds[frame][2][1] = zhi
box_dimensions[frame][0] = box_bounds[frame][0][1] - box_bounds[frame][0][0]
box_dimensions[frame][1] = box_bounds[frame][1][1] - box_bounds[frame][1][0]
box_dimensions[frame][2] = box_bounds[frame][2][1] - box_bounds[frame][2][0]
bin_volume[frame] = box_dimensions[frame][0] * box_dimensions[frame][1] * bin_width

# Note: num_atoms + 1 is used here so that the arrays are indexed by their LAMMPS atom ID.
r = np.zeros([alloc, num_atoms + 1, 3], np.float) # 3D array of x, y, z coordinates; r[frame][id][coordinate].
r_0 = np.zeros([num_atoms + 1, 3], np.float) # 2D array to store initial unscaled unwrapped x, y, z coordinates; r_0[id][coordinate].
ir = np.zeros([alloc, num_atoms + 1, 3], np.int) # 3D array of x, y, z image flags; r[frame][id][coordinate].

id2mol = np.zeros(num_atoms + 1, np.int) # Array to map from atom ID to molecule ID. Build this from the first frame, if available.
id2type = np.zeros(num_atoms + 1, np.int) # Array to map from atom ID to type. Build this from the first frame, if available.

# Separately do the first ATOMS section (so that variables can be initialized) and build the id2mol and id2type arrays (so that the main loop starts with reading in the header).
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

# Loop over the ATOMS lines.
for atom in range(num_atoms):
	line = f.readline()
	line = line.split()
	# Get the atom ID.
	my_id = int(line[id_index])
	# Get the x, y, z coordinates.
	r[frame][my_id][0] = float(line[x_index])
	r[frame][my_id][1] = float(line[y_index])
	r[frame][my_id][2] = float(line[z_index])

	# Get the x, y, z image flags.
	ir[frame][my_id][0] = int(line[ix_index])
	ir[frame][my_id][1] = int(line[iy_index])
	ir[frame][my_id][2] = int(line[iz_index])

	# Unscale, if necessary.
	if scaled:
		r[frame][my_id][0] = r[frame][my_id][0] * box_dimensions[frame][0] - box_dimensions[frame][0] / 2
		r[frame][my_id][1] = r[frame][my_id][1] * box_dimensions[frame][1] - box_dimensions[frame][1] / 2
		r[frame][my_id][2] = r[frame][my_id][2] * box_dimensions[frame][2] - box_dimensions[frame][2] / 2
		# Calculate and store initial unscaled unwrapped x, y, z coordinates.
		r_0[my_id][0] = r[frame][my_id][0] + ir[frame][my_id][0] * box_dimensions[frame][0]
		r_0[my_id][1] = r[frame][my_id][1] + ir[frame][my_id][1] * box_dimensions[frame][1]
		r_0[my_id][2] = r[frame][my_id][2] + ir[frame][my_id][2] * box_dimensions[frame][2]

	# If available, build the i2mol and id2type arrays.
	if mol_index is not None:
		id2mol[my_id] = int(line[mol_index])
	if type_index is not None:
		id2type[my_id] = int(line[type_index])

# Build the reverse of the id2mol array.
# This is a 2D array with rows of (potentially) varying length, so nest a NumPy array into a Python list.
if mol_index is not None:
	num_mols = id2mol.max()
	mol2ids = [[]]
	for molid in range(1, num_mols + 1):
		mol2ids.append(np.where(id2mol == molid)[0])
else:
	num_mols = None
	mol2ids = None

# Enforce that absolutely all atoms are wrapped within the box at the given frame.
r = wrap(frame, r, box_bounds)

# Loop over the number of frames num_frames. If num_frames is infinite, loop over all frames in the file.
frame = 1
while (frame < num_frames):
	# Print frame.
	# Try to read in a new header.
	try:
		my_timestep, my_num_atoms, my_xlo, my_xhi, my_ylo, my_yhi, my_zlo, my_zhi = read_header(f)
	except:
		# print >> sys.stderr, "Warning: end of file reached when reading frame ", frame, " in LAMMPS trajectory file ", fname, "."
		break

	# If the number of frames to read in are not known, more memory needs to be allocated for the arrays.
	if inf_frames:
		timestep = np.append(timestep, 0)
		box_bounds = np.concatenate((box_bounds, np.zeros([1, 3, 2], np.float)))
		box_dimensions = np.concatenate((box_dimensions, np.zeros([1, 3], np.float)))
		bin_volume = np.append(bin_volume, 0.)

		r = np.concatenate((r, np.zeros([1, num_atoms + 1, 3], np.float)))
		ir = np.concatenate((ir, np.zeros([1, num_atoms + 1, 3], np.float)))

	# Update the timestep and box size arrays.
	timestep[frame] = my_timestep
	box_bounds[frame][0][0] = my_xlo
	box_bounds[frame][0][1] = my_xhi
	box_bounds[frame][1][0] = my_ylo
	box_bounds[frame][1][1] = my_yhi
	box_bounds[frame][2][0] = my_zlo
	box_bounds[frame][2][1] = my_zhi
	box_dimensions[frame][0] = box_bounds[frame][0][1] - box_bounds[frame][0][0]
	box_dimensions[frame][1] = box_bounds[frame][1][1] - box_bounds[frame][1][0]
	box_dimensions[frame][2] = box_bounds[frame][2][1] - box_bounds[frame][2][0]
	bin_volume[frame] = box_dimensions[frame][0] * box_dimensions[frame][1] * bin_width

	f.readline() # ITEM: ATOMS
	# Loop over the ATOMS lines.
	for atom in range(num_atoms):
		line = f.readline()
		line = line.split()

		# Get the atom ID.
		my_id = int(line[id_index])

		# Get the x, y, z coordinates.
		r[frame][my_id][0] = float(line[x_index])
		r[frame][my_id][1] = float(line[y_index])
		r[frame][my_id][2] = float(line[z_index])

		# Get the x, y, z image flags.
		ir[frame][my_id][0] = int(line[ix_index])
		ir[frame][my_id][1] = int(line[iy_index])
		ir[frame][my_id][2] = int(line[iz_index])

		# Unscale, if necessary.
		if scaled:
			r[frame][my_id][0] = r[frame][my_id][0] * box_dimensions[frame][0] - box_dimensions[frame][0] / 2
			r[frame][my_id][1] = r[frame][my_id][1] * box_dimensions[frame][1] - box_dimensions[frame][1] / 2
			r[frame][my_id][2] = r[frame][my_id][2] * box_dimensions[frame][2] - box_dimensions[frame][2] / 2

	# Enforce that absolutely all atoms are wrapped within the box at the given frame.
	r = wrap(frame, r, box_bounds)

	frame += 1

num_frames = frame

# Enforce that the center of mass of the film and the center of the box coincide.
shift_cm = np.zeros([num_frames, 3], np.float)
shift_cm[0] = - np.sum(r_0, axis = 0) / num_atoms

if (((shift_cm[0][2] < 0) & (abs(shift_cm[0][2]) <= abs(- box_dimensions[0][2] / 2 - np.min(r[0, :, 2])))) | ((shift_cm[0][2] >= 0) & (abs(shift_cm[0][2]) <= abs(box_dimensions[0][2] / 2 - np.max(r[0, :, 2]))))):
	for atom in range(1, num_atoms + 1):
		r[0][atom] += shift_cm[0]
else:
	print "Warning: Shifting moves atoms outside of the box at frame 0. Wrapping of atoms that are moved outside of the box must be performed. Current version of the code does not do that."

for frame in range(1, num_frames):
	shift_cm[frame] = (shift_cm[0] / box_dimensions[0]) * box_dimensions[frame]
	if (((shift_cm[frame][2] < 0) & (abs(shift_cm[frame][2]) <= abs(- box_dimensions[frame][2] / 2 - np.min(r[frame, :, 2])))) | ((shift_cm[frame][2] >= 0) & (abs(shift_cm[frame][2]) <= abs(box_dimensions[frame][2] / 2 - np.max(r[frame, :, 2]))))):
		for atom in range(1, num_atoms + 1):
			r[frame][atom] += shift_cm[frame]
	else:
		print "Warning: Shifting moves atoms outside of the box at frame %d. Wrapping of atoms that are moved outside of the box must be performed. Current version of the code does not do that." % (frame)

# Shift z positions to shift interface to z = 0.
for frame in range(0, num_frames):
	for atom in range(1, num_atoms + 1):
		if (r[frame][atom][2] < 0):
			r[frame][atom][2] = r[frame][atom][2] + box_dimensions[frame][2] / 2
		elif (r[frame][atom][2] > 0):
			r[frame][atom][2] = r[frame][atom][2] - box_dimensions[frame][2] / 2

print "Calculating concentration profiles..."

# Initialize frames.
frames = len(r)
if (len(mol2ids[0]) == 0):
	del mol2ids[0]

# Set bins for each frame.
bins = [[] for _ in range(frames)]
for t in range(num_frames):
    count = m.ceil(- box_dimensions[t][2] / 2)
    while (count <= box_bounds[t][2][1] + bin_width):
        count_round = round(count, 2)
        bins[t].append(count)
        count = count + bin_width

# Initialize left and right.
left = np.zeros([frames, (len(max(bins, key = len)))])
right = np.zeros([frames, (len(max(bins, key = len)))])

# Count beads in each bin.
for time in range(frames):
    for i in bins[time]:
        Bin = bins[time].index(i)
        for atom in range(1, num_atoms + 1):
            # Beads are divided into left and right according to their position at t = 0.
            if (i < r[time][atom][2] < i + bin_width):
                if (ir[0][atom][2] == 0):
                    if (r[0][atom][2] < 0):
                        left[time][Bin] += 1
                    if (r[0][atom][2] > 0):
                        right[time][Bin] += 1
                else:
                    if (ir[0][atom][2] < 0):
                        left[time][Bin] += 1
                    if (ir[0][atom][2] > 0):
                        right[time][Bin] += 1

        # Calculate concentration profile.
        left[time][Bin] = left[time][Bin] / bin_volume[time]
        right[time][Bin] = right[time][Bin] / bin_volume[time]

# Write output to file.
OUT = open(file, 'w')
for time in range(0, frames):
    OUT.write("Frame Index\n")
    OUT.write("%7i\n" % (time))
    OUT.write("z Position,Left Concentration,Right Concentration\n")
    for b in range(0, len(bins[time])):
        # In the output file, concentrations are associated with the midpoints of their respective bins.
        OUT.write("%.3f,%.9f,%.9f\n" % (m.ceil(- box_dimensions[time][2] / 2) + (1. / 2. + b) * bin_width, left[time][b], right[time][b]))
