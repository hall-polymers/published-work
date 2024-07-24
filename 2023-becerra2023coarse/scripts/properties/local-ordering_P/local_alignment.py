import numpy as np
import math
from matplotlib import pyplot as plt
import time
from numba import jit
from pyquaternion import Quaternion
from scipy.spatial.transform import Rotation as Ro

# Adjustable Parameters
#-------------------------------------------------------------------------------------------------------------------
# Global formatting parameters:
pos_dimension = 3           # dimension of position vector
quat_dimension = 4          # dimension of the quaternion vector
skip_lines = 9              # number of lines before first line in input with particle properties
lines_before_dims = 5       # number of lines before box dimensions in input file
type_position = 1           # number of spaces in particle properties line before atom type
xcoord_position = 2         # number of spaces in particle properties line before first position coordinate
first_quat_position = 11    # number of spaces in particle properties line before first quaternion coordinate

# User-defined parameters
bin_size = 3.4                              # bin size
dot_prod_min = 0.8                          # minimum dot product between two vectors
atom_name = '3'                             # name of particle to be read
axis = [0, 0, 1]                            # alignment axis

base_dir = './'    # base directory
temperatures = ['1.5', '2.5']                    # all temperatures
case = 'sosc_100_6_1.5'                                                                # simulation case
outfile = base_dir + 'results/' + case + '.dat'                                  # output file name


# Average interval
use_all_frames = False  # set to be true if one desires to use every frame in the input
frames = [1, 5, 9]



#Functions
#-------------------------------------------------------------------------------------------------------------------
# reads positions of particular particles into 3D list
# inputs:   str file_name for address of trajectory file
#           str particle_type for ID of particular particle type
# outputs:  list of float arrays with entries in format [[x_1, y_1, z_1], ..., [x_n, y_n, z_n]]
#           with each entry corresponding to a frame and each list in each entry being coordinates of a particle
def read_particle_positions(file_name: str, particle_type: str) -> list:

    line_number = 0
    output = []
    frame = []

    #read data from file into list
    f = open(file_name, 'r')
    data = f.readlines()
    f.close()

    while(line_number < len(data)):
        if data[line_number].startswith('ITEM: TIMESTEP'):
            if (line_number != 0):
                output.append(frame)
                frame = []
            line_number += skip_lines
        else:
            split_line = data[line_number].split()
            if (split_line[type_position] == particle_type):
                temp = []
                for i in range(0, pos_dimension):
                    temp.append(float(split_line[xcoord_position + i]))
                frame.append(temp)
            line_number += 1
    output.append(frame)
    return output



# reads positions of particular particles into 3D list
# inputs:   str file_name for address of trajectory file
#           str particle_type for ID of particular particle type
# outputs:  list of float arrays with entries in format [[x_1, y_1, z_1], ..., [x_n, y_n, z_n]]
#           with each entry corresponding to a frame and each list in each entry being coordinates of a particle
def read_particle_quaternion(file_name: str, particle_type: str) -> list:

    line_number = 0
    output = []
    frame = []

    #read data from file into list
    f = open(file_name, 'r')
    data = f.readlines()
    f.close()

    while(line_number < len(data)):
        if data[line_number].startswith('ITEM: TIMESTEP'):
            if (line_number != 0):
                output.append(frame)
                frame = []
            line_number += skip_lines
        else:
            split_line = data[line_number].split()
            if (split_line[type_position] == particle_type):
                temp = []
                for i in range(0, quat_dimension):
                    temp.append(float(split_line[first_quat_position + i]))
                frame.append(temp)
            line_number += 1
    output.append(frame)
    return output



# Returns distances of all other vectors in array from given vector
# Inputs:   int num_positions representing total number of positions in array
#           list vector representing given vector
#           list position_array representing all vectors in a given frame
@jit(nopython=True)
def position_distances(num_positions, vector, position_array):
    distances = np.zeros(num_positions)
    for j in range(0, num_positions):
        dist = np.linalg.norm(vector - position_array[j])
        distances[j] = dist
    return distances

#@jit(nopython=True)
def quat_to_vect(quaternion_array, vector_to_rotate):
    quaternion = Quaternion([quaternion_array[1], quaternion_array[2], quaternion_array[3], quaternion_array[0]])
    vector = quaternion.rotate(vector_to_rotate)
    return vector



#Main Method
#----------------------------------------------------------------------------------------------------
f = open(outfile, 'w')

for temperature in temperatures:
    in_file = base_dir + 'traj_files/' + case + '/trajT' + temperature + '.lammpstrj'
    print(in_file)

    positions = np.array(read_particle_positions(in_file, atom_name))
    quaternions = np.array(read_particle_quaternion(in_file, atom_name))

    if use_all_frames:
        start_frame = 0
        end_frame = len(positions)

    total_particles = 0
    aligned_neighbors = 0

    for frame_number in frames:
        for i in range(0, len(positions[frame_number])):
            neighbors = []
            distances = position_distances(len(positions[frame_number]), positions[frame_number][i], positions[frame_number])
            for j in range(0, len(distances)):
                if (distances[j] < bin_size):
                    neighbors.append(j)

            base_vect = quat_to_vect(quaternions[frame_number][i], axis)
            for entry in neighbors:
                if(abs(np.dot(base_vect, quat_to_vect(quaternions[frame_number][entry], axis))) > dot_prod_min):
                    aligned_neighbors += 1
            aligned_neighbors -= 1
            total_particles += (len(neighbors) - 1)

        print('Frame number: ' + str(frame_number) + '\n\n')

    output = aligned_neighbors / total_particles
    f.write(temperature + ' ' + str(output) + '\n')

f.close()
