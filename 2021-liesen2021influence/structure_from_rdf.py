#!/usr/bin/env python3

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

# Script:  structure_from_rdf.py
# Usage: python3 structure_from_rdf.py --wdir "$working_dir" --fin "$input_fname" --rho "$density" --fout  "$desired_out_fname" 
#			          [dir w/ rdf data]  [g(r) data fname] [avg # density] [name for s(k) data]
# Author:  Nicholas Liesen
# Modified: Fall 2020


import sys
import os
import argparse
import numpy as np
#from scipy import integrate


def nanoVolCor(g_unc):  # Only relevant for nanocomposite

    phi = 0.008149895
    g_cor = g_unc*(1-phi)
    return(g_cor)


def get_spacing(x):
    """ This function quickly checks if spacing of data in a 1D numpy array (x) is uniform """

    # Find spacing
    spacings = x[1:]-x[0:-1]
    dx = spacings[0]

    # Check that spacing is uniform
    spacing_is_uniform = np.all(np.isclose(spacings, dx))
    if spacing_is_uniform:
        print("Spacing is {}".format(dx))
    else:
        sys.exit('Uniformly spaced data is required')
    
    return dx


def midpoint_int(y, x):
    """This function will integrate the given function using the midpoint method, after checking
    that spacing of data is uniform."""

    # Get spacing
    dx = get_spacing(x)

    # Find integral limits
    up_lim = x[-1] + dx/2.0
    low_lim = x[0] - dx/2.0

    # Check that lower limit is zero
    llim_is_zero = np.isclose(low_lim, 0.0)
    if llim_is_zero:
        print("Integrating from {} to {}".format(low_lim, up_lim))
    else:
        sys.exit('Lower limit is not zero in midpoint integration - will not continue with integral')

    # Perform summation
    midpoint_sum = np.sum(y*dx)
    return midpoint_sum


def compute_integral(h_array, r_array, q, up_bound):
    """This function takes in a 1D array containing h(r)=g(r)-1 and a second 1D array containing the
    corresponding r-values. These are used to calculate the integrand for s(q). This integrand is returned."""

    r_square = np.square(r_array)
    qr = q*r_array
    sin_qr = np.sin(qr)
    window_denom = np.pi*r_array/up_bound  # Numerator and denominator of Lorch window function -- correct low Q
    window_num = np.sin(window_denom)

    integrand=np.zeros(len(r_array))

    for ii in range(0, len(r_array)):
        #integrand[ii] = 4*np.pi*r_square[ii]*(sin_qr[ii]/qr[ii])*h_array[ii]
        window = window_num[ii]/window_denom[ii]
        integrand[ii] = 4*np.pi*r_array[ii]*(sin_qr[ii]/q)*h_array[ii]*window  # Corrected using "Lorch" window fcn

    integral = midpoint_int(integrand, r_array)

    return(integral)


# Add command-line arguments
parser = argparse.ArgumentParser(description="Generate S(k) via a 3D fourier transform of g(r)-1 \n\
in spherical coordinates, given a whitespace separated g(r) data file \n\ and an average number density.")
parser.add_argument('-w', '--wdir', type=str, required=True, help='Location of g(r) data file')
parser.add_argument('-i', '--fin', type=str, required=True, help='Name of g(r) data file')
parser.add_argument('-r', '--rho', type=float, required=True, help='Average number density of system')
parser.add_argument('-o', '--fout', type=str, required=True, help='Desired output filename for S(k) data file')
args = parser.parse_args()

# Process passed args
wdir = args.wdir
in_fname = args.fin
rho = args.rho
out_fname = args.fout

g_r = np.loadtxt(wdir+"/"+in_fname)
g_r[:,1] = g_r[:,1]-1.0  # Broadcast subtraction by 1 to create h(r)
h_r = g_r

"""
Minimum points needed to not "miss" sine wave is 3
<--------|-------->
-dt      0       dt
wavelength = 2dt, wavenumber = 2pi/(2dt)
resolution is based on this wavenumber (q)
so we can choose q-values in increments of 2pi/(2*Rmax)
i.e. pi/Rmax, 2pi/Rmax, 3pi/Rmax,...
selecting only values below our frequency limit of pi/dt
A common choice is Rmax = L/2, leading to the k-values
of (2pi/L)*n where n is the set of integers (1,2,3,...)
where k is below the frequency limit.
"""

r_dist = h_r[:,0]
h_vals = h_r[:,1]

# Check that spacing is uniform and return spacing
delta_r = get_spacing(r_dist)
npts = np.shape(r_dist)[0]

# Choice of wavevectors consistent with Amelie Frischknecht's box notes and
# Jeff's fftgofr.m script.
#sample_freq = 1.0/delta_r
#fmax = sample_freq/2.0
#qmax = 2.0*np.pi*fmax

last_r = r_dist[-1]

fmin = 1.0/(2.0*npts*delta_r)
qmin = 2.0*np.pi*fmin
q_vals = np.arange(1.0, npts+1.0)
q_vals = q_vals*qmin
#q_vals = np.arange(1.0*qmin, qmax, qmin)

qmax = q_vals[-1]
Rmax = last_r + delta_r/2.0  # Rmax ~ L/2
integral = np.zeros(len(q_vals))

jj=0
for q_mag in q_vals:  # Integrate for each choice of q-vector magnitude
    integral[jj] = compute_integral(h_vals, r_dist, q_mag, Rmax)
    jj = jj + 1

S_q = 1+rho*integral

Sq_vs_q = np.transpose(np.vstack((q_vals, S_q)))
np.savetxt(wdir+"/"+out_fname, Sq_vs_q, delimiter=' ')
