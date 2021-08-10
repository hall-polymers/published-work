# !/usr/bin/env pypy

#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

'''
Script to generate unconfined or confined random walks in the presence of walls or nanoparticles
the fundamental assumption is that when the RW would pass through a wall (or NP surface), it instead reflects off

For now, it is also assumed that the box is a unit cube with appropriate periodic boundary conditions
(or equivalently, we only start in the unit cube, but the system expands into all space with the appropriate symmetry)

input paramters are 
 M - number of chains 
 N - number of steps in the chain
 L - boxsize (set realative to the unit step size)
 constraint - 
   none: no constraint 
   lamellae: assumes that the z=0 and z=L planes are walls 
   cylinder: assumes that the wall is a cylinder surface along z-dir, and its center is (L/2, L/2) and radius is 2/L
   gyroid: assumes that the wall is the surface of a double gyroid unit cell

optional parameters
  x0, y0, z0 - the inital location of the RW, a negative number will choose random location between 0 and L FIXME: does not check that these choices are reasonable, so an infinite loop can happen if you force the system to start inside a particle
  outputchains - boolean to determine if printing to lammpstrj format

edited J.Brown 20140404 modified Kevin Shen 20170711
'''

from sys import stderr
from math import pi, sin, cos, sqrt, floor
from random import random, randint, seed
from scipy import optimize
import numpy as np
#seed(2)
#from random import SystemRandom
#sysrand = SystemRandom()
#random = sysrand.random

def crw(M, N, L, constraint, x0=-1, y0=-1, z0=-1):
	# set the starting locations for the RWs
	if x0 < 0:
		xstart = lambda: L*random()
	else:
		xstart = lambda: x0
	if y0 < 0:
		ystart = lambda: L*random()
	else:
		ystart = lambda: y0
	if z0 < 0:
		zstart = lambda: L*random()
	else:
		zstart = lambda: z0

	dumpSteps = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576]
	dumpArr = [0]* len(dumpSteps)

	# for gyroid surface function
	rescale = 2*pi/L
	fA = 0.35
	gy_crit = (1-fA)/0.067
	gy1 = lambda x,y,z : 10.0*(cos(x*rescale) * sin(y*rescale) + cos(y*rescale) * sin(z*rescale) + cos(z*rescale) * sin(x*rescale)) - 0.5*(cos(2*x*rescale)*cos(2*y*rescale) + cos(2*y*rescale)*cos(2*z*rescale) + cos(2*z*rescale)*cos(2*x*rescale)) - gy_crit
	gy2 = lambda x,y,z : -10.0*(cos(x*rescale) * sin(y*rescale) + cos(y*rescale) * sin(z*rescale) + cos(z*rescale) * sin(x*rescale)) - 0.5*(cos(2*x*rescale)*cos(2*y*rescale) + cos(2*y*rescale)*cos(2*z*rescale) + cos(2*z*rescale)*cos(2*x*rescale)) - gy_crit

	def gyroid(d, in_gy1, r0, r1, r2, dx, dy, dz):
		x = r0+dx*d
		y = r1+dy*d
		z = r2+dz*d
		if in_gy1:
			return round(gy1(x, y ,z), 10)
		return round(gy2(x, y ,z), 10)

	r = [0.0, 0.0, 0.0]
	D = 0
	
	if constraint in ["none", "lamellae"]: # bulk (no constraint) or reflective walls constraint
		# loop over all the chains
		for m in xrange(M):
			# place the starting location 
			ir = [0, 0, 0]
			r[0] = xstart()
			r[1] = ystart()
			r[2] = zstart()
		
			r0 = r[:]

			# initialize the previous step (this makes the first step free to go in any direction, so long as the required angle is < 90)
			dx_prev = 0
			dy_prev = 0
			dz_prev = 0

			# loop over all the links in the chain
			for n in xrange(1, N+1):
				# choose a random point on a unit sphere
				dz = 2*random()-1.0
				theta = 2*pi*random()
				ds = sqrt(1-dz*dz)
				dx = ds*cos(theta)
				dy = ds*sin(theta)

				#while the angle/dot product is too low, retry
				#while (dx_prev or dy_prev or dz_prev) and dx*dx_prev + dy*dy_prev + dz*dz_prev < 0.9:
				#while dx*dx_prev + dy*dy_prev + dz*dz_prev < -0.5:
				#	dz = 2*random()-1.0
				#	theta = 2*pi*random()
				#	ds = sqrt(1-dz*dz)
				#	dx = ds*cos(theta)
				#	dy = ds*sin(theta)
		
				r[0] += dx
				r[1] += dy
				r[2] += dz
				
				dx_prev = dx
				dy_prev = dy
				dz_prev = dz

				if constraint == "lamellae": # walls at z=0 and z=L
					#if we ended up past a wall, enforce the reflective BC
					if r[2] < 0 or r[2] > L:
						done = False
						while not done:
							if r[2] < 0:
								r[2] = -r[2]
							elif r[2] > L:
								r[2] = 2*L - r[2] 
							else:
								done = True

				# enforce periodic BCs
				ir[0] += floor(r[0]/L)
				r[0] = L*((r[0]/L)%1)
				ir[1] += floor(r[1]/L)
				r[1] = L*((r[1]/L)%1)
				ir[2] += floor(r[2]/L)
				r[2] = L*((r[2]/L)%1)
				#print ir, r

				if n in dumpSteps:
					dumpArr[dumpSteps.index(n)] += ((L*ir[0]+r[0]-r0[0])**2 + (L*ir[1]+r[1]-r0[1])**2 + (L*ir[2]+r[2]-r0[2])**2)/M
	
	elif constraint in ["cylinder"]: # cylindrical reflective wall
		# loop over all the chains
		for m in xrange(M):
			# place the starting location, as long as it's not outside the cylinder
			ir = [0, 0, 0]
			r[0] = xstart()
			r[1] = ystart()
			r[2] = zstart()
			while (r[0] - L/2)**2 + (r[1] - L/2)**2 > (L/2)**2:
				r[0] = xstart()
				r[1] = ystart()
				r[2] = zstart()
		
			r0 = r[:]
			
			# initialize the previous step (this makes the first step free to go in any direction, so long as the required angle is < 90)
			dx_prev = 0
			dy_prev = 0
			dz_prev = 0
			
			# loop over all the links in the chain
			for n in xrange(1, N+1):
				# choose a random point on a unit sphere
				dz = 2*random()-1.0
				theta = 2*pi*random()
				ds = sqrt(1-dz*dz)
				dx = ds*cos(theta)
				dy = ds*sin(theta)
				dx_0 = ds*cos(theta)
				dy_0 = ds*sin(theta)

				#while the angle/dot product is too low, retry
				#while (dx_prev or dy_prev or dz_prev) and dx*dx_prev + dy*dy_prev + dz*dz_prev < 0.9:
				#while dx*dx_prev + dy*dy_prev + dz*dz_prev < -0.5:
				#	dz = 2*random()-1.0
				#	theta = 2*pi*random()
				#	ds = sqrt(1-dz*dz)
				#	dx = ds*cos(theta)
				#	dy = ds*sin(theta)
	
				# # pull up the list of centers of spheres in the current 
				# # note that this only needs to be done once even if we have multiple reflections etc
				# octant = int(r[0]/(L/2)) + 2*int(r[1]/(L/2)) + 4*int(r[2]/(L/2))
				# spheres = NP_octant[octant]
		
				# check for reflections as long as there is distance remaining to be traveled
				dist_remaining = ds
				reflection = True
				# on_surface = False
				# loop over reflections
				# i = 0
				while reflection:
					reflection = False

					xc = L/2
					yc = L/2
					# compute the distance to the cylinder center, allowing quick disqualification
					centerdistsq = (r[0]-xc)**2 + (r[1]-yc)**2

					# make xy vector to be unit vector for calculation
					rsqrt_xy = sqrt(dx*dx+dy*dy)
					dx = dx/rsqrt_xy
					dy = dy/rsqrt_xy
					
					# if close enough from the center initially, we can ignore it, otherwise have to check reflections
					if centerdistsq < ((L/2)-1)**2:
						# default to infinite distance
						d = float("inf")
						continue

					# we're close enough to the cylindrical surface
					
					# find the intersection with the sphere, if there is one, using the quadratic formula
					# see: http://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
					a = dx**2 + dy**2
					b = 2*(dx*(r[0]-xc) + dy*(r[1]-yc))
					c = (r[0]-xc)**2 + (r[1]-yc)**2 - (L/2)**2
					radicand = b**2 - 4*a*c
	
					# if the term under the radical is negative, no intersection, otherwise test further for a reflection
					if radicand < 0:
						# default to infinite distance
						print >> stderr, "Error: A bead has been out of the cylinder, should be impossible."
						return

					# need to find the point of reflection and move from there
					d1 = (-b - sqrt(radicand))/(2*a)
					d2 = (-b + sqrt(radicand))/(2*a)
	
					# check if the intersection is outside of our range, if it is we can just move on
					# this is for the case where the floating number makes it slightly outside of the cylinder
					if d1 > 0 and d2 > 0:
						if d1 > d2:
							d = d1
						else:
							d = d2
					# normal cases
					elif d1 > 0 and d1 < dist_remaining:
						d = d1
					elif d2 > 0 and d2 < dist_remaining:
						d = d2
					else:
						# default to infinite distance
						d = float("inf")

					if d < dist_remaining:
						# found a reflection!
						dist_remaining = dist_remaining - d
						
						r[0] += d*dx
						r[1] += d*dy

						# still have to move l-d in the reflected direction (simplified by the fact that the current r is on the sphere)
						# see http://en.wikipedia.org/wiki/Reflection_%28mathematics%29
						drdotr = dx*(r[0]-xc) + dy*(r[1]-yc)
						dx = dx - 2*drdotr*(r[0]-xc)/((L/2)**2)
						dy = dy - 2*drdotr*(r[1]-yc)/((L/2)**2)
		
						# the second movement we don't explicitly perform, instead we just change dx, dy, and dz appropriately, and restart the loop
						reflection = True

	
				if dist_remaining == ds:
					r[0] += dx_0
					r[1] += dy_0
					r[2] += dz
				else:
					# move int he dx, dy, dz direction
					r[0] += dist_remaining*dx
					r[1] += dist_remaining*dy
					r[2] += dz

				# the vector that we check against for the nonreversing restriction is the last direction we moved (i.e. the reflected direction if there was a reflection)
				dx_prev = dx
				dy_prev = dy
				dz_prev = dz
	
				# enforce periodic BCs
				ir[0] += floor(r[0]/L)
				r[0] = L*((r[0]/L)%1)
				ir[1] += floor(r[1]/L)
				r[1] = L*((r[1]/L)%1)
				ir[2] += floor(r[2]/L)
				r[2] = L*((r[2]/L)%1)
				
				if ir[0] != 0 or ir[1] != 0:
					# image flag of x and y should not be other than 0
					print >> stderr, "Image Flag Warning: A bead has been out of the cylinder, should be impossible."
					return 
	
				if n in dumpSteps:
					dumpArr[dumpSteps.index(n)] += ((L*ir[0]+r[0]-r0[0])**2 + (L*ir[1]+r[1]-r0[1])**2 + (L*ir[2]+r[2]-r0[2])**2)/M
	
	elif constraint in ["gyroid"]: # gyroid reflective wall
		# loop over all the chains
		for m in xrange(M):
			# place the starting location, as long as it's not outside the cylinder
			ir = [0, 0, 0]
			r[0] = xstart()
			r[1] = ystart()
			r[2] = zstart()
			while max(gy1(r[0], r[1], r[2]), gy2(r[0], r[1], r[2])) < 0:
				r[0] = xstart()
				r[1] = ystart()
				r[2] = zstart()

			r0 = r[:]
			
			# initialize the previous step (this makes the first step free to go in any direction, so long as the required angle is < 90)
			dx_prev = 0
			dy_prev = 0
			dz_prev = 0

			# determine which network it is in
			in_gy1 = False
			if gy1(r[0], r[1], r[2]) > 0:
				in_gy1 = True
			
			# loop over all the links in the chain
			for n in xrange(1, N+1):
				# choose a random point on a unit sphere
				dz = 2*random()-1.0
				theta = 2*pi*random()
				ds = sqrt(1-dz*dz)
				dx = ds*cos(theta)
				dy = ds*sin(theta)

				#while the angle/dot product is too low, retry
				#while (dx_prev or dy_prev or dz_prev) and dx*dx_prev + dy*dy_prev + dz*dz_prev < 0.9:
				#while dx*dx_prev + dy*dy_prev + dz*dz_prev < -0.5:
				#	dz = 2*random()-1.0
				#	theta = 2*pi*random()
				#	ds = sqrt(1-dz*dz)
				#	dx = ds*cos(theta)
				#	dy = ds*sin(theta)

				# check for reflections as long as there is distance remaining to be traveled
				dist_remaining = 1.0
				reflection = True
				while reflection:
					reflection = False

					if gyroid(0, in_gy1, r[0], r[1], r[2], dx, dy, dz) < 0:
						print >> stderr, "Error: A bead has been out of the gyroid, should be impossible."
						return 

					elif gyroid(0, in_gy1, r[0], r[1], r[2], dx, dy, dz) > 0.1 and gyroid(dist_remaining, in_gy1, r[0], r[1], r[2], dx, dy, dz) > 0.1:
						d = float("inf")
						continue

					elif min(optimize.fmin_l_bfgs_b(gyroid, dist_remaining, args=(in_gy1, r[0], r[1], r[2], dx, dy, dz), bounds=[(0,dist_remaining)], approx_grad=True)[1], optimize.fmin_l_bfgs_b(gyroid, 0, args=(in_gy1, r[0], r[1], r[2], dx, dy, dz), bounds=[(0,dist_remaining)], approx_grad=True)[1]) < 0:
						init_position = 0
						final_position = dist_remaining
						
						if gyroid(0, in_gy1, r[0], r[1], r[2], dx, dy, dz) <= 0 and gyroid(dist_remaining, in_gy1, r[0], r[1], r[2], dx, dy, dz) < 0:
							init_position = optimize.minimize_scalar(lambda d: -gyroid(d, in_gy1, r[0], r[1], r[2], dx, dy, dz), bounds=(0, dist_remaining), method='bounded').x

						elif gyroid(dist_remaining, in_gy1, r[0], r[1], r[2], dx, dy, dz) >= 0:
							final_position = optimize.minimize_scalar(gyroid, bounds=(0, dist_remaining), method='bounded', args=(in_gy1, r[0], r[1], r[2], dx, dy, dz)).x

						d = optimize.brentq(gyroid, init_position, final_position, args=(in_gy1, r[0], r[1], r[2], dx, dy, dz))

						if d < dist_remaining:
							# found a reflection!
							dist_remaining = dist_remaining - d
							r[0] += d*dx
							r[1] += d*dy
							r[2] += d*dz

							# still have to move l-d in the reflected direction 
							gradient_x = (-10*sin(r[0]*rescale)*sin(r[1]*rescale)+10*cos(r[2]*rescale)*cos(r[0]*rescale)+sin(2*r[0]*rescale)*cos(2*r[1]*rescale)+cos(2*r[2]*rescale)*sin(2*r[0]*rescale))*rescale 
							gradient_y = (-10*sin(r[1]*rescale)*sin(r[2]*rescale)+10*cos(r[0]*rescale)*cos(r[1]*rescale)+sin(2*r[1]*rescale)*cos(2*r[2]*rescale)+cos(2*r[0]*rescale)*sin(2*r[1]*rescale))*rescale 
	 						gradient_z = (-10*sin(r[2]*rescale)*sin(r[0]*rescale)+10*cos(r[1]*rescale)*cos(r[2]*rescale)+sin(2*r[2]*rescale)*cos(2*r[0]*rescale)+cos(2*r[1]*rescale)*sin(2*r[2]*rescale))*rescale 
	 						if not in_gy1:
	 							gradient_x = (10*sin(r[0]*rescale)*sin(r[1]*rescale)-10*cos(r[2]*rescale)*cos(r[0]*rescale)+sin(2*r[0]*rescale)*cos(2*r[1]*rescale)+cos(2*r[2]*rescale)*sin(2*r[0]*rescale))*rescale 
								gradient_y = (10*sin(r[1]*rescale)*sin(r[2]*rescale)-10*cos(r[0]*rescale)*cos(r[1]*rescale)+sin(2*r[1]*rescale)*cos(2*r[2]*rescale)+cos(2*r[0]*rescale)*sin(2*r[1]*rescale))*rescale 
		 						gradient_z = (10*sin(r[2]*rescale)*sin(r[0]*rescale)-10*cos(r[1]*rescale)*cos(r[2]*rescale)+sin(2*r[2]*rescale)*cos(2*r[0]*rescale)+cos(2*r[1]*rescale)*sin(2*r[2]*rescale))*rescale
							vdota_over_adota = (dx*gradient_x+dy*gradient_y+dz*gradient_z)/(gradient_x*gradient_x+gradient_y*gradient_y+gradient_z*gradient_z)
							dx = dx - 2*vdota_over_adota*gradient_x
							dy = dy - 2*vdota_over_adota*gradient_y
							dz = dz - 2*vdota_over_adota*gradient_z

							# the second movement we don't explicitly perform, instead we just change dx, dy, and dz appropriately, and restart the loop
							reflection = True

				# move int he dx, dy, dz direction
				r[0] += dist_remaining*dx
				r[1] += dist_remaining*dy
				r[2] += dist_remaining*dz

				# the vector that we check against for the nonreversing restriction is the last direction we moved (i.e. the reflected direction if there was a reflection)
				dx_prev = dx
				dy_prev = dy
				dz_prev = dz

				# enforce periodic BCs
				ir[0] += floor(r[0]/L)
				r[0] = L*((r[0]/L)%1)
				ir[1] += floor(r[1]/L)
				r[1] = L*((r[1]/L)%1)
				ir[2] += floor(r[2]/L)
				r[2] = L*((r[2]/L)%1)
				
				if n in dumpSteps:
					dumpArr[dumpSteps.index(n)] += ((L*ir[0]+r[0]-r0[0])**2 + (L*ir[1]+r[1]-r0[1])**2 + (L*ir[2]+r[2]-r0[2])**2)/M
	
	return dumpSteps, dumpArr

# dumpSteps = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576]
repeat = 1
totResult = []
for i in range(repeat):
	result = crw(500, 1048576, 14, "lamellae")
	totResult.append(result[1])

avg = np.mean(totResult, axis=0)   #avg = sum(x)/len(totResult) for x in zip(*totResult) ]
std = np.std(totResult, axis=0)
#print "   step          msd          std"
#for i in range(len(result[0])):
#	print "%7i %12.4f %12.4f" % (result[0][i], avg[i], std[i])

for i in range(len(result[0])):
    print avg[i]
