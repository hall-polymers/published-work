############################
# function to find cluster data with output hard coded for Connor and Janani's deformed system
# FIXME: should probably use numpy arrays instead of lists?
# loops to calculate number of bridges, loops/cloops, and atoms in each clustter within each frame
# data is now printed to csv file and viewed with excel
# Date: Feb 2016
# Author: Jon Brown
# Edited by: Connor Barber

#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
############################

import sys, math, string
import numpy as np
from pppmd import read_lammpstrj
from math import ceil, floor, sqrt

# general parameters
nconf = 1

cluster_dist = 0.85
cluster_dist2 = cluster_dist*cluster_dist
cluster_types = [3, 4]
force_mol_into_cluster = False
verbose_output = True

# frame parameters
num_frames = float('inf')
skip_beginning = 0
skip_between = 2

fname=str(sys.argv[1])

print "Reading in the dump file", fname, "...",
sys.stdout.flush()
r, ir, timestep, box_bounds, id2type, id2mol, mol2ids = read_lammpstrj(fname, num_frames, skip_beginning, skip_between)
print "Done!"

# files for writing
cluster_data = open ('cb_clusters_'+fname[:fname.rfind('.')]+'_1.5.csv','w')
eigen_data = open ('eigen_clusters_'+fname[:fname.rfind('.')]+'_1.5.csv', 'w')
data_for_fspan = open ('data_4_fspan_'+fname[:fname.rfind('.')]+'_1.5.csv', 'w')
loop_bridge = open ('loop_bridge_'+fname[:fname.rfind('.')]+'_1.5.csv', 'w')
LMP = open ('LMP_'+fname[:fname.rfind('.')]+'_1.5.lammpstrj', 'w')
ACF = open ('autocorr_'+fname[:fname.rfind('.')]+'_1.5.csv', 'w')

cluster2ids = []
id2cluster = []

# lists to write from to files
Rg_list = []
size_list = []

loops_list = []
cloops_list = []
bridges_list = []

eig1_list = []
eig2_list = []
eig3_list = []
aniso_list = []

## # OUTPUT

#helper function to compute basic statistics for later
def stats(X, N):
	if not N:
		return 0, 0
	mean = float(sum(X))/N
	if N-1 > 0:
		stdev = sqrt( sum([(x-mean)**2 for x in X]) / (N-1) )
	else:
		stdev = float('nan')
	return mean, stdev


    	
###AutoCorrelation Function###
def autocorrelation(timestep, frame, sflag):
    global clustervec1, frame_0, sigma_start

    # variables
    clustervec1 = np.zeros(len(timestep), float)
    sigma_current = 0.0  #numerator of ACF
    index = 0
    fracperc = 0.3

    # establish where ions are at the beginning
    if frame == 0:
        frame_0 = []
        sigma_start = 0.0    #denominator of ACF
        for ii in range(len(sflag)):
            for jj in range(ii+1, len(sflag)):
                if sflag[ii] == sflag[jj]:
                    sigma_start += 1
                    frame_0.append(1)
                else:
                    frame_0.append(0)

	#print "frame+0 = ", frame_0
        return 1

    else:
        # autocorrelation numerator loop for frame
        for ii in range(len(sflag)):
            for jj in range(ii+1, len(sflag)):
                if sflag[ii] == sflag[jj]:
		    if frame_0[index] == 1:
                    	sigma_current += 1
		index += 1

    # array for all the numerators aligned with their specific frame
    clustervec1[frame] = sigma_current
    autocorr = sigma_current / sigma_start
    
    return autocorr


ACF.write("Timestep, Autocorrelation\n")

for frame in xrange(len(timestep)):
	print
	print
	print "Running cluster anaylysis for frame", frame, "at timestep", timestep[frame]
	print
	print "Binning location data...",
	sys.stdout.flush()
    # read in new box boundaries
    #box_bounds = boundary(fname, frame)
	# break the box into bins
	Lx = box_bounds[frame][0][1] - box_bounds[frame][0][0]
	Ly = box_bounds[frame][1][1] - box_bounds[frame][1][0]
	Lz = box_bounds[frame][2][1] - box_bounds[frame][2][0]

	xbins = int(floor(Lx/cluster_dist))
	ybins = int(floor(Ly/cluster_dist))
	zbins = int(floor(Lz/cluster_dist))

	bin2id = [[[ [] for k in xrange(zbins)] for j in xrange(ybins)] for i in xrange(xbins)]
	id2bin = [ [] for n in xrange(len(r[frame])) ] 
	
	#print box_bounds[frame]
	#print Lx, Ly, Lz
	#print xbins, ybins, zbins

	# map the atom locations into those bins
	for atom_id in xrange(1,len(r[frame])):
	    #print atom_id, id2type[atom_id]
            if id2type[atom_id] in cluster_types:
                i = int(floor( xbins*( (r[frame][atom_id][0]-box_bounds[frame][0][0]) / Lx ) ))
                j = int(floor( ybins*( (r[frame][atom_id][1]-box_bounds[frame][1][0]) / Ly ) ))
                k = int(floor( zbins*( (r[frame][atom_id][2]-box_bounds[frame][2][0]) / Lz ) ))
		
		#print r[frame][atom_id]
		#print i, j, k

		if i >= xbins:
			i = i - xbins
			r[frame][atom_id][0] = r[frame][atom_id][0] - Lx
		elif i <= -1:
			i = i + xbins
			r[frame][atom_id][0] = r[frame][atom_id][0] + Lx

		if j >= ybins:
			j = j - ybins
			r[frame][atom_id][1] = r[frame][atom_id][1] - Ly
		elif j <= -1:
			j = j + ybins
			r[frame][atom_id][1] = r[frame][atom_id][1] + Ly

		if k >= zbins:
			k = k - zbins
			r[frame][atom_id][2] = r[frame][atom_id][2] - Lz
		elif k <= -1:
			k = k + zbins
			r[frame][atom_id][2] = r[frame][atom_id][2] + Lz

                bin2id[i][j][k].append(atom_id)
                id2bin[atom_id] = [i,j,k]
       

	print "Done!"
	# now do the clustering
	
	print "Building neighbor list...",
	sys.stdout.flush()
	#initialize cluster data so that each bead is it's own cluster
	id2neighbors = [ [] for n in xrange(len(r[frame])) ]

	# optionally add the relevant beads of each molecule to the neighbor lists
	if force_mol_into_cluster:
		for mol in xrange(1,len(mol2ids)):
			for atom_i in mol2ids[mol]:
				if id2type[atom_i] not in cluster_types:
					continue

				for atom_j in mol2ids[mol]:
					if (id2type[atom_j] not in cluster_types) or (atom_i == atom_j):
						continue

					id2neighbors[atom_i].append(atom_j)

	# loop over beads and their local area to build a neighbor list
	for atom_id in xrange(1,len(r[frame])):
		if id2type[atom_id] in cluster_types:
			#print atom_id, id2bin[atom_id], r[frame][atom_id]

			# loop over the local area of each bead, and shift at the periodic BCs
			for i in xrange(id2bin[atom_id][0]-1,id2bin[atom_id][0]+2):
				if i == -1:
					i = xbins - 1
					xshift = -Lx
				elif i == xbins:
					i = 0
					xshift = Lx
				else:
					xshift = 0

				for j in xrange(id2bin[atom_id][1]-1,id2bin[atom_id][1]+2):
					if j == -1:
						j = ybins - 1
						yshift = -Ly
					elif j == ybins:
						j = 0
						yshift = Ly
					else:
						yshift = 0

					for k in xrange(id2bin[atom_id][2]-1,id2bin[atom_id][2]+2):
						if k == -1:
							k = zbins - 1
							zshift = -Lz
						elif k == zbins:
							k = 0
							zshift = Lz
						else:
							zshift = 0

						#print bin2id[i][j][k]

						# loop over the beads in this box and calculate the distance
						for test_id in bin2id[i][j][k]:

							# if it's already in the nhood, don't bother
							if (test_id==atom_id) or (test_id in id2neighbors[atom_id]):
								continue

							dx = r[frame][test_id][0] + xshift - r[frame][atom_id][0]
							dy = r[frame][test_id][1] + yshift - r[frame][atom_id][1]
							dz = r[frame][test_id][2] + zshift - r[frame][atom_id][2]

							dr2 = dx*dx + dy*dy + dz*dz

							#print r[frame][test_id]
							#print dr2

							# if close enough add to the neighbor list
							if dr2 < cluster_dist2:
								id2neighbors[atom_id].append(test_id)
								id2neighbors[test_id].append(atom_id)

							#print id2neighbors[atom_id]
	print "Done!"

	# now consolidate neighbor lists into a cluster list
	print "Consolidating neighbor lists into clusters...",
	sys.stdout.flush()
	clusters  = [ [] ]
	for atom_id in xrange(1,len(id2type)):
		if id2type[atom_id] in cluster_types:
			# check if we've already done this atom
			already_done = False
			for cluster_id in xrange(1, len(clusters)):
				if atom_id in clusters[cluster_id]:
					already_done = True
					break
			if already_done:
				continue


			# if we haven't already done this atom, make it the start of a new cluster
			cluster_id = len(clusters)
			tobeadded = [atom_id]
			clusters.append(tobeadded)

			# the way this loops works is it goes through the neighbor lists adding them into the cluster, every time it finds a new atom not already in the cluster, it also checks that new atom's neighbor list for other new atoms, repeating this process until no new atoms are found 
			done = False
			while not done:
				done  = True
				adding = tobeadded
				tobeadded = []
				for add_id in adding:
					for neighbor_id in id2neighbors[add_id]:
						if neighbor_id not in clusters[cluster_id]:
							tobeadded.append(neighbor_id)
							clusters[cluster_id].append(neighbor_id)
							done = False
						#id2neighbors[add_id] = []
			
	# sort by length
	clusters.sort(key=len)
	print "Done!"
	
	#for x in xrange(1,len(clusters)):
	#    print "Cluster: ", x
	#    print clusters[x]

	# make the reverse lookup FIXME: not actually using this for now
	for atom_id in xrange(0, len(clusters[cluster_id]) + 1):			
	    reverse_clusters = [0]*len(id2type)
	    for atom_id in xrange(1,len(id2type)):
		if id2type[atom_id] in cluster_types:
		    for cluster_id in xrange(1,len(clusters)):
			if atom_id in clusters[cluster_id]:
			    reverse_clusters[atom_id] = cluster_id
			    break
	 
	print len(reverse_clusters)
	# append to the frame by frame cluster data
	cluster2ids.append(clusters)
	id2cluster.append(reverse_clusters)
	#print cluster2ids
	#print id2cluster

    # init sflag, num, and kofion
	sflag = []
	kofion = []
	for atom_id in xrange(1, len(r[frame])):
	    if id2type[atom_id] in cluster_types:
		kofion.append(atom_id)    # lists each ion in system
		for i in xrange(1, len(clusters)):
		    if atom_id in clusters[i]:
	    		sflag.append(i)    # lists each ion, characterized by a number representing the cluster they are in (in order of lowest to highest cluster number)

	# list of clusters characterized by their size (number of ions in each)
	num = np.zeros(len(clusters), np.int)
	for i in xrange(1, len(clusters)):
	    num[i] = len(clusters[i])

	# init lists of coordinates of all ions in groups by axes
	xion = []
	yion = []
	zion = []
	for atom_id in kofion:
	    xion.append(r[frame][atom_id][0])
            yion.append(r[frame][atom_id][1])
            zion.append(r[frame][atom_id][2])
	
	# init variables for inertia tensor  -  MUST WRAP CLUSTER INTO SAME BOX FIRST
	ixx = np.zeros(len(clusters)+1, float)
	iyy = np.zeros(len(clusters)+1, float)
        izz = np.zeros(len(clusters)+1, float)
        ixy = np.zeros(len(clusters)+1, float)
        ixz = np.zeros(len(clusters)+1, float)
        iyz = np.zeros(len(clusters)+1, float)
	xcm = np.zeros(len(clusters)+1, float)
        ycm = np.zeros(len(clusters)+1, float)
        zcm = np.zeros(len(clusters)+1, float)    	
        xcmwrap = np.zeros(len(clusters)+1, float)
        ycmwrap = np.zeros(len(clusters)+1, float)
        zcmwrap = np.zeros(len(clusters)+1, float)
    
    	for i in xrange(1, len(kofion)):
            cl=sflag[i]
            xcm[cl]=xcm[cl]+xion[i]
            ycm[cl]=ycm[cl]+yion[i]
            zcm[cl]=zcm[cl]+zion[i]
        
    	for i in xrange(1, len(clusters)):
            xcm[i]=xcm[i]/float(num[i])
            ycm[i]=ycm[i]/float(num[i])
            zcm[i]=zcm[i]/float(num[i])
            xcmwrap[i]=round(xcm[i]/Lx)*Lx
            ycmwrap[i]=round(ycm[i]/Ly)*Ly
            zcmwrap[i]=round(zcm[i]/Lz)*Lz
            xcm[i]=xcm[i]-xcmwrap[i]
            ycm[i]=ycm[i]-ycmwrap[i]
            zcm[i]=zcm[i]-zcmwrap[i]

    	for i in xrange(1, len(kofion)):
            cl=sflag[i]
            xion[i]-=xcmwrap[cl]
            yion[i]-=ycmwrap[cl]
            zion[i]-=zcmwrap[cl]
                
	# init max/min values
	xmax = np.zeros(len(clusters), float)
	xmin = np.zeros(len(clusters), float)
        ymax = np.zeros(len(clusters), float)
        ymin = np.zeros(len(clusters), float)
        zmax = np.zeros(len(clusters), float)
        zmin = np.zeros(len(clusters), float)

	# determining minimum and maximum coordinate values for each cluster
        for i in xrange(1, len(kofion)):
            cl=sflag[i]
            if xmin[cl] == 0.0:
                xmin[cl]=xion[i]
                xmax[cl]=xion[i]
            elif xion[i]<xmin[cl]:
                xmin[cl]=xion[i]
            elif xion[i]>xmax[cl]:
                xmax[cl]=xion[i]
            if ymin[cl]==0.0:
                ymin[cl]=yion[i]
                ymax[cl]=yion[i]
            elif yion[i]<ymin[cl]:
                ymin[cl]=yion[i]
            elif yion[i]>ymax[cl]:
                ymax[cl]=yion[i]
            if zmin[cl]==0.0:
                zmin[cl]=zion[i]
                zmax[cl]=zion[i]
            elif zion[i]<zmin[cl]:
                zmin[cl]=zion[i]
            elif zion[i]>zmax[cl]:
                zmax[cl]=zion[i]

            #assume all masses are 1, otherwise have to multiply by mass here:
            ixx[cl]=ixx[cl]+(yion[i]-ycm[cl])**2+(zion[i]-zcm[cl])**2
            iyy[cl]=iyy[cl]+(xion[i]-xcm[cl])**2+(zion[i]-zcm[cl])**2
            izz[cl]=izz[cl]+(xion[i]-xcm[cl])**2+(yion[i]-ycm[cl])**2
            ixy[cl]=ixy[cl]-(xion[i]-xcm[cl])*(yion[i]-ycm[cl])
            ixz[cl]=ixz[cl]-(xion[i]-xcm[cl])*(zion[i]-zcm[cl])
            iyz[cl]=iyz[cl]-(yion[i]-ycm[cl])*(zion[i]-zcm[cl])

    	# anisotropy init
    	perc = 0
	aniso_2 = 0
	aniso = []
	eig1_tot = 0
	eig2_tot = 0
	eig3_tot = 0
	aniso_tot = 0
	aniso_count = 0

    	data_for_fspan.write( "Step, %i, number of clusters, %i\n" % (timestep[frame], len(clusters)-1))
	data_for_fspan.write( "Xbox, %2f, Ybox, %2f, Zbox, %2f\n" % (Lx, Ly, Lz)) 
    	data_for_fspan.write( "num ions in cluster, xcm, ycm, zcm, three principal moments of inertia, , ,figure axis, , ,long axis, , ,xspan,yspan,zspan\n" )
    
	eigen_data.write( "Timestep, %i, Number of Clusters, %i\n" % (timestep[frame], len(clusters)-1))
	eigen_data.write( "Cluster ID, eig_1, eig_2, eig_3, Anisotropy\n")
	# caculate eigenvalues for anisotropy
	for i in xrange(1,len(clusters)):
            if num[i] == 1:
                eig1 = 0
                eig2 = 0
                eig3 = 0
                faxis = np.array([0.0,0.0,0.0])
                longaxis = np.array([0.0,0.0,0.0])
            elif xmax[i]-xmin[i] > Lx or ymax[i]-ymin[i] > Ly or zmax[i]-zmin[i] > Lz:
                perc = perc+1
                eig1 = nan
                eig2 = nan
                eig3 = nan
                faxis = np.array([nan,nan,nan])
                longaxis = np.array([nan,nan,nan])
            else:
                eigenvalues, eigenvectors = np.linalg.eig([[ixx[i],ixy[i],ixz[i]],[ixy[i],iyy[i],iyz[i]],[ixz[i],iyz[i],izz[i]]])
	    #from linalg.eig help: column ``v[:,i]`` is the eigenvector corresponding to the eigenvalue ``w[i]``
                eig1 = max(eigenvalues)
                eig3 = min(eigenvalues)
                faxis = np.array([0.0,0.0,0.0])
                longaxis = np.array([0.0,0.0,0.0])
                if (eigenvalues[0]!=eig1 and eigenvalues[0]!=eig3):
                    eig2 = eigenvalues[0]
                elif (eigenvalues[1]!=eig1 and eigenvalues[1]!=eig3):
                    eig2 = eigenvalues[1]
                else:
                    eig2 = eigenvalues[2]
                if (eigenvalues[0] == eig1):
                    faxis = eigenvectors[:,0]
                elif (eigenvalues[1] == eig1):
                    faxis = eigenvectors[:,1]
                else:
                    faxis = eigenvectors[:,2]
                if (eigenvalues[0] == eig3):
                    longaxis = eigenvectors[:,0]
                elif (eigenvalues[1] == eig3):
                    longaxis = eigenvectors[:,1]
                else:
                    longaxis = eigenvectors[:,2]
	    # anisotropy function
	    if eig1 != 0 and eig2 != 0 and eig3 != 0:
		aniso_2 = 4 - 12*((eig1*eig2 + eig1*eig3 + eig2*eig3) / ((eig1 + eig2 + eig3)*(eig1 + eig2 + eig3)))
	   	aniso.append(aniso_2)
	        aniso_tot = aniso_tot + aniso_2
		aniso_count += 1

	    eig1_tot = eig1_tot + eig1
	    eig2_tot = eig2_tot + eig2
	    eig3_tot = eig3_tot + eig3

            data_for_fspan.write( "%i, %7.5f, %7.5f, %7.5f, %7.5f, %7.5f, %7.5f, %7.5f, %7.5f, %7.5f, %7.5f, %7.5f, %7.5f, %7.5f, %7.5f, %7.5f\n" % (num[i], xcm[i], ycm[i], zcm[i], eig1, eig2, eig3, faxis[0], faxis[1], faxis[2], longaxis[0] ,longaxis[1] ,longaxis[2] ,xmax[i]-xmin[i] ,ymax[i]-ymin[i] ,zmax[i]-zmin[i]) )
	
	    eigen_data.write( "%i, %f, %f, %f, %f\n" % (i, eig1, eig2, eig3, aniso_2))

	# average eigenvalues
	eig1_avg = eig1_tot / (len(clusters) - 1)
	eig1_list.append(eig1_avg)
	eig2_avg = eig2_tot / (len(clusters) - 1)
	eig2_list.append(eig2_avg)
	eig3_avg = eig3_tot / (len(clusters) - 1)
	eig3_list.append(eig3_avg)

	# calculating average anisotropy
	if aniso_count != 0:
	    aniso_avg = aniso_tot / (aniso_count)
	    aniso_list.append(aniso_avg)

        ## print "Reading cluster data file...."
        ## fspan(clusters, nconf, timestep)

        # init loops, close loops, and bridges counters
	numloops = 0
	numcloops = 0
	numbridges = 0
	loop = np.zeros(len(r[frame]))
	bridge = np.zeros(len(r[frame]))
	cloop = np.zeros(len(r[frame]), int)

	loop_bridge.write("Timestep, %i, Ions, %i\n" % (timestep[frame], len(kofion)))
	loop_bridge.write("Ion, Loop, Close Loop, Bridge\n")
 
	# function to calculate total number of loops, close loops, and bridges
	# neutral atoms will be given classifications of whether they are part of a loop or bridge; ions will have no specific classification
	for i in xrange(1, len(kofion)):
      	    if(id2type[kofion[i]] == 3):
                #find loops/bridges, assumes the atoms and ion lists are in order by connectivity along backbone
        	if id2mol[kofion[i]] == id2mol[kofion[i+1]]:
            	    if sflag[i] == sflag[i+1]:
                	numloops += 1
                	for k in xrange(kofion[i]+1,kofion[i+1]):
                    	    loop[k] = 1

        		# determining distance for classification of a close loop
	        	dx = r[frame][kofion[i]][0] - r[frame][kofion[i+1]][0]
                	dx = dx - round(dx/Lx)*Lx
                	dy = r[frame][kofion[i]][1] - r[frame][kofion[i+1]][1]
                	dy = dy - round(dy/Ly)*Ly
                	dz = r[frame][kofion[i]][2] - r[frame][kofion[i+1]][2]
                	dz = dz - round(dz/Lz)*Lz

                	if (dx*dx+dy*dy+dz*dz) > (4*cluster_dist2):
                    	    numcloops += 1
                    	    for k in xrange(kofion[i]+1,kofion[i+1]):
                                cloop[k] = 2
		    # if not a loop, then a bridge
            	    else:
                	numbridges += 1
                	for k in xrange(kofion[i]+1, kofion[i+1]):
                    	    bridge[k] = 3

	# for the neutral atoms before the first ion
	first_ion = kofion[0]
	last_ion = kofion[-1]
	if(id2type[kofion[0]] in cluster_types):
            #find loops/bridges, assumes the atoms and ion lists are in order by connectivity along backbone
            if id2mol[kofion[0]] == id2mol[kofion[1]]:
                if sflag[0] == sflag[1]:
		    for k in xrange(0, first_ion):
			loop[k] = 1
		else:
		    for k in xrange(0, first_ion):
			bridge[k] = 3
	# for the neutral atoms after the last ion
	if(id2type[kofion[-1]] in cluster_types):
            #find loops/bridges, assumes the atoms and ion lists are in order by connectivity along backbone
            if id2mol[kofion[-2]] == id2mol[kofion[-1]]:
                if sflag[-2] == sflag[-1]:                        
                    for k in xrange(last_ion+1, len(r[frame])):
                        loop[k] = 1
		    loop[-1] = 1
                else:
                    for k in xrange(last_ion+1, len(r[frame])):
                        bridge[k] = 3
		    bridge[-1] = 3


	loops_list.append(numloops)
	cloops_list.append(numcloops)
	bridges_list.append(numbridges)

	# init charge array
	qq = np.zeros(len(r[frame])) 

	# charge of each atom
	for atom in xrange(1, len(r[frame])):
	    if id2type[atom] == 3:
	        qq[atom] = -1
	    elif id2type[atom] == 4:
                qq[atom] = 1
	    else:
	        qq[atom] = 0


	## DATA for FSPAN ##
	# data will be extracted to make fspan calculations using 'cluster_fspan.py'
	index_largest = len(clusters)-1
	avg_ion = 0.0
	avg_ion2_el = 0.0
	for i in xrange(1, len(clusters)):
	    if (i != index_largest):
		avg_ion += num[i]
		avg_ion2_el += ((num[i])*(num[i]))

	if (len(clusters)-1 != 1):
	    avg_ion = avg_ion / (len(clusters)-2)
	    avg_ion2_el = avg_ion2_el / (len(clusters)-2)
	else:
	    avg_ion = 0
	    avg_ion2_el = 0

	avg_ion_ac = 0.0
	avg_ion2 = 0.0
	for i in xrange(1, len(clusters)):
	    avg_ion_ac += num[i]
	    avg_ion2 += ((num[i])*(num[i]))

	avg_ion_ac = avg_ion_ac / (len(clusters)-1)
	avg_ion2 = avg_ion2 / (len(clusters)-1)

    	data_for_fspan.write("Clusters crossing box, Loops, Bridges, Close Loops, Avg ions in cluster, Avg ions squared except largest, Avg ions in any cluster, Avg ions squared\n")
	data_for_fspan.write("%i, %i, %i, %i, %7.5f, %7.5f, %7.5f, %7.5f\n" % (perc, numloops, numbridges, numcloops, avg_ion, avg_ion2_el, avg_ion_ac, avg_ion2))

	## AUTOCORRELATION FUNCTION ##
        autocorr = autocorrelation(timestep, frame, sflag)
	print autocorr
        ACF.write( "%i, %4.6f\n" % (timestep[frame], autocorr))

	## PAIR POTENTIAL ##
	# calculate number of bridges/loops of a single chain to determine partition function
	# partition function will be used to ultimately calculate effective interparticle pair potential
	##
	##
	##
	
	## ORIENTATION VECTOR ##
	
	
	# Dot Product #
	#theta = angle between two 'chord' vectors
	
	
	

	# verbose output if the flag is set, makes a new file for each frame with cluster by cluster information 
	if verbose_output:
		print "Verbose output...",
		sys.stdout.flush()
		cluster_data.write( "Timestep, %i, Number of Clusters, %i\n" % (timestep[frame], len(clusters)-1))
		cluster_data.write( "Loops, %i, Close Loops, %i, Bridges, %i\n" % (numloops, numcloops, numbridges))
		cluster_data.write( "Cluster ID, Type 3, Type 4, Type 6, Rg, Size\n" )
		## f_out_frame = open(fname[:fname.rfind('.')] + "_clusters_step"+str(timestep[frame])+".txt",'w')
		## f_out_frame.write("cluster_id  type_3  type_4      Rg             size          loops  close_loops  bridges\n")

		Rg_tot = 0
		size_tot = 0
		for cluster_id in xrange(1,len(clusters)): 
			# init atom types
			type_3 = 0
			type_4 = 0
			type_6 = 0

			# init loop/bridge counters
			n_loop = 0
			n_cloop = 0
			n_bridge = 0

			# caculate type_3 and type_4 atoms in cluster
			for atom_id in clusters[cluster_id]:
				if id2type[atom_id] == 3:
					type_3 += 1	
				elif id2type[atom_id] == 4:
					type_4 += 1
				elif id2type[atom_id] == 6:
                                        type_6 += 1
		
			#first all the PEG beads into the clusters
			for mol in set( [ id2mol[atom_id] for atom_id in clusters[cluster_id] ] ):
				for atom_in_mol in mol2ids[mol]:
					if id2type[atom_in_mol] not in cluster_types:
						clusters[cluster_id].append(atom_in_mol)
						
			# now get unwrapped coordinates for the entire cluster
			x = [float('inf')]*len(clusters[cluster_id])
			y = [float('inf')]*len(clusters[cluster_id])
			z = [float('inf')]*len(clusters[cluster_id])
			# start by placing the first atom in the box
			x[0] = r[frame][clusters[cluster_id][0]][0]
			y[0] = r[frame][clusters[cluster_id][0]][1]
			z[0] = r[frame][clusters[cluster_id][0]][2]

			# now loop through the cluster and put them in the right box... FIXME: assumes we have no box spanning clusters and the size of the cluster is much smaller than the size of the box
			for n in xrange(1,len(clusters[cluster_id])):
				xw = r[frame][clusters[cluster_id][n]][0]
				yw = r[frame][clusters[cluster_id][n]][1]
				zw = r[frame][clusters[cluster_id][n]][2]
				
				dr2min = float('inf')
				for ix in xrange(-1,2):
					for iy in xrange(-1,2):
						for iz in xrange(-1,2):
							dx = xw + ix*Lx - x[0] 
							dy = yw + iy*Ly - y[0] 
							dz = zw + iz*Lz - z[0] 
			
							dr2 = dx**2 + dy**2 + dz**2

							if dr2 < dr2min:
								x_new = xw + ix*Lx
								y_new = yw + iy*Ly
								z_new = zw + iz*Lz
								dr2min = dr2

				x[n] = x_new
				y[n] = y_new
				z[n] = z_new

			# calculate Rg
			x_avg = sum(x)/len(clusters[cluster_id])
			y_avg = sum(y)/len(clusters[cluster_id])
			z_avg = sum(z)/len(clusters[cluster_id])
			Rg = sqrt( sum( [ (x[n]-x_avg)**2 + (y[n]-y_avg)**2 + (z[n]-z_avg)**2 for n in xrange(len(clusters[cluster_id])) ] )/len(clusters[cluster_id]) )
			
			Rg_tot += Rg

			# calculate size
			size = sum([max(x)-min(x), max(y)-min(y), max(z)-min(z)])/3
			size_tot += size

			cluster_data.write( " %i, %i, %i, %i, %3f, %3f\n" % (cluster_id, type_3, type_4, type_6, Rg, size))
			## f_out_frame.write("    {}         {}       {}    {}    {}      {}        {}          {}\n".format(cluster_id, type_3, type_4, Rg, size, numloops, numcloops, numbridges))

		## f_out_frame.close()

	###############################
	### PRINTING FORMAT TO READ ###
	###############################
	LMP.write("ITEM: TIMESTEP\n")
	LMP.write("%i\n" % (timestep[frame]))
	LMP.write("ITEM: NUMBER OF ATOMS\n")
        LMP.write("%i\n" % (len(r[frame])-1))
	LMP.write("ITEM: BOX BOUNDS pp pp pp\n")
	LMP.write("%4.4f %4.4f\n" % (box_bounds[frame][0][0], box_bounds[frame][0][1]))
        LMP.write("%4.4f %4.4f\n" % (box_bounds[frame][1][0], box_bounds[frame][1][1]))
        LMP.write("%4.4f %4.4f\n" % (box_bounds[frame][2][0], box_bounds[frame][2][1]))
	LMP.write("ITEM: ATOMS id mol type q xs ys zs ix iy iz vx vy vz\n")
	
	# Information for all the atoms
        for k in range(1, len(r[frame])):
	#scale to get original format
	    r[frame][k][0] = (r[frame][k][0]-box_bounds[frame][0][0])/(box_bounds[frame][0][1]-box_bounds[frame][0][0])
	    r[frame][k][1] = (r[frame][k][1]-box_bounds[frame][1][0])/(box_bounds[frame][1][1]-box_bounds[frame][1][0])
	    r[frame][k][2] = (r[frame][k][2]-box_bounds[frame][2][0])/(box_bounds[frame][2][1]-box_bounds[frame][2][0])

        #for poly vx  will give cluster number each atom is (if) in, vy gives whether the atom is part of a loop(1), close loop(2), or bridge(3) 
	    #if part of a close loop, print cloop[k]
	    if cloop[k] != 0:
		LMP.write( "%i %i %i %i %4.7f %4.7f %4.7f %i %i %i %i %i %i\n" % (k,id2mol[k],id2type[k],qq[k],r[frame][k][0],r[frame][k][1],r[frame][k][2],ir[frame][k][0],ir[frame][k][1],ir[frame][k][2],reverse_clusters[k],cloop[k],0))

	    #if part of a loop, print loop[k]
	    elif loop[k] != 0:
		LMP.write( "%i %i %i %i %4.7f %4.7f %4.7f %i %i %i %i %i %i\n" % (k,id2mol[k],id2type[k],qq[k],r[frame][k][0],r[frame][k][1],r[frame][k][2],ir[frame][k][0],ir[frame][k][1],ir[frame][k][2],reverse_clusters[k],loop[k],0))

	    #otherwise, print bridge[k]
	    else:
		LMP.write( "%i %i %i %i %4.7f %4.7f %4.7f %i %i %i %i %i %i\n" % (k,id2mol[k],id2type[k],qq[k],r[frame][k][0],r[frame][k][1],r[frame][k][2],ir[frame][k][0],ir[frame][k][1],ir[frame][k][2],reverse_clusters[k],bridge[k],0))



	# Unscale back
	for my_id in xrange(1,len(r[frame])):
	    r[frame][my_id][0] = r[frame][my_id][0]*(box_bounds[frame][0][1]-box_bounds[frame][0][0]) + box_bounds[frame][0][0]
            r[frame][my_id][1] = r[frame][my_id][1]*(box_bounds[frame][1][1]-box_bounds[frame][1][0]) + box_bounds[frame][1][0]
            r[frame][my_id][2] = r[frame][my_id][2]*(box_bounds[frame][2][1]-box_bounds[frame][2][0]) + box_bounds[frame][2][0]


	# Add Rg and size averages to list for all timesteps	
	Rg_avg = Rg_tot / (len(clusters)-1)
	Rg_list.append(Rg_avg)
	size_avg = size_tot / (len(clusters)-1)
	size_list.append(size_avg)

	cluster_data.write("Rg Avg, %2f, Size Avg, %2f\n\n" % (Rg_avg, size_avg))
	print "Done!"
# Finished gathering cluster data

# Write cluster and eigen data to CSV files
cluster_data.write("Timesteps,  Rg Avgs, Size Avgs, Loops, Close Loops, Bridges\n")
eigen_data.write("Timestep, Eig1 Avg, Eig2 Avg, Eig3 Avg, Aniso Avg\n")
for i in xrange(len(timestep)):
	cluster_data.write("%i, %2f, %2f, %i, %i, %i\n" % (timestep[i], Rg_list[i], size_list[i], loops_list[i], cloops_list[i], bridges_list[i]))
	eigen_data.write("%i, %2f, %2f, %2f, %2f\n" % (timestep[i], eig1_list[i], eig2_list[i], eig3_list[i], aniso_list[i]))


# Finish by closing out files written to
loop_bridge.close()
data_for_fspan.close()
eigen_data.close()
cluster_data.close()
LMP.close()
ACF.close()
## f_out.close()
