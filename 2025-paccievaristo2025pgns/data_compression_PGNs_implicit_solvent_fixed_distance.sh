#!/bin/bash

. ./params.sh
radius_nanoparticles=`echo $radius_nanoparticles`
# Value of radius_nanoparticles for all jobs to be submitted and run specified in and obtained from params.sh.
length_grafted_chains=`echo $length_grafted_chains`
# Value of length_grafted_chains for all jobs to be submitted and run specified in and obtained from params.sh.
target_grafting_density=`echo $target_grafting_density`
# Value of target_grafting_density for all jobs to be submitted and run specified in and obtained from params.sh.
epsilon_mon_mon_attractive=`echo $epsilon_mon_mon_attractive`
# Value of epsilon_mon_mon_attractive for all jobs to be submitted and run specified in and obtained from params.sh.

home=/fs/scratch/PAS2029/felipepacci/Projects/Solvated_Polymer_Grafted_Nanoparticles/Molecular_Dynamics/Implicit_Solvent/Fixed_Distance/${radius_nanoparticles}/${length_grafted_chains}/${target_grafting_density}/${epsilon_mon_mon_attractive}
cd $home/

# Importing of the distances between monomer-grafted nanoparticles listed in distances_monomer_grafted_nanoparticles.txt into distances_monomer_grafted_nanoparticles, one distance per entry.
mapfile -t distances_monomer_grafted_nanoparticles < distances_monomer_grafted_nanoparticles.txt

# Compression of all LAMMPS trajectory files.
for index in ${!distances_monomer_grafted_nanoparticles[@]}
do
    for file in ${home}/${distances_monomer_grafted_nanoparticles[index]}/{Soft_Pushoff,Equilibration_NVT}/*.lammpstrj
    do
	zstd -15 "${file}"
    done
done
