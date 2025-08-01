#!/bin/bash

module load python/3.12
# Python 3 is required to properly run make_data_PGNs_implicit_solvent_fixed_distance.py.

. ./params.sh
radius_nanoparticles=`echo $radius_nanoparticles`
# Value of radius_nanoparticles for all jobs to be submitted and run specified in and obtained from params.sh.
length_grafted_chains=`echo $length_grafted_chains`
# Value of length_grafted_chains for all jobs to be submitted and run specified in and obtained from params.sh.
target_grafting_density=`echo $target_grafting_density`
# Value of target_grafting_density for all jobs to be submitted and run specified in and obtained from params.sh.
epsilon_mon_mon_attractive=`echo $epsilon_mon_mon_attractive`
# Value of epsilon_mon_mon_attractive for all jobs to be submitted and run specified in and obtained from params.sh.

home=/fs/scratch/PAS2029/felipepacci/Projects/Solvated_Polymer_Grafted_Nanoparticles/Molecular_Dynamics/Implicit_Solvent/Fixed_Distance
cd $home/

# Importing of the value of radius_nanoparticles specified in params.sh into radius_nanoparticles.txt.
grep -m1 'radius_nanoparticles=' params.sh | sed -e 's#.*=\(\)#\1#' > radius_nanoparticles.txt

# Creation of a parent directory for all the simulations associated with the set of parameters specified in params.sh.
directory=./${radius_nanoparticles}
if test -d "$directory"
then
    cd ${radius_nanoparticles}/
else
    mkdir ${radius_nanoparticles}
    cd ${radius_nanoparticles}/
fi
directory=./${length_grafted_chains}
if test -d "$directory"
then
    cd ${length_grafted_chains}/
else
    mkdir ${length_grafted_chains}
    cd ${length_grafted_chains}/
fi
directory=./${target_grafting_density}
if test -d "$directory"
then
    cd ${target_grafting_density}/
else
    mkdir ${target_grafting_density}
    cd ${target_grafting_density}/
fi
directory=./${epsilon_mon_mon_attractive}
if test -d "$directory"
then
    echo "The simulations associated with the set of parameters specified in params.sh have already been performed."
    exit 1
else
    mkdir ${epsilon_mon_mon_attractive}
    cd ${epsilon_mon_mon_attractive}/
    find $home/ -maxdepth 1 -type f | xargs mv -t ./
fi

min_distance_monomer_grafted_nanoparticles_cell=$(echo "2. * $radius_nanoparticles" | bc)
# Minimum possible distance between monomer-grafted nanoparticles in a unit cell.
min_sampled_distance_monomer_grafted_nanoparticles_cell=$(echo "$min_distance_monomer_grafted_nanoparticles_cell + 3." | bc)
min_sampled_distance_monomer_grafted_nanoparticles_cell=$(printf "%.6f" $min_sampled_distance_monomer_grafted_nanoparticles_cell)
# Minimum sampled distance between monomer-grafted nanoparticles in a unit cell set as a 6-decimal place number.
max_sampled_distance_monomer_grafted_nanoparticles_cell=$(echo "$min_sampled_distance_monomer_grafted_nanoparticles_cell + .9 * $length_grafted_chains" | bc)
# The formula above is empirically determined to safely estimate the shortest distance between monomer-grafted nanoparticles at which the force between the polymer-grafted nanoparticles is effectively zero.
max_sampled_distance_monomer_grafted_nanoparticles_cell=$(printf "%.6f" $max_sampled_distance_monomer_grafted_nanoparticles_cell)
# Maximum sampled distance between monomer-grafted nanoparticles in a unit cell set as a 6-decimal place number.
interval_sampled_distance_monomer_grafted_nanoparticles_cell=2.
# Interval between sampled distances between monomer-grafted nanoparticles in a unit cell.

home=/fs/scratch/PAS2029/felipepacci/Projects/Solvated_Polymer_Grafted_Nanoparticles/Molecular_Dynamics/Implicit_Solvent/Fixed_Distance/${radius_nanoparticles}/${length_grafted_chains}/${target_grafting_density}/${epsilon_mon_mon_attractive}
cd $home/

# Deletion of the file distances_monomer_grafted_nanoparticles.txt if it already exists.
file=./distances_monomer_grafted_nanoparticles.txt
if test -f "$file"
then
    rm $file
fi
# Creation of the file distances_monomer_grafted_nanoparticles.txt containing all the distances between monomer-grafted nanoparticles to be sampled, from min_sampled_distance_monomer_grafted_nanoparticles_cell to max_sampled_distance_monomer_grafted_nanoparticles_cell, inclusive, in intervals of interval_sampled_distance_monomer_grafted_nanoparticles_cell σ, one sampled distance per line.
distance_monomer_grafted_nanoparticles=${min_sampled_distance_monomer_grafted_nanoparticles_cell}
echo "$distance_monomer_grafted_nanoparticles" >> distances_monomer_grafted_nanoparticles.txt
while (( $(echo "$distance_monomer_grafted_nanoparticles <= $max_sampled_distance_monomer_grafted_nanoparticles_cell - $interval_sampled_distance_monomer_grafted_nanoparticles_cell" | bc -l) ))
do
    distance_monomer_grafted_nanoparticles=$(echo "$distance_monomer_grafted_nanoparticles + $interval_sampled_distance_monomer_grafted_nanoparticles_cell" | bc)
    echo "$distance_monomer_grafted_nanoparticles" >> distances_monomer_grafted_nanoparticles.txt
done
# Importing of the distances between monomer-grafted nanoparticles listed in distances_monomer_grafted_nanoparticles.txt into distances_monomer_grafted_nanoparticles, one distance per entry.
mapfile -t distances_monomer_grafted_nanoparticles < distances_monomer_grafted_nanoparticles.txt

name='PGNs_implicit_solvent_fixed_distance_soft_pushoff_'${radius_nanoparticles}'_'${length_grafted_chains}'_'${target_grafting_density}'_'${epsilon_mon_mon_attractive}
# Creation of directories, set up of a Soft Pushoff run, and submission of a Soft Pushoff job for every distance between monomer-grafted nanoparticles to be sampled.
for index in ${!distances_monomer_grafted_nanoparticles[@]}
do
    mkdir ${distances_monomer_grafted_nanoparticles[index]}
    cd ${distances_monomer_grafted_nanoparticles[index]}/
    mkdir Soft_Pushoff
    mkdir Equilibration_NVT
    cd Soft_Pushoff/
    cp $home/{time.txt,radius_nanoparticles.txt,make_data_PGNs_implicit_solvent_fixed_distance.py,job_submission_PGNs_implicit_solvent_fixed_distance_soft_pushoff_cardinal.sh} ./
    file=./ctrl.py
    if test -f "$file"
    then
	rm $file
    fi
    file=./params.py
    if test -f "$file"
    then
	rm $file
    fi
    echo "current_seed = ${distances_monomer_grafted_nanoparticles[index]}" > ctrl.py
    echo "radius_nanoparticles = ${radius_nanoparticles}" > params.py
    echo "# Radius of a nanoparticle." >> params.py
    echo "length_grafted_chains = ${length_grafted_chains}" >> params.py
    echo "# Length of a grafted chain (number of monomers/grafted chain)." >> params.py
    echo "target_grafting_density = ${target_grafting_density}" >> params.py
    echo "# Target grafting density." >> params.py
    echo "distance_monomer_grafted_nanoparticles_cell = ${distances_monomer_grafted_nanoparticles[index]}" >> params.py
    echo "# Distance between monomer-grafted nanoparticles in a unit cell." >> params.py
    python make_data_PGNs_implicit_solvent_fixed_distance.py
    if (( $(echo "$radius_nanoparticles == 2.500000" | bc -l) ))
    then
	number_tasks=92
    fi
    if (( $(echo "$radius_nanoparticles == 5.000000" | bc -l) ))
    then
	number_tasks=32
    fi
    if (( $(echo "$radius_nanoparticles == 7.500000" | bc -l) ))
    then
	number_tasks=44
    fi
    if (( $(echo "$radius_nanoparticles == 10.000000" | bc -l) ))
    then
	number_tasks=68
    fi
    sbatch --nodes=1 --ntasks=$number_tasks --partition=cpu --job-name=${name}_${distances_monomer_grafted_nanoparticles[index]} --output=${name}_${distances_monomer_grafted_nanoparticles[index]}.out --error=${name}_${distances_monomer_grafted_nanoparticles[index]}.err job_submission_PGNs_implicit_solvent_fixed_distance_soft_pushoff_cardinal.sh
    cd $home/
done
