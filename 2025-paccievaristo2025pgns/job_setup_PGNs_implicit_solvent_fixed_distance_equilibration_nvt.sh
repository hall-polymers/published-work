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

# Importing of the value of epsilon_mon_mon_attractive specified in params.sh into eps.txt.
grep -m1 'epsilon_mon_mon_attractive=' params.sh | sed -e 's#.*=\(\)#\1#' > eps.txt
# Importing of the distances between monomer-grafted nanoparticles listed in distances_monomer_grafted_nanoparticles.txt into distances_monomer_grafted_nanoparticles, one distance per entry.
mapfile -t distances_monomer_grafted_nanoparticles < distances_monomer_grafted_nanoparticles.txt

name='PGNs_implicit_solvent_fixed_distance_equilibration_nvt_'${radius_nanoparticles}'_'${length_grafted_chains}'_'${target_grafting_density}'_'${epsilon_mon_mon_attractive}
# Set up of an NVT Equilibration run and submission of an NVT Equilibration job for every distance between monomer-grafted nanoparticles to be sampled.
for index in ${!distances_monomer_grafted_nanoparticles[@]}
do
    cd ${distances_monomer_grafted_nanoparticles[index]}/
    cd Equilibration_NVT/
    cp $home/${distances_monomer_grafted_nanoparticles[index]}/Soft_Pushoff/{radius_nanoparticles.txt,*.restart.*} ./
    cp $home/{time.txt,eps.txt,job_submission_PGNs_implicit_solvent_fixed_distance_equilibration_nvt_cardinal.sh} ./
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
    sbatch --nodes=1 --ntasks=$number_tasks --partition=cpu --job-name=${name}_${distances_monomer_grafted_nanoparticles[index]} --output=${name}_${distances_monomer_grafted_nanoparticles[index]}.out --error=${name}_${distances_monomer_grafted_nanoparticles[index]}.err job_submission_PGNs_implicit_solvent_fixed_distance_equilibration_nvt_cardinal.sh
    cd $home/
done
