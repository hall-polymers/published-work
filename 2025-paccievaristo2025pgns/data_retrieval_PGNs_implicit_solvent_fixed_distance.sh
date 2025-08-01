#!/opt/homebrew/Cellar/bash/5.2.37/bin/bash
# The shebang must point to an installation of Bash 4.x or higher for the script to run.

. ./params.sh
radius_nanoparticles=`echo $radius_nanoparticles`
# Value of radius_nanoparticles specified in and obtained from params.sh.
length_grafted_chains=`echo $length_grafted_chains`
# Value of length_grafted_chains specified in and obtained from params.sh.
target_grafting_density=`echo $target_grafting_density`
# Value of target_grafting_density specified in and obtained from params.sh.
epsilon_mon_mon_attractive=`echo $epsilon_mon_mon_attractive`
# Value of epsilon_mon_mon_attractive specified in and obtained from params.sh.

origin=/fs/scratch/PAS2029/felipepacci/Projects/Solvated_Polymer_Grafted_Nanoparticles/Molecular_Dynamics/Implicit_Solvent/Fixed_Distance/${radius_nanoparticles}/${length_grafted_chains}/${target_grafting_density}/${epsilon_mon_mon_attractive}
destination=/Users/felipe.pace.evaristo/Library/CloudStorage/Dropbox/Mac/Desktop/Projects/Solvated\ Polymer-Grafted\ Nanoparticles/Molecular\ Dynamics/Implicit\ Solvent/Fixed\ Distance
cd "$destination"/
# Creation of a parent directory for the data from all the simulations associated with the set of parameters specified in params.sh.
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
    find "$destination" -maxdepth 1 -type f -print0 | xargs -0 gcp -t ./
fi

destination=/Users/felipe.pace.evaristo/Library/CloudStorage/Dropbox/Mac/Desktop/Projects/Solvated\ Polymer-Grafted\ Nanoparticles/Molecular\ Dynamics/Implicit\ Solvent/Fixed\ Distance/${radius_nanoparticles}/${length_grafted_chains}/${target_grafting_density}/${epsilon_mon_mon_attractive}
cd "$destination"/

# Exporting of the simulation parameters listed in params.sh into params.txt.
echo "$radius_nanoparticles $length_grafted_chains $target_grafting_density $epsilon_mon_mon_attractive" > params.txt

# Retrieval of distances_monomer_grafted_nanoparticles.txt and importing of the distances between monomer-grafted nanoparticles listed in distances_monomer_grafted_nanoparticles.txt into distances_monomer_grafted_nanoparticles, one distance per entry.
scp felipepacci@cardinal.osc.edu:$origin/distances_monomer_grafted_nanoparticles.txt ./
mapfile -t distances_monomer_grafted_nanoparticles < distances_monomer_grafted_nanoparticles.txt
# mapfile is only available in Bash 4.x or higher.

parameters=${radius_nanoparticles}'_'${length_grafted_chains}'_'${target_grafting_density}'_'${epsilon_mon_mon_attractive}
# Creation of local copies of all .txt files generated for all distances between monomer-grafted nanoparticles sampled.
for index in ${!distances_monomer_grafted_nanoparticles[@]}
do
    mkdir ${distances_monomer_grafted_nanoparticles[index]}
    cd ${distances_monomer_grafted_nanoparticles[index]}/
    for file in *.dat
    do
	scp felipepacci@cardinal.osc.edu:$origin/${distances_monomer_grafted_nanoparticles[index]}/Soft_Pushoff/$file ./
    done
    for file in *.txt
    do
	scp felipepacci@cardinal.osc.edu:$origin/${distances_monomer_grafted_nanoparticles[index]}/Soft_Pushoff/$file ./
	scp felipepacci@cardinal.osc.edu:$origin/${distances_monomer_grafted_nanoparticles[index]}/Equilibration_NVT/$file ./
    done
#    for file in *.lammpstrj.zst
#    do
#	scp felipepacci@cardinal.osc.edu:$origin/${distances_monomer_grafted_nanoparticles[index]}/Soft_Pushoff/$file ./
#	scp felipepacci@cardinal.osc.edu:$origin/${distances_monomer_grafted_nanoparticles[index]}/Equilibration_NVT/$file ./
#    done
    for file in *.dat
    do
	mv "$file" "${file%.dat}_${parameters}_${distances_monomer_grafted_nanoparticles[index]}.dat"
    done
    for file in *.txt
    do
	mv "$file" "${file%.txt}_${parameters}_${distances_monomer_grafted_nanoparticles[index]}.txt"
    done
#    for file in *.lammpstrj.zst
#    do
#	mv "$file" "${file%.lammpstrj.zst}_${parameters}_${distances_monomer_grafted_nanoparticles[index]}.lammpstrj.zst"
#    done
    for file in position*
    do
	sed -i '' "2s/.*/# Timestep Time xu yu zu/" "$file"
    done
    "$destination"/reformatting_average_radii_of_gyration.sh < average_radii_of_gyration_soft_pushoff_${parameters}_${distances_monomer_grafted_nanoparticles[index]}.txt | tee average_radii_of_gyration_soft_pushoff_${parameters}_${distances_monomer_grafted_nanoparticles[index]}_temp.txt > /dev/null
    mv average_radii_of_gyration_soft_pushoff_${parameters}_${distances_monomer_grafted_nanoparticles[index]}_temp.txt average_radii_of_gyration_soft_pushoff_${parameters}_${distances_monomer_grafted_nanoparticles[index]}.txt
    "$destination"/reformatting_average_radii_of_gyration.sh < average_radii_of_gyration_equilibration_nvt_${parameters}_${distances_monomer_grafted_nanoparticles[index]}.txt | tee average_radii_of_gyration_equilibration_nvt_${parameters}_${distances_monomer_grafted_nanoparticles[index]}_temp.txt > /dev/null
    mv average_radii_of_gyration_equilibration_nvt_${parameters}_${distances_monomer_grafted_nanoparticles[index]}_temp.txt average_radii_of_gyration_equilibration_nvt_${parameters}_${distances_monomer_grafted_nanoparticles[index]}.txt
    rm time*.txt
    rm eps*.txt
    cd "$destination"/
done
