#!/bin/bash

. ./params.sh
radius_nanoparticles=`echo $radius_nanoparticles`
# Value of radius_nanoparticles specified in and obtained from params.sh.
length_grafted_chains=`echo $length_grafted_chains`
# Value of length_grafted_chains specified in and obtained from params.sh.
target_grafting_density=`echo $target_grafting_density`
# Value of target_grafting_density specified in and obtained from params.sh.
epsilon_mon_mon_attractive=`echo $epsilon_mon_mon_attractive`
# Value of epsilon_mon_mon_attractive specified in and obtained from params.sh.

home=/fs/scratch/PAS2029/felipepacci/Projects/Solvated_Polymer_Grafted_Nanoparticles/Molecular_Dynamics/Implicit_Solvent/Fixed_Distance/${radius_nanoparticles}/${length_grafted_chains}/${target_grafting_density}/${epsilon_mon_mon_attractive}
cd $home/

echo 'Existence of Non-Empty Output Files'
switch="false"
find . -type f -name "*.out" |
while read -r output_file
do
    if [ ! -e "$output_file" ] || [ ! -s "$output_file" ]
    then
        switch="true"
        echo "The output file $output_file does not exist or exists but is empty."
    fi
done
if [ "$switch" == "false" ]
then
    echo -e ''
fi

echo 'Content of Error Files'
switch="false"
find . -type f -name "*.err" |
while read -r error_file
do
    if [ -s "$error_file" ]
    then
        switch="true"
        echo $error_file
        awk '!seen[$0]++' "$error_file"
    fi
done
if [ "$switch" == "false" ]
then
    echo -e ''
fi

echo 'List of Warnings'
grep --include=\*.{out,err} -Rnw './' -e 'WARNING'

echo -e ''
echo 'List of Errors'
grep --include=\*.{out,err} -Rnw './' -e 'ERROR'

echo -e ''
echo 'List of Dangerous Builds'
grep --include=\*.out -Rnw './' -e 'Dangerous builds ='
