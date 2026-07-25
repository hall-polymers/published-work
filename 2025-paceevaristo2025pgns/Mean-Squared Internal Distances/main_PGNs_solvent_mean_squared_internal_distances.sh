#!/bin/bash

source ~/venv/bin/activate
# Activation of virtual environment containing pandas required for script to run.

home=.
cd $home/
file=params.sh
if test -f "$file"
then
    echo "Verify that the parameter values in params.sh indeed correspond to the simulation of interest."
    echo ""
    . ./params.sh
    radius_nanoparticles=`echo $radius_nanoparticles`
    # Value of radius_nanoparticles specified in and obtained from params.sh.
    length_grafted_chains=`echo $length_grafted_chains`
    # Value of length_grafted_chains specified in and obtained from params.sh.
    target_grafting_density=`echo $target_grafting_density`
    # Value of target_grafting_density specified in and obtained from params.sh.
    epsilon_mon_mon_attractive=`echo $epsilon_mon_mon_attractive`
    # Value of epsilon_mon_mon_attractive specified in and obtained from params.sh.
    read -p "Enter the distance between monomer-grafted nanoparticles in σ that corresponds to the simulation of interest: " distance_monomer_grafted_nanoparticles
    # Value of distance_monomer_grafted_nanoparticles provided by user.
else
    echo "Script requires existence and presence of params.sh in its directory."
    exit 1
fi

home=`pwd`
simulation="$PWD/.."/${radius_nanoparticles}/${length_grafted_chains}/${target_grafting_density}/${epsilon_mon_mon_attractive}/${distance_monomer_grafted_nanoparticles}
input="$home"/Input
output="$home"/Output
cd "$simulation"/
flag=false
counter=0
for file in *.dat
do
    let counter+=1
    template_file=$file
    cp "$simulation"/$template_file "$input"/
done
if (($counter!=1))
then
    echo "There is none or more than one template file in the simulation directory and there should only be one."
    flag=true
fi
counter=0
for file in *.lammpstrj
do
    if echo "$file" | grep -q "equilibration"
    then
	if echo "$file" | grep -q "truncated"
	then
	    :
	else
	    let counter+=1
	    read -p "Enter the number of final frames in the equilibration LAMMPS trajectory file to export to the truncated LAMMPS trajectory file: " number_frames
	    # Value of number_frames provided by user.
	    number_atoms=`sed '4q;d' $file`
	    index_last_line=`wc -l < $file`
	    # More than one truncated trajectory file is generated because a single one is oftentimes so large it causes segmentation faults when used as an input file.
	    index_first_line=$(($index_last_line+1-($number_frames/2)*($number_atoms+9)))
	    trajectory_final_file="${file%.lammpstrj}_truncated_final.lammpstrj"
	    sed -n "${index_first_line},${index_last_line}p" $file > "$input"/$trajectory_final_file
	    index_first_line=$(($index_last_line+1-$number_frames*($number_atoms+9)))
	    index_last_line=$(($index_last_line-($number_frames/2)*($number_atoms+9)))
	    trajectory_initial_file="${file%.lammpstrj}_truncated_initial.lammpstrj"
	    sed -n "${index_first_line},${index_last_line}p" $file > "$input"/$trajectory_initial_file
	fi
    fi
done
if (($counter!=1))
then
    echo "There is none or more than one equilibration LAMMPS trajectory file in the simulation directory and there should only be one."
    flag=true
fi
if $flag
then
    rm "$input"/*.dat
    rm "$input"/*.lammpstrj
    exit 1
fi

mkdir "$output"/Post-Processed
post_processed="$output"/Post-Processed
echo ""
mv "$input"/{frame_file_generation.py,$template_file,$trajectory_initial_file,$trajectory_final_file} "$output"/
wait
python3 "$output"/frame_file_generation.py "$output"/$template_file "$output"/$trajectory_initial_file
wait
echo ""
wait
python3 "$output"/frame_file_generation.py "$output"/$template_file "$output"/$trajectory_final_file
wait
mv "$output"/{frame_file_generation.py,$template_file,$trajectory_initial_file,$trajectory_final_file} "$input"/
cd "$output"/
echo ""
echo "Generating frame files with unscaled and unwrapped coordinates from frame files with scaled and wrapped coordinates."
ls *.dat | sed 's/^\([^0-9]*\)\([0-9]*\)/\1 \2/' | sort -k2,2n | tr -d ' ' |
while read file
do
    timestep=`awk '/Timestep/{print $NF}' $file`
    echo "Timestep: $timestep"
    scaled_wrapped_coordinates_frame_file=$file
    unscaled_unwrapped_coordinates_frame_file="${file%.dat}_reworked.dat"
    wait
    python3 "$input"/frame_file_coordinates_unscaling_unwrapping.py "$output"/$scaled_wrapped_coordinates_frame_file "$post_processed"/$unscaled_unwrapped_coordinates_frame_file
    wait
done
cd "$post_processed"/
echo ""
echo "Generating data files with mean-squared internal distances and average mean-squared internal distances of grafted chains in frame files with unscaled and unwrapped coordinates."
ls *reworked.dat | sed 's/^\([^0-9]*\)\([0-9]*\)/\1 \2/' | sort -k2,2n | tr -d ' ' |
while read file
do
    timestep=`awk '/Timestep/{print $NF}' $file`
    echo "Timestep: $timestep"
    unscaled_unwrapped_coordinates_frame_file=$file
    mean_squared_internal_distances_radial_data_file=mean_squared_internal_distances_radial_"$timestep".txt
    mean_squared_internal_distances_antiradial_data_file=mean_squared_internal_distances_antiradial_"$timestep".txt
    average_mean_squared_internal_distances_radial_data_file=average_mean_squared_internal_distances_radial_"$timestep".txt
    average_mean_squared_internal_distances_antiradial_data_file=average_mean_squared_internal_distances_antiradial_"$timestep".txt
    wait
    python3 "$input"/mean_squared_internal_distances_radial_calculation_PGNs_solvent.py "$post_processed"/$unscaled_unwrapped_coordinates_frame_file "$post_processed"/$mean_squared_internal_distances_radial_data_file "$post_processed"/$average_mean_squared_internal_distances_radial_data_file
    wait
    python3 "$input"/mean_squared_internal_distances_antiradial_calculation_PGNs_solvent.py "$post_processed"/$unscaled_unwrapped_coordinates_frame_file "$post_processed"/$mean_squared_internal_distances_antiradial_data_file "$post_processed"/$average_mean_squared_internal_distances_antiradial_data_file
    wait
done
echo ""
echo "Generating data files with combined average mean-squared internal distances and ensemble averaged mean-squared internal distances."
for direction_mean_squared_internal_distances in "radial" "antiradial"
do
    wait
    python3 "$input"/mean_squared_internal_distances_concatenation_averaging_PGNs_solvent.py "$direction_mean_squared_internal_distances" "$post_processed"/ combined_average_mean_squared_internal_distances_"$direction_mean_squared_internal_distances".txt ensemble_averaged_mean_squared_internal_distances_"$direction_mean_squared_internal_distances".txt
    wait
done
find "$post_processed"/ -name '*.*' -exec mv {} "$simulation"/ \;
rm -r "$post_processed"/
find "$output"/ -name '*.*' -exec mv {} "$simulation"/ \;
rm "$input"/*.dat
mv "$input"/*.lammpstrj "$simulation"/
