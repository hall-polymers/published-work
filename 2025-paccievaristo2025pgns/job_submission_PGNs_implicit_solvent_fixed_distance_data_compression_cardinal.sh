#!/bin/bash
#SBATCH --job-name=name
#SBATCH --account=PAS2029
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1 --partition=cpu
#SBATCH --output=./name.out
#SBATCH --error=./name.err
name="PGNs_implicit_solvent_fixed_distance_data_compression"
cd $SLURM_SUBMIT_DIR
./data_compression_PGNs_implicit_solvent_fixed_distance.sh
