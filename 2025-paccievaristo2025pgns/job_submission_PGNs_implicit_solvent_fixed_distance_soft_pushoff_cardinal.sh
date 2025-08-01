#!/bin/bash
#SBATCH --job-name=name
#SBATCH --account=PAS2029
#SBATCH --time=1-00:00:00
#SBATCH --output=./name.out
#SBATCH --error=./name.err
name="PGNs_implicit_solvent_fixed_distance_soft_pushoff"

echo "SLURM_JOB_NAME     = " $SLURM_JOB_NAME
echo "SLURM_JOB_ID       = " $SLURM_JOB_ID
echo "SLURM_SUBMIT_DIR   = " $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST = " $SLURM_JOB_NODELIST

LAMMPS=/fs/scratch/PAS2029/felipepacci/lammps-29Aug2024/src

cp ${LAMMPS}/lmp_mpi $SLURM_SUBMIT_DIR/
cd $SLURM_SUBMIT_DIR/

srun lmp_mpi -in ${name}".lmp" | grep -v "WARNING: FENE bond too long"

echo "Job completed."

rm lmp_mpi

for file in *.err
do
    if grep -Fq '[][mvp_generate_implicit_cpu_mapping] WARNING: You appear to be running at full subscription for this job. UCX spawns an additional thread for each process which may result in oversubscribed cores and poor performance. Please consider reserving at least 2 cores per node for the additional threads, enabling SMT, or setting MVP_THREADS_PER_PROCESS=2 to ensure that sufficient resources are available.' "$file"
    then
	sed -i '/\[\]\[mvp_generate_implicit_cpu_mapping\] WARNING: You appear to be running at full subscription for this job. UCX spawns an additional thread for each process which may result in oversubscribed cores and poor performance. Please consider reserving at least 2 cores per node for the additional threads, enabling SMT, or setting MVP_THREADS_PER_PROCESS=2 to ensure that sufficient resources are available./d' "$file"
    fi
done
