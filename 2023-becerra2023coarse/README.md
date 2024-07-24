## Coarse-grained modeling of polymers with end-on and side-on liquid crystal moieties: Effect of architecture
### To view this publication, click [here](https://doi.org/10.1063/5.0152817)

sclcp_lammps_implementation/sclcp_lammps/src/lammps/: Necessary files to make the SCLCP model work with lammps-16Mar2018.
Here, angle_orient.cpp, angle_orient.h, angle_wlc_twist.cpp, angle_wlc_twist.h, and atom.cpp were taken from prior work, 
which you can find [here](https://doi.org/10.1063/1.5092976).
The fix_ext_alignment.cpp and fix_ext_alignment.h files were developed in this work

scripts/initial_configurations/init_conf_X-6-X/

- [lcp.py](https://github.com/hall-polymers/published-work/blob/master/2023-becerra2023coarse/scripts/initial_configurations/init_conf_X-6-X/lcp.py) -- Python script to create SCLCP initial configuration.
  It produces in.lammps (LAMMPS initial configuration) and traj.lammpstrj (LAMMPS trajectory file for visualization)

scripts/lammps_input/

- [1-eosc_step1.in](https://github.com/hall-polymers/published-work/blob/master/2023-becerra2023coarse/scripts/lammps_input/1-eosc_step1.in) -- LAMMPS input script to push off overlapped chains and start equilibration of the EOSC system

- [2-eosc-preproduction.in](https://github.com/hall-polymers/published-work/blob/master/2023-becerra2023coarse/scripts/lammps_input/2-eosc-preproduction.in) -- LAMMPS input script for the final equilibration step of the EOSC system

- [3-eosc_productionT6-T5.in](https://github.com/hall-polymers/published-work/blob/master/2023-becerra2023coarse/scripts/lammps_input/3-eosc_productionT6-T5.in) -- LAMMPS input script for the production run of the EOSC system.
  First the external alignment field is applied to the LC moieties while cooling the system from T=6 to T=5. Second, the system is again relaxed at T=5

- [1-sosc_step1.in](https://github.com/hall-polymers/published-work/blob/master/2023-becerra2023coarse/scripts/lammps_input/1-sosc_step1.in) -- LAMMPS input script to push off overlapped chains and start equilibration of the SOSC system

- [2-sosc-preproduction.in](https://github.com/hall-polymers/published-work/blob/master/2023-becerra2023coarse/scripts/lammps_input/2-sosc-preproduction.in) -- LAMMPS input script for the final equilibration step of the SOSC system

- [3-sosc_productionT6-T5.in](https://github.com/hall-polymers/published-work/blob/master/2023-becerra2023coarse/scripts/lammps_input/3-sosc_productionT6-T5.in) -- LAMMPS input script for the production run of the SOSC system.
  First the external alignment field is applied to the LC moieties while cooling the system from T=6 to T=5. Second, the system is again relaxed at T=5

- [in.variables1.5](https://github.com/hall-polymers/published-work/blob/master/2023-becerra2023coarse/scripts/lammps_input/in.variables1.5) -- Necessary force field parameters for the LAMMPS input script

scripts/properties/local-ordering_P/

- [local_alignment.py](https://github.com/hall-polymers/published-work/blob/master/2023-becerra2023coarse/scripts/properties/local-ordering_P/local_alignment.py) -- Python script to calculate local alignment P.
  Some trajectory files are given in the folder traj_files to test the script

scripts/properties/orientation_S2/1.configurations/1.lcp_oblate_100_6_1.5/: trajectory file (trajlammpstrj.tar.gz) to test the scripts. 
Type ./automate.sh in the terminal

scripts/properties/orientation_S2/2.templates/

- [orientation.py](https://github.com/hall-polymers/published-work/blob/master/2023-becerra2023coarse/scripts/properties/orientation_S2/2.templates/orientation.py) -- Python script to measure global order parameter S2
