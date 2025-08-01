## Coarse-Grained Simulations of Pairs of Polymer-Grafted Nanoparticles in Implicit Solvent
### To view this publication, click [here]().

### Simulation Scripts

- [params.sh](https://github.com/hall-polymers/published-work/blob/master/2025-paccievaristo2025pgns/params.sh) - Data file containing the parameters that define the system (i.e., nanoparticle radius, chain length, grafting density, and solvent quality)

- [make_data_PGNs_implicit_solvent_fixed_distance.py](https://github.com/hall-polymers/published-work/blob/master/2025-paccievaristo2025pgns/make_data_PGNs_implicit_solvent_fixed_distance.py) - Python script to create an initial configuration for the system defined by the parameters in params.sh

- [PGNs_implicit_solvent_fixed_distance_soft_pushoff.lmp](https://github.com/hall-polymers/published-work/blob/master/2025-paccievaristo2025pgns/PGNs_implicit_solvent_fixed_distance_soft_pushoff.lmp) - LAMMPS input script to perform a soft pushoff run on the initial configuration of the system created by make_data_PGNs_implicit_solvent_fixed_distance.py

- [PGNs_implicit_solvent_fixed_distance_equilibration_nvt.lmp](https://github.com/hall-polymers/published-work/blob/master/2025-paccievaristo2025pgns/PGNs_implicit_solvent_fixed_distance_equilibration_nvt.lmp) - LAMMPS input script to perform an equilibration run on the configuration of the system output by PGNs_implicit_solvent_fixed_distance_soft_pushoff.lmp

### Analysis Scripts

- [PGNs_implicit_solvent.m](https://github.com/hall-polymers/published-work/blob/master/2025-paccievaristo2025pgns/PGNs_implicit_solvent.m) - Plots the properties calculated by PGNs_implicit_solvent_fixed_distance_equilibration_nvt.lmp and calculates and plots the potential of mean force (PMF) of the system (Figure 6a)
