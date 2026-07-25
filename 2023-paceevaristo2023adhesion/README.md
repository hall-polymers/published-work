## Ionomer Interfacial Behavior from Molecular Dynamics Simulations: Impact of Ion Content on Interfacial Structure and Mixing
### To view this publication, click [here](https://pubs.acs.org/doi/10.1021/acs.macromol.3c01186).

### Simulation Scripts

- [mkinput_pendant_3.py](https://github.com/hall-polymers/published-work/blob/master/2023-paccievaristo2023adhesion/mkinput_pendant_3.py)  - Python script to create an initial configuration for the Nbb3-50%Na system

- [SPO1Nbb350.imers](https://github.com/hall-polymers/published-work/blob/master/2023-paccievaristo2023adhesion/SPO1Nbb350.imers) - LAMMPS input script to perform a soft pushoff run on the initial configuration of the Nbb3-50%Na system created by mkinput_pendant_3.py

- [equil350.imers](https://github.com/hall-polymers/published-work/blob/master/2023-paccievaristo2023adhesion/equil350.imers) - LAMMPS input script to perform an equilibration run on the configuration of the Nbb3-50%Na system output by SPO1Nbb350.imers

- [run350(1).imers](https://github.com/hall-polymers/published-work/blob/master/2023-paccievaristo2023adhesion/run350(1).imers) - LAMMPS input script to perform a data collection run on the freestanding ionomer film of the Nbb3-50%Na system output by equil350.imers

- [run350(2).imers](https://github.com/hall-polymers/published-work/blob/master/2023-paccievaristo2023adhesion/run350(2).imers) - LAMMPS input script to perform an MD simulation run on the adhering ionomer film surfaces of the Nbb3-50%Na system created from the output of run350(1).imers

### Analysis Scripts

- [Density_SingleFilm.py](https://github.com/hall-polymers/published-work/blob/master/2023-paccievaristo2023adhesion/Density_SingleFilm.py) - Calculates the density profiles by bead type of a single freestanding film as a function of distance from the center of the film (Figure 2)

- [Orientation_SingleFilm.py](https://github.com/hall-polymers/published-work/blob/master/2023-paccievaristo2023adhesion/Orientation_SingleFilm.py) - Calculates the bond orientation parameter of backbone bond vectors in a single freestanding film as a function of distance from the center of the film (Figure 3)

- [Concentration_TwoFilms.py](https://github.com/hall-polymers/published-work/blob/master/2023-paccievaristo2023adhesion/Concentration_TwoFilms.py) - Calculates the concentration profile of all beads as a function of distance from the healing interface (Figure 6)

- [Overlap_Parameter_and_Interdiffusion_Distance_3_50.xlsx](https://github.com/hall-polymers/published-work/blob/master/2023-paccievaristo2023adhesion/Overlap_Parameter_and_Interdiffusion_Distance_3_50.xlsx) - Calculates overlap parameter and interdiffusion distance from the counterion and noncounterion density profiles of the Nbb3-50%Na system (Figure 7 and Figure 8)
