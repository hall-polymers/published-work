## Ion Correlation Effects in Salt-Doped Block Copolymers
### To view this publication, click [here](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.127801). 

- [mkinput_lamellae_ions.py](https://github.com/hall-polymers/published-work/blob/master/2018-brown2018ion/mkinput_lamellae_ions.py) -- generate initial configuration of diblock copolymers in lamellae with penetrants (LAMMPS data file)

- [in.diblock](https://github.com/hall-polymers/published-work/blob/master/2018-brown2018ion/in.diblock) -- LAMMPS input script to equilibrate diblock copolymers with monomer penetrants

- [pair_bornsolv.cpp](https://github.com/hall-polymers/published-work/blob/master/2018-brown2018ion/pair_bornsolv.cpp) -- LAMMPS cpp source file for using the customized pair potential. LAMMPS must be compiled with the file to use the bornsolv pair potential in the in.diblock file