## Role of Solvation on Diffusion of Ions in Diblock Copolymers: Understanding the Molecular Weight Effect through Modeling
### To view this publication, click [here](https://pubs.acs.org/doi/10.1021/jacs.9b07227). 

- [denprof_log_mod.py](https://github.com/hall-polymers/published-work/blob/master/2019-seo2019role/denprof_log_mod.py) - python script to create a density profile of block copolymers along z direction

- [mkinput_lamellae_ions.py](https://github.com/hall-polymers/published-work/blob/master/2019-seo2019role/mkinput_lamellae_ions.py) - python script to generate input data file of block copolymers with ions in a lamellar structure

- [msdxy.py](https://github.com/hall-polymers/published-work/blob/master/2019-seo2019role/msdxy.py) - python script to calculate mean square displacement of ions along x and y directions (parallel to the block copolymer interfaces)

- [pair_bornsolv.cpp](https://github.com/hall-polymers/published-work/blob/master/2019-seo2019role/pair_bornsolv.cpp) - LAMMPS C++ script to register 1/r^4 solvation potential to LAMMPS

- [in.diblock](https://github.com/hall-polymers/published-work/blob/master/2019-seo2019role/in.diblock) - LAMMPS script to perform MD simulations of block copolymers with ions (with a 1/r^4 solvation potential)