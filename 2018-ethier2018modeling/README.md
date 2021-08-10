## Modeling individual and pairs of adsorbed polymer-grafted nanoparticles: structure and entanglements
### To view this publication, click [here](https://pubs.rsc.org/en/content/articlelanding/2018/sm/c7sm02116j#!divAbstract). 

- [input_2PGN6_160.py](https://github.com/hall-polymers/published-work/blob/master/2018-ethier2018modeling/input_2PGN6_160.py) -- generate a random initial configuration for two polymer-grafted nanoparticles (at a specified graft density and chain length) adsorbed on a surface

- [in_nvt.nano](https://github.com/hall-polymers/published-work/blob/master/2018-ethier2018modeling/in_nvt.nano) -- LAMMPS input script to push off overlapped chains and equilibrate on a 9/3 LJ surface

- [radavght3.py](https://github.com/hall-polymers/published-work/blob/master/2018-ethier2018modeling/radavght3.py) -- python script to record the height at which a density of 0.25 is reached in the simulation box (radially averaged height profile)

- [plotcontour.py](https://github.com/hall-polymers/published-work/blob/master/2018-ethier2018modeling/plotcontour.py) -- python script to plot the results of radavght3.py from an excel file (3 columns: x coordinate, y coordinate, and height)

- [readentangle_numtwokinks.py](https://github.com/hall-polymers/published-work/blob/master/2018-ethier2018modeling/readentangle_numtwokinks.py) -- python script to calculate interparticle double kinks (interparticle entanglements) (Z1 output files REQUIRED)



- [denprof_bygroup.py](https://github.com/hall-polymers/published-work/blob/master/2015-seo2015effect/denprof_bygroup.py) – python script to obtain the density profiles (as a function of z) of lamellar structured tapered polymers as grouped by number fraction of one monomer

- [in.diblock](https://github.com/hall-polymers/published-work/blob/master/2015-seo2015effect/in.diblock) – LAMMPS script to equilibrate the system and obtain the outputs

- [mkinput_lamellae_taper.py](https://github.com/hall-polymers/published-work/blob/master/2015-seo2015effect/mkinput_lamellae_taper.py) – python script to create an input data of tapered polymers initially in lamellar structure

- [msd.py](https://github.com/hall-polymers/published-work/blob/master/2015-seo2015effect/msd.py) – python script to calculate msd

- [pppmd.py](https://github.com/hall-polymers/published-work/blob/master/2015-seo2015effect/pppmd.py) – python script to calculate the end-to-end vector autocorrelation function