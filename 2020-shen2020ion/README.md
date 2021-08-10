## Ion Conductivity and Correlations in Model Salt-Doped Polymers: Effects of Interaction Strength and Concentration
### To view this publication, click [here](https://pubs.acs.org/doi/10.1021/acs.macromol.0c00216). 

- [mkinput_homo_ions.py](https://github.com/hall-polymers/published-work/blob/master/2020-shen2020ion/mkinput_homo_ions.py) - python script to generate input data file of a homopolymer melt with ions

- [mkinput_lamellae_ions.py](https://github.com/hall-polymers/published-work/blob/master/2020-shen2020ion/mkinput_lamellae_ions.py) - python script to generate input data file of block copolymers with ions in a lamellar structure

- [in.nvt](https://github.com/hall-polymers/published-work/blob/master/2020-shen2020ion/in.nvt) - LAMMPS script to perform MD simulations (with a 1/r^4 solvation potential) with an applied electric field

- [msd_ions_EMD.py](https://github.com/hall-polymers/published-work/blob/master/2020-shen2020ion/msd_ions_emd.py) - python script to calculate collective and mean square displacement of ions from EMD simulations

- [msd_ions_NEMD.py](https://github.com/hall-polymers/published-work/blob/master/2020-shen2020ion/msd_ions_nemd.py) - python script to calculate mean square displacement of ions in directions parallel and perpendicular to the electric field from NEMD simulations

- [extract_msds.py](https://github.com/hall-polymers/published-work/blob/master/2020-shen2020ion/extract_msds.py) - python script to extract data from txt files created from msd_ions_NEMD.py, calculate the slopes for ion diffusion constant and mobility, and make corresponding plots

- [transport_analysis.py](https://github.com/hall-polymers/published-work/blob/master/2020-shen2020ion/transport_analysis.py) - extract data from txt files created from msd_ions_NEMD.py, calculate cations' and anions' diffusion constant, mobility, degree of uncorrelated ion motion, and make corresponding plots

- [Fs_direct_avg.py](https://github.com/hall-polymers/published-work/blob/master/2020-shen2020ion/fs_direct_avg.py) - python script to calculate self-intermediate scattering function

- [clustrautocorr.py](https://github.com/hall-polymers/published-work/blob/master/2020-shen2020ion/clustrautocorr.py) - python script to calculate ion cluster autocorrelation function

- [stretched_exponential_fit.py](https://github.com/hall-polymers/published-work/blob/master/2020-shen2020ion/stretched_exponential_fit.py) - python script to fit polymer self-intermediate scattering function and ion cluster autocorrelation function with a stretched exponential function to get polymer relaxation rate and cluster relaxation rate, respectively

- [ip_iccorr.py](https://github.com/hall-polymers/published-work/blob/master/2020-shen2020ion/ip_iccorr.py) - python script to calculate ion pair/cage relaxation rate (data and description are in the Supporting Information)