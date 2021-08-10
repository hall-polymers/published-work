#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#!/usr/bin/python
 
# Script:  transport_analysis.py
# Purpose: extract data from txt files created from msd_ions_NEMD.py, calculate cations' and anions' diffusion constant, mobility, degree of uncorrelated ion motion, and make corresponding plots
# Syntax:  transport_analysis.py 
# Author:  Kevin Shen Feb. 2019

import os
import string
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit

# Note: this autovivification method is borrowed from https://stackoverflow.com/questions/635483/what-is-the-best-way-to-implement-nested-dictionaries
class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)() # retain local pointer to value
        return value                     # faster to return than dict lookup

def calculate_alpha(num_frames, estrength, D_dir, output_file='blockavg_results', makeplot=True):
	
	def read_txtfiles(directory="."):
		df = Vividict()
		for root, dirs, files in os.walk(directory):
			for file in files:
				if file.endswith(".txt"):
					rf = file.replace("avged", "_").replace("tot", "_").replace(".","_") # FIXME: could just name the files better while outputing msd data
					if string.split(rf,"_")[-2] == 'x' or string.split(rf,"_")[-2] == D_dir:
						efield_dir = string.split(rf,"_")[6]
						num_blocks = int(string.split(rf,"_")[3])
						nframe = int(string.split(rf,"_")[5])

						if num_frames == nframe:
							df_ind = pd.read_csv(os.path.join(root, file), delim_whitespace=True) # Equivalent to setting sep='\s+'
							df_ind = df_ind.ix[np.log2(df_ind["time"]/100).apply(float.is_integer)| (np.log2(df_ind["time"]/100/3)).apply(float.is_integer)| (np.log2(df_ind["time"]/100/5)).apply(float.is_integer)| (np.log2(df_ind["time"]/100/7)).apply(float.is_integer)]
							df[num_blocks][num_frames][efield_dir] = df_ind
					# print os.path.join(root, file)
		return df

	def find_nearest(array, value):
		array = np.asarray(array)
		nearest_index = (np.abs(array - value)).argmin()
		if array[nearest_index] > value:
			nearest_index -= 1
		return nearest_index

	def plot_setup(plot, title, x_label, y_label, font_size=15):
		# The title
		if title is not None:
			plot.set_title(title, fontweight="bold")
		# The axes
		plot.minorticks_on()
		plot.tick_params(axis='both',which='minor', direction="in",length=5,width=2,labelsize=font_size)
		plot.tick_params(axis='both',which='major', direction="in",length=8,width=2,labelsize=font_size)
		plot.set_xlabel(x_label, fontsize=font_size)
		plot.set_ylabel(y_label, fontsize=font_size)
		# The spines
		plt.setp(plot.spines.values(), linewidth=2)
	
	df = read_txtfiles()
	dif_cat_mean = np.zeros(len(df))
	dif_an_mean = np.zeros(len(df))
	mobil_cat_mean = np.zeros(len(df))
	mobil_an_mean = np.zeros(len(df))
	transf_dif_mean = np.zeros(len(df))
	transf_mobil_mean = np.zeros(len(df))
	alpha_mean = np.zeros(len(df))
	dif_cat_std = np.zeros(len(df))
	dif_an_std = np.zeros(len(df))
	mobil_cat_std = np.zeros(len(df))
	mobil_an_std = np.zeros(len(df))
	transf_dif_std = np.zeros(len(df))
	transf_mobil_std = np.zeros(len(df))
	alpha_std = np.zeros(len(df))
	block_list = np.zeros(len(df))
	for index, num_blocks in enumerate(df):
		time = list(df[num_blocks][num_frames]['x']['time'])
		cut_point = find_nearest(time, num_frames//(num_blocks*10)*100)
		
		nBlock = int(num_blocks)
		dif_cat = np.zeros(nBlock)
		dif_an = np.zeros(nBlock)
		mobil_cat = np.zeros(nBlock)
		mobil_an = np.zeros(nBlock)
		transf_dif = np.zeros(nBlock)
		transf_mobil = np.zeros(nBlock)
		alpha = np.zeros(nBlock)
		for n in range(nBlock):
			cation = np.array(df[num_blocks][num_frames]['x']['blk{}_cat'.format(n)]-df[num_blocks][num_frames][D_dir]['blk{}_cat'.format(n)]/len(D_dir))
			anion = np.array(df[num_blocks][num_frames]['x']['blk{}_an'.format(n)]-df[num_blocks][num_frames][D_dir]['blk{}_an'.format(n)]/len(D_dir))
			fit_func = lambda x, a, b: a*np.square(x)+b
			mobil_slope_cat = curve_fit(fit_func, df[num_blocks][num_frames]['x']['time'][cut_point:], cation[cut_point:])[0][0]
			mobil_slope_an = curve_fit(fit_func, df[num_blocks][num_frames]['x']['time'][cut_point:], anion[cut_point:])[0][0]
			msd_slope_cat = LinearRegression().fit(df[num_blocks][num_frames][D_dir][['time']][cut_point:], df[num_blocks][num_frames][D_dir]['blk{}_cat'.format(n)][cut_point:]).coef_
			msd_slope_an = LinearRegression().fit(df[num_blocks][num_frames][D_dir][['time']][cut_point:], df[num_blocks][num_frames][D_dir]['blk{}_an'.format(n)][cut_point:]).coef_
			dif_cat[n] = msd_slope_cat/(2*len(D_dir))
			dif_an[n] = msd_slope_an/(2*len(D_dir))
			mobil_cat[n] = np.sqrt(mobil_slope_cat)/estrength
			mobil_an[n] = np.sqrt(mobil_slope_an)/estrength
		transf_dif = dif_cat/(dif_cat+dif_an)
		transf_mobil = mobil_cat/(mobil_cat+mobil_an)
		alpha = (mobil_cat+mobil_an)/(dif_cat+dif_an)

		dif_cat_mean[index] = dif_cat.mean()
		dif_an_mean[index] = dif_an.mean()
		mobil_cat_mean[index] = mobil_cat.mean()
		mobil_an_mean[index] = mobil_an.mean()
		transf_dif_mean[index] = transf_dif.mean()
		transf_mobil_mean[index] = transf_mobil.mean()
		alpha_mean[index] = alpha.mean()
		dif_cat_std[index] = dif_cat.std()/np.sqrt(nBlock)
		dif_an_std[index] = dif_an.std()/np.sqrt(nBlock)
		mobil_cat_std[index] = mobil_cat.std()/np.sqrt(nBlock)
		mobil_an_std[index] = mobil_an.std()/np.sqrt(nBlock)
		transf_dif_std[index] = transf_dif.std()/np.sqrt(nBlock)
		transf_mobil_std[index] = transf_mobil.std()/np.sqrt(nBlock)
		alpha_std[index] = alpha.std()/np.sqrt(nBlock)
		block_list[index] = nBlock

		print 'Mean:', num_frames, num_blocks, time[cut_point], time[-1], dif_cat_mean[index], dif_an_mean[index], mobil_cat_mean[index], mobil_an_mean[index], transf_dif_mean[index], transf_mobil_mean[index], alpha_mean[index]
		print 'Stde:', num_frames, num_blocks, time[cut_point], time[-1], dif_cat_std[index], dif_an_std[index], mobil_cat_std[index], mobil_an_std[index], transf_dif_std[index], transf_mobil_std[index], alpha_std[index]
	
	sort = block_list.argsort()
	if makeplot:
		fig, axs = plt.subplots(4, 2, figsize=(16, 16), dpi=120)
		axs[0][0].errorbar(block_list[sort], dif_cat_mean[sort], dif_cat_std[sort], label='1.0C', marker='o', color='C1')
		axs[0][1].errorbar(block_list[sort], mobil_cat_mean[sort], mobil_cat_std[sort], label='1.0C', marker='o', color='C2')
		axs[1][0].errorbar(block_list[sort], dif_an_mean[sort], dif_an_std[sort], label='1.0C', marker='o', color='C3')
		axs[1][1].errorbar(block_list[sort], mobil_an_mean[sort], mobil_an_std[sort], label='1.0C', marker='o', color='C4')
		axs[2][0].errorbar(block_list[sort], transf_dif_mean[sort], transf_dif_std[sort], label='1.0C', marker='o', color='C5')
		axs[2][1].errorbar(block_list[sort], transf_mobil_mean[sort], transf_mobil_std[sort], label='1.0C', marker='o', color='C6')
		axs[3][0].errorbar(block_list[sort], alpha_mean[sort], alpha_std[sort], label='1.0C', marker='o', color='C7')
		plot_setup(axs[0][0], None, 'Number of blocks', r'Cation diffusion constant $(\sigma^2/\tau)$')
		plot_setup(axs[0][1], None, 'Number of blocks', r'Cation mobility $(\sigma^2/\tau)$')
		plot_setup(axs[1][0], None, 'Number of blocks', r'Anion diffusion constant $(\sigma^2/\tau)$')
		plot_setup(axs[1][1], None, 'Number of blocks', r'Anion mobility $(\sigma^2/\tau)$')
		plot_setup(axs[2][0], None, 'Number of blocks', r'Transference # by diffusion')
		plot_setup(axs[2][1], None, 'Number of blocks', r'Transference # by mobility')
		plot_setup(axs[3][0], None, 'Number of blocks', r'alpha')

		for num_blocks in df:
			axs[3][1].plot(df[num_blocks][num_frames]['x']['time'], df[num_blocks][num_frames]['x']['dr_cat'], label='1.0C')
			axs[3][1].plot(df[num_blocks][num_frames]['x']['time'], df[num_blocks][num_frames]['x']['dr_an'], label='1.0C')
		plot_setup(axs[3][1], None, r'time $(\tau)$', r'Ion xcoms $(\sigma)$')

		fig.savefig("{}.png".format(output_file))

def main():
	diffusion_dir = 'yz'
	calculate_alpha(6001, 0.15, diffusion_dir)

if __name__ == "__main__":
	main()
