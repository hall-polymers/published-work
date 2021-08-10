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
 
# Script:  extract_msds.py
# Purpose: extract data from txt files created from msd_ions_NEMD.py, calculate the slopes for ion diffusion constant and mobility, and make corresponding plots
# Syntax:  extract_msds.py 
# Author:  Kevin Shen Feb. 2019

import os, string, json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit

# NOTE: this script works only for files with file names: ion_whole_msdavged%d_tot%d_x and ion_whole_msdavged%d_tot%d_yz, FIXME: make a more general one
# NOTE: this script works only for the following dir/file structure:
'''
root/
|__ sys1/ (example structure for NEMD)
	|	condition1/
	|	|-- extract_msds.py
	|	|-- ion_whole_msdavged12_tot6001_x
	|	|__ ion_whole_msdavged12_tot6001_yz
	|__ condition1/
		|-- extract_msds.py
		|-- ion_whole_msdavged6_tot6001_x
		|__ ion_whole_msdavged6_tot6001_yz
'''
# Note: this autovivification method is from https://stackoverflow.com/questions/635483/what-is-the-best-way-to-implement-nested-dictionaries
class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

def read_txtfiles(nframe, D_dir, directory="."):
	df = Vividict()
	for root, dirs, files in os.walk(directory):
		for file in files:
			if file.endswith(".txt"):
				rf = file.replace("avged", "_").replace("tot", "_").replace(".","_") # FIXME: could just name the files better while outputing msd data
				if string.split(rf,"_")[-2] == 'x' or string.split(rf,"_")[-2] == D_dir:
					efield_dir = string.split(rf,"_")[6]
					num_blocks = int(string.split(rf,"_")[3])
					num_frames = int(string.split(rf,"_")[5])
					if num_frames == nframe:
						df[num_blocks][num_frames][efield_dir] = pd.read_csv(os.path.join(root, file), delim_whitespace=True) # Equivalent to setting sep='\s+'
	return df

def plot_setup(plot, title, x_label, y_label, font_size=20):
	# The title
	if title is not None:
		plot.set_title(title, fontweight="bold")

	# The axes
	plot.minorticks_on()
	plot.tick_params(axis='both',which='minor', direction="in",length=5,width=2,labelsize=font_size)
	plot.tick_params(axis='both',which='major', direction="in",length=8,width=2,labelsize=font_size)
	plot.set_xlabel(x_label, fontsize=font_size)
	plot.set_ylabel(y_label, fontsize=font_size)
	plot.set_xscale("log")

	# The spines
	plt.setp(plot.spines.values(), linewidth=2)

def extract_delta_msd(df, num_frames, D_dir):
	delta_msd = Vividict()
	delta_msd._name = "delta_msd"
	for num_blocks in df:
		for num_frames in df[num_blocks]:
			delta_msd_df = pd.DataFrame({"time": df[num_blocks][num_frames]["x"]["time"]})
			delta_msd_df['delta_msd_cation'] = np.array(df[num_blocks][num_frames]['x']['avg_cat']-df[num_blocks][num_frames][D_dir]['avg_cat']/len(D_dir))
			delta_msd_df['delta_msd_anion'] = np.array(df[num_blocks][num_frames]['x']['avg_an']-df[num_blocks][num_frames][D_dir]['avg_an']/len(D_dir))
			delta_msd[num_blocks][num_frames] = delta_msd_df
	return delta_msd

def read_delta_msd(delta_msd_df):
	fig, axs = plt.subplots(3, 2, figsize=(30, 20), dpi=200)
	for num_blocks in delta_msd_df:
		for num_frames in delta_msd_df[num_blocks]:
			df = delta_msd_df[num_blocks][num_frames]
			for col in df.columns[1:]: # FIXME: this can be written in lambda instead of an actual for loop?
				df = df.ix[df[col] > 0]
			df = df.ix[np.log2(df["time"]/100).apply(float.is_integer)| (np.log2(df["time"]/100/3)).apply(float.is_integer)| (np.log2(df["time"]/100/5)).apply(float.is_integer)| (np.log2(df["time"]/100/7)).apply(float.is_integer)]

			df_log = np.log10(df)

			X = df["time"]
			X_log = df_log[["time"]]
			ndata_least = 5 # must be an odd number
			time = np.array(df["time"][ndata_least/2:-(ndata_least/2)])
			len_data = df.shape[0]-ndata_least+1
			slope_inst = np.zeros((df.shape[1]-1, len_data))
			slope_avg = np.zeros((df.shape[1]-1, len_data))
			slope_linear = np.zeros((df.shape[1]-1, len_data))

			for col in df.columns[1:]:
				y = df[col]
				y_log = df_log[col]
				col_loc = df.columns[1:].get_loc(col)

				for row_loc in range(len_data):
					# calculate instantaneous log-log slope
					reg = LinearRegression().fit(X_log[row_loc:row_loc+ndata_least], y_log[row_loc:row_loc+ndata_least])
					slope_inst[col_loc][row_loc] = reg.coef_
					
					# calculate average log-log slope
					reg2 = LinearRegression().fit(X_log[row_loc:], y_log[row_loc:])
					slope_avg[col_loc][row_loc] = reg2.coef_

					# calculate average linear-linear slope
					fit_func = lambda x, a, b: a*np.square(x)+b
					popt, pcov = curve_fit(fit_func, X[row_loc:], y[row_loc:])
					slope_linear[col_loc][row_loc] = popt[0]

				axs[2, 0].plot(time, slope_linear[col_loc], marker=".")
				axs[2, 1].plot(X, y, marker=".")
				axs[0, 1].plot(time, slope_inst[col_loc]-2, marker=".")
				axs[1, 1].plot(time, slope_avg[col_loc]-2, marker=".")
			axs[0, 0].plot(time, np.mean(slope_inst, axis=0), label="{}frames_{}blocks".format(num_frames, num_blocks), marker="o", linewidth=2)
			axs[1, 0].plot(time, np.mean(slope_avg, axis=0), marker="o", linewidth=2)
			# axs[2, 0].plot(time, np.mean(slope_linear, axis=0), marker="o")

			# dev_fun = lambda x, mean: np.mean(x, axis=0)-mean
			# axs[0, 1].plot(time, dev_fun(slope_inst, 2), marker=".")
			# axs[1, 1].plot(time, dev_fun(slope_avg, 2), marker=".")
			axs[0, 1].plot(time, np.zeros(len_data), linewidth=2, color="black")
			axs[1, 1].plot(time, np.zeros(len_data), linewidth=2, color="black")
		
	plot_setup(axs[0, 0], None, r'time $(\tau)$', r'instant log-log slope $(\sigma^2/\tau)$')
	plot_setup(axs[1, 0], None, r'time $(\tau)$', r'average log-log slope $(\sigma^2/\tau)$')
	plot_setup(axs[2, 0], None, r'time $(\tau)$', r'average slope $(\sigma^2/\tau)$')
	plot_setup(axs[0, 1], None, r'time $(\tau)$', r'level of deviation from 2 $(\sigma^2/\tau)$')
	plot_setup(axs[1, 1], None, r'time $(\tau)$', r'level of deviation from 2 $(\sigma^2/\tau)$')
	plot_setup(axs[2, 1], None, r'time $(\tau)$', r'Mean squared displacement $(\sigma^2)$')
	axs[2, 1].set_yscale("log")
	axs[2, 1].grid(which='major', linestyle='-')

	fig.legend(loc="center right", fontsize=18)

	fig.savefig("drift_vel_overall.png")

def read_msd(msd_df):
	fig2, axs = plt.subplots(3, 2, figsize=(30, 20), dpi=200)
	for num_blocks in msd_df:
		for num_frames in msd_df[num_blocks]:
			df = msd_df[num_blocks][num_frames]
			for col in df.columns[1:]: # FIXME: this can be written in lambda instead of an actual for loop?
				df = df.ix[df[col] > 0]
			df = df.ix[np.log2(df["time"]/100).apply(float.is_integer)| (np.log2(df["time"]/100/3)).apply(float.is_integer)| (np.log2(df["time"]/100/5)).apply(float.is_integer)| (np.log2(df["time"]/100/7)).apply(float.is_integer)]
			df_log = np.log10(df)

			X = df[["time"]]
			X_log = df_log[["time"]]
			ndata_least = 5 # must be an odd number
			time = np.array(df["time"][ndata_least/2:-(ndata_least/2)])

			len_data = df.shape[0]-ndata_least+1
			slope_inst = np.zeros((df.shape[1]-1, len_data))
			slope_avg = np.zeros((df.shape[1]-1, len_data))
			slope_linear = np.zeros((df.shape[1]-1, len_data))

			for col in df.columns[1:]:
				y = df[col]
				y_log = df_log[col]
				col_loc = df.columns[1:].get_loc(col)

				for row_loc in range(len_data):
					# calculate instantaneous log-log slope
					reg = LinearRegression().fit(X_log[row_loc:row_loc+ndata_least], y_log[row_loc:row_loc+ndata_least])
					slope_inst[col_loc][row_loc] = reg.coef_
					
					# calculate average log-log slope
					reg2 = LinearRegression().fit(X_log[row_loc:], y_log[row_loc:])
					slope_avg[col_loc][row_loc] = reg2.coef_

					# calculate average linear-linear slope (at least "ndata_least" points at the end)
					reg3 = LinearRegression().fit(X[row_loc:], y[row_loc:])
					slope_linear[col_loc][row_loc] = reg3.coef_

				axs[2, 0].plot(time, slope_linear[col_loc], marker=".")
				axs[2, 1].plot(X, y, marker=".")
				axs[0, 1].plot(time, slope_inst[col_loc]-1, marker=".")
				axs[1, 1].plot(time, slope_avg[col_loc]-1, marker=".")
			axs[0, 0].plot(time, np.mean(slope_inst, axis=0), label="{}frames_{}blocks".format(num_frames, num_blocks), marker="o", linewidth=2)
			axs[1, 0].plot(time, np.mean(slope_avg, axis=0), marker="o", linewidth=2)
			# axs[2, 0].plot(time, np.mean(slope_linear, axis=0), marker="o")

			# dev_fun = lambda x, mean: np.mean(x, axis=0)-mean
			# axs[0, 1].plot(time, dev_fun(slope_inst, 1), marker=".")
			# axs[1, 1].plot(time, dev_fun(slope_avg, 1), marker=".")
			axs[0, 1].plot(time, np.zeros(len_data), linewidth=2, color="black")
			axs[1, 1].plot(time, np.zeros(len_data), linewidth=2, color="black")
		
	plot_setup(axs[0, 0], None, r'time $(\tau)$', r'instant log-log slope $(\sigma^2/\tau)$')
	plot_setup(axs[1, 0], None, r'time $(\tau)$', r'average log-log slope $(\sigma^2/\tau)$')
	plot_setup(axs[2, 0], None, r'time $(\tau)$', r'average slope $(\sigma^2/\tau)$')
	plot_setup(axs[0, 1], None, r'time $(\tau)$', r'level of deviation from 1 $(\sigma^2/\tau)$')
	plot_setup(axs[1, 1], None, r'time $(\tau)$', r'level of deviation from 1 $(\sigma^2/\tau)$')
	plot_setup(axs[2, 1], None, r'time $(\tau)$', r'Mean squared displacement $(\sigma^2)$')
	axs[2, 1].set_yscale("log")
	axs[2, 1].grid(which='major', linestyle='-')

	fig2.legend(loc="center right", fontsize=18)

	fig2.savefig("diffusion_overall.png")

def output_result(vividict, output_dir):
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	for num_blocks in vividict:
		for num_frames in vividict[num_blocks]:
			file = "{}_{}frames_{}blocks.csv".format(vividict._name, num_frames, num_blocks)
			vividict[num_blocks][num_frames].to_csv(os.path.join(output_dir,file), index=False)

def extract_msd(df, num_frames, D_dir):
	msd = Vividict()
	msd._name = "msd"
	for num_blocks in df:
		for num_frames in df[num_blocks]:
			msd_df = pd.DataFrame({"time": df[num_blocks][num_frames]["x"]["time"]})
			msd_df['avg_cat'] = np.array(df[num_blocks][num_frames][D_dir]['avg_cat']/len(D_dir))
			msd_df['avg_an'] = np.array(df[num_blocks][num_frames][D_dir]['avg_an']/len(D_dir))
			msd[num_blocks][num_frames] = msd_df
	return msd


def main():
	num_frames = 6001
	diffusion_dir = 'yz'

	df = read_txtfiles(num_frames, diffusion_dir)
	delta_msd = extract_delta_msd(df, num_frames, diffusion_dir)
	read_delta_msd(delta_msd)
	# output_result(delta_msd, "analysis")

	msd = extract_msd(df, num_frames, diffusion_dir)
	read_msd(msd)
	# output_result(msd, "analysis")


if __name__ == "__main__":
	main()