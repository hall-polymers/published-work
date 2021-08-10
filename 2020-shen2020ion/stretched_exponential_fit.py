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
 
# Script:  stretched_exponential_fit.py
# Purpose: fit F(k,t) or ACF_cluster results
# Syntax:  stretched_exponential_fit.py 
# Author:  Kevin Shen Apr. 2019

import os, string, json, math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution

# NOTE: this script works only for files with file names specified in the read_txtfiles function
# NOTE: this script works only for the following dir/file structure:
'''
root/
|
|-- stretched_exponential_fit.py
|
|-- sys1/ 
|   |-- fs_2sigma.txt
|
|-- sys2/ 
    |-- fs_2sigma.txt

'''
# Note: this autovivification method is from https://stackoverflow.com/questions/635483/what-is-the-best-way-to-implement-nested-dictionaries
class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

def read_txtfiles(directory="."):
	df = Vividict()
	for root, dirs, files in os.walk(directory):
		for file in files:
			if file.endswith("fs_2sigma.txt"):
				sys = root.split('/')[1]
				df[sys] = pd.read_csv(open(os.path.join(root, file)), delim_whitespace=True) # Equivalent to setting sep='\s+', delim_whitespace=True
				print os.path.join(root, file)
	return df

def Exponential_Fit(x, T, B):
    return np.exp(-(x/T)**B)

def Exponential_Fit_difevo(parameters, *data):
    T,B = parameters
    x,y,z = data 
    
    result = 0
    for i in range(len(x)):
    	result += ((np.exp(-(x[i]/T)**B) - y[i])*np.reciprocal(z[i]))**2
    
    return result**0.5


def extract_exponential_fit():
	df = read_txtfiles()
	fig, axs = plt.subplots(1, 1, figsize=(20, 10), dpi=200)

	print 'Exponential Fit:'
	for sys in sorted(df):
		x = np.array(df[sys]['timesteps'])
		y = np.array(df[sys]['avg'])
		z = np.array(df[sys]['std'])
		df[sys] = df[sys][df[sys]['timesteps']>0]
		df[sys] = df[sys][df[sys]['std']>0] # to get rid of /0 issue
		x_fit = np.array(df[sys]['timesteps'])
		y_fit = np.array(df[sys]['avg'])
		z_fit = np.array(df[sys]['std']) 
		
		args = (x_fit, y_fit, z_fit)

		# Gaussian fit
		try:
			result = differential_evolution(Exponential_Fit_difevo, bounds=[(0, 1000),(0, 1)], args=args) #, tol=0, maxiter=10000
			relaxation_time = result.x[0]/result.x[1]*math.gamma(1/result.x[1])
			print sys, result.x, relaxation_time
			axs.scatter(x, y, label=r'{0} ($\tau={1:.2f}, \beta={2:.2f}, \tau_R={3:.2f}$)'.format(sys, result.x[0], result.x[1], relaxation_time))
			axs.plot(x_fit, Exponential_Fit(x_fit, *result.x))
			for i in Exponential_Fit(x_fit, *result.x):
				print i

		except:
			print 'Exponential does not fit {}'.format(sys)
	
	axs.set_title('Exponential Fit', fontweight='bold')
	axs.legend()
	axs.set_yscale('log')
	axs.set_xscale('log')
	axs.set_xlabel(r'Time', fontweight='bold')
	axs.set_ylabel('Average', fontweight='bold')
	axs.set_xlim([0.01,np.max(x)])
	axs.set_ylim([0.01,1])
	fig.savefig('fs_2sigma_std.png')


def main():
	extract_exponential_fit()


if __name__ == "__main__":
	main()