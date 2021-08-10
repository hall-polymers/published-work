#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


#!/usr/bin/env python

import sys
sys.path.append('/Users/brownj/Research/Scripts')
import numpy as np
import pppmd



r, ir, timestep, boxbds, id2type, id2mol, mol2ids = pppmd.read_lammpstrj('dynamicnpt.lammpstrj',skip_between=0)

e2e_mols = []
for mol in range(len(mol2ids)):
	if len(mol2ids[mol]) > 1:
		e2e_mols.append(mol2ids[mol])

e2e_autocorr = pppmd.end2end_autocorr(r, ir, boxbds, e2e_mols)

print "t Reet_dot_Ree0"
for t in range(len(timestep)):
	print timestep[t]-timestep[0], e2e_autocorr[t]
