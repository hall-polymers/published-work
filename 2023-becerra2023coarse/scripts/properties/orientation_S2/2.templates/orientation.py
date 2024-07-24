#####################################################
# Filename: orientation.py                          #
# Author: Diego Becerra, 2022                       #
#                                                   #
# Characteristics: Measure director vector and      #
# order parameter                                   #
#                                                   #
# Updated to May, 2022                              #
#####################################################
#!/usr/bin/env python

"""
General functions
"""
import os
import sys
import math
import numpy as np
from numpy.linalg import eig
import random
from scipy import stats
from scipy import constants
import seaborn as sns
from pyquaternion import Quaternion
from scipy.spatial.transform import Rotation as Ro

class Atom:
    def __init__(self, id, type, x, y, z, mol, q,
                 ellipsoidflag, density, nx, ny, nz):
        self.id = id
        self.type = type
        self.x = x
        self.y = y
        self.z = z
        self.mol = mol
        self.q = q
        self.ellipsoidflag = ellipsoidflag
        self.density = density
        self.nx = nx
        self.ny = ny
        self.nz = nz

class Ellipsoid:
    def __init__(self, id, shapex, shapey, shapez,
                 quatw, quati, quatj, quatk):
        self.id = id
        self.shapex = shapex
        self.shapey = shapey
        self.shapez = shapez
        self.quatw = quatw
        self.quati = quati
        self.quatj = quatj
        self.quatk = quatk

class Bond:
    def __init__(self, type, a, b):
        self.type = type
        self.a = a
        self.b = b

class Angle:
    def __init__(self, type, a, b, c):
        self.type = type
        self.a = a
        self.b = b
        self.c = c

class LmpConf:
    def __init__(self):
        self.xlo = None
        self.xhi = None
        self.ylo = None
        self.yhi = None
        self.zlo = None
        self.yhi = None
        self.nAtoms = None
        self.nEllipsoids = None
        self.nBonds = None
        self.nAngles = None
        self.nAtomTypes = None
        self.nBondTypes = None
        self.nAngleTypes = None
        self.mass = []
        self.atoms = []
        self.ellipsoids = []
        self.bonds = []
        self.angles = []

    def read_configuration_file(self, confFile):
        try:
            file = open(confFile, "r")
            line = file.readline()
            while line:
                if "atoms" in line:
                    self.nAtoms = int(line.strip().split()[0])
                if "ellipsoids" in line:
                    self.nEllipsoids = int(line.strip().split()[0])
                if "bonds" in line:
                    self.nBonds = int(line.strip().split()[0])
                if "angles" in line:
                    self.nAngles = int(line.strip().split()[0])
                if "dihedrals" in line:
                    self.nDihedrals = int(line.strip().split()[0])
                if "atom types" in line:
                    self.nAtomTypes = int(line.strip().split()[0])
                if "bond types" in line:
                    self.nBondTypes = int(line.strip().split()[0])
                if "angle types" in line:
                    self.nAngleTypes = int(line.strip().split()[0])
                if "dihedral types" in line:
                    self.nDihedralTypes = int(line.strip().split()[0])
                if "xlo xhi" in line:
                    self.xlo = float(line.strip().split()[0])
                    self.xhi = float(line.strip().split()[1])
                if "ylo yhi" in line:
                    self.ylo = float(line.strip().split()[0])
                    self.yhi = float(line.strip().split()[1])
                if "zlo zhi" in line:
                    self.zlo = float(line.strip().split()[0])
                    self.zhi = float(line.strip().split()[1])
                if "Masses" in line :
                    file.readline()
                    for i in range(self.nAtomTypes):
                        l = file.readline().strip().split()
                        self.mass.append(float(l[1]))
                if "Atoms" in line:
                    file.readline()
                    for i in range(self.nAtoms):
                        l = file.readline().strip().split()
                        if len(l) == 9:
                            nx = 0
                            ny = 0
                            nz = 0
                            self.atoms.append(Atom(
                                int(l[0]),
                                int(l[1]),
                                float(l[2]),
                                float(l[3]),
                                float(l[4]),
                                int(l[5]),
                                float(l[6]),
                                int(l[7]),
                                float(l[8]),
                                nx, ny, nz,
                            )
                            )
                        elif len(l) == 12:
                            self.atoms.append(Atom(
                                int(l[0]),
                                int(l[1]),
                                float(l[2]),
                                float(l[3]),
                                float(l[4]),
                                int(l[5]),
                                float(l[6]),
                                int(l[7]),
                                float(l[8]),
                                int(l[9]),
                                int(l[10]),
                                int(l[11]),
                            )
                            )
                        else:
                            print ("Inconsistent number of entries in Atoms line %d for molecular atom_style" % i+1)
                            sys.exit(1)
                if "Ellipsoids" in line:
                    file.readline()
                    for i in range(self.nEllipsoids):
                        l = file.readline().strip().split()
                        self.ellipsoids.append(Ellipsoid(
                            int(l[0]),
                            float(l[1]),
                            float(l[2]),
                            float(l[3]),
                            float(l[4]),
                            float(l[5]),
                            float(l[6]),
                            float(l[7]),
                        )
                        )
                if "Bonds" in line:
                    file.readline()
                    for i in range(self.nBonds):
                        l = file.readline().strip().split()
                        self.bonds.append(Bond(
                            int(l[1]),
                            int(l[2]),
                            int(l[3]),
                        )
                        )
                if "Angles" in line:
                    file.readline()
                    for i in range(self.nAngles):
                        l = file.readline().strip().split()
                        self.angles.append(Angle(
                            int(l[1]),
                            int(l[2]),
                            int(l[3]),
                            int(l[4]),
                        )
                        )
                line = file.readline()
        except IOError:
            print ("Count not open %s" % confFile)
            sys.exit(2)
        file.close()

    def orientation(self):
        Qavg = 0
        N = 0
        for i in range(self.nAtoms):
            a = self.atoms[i]
            e = self.ellipsoids[i]
            if a.type == 3 and a.id == e.id:
                N += 1
                my_quaternionslc = (e.quatw, e.quati, e.quatj, e.quatk)
                r = Ro.from_quat([
                    my_quaternionslc[1],
                    my_quaternionslc[2],
                    my_quaternionslc[3],
                    my_quaternionslc[0],
                ]
                )
                lc_0 = [0, 0, 1]
                lc_vec = Quaternion(matrix = r.as_matrix()).rotate(lc_0)
                uiui = np.array([
                    [lc_vec[0] * lc_vec[0],
                     lc_vec[0] * lc_vec[1],
                     lc_vec[0] * lc_vec[2]],
                    [lc_vec[1] * lc_vec[0],
                     lc_vec[1] * lc_vec[1],
                     lc_vec[1] * lc_vec[2]],
                    [lc_vec[2] * lc_vec[0],
                     lc_vec[2] * lc_vec[1],
                     lc_vec[2] * lc_vec[2]],
                ]
                )
                I = np.array([
                    [1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1],
                ]
                )
                Qi = 1.5 * uiui - 0.5 * I
                Qavg += Qi
        Q = 1 / N * Qavg
        w,v = eig(Q)
        print('eigen-value:', w)
        print('eigen-vector', v)
        max_index = np.argmax(w)

        file_status = open('status','a')
        if os.stat("status").st_size == 0:
            file_status.write("# eigenvalue_max eigenvalue_max^2 eigenvector0 eigenvector1 eigenvector2\n")
        file_status.write(
            str(w[max_index])
            + "\t"
            + str(w[max_index] ** 2)
            + "\t"
            + str(v[0, max_index])
            + "\t"
            + str(v[1, max_index])
            + "\t"
            + str(v[2, max_index])
            + "\n"
        )
        file_status.close()

        d = (v[0, max_index], v[1, max_index], v[2, max_index])
        Savg = 0
        N = 0
        for i in range(self.nAtoms):
            a = self.atoms[i]
            e = self.ellipsoids[i]
            if a.type == 3 and a.id == e.id:
                N += 1
                my_quaternionslc = (e.quatw, e.quati, e.quatj, e.quatk)
                r = Ro.from_quat([
                    my_quaternionslc[1],
                    my_quaternionslc[2],
                    my_quaternionslc[3],
                    my_quaternionslc[0],
                ]
                )
                lc_0 = [0, 0, 1]
                lc_vec = Quaternion(matrix = r.as_matrix()).rotate(lc_0)
                Si = 1.5 * (np.dot(d, lc_vec) ** 2) - 0.5
                Savg += Si
        S = 1 / N * Savg
        print (S)

def get_command_line_args(args):
    if len(args) != 2:
        print ("Usage: %s <input files> <output files>"  % args[0])
        sys.exit(1)
    return args

def main():
    args = get_command_line_args(sys.argv)
    lmpconf = LmpConf()
    lmpconf.read_configuration_file(args[1])
    lmpconf.orientation()

if __name__ == "__main__":
    main()
