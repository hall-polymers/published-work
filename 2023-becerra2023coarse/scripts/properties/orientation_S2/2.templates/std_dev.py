#####################################################
# Filename: orientation.py                          #
# Author: Diego Becerra, 2022                       #
#                                                   #
# Characteristics: Calculates the std dev of Q      #
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

class Status:
    def __init__(self, evalue_max, evalue_max2, evector0, evector1, evector2):
        self.evalue_max = evalue_max
        self.evalue_max2 = evalue_max2
        self.evector0 = evector0
        self.evector1 = evector1
        self.evector2 = evector2

class LmpConf:
    def __init__(self):
        self.nStatus = None
        self.status = []

    def read_configuration_file(self, confFile):
        try:
            file = open(confFile, "r")
            line = file.readline()
            while line:
                if "# eigenvalue_max" in line:
                    #file.readline()
                    for i in range(0, 10):
                        l = file.readline().strip().split()
                        if len(l) == 5:
                            self.status.append(Status(
                                float(l[0]),
                                float(l[1]),
                                float(l[2]),
                                float(l[3]),
                                float(l[4]),
                            )
                            )
                        else:
                            print ("Inconsistent number of entries in Status line %d" % i + 1)
                            sys.exit(1)
                line = file.readline()
        except IOError:
            print ("Count not open %s" % confFile)
            sys.exit(2)
        file.close()

    def stddev(self):
        self.evalue = []
        self.evalue2 = []
        for i in range(0, len(self.status)):
            s=self.status[i]
            self.evalue.append(s.evalue_max)
            self.evalue2.append(s.evalue_max2)
        evaluemax_avg = np.sum(self.evalue) / len(self.status)
        std_dev = 1 / np.sqrt(len(self.status)) \
                  * np.sqrt(np.sum(self.evalue2) / len(self.status) \
                  - (np.sum(self.evalue) / len(self.status)) ** 2)
        # print (std_dev)

        file_status = open('status','a')
        file_status.write("\n")
        file_status.write("\n")
        file_status.write(r"#  eigenvalue_max_avg std_dev (1/\sqrt(N)*\sqrt([x^2]-[x]^2))\n")
        file_status.write(
           "\n"
           + str(evaluemax_avg)
           + "\t"
           + str(std_dev)
        )
        file_status.close()

def get_command_line_args(args):
    if len(args) != 2:
        print ("Usage: %s <input files> <output files>"  % args[0])
        sys.exit(1)
    return args

def main():
    args = get_command_line_args(sys.argv)
    lmpconf = LmpConf()
    lmpconf.read_configuration_file(args[1])
    lmpconf.stddev()

if __name__ == "__main__":
    main()
