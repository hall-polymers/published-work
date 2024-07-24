#!/usr/bin/env python

import sys
#from lmp_io import *
class DumpFrame:
    def __init__(self):
        self.timestep = None
        self.numberAtoms = None
        self.boxBounds = None
        self.xlo = None
        self.xhi = None
        self.ylo = None
        self.yhi = None
        self.zlo = None
        self.zhi = None
        self.atomInfo = None
        self.atomData = []

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
    def __init__(self,type,a,b):
        self.type = type
        self.a = a
        self.b = b

class Angle:
    def __init__(self,type,a,b,c):
        self.type = type
        self.a = a
        self.b = b
        self.c = c

class Dihedral:
    def __init__(self,type,a,b,c,d):
        self.type = type
        self.a = a
        self.b = b
        self.c = c
        self.d = d

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
        self.nDihedrals = None
        self.nAtomTypes = None
        self.nBondTypes = None
        self.nAngleTypes = None
        self.nDihedralTypes = None
        self.mass = []
        self.atoms = []
        self.ellipsoids = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.blen = []
        self.bhlf = []
        self.blni = []

    def initialize_box(self):
        self.blen.append(abs(self.xhi - self.xlo))
        self.blen.append(abs(self.yhi - self.ylo))
        self.blen.append(abs(self.zhi - self.zlo))
        for i in range(3):
            self.bhlf.append(self.blen[i] / 2.0)
            self.blni.append(1.0 / self.blen[i])

    def apply_pbcs(self, dr, dim):
        if dr > self.bhlf[dim]:
            dr = dr - self.blen[dim] * float(int(dr * self.blni[dim] + 0.5))
        elif dr < - self.bhlf[dim]:
            dr = dr - self.blen[dim] * float(int(dr * self.blni[dim] - 0.5))
        return dr

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
                if "Pair Coeffs" in line:
                    print ("This script is not designed to read Pair Coeffs data!")
                    sys.exit(1)
                if "Bond Coeffs" in line:
                    file.readline()
                    for i in range(self.nBondTypes):
                        l = file.readline()
                        self.bonddata.append(l)
                if "Angle Coeffs" in line:
                    file.readline()
                    for i in range(self.nAngleTypes):
                        l = file.readline()
                        self.angledata.append(l)
                if "Dihedral Coeffs" in line:
                    file.readline()
                    for i in range(self.nDihedralTypes):
                        l = file.readline()
                        self.dihedraldata.append(l)
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
                if "Dihedrals" in line:
                    file.readline()
                    for i in range(self.nDihedrals):
                        l = file.readline().strip().split()
                        self.dihedrals.append(Dihedral(
                            int(l[1]),
                            int(l[2]),
                            int(l[3]),
                            int(l[4]),
                            int(l[5])
                        )
                        )
                line = file.readline()
        except IOError:
            print ("Count not open %s" % confFile)
            sys.exit(2)
        self.initialize_box()
        file.close()

    def read_dump_file(self, dumpFile):
        try:
            file = open(dumpFile, "r")
        except IOError:
            print ("%s cannot be opened!" % dumpFile)
        else:
            localDumps = []
            line = file.readline()
            while line:
                if "ITEM: TIMESTEP" in line:
                    frame = DumpFrame()
                    frame.timestep = int(file.readline().strip().split()[0])
                    self.timestep =  frame.timestep
                if "ITEM: NUMBER OF ATOMS" in line:
                    frame.numberAtoms = int(file.readline().strip().split()[0])
                    number = frame.numberAtoms
                if "ITEM: BOX BOUNDS" in line:
                    line = line.strip().split()
                    b1 = line[-3]
                    b2 = line[-2]
                    b3 = line[-1]
                    frame.boxBounds = "%s %s %s" % (b1, b2, b3)
                    line = file.readline().strip().split()
                    frame.xlo = float(line[0])
                    frame.xhi = float(line[1])
                    line = file.readline().strip().split()
                    frame.ylo = float(line[0])
                    frame.yhi = float(line[1])
                    line = file.readline().strip().split()
                    frame.zlo = float(line[0])
                    frame.zhi = float(line[1])
                if "ITEM: ATOMS" in line:
                    frame.atomInfo = line.strip().replace(
                    "ITEM: ATOMS ","").split()
                    for i in range(frame.numberAtoms):
                        frame.atomData.append(file.readline().strip().split())
                    localDumps.append(frame)
                line = file.readline()
            file.close()
            return localDumps[<traj>]

    def update_configuration_file(self,dumpData,confFile):
        # First, build a map from dump to conf
        siteMap = {}
        # Sense strand
        for i in range(dumpData.numberAtoms+1):
            siteMap[i] = i
        # First, reorder the data in dumpData
        info = dumpData.atomInfo
        infoDict = {}
        for i in range(len(info)):
            infoDict[info[i]]=  int(i)
        final = dumpData.numberAtoms
        for i in range (final):
            data = dumpData.atomData[i]
            index = int(data[infoDict["id"]])
            if index in siteMap:
                idx = siteMap[index]
                self.atoms[idx - 1].x = float(data[infoDict["xu"]])
                self.atoms[idx - 1].y = float(data[infoDict["yu"]])
                self.atoms[idx - 1].z = float(data[infoDict["zu"]])
                self.ellipsoids[idx-1].quatw = float(data[infoDict["c_orient[4]"]])
                self.ellipsoids[idx - 1].quati = float(data[infoDict["c_orient[1]"]])
                self.ellipsoids[idx - 1].quatj = float(data[infoDict["c_orient[2]"]])
                self.ellipsoids[idx - 1].quatk = float(data[infoDict["c_orient[3]"]])

    def write_configuration_file(self,outFile):
        file = open(outFile,"w")
        file.write("LAMMPS configuration file modified using lmp_io\n\n")
        file.write("\t%ld atoms\n" % self.nAtoms)
        file.write("\t%ld ellipsoids\n" % self.nEllipsoids)
        file.write("\t%ld bonds\n" % self.nBonds)
        file.write("\t%ld angles\n\n" % self.nAngles)
        #file.write("\t%ld dihedrals\n\n" % self.nDihedrals)
        file.write("\t%ld atom types\n" % self.nAtomTypes)
        file.write("\t%ld bond types\n" % self.nBondTypes)
        file.write("\t%ld angle types\n\n" % self.nAngleTypes)
        #file.write("\t%ld dihedral types\n\n" % self.nDihedralTypes)
        file.write("\t%lf\t%lf xlo xhi\n" % (self.xlo, self.xhi))
        file.write("\t%lf\t%lf ylo yhi\n" % (self.ylo, self.yhi))
        file.write("\t%lf\t%lf zlo zhi\n\n" % (self.zlo, self.zhi))
        file.write("Masses\n\n")
        for i in range(self.nAtomTypes):
            file.write("\t%d\t%lf\n" % (i + 1, (self.mass[i])))
        file.write("\n")
        file.write("Atoms\n\n")
        for i in range(len(self.atoms)):
            a = self.atoms[i]
            file.write("\t%ld\t%ld\t%lf\t%lf\t%lf\t%ld\t%lf\t%ld\t%f\t%ld\t%ld\t%ld\n" %
            (a.id, a.type, self.atoms[i].x, self.atoms[i].y, self.atoms[i].z,
            a.mol, a.q, a.ellipsoidflag, a.density, 0, 0, 0)
            )
        file.write("\n")
        file.write("Ellipsoids\n\n")
        for i in range(len(self.ellipsoids)):
            e = self.ellipsoids[i]
            file.write("\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n" %
            (e.id,
            self.ellipsoids[i].shapex,
            self.ellipsoids[i].shapey,
            self.ellipsoids[i].shapez,
            self.ellipsoids[i].quatw,
            self.ellipsoids[i].quati,
            self.ellipsoids[i].quatj,
            self.ellipsoids[i].quatk)
            )
        file.write("\n")
        if len(self.bonds) > 0:
           file.write("Bonds\n\n")
           for i in range(len(self.bonds)):
               b = self.bonds[i]
               file.write("\t%ld\t%ld\t%ld\t%ld\n" %
               (i + 1, b.type, b.a, b.b))
           file.write("\n")
        if len(self.angles) > 0:
           file.write("Angles\n\n")
           for i in range(len(self.angles)):
               ang = self.angles[i]
               file.write("\t%ld\t%ld\t%ld\t%ld\t%ld\n" %
               (i+1, ang.type, ang.a, ang.b, ang.c))
           file.write("\n")
        # if len(self.dihedrals) > 0:
        #    file.write("Dihedrals\n\n")
        #    for i in range(len(self.dihedrals)):
        #        d = self.dihedrals[i]
        #        file.write("\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n" % (i+1,d.type,d.a,d.b,d.c,d.d))
        #    file.close()

def get_command_line_args(args):
    if len(args) != 4:
        print ("Usage: %s <confFile> <dumpFile> <modifiedConfFile" % sys.argv[0])
        sys.exit(1)
    return args

def main():
    args = get_command_line_args(sys.argv)
    confFile = args[1]
    dumpFile = args[2]
    outFile = args[3]
    lmpconf = LmpConf()
    lmpconf.read_configuration_file(confFile)
    dumpData = lmpconf.read_dump_file(args[2])
    lmpconf.update_configuration_file(dumpData,confFile)
    lmpconf.write_configuration_file(outFile)

if __name__ == "__main__":
    main()
