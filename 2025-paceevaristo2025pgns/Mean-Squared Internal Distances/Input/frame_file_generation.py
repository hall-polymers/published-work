#!/usr/bin/python
# -*- coding: utf-8 -*-

# %% Overview
# This Python 3 script generates a separate frame file for a select frame of a
# LAMMPS trajectory file in the format of a provided template file.

# %% Modification History
# This script has been edited by:
# * Diego Becerra (06/05/2022)
# * Felipe Fabricio Pacci Evaristo (12/31/2023): Revision of entire script.
# Optimization of redundant processes. Correction of minor bugs. Adequacy to
# PEP 8 coding conventions. Addition of thorough commentary.

# %% Library Imports
from sys import argv, exit


# %% Functions
# Retrieval of the Arguments used with this Script
def get_command_arguments(arguments):
    if len(arguments) != 3:
        print("Unexpected number or arguments.\n"
              "Usage: %s "
              "<template_file> "
              "<trajectory_file>"
              % argv[0])
        exit(1)
    return arguments


# Calculation of the Relevant Parameters that define the Simulation Box
def calculate_box_parameters(frame):
    frame.boxEdgeLength.append(abs(frame.xhi - frame.xlo))
    frame.boxEdgeLength.append(abs(frame.yhi - frame.ylo))
    frame.boxEdgeLength.append(abs(frame.zhi - frame.zlo))


# Import of the LAMMPS Trajectory File
def read_trajectory_file(trajectory_file):
    try:
        file = open(trajectory_file, "r")
    except IOError:
        print("Trajectory file %s could not be opened." % trajectory_file)
        exit(1)
    else:
        frame = Frame()
        frames = []
        line = file.readline()
        while line:
            if "ITEM: TIMESTEP" in line:
                frame = Frame()
                frame.timestep = int(file.readline().strip().split()[0])
            if "ITEM: NUMBER OF ATOMS" in line:
                frame.numberAtoms = int(file.readline().strip().split()[0])
            if "ITEM: BOX BOUNDS" in line:
                line = line.strip().split()
                x_type = line[-3]
                y_type = line[-2]
                z_type = line[-1]
                frame.typesBoundaries = \
                    "%s %s %s" % (x_type, y_type, z_type)
                subline = file.readline().strip().split()
                frame.xlo = float(subline[0])
                frame.xhi = float(subline[1])
                subline = file.readline().strip().split()
                frame.ylo = float(subline[0])
                frame.yhi = float(subline[1])
                subline = file.readline().strip().split()
                frame.zlo = float(subline[0])
                frame.zhi = float(subline[1])
            if "ITEM: ATOMS" in line:
                frame.atomInfo = \
                    line.strip().replace("ITEM: ATOMS ", "").split()
                for index in range(frame.numberAtoms):
                    frame.atomData.append(file.readline().strip().split())
                calculate_box_parameters(frame)
                frames.append(frame)
            line = file.readline()
    file.close()
    return frames


def main():
    arguments = get_command_arguments(argv)
    template_file = arguments[1]
    trajectory_file = arguments[2]
    template = Template()
    template.read_template_file(template_file)
    frames = read_trajectory_file(trajectory_file)
    print("Generating frame files with scaled and wrapped coordinates.")
    for frame in frames:
        print("Timestep: %d" % frame.timestep)
        template.update_template_file(frame)
        timestep_index = template_file.find(".dat")
        frame_file = (
                template_file[:timestep_index] +
                "_" + str(frame.timestep) +
                template_file[timestep_index:]
        )
        template.write_frame_file(frame, frame_file)


# %% Classes
class Atom:
    def __init__(self, id_, mol, type_, xs, ys, zs, ix, iy, iz):
        self.id = id_
        self.mol = mol
        self.type = type_
        self.xs = xs
        self.ys = ys
        self.zs = zs
        self.ix = ix
        self.iy = iy
        self.iz = iz


class Bond:
    def __init__(self, id_, type_, a, b):
        self.id = id_
        self.type = type_
        self.a = a
        self.b = b


class Frame:
    def __init__(self):
        self.timestep = None
        self.numberAtoms = None
        self.typesBoundaries = None
        self.xlo = None
        self.xhi = None
        self.ylo = None
        self.yhi = None
        self.zlo = None
        self.zhi = None
        self.atomInfo = None
        self.atomData = []
        self.boxEdgeLength = []


class Template:
    def __init__(self):
        self.numberAtoms = None
        self.numberBonds = None
        self.numberAtomTypes = None
        self.numberBondTypes = None
        self.xlo = None
        self.xhi = None
        self.ylo = None
        self.yhi = None
        self.zlo = None
        self.zhi = None
        self.masses = []
        self.atoms = []
        self.bonds = []
        self.boxEdgeLength = []

    # Import of the Template File
    def read_template_file(self, template_file):
        try:
            file = open(template_file, "r")
        except IOError:
            print("Template file %s could not be opened." % template_file)
            exit(1)
        else:
            line = file.readline()
            while line:
                if "atoms" in line:
                    self.numberAtoms = int(line.strip().split()[0])
                if "bonds" in line:
                    self.numberBonds = int(line.strip().split()[0])
                if "atom types" in line:
                    self.numberAtomTypes = int(line.strip().split()[0])
                if "bond types" in line:
                    self.numberBondTypes = int(line.strip().split()[0])
                if "xlo xhi" in line:
                    self.xlo = float(line.strip().split()[0])
                    self.xhi = float(line.strip().split()[1])
                if "ylo yhi" in line:
                    self.ylo = float(line.strip().split()[0])
                    self.yhi = float(line.strip().split()[1])
                if "zlo zhi" in line:
                    self.zlo = float(line.strip().split()[0])
                    self.zhi = float(line.strip().split()[1])
                if "Masses" in line:
                    file.readline()
                    for index in range(self.numberAtomTypes):
                        subline = file.readline().strip().split()
                        self.masses.append(float(subline[1]))
                if "Atoms" in line:
                    file.readline()
                    for index in range(self.numberAtoms):
                        subline = file.readline().strip().split()
                        if len(subline) == 9:
                            self.atoms.append(Atom(
                                int(subline[0]),
                                int(subline[1]),
                                int(subline[2]),
                                float(subline[3]),
                                float(subline[4]),
                                float(subline[5]),
                                int(subline[6]),
                                int(subline[7]),
                                int(subline[8])
                            )
                            )
                        else:
                            print("Unexpected number of entries "
                                  "in Atoms line %d."
                                  % (index + 1))
                            exit(2)
                if "Bonds" in line:
                    file.readline()
                    for index in range(self.numberBonds):
                        subline = file.readline().strip().split()
                        self.bonds.append(Bond(
                            int(subline[0]),
                            int(subline[1]),
                            int(subline[2]),
                            int(subline[3])
                        )
                        )
                line = file.readline()
            calculate_box_parameters(self)
        file.close()

    # Update of the Template File with the Data in the Selected Frame of the
    # LAMMPS Trajectory File
    def update_template_file(self, frame):
        atom_info = frame.atomInfo
        atom_info_dictionary = {}
        for index in range(len(atom_info)):
            atom_info_dictionary[atom_info[index]] = index
        for index in range(frame.numberAtoms):
            atom_data = frame.atomData[index]
            atom_index = int(atom_data[atom_info_dictionary["id"]])
            self.atoms[atom_index - 1].xs = \
                float(atom_data[atom_info_dictionary["xs"]])
            self.atoms[atom_index - 1].ys = \
                float(atom_data[atom_info_dictionary["ys"]])
            self.atoms[atom_index - 1].zs = \
                float(atom_data[atom_info_dictionary["zs"]])
            self.atoms[atom_index - 1].ix = \
                float(atom_data[atom_info_dictionary["ix"]])
            self.atoms[atom_index - 1].iy = \
                float(atom_data[atom_info_dictionary["iy"]])
            self.atoms[atom_index - 1].iz = \
                float(atom_data[atom_info_dictionary["iz"]])

    # Export of the Frame File
    def write_frame_file(self, frame, frame_file):
        file = open(frame_file, "w")
        file.write("Frame at Timestep %d\n" % frame.timestep)
        file.write("\n")
        file.write("%6i atoms\n" % self.numberAtoms)
        file.write("%6i bonds\n" % self.numberBonds)
        file.write("\n")
        file.write("%6i atom types\n" % self.numberAtomTypes)
        file.write("%6i bond types\n" % self.numberBondTypes)
        file.write("\n")
        file.write("%18.16f %18.16f xlo xhi\n" % (frame.xlo, frame.xhi))
        file.write("%18.16f %18.16f ylo yhi\n" % (frame.ylo, frame.yhi))
        file.write("%18.16f %18.16f zlo zhi\n" % (frame.zlo, frame.zhi))
        file.write("\n")
        file.write("Masses\n")
        file.write("\n")
        for index in range(self.numberAtomTypes):
            file.write("%1i %5.1f\n" % (index + 1, self.masses[index]))
        file.write("\n")
        file.write("Atoms\n")
        file.write("\n")
        for index in range(len(self.atoms)):
            atom_data = self.atoms[index]
            file.write(
                "%6i %4i %1i %11.6f %11.6f %11.6f %2i %2i %2i\n"
                % (
                    atom_data.id,
                    atom_data.mol,
                    atom_data.type,
                    atom_data.xs,
                    atom_data.ys,
                    atom_data.zs,
                    atom_data.ix,
                    atom_data.iy,
                    atom_data.iz
                )
            )
        file.write("\n")
        if len(self.bonds) > 0:
            file.write("Bonds\n")
            file.write("\n")
            for index in range(len(self.bonds)):
                bond_data = self.bonds[index]
                file.write(
                    "%6i %1i %6i %6i\n"
                    % (
                        bond_data.id,
                        bond_data.type,
                        bond_data.a,
                        bond_data.b
                    )
                )
        file.close()


if __name__ == "__main__":
    main()
