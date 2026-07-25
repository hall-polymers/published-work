#!/usr/bin/python
# -*- coding: utf-8 -*-

# %% Overview
# This Python 3 script generates a frame file with unscaled and unwrapped
# coordinates from a frame file with scaled and wrapped coordinates generated
# by frame_file_generation.py.

# %% Modification History
# This script has been edited by:
# * Felipe Fabricio Pacci Evaristo (12/31/2023)

# %% Library Imports
from sys import argv, exit


# %% Functions
# Retrieval of the Arguments used with this Script
def get_command_arguments(arguments):
    if len(arguments) != 3:
        print("Unexpected number or arguments.\n"
              "Usage: %s "
              "<scaled_wrapped_coordinates_frame_file> "
              "<unscaled_unwrapped_coordinates_frame_file>"
              % argv[0])
        exit(1)
    return arguments


def main():
    arguments = get_command_arguments(argv)
    scaled_wrapped_coordinates_frame_file = arguments[1]
    unscaled_unwrapped_coordinates_frame_file = arguments[2]
    frame = Frame()
    frame.read_scaled_wrapped_coordinates_frame_file(
        scaled_wrapped_coordinates_frame_file)
    frame.unscale_unwrap_coordinates()
    frame.write_unscaled_unwrapped_coordinates_frame_file(
        unscaled_unwrapped_coordinates_frame_file)


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
        self.xu = None
        self.yu = None
        self.zu = None


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

    # Calculation of the Relevant Parameters that define the Simulation Box
    def calculate_box_parameters(self):
        self.boxEdgeLength.append(abs(self.xhi - self.xlo))
        self.boxEdgeLength.append(abs(self.yhi - self.ylo))
        self.boxEdgeLength.append(abs(self.zhi - self.zlo))

    # Import of the Frame File with Scaled and Wrapped Coordinates
    def read_scaled_wrapped_coordinates_frame_file(
            self, scaled_wrapped_coordinates_frame_file):
        try:
            file = open(scaled_wrapped_coordinates_frame_file, "r")
        except IOError:
            print("Scaled wrapped coordinates frame file %s "
                  "could not be opened."
                  % scaled_wrapped_coordinates_frame_file)
            exit(1)
        else:
            line = file.readline()
            while line:
                if "Timestep" in line:
                    self.timestep = int(line.strip().split()[-1])
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
            self.calculate_box_parameters()
        file.close()

    # Unscaling and Unwrapping of the Coordinates of the Imported Frame File
    # with Scaled and Wrapped Coordinates
    def unscale_unwrap_coordinates(self):
        for index in range(self.numberAtoms):
            self.atoms[index].xu = \
                (self.atoms[index].xs + self.atoms[index].ix) * \
                self.boxEdgeLength[0] + \
                self.xlo
            self.atoms[index].yu = \
                (self.atoms[index].ys + self.atoms[index].iy) * \
                self.boxEdgeLength[1] + \
                self.ylo
            self.atoms[index].zu = \
                (self.atoms[index].zs + self.atoms[index].iz) * \
                self.boxEdgeLength[2] + \
                self.zlo

    # Export of the Frame File with Unscaled and Unwrapped Coordinates
    def write_unscaled_unwrapped_coordinates_frame_file(
            self, unscaled_unwrapped_coordinates_frame_file):
        file = open(unscaled_unwrapped_coordinates_frame_file, "w")
        file.write("Frame at Timestep %d\n" % self.timestep)
        file.write("\n")
        file.write("%6i atoms\n" % self.numberAtoms)
        file.write("%6i bonds\n" % self.numberBonds)
        file.write("\n")
        file.write("%6i atom types\n" % self.numberAtomTypes)
        file.write("%6i bond types\n" % self.numberBondTypes)
        file.write("\n")
        file.write("%18.16f %18.16f xlo xhi\n" % (self.xlo, self.xhi))
        file.write("%18.16f %18.16f ylo yhi\n" % (self.ylo, self.yhi))
        file.write("%18.16f %18.16f zlo zhi\n" % (self.zlo, self.zhi))
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
                "%6i %4i %1i %11.6f %11.6f %11.6f\n"
                % (
                    atom_data.id,
                    atom_data.mol,
                    atom_data.type,
                    atom_data.xu,
                    atom_data.yu,
                    atom_data.zu
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
