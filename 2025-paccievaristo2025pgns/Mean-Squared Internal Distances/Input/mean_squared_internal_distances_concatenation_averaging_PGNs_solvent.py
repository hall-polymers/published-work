#!/usr/bin/python
# -*- coding: utf-8 -*-

# %% Overview
# This Python 3 script concatenates and averages all data files with average
# mean-squared internal distances associated with a simulation exported by
# mean_squared_internal_distances_calculation_PGNs_solvent.py.

# %% Modification History
# This script has been edited by:
# * Felipe Fabricio Pacci Evaristo (12/31/2023)

# %% Library Imports
from os import listdir
from sys import argv, exit
from pandas import DataFrame, read_csv, concat


# %% Functions
# Retrieval of the Arguments used with this Script
def get_command_arguments(arguments):
    if len(arguments) != 5:
        print("Unexpected number or arguments.\n"
              "Usage: %s "
              "<average_mean_squared_internal_distances_data_files_direction> "
              "<average_mean_squared_internal_distances_data_files_path> "
              "<combined_average_mean_squared_internal_distances_data_file> "
              "<ensemble_averaged_mean_squared_internal_distances_data_file>"
              % argv[0])
        exit(1)
    return arguments


# Export of a Data File
def write_data_file(data_frame, data_file_path, data_file):
    with open(data_file_path + data_file, "w") as file:
        string = data_frame.to_string()
        file.write(string)


def main():
    arguments = get_command_arguments(argv)
    average_mean_squared_internal_distances_data_files_direction = arguments[1]
    average_mean_squared_internal_distances_data_files_path = arguments[2]
    combined_average_mean_squared_internal_distances_data_file = arguments[3]
    ensemble_averaged_mean_squared_internal_distances_data_file = arguments[4]
    separate_files = []
    for file in \
            listdir(average_mean_squared_internal_distances_data_files_path):
        if file.startswith(
                "average_mean_squared_internal_distances" +
                "_" +
                average_mean_squared_internal_distances_data_files_direction):
            separate_files.append(file)
    separate_files.sort()
    combined_average_mean_squared_internal_distances = DataFrame()
    for separate_file in separate_files:
        file = open(separate_file, "r")
        timestep = file.readline().split()[-1]
        file.close()
        combined_average_mean_squared_internal_distances = \
            concat([combined_average_mean_squared_internal_distances,
                    read_csv(separate_file,
                             sep=r"\s+",
                             header=0,
                             names=[str(timestep)],
                             usecols=[1])],
                   axis=1)
    ensemble_averaged_mean_squared_internal_distances_series = \
        combined_average_mean_squared_internal_distances.mean(axis=1)
    sd_ensemble_averaged_mean_squared_internal_distances_series = \
        combined_average_mean_squared_internal_distances.std(axis=1)
    ensemble_averaged_mean_squared_internal_distances = \
        ensemble_averaged_mean_squared_internal_distances_series.to_frame()
    sd_ensemble_averaged_mean_squared_internal_distances = \
        sd_ensemble_averaged_mean_squared_internal_distances_series.to_frame()
    ensemble_averaged_mean_squared_internal_distances = \
        concat([ensemble_averaged_mean_squared_internal_distances,
                sd_ensemble_averaged_mean_squared_internal_distances],
               axis=1)
    ensemble_averaged_mean_squared_internal_distances.columns = \
        ["Average", "Standard Deviation"]
    write_data_file(
        combined_average_mean_squared_internal_distances,
        average_mean_squared_internal_distances_data_files_path,
        combined_average_mean_squared_internal_distances_data_file)
    write_data_file(
        ensemble_averaged_mean_squared_internal_distances,
        average_mean_squared_internal_distances_data_files_path,
        ensemble_averaged_mean_squared_internal_distances_data_file)


if __name__ == "__main__":
    main()
