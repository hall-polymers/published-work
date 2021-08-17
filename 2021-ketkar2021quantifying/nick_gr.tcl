#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#!/usr/bin/tclsh

# This script calculates g(r) for all possible pairs of particles

# Usage: tclsh nick_gr.tcl "$path_to_file"  "$dump_file_name"   "$spacing"     "$maximum_r"    "$output_file_prefix"
#                          [lindex $argv 0] [lindex $argv 1] [lindex $argv 2] [lindex $argv 3] [lindex $argv 4]
#
# Best to call via tclsh, your tcl interpreter may be in some other directory in your $PATH

# package require pbctools
set prefix [lindex $argv 4];  # Prefix for output file name

# Read in molecule
# See: https://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/18844.html
set path [lindex $argv 0];  # The first argument passed should be the relevant path
set dumpf [lindex $argv 1];  # The second arg is the dump file name 
mol new "${path}/${dumpf}" type lammpstrj waitfor all first 0 last -1 step 1 autobonds off 

set dr [lindex $argv 2];  # 3rd arg is spacing
set maxr [lindex $argv 3]

set sel1 [atomselect top "type 3"]
set sel2 [atomselect top "type 4"]
set num1 [$sel1 num]
set num2 [$sel2 num]
puts "# Atoms in selection 1 are $num1"
puts "# Atoms in selection 2 are $num2"
#set rdfcn [measure rdf $sel1 $sel2 delta $dr rmax $maxr usepbc 1 first 0 last -1 step 1]; # GPU code
set rdfcn [measure gofr $sel1 $sel2 delta $dr rmax $maxr usepbc 1 first 0 last -1 step 1]; # CPU code
# Note: Normalization differs b/twn the CPU code "gofr" and GPU code "rdfcn". Normalization in the
#       GPU code normalization is based on average density of particles found inside the various spheres of
#       radius maxr (centered on each origin atom) rather than the entire box. CPU code normalizes using
#       avg density from the entire box.
set r_vals [lindex $rdfcn 0]; # Pull out nested list of positions
set g_vals [lindex $rdfcn 1]; # Pull out nested list of g-values
set n_vals [lindex $rdfcn 2]; # Pull out nested list of n-values

# Writing to systematically named files...
set fname "${prefix}.txt"
set rdf_out [open "${path}/${fname}" "w"]

foreach x_vals $r_vals y_vals $g_vals z_vals $n_vals {
    puts $rdf_out "${x_vals}\t${y_vals}\t${z_vals}"
}
close $rdf_out

# Calculate Average Density
# set rho_fname "${prefix}_density.dat"
# puts "Computing density for selection"
# set average_rho [expr $num2 / ${average_V}];  # Calculate average density
# set rho_out [open "$rho_fname" "w+"]; # Writes, but if file doesn't exist it creates an empty file
# puts $rho_out "${average_rho}";  # tcl always adds newline char on its own
# close $rho_out

exit

