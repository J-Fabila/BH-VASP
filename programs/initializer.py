import subprocess
import glob
import os
import math
import random
import re
import sys
import datetime
from subprocess import*

# Variable declarations
init = ""
input_dir = ""
programs_dir = ""
substrate_nat = ""
cluster_nat = ""
cluster_ntyp = ""
temperature_K = 0
step_width = ""
iterations = ""
x_range = ""
y_range = ""
z_range = ""
z_vacuum = ""

time_initial = 0
time_final = 0

# Initial time
time_initial = datetime.datetime.now()

file = str(sys.argv[1])

# Variable initialization with information from command line
f = open(file, 'r')
lines = f.readlines()
f.close()

for line in lines:
	l = line.strip()
	if l.startswith("initialization_file"):
		init = l.split("initialization_file = ", 1)[1]
	if l.startswith("input_dir"):
		input_dir = l.split("input_dir = ", 1)[1]
	if l.startswith("programs_dir"):
		programs_dir = l.split("programs_dir = ", 1)[1]
	if l.startswith("substrate_nat"):
		substrate_nat = l.split("substrate_nat = ", 1)[1]
	if l.startswith("cluster_ntyp"):
		cluster_ntyp = l.split("cluster_ntyp = ", 1)[1]
	if l.startswith("temperature_K"):
		temperature_K = float(l.split("temperature_K = ", 1)[1])
	if l.startswith("step_width"):
		step_width = l.split("step_width = ", 1)[1]
	if l.startswith("iterations"):
		iterations = int(l.split("iterations = ", 1)[1])
	if l.startswith("x_range"):
		x_range = l.split("x_range = ", 1)[1]
	if l.startswith("y_range"):
		y_range = l.split("y_range = ", 1)[1]
	if l.startswith("z_range"):
		z_range = l.split("z_range = ", 1)[1]
	if l.startswith("z_vacuum"):
		z_vacuum = l.split("z_vacuum = ", 1)[1]

sys.path.append(programs_dir)

from input_parser import*
# Cluster total number of atoms
cluster_nat = str(get_nAtoms_str(cluster_ntyp))

# Temperature
kBT = float(0.00008617 * temperature_K)

# Output folder name
output_dir_name = get_output_name(cluster_ntyp)

# Output folder is created
call(['cp', '-r', input_dir, './' + output_dir_name])




print("=========================================================================================================")
print("Basin Hopping Monte Carlo global optimization program coupled to (DFT-Quantum Espresso) ")
print("for the analysis of supported clusters, and (eventually) surface/bulk alloys ")
print("Instituto de Fisica - UNAM (2016) ")
print("=========================================================================================================")
print(" ")
print(" --> Starting DFT relaxation of initial/random configuration ")
print("Iteration 1 of " + output_dir_name)
print(" ")


