#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that show abundances Vs position in log-log for a given list of species

import os, pdb
import numpy as np
import sys # to be able to retrieve arguments of the script
import pylab as pl
from matplotlib.ticker import FormatStrFormatter

###############################################
## Beginning of the program
###############################################

AB_FOLDER = 'ab' # Name of the sub-folder containing abundances informations for each species
OUTPUT_EXTENSION = 'pdf' # default value in bitmap, because vectoriel can take time and space if there is a lot of data

# By default, 1, but represent the column containing abundances in each file. 
index_time = 1

species_name = None

isProblem = False
problem_message = """AIM : Display in log-log the spatial distribution
of abundances for a set of species.
The script must be launched in a folder that contain a simulation.

The script can take various arguments :
(no spaces between the key and the values, only separated by '=')
 * xmax=20 : the end of the output [a.u.]
 * xmin=5000 : the beginning of the output [a.u.]
 * t=%d : index of the desired time point
 * species=CO,H20 : the list of species we want to display /!\ no space !
 * ext=%s : The extension for the output files
 * help : display a little help message on HOW to use various options.

EXAMPLE:
 > nautilus-plot-abundances.py species=CO,H20,C2H6
 > nautilus-plot-abundances.py species=CO,H20 xmin=1. xmax=1e6""" % (index_time, OUTPUT_EXTENSION)


value_message = "/!\ Warning: %s does not need any value, but you defined '%s=%s' ; value ignored."

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
    value = None
  if (key == 'xmin'):
    x_min = float(value)
  elif (key == 'xmax'):
    x_max = float(value)
  elif (key == 't'):
    index_time = int(value)
  elif (key == 'ext'):
    OUTPUT_EXTENSION = value
  elif (key == 'species'):
    species_name = value.split(',')
  elif (key == 'help'):
    isProblem = True
    if (value != None):
      print(value_message % (key, key, value))
  else:
    print("the key '%s' does not match" % key)
    isProblem = True

if (species_name == None):
  isProblem = True

if isProblem:
  print(problem_message)
  exit()

####################
# We declare arrays to store read informations
####################
ref_pos = [] # position in a.u.
abundances = {} # Array containing Species abundances (relative to H) [number ratio].

# For each species, we read the corresponding file and store its values
for name in species_name:
  filename = "%s.ab" % name
  filePath = os.path.join(AB_FOLDER,filename)
  if (not(os.path.isfile(filePath))):
    raise NameError("The file %s doesn't exists in %s/" % (filename, AB_FOLDER))
  
  try:
    tmp = np.loadtxt(filePath, skiprows=1, dtype=float, unpack=True)[:,index_time]
    time_value = tmp[0]
    tmp_abundance = tmp[1:]
  except IndexError:
    print("Spatial point out of bounds. Please lower your value of 't'.")
    exit()
  
  filename = "%s.ab" % 'space'
  filePath = os.path.join(AB_FOLDER,filename)
  if (not(os.path.isfile(filePath))):
    raise NameError("The file %s doesn't exists in %s/" % (filename, AB_FOLDER))
  
  try:
    ref_pos = np.loadtxt(filePath, skiprows=1, dtype=float, unpack=True)
  except IndexError:
    print("Spatial point out of bounds. Please lower your value of 't'.")
    exit()
  
  abundances[name] = tmp_abundance


# We get the index for the x_max value
if not('x_max' in locals()):
  x_max = ref_pos[-1]

# We get the index for the x_min value
if not('x_min' in locals()):
  x_min = ref_pos[0]
  
point_selection = (ref_pos >= x_min) & (ref_pos <= x_max)

# We generate a list of colors
#~ colors = [ '#'+li for li in autiwa.colorList(nb_planete)]

fig = pl.figure()
pl.clf()
fig.subplots_adjust(left=0.12, bottom=0.1, right=0.96, top=0.95, wspace=0.26, hspace=0.26)
# We create subplots. add_subplot(2, 3, 1) means we have 2 lines, 3 columns, 
# and that the active plot is the first, starting from top left (for 6 plots in total)
plot_ab = fig.add_subplot(1, 1, 1)

plot = plot_ab.loglog

for species in species_name:
  plot(ref_pos[point_selection], abundances[species][point_selection], label=species)

plot_ab.set_xlabel("Position [a.u.]")
plot_ab.set_ylabel("Abundance")
plot_ab.set_title("Time= {:5.3e} yr".format(time_value))
plot_ab.grid(True)

timeFormat = FormatStrFormatter("%.3g")
plot_ab.xaxis.set_major_formatter(timeFormat)

plot_ab.legend(loc="best")


nom_fichier_plot = "abundances_1Dstatic_timeidx_{}".format(int(index_time))
fig.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)

#~ pdb.set_trace()

pl.show()



