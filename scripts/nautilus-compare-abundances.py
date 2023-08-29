#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that show abundances Vs time in log-log for a given list of species

import os, pdb
import numpy as np
import sys # to be able to retrieve arguments of the script
import pylab as pl
from matplotlib.ticker import FormatStrFormatter

import colorsys

def get_N_HexCol(number):
  """Generate equally spaced colors using HUE color space 
  (that I don't understand though, since I picked two lines)."""
  ref_color = ['#ff0000', '#00ff00', '#0000ff']
  nb_ref = len(ref_color)
  if (number <= nb_ref):
    RGB_tuples = ref_color[0:number]
  else:
    RGB_tuples = ref_color
    
    remaining = number - nb_ref
    HSV_tuples = [(x*1.0/remaining, 0.5, 0.5) for x in range(remaining)]
    RGB_tuples.extend(map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples))

  return RGB_tuples

###############################################
## Beginning of the program
###############################################
linestyles = ['-', '--', '-.', ':', '.', ',', 'o', 'v', '^', '<', '>', '1', '2']

AB_FOLDER = 'ab' # Name of the sub-folder containing abundances informations for each species
OUTPUT_EXTENSION = 'pdf' # default value in bitmap, because vectoriel can take time and space if there is a lot of data
nom_fichier_plot = "compare_abundances"

directories = [dir for dir in os.listdir(".") if os.path.isdir(dir)] # We list folders of the current directory, each one assumed to be a simulation
directories.sort()

# By default, 1, but represent the column containing abundances in each file. 
## If one, this is either the first 1D point, or the only one, if this is not a 1D simulation
index_1D = 1

isProblem = False
problem_message = """AIM : Display in log-log the evolution of abundances for a set of species. 
Must be launched in a folder that contains sub-folders, each one containing a simulation.

The script can take various arguments:
(no spaces between the key and the values, only separated by '=')
 * tmax=1.e6 : the end of the output [year]
 * tmin=5e5 : the beginning of the output [year]
 * x=%d : index of the desired spatial point
 * species=CO,H20 : the list of species we want to display /!\ no space !
 * dir=simu0001,simu0002 : the list of sub-folder we want to display /!\ no space !
 * ext=%s : The extension for the output files
 * help : display a little help message on HOW to use various options.

EXAMPLE:
 > nautilus-compare-abundances.py species=CO,H20,C2H6
 > nautilus-compare-abundances.py species=CO,H20,C2H6 dir=simu1,simu2
 > nautilus-compare-abundances.py species=CO,H20,C2H6 x=2
 > nautilus-compare-abundances.py species=CO,H20 tmin=1. tmax=1e6""" % (index_1D, OUTPUT_EXTENSION)


value_message = "/!\ Warning: %s does not need any value, but you defined '%s=%s' ; value ignored."

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
    value = None
  if (key == 'tmin'):
    t_min = float(value)
  elif (key == 'tmax'):
    t_max = float(value)
  elif (key == 'x'):
    index_1D = int(value)
  elif (key == 'ext'):
    OUTPUT_EXTENSION = value
  elif (key == 'species'):
    species_name = value.split(',')
  elif (key == 'dir'):
    directories = value.split(',')
  elif (key == 'help'):
    isProblem = True
    if (value != None):
      print(value_message % (key, key, value))
  else:
    print("the key '%s' does not match" % key)
    isProblem = True

if isProblem:
  print(problem_message)
  exit()

if (len(directories)> len(linestyles)):
  print("Error: Not enough linestyles to plot the various simulations")
  exit()
  
if not('species_name' in locals()):
  tmp = raw_input("Enter species names, without space, separated by ',' (e.g.: 'CO,H2O'):\n")
  species_name = tmp.split(',')

####################
# We declare arrays to store read informations
####################
ref_time = [] # time in years, one sub-list per simulation
abundances = [] # array of dictionnaries, each one containing Species abundances (relative to H) [number ratio].

nb_species = len(species_name)

for (index, folder_name) in enumerate(directories):
  os.chdir(folder_name)
  
  abundances.append({})
  
  # For each species, we read the corresponding file and store its values
  for name in species_name:
    filename = "%s.ab" % name
    filePath = os.path.join(AB_FOLDER,filename)
    if (not(os.path.isfile(filePath))):
      raise NameError("The file %s doesn't exists in %s/" % (filename, AB_FOLDER))
    
    try:
      (tmp_time, tmp_abundance) = np.loadtxt(filePath, skiprows=1, usecols=(0,index_1D), dtype=float, unpack=True)
    except IndexError:
      print("Spatial point out of bounds. Please lower your value of 'x'.")
      exit()

    abundances[-1][name] = tmp_abundance
  ref_time.append(tmp_time) # only one array per simulation, regardless of the number of species
  
  os.chdir("..")

# We get the index for the t_max value
if not('t_max' in locals()):
  tmp = [t[-1] for t in ref_time]
  t_max = max(tmp)

# We get the index for the t_min value
if not('t_min' in locals()):
  tmp = [t[0] for t in ref_time]
  t_min = min(tmp)
  

point_selection = []
for time in ref_time:
  tmp_selection = (time >= t_min) & (time <= t_max)
  point_selection.append(tmp_selection)
  


# We generate a list of colors
colors = get_N_HexCol(nb_species)

fig = pl.figure()
pl.clf()
fig.subplots_adjust(left=0.12, bottom=0.1, right=0.96, top=0.95, wspace=0.26, hspace=0.26)
# We create subplots. add_subplot(2, 3, 1) means we have 2 lines, 3 columns, 
# and that the active plot is the first, starting from top left (for 6 plots in total)
plot_ab = fig.add_subplot(1, 1, 1)

plot = plot_ab.loglog
for (folder_name, time, ab, linestyle, selection) in zip(directories, ref_time, abundances, linestyles, point_selection):
  
  for (species,color) in zip(species_name, colors):
    plot(time[selection], ab[species][selection], color=color, linestyle=linestyle, label=species)

handles = [] # Style of line to fill legend
labels = [] # corresponding name to fill legend

for (species, color) in zip(species_name, colors):
  tmpLine = pl.Line2D((0,1),(0,0), color=color, linestyle='-') # Straight line but various colors
  handles.append(tmpLine)
  labels.append(species)

for (folder_name, linestyle) in zip(directories, linestyles):
  tmpLine = pl.Line2D((0,1),(0,0), color='k', linestyle=linestyle) # Black line, but various linestyles
  handles.append(tmpLine)
  labels.append(folder_name)

plot_ab.set_xlabel("Time [years]")
plot_ab.set_ylabel("Abundance")
plot_ab.grid(True)

timeFormat = FormatStrFormatter("%.3g")
plot_ab.xaxis.set_major_formatter(timeFormat)

plot_ab.legend(handles, labels, loc="best")


fig.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)

pl.show()



