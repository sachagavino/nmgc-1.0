#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that show abundances Vs time in log-log for a given list of species

import os, pdb
import numpy as np
import sys # to be able to retrieve arguments of the script
import pylab as pl
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter

###############################################
## Beginning of the program
###############################################

OUTPUT_EXTENSION = 'pdf' # default value in bitmap, because vectoriel can take time and space if there is a lot of data

# By default, 1, but represent the column containing abundances in each file. 
## If one, this is either the first 1D point, or the only one, if this is not a 1D simulation

species_name = None

isProblem = False
problem_message = """AIM : For a given species, display the most important reaction
both for production and destruction. Display the evolution of theses importances
with time.

The script can take various arguments :
(no spaces between the key and the values, only separated by '=')
 * species=CO : the species name you want to plot the most important reactions
 * ext=%s : The extension for the output files
 * help : display a little help message on HOW to use various options.

EXAMPLE:
 > nautilus-trace-species.py species=CO""" % (OUTPUT_EXTENSION)


value_message = "/!\ Warning: %s does not need any value, but you defined '%s=%s' ; value ignored."

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
    value = None
  if (key == 'ext'):
    OUTPUT_EXTENSION = value
  elif (key == 'species'):
    species_name = value
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

filename = 'trace_dest_%s.percentage' % species_name
if not(os.path.isfile(filename)):
  print("Error: The file %s doesn't exists." % filename)
  print("Please run nautilus_trace_major for the given reaction beforehand.")
  exit()
  
time = []
production_fraction = []
destruction_fraction = []

# Read production reactions
object_file = open('trace_prod_%s.percentage' % species_name, 'r')

# First header
dummy = object_file.readline()

# reaction IDs
line = object_file.readline()
reaction_prod_IDs = [word for word in line.split()]

nb_prod_reactions = len(reaction_prod_IDs)
for i in range(nb_prod_reactions):
  production_fraction.append([])
# 2nd header
dummy = object_file.readline()

# Actual datas
lines = object_file.readlines()

for line in lines:
  words = line.split()
  
  time.append(float(words[0]))
  for i in range(nb_prod_reactions):
    production_fraction[i].append(float(words[i+1]))

object_file.close()

# Read destruction reactions

time = [] # We reset the time array because we assume the same times are valid both for production AND destruction !
object_file = open('trace_dest_%s.percentage' % species_name, 'r')

# First header
dummy = object_file.readline()

# reaction IDs
line = object_file.readline()
reaction_dest_IDs = [word for word in line.split()]

nb_dest_reactions = len(reaction_dest_IDs)
for i in range(nb_dest_reactions):
  destruction_fraction.append([])
# 2nd header
dummy = object_file.readline()

# Actual datas
lines = object_file.readlines()

for line in lines:
  words = line.split()
  
  time.append(float(words[0]))
  for i in range(nb_dest_reactions):
    destruction_fraction[i].append(float(words[i+1]))

object_file.close()

#~ pdb.set_trace()

fig = pl.figure()
pl.clf()
fig.suptitle("Main production and destruction reactions for %s" % species_name)

#~ fig.subplots_adjust(left=0.12, bottom=0.1, right=0.96, top=0.95, wspace=0.26, hspace=0.26)
# We create subplots. add_subplot(2, 3, 1) means we have 2 lines, 3 columns, 
# and that the active plot is the first, starting from top left (for 6 plots in total)

# Plot production
plot_prod = fig.add_subplot(2, 1, 1)

plot_prod.set_xscale("log")
plot = plot_prod.plot

for reaction in range(nb_prod_reactions):
  plot(time, production_fraction[reaction], label="ID=%s" % reaction_prod_IDs[reaction])

plot_prod.set_xlabel("Time [years]")
plot_prod.set_ylabel("Production fraction [%]")
plot_prod.grid(True)

#~ timeFormat = FormatStrFormatter("%.3g")
#~ plot_prod.xaxis.set_major_formatter(timeFormat)

plot_prod.legend()

# Plot destruction
plot_dest = fig.add_subplot(2, 1, 2)

plot_dest.set_xscale("log")
plot = plot_dest.plot

for reaction in range(nb_dest_reactions):
  plot(time, destruction_fraction[reaction], label="ID=%s" % reaction_dest_IDs[reaction])

plot_dest.set_xlabel("Time [years]")
plot_dest.set_ylabel("Destruction fraction [%]")
plot_dest.grid(True)

plot_dest.legend()


#~ percentage_format = FormatStrFormatter("%3.1f")
percentage_format = ScalarFormatter(useOffset=False)

plot_prod.yaxis.set_major_formatter(percentage_format)
plot_dest.yaxis.set_major_formatter(percentage_format)



nom_fichier_plot = "major_reaction_%s" % species_name
fig.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)

#~ pdb.set_trace()

pl.show()



