#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import git_infos
import os
import sys
import glob
import string

LOG_NAME = 'compilation.log'
SIMULATION_FOLDER = "tests"
GNUPLOT_EXTENSION = "gnuplot"

VISUALISEUR = {"svg":"gthumb", "png":"gthumb", "pdf":"gv", "jpg":"gthumb"}

# Parameters
debug = False
gdb = False
profiling = False
force = False # To force the compilation of every module
indice_script = None # to allow default actions if nothing is given in parameters
generate_all = False

isProblem = False
problem_message = """This script will compile then execute the fortran program for 
unitary tests. Theses tests are divided into routines that generate a data file and a corresponding gnuplot file. 
Gnuplot generate *.pdf (usually) that are automatically displayed by this script. But everything is available 
in the sub-folder 'tests' of the Git repository (that is automatically created)

The script can take various arguments:
(no spaces between the key and the values, only separated by '=')
 * help : display a little help message on HOW to use various options
 * force : To force the compilation of every module even those not modified
 * debug : [%s] activate debug options
 * script=1 to avoid be asked what test we want. (index refer to a test, 'all' if we want to execute them all)
 * gdb : [%s] activate options for gdb
 * profiling : [%s] activate options for profiling
 
 Example : 
 unitary_tests.py gdb""" % (debug, gdb, profiling)

def run(commande):
  """lance une commande qui sera typiquement soit une liste, soit une 
  commande seule. La fonction renvoit un tuple avec la sortie, 
  l'erreur et le code de retour"""
  if (type(commande)==list):
    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  elif (type(commande)==str):
    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  else:
    raise TypeError("The command is neither a string nor a list.")
  (process_stdout, process_stderr) = process.communicate()
  returncode = process.poll()
  # there is .poll() or .wait() but I don't remember the difference. For some kind of things, one of the two was not working
  return (process_stdout, process_stderr, returncode)
  
def enumeration(liste):
  """liste, un élément par ligne, les éléments de l'argument
  avec les indices associés"""
  output_list = []
  for (indice,element) in enumerate(liste):
    txt = "%2i : %-34s" % (indice, os.path.splitext(element)[0])
    output_list.append(txt)
  
  print(" ".join(output_list))

def getOutput(scriptname):
  """function that return the name of the output file of the given gnuplot script, if there is any"""
  
  f=open(scriptname,'r')
  lines = f.readlines()
  f.close()
  
  output = None
  
  for line in reversed(lines):
    words = line.split()
    if (len(words) > 1):
      if (words[0] == "set") and (words[1] == "output"):
        output = words[2][1:-1]
        break
        
  # in case of output file printed not in the directory of the .dat file
  temp = output.split("/")
  output = temp[-1]
  
  return output
  
# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 'debug'):
    debug = True
  elif (key == 'force'):
    force = True
  elif (key == 'gdb'):
    gdb = True
  elif (key == 'script'):
    if (value == 'all'):
      generate_all = True
      indice_script = 0
    else:
      generate_all = False
      indice_script = int(value)
  elif (key == 'profiling'):
    profiling = True
  elif (key == 'help'):
    isProblem = True
  else:
    print("the key '"+key+"' does not match")
    isProblem = True



if isProblem:
  print(problem_message)
  exit()

compilation_command = "./Makefile.py name=unitary_tests.f90"

if debug:
  compilation_command += " debug"

if force:
  compilation_command += " force"

if profiling:
  compilation_command += " profiling"

print(compilation_command)
(process_stdout, process_stderr, returncode) = run(compilation_command)
print(process_stdout)

if (returncode != 0):
  print(process_stderr)

# We run the unitary tests
os.chdir(SIMULATION_FOLDER)

# We delete the old gnuplot scripts
run("rm *.gnuplot")

test_command = "../unitary_tests"
print(test_command)
(process_stdout, process_stderr, returncode) = run(test_command)
print(process_stdout)

if (returncode != 0):
  print(process_stderr)

# We run the gnuplot scripts we wants
scripts = glob.glob("*."+GNUPLOT_EXTENSION)
scripts.sort()
nb_scripts = len(scripts)

# If the action is not given in parameters, we ask the user through standard input
if (indice_script == None):
  enumeration(scripts)

  while not(0<=indice_script<nb_scripts):
    try:
      txt_input = raw_input("What test do you want to display? (0-"+str(nb_scripts-1)+" ; 'all' treat them all ; 'l' display list again)\n")
      
      indice_script = int(txt_input)
    except ValueError:
      if (txt_input == 'l'):
        enumeration(scripts)
      elif (txt_input == 'all'):
        generate_all = True
        indice_script = 0
      else:
        print("The parameter must be between 0 and %i" % (nb_scripts-1))

if not(generate_all):
  indexes = [indice_script]
else:
  indexes = range(nb_scripts)

# We run gnuplot for every script needed. 
for index in indexes:
  script = scripts[index]
  sys.stdout.write("Running %s...\r" % script)
  sys.stdout.flush()
  (stdout, stderr, returncode) = run("gnuplot "+script)
  sys.stdout.write("Running %s...OK\n" % script)
  sys.stdout.flush()

# We only display the output file for the last script (in order to display only one graphic if we choose to run them all
output_file = getOutput(script)
if (output_file != None):
  output_extension = os.path.splitext(output_file)[1][1:] # [1:] is here to get rid of the prefixed dot  "." of the extension

  run(VISUALISEUR[output_extension]+" "+output_file)
else:
  print("Warning: The current gnuplot file has no output file")
